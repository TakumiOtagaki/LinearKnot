#include "scfg2/proposal/hotspot_wavefront.hh"

#include "scfg2/types.hh"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::proposal {
namespace {

constexpr std::uint32_t kNoTrace = std::numeric_limits<std::uint32_t>::max();

struct TraceNode {
  cand_pos_t i = -1;
  cand_pos_t j = -1;
  std::uint32_t parent = kNoTrace;
};

struct Candidate {
  cand_pos_t i = -1;
  cand_pos_t j = -1;
  double physical_score = 0.0;
  double priority = 0.0;
  std::size_t pair_count = 0;
  std::uint32_t parent_trace = kNoTrace;
  cand_pos_t seed_outer_i = -1;
  cand_pos_t seed_outer_j = -1;
  bool is_seed = false;
  std::uint64_t origin = 0;
};

struct State {
  cand_pos_t i = -1;
  cand_pos_t j = -1;
  double physical_score = 0.0;
  std::size_t pair_count = 0;
  std::uint32_t trace = kNoTrace;
  std::uint64_t origin = 0;
};

struct Completion {
  double physical_score = 0.0;
  std::uint32_t trace = kNoTrace;
  std::uint64_t origin = 0;
};

struct CompletionMinHeap {
  bool operator()(const Completion &lhs, const Completion &rhs) const {
    if (lhs.physical_score != rhs.physical_score) {
      return lhs.physical_score > rhs.physical_score;
    }
    return lhs.origin > rhs.origin;
  }
};

std::uint64_t mix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

std::uint64_t combine_hash(std::uint64_t seed, std::uint64_t value) {
  return mix64(seed ^ mix64(value));
}

double signed_unit(std::uint64_t value) {
  constexpr double denominator =
      static_cast<double>(std::uint64_t{1} << 53U);
  const double unit =
      static_cast<double>(mix64(value) >> 11U) / denominator;
  return 2.0 * unit - 1.0;
}

double diversity_jitter(std::uint64_t rank_hash,
                        const WavefrontHotspotConfig &config) {
  if (config.diversity_mode == HotspotDiversityMode::None) return 0.0;
  return config.diversity_priority_jitter_amplitude *
         signed_unit(rank_hash);
}

int base_code(char base) {
  switch (base) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'U':
      return 3;
    default:
      return -1;
  }
}

bool normalize_sequence(const std::string &sequence,
                        std::string *normalized,
                        std::string *error) {
  normalized->clear();
  normalized->reserve(sequence.size());
  for (char raw : sequence) {
    char base = static_cast<char>(
        std::toupper(static_cast<unsigned char>(raw)));
    if (base == 'T') base = 'U';
    if (base_code(base) < 0) {
      if (error) *error = "sequence contains a non-ACGUT symbol";
      return false;
    }
    normalized->push_back(base);
  }
  if (normalized->empty()) {
    if (error) *error = "sequence must not be empty";
    return false;
  }
  return true;
}

double pair_score(char left,
                  char right,
                  const WavefrontHotspotConfig &config) {
  if ((left == 'G' && right == 'C') ||
      (left == 'C' && right == 'G')) {
    return config.gc_pair_score;
  }
  if ((left == 'A' && right == 'U') ||
      (left == 'U' && right == 'A')) {
    return config.au_pair_score;
  }
  if (config.allow_gu &&
      ((left == 'G' && right == 'U') ||
       (left == 'U' && right == 'G'))) {
    return config.gu_pair_score;
  }
  return -std::numeric_limits<double>::infinity();
}

bool can_pair(char left,
              char right,
              const WavefrontHotspotConfig &config) {
  return std::isfinite(pair_score(left, right, config));
}

std::size_t encode_kmer(const std::string &sequence,
                        std::size_t start,
                        std::size_t length) {
  std::size_t code = 0;
  for (std::size_t offset = 0; offset < length; ++offset) {
    code = (code << 2U) |
           static_cast<std::size_t>(base_code(sequence[start + offset]));
  }
  return code;
}

void append_compatible_left_codes(char right,
                                  bool allow_gu,
                                  std::vector<int> *codes) {
  codes->clear();
  switch (right) {
    case 'A':
      codes->push_back(base_code('U'));
      break;
    case 'C':
      codes->push_back(base_code('G'));
      break;
    case 'G':
      codes->push_back(base_code('C'));
      if (allow_gu) codes->push_back(base_code('U'));
      break;
    case 'U':
      codes->push_back(base_code('A'));
      if (allow_gu) codes->push_back(base_code('G'));
      break;
  }
}

void compatible_left_kmer_codes(const std::string &sequence,
                                std::size_t right_end,
                                std::size_t length,
                                bool allow_gu,
                                std::vector<std::size_t> *out) {
  out->assign(1, 0);
  std::vector<int> allowed;
  for (std::size_t offset = 0; offset < length; ++offset) {
    append_compatible_left_codes(
        sequence[right_end - offset], allow_gu, &allowed);
    std::vector<std::size_t> next;
    next.reserve(out->size() * allowed.size());
    for (std::size_t prefix : *out) {
      for (int code : allowed) {
        next.push_back((prefix << 2U) | static_cast<std::size_t>(code));
      }
    }
    *out = std::move(next);
  }
}

std::vector<std::size_t> selected_global_ranks(std::size_t total,
                                               std::size_t cap,
                                               std::uint64_t salt) {
  std::vector<std::size_t> ranks;
  if (total <= cap) {
    ranks.reserve(total);
    for (std::size_t rank = 0; rank < total; ++rank) {
      ranks.push_back(rank);
    }
    return ranks;
  }
  ranks.reserve(cap);
  const auto append_unique = [&](std::size_t rank,
                                 std::vector<std::size_t> *values) {
    if (values->size() >= cap) return;
    if (std::find(values->begin(), values->end(), rank) == values->end()) {
      values->push_back(rank);
    }
  };
  append_unique(0, &ranks);
  append_unique(total - 1, &ranks);
  for (std::size_t attempt = 0;
       ranks.size() < cap && attempt < cap * 8;
       ++attempt) {
    append_unique(
        static_cast<std::size_t>(
            mix64(combine_hash(salt, attempt)) % total),
        &ranks);
  }
  for (std::size_t slot = 1; ranks.size() < cap; ++slot) {
    append_unique((slot * (total - 1)) / cap, &ranks);
  }
  std::sort(ranks.begin(), ranks.end());
  return ranks;
}

std::string pair_map_key(const std::vector<api::BasePair> &pairs) {
  std::string key;
  key.reserve(pairs.size() * 16);
  for (const auto &pair : pairs) {
    key += std::to_string(pair.i);
    key.push_back(':');
    key += std::to_string(pair.j);
    key.push_back(';');
  }
  return key;
}

// Keep existing decimal pair-map identity and order. Map nodes/keys remain
// valid throughout sorting; only values move, after the pointer sort finishes.
using ProposalMap = std::unordered_map<std::string, HotspotProposal>;

std::vector<HotspotProposal> ranked_proposals(ProposalMap &proposals,
                                             std::size_t spare = 0) {
  std::vector<ProposalMap::value_type *> ranked;
  ranked.reserve(proposals.size());
  for (auto &entry : proposals) ranked.push_back(&entry);
  std::sort(ranked.begin(), ranked.end(), [](const auto *lhs, const auto *rhs) {
    const auto &a = lhs->second;
    const auto &b = rhs->second;
    if (a.heuristic_score != b.heuristic_score)
      return a.heuristic_score > b.heuristic_score;
    if (lhs->first != rhs->first) return lhs->first < rhs->first;
    if (a.pass != b.pass) return a.pass < b.pass;
    return a.origin < b.origin;
  });
  std::vector<HotspotProposal> out;
  out.reserve(proposals.size() + spare);
  for (auto *entry : ranked) out.push_back(std::move(entry->second));
  return out;
}

std::vector<api::BasePair> traceback_pairs(
    std::uint32_t trace,
    const std::vector<TraceNode> &arena) {
  std::vector<api::BasePair> pairs;
  while (trace != kNoTrace) {
    const auto &node = arena.at(trace);
    pairs.push_back(api::BasePair{
        static_cast<int>(node.i), static_cast<int>(node.j)});
    trace = node.parent;
  }
  // Every seed/extension increases i strictly; following parents reverses it.
  std::reverse(pairs.begin(), pairs.end());
  return pairs;
}

std::uint32_t append_trace(std::vector<TraceNode> *arena,
                           cand_pos_t i,
                           cand_pos_t j,
                           std::uint32_t parent) {
  arena->push_back(TraceNode{i, j, parent});
  return static_cast<std::uint32_t>(arena->size() - 1);
}

bool candidate_better(const Candidate &lhs, const Candidate &rhs) {
  if (lhs.priority != rhs.priority) return lhs.priority > rhs.priority;
  if (lhs.physical_score != rhs.physical_score) {
    return lhs.physical_score > rhs.physical_score;
  }
  if (lhs.pair_count != rhs.pair_count) {
    return lhs.pair_count > rhs.pair_count;
  }
  return lhs.origin < rhs.origin;
}

std::vector<HotspotProposal> run_one_pass(
    const std::string &sequence,
    const WavefrontHotspotConfig &config,
    std::size_t pass,
    WavefrontHotspotStats *stats) {
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto k = static_cast<cand_pos_t>(config.seed_pairs);
  const std::size_t code_count =
      std::size_t{1} << (2U * config.seed_pairs);
  std::vector<std::vector<cand_pos_t>> occurrences(code_count);
  for (cand_pos_t start = 0; start + k <= n; ++start) {
    occurrences[encode_kmer(
                    sequence,
                    static_cast<std::size_t>(start),
                    config.seed_pairs)]
        .push_back(start);
  }

  std::vector<std::size_t> eligible_count(code_count, 0);
  std::vector<std::vector<Candidate>> buckets(
      static_cast<std::size_t>(n));
  std::vector<std::size_t> compatible_codes;

  for (cand_pos_t outer_j = k - 1; outer_j < n; ++outer_j) {
    const cand_pos_t max_left =
        outer_j - 2 * k + 1 -
        static_cast<cand_pos_t>(config.minimum_loop_unpaired);
    if (max_left < 0) continue;
    compatible_left_kmer_codes(
        sequence,
        static_cast<std::size_t>(outer_j),
        config.seed_pairs,
        config.allow_gu,
        &compatible_codes);

    std::size_t total = 0;
    for (std::size_t code : compatible_codes) {
      auto &count = eligible_count[code];
      const auto &positions = occurrences[code];
      while (count < positions.size() && positions[count] <= max_left) {
        ++count;
      }
      total += count;
    }
    if (total == 0) continue;

    const std::size_t diversity_pass =
        config.diversity_mode == HotspotDiversityMode::None ? 0 : pass;
    std::uint64_t salt = combine_hash(config.seed, diversity_pass);
    salt = combine_hash(salt, static_cast<std::uint64_t>(outer_j));
    const auto ranks =
        selected_global_ranks(total, config.anchor_cap, salt);
    for (std::size_t global_rank : ranks) {
      std::size_t code_index = 0;
      while (global_rank >= eligible_count[compatible_codes[code_index]]) {
        global_rank -= eligible_count[compatible_codes[code_index]];
        ++code_index;
      }
      const cand_pos_t outer_i =
          occurrences[compatible_codes[code_index]][global_rank];
      double score = config.stack_bonus *
                     static_cast<double>(config.seed_pairs - 1);
      for (cand_pos_t offset = 0; offset < k; ++offset) {
        score += pair_score(
            sequence[static_cast<std::size_t>(outer_i + offset)],
            sequence[static_cast<std::size_t>(outer_j - offset)],
            config);
      }
      const cand_pos_t inner_i = outer_i + k - 1;
      const cand_pos_t inner_j = outer_j - k + 1;
      std::uint64_t origin = combine_hash(config.seed, pass);
      origin = combine_hash(origin, static_cast<std::uint64_t>(outer_i));
      origin = combine_hash(origin, static_cast<std::uint64_t>(outer_j));
      const double priority_jitter = diversity_jitter(origin, config);
      buckets[static_cast<std::size_t>(inner_j)].push_back(Candidate{
          inner_i,
          inner_j,
          score,
          score + priority_jitter,
          config.seed_pairs,
          kNoTrace,
          outer_i,
          outer_j,
          true,
          origin,
      });
      ++stats->anchors_retained;
    }
  }

  std::vector<TraceNode> arena;
  arena.reserve(sequence.size() * config.beam_width);
  const std::size_t completion_cap =
      std::max<std::size_t>(config.max_results, config.max_results * 8);
  std::priority_queue<
      Completion, std::vector<Completion>, CompletionMinHeap>
      completions;

  constexpr std::array<std::pair<cand_pos_t, cand_pos_t>, 4> moves = {
      std::pair<cand_pos_t, cand_pos_t>{1, 1},
      std::pair<cand_pos_t, cand_pos_t>{2, 1},
      std::pair<cand_pos_t, cand_pos_t>{1, 2},
      std::pair<cand_pos_t, cand_pos_t>{2, 2},
  };

  for (cand_pos_t j = n - 1; j >= 0; --j) {
    ++stats->boundaries;
    auto &incoming = buckets[static_cast<std::size_t>(j)];
    std::unordered_map<cand_pos_t, std::size_t> by_left;
    std::vector<Candidate> merged;
    merged.reserve(incoming.size());
    for (const auto &candidate : incoming) {
      const auto found = by_left.find(candidate.i);
      if (found == by_left.end()) {
        by_left.emplace(candidate.i, merged.size());
        merged.push_back(candidate);
      } else if (candidate_better(candidate, merged[found->second])) {
        merged[found->second] = candidate;
      }
    }
    std::sort(merged.begin(), merged.end(), candidate_better);
    if (merged.size() > config.beam_width) {
      merged.resize(config.beam_width);
    }

    std::vector<State> retained;
    retained.reserve(merged.size());
    for (const auto &candidate : merged) {
      std::uint32_t trace = candidate.parent_trace;
      if (candidate.is_seed) {
        trace = kNoTrace;
        for (cand_pos_t offset = 0; offset < k; ++offset) {
          trace = append_trace(
              &arena,
              candidate.seed_outer_i + offset + 1,
              candidate.seed_outer_j - offset + 1,
              trace);
        }
      } else {
        trace = append_trace(
            &arena, candidate.i + 1, candidate.j + 1, trace);
      }
      retained.push_back(State{
          candidate.i,
          candidate.j,
          candidate.physical_score,
          candidate.pair_count,
          trace,
          candidate.origin,
      });
      ++stats->states_retained;
      if (candidate.pair_count >= config.minimum_pairs &&
          candidate.physical_score >= config.minimum_score) {
        completions.push(Completion{
            candidate.physical_score, trace, candidate.origin});
        ++stats->completed_candidates;
        if (completions.size() > completion_cap) completions.pop();
      }
    }

    for (const auto &state : retained) {
      for (std::size_t move_index = 0; move_index < moves.size();
           ++move_index) {
        const cand_pos_t next_i = state.i + moves[move_index].first;
        const cand_pos_t next_j = state.j - moves[move_index].second;
        if (next_i < 0 || next_j < 0 || next_i >= next_j ||
            next_j - next_i - 1 <
                static_cast<cand_pos_t>(config.minimum_loop_unpaired) ||
            !can_pair(
                sequence[static_cast<std::size_t>(next_i)],
                sequence[static_cast<std::size_t>(next_j)],
                config)) {
          continue;
        }
        double delta = pair_score(
            sequence[static_cast<std::size_t>(next_i)],
            sequence[static_cast<std::size_t>(next_j)],
            config);
        if (move_index == 0) {
          delta += config.stack_bonus;
        } else if (move_index == 1 || move_index == 2) {
          delta -= config.bulge_penalty;
        } else {
          delta -= config.internal_loop_penalty;
        }
        std::uint64_t rank_hash = combine_hash(state.origin, move_index);
        rank_hash = combine_hash(rank_hash, static_cast<std::uint64_t>(next_i));
        rank_hash = combine_hash(rank_hash, static_cast<std::uint64_t>(next_j));
        const double physical = state.physical_score + delta;
        buckets[static_cast<std::size_t>(next_j)].push_back(Candidate{
            next_i,
            next_j,
            physical,
            physical + diversity_jitter(rank_hash, config),
            state.pair_count + 1,
            state.trace,
            -1,
            -1,
            false,
            state.origin,
        });
        ++stats->successor_candidates;
      }
    }
    std::vector<Candidate>().swap(incoming);
  }

  std::unordered_map<std::string, HotspotProposal> best_by_pairs;
  while (!completions.empty()) {
    const auto completion = completions.top();
    completions.pop();
    auto pairs = traceback_pairs(completion.trace, arena);
    const std::string key = pair_map_key(pairs);
    HotspotProposal proposal{
        completion.physical_score,
        std::move(pairs),
        pass,
        completion.origin,
    };
    const auto found = best_by_pairs.find(key);
    if (found == best_by_pairs.end() ||
        proposal.heuristic_score > found->second.heuristic_score) {
      best_by_pairs[key] = std::move(proposal);
    }
  }

  auto out = ranked_proposals(best_by_pairs);
  if (out.size() > config.max_results) out.resize(config.max_results);
  return out;
}

}  // namespace

WavefrontHotspotResult propose_hotspots_wavefront(
    const std::string &sequence,
    const WavefrontHotspotConfig &config) {
  WavefrontHotspotResult out;
  if (config.anchor_cap == 0 || config.beam_width == 0 ||
      config.passes == 0 || config.max_results == 0) {
    out.error =
        "anchor cap, beam width, passes, and max results must be positive";
    return out;
  }
  if (config.seed_pairs != 3) {
    out.error = "the initial wavefront proposer requires seed_pairs == 3";
    return out;
  }
  if (config.minimum_pairs < config.seed_pairs) {
    out.error = "minimum_pairs must be at least seed_pairs";
    return out;
  }
  if (!std::isfinite(config.diversity_priority_jitter_amplitude) ||
      config.diversity_priority_jitter_amplitude < 0.0) {
    out.error =
        "diversity priority-jitter amplitude must be finite and nonnegative";
    return out;
  }

  std::string normalized;
  if (!normalize_sequence(sequence, &normalized, &out.error)) return out;

  std::unordered_map<std::string, HotspotProposal> best_by_pairs;
  for (std::size_t pass = 0; pass < config.passes; ++pass) {
    auto proposals =
        run_one_pass(normalized, config, pass, &out.stats);
    for (auto &proposal : proposals) {
      const std::string key = pair_map_key(proposal.scaffold_pairs);
      const auto found = best_by_pairs.find(key);
      if (found == best_by_pairs.end() ||
          proposal.heuristic_score > found->second.heuristic_score) {
        best_by_pairs[key] = std::move(proposal);
      }
    }
  }

  out.proposals = ranked_proposals(best_by_pairs, 1);

  // The empty scaffold is a required baseline and is kept independently of
  // hotspot score.
  const std::size_t nonempty_cap =
      config.max_results > 0 ? config.max_results - 1 : 0;
  if (out.proposals.size() > nonempty_cap) {
    out.proposals.resize(nonempty_cap);
  }
  out.proposals.push_back(HotspotProposal{});
  out.ok = true;
  return out;
}

}  // namespace cparty::scfg2::proposal

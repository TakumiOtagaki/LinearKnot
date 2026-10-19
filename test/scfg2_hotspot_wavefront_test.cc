#include "scfg2/proposal/hotspot_wavefront.hh"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

using cparty::scfg2::proposal::HotspotProposal;
using cparty::scfg2::api::BasePair;

bool same_pairs(const std::vector<BasePair> &lhs,
                const std::vector<BasePair> &rhs) {
  if (lhs.size() != rhs.size()) return false;
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    if (lhs[index].i != rhs[index].i || lhs[index].j != rhs[index].j) {
      return false;
    }
  }
  return true;
}

bool same_proposals(const std::vector<HotspotProposal> &lhs,
                    const std::vector<HotspotProposal> &rhs) {
  if (lhs.size() != rhs.size()) return false;
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    if (lhs[index].heuristic_score != rhs[index].heuristic_score ||
        lhs[index].pass != rhs[index].pass ||
        lhs[index].origin != rhs[index].origin ||
        !same_pairs(
            lhs[index].scaffold_pairs, rhs[index].scaffold_pairs)) {
      return false;
    }
  }
  return true;
}

bool canonical(char left, char right, bool allow_gu) {
  if ((left == 'G' && right == 'C') ||
      (left == 'C' && right == 'G') ||
      (left == 'A' && right == 'U') ||
      (left == 'U' && right == 'A')) {
    return true;
  }
  return allow_gu &&
         ((left == 'G' && right == 'U') ||
          (left == 'U' && right == 'G'));
}

bool valid_nested_proposal(const std::string &sequence,
                           const HotspotProposal &proposal,
                           bool allow_gu) {
  std::set<int> endpoints;
  for (const auto &pair : proposal.scaffold_pairs) {
    if (!(1 <= pair.i && pair.i < pair.j &&
          pair.j <= static_cast<int>(sequence.size())) ||
        !endpoints.insert(pair.i).second ||
        !endpoints.insert(pair.j).second ||
        !canonical(
            sequence[static_cast<std::size_t>(pair.i - 1)],
            sequence[static_cast<std::size_t>(pair.j - 1)],
            allow_gu)) {
      return false;
    }
  }
  for (std::size_t lhs = 0; lhs < proposal.scaffold_pairs.size(); ++lhs) {
    const auto &a = proposal.scaffold_pairs[lhs];
    for (std::size_t rhs = lhs + 1;
         rhs < proposal.scaffold_pairs.size();
         ++rhs) {
      const auto &b = proposal.scaffold_pairs[rhs];
      if ((a.i < b.i && b.i < a.j && a.j < b.j) ||
          (b.i < a.i && a.i < b.j && b.j < a.j)) {
        return false;
      }
    }
  }
  return true;
}

bool contains_pairs(const std::vector<HotspotProposal> &proposals,
                    const std::vector<BasePair> &expected) {
  return std::any_of(
      proposals.begin(), proposals.end(), [&](const auto &proposal) {
        return same_pairs(proposal.scaffold_pairs, expected);
      });
}

std::size_t brute_force_seed_count(
    const std::string &sequence,
    const cparty::scfg2::proposal::WavefrontHotspotConfig &config) {
  const auto n = static_cast<int>(sequence.size());
  const auto k = static_cast<int>(config.seed_pairs);
  std::size_t count = 0;
  for (int outer_j = k - 1; outer_j < n; ++outer_j) {
    for (int outer_i = 0; outer_i + k <= outer_j - k + 1; ++outer_i) {
      const int inner_i = outer_i + k - 1;
      const int inner_j = outer_j - k + 1;
      if (inner_j - inner_i - 1 <
          static_cast<int>(config.minimum_loop_unpaired)) {
        continue;
      }
      bool valid = true;
      for (int offset = 0; offset < k; ++offset) {
        valid = valid &&
                canonical(
                    sequence[static_cast<std::size_t>(outer_i + offset)],
                    sequence[static_cast<std::size_t>(outer_j - offset)],
                    config.allow_gu);
      }
      if (valid) ++count;
    }
  }
  return count;
}

}  // namespace

int main() {
  cparty::scfg2::proposal::WavefrontHotspotConfig config;
  config.anchor_cap = 8;
  config.beam_width = 16;
  config.passes = 3;
  config.max_results = 12;
  config.seed = 20260727;

  const std::string sequence = "GGGAAACCC";
  const auto first =
      cparty::scfg2::proposal::propose_hotspots_wavefront(sequence, config);
  const auto second =
      cparty::scfg2::proposal::propose_hotspots_wavefront(sequence, config);
  const std::vector<BasePair> expected = {{1, 9}, {2, 8}, {3, 7}};
  if (!first.ok || !second.ok ||
      !contains_pairs(first.proposals, expected) ||
      !same_proposals(first.proposals, second.proposals) ||
      first.proposals.empty() ||
      !first.proposals.back().scaffold_pairs.empty()) {
    std::cerr << "basic/deterministic wavefront proposal failed: "
              << first.error << '\n';
    return EXIT_FAILURE;
  }

  for (const auto &proposal : first.proposals) {
    if (!valid_nested_proposal(sequence, proposal, config.allow_gu)) {
      std::cerr << "invalid scaffold proposal\n";
      return EXIT_FAILURE;
    }
  }

  auto oracle_config = config;
  oracle_config.anchor_cap = 1024;
  oracle_config.beam_width = 1024;
  oracle_config.passes = 1;
  oracle_config.diversity_mode =
      cparty::scfg2::proposal::HotspotDiversityMode::None;
  const std::string oracle_sequence = "GGGCGAAAUCGCCC";
  const auto oracle =
      cparty::scfg2::proposal::propose_hotspots_wavefront(
          oracle_sequence, oracle_config);
  if (!oracle.ok ||
      oracle.stats.anchors_retained !=
          brute_force_seed_count(oracle_sequence, oracle_config)) {
    std::cerr << "small-n anchor oracle failed: expected="
              << brute_force_seed_count(oracle_sequence, oracle_config)
              << " observed=" << oracle.stats.anchors_retained << '\n';
    return EXIT_FAILURE;
  }

  std::string repetitive(36, 'G');
  repetitive.append(36, 'C');
  const auto bounded =
      cparty::scfg2::proposal::propose_hotspots_wavefront(repetitive, config);
  const std::size_t n = repetitive.size();
  if (!bounded.ok ||
      bounded.proposals.size() > config.max_results ||
      bounded.stats.boundaries != config.passes * n ||
      bounded.stats.anchors_retained >
          config.passes * n * config.anchor_cap ||
      bounded.stats.states_retained >
          config.passes * n * config.beam_width ||
      bounded.stats.successor_candidates >
          4 * bounded.stats.states_retained) {
    std::cerr << "bounded-complexity contract failed\n";
    return EXIT_FAILURE;
  }
  for (const auto &proposal : bounded.proposals) {
    if (!valid_nested_proposal(repetitive, proposal, config.allow_gu)) {
      std::cerr << "repetitive-input proposal is invalid\n";
      return EXIT_FAILURE;
    }
  }

  auto single_anchor_config = config;
  single_anchor_config.anchor_cap = 1;
  single_anchor_config.passes = 1;
  const auto single_anchor =
      cparty::scfg2::proposal::propose_hotspots_wavefront(
          repetitive, single_anchor_config);
  if (!single_anchor.ok ||
      single_anchor.stats.anchors_retained > repetitive.size()) {
    std::cerr << "single-anchor cap contract failed\n";
    return EXIT_FAILURE;
  }

  auto no_diversity_one_pass = config;
  no_diversity_one_pass.passes = 1;
  no_diversity_one_pass.diversity_mode =
      cparty::scfg2::proposal::HotspotDiversityMode::None;
  auto no_diversity_many_passes = no_diversity_one_pass;
  no_diversity_many_passes.passes = 4;
  const auto one_pass =
      cparty::scfg2::proposal::propose_hotspots_wavefront(
          repetitive, no_diversity_one_pass);
  const auto redundant_passes =
      cparty::scfg2::proposal::propose_hotspots_wavefront(
          repetitive, no_diversity_many_passes);
  if (!one_pass.ok || !redundant_passes.ok ||
      !same_proposals(one_pass.proposals, redundant_passes.proposals)) {
    std::cerr << "disabled-diversity pass contract failed\n";
    return EXIT_FAILURE;
  }

  const std::string gu_sequence = "GGGAAAUUU";
  auto no_gu_config = config;
  no_gu_config.allow_gu = false;
  const auto without_gu =
      cparty::scfg2::proposal::propose_hotspots_wavefront(
          gu_sequence, no_gu_config);
  const auto with_gu =
      cparty::scfg2::proposal::propose_hotspots_wavefront(
          gu_sequence, config);
  if (!without_gu.ok || !with_gu.ok ||
      contains_pairs(without_gu.proposals, expected) ||
      !contains_pairs(with_gu.proposals, expected)) {
    std::cerr << "GU policy failed\n";
    return EXIT_FAILURE;
  }

  auto invalid_config = config;
  invalid_config.anchor_cap = 0;
  const auto invalid =
      cparty::scfg2::proposal::propose_hotspots_wavefront(
          sequence, invalid_config);
  if (invalid.ok || invalid.error.empty()) {
    std::cerr << "invalid-config contract failed\n";
    return EXIT_FAILURE;
  }

  if (const char *long_n_env = std::getenv("SCFG2_HOTSPOT_LONG_N")) {
    const auto requested =
        static_cast<std::size_t>(std::strtoull(long_n_env, nullptr, 10));
    if (requested < 12) {
      std::cerr << "SCFG2_HOTSPOT_LONG_N must be at least 12\n";
      return EXIT_FAILURE;
    }
    std::string long_sequence(requested / 2, 'G');
    long_sequence.append(requested - requested / 2, 'C');
    auto long_config = config;
    long_config.passes = 1;
    long_config.max_results = 8;
    const auto long_result =
        cparty::scfg2::proposal::propose_hotspots_wavefront(
            long_sequence, long_config);
    if (!long_result.ok ||
        long_result.stats.anchors_retained >
            requested * long_config.anchor_cap ||
        long_result.stats.states_retained >
            requested * long_config.beam_width ||
        long_result.stats.successor_candidates >
            4 * long_result.stats.states_retained) {
      std::cerr << "long-input bounded-complexity contract failed\n";
      return EXIT_FAILURE;
    }
    std::cout << " long_n=" << requested
              << " long_anchors=" << long_result.stats.anchors_retained
              << " long_retained=" << long_result.stats.states_retained
              << '\n';
  }

  std::cout << "scfg2_hotspot_wavefront=ok"
            << " anchors=" << bounded.stats.anchors_retained
            << " retained=" << bounded.stats.states_retained
            << " successors=" << bounded.stats.successor_candidates
            << '\n';
  return EXIT_SUCCESS;
}

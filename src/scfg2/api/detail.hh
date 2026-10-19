#ifndef CPARTY_SCFG2_API_DETAIL_HH
#define CPARTY_SCFG2_API_DETAIL_HH

#include "scfg2/api/exact.hh"
#include "scfg2/geometry/scaffold_input.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

namespace cparty::scfg2::api::detail {

inline std::vector<engine::TracebackPair> to_traceback_pairs(
    const std::vector<BasePair> &pairs) {
  std::vector<engine::TracebackPair> out;
  out.reserve(pairs.size());
  for (const auto &pair : pairs) {
    out.push_back(engine::TracebackPair{
        static_cast<cand_pos_t>(std::min(pair.i, pair.j)),
        static_cast<cand_pos_t>(std::max(pair.i, pair.j)),
        engine::TracebackPairLane::Square});
  }
  return out;
}

inline std::vector<BasePair> to_base_pairs(
    const std::vector<engine::TracebackPair> &pairs) {
  std::vector<BasePair> out;
  out.reserve(pairs.size());
  for (const auto &pair : pairs) {
    out.push_back(BasePair{static_cast<int>(pair.i), static_cast<int>(pair.j)});
  }
  std::sort(out.begin(), out.end(), [](const BasePair &lhs, const BasePair &rhs) {
    return lhs.i != rhs.i ? lhs.i < rhs.i : lhs.j < rhs.j;
  });
  return out;
}

// Pair coordinates lie in 1..n. Two stable counting passes preserve the
// public (i,j) order in O(number_of_pairs+n), including duplicate coordinates.
inline std::vector<BasePair> to_base_pairs_linear(
    const std::vector<engine::TracebackPair> &pairs, int n) {
  if (n < 0) throw std::invalid_argument("negative pair output sequence length");
  std::vector<BasePair> out;
  out.reserve(pairs.size());
  for (const auto &pair : pairs) {
    if (pair.i < 1 || pair.i > n || pair.j < 1 || pair.j > n)
      throw std::invalid_argument("pair output coordinate is outside the sequence");
    out.push_back({pair.i, pair.j});
  }
  if (out.size() < 2) return out;
  std::vector<BasePair> scratch(out.size());
  for (const bool by_i : {false, true}) {
    std::vector<std::size_t> offsets(static_cast<std::size_t>(n) + 2);
    for (const auto &pair : out) ++offsets[(by_i ? pair.i : pair.j) + 1];
    for (std::size_t k = 1; k < offsets.size(); ++k) offsets[k] += offsets[k - 1];
    for (const auto &pair : out) scratch[offsets[by_i ? pair.i : pair.j]++] = pair;
    out.swap(scratch);
  }
  return out;
}

template <class PairabilityContext>
std::string render_generated_structure(
    int n,
    const std::vector<engine::TracebackPair> &selected_pairs,
    const PairabilityContext &ctx) {
  std::vector<engine::TracebackPair> generated;
  generated.reserve(selected_pairs.size());
  for (const auto &pair : selected_pairs) {
    if (!engine::is_fixed_scaffold_pair(pair, ctx)) {
      generated.push_back(pair);
    }
  }
  engine::stable_unique_traceback_pairs(&generated);
  const auto lanes = engine::render_traceback_pairs(
      static_cast<cand_pos_t>(n), generated);
  std::string out = lanes.round;
  for (std::size_t pos = 0; pos < out.size(); ++pos) {
    if (lanes.square[pos] != '.') out[pos] = lanes.square[pos];
  }
  return out;
}

inline std::string compose_union_structure(const std::string &scaffold,
                                           const std::string &generated) {
  if (scaffold.size() != generated.size()) {
    throw std::invalid_argument(
        "scaffold and generated structures have different lengths");
  }
  std::string out = scaffold;
  for (std::size_t pos = 0; pos < generated.size(); ++pos) {
    if (generated[pos] == '.') continue;
    if (out[pos] != '.' && out[pos] != '_') {
      throw std::invalid_argument(
          "scaffold and generated structures reuse one endpoint");
    }
    out[pos] = generated[pos];
  }
  return out;
}

inline bool validate_pair_map(const std::vector<BasePair> &pairs,
                              int n,
                              std::string *error) {
  std::vector<int> partner(static_cast<std::size_t>(n + 1), -1);
  for (const auto &raw : pairs) {
    const int i = std::min(raw.i, raw.j);
    const int j = std::max(raw.i, raw.j);
    if (!(1 <= i && i < j && j <= n)) {
      if (error) {
        *error = "target base pair is outside the 1-based sequence domain";
      }
      return false;
    }
    if (partner[static_cast<std::size_t>(i)] >= 0 ||
        partner[static_cast<std::size_t>(j)] >= 0) {
      if (error) {
        *error = "target base-pair map uses one endpoint more than once";
      }
      return false;
    }
    partner[static_cast<std::size_t>(i)] = j;
    partner[static_cast<std::size_t>(j)] = i;
  }
  return true;
}

inline bool parse_dot_bracket_pairs(const std::string &structure,
                                    std::vector<BasePair> *pairs,
                                    std::string *error) {
  if (!pairs) {
    if (error) *error = "null target-pair output";
    return false;
  }
  pairs->clear();
  constexpr std::array<char, 4> opens = {'(', '[', '{', '<'};
  constexpr std::array<char, 4> closes = {')', ']', '}', '>'};
  std::array<std::vector<int>, 4> stacks;
  for (std::size_t pos = 0; pos < structure.size(); ++pos) {
    const char symbol = structure[pos];
    bool handled = false;
    for (std::size_t kind = 0; kind < opens.size(); ++kind) {
      if (symbol == opens[kind]) {
        stacks[kind].push_back(static_cast<int>(pos + 1));
        handled = true;
        break;
      }
      if (symbol == closes[kind]) {
        if (stacks[kind].empty()) {
          if (error) *error = "unbalanced target dot-bracket structure";
          return false;
        }
        pairs->push_back(
            BasePair{stacks[kind].back(), static_cast<int>(pos + 1)});
        stacks[kind].pop_back();
        handled = true;
        break;
      }
    }
    if (!handled && symbol != '.' && symbol != '_' && symbol != 'x' &&
        symbol != 'X') {
      if (error) *error = "unsupported target dot-bracket symbol";
      return false;
    }
  }
  for (const auto &stack : stacks) {
    if (!stack.empty()) {
      if (error) *error = "unbalanced target dot-bracket structure";
      return false;
    }
  }
  return true;
}

inline bool validate_session_inputs(const std::string &sequence,
                                    const std::string &scaffold,
                                    std::string *error) {
  if (sequence.empty()) {
    if (error) *error = "sequence must not be empty";
    return false;
  }
  if (sequence.size() != scaffold.size()) {
    if (error) *error = "scaffold structure length differs from sequence length";
    return false;
  }
  return geometry::validate_scaffold_input(scaffold, error);
}

}  // namespace cparty::scfg2::api::detail

#endif  // CPARTY_SCFG2_API_DETAIL_HH

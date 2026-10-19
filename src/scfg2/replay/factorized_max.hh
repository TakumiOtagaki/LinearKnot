#ifndef CPARTY_SCFG2_REPLAY_FACTORIZED_MAX_HH
#define CPARTY_SCFG2_REPLAY_FACTORIZED_MAX_HH

#include "scfg2/semirings/factorized_logsumexp.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cparty::scfg2::replay {

struct FactorizedMaximum {
  std::vector<double> scores;
  std::vector<std::uint32_t> choices;
  std::size_t terms_evaluated = 0;
  std::size_t factors_evaluated = 0;
};

// Derived consumer of the frozen inside forest. Its terms already include
// every forward admission/selection decision. Local thermodynamic weights
// decide that support, while this max-plus sweep uses additive node rewards.
// Child IDs precede parent IDs. Each term and factor is visited once.
inline FactorizedMaximum maximize_factorized_forest(
    const semirings::FactorizedInsideGraph &graph,
    const std::vector<double> &rewards) {
  if (rewards.size() != graph.states().size())
    throw std::invalid_argument("forest reward count differs from state count");
  FactorizedMaximum out;
  out.scores.resize(rewards.size(), -std::numeric_limits<double>::infinity());
  out.choices.resize(rewards.size(), semirings::invalid_outside_state);
  for (std::size_t node = 0; node < rewards.size(); ++node) {
    if (!std::isfinite(rewards[node]))
      throw std::invalid_argument("forest reward is non-finite");
    const auto &state = graph.states()[node];
    for (std::size_t t = state.first_term;
         t < std::size_t(state.first_term) + state.term_count; ++t) {
      const auto &term = graph.terms()[t];
      double score = rewards[node];
      ++out.terms_evaluated;
      for (std::size_t f = term.first_factor;
           f < std::size_t(term.first_factor) + term.factor_count; ++f) {
        const auto child = graph.factors()[f];
        if (child >= node) throw std::logic_error("forest is not topological");
        score += out.scores[child];
        ++out.factors_evaluated;
      }
      // Exact ties keep the first term in the frozen forward order.
      if (score > out.scores[node]) {
        out.scores[node] = score;
        out.choices[node] = static_cast<std::uint32_t>(t);
      }
    }
    if (!std::isfinite(out.scores[node]))
      throw std::runtime_error("forest maximum is non-finite");
  }
  return out;
}

// Saturating occurrence counts for the chosen derivation. Reverse ID order
// follows shared subtrees without expanding a DAG into an exponential tree.
inline std::vector<unsigned char> maximum_state_occurrences(
    const semirings::FactorizedInsideGraph &graph,
    const FactorizedMaximum &maximum, std::uint32_t root) {
  if (root >= graph.states().size() ||
      maximum.choices.size() != graph.states().size())
    throw std::invalid_argument("forest maximum root/choices are invalid");
  std::vector<unsigned char> counts(graph.states().size());
  counts[root] = 1;
  for (std::size_t cursor = std::size_t(root) + 1; cursor > 0; --cursor) {
    const auto node = cursor - 1;
    if (!counts[node]) continue;
    const auto &term = graph.terms().at(maximum.choices[node]);
    for (std::size_t f = term.first_factor;
         f < std::size_t(term.first_factor) + term.factor_count; ++f) {
      const auto child = graph.factors()[f];
      counts[child] = static_cast<unsigned char>(
          std::min(2u, unsigned(counts[child]) + counts[node]));
    }
  }
  return counts;
}

}  // namespace cparty::scfg2::replay

#endif

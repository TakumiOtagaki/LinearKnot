#ifndef CPARTY_SCFG2_ENGINE_OUTSIDE_MARGINALS_HH
#define CPARTY_SCFG2_ENGINE_OUTSIDE_MARGINALS_HH

#include "scfg2/types.hh"
#include "scfg2/grammar/nonterminal.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cparty::scfg2::engine {

struct LanePairMarginal {
  cand_pos_t i = -1;
  cand_pos_t j = -1;
  double round_probability = 0.0;
  double square_probability = 0.0;
  double union_probability = 0.0;
};

struct OutsideMarginalSummary {
  bool ok = false;
  std::string error;
  std::vector<LanePairMarginal> pairs;
  // One-based: index zero is unused.
  std::vector<double> unpaired_probability;
};

template <class InsideChart, class ItemPosterior, class PairabilityContext>
OutsideMarginalSummary collect_outside_marginals_by_item(
    const InsideChart &chart,
    ItemPosterior &&posterior,
    cand_pos_t n,
    const PairabilityContext &ctx,
    double tolerance = 1e-9) {
  std::vector<LanePairMarginal> contributions;
  for (cand_pos_t i = 1; i <= n; ++i) {
    const cand_pos_t j = ctx.pair_round(i);
    if (j > i) {
      contributions.push_back(LanePairMarginal{i, j, 1.0, 0.0, 1.0});
    }
  }
  for (const auto &[item, value] : chart.values()) {
    const bool round_owner =
        item.nonterminal == grammar::NonTerminal::V;
    const bool square_owner =
        item.nonterminal == grammar::NonTerminal::VP;
    if (!round_owner && !square_owner) continue;
    if (ctx.pair_partner(item.i) == item.j ||
        ctx.pair_partner(item.j) == item.i) {
      continue;
    }
    const std::optional<double> probability = posterior(item, value);
    if (!probability.has_value()) continue;
    contributions.push_back(LanePairMarginal{item.i, item.j,
        round_owner ? *probability : 0.0, square_owner ? *probability : 0.0, *probability});
  }

  // Two stable counting passes on bounded nucleotide coordinates give the
  // public lexicographic pair order in O(number_of_owners + n), without an
  // ordered-map lookup or a log(n) sort per retained owner.
  std::vector<LanePairMarginal> scratch(contributions.size());
  for (const bool by_i : {false, true}) {
    std::vector<std::size_t> offsets(static_cast<std::size_t>(n) + 2, 0);
    for (const auto &p : contributions) ++offsets[(by_i ? p.i : p.j) + 1];
    for (std::size_t k = 1; k < offsets.size(); ++k) offsets[k] += offsets[k - 1];
    for (const auto &p : contributions) scratch[offsets[by_i ? p.i : p.j]++] = p;
    contributions.swap(scratch);
  }
  std::vector<LanePairMarginal> by_pair;
  for (const auto &p : contributions) {
    if (!by_pair.empty() && by_pair.back().i == p.i && by_pair.back().j == p.j) {
      auto &merged = by_pair.back();
      merged.round_probability += p.round_probability;
      merged.square_probability += p.square_probability;
      merged.union_probability = merged.round_probability + merged.square_probability;
    } else {
      by_pair.push_back(p);
    }
  }

  OutsideMarginalSummary out;
  out.unpaired_probability.assign(static_cast<std::size_t>(n) + 1, 1.0);
  std::vector<double> paired(static_cast<std::size_t>(n) + 1, 0.0);
  for (auto &marginal : by_pair) {
    const auto valid_probability = [&](double value) {
      return std::isfinite(value) && value >= -tolerance &&
             value <= 1.0 + tolerance;
    };
    if (!valid_probability(marginal.round_probability) ||
        !valid_probability(marginal.square_probability) ||
        !valid_probability(marginal.union_probability)) {
      out.error = "outside pair probability is outside [0,1]";
      return out;
    }
    marginal.round_probability =
        std::clamp(marginal.round_probability, 0.0, 1.0);
    marginal.square_probability =
        std::clamp(marginal.square_probability, 0.0, 1.0);
    marginal.union_probability =
        std::clamp(marginal.union_probability, 0.0, 1.0);
    paired[static_cast<std::size_t>(marginal.i)] +=
        marginal.union_probability;
    paired[static_cast<std::size_t>(marginal.j)] +=
        marginal.union_probability;
    out.pairs.push_back(marginal);
  }
  for (cand_pos_t i = 1; i <= n; ++i) {
    const double probability = paired[static_cast<std::size_t>(i)];
    if (!std::isfinite(probability) || probability > 1.0 + tolerance) {
      out.error = "outside residue paired probability exceeds one";
      return out;
    }
    out.unpaired_probability[static_cast<std::size_t>(i)] =
        std::clamp(1.0 - probability, 0.0, 1.0);
  }
  out.ok = true;
  return out;
}

template <class TapeChart, class PairabilityContext>
OutsideMarginalSummary collect_outside_marginals(
    const TapeChart &chart, const std::vector<double> &adjoints,
    cand_pos_t n, const PairabilityContext &ctx, double tolerance = 1e-9) {
  return collect_outside_marginals_by_item(chart, [&](const auto &, const auto &value) -> std::optional<double> {
    if (value.node == semirings::invalid_tape_node || value.node >= adjoints.size()) return std::nullopt;
    return adjoints[value.node];
  }, n, ctx, tolerance);
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_OUTSIDE_MARGINALS_HH

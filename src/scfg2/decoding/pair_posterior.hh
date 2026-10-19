#ifndef CPARTY_SCFG2_DECODING_PAIR_POSTERIOR_HH
#define CPARTY_SCFG2_DECODING_PAIR_POSTERIOR_HH

#include "scfg2/engine/outside_marginals.hh"
#include "scfg2/engine/traceback_pairs.hh"
#include "scfg2/replay/factorized_max.hh"

#include <cstdint>
#include <numeric>
#include <optional>
#include <unordered_map>

namespace cparty::scfg2::decoding {

enum class PairPosteriorObjective { Mea, GammaCentroid };

struct DecodedPairPosterior {
  double expected_accuracy = 0.0;
  double expected_true_positive_pairs = 0.0;
  double expected_true_negative_pairs = 0.0;
  std::vector<engine::TracebackPair> pairs;  // Generated, ordered by opening.
  std::string generated_structure;
  std::size_t terms_evaluated = 0;
  std::size_t factors_evaluated = 0;
};

inline std::uint64_t posterior_pair_key(cand_pos_t i, cand_pos_t j) {
  return (std::uint64_t(static_cast<std::uint32_t>(i)) << 32) |
         static_cast<std::uint32_t>(j);
}

// Both objectives use one pair-owner mapping and one additive max surface.
// MEA:       2*gamma*TP_pair + correctly-unpaired bases.
// Centroid:  gamma*TP_pair + TN_pair, over all unordered positions i<j.
// Thermal scores decide support, never the posterior objective. Negative
// connector rewards remain available.
inline double pair_posterior_reward(PairPosteriorObjective objective,
                                    double p, double qi, double qj, double gamma) {
  return objective == PairPosteriorObjective::GammaCentroid
      ? (gamma + 1.0) * p - 1.0 : 2.0 * gamma * p - qi - qj;
}

struct PosteriorPairFeature {
  engine::TracebackPair pair;
  double probability;
};

// Consumer rewards live here, outside the local thermodynamic weight layer.
class PairPosteriorScoring {
 public:
  PairPosteriorScoring(const engine::OutsideMarginalSummary &marginal,
                      cand_pos_t n, double gamma, PairPosteriorObjective objective)
      : marginal_(marginal), gamma_(gamma), objective_(objective) {
    if (!marginal.ok || marginal.unpaired_probability.size() != std::size_t(n) + 1 ||
        !std::isfinite(gamma) || gamma < 0.0)
      throw std::invalid_argument("invalid decoding posterior or gamma");
    probability_.reserve(marginal.pairs.size());
    for (const auto &pair : marginal.pairs)
      probability_.emplace(posterior_pair_key(pair.i, pair.j), pair.union_probability);
  }

  template <class Context>
  std::optional<PosteriorPairFeature> feature(const deductions::ItemKey &item,
                                             const Context &ctx,
                                             bool pending_candidate = false) const {
    const bool round = item.nonterminal == grammar::NonTerminal::V;
    const bool square = item.nonterminal == grammar::NonTerminal::VP;
    if ((!round && !square) || ctx.pair_partner(item.i) == item.j) return {};
    const auto found = probability_.find(posterior_pair_key(item.i, item.j));
    // Replay also evaluates candidates that alpha subsequently prunes. They
    // have no marginal entry and cannot reach a retained consumer. Selected
    // features and frozen-chart owners must have an observed marginal.
    if (found == probability_.end() && !pending_candidate)
      throw std::runtime_error("selected generated pair has no posterior entry");
    return PosteriorPairFeature{
        {item.i, item.j, round ? engine::TracebackPairLane::Round
                              : engine::TracebackPairLane::Square},
        found == probability_.end() ? 0.0 : found->second};
  }

  double reward(const PosteriorPairFeature &feature) const {
    const auto &q = marginal_.unpaired_probability;
    return pair_posterior_reward(objective_, feature.probability,
                                q[feature.pair.i], q[feature.pair.j], gamma_);
  }
 private:
  const engine::OutsideMarginalSummary &marginal_;
  double gamma_;
  PairPosteriorObjective objective_;
  std::unordered_map<std::uint64_t, double> probability_;
};

template <class Context>
DecodedPairPosterior finish_pair_posterior(
    const std::vector<PosteriorPairFeature> &selected,
    const engine::OutsideMarginalSummary &marginal, cand_pos_t n, const Context &ctx,
    double gamma, PairPosteriorObjective objective, double maximum_score,
    std::size_t terms_evaluated, std::size_t factors_evaluated) {
  const auto &q = marginal.unpaired_probability;
  std::vector<cand_pos_t> partner(std::size_t(n) + 1);
  std::vector<engine::TracebackPairLane> lane(std::size_t(n) + 1);
  std::vector<double> selected_probability(std::size_t(n) + 1);
  for (cand_pos_t i = 1; i <= n; ++i) {
    const auto mate = ctx.pair_partner(i);
    if (mate > 0) partner[i] = mate;
  }
  for (const auto &feature : selected) {
    const auto &pair = feature.pair;
    if (partner[pair.i] || partner[pair.j])
      throw std::runtime_error("posterior decoding derivation reuses a paired endpoint");
    partner[pair.i] = pair.j;
    partner[pair.j] = pair.i;
    lane[pair.i] = pair.lane;
    selected_probability[pair.i] = feature.probability;
  }

  DecodedPairPosterior out;
  double correct_unpaired = 0.0;
  double selected_pair_reward = 0.0;
  double predicted_pair_count = 0.0;
  out.terms_evaluated = terms_evaluated;
  out.factors_evaluated = factors_evaluated;
  out.generated_structure.assign(std::size_t(n), '.');
  // Opening-index traversal renders and orders output in O(n), avoiding
  // repeated pair-set scans or an O(P log P) traceback sort.
  for (cand_pos_t i = 1; i <= n; ++i) {
    const auto j = partner[i];
    if (!j) correct_unpaired += q[i];
    if (j <= i) continue;
    ++predicted_pair_count;
    if (ctx.pair_partner(i) == j) {
      out.expected_true_positive_pairs += 1.0;
    } else {
      const double p = selected_probability[i];
      out.expected_true_positive_pairs += p;
      selected_pair_reward += pair_posterior_reward(objective, p, q[i], q[j], gamma);
      out.pairs.push_back({i, j, lane[i]});
      const bool round = lane[i] == engine::TracebackPairLane::Round;
      out.generated_structure[i - 1] = round ? '(' : '[';
      out.generated_structure[j - 1] = round ? ')' : ']';
    }
  }
  // E[TN] = total position pairs - E[true pairs] - E[false positives].
  // Sparse marginal mass gives E[true pairs], so no n-by-n TN scan is needed.
  double true_pair_mass = 0.0;
  for (const auto &pair : marginal.pairs) true_pair_mass += pair.union_probability;
  out.expected_true_negative_pairs = double(n) * (double(n) - 1.0) / 2.0
      - true_pair_mass - (predicted_pair_count - out.expected_true_positive_pairs);
  out.expected_accuracy = objective == PairPosteriorObjective::GammaCentroid
      ? gamma * out.expected_true_positive_pairs + out.expected_true_negative_pairs
      : 2.0 * gamma * out.expected_true_positive_pairs + correct_unpaired;
  // Validate the reduced pair objective before the potentially large TN
  // constant is added. That constant must not hide a traceback discrepancy.
  if (!std::isfinite(out.expected_accuracy) ||
      std::fabs(selected_pair_reward - maximum_score) >
          1e-8 * std::max(1.0, std::fabs(maximum_score)))
    throw std::runtime_error("decoded pair objective differs from forest maximum");
  return out;
}

template <class Chart, class Context>
DecodedPairPosterior decode_pair_posterior(const semirings::FactorizedInsideGraph &graph,
                     std::uint32_t root, const Chart &chart,
                     const engine::OutsideMarginalSummary &marginal,
                     cand_pos_t n, const Context &ctx, double gamma,
                     PairPosteriorObjective objective) {
  const PairPosteriorScoring scoring(marginal, n, gamma, objective);
  struct Feature { std::uint32_t node; PosteriorPairFeature pair; };
  std::vector<Feature> features;
  std::vector<double> reward(graph.states().size());
  for (const auto &[item, value] : chart.values()) {
    if (value.node >= reward.size()) continue;
    if (const auto feature = scoring.feature(item, ctx)) {
      reward[value.node] += scoring.reward(*feature);
      features.push_back({value.node, *feature});
    }
  }
  const auto maximum = replay::maximize_factorized_forest(graph, reward);
  const auto counts = replay::maximum_state_occurrences(graph, maximum, root);
  std::vector<PosteriorPairFeature> selected;
  for (const auto &feature : features) {
    if (!counts[feature.node]) continue;
    if (counts[feature.node] > 1)
      throw std::runtime_error("posterior decoding derivation repeats a generated-pair state");
    selected.push_back(feature.pair);
  }
  return finish_pair_posterior(selected, marginal, n, ctx, gamma, objective,
      maximum.scores.at(root), maximum.terms_evaluated, maximum.factors_evaluated);
}

}  // namespace cparty::scfg2::decoding

#endif

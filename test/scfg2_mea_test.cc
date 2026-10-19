// Independent exhaustive grammar oracle for MEA and Hamada gamma-centroid. Enumeration
// uses rule-local thermal weights only to form the posterior; it scores each
// output pair set directly, without the production forest reward recurrence.
#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/mea.hh"
#include "scfg2/api/centroid.hh"
#include "scfg2/api/detail.hh"
#include "scfg2/replay/factorized_max.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/beam/prepare_ranking.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <unordered_map>

namespace {
using namespace cparty::scfg2;
using Pair = std::pair<int, int>;
using Pairs = std::vector<Pair>;
using Groups = std::map<Pairs, double>;

bool require(bool condition, const std::string &message) {
  if (!condition) std::cerr << message << '\n';
  return condition;
}

bool core() {
  using S = semirings::FactorizedLogSumExpSemiring;
  semirings::FactorizedInsideGraph graph;
  S::Scope scope(graph);
  const auto x = S::from_weight(100.0).checkpointed();
  const auto y = S::from_weight(1.0).checkpointed();
  const auto z = S::from_weight(1.0).checkpointed();
  const auto root = S::plus(x, S::times(y, z)).checkpointed();
  std::vector<double> rewards(graph.states().size());
  rewards[x.node] = 0.5;
  rewards[y.node] = -1.0;  // Negative connectors must remain available.
  rewards[z.node] = 3.0;
  const auto maximum = replay::maximize_factorized_forest(graph, rewards);
  const auto counts = replay::maximum_state_occurrences(graph, maximum, root.node);
  return require(maximum.scores[root.node] == 2.0 && !counts[x.node] &&
                     counts[y.node] == 1 && counts[z.node] == 1 &&
                     maximum.terms_evaluated == graph.terms().size() &&
                     maximum.factors_evaluated == graph.factors().size(),
                 "additive forest maximum follows thermal weights or drops a negative connector");
}

Groups enumerate(replay::WFinalExactInsideContext &ctx, int n) {
  std::unordered_map<deductions::ItemKey, Groups, deductions::ItemKeyHash,
                     deductions::ItemKeyEq> chart;
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weights(&ctx);
  for (const auto &item : engine::build_pkfree_basic_exact_schedule(n, ctx)) {
    Groups groups;
    for (const auto &d : deductions::generate_exact_basic_deductions_for_item(item, ctx)) {
      if (!engine::legacy_span_schedule_supports(d)) continue;
      const double local = weights.weight_for(d);
      if (!(local > 0.0)) continue;
      Pairs own;
      for (const auto &p : engine::deduction_output_pairs(d, ctx)) {
        if (!engine::is_fixed_scaffold_pair(p, ctx) &&
            (d.parent.nonterminal == grammar::NonTerminal::V ||
             d.parent.nonterminal == grammar::NonTerminal::VP))
          own.emplace_back(p.i, p.j);
      }
      Groups partial{{own, local}};
      for (std::size_t c = 0; c < d.child_count; ++c) {
        auto key = d.children[c].item;
        Groups next;
        for (const auto &[left, lw] : partial) {
          for (const auto &[right, rw] : chart[key]) {
            Pairs pairs = left;
            pairs.insert(pairs.end(), right.begin(), right.end());
            std::sort(pairs.begin(), pairs.end());
            pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
            next[pairs] += lw * rw;
          }
        }
        partial = std::move(next);
      }
      for (const auto &[pairs, weight] : partial) groups[pairs] += weight;
    }
    chart[item] = std::move(groups);
  }
  return chart[{grammar::NonTerminal::W, 1, n, -1, -1}];
}

// Exhaustively form pair sets for the *actual* finite retained forest. This
// enumerates alternatives/products, rather than using additive max rewards.
std::set<Pairs> finite_outputs(std::string sequence, std::string scaffold,
                               const api::BeamConfig &config) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalLogWeightContext ctx(owner, config.wmb_closure_policy);
  const int n = sequence.size();
  const auto prepared = beam::prepare_ranking(ctx, n,
      replay::WFinalVpWmbpOnlineWidths::uniform(config.beam_size), config.beam_size,
      config.ranking, config.outer_pilot_exponent);
  using S = semirings::FactorizedLogSumExpSemiring;
  semirings::FactorizedInsideGraph graph;
  const auto chart = [&] {
    S::Scope scope(graph);
    return replay::run_w_final_vp_wmbp_online<S>(
        ctx, n, prepared.widths, prepared.ranking).chart;
  }();
  std::unordered_map<std::uint32_t, Pairs> owners;
  for (const auto &[item, value] : chart.values()) {
    const auto nt = item.nonterminal;
    if ((nt == grammar::NonTerminal::V || nt == grammar::NonTerminal::VP) &&
        ctx.pair_partner(item.i) != item.j)
      owners[value.node].emplace_back(item.i, item.j);
  }
  std::vector<std::set<Pairs>> outputs(graph.states().size());
  for (std::size_t node = 0; node < outputs.size(); ++node) {
    const auto &state = graph.states()[node];
    for (std::size_t t = state.first_term;
         t < std::size_t(state.first_term) + state.term_count; ++t) {
      std::set<Pairs> partial{owners[node]};
      const auto &term = graph.terms()[t];
      for (std::size_t f = term.first_factor;
           f < std::size_t(term.first_factor) + term.factor_count; ++f) {
        std::set<Pairs> next;
        for (const auto &left : partial) {
          for (const auto &right : outputs[graph.factors()[f]]) {
            Pairs joined = left;
            joined.insert(joined.end(), right.begin(), right.end());
            std::sort(joined.begin(), joined.end());
            joined.erase(std::unique(joined.begin(), joined.end()), joined.end());
            next.insert(std::move(joined));
          }
        }
        partial = std::move(next);
      }
      outputs[node].insert(partial.begin(), partial.end());
    }
  }
  return outputs[chart.get({grammar::NonTerminal::W, 1, n, -1, -1}).node];
}

bool check(std::string sequence, std::string scaffold, bool &saw_pk,
           bool &saw_centroid_pk, bool &saw_difference,
           api::BeamRanking ranking) {
  const int n = sequence.size();
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto groups = enumerate(ctx, n);
  double z = 0.0;
  for (const auto &[pairs, weight] : groups) z += weight;
  if (!require(z > 0.0, "enumeration has no root")) return false;
  std::map<Pair, double> p;
  std::vector<double> q(n + 1, 1.0);
  int fixed = 0;
  for (int i = 1; i <= n; ++i) {
    if (ctx.pair_partner(i) > 0) q[i] = 0.0;
    if (ctx.pair_partner(i) > i) ++fixed;
  }
  for (const auto &[pairs, weight] : groups)
    for (const auto &pair : pairs) p[pair] += weight / z;
  for (const auto &[pair, probability] : p) {
    q[pair.first] -= probability;
    q[pair.second] -= probability;
  }
  const auto score = [&](const Pairs &pairs, double gamma) {
    std::vector<bool> paired(n + 1);
    double value = 2.0 * gamma * fixed;
    for (const auto &pair : pairs) {
      value += 2.0 * gamma * p.at(pair);
      paired[pair.first] = paired[pair.second] = true;
    }
    for (int i = 1; i <= n; ++i)
      if (!paired[i]) value += q[i];
    return value;
  };
  api::BeamConfig finite_config{3, 30, 30};
  finite_config.ranking = ranking;
  const auto retained = finite_outputs(sequence, scaffold, finite_config);
  api::BeamSession session(sequence, scaffold, finite_config);
  const auto outside = session.outside();
  if (!require(outside.ok && outside.reachable, "finite posterior failed: " + outside.error)) return false;
  const auto finite_score = [&](const Pairs &pairs, double gamma) {
    double value = 0.0;
    std::vector<bool> paired(n + 1);
    for (const auto &prob : outside.pair_probabilities) {
      if (ctx.pair_partner(prob.i) == prob.j ||
          std::binary_search(pairs.begin(), pairs.end(), Pair{prob.i, prob.j})) {
        value += 2.0 * gamma * prob.union_probability;
        paired[prob.i] = paired[prob.j] = true;
      }
    }
    for (int i = 1; i <= n; ++i) if (!paired[i]) value += outside.unpaired_probabilities[i];
    return value;
  };
  auto exact_posterior = p;
  std::map<Pair, double> finite_posterior;
  for (int i = 1; i <= n; ++i)
    if (ctx.pair_partner(i) > i) exact_posterior[{i, ctx.pair_partner(i)}] = 1.0;
  for (const auto &prob : outside.pair_probabilities)
    finite_posterior[{prob.i, prob.j}] = prob.union_probability;
  // Directly accumulate TP/TN over *every* position pair. No reduced reward
  // formula or sparse TN reconstruction from production is used here.
  const auto pair_counts = [&](const Pairs &pairs, const std::map<Pair, double> &posterior) {
    double tp = 0.0, tn = 0.0;
    for (int i = 1; i < n; ++i) for (int j = i + 1; j <= n; ++j) {
      const auto it = posterior.find({i, j});
      const double probability = it == posterior.end() ? 0.0 : it->second;
      const bool predicted = ctx.pair_partner(i) == j ||
          std::binary_search(pairs.begin(), pairs.end(), Pair{i, j});
      if (predicted) tp += probability;
      else tn += 1.0 - probability;
    }
    return std::pair<double, double>{tp, tn};
  };
  const auto centroid_score = [&](const Pairs &pairs, double gamma,
                                   const std::map<Pair, double> &posterior) {
    const auto [tp, tn] = pair_counts(pairs, posterior);
    return gamma * tp + tn;
  };
  for (double gamma : {0.0, 0.5, 1.0, 3.0, 10.0, 1000.0}) {
    double best = -1.0;
    for (const auto &[pairs, weight] : groups) best = std::max(best, score(pairs, gamma));
    const auto exact = api::exact_mea(sequence, scaffold, gamma);
    const auto saturated = api::beam_mea(sequence, scaffold, gamma,
                                        api::BeamConfig{std::size_t(n + 1), std::size_t(n), std::size_t(n)});
    for (const auto *result : {&exact, &saturated}) {
      if (!require(result->ok && result->reachable, "MEA failed: " + result->error)) return false;
      Pairs pairs;
      for (const auto &pair : result->generated_pairs) pairs.emplace_back(pair.i, pair.j);
      if (!require(groups.count(pairs) && std::fabs(result->expected_accuracy - best) < 1e-8 &&
                       std::fabs(score(pairs, gamma) - best) < 1e-8 &&
                       std::fabs(result->log_partition_weight - std::log(z)) < 1e-9,
                   "MEA differs from exhaustive objective: " + sequence + " / " + scaffold)) return false;
      saw_pk = saw_pk || result->generated_structure.find('[') != std::string::npos;
    }
    const auto finite = api::beam_mea(sequence, scaffold, gamma, finite_config);
    if (!require(finite.ok && finite.reachable, "finite MEA failed: " + finite.error)) return false;
    if (!require(std::fabs(outside.log_partition_weight - finite.log_partition_weight) < 1e-9,
                 "MEA changed the inside-ranked beam")) return false;
    if (!require(api::ExactSession(sequence, scaffold).probability(finite.generated_pairs).reachable,
                 "finite MEA is outside the exact grammar")) return false;
    Pairs chosen;
    for (const auto &pair : finite.generated_pairs) chosen.emplace_back(pair.i, pair.j);
    double finite_best = -1.0;
    for (const auto &pairs : retained) finite_best = std::max(finite_best, finite_score(pairs, gamma));
    if (!require(retained.count(chosen) &&
                     std::fabs(finite_score(chosen, gamma) - finite.expected_accuracy) < 1e-8 &&
                     std::fabs(finite_best - finite.expected_accuracy) < 1e-8,
                 "finite MEA differs from exhaustive retained-forest optimum")) return false;
    double centroid_best = -1.0;
    for (const auto &[pairs, weight] : groups)
      centroid_best = std::max(centroid_best, centroid_score(pairs, gamma, exact_posterior));
    const auto centroid_exact = api::exact_centroid(sequence, scaffold, gamma);
    const auto centroid_saturated = api::beam_centroid(sequence, scaffold, gamma,
        api::BeamConfig{std::size_t(n + 1), std::size_t(n), std::size_t(n)});
    for (const auto *result : {&centroid_exact, &centroid_saturated}) {
      if (!require(result->ok && result->reachable, "centroid failed: " + result->error)) return false;
      Pairs pairs;
      for (const auto &pair : result->generated_pairs) pairs.emplace_back(pair.i, pair.j);
      const auto [tp, tn] = pair_counts(pairs, exact_posterior);
      if (!require(groups.count(pairs) &&
                       std::fabs(result->expected_accuracy - centroid_best) < 1e-8 &&
                       std::fabs(gamma * tp + tn - centroid_best) < 1e-8 &&
                       std::fabs(result->expected_true_positive_pairs - tp) < 1e-9 &&
                       std::fabs(result->expected_true_negative_pairs - tn) < 1e-9 &&
                       std::fabs(result->log_partition_weight - std::log(z)) < 1e-9,
                   "centroid differs from direct exhaustive gamma*TP+TN")) return false;
      saw_centroid_pk = saw_centroid_pk || result->generated_structure.find('[') != std::string::npos;
      if (best - score(pairs, gamma) > 1e-8) {
        if (!saw_difference)
          std::cout << "different_objectives_witness sequence=" << sequence << " scaffold=" << scaffold
                    << " gamma=" << gamma << " mea=" << exact.union_structure
                    << " centroid=" << result->union_structure << '\n';
        saw_difference = true;
      }
    }
    const auto centroid_finite = api::beam_centroid(sequence, scaffold, gamma, finite_config);
    if (!require(centroid_finite.ok && centroid_finite.reachable,
                 "finite centroid failed: " + centroid_finite.error)) return false;
    Pairs centroid_chosen;
    for (const auto &pair : centroid_finite.generated_pairs) centroid_chosen.emplace_back(pair.i, pair.j);
    double finite_centroid_best = -1.0;
    for (const auto &pairs : retained)
      finite_centroid_best = std::max(finite_centroid_best, centroid_score(pairs, gamma, finite_posterior));
    const auto [finite_tp, finite_tn] = pair_counts(centroid_chosen, finite_posterior);
    if (!require(retained.count(centroid_chosen) &&
                     std::fabs(centroid_finite.expected_accuracy - finite_centroid_best) < 1e-8 &&
                     std::fabs(gamma * finite_tp + finite_tn - finite_centroid_best) < 1e-8 &&
                     std::fabs(centroid_finite.log_partition_weight - finite.log_partition_weight) < 1e-9,
                 "finite centroid differs from exhaustive retained gamma*TP+TN or changed support")) return false;
  }
  return true;
}
}  // namespace

int main(int argc, char **argv) {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  if (!core()) return EXIT_FAILURE;
  bool saw_pk = false, saw_centroid_pk = false, saw_difference = false;
  for (const auto &fixture : std::vector<std::pair<std::string, std::string>>{
      {"AAAAAA", "......"}, {"GCGCGCGCGCGCGC", "(.(......)...)"},
      {"GGGAAACCCAGGGAAACCC", "(((...))).........."},
      {"GGGAAAUCCCAGG", "............."}})
    for (auto ranking : {api::BeamRanking::Baseline, api::BeamRanking::Outer})
      if (!check(fixture.first, fixture.second, saw_pk, saw_centroid_pk, saw_difference, ranking)) return EXIT_FAILURE;
  // Every free G/C partner crosses this fixed stem; other bases are A.
  // High gamma selects a genuine PK even when its thermodynamic mass is small.
  std::string crossing_sequence(24, 'A'), crossing_scaffold(24, '.');
  for (int k = 0; k < 3; ++k) {
    crossing_sequence[3 + k] = crossing_sequence[9 + k] = 'G';
    crossing_sequence[14 - k] = crossing_sequence[20 - k] = 'C';
    crossing_scaffold[3 + k] = '(';
    crossing_scaffold[14 - k] = ')';
  }
  for (auto ranking : {api::BeamRanking::Baseline, api::BeamRanking::Outer})
    if (!check(crossing_sequence, crossing_scaffold, saw_pk, saw_centroid_pk, saw_difference, ranking)) return EXIT_FAILURE;
  std::mt19937 rng(20261008);
  const int cases = argc > 1 ? std::atoi(argv[1]) : 8;
  for (int c = 0; c < cases; ++c) {
    std::string seq(12 + c % 9, 'A');
    for (auto &base : seq) base = "ACGU"[rng() % 4];
    std::string g(seq.size(), '.');
    if (c % 2) { g[2] = '('; g[seq.size() - 3] = ')'; seq[seq.size() - 3] = 'C'; seq[2] = 'G'; }
    for (auto ranking : {api::BeamRanking::Baseline, api::BeamRanking::Outer})
      if (!check(seq, g, saw_pk, saw_centroid_pk, saw_difference, ranking)) return EXIT_FAILURE;
  }
  if (!require(saw_pk && saw_centroid_pk && saw_difference,
               "posterior decoding witnesses missed PK or distinct objective choices")) return EXIT_FAILURE;
  for (double gamma : {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    if (!require(!api::beam_mea("AAAAAA", "......", gamma).ok &&
                     !api::beam_centroid("AAAAAA", "......", gamma).ok &&
                     !api::exact_centroid("AAAAAA", "......", gamma).ok, "invalid gamma accepted")) return EXIT_FAILURE;
  if (!require(!api::beam_mea("", "").ok && !api::beam_mea("AAAAAA", "...").ok &&
                   !api::beam_mea("AAAAAA", "......", 1.0, api::BeamConfig{0}).ok &&
                   !api::beam_centroid("", "").ok && !api::beam_centroid("AAAAAA", "...").ok &&
                   !api::beam_centroid("AAAAAA", "......", 1.0, api::BeamConfig{0}).ok,
               "invalid MEA inputs accepted")) return EXIT_FAILURE;
  std::cout << "scfg2_mea_and_centroid=ok cases=" << cases + 5 << " rankings=2 gamma_values=6 crossing_witness=1 different_objectives_witness=1\n";
}

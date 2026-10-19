#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/detail.hh"
#include "scfg2/api/beam_detail.hh"
#include "scfg2/beam/retained_outside.hh"
#include "scfg2/decoding/inside_pair_max.hh"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace cparty::scfg2;
using S = semirings::InsidePairMaxSemiring;
using F = semirings::FactorizedLogSumExpSemiring;

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

void algebra() {
  auto storage = std::make_shared<semirings::PairMaxStorage>();
  {
    S::Scope scope(storage);
    deductions::ItemKey a{grammar::NonTerminal::V, 1, 10, -1, -1};
    auto x = S::with_pair_reward(S::from_weight(100.0), a, -1.0).checkpointed();
    auto y = S::with_pair_reward(S::one(), a, 3.0).checkpointed();
    const auto z = S::plus(x, y).checkpointed();
    require(z.score == 3.0 && z.log_value == std::log(101.0), "posterior score follows alpha");
    require(S::better(x, y) && S::priority(x) == std::log(100.0), "scheduler uses posterior score");
    const auto tie = S::plus(y, S::with_pair_reward(S::one(), a, 3.0)).checkpointed();
    require(tie.node == y.node, "first tie alternative was lost");
    struct Context { int pair_partner(int) const { return 0; } } ctx;
    engine::OutsideMarginalSummary posterior;
    posterior.ok = true;
    posterior.unpaired_probability.assign(11, 0.0);
    posterior.pairs.push_back({1, 10, 1.0, 0.0, 1.0});
    bool repeated = false;
    try {
      decoding::decode_inside_pair_max(S::times(y, y), posterior, 10, ctx, 1.0,
                                       decoding::PairPosteriorObjective::GammaCentroid);
    } catch (const std::runtime_error &) { repeated = true; }
    require(repeated, "shared generated feature occurrence was silently skipped");
    for (int i = 0; i < 100000; ++i)
      x = S::with_pair_reward(S::times(x, S::one()), a, 0.0).checkpointed();
    require(storage->nodes >= 100001, "deep chain witness is empty");
  }
  require(storage->nodes == 0 && storage->factors == 0 && storage->products == 0 &&
          storage->product_factors == 0, "chosen trace storage leaked");
}

std::size_t compare(std::string sequence, std::string scaffold, api::BeamConfig config) {
  weights::RuntimeSnapshot owner(sequence, scaffold);
  replay::WFinalLogWeightContext ctx(owner, config.wmb_closure_policy);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
  const auto prepared = api::detail::prepare_beam_ranking(ctx, n, config);
  semirings::FactorizedInsideGraph graph;
  const auto frozen = [&] {
    F::Scope scope(graph);
    return replay::run_w_final_vp_wmbp_online<F>(ctx, n, prepared.widths, prepared.ranking);
  }();
  const auto marginal = engine::collect_outside_marginals(
      frozen.chart, graph.posteriors(frozen.chart.get(root).node), n, ctx);
  require(marginal.ok, "frozen posterior failed");

  using L = semirings::LogSumExpSemiring;
  beam::RetainedRuleAdmission<L> rules;
  const auto numeric = replay::run_w_final_vp_wmbp_online_observed<L>(
      ctx, n, prepared.widths, rules, prepared.ranking);
  rules.finish(numeric.chart);
  beam::RetainedOutsideProgram<L, replay::WFinalLogWeightContext> outside(
      n, numeric.chart, ctx, rules, prepared.widths, prepared.ranking);
  const auto reverse = outside.outside(root);
  const auto regenerated = engine::collect_outside_marginals_by_item(numeric.chart,
      [&](const auto &item, double alpha) {
        return std::exp(alpha + reverse.beta.get(item) - numeric.chart.get(root));
      }, n, ctx);
  require(regenerated.ok, "regenerated posterior failed");
  require(numeric.chart.values().size() == frozen.chart.values().size(), "numeric support differs");
  std::size_t comparisons = 0;
  std::vector<double> gammas{0.0, 0.5, 1.0, 3.0, 10.0, 1000.0};
  for (const auto &pair : marginal.pairs) {
    if (pair.union_probability > 0.0 && pair.union_probability < 1.0) {
      gammas.push_back(1.0 / pair.union_probability - 1.0);
      break; // A zero-reward centroid threshold, in addition to gamma=0.
    }
  }
  auto storage = std::make_shared<semirings::PairMaxStorage>();
  for (const auto objective : {decoding::PairPosteriorObjective::Mea,
                               decoding::PairPosteriorObjective::GammaCentroid}) {
    for (double gamma : gammas) {
      for (const auto *posterior : {&marginal, &regenerated}) {
        const decoding::PairPosteriorScoring scoring(*posterior, n, gamma, objective);
        const auto compact = [&] {
          S::Scope scope(storage);
          return replay::run_w_final_vp_wmbp_online_observed_transformed<S>(
              ctx, n, prepared.widths,
              [](const auto &, const auto &, auto, const auto &, auto) {},
              [&](const deductions::Deduction &d, auto value) {
                const auto feature = scoring.feature(d.parent, ctx, true);
                return feature ? S::with_pair_reward(value, d.parent, scoring.reward(*feature)) : value;
              }, prepared.ranking);
        }();
        require(compact.chart.values().size() == numeric.chart.values().size(), "compact support differs");
        for (const auto &[item, alpha] : numeric.chart.values())
          require(compact.chart.get(item).log_value == alpha, "compact alpha differs");
        const auto expected = decoding::decode_pair_posterior(graph, frozen.chart.get(root).node,
            frozen.chart, *posterior, n, ctx, gamma, objective);
        const auto actual = decoding::decode_inside_pair_max(compact.chart.get(root),
            *posterior, n, ctx, gamma, objective);
        if (actual.generated_structure != expected.generated_structure ||
            std::fabs(actual.expected_accuracy - expected.expected_accuracy) > 1e-8) {
          std::cerr << sequence << '/' << scaffold << " b=" << config.beam_size
                    << " gamma=" << gamma << " objective=" << int(objective)
                    << "\nold=" << expected.generated_structure << " " << expected.expected_accuracy
                    << "\nnew=" << actual.generated_structure << " " << actual.expected_accuracy << '\n';
          throw std::runtime_error("compact posterior decoding differs from frozen forest");
        }
        require(storage->max_arity <= 8, "production selected factor arity is unbounded");
        const auto budget = std::size_t(n) * (config.beam_size + 2);
        require(storage->peak_nodes + storage->peak_factors + storage->peak_products +
                storage->peak_product_factors <= 32 * budget, "selected storage exceeds linear witness budget");
        ++comparisons;
      }
      require(storage->nodes == 0 && storage->factors == 0 && storage->products == 0 &&
              storage->product_factors == 0, "decoder storage survived its chart");
    }
  }
  std::cout << "n=" << n << " b=" << config.beam_size << " peak_nodes=" << storage->peak_nodes
            << " peak_factors=" << storage->peak_factors << " peak_products=" << storage->peak_products
            << " arity=" << storage->max_arity << '\n';
  return comparisons;
}
} // namespace

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", SCFG2_TEST_PARAM_FILE, 1);
    algebra();
    std::size_t comparisons = 0;
    for (const auto &fixture : std::vector<std::pair<std::string, std::string>>{
        {"GCGCAAAAGCGCGC", ".............."},
        {"GGGGAACCCCUUUU", "((....))......"},
        {"GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
         "....(...)((....)..)((.(.((...))))(.(....)((...))))"}}) {
      require(fixture.first.size() == fixture.second.size(), "bad fixture");
      for (std::size_t b : {std::size_t(3), std::size_t(20)}) {
        if (b > fixture.first.size()) continue;
        for (auto ranking : {api::BeamRanking::Baseline, api::BeamRanking::Outer}) {
          for (auto policy : {WmbClosurePolicy::LegacyCompatible, WmbClosurePolicy::RequireWeakClosure}) {
            api::BeamConfig config{b, 30, 30};
            config.ranking = ranking;
            config.wmb_closure_policy = policy;
            comparisons += compare(fixture.first, fixture.second, config);
          }
        }
      }
    }
    std::mt19937 rng(20261014);
    for (int sample = 0; sample < 12; ++sample) {
      std::string scaffold = sample % 2 ? "(.(......)...).................."
                                        : "((....))........((....))........";
      std::string sequence(scaffold.size(), 'A');
      std::vector<int> openings;
      for (std::size_t i = 0; i < scaffold.size(); ++i) {
        sequence[i] = "ACGU"[rng() % 4];
        if (scaffold[i] == '(') openings.push_back(i);
        if (scaffold[i] == ')') {
          sequence[openings.back()] = 'G'; sequence[i] = 'C'; openings.pop_back();
        }
      }
      for (auto ranking : {api::BeamRanking::Baseline, api::BeamRanking::Outer}) {
        api::BeamConfig config{sample % 2 ? std::size_t(3) : std::size_t(8), 30, 30};
        config.ranking = ranking;
        config.outer_pilot_exponent = sample % 2 ? api::OuterPilotExponent::Half
                                                : api::OuterPilotExponent::TwoThirds;
        config.wmb_closure_policy = sample % 3 ? WmbClosurePolicy::RequireWeakClosure
                                              : WmbClosurePolicy::LegacyCompatible;
        comparisons += compare(sequence, scaffold, config);
      }
    }
    std::cout << "scfg2_inside_pair_max=ok comparisons=" << comparisons << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

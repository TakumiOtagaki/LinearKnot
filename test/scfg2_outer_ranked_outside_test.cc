#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/semirings/logsumexp_tape.hh"
#include "scfg2/semirings/value_boundary.hh"

// Distinguish public assignments from aliases of private carries in the
// independent tape reference, as in the baseline retained-outside witness.
namespace cparty::scfg2::semirings {
template <> struct ValueBoundary<LogSumExpTapeValue> {
  static LogSumExpTapeValue checkpoint(LogSumExpTapeValue value) {
    return TapeLogSumExpSemiring::is_zero(value) ? value :
        TapeLogSumExpSemiring::times(value, TapeLogSumExpSemiring::one());
  }
  static bool discard(const LogSumExpTapeValue &) { return false; }
};
}

#include "scfg2/beam/prepare_ranking.hh"
#include "scfg2/beam/retained_outside.hh"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>

namespace {
using namespace cparty::scfg2;
using S = semirings::LogSumExpSemiring;
using T = semirings::TapeLogSumExpSemiring;
using Key = deductions::ItemKey;
using NT = grammar::NonTerminal;
std::size_t cases = 0, parents = 0, positive_posteriors = 0;

void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
bool same_log(double a, double b) {
  return a == b || std::fabs(a - b) <= 2e-9;
}

void preparation_boundaries() {
  std::string seq = "AAAAAA", g = "......";
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
  replay::WFinalLogWeightContext ctx(owner);
  const auto run = [&](std::size_t b, beam::RankingMode mode, beam::PilotExponent q) {
    return beam::prepare_ranking(ctx, 6, replay::WFinalVpWmbpOnlineWidths::uniform(b), b, mode, q);
  };
  require(!run(2, beam::RankingMode::Baseline, beam::PilotExponent::TwoThirds).ranking,
          "baseline unexpectedly builds a guide");
  const auto keep_all = run(std::numeric_limits<std::size_t>::max(),
                            beam::RankingMode::Outer, beam::PilotExponent::Half);
  require(!keep_all.ranking && keep_all.pilot_width == 0,
          "KeepAll needlessly builds an unbounded pilot");
  for (int which = 0; which < 4; ++which) {
    bool rejected = false;
    try {
      auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(which == 3
          ? std::numeric_limits<std::size_t>::max() : 2);
      if (which == 3) widths.wmb = 1;
      (void)beam::prepare_ranking(ctx, 6, widths, which == 0 ? 0 : widths.v,
          which == 1 ? static_cast<beam::RankingMode>(99) : beam::RankingMode::Outer,
          which == 2 ? static_cast<beam::PilotExponent>(99) : beam::PilotExponent::TwoThirds);
    } catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "invalid outer preparation was accepted");
  }
}

void view_and_fifo_witness() {
  engine::Chart<double> chart(S::zero());
  const Key high{NT::WMB, 1, 12, -1, -1}, guided{NT::WMB, 2, 12, -1, -1};
  const Key first{NT::BE, 1, 10, 2, 8}, second{NT::BE, 1, 10, 3, 7};
  chart.set(high, 10.0); chart.set(guided, 0.0);
  chart.set({NT::W, 1, 1, -1, -1}, 0.0);
  chart.set(first, 0.0); chart.set(second, 0.0);
  chart.set({NT::WI, 3, 5, -1, -1}, 0.0);
  chart.set({NT::WI, 4, 5, -1, -1}, 0.0);
  chart.set({NT::VP, 6, 12, -1, -1}, 0.0);
  const beam::PublicItemRanking rank = [&](const Key &key, double raw, double existing) {
    if (deductions::ItemKeyEq{}(key, high)) return raw - 100.0;
    if (deductions::ItemKeyEq{}(key, second)) return raw + 10.0;
    return existing;
  };
  const auto views = beam::RetainedChartViews::build<S>(12, chart, 1, rank);
  require(views.internal_starts(NT::WMB, 12) == std::vector<cand_pos_t>{1} &&
          views.root_starts(NT::WMB, 12) == std::vector<cand_pos_t>{2},
          "guided WMB root view collapsed into raw internal view");
  require(deductions::ItemKeyEq{}(views.items(NT::BE, 10).front(), second),
          "retained BE order ignored the forward ranking");
  const auto baseline = beam::WmbpWiPrefixAdmission::from_retained_chart<S>(12, chart, 1, 1);
  const auto outer = beam::WmbpWiPrefixAdmission::from_retained_chart<S>(12, chart, 1, 1, rank);
  const auto source = [&](const Key &be, cand_pos_t wi_start) {
    deductions::Deduction d;
    d.parent = {NT::WMBP, 1, 12, -1, -1};
    d.rule = grammar::RuleId::WMBP_SPLIT_BE_WI_VP;
    d.child_count = 3;
    d.children[0].item = be;
    d.children[1].item = {NT::WI, wi_start, 5, -1, -1};
    d.children[2].item = {NT::VP, 6, 12, -1, -1};
    return d;
  };
  require(baseline.allows<S>(source(first, 3), chart) &&
          !baseline.allows<S>(source(second, 4), chart) &&
          outer.allows<S>(source(second, 4), chart) &&
          !outer.allows<S>(source(first, 3), chart),
          "BE FIFO boundary did not preserve the guided commit order");
}

void check(std::string seq, std::string scaffold, std::size_t b,
           beam::PilotExponent exponent, WmbClosurePolicy policy) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, scaffold);
  replay::WFinalLogWeightContext ctx(owner, policy);
  const auto n = static_cast<cand_pos_t>(seq.size());
  auto prepared = beam::prepare_ranking(ctx, n,
      replay::WFinalVpWmbpOnlineWidths::uniform(b), b,
      beam::RankingMode::Outer, exponent);
  require(prepared.ranking && prepared.guide_items > 0 && prepared.widths.forced_items,
          "guided preparation omitted scores or shared forced keys");
  beam::RetainedRuleAdmission<S> rules;
  const auto run = replay::run_w_final_vp_wmbp_online_observed<S>(
      ctx, n, prepared.widths, rules, prepared.ranking);
  rules.finish(run.chart);
  beam::RetainedOutsideProgram<S, replay::WFinalLogWeightContext> program(
      n, run.chart, ctx, rules, prepared.widths, prepared.ranking);
  engine::AllLocalWeightModel<replay::WFinalLogWeightContext> weight(&ctx);
  for (cand_pos_t j = n; j >= 1; --j) {
    program.prepare_column(j);
    for (auto nt : beam::vp_wmbp_online_column_phases) {
      for (const auto &parent : program.views().items(nt, j)) {
        auto alpha = program.private_contribution(parent);
        for (const auto &d : program.deductions_for(parent)) {
          auto term = S::from_weight(weight.local_weight_for(d));
          for (std::size_t c = 0; c < d.child_count; ++c)
            term = S::times(term, run.chart.get(engine::RuntimeDeductionChildItem{}(d, c)));
          alpha = S::plus(alpha, term);
        }
        if (!same_log(alpha, run.chart.get(parent))) {
          std::cerr << "guided alpha mismatch case=" << cases << " family="
                    << grammar::nonterminal_name(nt) << " i=" << parent.i
                    << " j=" << parent.j << " rebuilt=" << alpha
                    << " forward=" << run.chart.get(parent) << '\n';
          throw std::runtime_error("guided retained source reconstruction failed");
        }
        ++parents;
      }
    }
  }
  const Key root{NT::W, 1, n, -1, -1};
  const auto reverse = program.outside(root);
  semirings::LogSumExpTape tape;
  const auto reference = [&] {
    T::Scope scope(tape);
    return replay::run_w_final_vp_wmbp_online<T>(
        ctx, n, prepared.widths, prepared.ranking);
  }();
  require(run.chart.values().size() == reference.chart.values().size(),
          "guided tape/numeric chart support differs");
  for (const auto &[key, value] : reference.chart.values())
    require(same_log(value.log_value, run.chart.get(key)), "guided tape/numeric alpha differs");
  require(std::isfinite(run.chart.get(root)), "guided witness lost its root");
  const auto expected = tape.reverse(reference.chart.get(root).node);
  for (const auto &[key, value] : reference.chart.values()) {
    const double a = expected.at(value.node);
    const double z = std::exp(reverse.beta.get(key) + run.chart.get(key) - run.chart.get(root));
    if (std::fabs(a - z) > 1e-7 * std::max(a, z)) {
      std::cerr << "guided beta mismatch case=" << cases << " family="
                << grammar::nonterminal_name(key.nonterminal) << " i=" << key.i
                << " j=" << key.j << " tape=" << a << " regenerated=" << z << '\n';
      throw std::runtime_error("guided reverse differs from the arithmetic tape");
    }
    positive_posteriors += a > 0.0;
  }
  const auto targeted = replay::run_w_final_vp_wmbp_online_target_probability(
      ctx, n, prepared.widths, std::vector<engine::TracebackPair>{}, prepared.ranking);
  require(same_log(targeted.log_partition, run.chart.get(root)),
          "guided targeted denominator uses another beam");
  const auto best = replay::run_w_final_vp_wmbp_online_viterbi_with_semiring<
      semirings::LogMaxSemiring>(ctx, n, prepared.widths, prepared.ranking);
  require(std::isfinite(best.dp.chart.get(root)) &&
          best.dp.chart.values().size() == best.choices.values().size(),
          "guided Viterbi lost its root or retained choice ownership");
  ++cases;
}
} // namespace

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    preparation_boundaries();
    view_and_fifo_witness();
    std::vector<std::pair<std::string, std::string>> inputs{
      {"GGAAACGAAACC", "((...)(...))"},
      {"GCGCGCGCGCGCGCGCGCGCGCGC", "........................"},
      {"GCGCGCGCGCGCGC", "(.(......)...)"},
      {"AGGUACCCUACAAUUACGGCAUCAGUCUAAAAACUCUGGGAUGGGCA",
       "(.(....).....)..(..(....)........(....)).(....)"}};
    std::ifstream fixture("test/fixtures/scfg2_wmbp_leftward_be_20261008/stem20.txt");
    std::string seq, scaffold;
    require(bool(std::getline(fixture, seq) && std::getline(fixture, scaffold)), "missing BE fixture");
    inputs.emplace_back(seq, scaffold);
    for (auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible})
      for (const auto &[s, g] : inputs)
        for (std::size_t b : {1, 3, 20})
          for (auto q : {beam::PilotExponent::Half, beam::PilotExponent::TwoThirds})
            check(s, g, b, q, policy);
    require(positive_posteriors > 0, "guided outside reference is vacuous");
    std::cout << "scfg2_outer_ranked_outside=ok cases=" << cases
              << " parents=" << parents << " positive=" << positive_posteriors << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

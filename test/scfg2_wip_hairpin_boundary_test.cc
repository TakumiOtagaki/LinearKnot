#include "scfg2/weights/runtime_snapshot.hh"
// Independent energy witnesses for the WIP split/start TURN boundary.
// DP09, 37 C, dangles=0: H3=3.69, H4=2.88, b'=0.56, c'=0.12 kcal/mol.
#include "scfg2/api/beam.hh"
#include "scfg2/api/exact.hh"
#include "scfg2/api/profiles.hh"
#include "scfg2/beam/retained_outside.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {
using namespace cparty::scfg2;
using S = semirings::LogSumExpSemiring;
using M = semirings::LogMaxSemiring;
using NT = grammar::NonTerminal;
using R = grammar::RuleId;
using Key = deductions::ItemKey;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

bool close(double a, double b) { return std::fabs(a - b) < 1e-10; }

void check_local(std::string seq, std::string g, cand_pos_t i, cand_pos_t j,
                 double energy, R rule, cand_pos_t k) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g, cparty::scfg2::weights::RuntimeOptions{0});
  replay::WFinalLogWeightContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(seq.size());
  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  const auto allow = [](const auto &) { return true; };
  const auto ignore = [](const auto &, const auto &, const auto &, const auto &, auto) {};
  const auto inside = replay::run_w_final_exact_dp_schedule<S>(schedule, ctx, allow, ignore);
  const auto best = replay::run_w_final_exact_dp_schedule<M>(schedule, ctx, allow, ignore);
  const Key parent{NT::WIP, i, j, -1, -1};
  const std::string label = "WIP(" + std::to_string(i) + "," + std::to_string(j) + ") " + seq;
  require(close(owner.partition_free_energy_from_log_weight(inside.get(parent)), energy),
          label + " exact inside energy != " + std::to_string(energy));
  require(close(owner.partition_free_energy_from_log_weight(best.get(parent)), energy),
          label + " exact Viterbi energy mismatch");

  // One legal derivation in each witness. Test finite sparse starts/joins and
  // KeepAll's dense prefix separately, including regeneration for outside.
  for (const auto widths : {replay::WFinalVpWmbpOnlineWidths::uniform(20),
                            replay::WFinalVpWmbpOnlineWidths{}}) {
    beam::RetainedRuleAdmission<S> admission;
    const auto run = replay::run_w_final_vp_wmbp_online_observed<S>(ctx, n, widths, admission);
    admission.finish(run.chart);
    require(close(run.chart.get(parent), inside.get(parent)), label + " beam inside mismatch");
    const auto max_run = replay::run_w_final_vp_wmbp_online<M>(ctx, n, widths);
    require(close(max_run.chart.get(parent), best.get(parent)), label + " beam Viterbi mismatch");

    beam::RetainedOutsideProgram<S, replay::WFinalLogWeightContext> program(
        n, run.chart, ctx, admission, widths);
    program.prepare_column(j);
    const auto sources = program.deductions_for(parent);
    require(sources.size() == 1 && sources[0].rule == rule &&
                (k < 0 || sources[0].split.k == k),
            label + " regenerated the wrong start/split");
    const auto reverse = program.outside(parent);
    const Key child{NT::V, k < 0 ? i : k, j, -1, -1};
    const double posterior = std::exp(reverse.beta.get(child) + run.chart.get(child) - run.chart.get(parent));
    require(close(posterior, 1.0), label + " lost its hairpin in outside regeneration");
  }
}

void check_too_short() {
  std::string seq = "AAGAAC", g = "......";
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g, cparty::scfg2::weights::RuntimeOptions{0});
  replay::WFinalLogWeightContext ctx(owner);
  const Key parent{NT::WIP, 2, 6, -1, -1};
  // k=3 would give a two-base hairpin. Enlarging the split range by one
  // must not admit it, for either V or WMB start/binary rules.
  for (const auto rule : {R::WIP_BASEPAIR_V, R::WIP_BASEPAIR_WMB,
                           R::WIP_SPLIT_V, R::WIP_SPLIT_WMB}) {
    for (const auto &source : deductions::generate_exact_basic_deductions(parent, rule, ctx))
      require(source.split.k <= 2, "WIP admitted a sub-TURN branch");
  }
  require(!deductions::exact_wip_start_split_allows(parent, 3, ctx),
          "WIP start guard admitted a sub-TURN branch");
}

void check_root_consumers() {
  const std::string seq = "GAGAAACAAGAAGAAAAAAAAACCAAAAC";
  const std::string g = "(.(...)..(............)).....";
  const std::string target = "............(...............)";
  auto sequence = seq, scaffold = g;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold, cparty::scfg2::weights::RuntimeOptions{0});
  const std::vector<engine::TracebackPair> pairs{
      {13, 29, engine::TracebackPairLane::Square}};
  const auto fixed = replay::run_w_final_exact_output_constrained_viterbi_traceback<true>(owner, pairs);
  require(close(owner.partition_free_energy_from_log_weight(fixed.run.chart.get(fixed.root)), 13.58),
          "minimal WIP hairpin did not reach the completed 29-nt target at 13.58 kcal/mol");
  require(engine::same_traceback_pair_set(fixed.generated_pairs, pairs),
          "fixed-output traceback lost the target pair");

  // Public consumers retain their usual dangles=2 setting.
  api::ExactSession exact(seq, g);
  api::BeamSession finite(seq, g, 100);
  const auto evaluated = exact.evaluate_energy(target);
  const auto exact_p = exact.probability(target);
  const auto beam_p = finite.probability(target);
  require(evaluated.ok && evaluated.reachable && exact_p.ok && exact_p.reachable &&
              beam_p.ok && beam_p.target_included && exact_p.probability > 0.0 &&
              close(exact_p.log_probability, beam_p.log_probability),
          "completed WIP hairpin target missing from public exact/beam consumers");
  const auto exact_out = exact.outside();
  const auto beam_out = finite.outside();
  require(exact_out.ok && beam_out.ok && close(exact_out.log_partition_weight, beam_out.log_partition_weight),
          "completed witness outside log Z mismatch");
  const auto exact_profile = api::exact_structural_profiles(seq, g, {1, 2, 4});
  const auto beam_profile = api::beam_structural_profiles(seq, g, {1, 2, 4}, api::BeamConfig{100});
  require(exact_profile.ok && beam_profile.ok, "completed witness profiles failed");
  for (std::size_t pos = 1; pos <= seq.size(); ++pos) {
    require(close(exact_out.unpaired_probabilities[pos], beam_out.unpaired_probabilities[pos]),
            "completed witness outside marginal mismatch");
    for (std::size_t c = 0; c < api::structural_context_count; ++c)
      require(close(exact_profile.contexts[pos][c], beam_profile.contexts[pos][c]),
              "completed witness profile context mismatch");
    require(close(1.0 - beam_profile.contexts[pos][0], beam_out.unpaired_probabilities[pos]),
            "completed witness profile/pair emission mismatch");
  }
  for (std::size_t w = 0; w < exact_profile.accessibility.size(); ++w)
    for (std::size_t pos = 1; pos < exact_profile.accessibility[w].size(); ++pos)
      require(close(exact_profile.accessibility[w][pos], beam_profile.accessibility[w][pos]),
              "completed witness accessibility mismatch");
}
}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_DirksPierce09.par", 1);
  check_local("AAGAAAC", ".......", 2, 7, 4.37, R::WIP_BASEPAIR_V, 3);
  check_local("AAGAAAC", ".......", 3, 7, 4.25, R::WIP_BASE_V, -1);
  check_local("AAGAAAAC", "........", 2, 8, 3.56, R::WIP_BASEPAIR_V, 3);
  check_local("AAGAAACAGAAAC", "..(...).(...)", 3, 13, 8.62, R::WIP_SPLIT_V, 9);
  check_too_short();
  check_root_consumers();
  std::cout << "scfg2_wip_hairpin_boundary=ok\n";
}

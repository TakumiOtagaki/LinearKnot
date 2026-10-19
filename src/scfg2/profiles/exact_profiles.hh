#ifndef CPARTY_SCFG2_PROFILES_EXACT_PROFILES_HH
#define CPARTY_SCFG2_PROFILES_EXACT_PROFILES_HH

#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/profiles/loop_context.hh"
#include "scfg2/profiles/profile_assembly.hh"
#include "scfg2/replay/w_final_exact_outside.hh"

#include <cmath>
#include <string>
#include <vector>

namespace cparty::scfg2::profiles {

// Exact PK-CapR/Raccess on numeric alpha/beta. Regenerated emissions go
// directly to output-sized difference arrays; only state chain steps survive.
inline ProfileResult exact_profiles(const std::string &sequence,
                                    const std::string &scaffold,
                                    const std::vector<cand_pos_t> &window_lengths,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure,
    bool by_loop = false) {
  using Semiring = semirings::LogSumExpSemiring;
  std::string seq = sequence;
  std::string g = scaffold;
  weights::RuntimeSnapshot pf(seq, g);
  replay::WFinalLogWeightContext ctx(pf, policy);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto is_unpaired = [&](cand_pos_t pos) { return ctx.is_unpaired(pos); };

  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  const auto chart = replay::run_w_final_exact_dp_schedule<Semiring>(schedule, ctx,
      [](const deductions::Deduction &) { return true; },
      [](const auto &...) {});
  const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
  const double root_value = chart.get(root);
  if (!std::isfinite(root_value)) {
    ProfileResult out;
    out.error = "scaffold admits no structure";
    return out;
  }
  engine::AllLocalWeightModel<replay::WFinalLogWeightContext> weight(&ctx);
  ProfileAccumulator acc(n, window_lengths, by_loop);
  replay::run_w_final_exact_regenerated_outside(schedule, ctx, chart, root,
      [&](const deductions::Deduction &deduction, double parent_beta) {
    const auto local = weight.local_weight_for(deduction);
    if (weights::is_zero(local)) return;
    const auto emission = unpaired_emission(deduction, is_unpaired);
    if (emission.count == 0) return;
    const double contribution = replay::runtime_deduction_contribution<Semiring>(
        chart, deduction, local);
    const double probability = std::exp(parent_beta + contribution - root_value);
    if (!emission.chain_step) {
      acc.add_whole_runs(emission, probability);
      return;
    }
    const auto &parent = deduction.parent;
    auto &step = acc.chain_step(static_cast<int>(parent.nonterminal), parent.i,
                                parent.j, emission.runs[0].context);
    step.log_entry = contribution;
    step.log_step = weights::log_weight(local);
    step.has_prev = deduction.child_count > 0;
    step.emit_probability = probability;
    step.log_state = chart.get(parent);
    step.state_adjoint = std::exp(parent_beta + step.log_state - root_value);
  });
  auto out = acc.finish(root_value, window_lengths);
  out.ensemble_free_energy =
      pf.partition_free_energy_from_log_weight(root_value);
  out.rt_kcal_mol = pf.partition_free_energy_from_log_weight(-1.0) -
                    pf.partition_free_energy_from_log_weight(0.0);
  return out;
}

}  // namespace cparty::scfg2::profiles

#endif  // CPARTY_SCFG2_PROFILES_EXACT_PROFILES_HH

#ifndef CPARTY_SCFG2_ENGINE_BASIC_EXACT_INSIDE_RUNNER_HH
#define CPARTY_SCFG2_ENGINE_BASIC_EXACT_INSIDE_RUNNER_HH

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/engine/local_weight_model.hh"
#include "scfg2/engine/run_schedule.hh"
#include "scfg2/semirings/semiring.hh"

namespace cparty::scfg2::engine {

template <class Semiring, class PairabilityContext>
Chart<typename Semiring::value_type> run_basic_exact_inside(
    cand_pos_t n,
    const LocalWeightModel &weight_model,
    const PairabilityContext &pairability_ctx) {
  const auto schedule = build_pkfree_basic_exact_schedule(n, pairability_ctx);
  return run_schedule<Semiring>(
      schedule,
      [&](const deductions::ItemKey &item) {
        return deductions::generate_exact_basic_deductions_for_item(item, pairability_ctx);
      },
      [&](const deductions::Deduction &deduction) {
        return weight_model.weight_for(deduction);
      });
}

template <class PairabilityContext>
ViterbiRunResult<typename semirings::MaxProductSemiring::value_type>
run_basic_exact_viterbi(cand_pos_t n,
                        const LocalWeightModel &weight_model,
                        const PairabilityContext &pairability_ctx) {
  const auto schedule = build_pkfree_basic_exact_schedule(n, pairability_ctx);
  return run_schedule_with_choices<semirings::MaxProductSemiring>(
      schedule,
      [&](const deductions::ItemKey &item) {
        return deductions::generate_exact_basic_deductions_for_item(item, pairability_ctx);
      },
      [&](const deductions::Deduction &deduction) {
        return weight_model.weight_for(deduction);
      });
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_BASIC_EXACT_INSIDE_RUNNER_HH

#ifndef CPARTY_SCFG2_REPLAY_W_FINAL_EXACT_OUTSIDE_HH
#define CPARTY_SCFG2_REPLAY_W_FINAL_EXACT_OUTSIDE_HH

#include "scfg2/engine/outside.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <utility>

namespace cparty::scfg2::replay {

// Derived reverse consumer of the exact forward deduction surface. Preserve
// the same span invariant, local log weights and boundary child identities;
// only the current parent's candidates are materialized, never all edges.
template <class Context, class Observer = engine::IgnoreOutsideDeduction>
engine::Chart<double> run_w_final_exact_regenerated_outside(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx, const engine::Chart<double> &inside,
    const deductions::ItemKey &root, Observer observe = {}) {
  engine::AllLocalWeightModel<Context> weight(&ctx);
  return engine::run_regenerated_outside<semirings::LogSumExpSemiring>(
      schedule, inside, root,
      [&](const deductions::ItemKey &parent) {
        auto candidates = deductions::generate_exact_basic_deductions_for_item(parent, ctx);
        for (const auto &d : candidates) engine::require_legacy_span_schedule_support(d);
        return candidates;
      },
      [&](const deductions::Deduction &d) { return weight.local_weight_for(d); },
      engine::RuntimeDeductionChildItem{}, std::move(observe), inside.dense_length());
}

}  // namespace cparty::scfg2::replay
#endif

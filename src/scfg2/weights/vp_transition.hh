#ifndef CPARTY_SCFG2_WEIGHTS_VP_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_VP_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_vp.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_vp(const deductions::Deduction &deduction,
                                   const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::VP_WI_CASE1:
    case grammar::RuleId::VP_WI_CASE2:
    case grammar::RuleId::VP_WI_CASE3:
      return ctx.vp_wi_case_weight();
    case grammar::RuleId::VP_STACK:
      return ctx.vp_stack_weight(deduction.parent.i, deduction.parent.j);
    case grammar::RuleId::VP_INTERNAL_LOOP:
      if (deduction.child_count != 1) {
        return impossible_weight<Context>();
      }
      return ctx.vp_internal_loop_weight(deduction.parent.i,
                                         deduction.parent.j,
                                         deduction.children[0].item.i,
                                         deduction.children[0].item.j);
    case grammar::RuleId::VP_WIP_VP_LEFT:
    case grammar::RuleId::VP_VP_WIP_RIGHT:
    case grammar::RuleId::VP_WIP_VPR:
    case grammar::RuleId::VP_VPL_WIP:
      return ctx.vp_split_weight();
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_VP_TRANSITION_HH

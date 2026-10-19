#ifndef CPARTY_SCFG2_WEIGHTS_WM_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_WM_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_wm.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_wm(const deductions::Deduction &deduction,
                                   const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::WM_START_V:
      return ctx.start_v_weight(deduction.parent.i, deduction.parent.j, deduction.split.k);
    case grammar::RuleId::WM_START_WMB:
      return ctx.start_wmb_weight(deduction.parent.i, deduction.parent.j, deduction.split.k);
    case grammar::RuleId::WM_SPLIT_V:
      return ctx.split_v_weight(deduction.split.k, deduction.parent.j);
    case grammar::RuleId::WM_SPLIT_WMB:
      return ctx.split_wmb_weight(deduction.split.k, deduction.parent.j);
    case grammar::RuleId::WM_EXTEND_UNPAIRED:
      return ctx.extend_unpaired_weight();
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_WM_TRANSITION_HH

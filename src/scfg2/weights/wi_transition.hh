#ifndef CPARTY_SCFG2_WEIGHTS_WI_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_WI_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_wi.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_wi(const deductions::Deduction &deduction,
                                   const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::WI_BASE_SINGLE:
      return ctx.wi_base_single_weight();
    case grammar::RuleId::WI_SPLIT_V:
      return ctx.wi_split_v_weight();
    case grammar::RuleId::WI_SPLIT_WMB:
      return ctx.wi_split_wmb_weight();
    case grammar::RuleId::WI_EXTEND_UNPAIRED:
      return ctx.wi_extend_unpaired_weight();
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_WI_TRANSITION_HH

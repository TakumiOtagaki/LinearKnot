#ifndef CPARTY_SCFG2_WEIGHTS_WIP_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_WIP_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_wip.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_wip(const deductions::Deduction &deduction,
                                    const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::WIP_BASE_V:
      return ctx.wip_base_v_weight();
    case grammar::RuleId::WIP_BASE_WMB:
      return ctx.wip_base_wmb_weight();
    case grammar::RuleId::WIP_SPLIT_V:
      return ctx.wip_split_v_weight();
    case grammar::RuleId::WIP_SPLIT_WMB:
      return ctx.wip_split_wmb_weight();
    case grammar::RuleId::WIP_BASEPAIR_V:
      return ctx.wip_basepair_v_weight(deduction.split.k - deduction.parent.i);
    case grammar::RuleId::WIP_BASEPAIR_WMB:
      return ctx.wip_basepair_wmb_weight(deduction.split.k - deduction.parent.i);
    case grammar::RuleId::WIP_EXTEND_UNPAIRED:
      return ctx.wip_extend_unpaired_weight();
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_WIP_TRANSITION_HH

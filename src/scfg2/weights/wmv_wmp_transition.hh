#ifndef CPARTY_SCFG2_WEIGHTS_WMV_WMP_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_WMV_WMP_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_wmv_wmp.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_wmv_wmp(const deductions::Deduction &deduction,
                                        const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::WMv_STEM_V:
      return ctx.wmv_stem_weight(deduction.parent.i, deduction.parent.j);
    case grammar::RuleId::WMp_STEM_WMB:
      return ctx.wmp_stem_weight(deduction.parent.i, deduction.parent.j);
    case grammar::RuleId::WMv_EXTEND_UNPAIRED:
      return ctx.wmv_extend_unpaired_weight();
    case grammar::RuleId::WMp_EXTEND_UNPAIRED:
      return ctx.wmp_extend_unpaired_weight();
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_WMV_WMP_TRANSITION_HH

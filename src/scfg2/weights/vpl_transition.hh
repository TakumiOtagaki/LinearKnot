#ifndef CPARTY_SCFG2_WEIGHTS_VPL_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_VPL_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_vpl.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_vpl(const deductions::Deduction &deduction,
                                    const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::VPL_SPLIT_VP:
      return ctx.vpl_split_vp_weight(deduction.split.k - deduction.parent.i);
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_VPL_TRANSITION_HH

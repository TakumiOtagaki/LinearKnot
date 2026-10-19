#ifndef CPARTY_SCFG2_WEIGHTS_WMBW_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_WMBW_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_wmbw.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_wmbw(const deductions::Deduction &deduction,
                                     const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::WMBW_SPLIT_WMBP_WI:
      return ctx.wmbw_split_wmbp_wi_weight();
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_WMBW_TRANSITION_HH

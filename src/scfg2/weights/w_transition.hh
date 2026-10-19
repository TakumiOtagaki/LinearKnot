#ifndef CPARTY_SCFG2_WEIGHTS_W_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_W_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_w.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_w(const deductions::Deduction &deduction,
                                  const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::W_EXTEND_UNPAIRED:
      return ctx.w_extend_unpaired_weight();
    case grammar::RuleId::W_SPLIT_V:
      return ctx.w_split_v_weight(deduction.split.k, deduction.parent.j);
    case grammar::RuleId::W_SPLIT_WMB:
      return ctx.w_split_wmb_weight();
    case grammar::RuleId::W_ROOT_WMB_DECOMPOSITION:
      return ctx.w_root_wmb_decomposition_weight();
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_W_TRANSITION_HH

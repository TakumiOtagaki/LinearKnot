#ifndef CPARTY_SCFG2_WEIGHTS_WMB_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_WMB_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/deductions/wmb_geometry.hh"
#include "scfg2/weights/context_wmb.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_wmb(const deductions::Deduction &deduction,
                                    const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::WMB_SPLIT_BE_WMBP_WI: {
      // Weight-layer anchors keep only the BE endpoints needed by the local
      // multiplier boundary; the full anchored BE support stays outside here.
      const WmbSplitBeWmbpWiWeightAnchors anchors{
          .be_i = deduction.split.p,
          .be_j = deduction.parent.j,
          .be_ip = deduction.split.q,
          .be_jp = deduction.split.l,
      };
      return ctx.wmb_split_be_wmbp_wi_weight(anchors);
    }
    case grammar::RuleId::WMB_DIRECT_WMBP:
      return ctx.wmb_direct_wmbp_weight();
    case grammar::RuleId::WMB_EMPTY:
      return ctx.wmb_empty_weight();
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_WMB_TRANSITION_HH

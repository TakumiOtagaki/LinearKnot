#ifndef CPARTY_SCFG2_WEIGHTS_WMBP_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_WMBP_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/deductions/wmb_geometry.hh"
#include "scfg2/weights/context_wmbp.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_wmbp(const deductions::Deduction &deduction,
                                     const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP: {
      // Shared-BE rules expose the anchor geometry to the weight layer even
      // when the exact local factor currently remains a compatibility constant.
      const WmbpSharedBeVpAnchors anchors{
          .split_k = deduction.split.k,
          .be_i = deduction.split.p,
          .be_j = deduction.split.q,
          .be_ip = deduction.split.l,
      };
      return ctx.wmbp_split_be_wmbp_vp_weight(anchors);
    }
    case grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP: {
      const WmbpSharedBeVpAnchors anchors{
          .split_k = deduction.split.k,
          .be_i = deduction.split.p,
          .be_j = deduction.split.q,
          .be_ip = deduction.split.l,
      };
      return ctx.wmbp_split_be_wmbw_vp_weight(anchors);
    }
    case grammar::RuleId::WMBP_DIRECT_VP:
      return ctx.wmbp_direct_vp_weight(deduction.parent.i, deduction.parent.j);
    case grammar::RuleId::WMBP_SPLIT_BE_WI_VP: {
      const deductions::WmbpWiVpAnchors anchors{
          .parent_i = deduction.parent.i,
          .be_ip = deduction.split.k,
          .be_jp = deduction.split.l,
          .vp_i = deduction.split.p,
          .vp_j = deduction.parent.j,
      };
      return ctx.wmbp_split_be_wi_vp_weight(anchors);
    }
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_WMBP_TRANSITION_HH

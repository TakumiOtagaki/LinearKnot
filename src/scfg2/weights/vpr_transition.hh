#ifndef CPARTY_SCFG2_WEIGHTS_VPR_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_VPR_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_vpr.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_vpr(const deductions::Deduction &deduction,
                                    const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::VPR_SPLIT_VP_WIP:
      return ctx.vpr_split_vp_wip_weight();
    case grammar::RuleId::VPR_SPLIT_VP_BASEPAIR:
      // Right unpaired run [k+1,j] (CParty SI eq. vi; TheCOBRALab/CParty#5).
      return ctx.vpr_split_vp_basepair_weight(deduction.parent.j - deduction.split.k);
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_VPR_TRANSITION_HH

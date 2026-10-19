#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_VPR_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_VPR_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class VPRTransitionContext {
 public:
  virtual ~VPRTransitionContext() = default;

  virtual double vpr_split_vp_wip_weight() const = 0;
  virtual double vpr_split_vp_basepair_weight(cand_pos_t right_unpaired) const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_VPR_HH

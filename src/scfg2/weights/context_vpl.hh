#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_VPL_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_VPL_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class VPLTransitionContext {
 public:
  virtual ~VPLTransitionContext() = default;

  virtual double vpl_split_vp_weight(cand_pos_t left_unpaired) const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_VPL_HH

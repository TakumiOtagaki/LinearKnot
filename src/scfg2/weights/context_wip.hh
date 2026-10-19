#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_WIP_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_WIP_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class WIPTransitionContext {
 public:
  virtual ~WIPTransitionContext() = default;

  virtual double wip_base_v_weight() const = 0;
  virtual double wip_base_wmb_weight() const = 0;
  virtual double wip_split_v_weight() const = 0;
  virtual double wip_split_wmb_weight() const = 0;
  virtual double wip_basepair_v_weight(cand_pos_t left_unpaired) const = 0;
  virtual double wip_basepair_wmb_weight(cand_pos_t left_unpaired) const = 0;
  virtual double wip_extend_unpaired_weight() const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_WIP_HH

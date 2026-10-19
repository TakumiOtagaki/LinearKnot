#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_WI_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_WI_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class WITransitionContext {
 public:
  virtual ~WITransitionContext() = default;

  virtual double wi_base_single_weight() const = 0;
  virtual double wi_split_v_weight() const = 0;
  virtual double wi_split_wmb_weight() const = 0;
  virtual double wi_extend_unpaired_weight() const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_WI_HH

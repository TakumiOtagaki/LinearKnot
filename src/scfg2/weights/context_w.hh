#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_W_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_W_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class WTransitionContext {
 public:
  virtual ~WTransitionContext() = default;

  virtual double w_extend_unpaired_weight() const = 0;
  virtual double w_split_v_weight(cand_pos_t i, cand_pos_t j) const = 0;
  virtual double w_split_wmb_weight() const = 0;
  virtual double w_root_wmb_decomposition_weight() const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_W_HH

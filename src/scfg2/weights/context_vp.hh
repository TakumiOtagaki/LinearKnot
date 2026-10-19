#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_VP_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_VP_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class VPTransitionContext {
 public:
  virtual ~VPTransitionContext() = default;

  virtual double vp_wi_case_weight() const = 0;
  virtual double vp_stack_weight(cand_pos_t i, cand_pos_t j) const = 0;
  virtual double vp_internal_loop_weight(cand_pos_t i,
                                         cand_pos_t j,
                                         cand_pos_t k,
                                         cand_pos_t l) const = 0;
  virtual double vp_split_weight() const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_VP_HH

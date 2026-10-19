#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_V_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_V_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class VTransitionContext {
 public:
  virtual ~VTransitionContext() = default;

  virtual double hairpin_weight(cand_pos_t i, cand_pos_t j) const = 0;
  virtual double internal_weight(cand_pos_t i,
                                 cand_pos_t j,
                                 cand_pos_t k,
                                 cand_pos_t l) const = 0;
  virtual double vm_closing_weight(cand_pos_t i, cand_pos_t j) const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_V_HH

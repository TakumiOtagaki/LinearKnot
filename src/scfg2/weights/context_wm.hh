#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_WM_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_WM_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class WMTransitionContext {
 public:
  virtual ~WMTransitionContext() = default;

  virtual double start_v_weight(cand_pos_t i, cand_pos_t j, cand_pos_t k) const = 0;
  virtual double start_wmb_weight(cand_pos_t i, cand_pos_t j, cand_pos_t k) const = 0;
  virtual double split_v_weight(cand_pos_t k, cand_pos_t j) const = 0;
  virtual double split_wmb_weight(cand_pos_t k, cand_pos_t j) const = 0;
  virtual double extend_unpaired_weight() const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_WM_HH

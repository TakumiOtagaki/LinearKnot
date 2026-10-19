#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_VM_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_VM_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class VMTransitionContext {
 public:
  virtual ~VMTransitionContext() = default;

  virtual double split_wm_wmv_weight(cand_pos_t i, cand_pos_t j) const = 0;
  virtual double split_wm_wmp_weight(cand_pos_t i, cand_pos_t j) const = 0;
  virtual double split_wmp_base_weight(cand_pos_t i,
                                       cand_pos_t j,
                                       cand_pos_t left_unpaired) const = 0;
  virtual double vm_base_unpaired_step_weight() const = 0;
  virtual double vm_base_gap_weight(cand_pos_t left_unpaired) const = 0;
  virtual double vm_wrapper_weight() const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_VM_HH

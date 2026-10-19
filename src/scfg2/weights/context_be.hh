#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_BE_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_BE_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class BETransitionContext {
 public:
  virtual ~BETransitionContext() = default;

  virtual double be_base_samepair_weight() const = 0;
  virtual double be_stack_weight(cand_pos_t i, cand_pos_t j) const = 0;
  virtual double be_internal_loop_weight(cand_pos_t i,
                                         cand_pos_t j,
                                         cand_pos_t k,
                                         cand_pos_t l) const = 0;
  virtual double be_split_wip_wip_weight() const = 0;
  virtual double be_split_wip_basepair_weight(cand_pos_t right_unpaired) const = 0;
  virtual double be_split_basepair_wip_weight(cand_pos_t left_unpaired) const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_BE_HH

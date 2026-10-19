#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_WMV_WMP_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_WMV_WMP_HH

#include "scfg2/types.hh"

namespace cparty::scfg2::weights {

class WMvWMpTransitionContext {
 public:
  virtual ~WMvWMpTransitionContext() = default;

  virtual double wmv_stem_weight(cand_pos_t i, cand_pos_t j) const = 0;
  virtual double wmp_stem_weight(cand_pos_t i, cand_pos_t j) const = 0;
  virtual double wmv_extend_unpaired_weight() const = 0;
  virtual double wmp_extend_unpaired_weight() const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_WMV_WMP_HH

#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_WMB_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_WMB_HH

#include "scfg2/deductions/wmb_geometry.hh"

namespace cparty::scfg2::weights {

struct WmbSplitBeWmbpWiWeightAnchors {
  cand_pos_t be_i = -1;
  cand_pos_t be_j = -1;
  cand_pos_t be_ip = -1;
  cand_pos_t be_jp = -1;
};

class WMBTransitionContext {
 public:
  virtual ~WMBTransitionContext() = default;

  virtual double wmb_split_be_wmbp_wi_weight(
      const WmbSplitBeWmbpWiWeightAnchors &anchors) const = 0;
  virtual double wmb_direct_wmbp_weight() const = 0;
  virtual double wmb_empty_weight() const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_WMB_HH

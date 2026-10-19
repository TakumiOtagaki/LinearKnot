#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_WMBP_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_WMBP_HH

#include "scfg2/deductions/wmb_geometry.hh"

namespace cparty::scfg2::weights {

struct WmbpSharedBeVpAnchors {
  cand_pos_t split_k = -1;
  cand_pos_t be_i = -1;
  cand_pos_t be_j = -1;
  cand_pos_t be_ip = -1;
};

class WMBPTransitionContext {
 public:
  virtual ~WMBPTransitionContext() = default;

  virtual double wmbp_split_be_wmbp_vp_weight(
      const WmbpSharedBeVpAnchors &anchors) const = 0;
  virtual double wmbp_split_be_wmbw_vp_weight(
      const WmbpSharedBeVpAnchors &anchors) const = 0;
  virtual double wmbp_direct_vp_weight(cand_pos_t i, cand_pos_t j) const = 0;
  virtual double wmbp_split_be_wi_vp_weight(
      const deductions::WmbpWiVpAnchors &anchors) const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_WMBP_HH

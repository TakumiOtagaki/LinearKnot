#ifndef CPARTY_SCFG2_WEIGHTS_CONTEXT_WMBW_HH
#define CPARTY_SCFG2_WEIGHTS_CONTEXT_WMBW_HH

namespace cparty::scfg2::weights {

class WMBWTransitionContext {
 public:
  virtual ~WMBWTransitionContext() = default;

  virtual double wmbw_split_wmbp_wi_weight() const = 0;
};

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_CONTEXT_WMBW_HH

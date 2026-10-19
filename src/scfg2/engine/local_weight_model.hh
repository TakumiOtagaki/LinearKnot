#ifndef CPARTY_SCFG2_ENGINE_LOCAL_WEIGHT_MODEL_HH
#define CPARTY_SCFG2_ENGINE_LOCAL_WEIGHT_MODEL_HH

#include "scfg2/deductions/deduction.hh"

namespace cparty::scfg2::engine {

class LocalWeightModel {
 public:
  virtual ~LocalWeightModel() = default;
  virtual double weight_for(const deductions::Deduction &deduction) const = 0;
};

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_LOCAL_WEIGHT_MODEL_HH

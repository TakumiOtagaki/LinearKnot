#ifndef CPARTY_SCFG2_API_TYPES_HH
#define CPARTY_SCFG2_API_TYPES_HH

#include "scfg2/grammar/model_policy.hh"

namespace cparty::scfg2::api {

using ::cparty::scfg2::WmbClosurePolicy;

// Public base-pair coordinates are 1-based and inclusive. Pair-list inputs
// accept either endpoint order; output pairs are normalized to i < j.
struct BasePair {
  int i = -1;
  int j = -1;
};

struct PairProbability {
  int i = -1;
  int j = -1;
  double round_probability = 0.0;
  double square_probability = 0.0;
  double union_probability = 0.0;
};

}  // namespace cparty::scfg2::api

#endif  // CPARTY_SCFG2_API_TYPES_HH

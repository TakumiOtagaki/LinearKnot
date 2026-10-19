#ifndef CPARTY_SCFG2_API_POSTERIOR_DECODING_HH
#define CPARTY_SCFG2_API_POSTERIOR_DECODING_HH

#include "scfg2/api/beam.hh"

namespace cparty::scfg2::api {

struct PosteriorDecodingResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool reachable = false;
  std::string error;
  std::size_t beam_size = 0;  // zero for the exact computation
  double gamma = 1.0;
  // Expected gain for the selected objective.
  double expected_accuracy = std::numeric_limits<double>::quiet_NaN();
  double log_partition_weight = -std::numeric_limits<double>::infinity();
  std::string scaffold_structure;
  std::string generated_structure;
  std::string union_structure;
  std::vector<BasePair> generated_pairs;
  // Beam: peak live selected-DAG nodes; terms/factors visited in chosen root.
  // Exact: all frozen-forest states/terms/factors in the maximum sweep.
  std::size_t decoding_states = 0;
  std::size_t decoding_terms = 0;
  std::size_t decoding_factors = 0;
  // Total includes inside and outside. Decoding timing starts after marginals
  // are available and includes reward assembly, maximization and traceback.
  double elapsed_seconds = 0.0;
  double decoding_seconds = 0.0;
  double expected_true_positive_pairs = 0.0;
  // Over every unordered position pair i<j, including impossible pairs.
  double expected_true_negative_pairs = 0.0;
};

}  // namespace cparty::scfg2::api

#endif

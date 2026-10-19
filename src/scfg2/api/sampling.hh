#ifndef CPARTY_SCFG2_API_SAMPLING_HH
#define CPARTY_SCFG2_API_SAMPLING_HH

#include "scfg2/api/beam.hh"

#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace cparty::scfg2::api {

// Structures G' drawn from the conditional ensemble given the scaffold G
// (exact) or from the retained beam DAG (`linearknot sample`), with probability
// weight / Z of that DAG.
struct SampleResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  std::string error;
  std::size_t beam_size = 0;  // zero for the exact computation
  double log_partition_weight = -std::numeric_limits<double>::infinity();
  double ensemble_free_energy = std::numeric_limits<double>::quiet_NaN();
  // One entry per sample, in draw order.
  std::vector<std::vector<BasePair>> generated_pairs;
  std::vector<std::string> union_structures;
  double elapsed_seconds = 0.0;
};

struct SampledStructure {
  std::vector<BasePair> generated_pairs;
  std::string union_structure;
};

struct SampleStreamResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  std::string error;
  std::size_t beam_size = 0;
  double log_partition_weight = -std::numeric_limits<double>::infinity();
  double ensemble_free_energy = std::numeric_limits<double>::quiet_NaN();
  // Includes the structure passed to a consumer that returns false/throws.
  std::size_t samples_emitted = 0;
  bool stopped_by_consumer = false;
  double elapsed_seconds = 0.0;  // Includes consumer work.
};

using SampleConsumer = std::function<bool(SampledStructure &&)>;

// Synchronous, one structure at a time. Move from the argument to retain it;
// otherwise it is released on return. false stops after this structure.
// A later error cannot retract an already-delivered prefix; inspect ok/error.
// Forest/sampler storage remains; this removes completed-output retention.
SampleStreamResult exact_sample_each(const std::string &sequence,
                                     const std::string &scaffold_structure,
                                     std::size_t count, std::uint64_t seed,
                                     const SampleConsumer &consumer,
                                     WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);
SampleStreamResult beam_sample_each(const std::string &sequence,
                                    const std::string &scaffold_structure,
                                    std::size_t count, std::uint64_t seed,
                                    const SampleConsumer &consumer,
                                    BeamConfig config = {});

SampleResult exact_sample(const std::string &sequence,
                          const std::string &scaffold_structure,
                          std::size_t count,
                          std::uint64_t seed,
                          WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);

SampleResult beam_sample(const std::string &sequence,
                         const std::string &scaffold_structure,
                         std::size_t count,
                         std::uint64_t seed,
                         BeamConfig config = {});

}  // namespace cparty::scfg2::api

#endif  // CPARTY_SCFG2_API_SAMPLING_HH

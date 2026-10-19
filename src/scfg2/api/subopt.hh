#ifndef CPARTY_SCFG2_API_SUBOPT_HH
#define CPARTY_SCFG2_API_SUBOPT_HH

#include "scfg2/api/beam.hh"

#include <functional>

namespace cparty::scfg2::api {

struct SuboptConfig {
  std::size_t top_k = 10;  // positive; fewer results if support/band exhausts
  std::optional<double> delta_energy;  // finite, nonnegative kcal/mol
};

enum class SuboptStopReason { TopK, DeltaEnergy, Exhausted, Consumer };
inline const char *subopt_stop_reason_name(SuboptStopReason reason) {
  switch (reason) {
    case SuboptStopReason::TopK: return "top_k";
    case SuboptStopReason::DeltaEnergy: return "delta_energy";
    case SuboptStopReason::Exhausted: return "exhausted";
    case SuboptStopReason::Consumer: return "consumer";
  }
  return "unknown";
}

struct SuboptStructure {
  double log_weight = -std::numeric_limits<double>::infinity();
  double energy = std::numeric_limits<double>::quiet_NaN();
  std::vector<BasePair> generated_pairs;  // ordered, one-based
  std::string generated_structure;
  std::string union_structure;
};

struct SuboptResult {
  bool ok = false;
  bool reachable = false;
  std::string error;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  std::size_t beam_size = 0;  // zero means exact
  SuboptStopReason stop_reason = SuboptStopReason::Exhausted;
  double log_partition_weight = -std::numeric_limits<double>::infinity();
  double mfe_energy = std::numeric_limits<double>::quiet_NaN();
  std::vector<SuboptStructure> structures;  // nondecreasing energy
  std::size_t graph_states = 0, graph_terms = 0, graph_factors = 0;
  std::size_t candidates_created = 0, derivations_cached = 0;
  double elapsed_seconds = 0.0;
};

struct SuboptStreamResult {
  bool ok = false;
  bool reachable = false;
  std::string error;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  std::size_t beam_size = 0;
  SuboptStopReason stop_reason = SuboptStopReason::Exhausted;
  double log_partition_weight = -std::numeric_limits<double>::infinity();
  double mfe_energy = std::numeric_limits<double>::quiet_NaN();
  // Includes the structure passed to a consumer that returns false/throws.
  std::size_t structures_emitted = 0;
  std::size_t graph_states = 0, graph_terms = 0, graph_factors = 0;
  std::size_t candidates_created = 0, derivations_cached = 0;
  double elapsed_seconds = 0.0;  // Includes consumer work.
};

using SuboptConsumer = std::function<bool(SuboptStructure &&)>;

// Synchronous ordered delivery; false stops after this structure, without
// requesting the next rank. Move the argument to retain it. Earlier callbacks
// cannot be retracted on a later error. Inspect ok/error and emitted count.
// Forest, rank/heap caches and exact ambiguity keys still consume memory;
// this removes the completed-output vector, not all K-dependent storage.
SuboptStreamResult exact_subopt_each(const std::string &sequence, const std::string &scaffold,
                                     const SuboptConsumer &consumer, SuboptConfig config = {},
                                     WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);
SuboptStreamResult beam_subopt_each(const std::string &sequence, const std::string &scaffold,
                                    const SuboptConsumer &consumer, SuboptConfig config = {},
                                    BeamConfig beam = {});

// One common lazy k-best consumer of the shared deduction forest. Exact
// ranks the full grammar; beam ranks the frozen INSIDE-ranked forest, not
// the separately LogMax-pruned BeamSession::viterbi support. Delta is relative
// to this same forest's optimum (inclusive, 1e-9 kcal/mol rounding tolerance).
// Equal-energy ties are deterministic within a forest, with no lexicographic
// promise across different forests. Duplicate physical pair maps are errors:
// no deduplication hides a violation of the grammar's unambiguity contract.
SuboptResult exact_subopt(const std::string &sequence, const std::string &scaffold,
                          SuboptConfig config = {},
                          WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);
SuboptResult beam_subopt(const std::string &sequence, const std::string &scaffold,
                         SuboptConfig config = {}, BeamConfig beam = {});

}  // namespace cparty::scfg2::api
#endif

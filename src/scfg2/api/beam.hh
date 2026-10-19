#ifndef CPARTY_SCFG2_API_BEAM_HH
#define CPARTY_SCFG2_API_BEAM_HH

#include "scfg2/api/exact.hh"
#include "scfg2/beam/ranking_config.hh"
#include "scfg2/beam/memory_observation.hh"

#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace cparty::scfg2::beam { struct PreparedRanking; }

namespace cparty::scfg2::api {

using BeamRanking = beam::RankingMode;
using OuterPilotExponent = beam::PilotExponent;

struct BeamConfig {
  std::size_t beam_size = 100;
  std::size_t max_hairpin_unpaired = 30;
  std::size_t max_case3_ancestor_depth = 30;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  BeamRanking ranking = BeamRanking::Outer;
  OuterPilotExponent outer_pilot_exponent = OuterPilotExponent::TwoThirds;
  // Optional serial diagnostic callback. Cached calls need not repeat events.
  // Phase-boundary events are always emitted when enabled. A positive interval
  // additionally samples the forward scratch containers every N columns.
  beam::MemoryObserver memory_observer = {};
  std::size_t memory_column_interval = 0;
};

struct BeamFamilyWork {
  std::string nonterminal;
  // Cumulative evaluations before pruning, not live/retained state counts.
  std::size_t item_visits = 0;
  std::size_t positive_items = 0;
  std::size_t frontier_candidates = 0;
  std::size_t frontier_selected = 0;
};

struct BeamWorkStats {
  // Cumulative main-forward work; excludes forcing/pilot/outer. For actual
  // live state counts use MemoryObservation at main/forward_end.
  std::size_t total_item_visits = 0;
  std::size_t total_positive_items = 0;
  std::size_t total_frontier_candidates = 0;
  std::size_t total_frontier_selected = 0;
  std::size_t vp_candidates = 0;
  std::size_t vp_selected = 0;
  std::size_t wmbp_candidates = 0;
  std::size_t wmbp_selected = 0;
  std::size_t case3_candidates = 0;
  std::size_t case3_selected = 0;
  std::size_t binary_join_probes = 0;
  std::size_t extension_probes = 0;
  std::size_t v_parent_probes = 0;
  std::size_t vp_loop_parent_probes = 0;
  // Unpaired positions with more branching scaffold ancestors than
  // BeamConfig::max_case3_ancestor_depth. Zero means that limit removed no
  // CASE3 crossing candidate from this run.
  std::size_t case3_truncated_endpoints = 0;
  std::vector<BeamFamilyWork> families;
};

struct BeamInsideResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool reachable = false;
  std::string error;
  std::size_t beam_size = 0;
  double log_partition_weight =
      -std::numeric_limits<double>::infinity();
  // Convenience projection. It may underflow to zero on long sequences;
  // log_partition_weight is the canonical value.
  double partition_weight = 0.0;
  double ensemble_free_energy = std::numeric_limits<double>::quiet_NaN();
  double elapsed_seconds = 0.0;
  BeamWorkStats work;
};

struct BeamOutsideResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool reachable = false;
  std::string error;
  std::size_t beam_size = 0;
  double log_partition_weight =
      -std::numeric_limits<double>::infinity();
  double partition_weight = 0.0;
  double ensemble_free_energy = std::numeric_limits<double>::quiet_NaN();
  std::vector<PairProbability> pair_probabilities;
  // One-based: index zero is unused.
  std::vector<double> unpaired_probabilities;
  // Whole numeric forward plus regenerated beta/marginal run.
  double elapsed_seconds = 0.0;
  // Deprecated diagnostic; production beam outside no longer creates a tape.
  std::size_t tape_nodes = 0;
  std::size_t reverse_records = 0;
  // Retained backend: no saved terms/factors; byte estimate unavailable (0).
  std::size_t reverse_terms = 0;
  std::size_t reverse_factors = 0;
  std::size_t reverse_storage_bytes = 0;
  BeamWorkStats work;
};

struct BeamViterbiResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool reachable = false;
  std::string error;
  std::size_t beam_size = 0;
  double log_weight = -std::numeric_limits<double>::infinity();
  // Convenience projection. It may underflow to zero on long sequences.
  double weight = 0.0;
  double energy = std::numeric_limits<double>::quiet_NaN();
  std::vector<BasePair> generated_pairs;
  std::string scaffold_structure;
  std::string generated_structure;
  std::string union_structure;
  BeamWorkStats work;
};

struct BeamProbabilityResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool denominator_reachable = false;
  bool target_included = false;
  std::string error;
  std::size_t beam_size = 0;
  double log_partition_weight =
      -std::numeric_limits<double>::infinity();
  double log_target_weight = -std::numeric_limits<double>::infinity();
  double log_probability = -std::numeric_limits<double>::infinity();
  double probability = 0.0;
  // This is the free energy of all retained derivations producing the target
  // pair map, not a promise that the grammar has only one such derivation.
  double target_free_energy = std::numeric_limits<double>::quiet_NaN();
  std::vector<BasePair> target_generated_pairs;
};

// BeamSession is intentionally serial. Thermodynamic parameter loading uses
// Vienna's process-global state, so callers must not invoke multiple sessions
// concurrently.
class BeamSession {
 public:
  BeamSession(std::string sequence,
              std::string scaffold_structure,
              BeamConfig config = {});
  BeamSession(std::string sequence,
              std::string scaffold_structure,
              std::size_t beam_size);

  const std::string &sequence() const;
  const std::string &scaffold_structure() const;
  const BeamConfig &config() const;

  // LogSumExp inside on the configured finite beam DAG. Cached per session.
  BeamInsideResult inside();
  BeamOutsideResult outside() const;

  // LogMax prediction with the configured ranking. Outer uses the same
  // LogSumExp pilot guide as inside, with LogMax main priorities.
  BeamViterbiResult viterbi() const;

  // Fail-closed probability on one LogSumExp-ranked retained DAG. If the target
  // was pruned, ok remains true and target_included is false.
  BeamProbabilityResult probability(
      const std::vector<BasePair> &target_pairs) const;
  BeamProbabilityResult probability(
      const std::string &target_structure) const;

 private:
  bool validate_inputs(std::string *error) const;

  std::string sequence_;
  std::string scaffold_structure_;
  BeamConfig config_;
  std::optional<BeamInsideResult> inside_cache_;
  // Immutable numeric guide and forced keys; shared across session calls.
  mutable std::shared_ptr<const beam::PreparedRanking> ranking_cache_;
};

BeamInsideResult beam_inside(const std::string &sequence,
                             const std::string &scaffold_structure,
                             BeamConfig config = {});
BeamOutsideResult beam_outside(const std::string &sequence,
                               const std::string &scaffold_structure,
                               BeamConfig config = {});
BeamViterbiResult beam_viterbi(const std::string &sequence,
                               const std::string &scaffold_structure,
                               BeamConfig config = {});
BeamProbabilityResult beam_probability(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::vector<BasePair> &target_pairs,
    BeamConfig config = {});
BeamProbabilityResult beam_probability(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::string &target_structure,
    BeamConfig config = {});

}  // namespace cparty::scfg2::api

#endif  // CPARTY_SCFG2_API_BEAM_HH

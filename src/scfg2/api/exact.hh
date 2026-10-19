#ifndef CPARTY_SCFG2_API_EXACT_HH
#define CPARTY_SCFG2_API_EXACT_HH

#include "scfg2/api/types.hh"

#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace cparty::scfg2::api {

struct InsideResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool reachable = false;
  std::string error;
  // Compatibility projection: may be zero or infinity for a valid log result.
  double partition_weight = 0.0;
  double log_partition_weight = -std::numeric_limits<double>::infinity();
  double ensemble_free_energy = std::numeric_limits<double>::quiet_NaN();
};

struct OutsideResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool reachable = false;
  std::string error;
  double log_partition_weight =
      -std::numeric_limits<double>::infinity();
  double partition_weight = 0.0;
  double ensemble_free_energy = std::numeric_limits<double>::quiet_NaN();
  std::vector<PairProbability> pair_probabilities;
  // One-based: index zero is unused.
  std::vector<double> unpaired_probabilities;
};

struct ViterbiResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool reachable = false;
  std::string error;
  double weight = 0.0;
  double log_weight = -std::numeric_limits<double>::infinity();
  double energy = std::numeric_limits<double>::quiet_NaN();
  std::vector<BasePair> generated_pairs;
  std::string scaffold_structure;
  std::string generated_structure;
  std::string union_structure;
};

struct EnergyEvaluationResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool reachable = false;
  std::string error;
  double weight = 0.0;
  double log_weight = -std::numeric_limits<double>::infinity();
  double energy = std::numeric_limits<double>::quiet_NaN();
  std::vector<BasePair> generated_pairs;
  std::string scaffold_structure;
  std::string generated_structure;
  std::string union_structure;
};

struct ProbabilityResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  bool reachable = false;
  std::string error;
  double partition_weight = 0.0;
  double target_weight = 0.0;
  double log_partition_weight = -std::numeric_limits<double>::infinity();
  double log_target_weight = -std::numeric_limits<double>::infinity();
  double energy = std::numeric_limits<double>::quiet_NaN();
  double probability = 0.0;
  double log_probability = -std::numeric_limits<double>::infinity();
  std::string scaffold_structure;
  std::string generated_structure;
  std::string union_structure;
};

// ExactSession is intentionally serial. Thermodynamic parameter loading uses
// Vienna's process-global state, so callers must not invoke methods concurrently
// on one or multiple sessions.
class ExactSession {
 public:
  ExactSession(std::string sequence, std::string scaffold_structure,
               WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);

  const std::string &sequence() const;
  const std::string &scaffold_structure() const;

  WmbClosurePolicy wmb_closure_policy() const;

  // The first call computes Z; subsequent calls return the cached value.
  InsideResult inside();
  OutsideResult outside() const;
  ViterbiResult viterbi() const;
  EnergyEvaluationResult evaluate_energy(
      const std::vector<BasePair> &target_pairs) const;
  EnergyEvaluationResult evaluate_energy(
      const std::string &target_structure) const;
  ProbabilityResult probability(const std::vector<BasePair> &target_pairs);
  ProbabilityResult probability(const std::string &target_structure);

 private:
  bool validate_inputs(std::string *error) const;
  ProbabilityResult probability_from_evaluation(
      const EnergyEvaluationResult &evaluated);

  std::string sequence_;
  std::string scaffold_structure_;
  WmbClosurePolicy wmb_closure_policy_;
  std::optional<InsideResult> inside_cache_;
};

InsideResult inside(const std::string &sequence,
                    const std::string &scaffold_structure,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);
OutsideResult outside(const std::string &sequence,
                      const std::string &scaffold_structure,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);
ViterbiResult viterbi(const std::string &sequence,
                      const std::string &scaffold_structure,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);
EnergyEvaluationResult evaluate_energy(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::vector<BasePair> &target_pairs,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);
EnergyEvaluationResult evaluate_energy(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::string &target_structure,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);
ProbabilityResult probability(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::vector<BasePair> &target_pairs,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);
ProbabilityResult probability(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::string &target_structure,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);

}  // namespace cparty::scfg2::api

#endif  // CPARTY_SCFG2_API_EXACT_HH

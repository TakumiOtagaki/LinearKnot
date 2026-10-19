#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/exact.hh"

#include "scfg2/api/detail.hh"
#include "scfg2/engine/outside_marginals.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/replay/w_final_exact_outside.hh"

#include <cmath>
#include <exception>
#include <utility>

namespace cparty::scfg2::api {
using detail::compose_union_structure;
using detail::parse_dot_bracket_pairs;
using detail::render_generated_structure;
using detail::to_base_pairs;
using detail::to_traceback_pairs;
using detail::validate_pair_map;

ExactSession::ExactSession(std::string sequence, std::string scaffold_structure,
                           WmbClosurePolicy policy)
    : sequence_(std::move(sequence)),
      scaffold_structure_(std::move(scaffold_structure)),
      wmb_closure_policy_(policy) {}

WmbClosurePolicy ExactSession::wmb_closure_policy() const {
  return wmb_closure_policy_;
}

const std::string &ExactSession::sequence() const { return sequence_; }

const std::string &ExactSession::scaffold_structure() const {
  return scaffold_structure_;
}

InsideResult ExactSession::inside() {
  if (inside_cache_.has_value()) return *inside_cache_;
  InsideResult out;
  out.wmb_closure_policy = wmb_closure_policy_;
  if (!validate_inputs(&out.error)) {
    inside_cache_ = out;
    return out;
  }
  try {
    std::string sequence = sequence_;
    std::string scaffold = scaffold_structure_;
    weights::RuntimeSnapshot pf(sequence, scaffold);
    replay::WFinalLogWeightContext ctx(pf, wmb_closure_policy_);
    const auto n = static_cast<cand_pos_t>(sequence_.size());
    const auto chart = replay::run_w_final_exact_dp_schedule<semirings::LogSumExpSemiring>(
        engine::build_pkfree_basic_exact_schedule(n, ctx), ctx,
        [](const auto &) { return true; },
        [](const auto &, const auto &, const auto &, const auto &, auto) {});
    out.log_partition_weight = chart.get({grammar::NonTerminal::W, 1, n, -1, -1});
    out.reachable = std::isfinite(out.log_partition_weight);
    if (!out.reachable) {
      out.error = "exact inside returned an unreachable or non-finite log weight";
    } else {
      out.partition_weight = std::exp(out.log_partition_weight);
      out.ensemble_free_energy = pf.partition_free_energy_from_log_weight(out.log_partition_weight);
      out.ok = std::isfinite(out.ensemble_free_energy);
      if (!out.ok) out.error = "exact inside returned a non-finite free energy";
    }
  } catch (const std::exception &e) {
    out.error = e.what();
  }
  inside_cache_ = out;
  return out;
}

OutsideResult ExactSession::outside() const {
  OutsideResult out;
  out.wmb_closure_policy = wmb_closure_policy_;
  if (!validate_inputs(&out.error)) return out;
  try {
    std::string sequence = sequence_;
    std::string scaffold = scaffold_structure_;
    weights::RuntimeSnapshot pf(sequence, scaffold);
    replay::WFinalLogWeightContext ctx(pf, wmb_closure_policy_);
    const auto n = static_cast<cand_pos_t>(sequence_.size());
    const auto schedule =
        engine::build_pkfree_basic_exact_schedule(n, ctx);
    using S = semirings::LogSumExpSemiring;
    const auto chart = replay::run_w_final_exact_dp_schedule<S>(schedule, ctx,
        [](const deductions::Deduction &) { return true; },
        [](const auto &...) {});
    const deductions::ItemKey root{
        grammar::NonTerminal::W, 1, n, -1, -1};
    const auto root_value = chart.get(root);
    out.log_partition_weight = root_value;
    out.reachable = std::isfinite(out.log_partition_weight);
    if (!out.reachable) {
      out.ok = true;
      return out;
    }
    const auto beta = replay::run_w_final_exact_regenerated_outside(schedule, ctx, chart, root);
    const auto marginal = engine::collect_outside_marginals_by_item(chart,
        [&](const auto &item, double alpha) -> std::optional<double> {
          // Exact charts preserve zero entries, while a zero tape value had
          // no node and contributed no public pair row. Keep that contract.
          if (S::is_zero(alpha)) return std::nullopt;
          return std::exp(alpha + beta.get(item) - root_value);
        }, n, ctx);
    if (!marginal.ok) {
      out.error = marginal.error;
      return out;
    }
    out.partition_weight = std::exp(out.log_partition_weight);
    out.ensemble_free_energy =
        pf.partition_free_energy_from_log_weight(out.log_partition_weight);
    out.unpaired_probabilities = marginal.unpaired_probability;
    out.pair_probabilities.reserve(marginal.pairs.size());
    for (const auto &pair : marginal.pairs) {
      out.pair_probabilities.push_back(PairProbability{
          pair.i,
          pair.j,
          pair.round_probability,
          pair.square_probability,
          pair.union_probability});
    }
    out.ok = std::isfinite(out.ensemble_free_energy);
    if (!out.ok) out.error = "exact outside returned non-finite free energy";
  } catch (const std::exception &e) {
    out.error = e.what();
  }
  return out;
}

ViterbiResult ExactSession::viterbi() const {
  ViterbiResult out;
  out.wmb_closure_policy = wmb_closure_policy_;
  if (!validate_inputs(&out.error)) return out;
  try {
    std::string sequence = sequence_;
    std::string scaffold = scaffold_structure_;
    weights::RuntimeSnapshot pf(sequence, scaffold);
    const auto result = replay::run_w_final_exact_viterbi_traceback<true>(pf, wmb_closure_policy_);
    out.log_weight = result.run.chart.get(result.root);
    out.weight = std::exp(out.log_weight);
    out.reachable = std::isfinite(out.log_weight);
    if (!out.reachable) {
      out.error = "exact Viterbi returned an unreachable or non-finite log weight";
      return out;
    }
    out.energy = pf.partition_free_energy_from_log_weight(out.log_weight);
    out.generated_pairs = to_base_pairs(result.generated_pairs);
    out.scaffold_structure = scaffold_structure_;
    replay::WFinalLogWeightContext ctx(pf, wmb_closure_policy_);
    out.generated_structure = render_generated_structure(
        static_cast<int>(sequence_.size()), result.pairs, ctx);
    out.union_structure = compose_union_structure(
        out.scaffold_structure, out.generated_structure);
    out.ok = std::isfinite(out.energy);
    if (!out.ok) out.error = "exact Viterbi returned a non-finite energy";
  } catch (const std::exception &e) {
    out.error = e.what();
  }
  return out;
}

EnergyEvaluationResult ExactSession::evaluate_energy(
    const std::vector<BasePair> &target_pairs) const {
  EnergyEvaluationResult out;
  out.wmb_closure_policy = wmb_closure_policy_;
  if (!validate_inputs(&out.error) ||
      !validate_pair_map(target_pairs, static_cast<int>(sequence_.size()),
                         &out.error)) {
    return out;
  }
  try {
    std::string sequence = sequence_;
    std::string scaffold = scaffold_structure_;
    weights::RuntimeSnapshot pf(sequence, scaffold);
    const auto result =
        replay::run_w_final_exact_output_constrained_viterbi_traceback<true>(
            pf, to_traceback_pairs(target_pairs), wmb_closure_policy_);
    out.log_weight = result.run.chart.get(result.root);
    out.weight = std::exp(out.log_weight);
    out.scaffold_structure = scaffold_structure_;
    replay::WFinalLogWeightContext ctx(pf, wmb_closure_policy_);
    out.generated_structure = render_generated_structure(
        static_cast<int>(sequence_.size()), result.pairs, ctx);
    out.generated_pairs = to_base_pairs(result.generated_pairs);
    out.union_structure = compose_union_structure(
        out.scaffold_structure, out.generated_structure);
    out.ok = true;
    out.reachable = std::isfinite(out.log_weight);
    if (out.reachable) {
      out.energy = pf.partition_free_energy_from_log_weight(out.log_weight);
      if (!std::isfinite(out.energy)) {
        out.ok = false;
        out.reachable = false;
        out.error = "fixed-output evaluator returned a non-finite energy";
      }
    }
  } catch (const std::exception &e) {
    out.error = e.what();
  }
  return out;
}

EnergyEvaluationResult ExactSession::evaluate_energy(
    const std::string &target_structure) const {
  EnergyEvaluationResult out;
  out.wmb_closure_policy = wmb_closure_policy_;
  if (target_structure.size() != sequence_.size()) {
    out.error = "target structure length differs from sequence length";
    return out;
  }
  std::vector<BasePair> target_pairs;
  if (!parse_dot_bracket_pairs(target_structure, &target_pairs, &out.error)) {
    return out;
  }
  return evaluate_energy(target_pairs);
}

ProbabilityResult ExactSession::probability(
    const std::vector<BasePair> &target_pairs) {
  return probability_from_evaluation(evaluate_energy(target_pairs));
}

ProbabilityResult ExactSession::probability(
    const std::string &target_structure) {
  return probability_from_evaluation(evaluate_energy(target_structure));
}

bool ExactSession::validate_inputs(std::string *error) const {
  return detail::validate_session_inputs(sequence_, scaffold_structure_, error);
}

ProbabilityResult ExactSession::probability_from_evaluation(
    const EnergyEvaluationResult &evaluated) {
  ProbabilityResult out;
  out.wmb_closure_policy = wmb_closure_policy_;
  out.reachable = evaluated.reachable;
  out.target_weight = evaluated.weight;
  out.log_target_weight = evaluated.log_weight;
  out.energy = evaluated.energy;
  out.scaffold_structure = evaluated.scaffold_structure;
  out.generated_structure = evaluated.generated_structure;
  out.union_structure = evaluated.union_structure;
  if (!evaluated.ok) {
    out.error = evaluated.error;
    return out;
  }
  const InsideResult partition = inside();
  out.partition_weight = partition.partition_weight;
  out.log_partition_weight = partition.log_partition_weight;
  if (!partition.ok) {
    out.error = partition.error;
    return out;
  }
  if (!evaluated.reachable) {
    out.ok = true;
    return out;
  }
  out.log_probability =
      out.log_target_weight - out.log_partition_weight;
  out.probability = std::exp(out.log_probability);
  if (!std::isfinite(out.log_probability) || !std::isfinite(out.probability) ||
      out.probability > 1.0 + 1e-12) {
    out.error = "fixed-output probability is non-finite or exceeds one";
    return out;
  }
  out.ok = true;
  return out;
}

InsideResult inside(const std::string &sequence,
                    const std::string &scaffold_structure, WmbClosurePolicy policy) {
  return ExactSession(sequence, scaffold_structure, policy).inside();
}

OutsideResult outside(const std::string &sequence,
                      const std::string &scaffold_structure, WmbClosurePolicy policy) {
  return ExactSession(sequence, scaffold_structure, policy).outside();
}

ViterbiResult viterbi(const std::string &sequence,
                      const std::string &scaffold_structure, WmbClosurePolicy policy) {
  return ExactSession(sequence, scaffold_structure, policy).viterbi();
}

EnergyEvaluationResult evaluate_energy(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::vector<BasePair> &target_pairs, WmbClosurePolicy policy) {
  return ExactSession(sequence, scaffold_structure, policy).evaluate_energy(target_pairs);
}

EnergyEvaluationResult evaluate_energy(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::string &target_structure, WmbClosurePolicy policy) {
  return ExactSession(sequence, scaffold_structure, policy)
      .evaluate_energy(target_structure);
}

ProbabilityResult probability(const std::string &sequence,
                              const std::string &scaffold_structure,
                              const std::vector<BasePair> &target_pairs, WmbClosurePolicy policy) {
  return ExactSession(sequence, scaffold_structure, policy).probability(target_pairs);
}

ProbabilityResult probability(const std::string &sequence,
                              const std::string &scaffold_structure,
                              const std::string &target_structure, WmbClosurePolicy policy) {
  return ExactSession(sequence, scaffold_structure, policy).probability(target_structure);
}

}  // namespace cparty::scfg2::api

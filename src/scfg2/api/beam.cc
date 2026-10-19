#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/beam.hh"

#include "scfg2/api/detail.hh"
#include "scfg2/api/beam_detail.hh"
#include "scfg2/beam/retained_outside.hh"
#include "scfg2/engine/outside_marginals.hh"
#include "scfg2/engine/chosen_pair_output.hh"
#include "scfg2/engine/compact_choice_map.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <utility>

namespace cparty::scfg2::api {
namespace {

const beam::PreparedRanking &prepared_ranking(
    std::shared_ptr<const beam::PreparedRanking> &cache,
    replay::WFinalLogWeightContext &ctx, cand_pos_t n, const BeamConfig &config) {
  if (!cache)
    cache = std::make_shared<const beam::PreparedRanking>(
        detail::prepare_beam_ranking(ctx, n, config));
  return *cache;
}

deductions::ItemKey root_item(cand_pos_t n) {
  return deductions::ItemKey{
      grammar::NonTerminal::W, 1, n, -1, -1};
}

double projected_weight(double log_weight) {
  return std::isfinite(log_weight) ? std::exp(log_weight) : 0.0;
}

template <class Run>
BeamWorkStats summarize_work(const Run &run) {
  BeamWorkStats out;
  out.vp_candidates = run.vp_candidates;
  out.vp_selected = run.vp_selected;
  out.wmbp_candidates = run.wmbp_candidates;
  out.wmbp_selected = run.wmbp_selected;
  out.case3_candidates =
      run.case3_producer_candidates + run.case3_tail_candidates +
      run.case3_left_candidates;
  out.case3_selected =
      run.case3_producer_selected + run.case3_tail_selected +
      run.case3_left_selected;
  out.binary_join_probes =
      run.vm_binary_join_probes + run.wm_binary_join_probes +
      run.wi_binary_join_probes + run.wip_binary_join_probes +
      run.wmbw_binary_join_probes + run.vpr_binary_join_probes +
      run.w_binary_join_probes;
  out.extension_probes =
      run.wm_extension_probes + run.wi_extension_probes +
      run.wip_extension_probes + run.w_extension_probes;
  out.v_parent_probes =
      run.v_hairpin_parent_probes + run.v_internal_parent_probes +
      run.v_vm_parent_probes;
  out.vp_loop_parent_probes = run.vp_loop_parent_probes;
  out.case3_truncated_endpoints = run.case3_truncated_endpoints;
  out.families.reserve(replay::w_final_online_nonterminal_count);
  for (std::size_t index = 0;
       index < replay::w_final_online_nonterminal_count;
       ++index) {
    const auto nonterminal =
        static_cast<grammar::NonTerminal>(index);
    const auto visits = run.item_visits[index];
    const auto positive = run.positive_items[index];
    const auto candidates = run.frontier_candidates[index];
    const auto selected = run.frontier_selected[index];
    out.total_item_visits += visits;
    out.total_positive_items += positive;
    out.total_frontier_candidates += candidates;
    out.total_frontier_selected += selected;
    out.families.push_back(BeamFamilyWork{
        std::string(grammar::nonterminal_name(nonterminal)),
        visits,
        positive,
        candidates,
        selected,
    });
  }
  return out;
}

}  // namespace

BeamSession::BeamSession(std::string sequence,
                         std::string scaffold_structure,
                         BeamConfig config)
    : sequence_(std::move(sequence)),
      scaffold_structure_(std::move(scaffold_structure)),
      config_(config) {}

BeamSession::BeamSession(std::string sequence,
                         std::string scaffold_structure,
                         std::size_t beam_size)
    : BeamSession(std::move(sequence),
                  std::move(scaffold_structure),
                  BeamConfig{beam_size, 30, 30}) {}

const std::string &BeamSession::sequence() const { return sequence_; }

const std::string &BeamSession::scaffold_structure() const {
  return scaffold_structure_;
}

const BeamConfig &BeamSession::config() const { return config_; }

BeamInsideResult BeamSession::inside() {
  if (inside_cache_.has_value()) return *inside_cache_;
  const auto started = std::chrono::steady_clock::now();
  BeamInsideResult out;
  out.wmb_closure_policy = config_.wmb_closure_policy;
  out.beam_size = config_.beam_size;
  if (!validate_inputs(&out.error)) {
    inside_cache_ = out;
    return out;
  }
  try {
    std::string sequence = sequence_;
    std::string scaffold = scaffold_structure_;
    weights::RuntimeSnapshot pf(sequence, scaffold);
    replay::WFinalLogWeightContext ctx(pf, config_.wmb_closure_policy);
    const auto n = static_cast<cand_pos_t>(sequence_.size());
    const auto &prepared = prepared_ranking(ranking_cache_, ctx, n, config_);
    const auto run = replay::run_w_final_vp_wmbp_online<semirings::LogSumExpSemiring>(
        ctx, n, prepared.widths, prepared.ranking);
    out.work = summarize_work(run);
    out.log_partition_weight = run.chart.get(root_item(n));
    out.reachable = std::isfinite(out.log_partition_weight);
    out.ok = true;
    if (out.reachable) {
      out.partition_weight = projected_weight(out.log_partition_weight);
      out.ensemble_free_energy =
          pf.partition_free_energy_from_log_weight(out.log_partition_weight);
      if (!std::isfinite(out.ensemble_free_energy)) {
        out.ok = false;
        out.error = "beam inside returned a non-finite free energy";
      }
    }
  } catch (const std::exception &e) {
    out.ok = false;
    out.reachable = false;
    out.error = e.what();
  }
  out.elapsed_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  inside_cache_ = out;
  return out;
}

BeamOutsideResult BeamSession::outside() const {
  const auto started = std::chrono::steady_clock::now();
  BeamOutsideResult out;
  out.wmb_closure_policy = config_.wmb_closure_policy;
  const auto finish_timing = [&] {
    out.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
  };
  out.beam_size = config_.beam_size;
  if (!validate_inputs(&out.error)) {
    finish_timing();
    return out;
  }
  try {
    std::string sequence = sequence_;
    std::string scaffold = scaffold_structure_;
    weights::RuntimeSnapshot pf(sequence, scaffold);
    replay::WFinalLogWeightContext ctx(pf, config_.wmb_closure_policy);
    const auto n = static_cast<cand_pos_t>(sequence_.size());
    using S = semirings::LogSumExpSemiring;
    const auto &prepared = prepared_ranking(ranking_cache_, ctx, n, config_);
    const auto &widths = prepared.widths;
    beam::RetainedRuleAdmission<S> rules;
    const auto run = replay::run_w_final_vp_wmbp_online_observed<S>(ctx, n, widths, rules, prepared.ranking);
    rules.finish(run.chart);
    beam::RetainedOutsideProgram<S, replay::WFinalLogWeightContext> program(
        n, run.chart, ctx, rules, widths, prepared.ranking);
    out.reverse_records = run.chart.values().size() + program.private_cell_count();
    out.work = summarize_work(run);
    const auto root = root_item(n);
    const auto root_value = run.chart.get(root);
    out.log_partition_weight = root_value;
    out.reachable = std::isfinite(out.log_partition_weight);
    if (!out.reachable) {
      out.ok = true;
      finish_timing();
      return out;
    }
    const auto reverse = program.outside(root);
    const auto marginal = engine::collect_outside_marginals_by_item(run.chart,
        [&](const auto &item, double alpha) {
          return std::exp(alpha + reverse.beta.get(item) - root_value);
        }, n, ctx);
    if (!marginal.ok) {
      out.error = marginal.error;
      finish_timing();
      return out;
    }
    out.partition_weight = projected_weight(out.log_partition_weight);
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
    if (!out.ok) out.error = "beam outside returned non-finite free energy";
  } catch (const std::exception &e) {
    out.error = e.what();
  }
  finish_timing();
  return out;
}

BeamViterbiResult BeamSession::viterbi() const {
  BeamViterbiResult out;
  out.wmb_closure_policy = config_.wmb_closure_policy;
  out.beam_size = config_.beam_size;
  out.scaffold_structure = scaffold_structure_;
  if (!validate_inputs(&out.error)) return out;
  try {
    std::string sequence = sequence_;
    std::string scaffold = scaffold_structure_;
    weights::RuntimeSnapshot pf(sequence, scaffold);
    replay::WFinalLogWeightContext ctx(pf, config_.wmb_closure_policy);
    const auto n = static_cast<cand_pos_t>(sequence_.size());
    const auto &prepared = prepared_ranking(ranking_cache_, ctx, n, config_);
    auto online = replay::run_w_final_vp_wmbp_online_viterbi_stored<
        semirings::LogMaxSemiring, engine::CompactChoiceMap<double>>(
            ctx, n, prepared.widths, prepared.ranking);
    out.work = summarize_work(online.dp);
    const auto root = root_item(n);
    out.log_weight = online.dp.chart.get(root);
    out.reachable = std::isfinite(out.log_weight);
    out.ok = true;
    if (!out.reachable) return out;

    auto result = engine::collect_chosen_pair_output(root, online.choices, ctx);
    out.weight = projected_weight(out.log_weight);
    out.energy = pf.partition_free_energy_from_log_weight(out.log_weight);
    out.generated_pairs = detail::to_base_pairs_linear(result.generated_pairs, n);
    out.generated_structure = detail::render_generated_structure(
        static_cast<int>(sequence_.size()), result.pairs, ctx);
    out.union_structure = detail::compose_union_structure(
        out.scaffold_structure, out.generated_structure);
    if (!std::isfinite(out.energy)) {
      out.ok = false;
      out.reachable = false;
      out.error = "beam Viterbi returned a non-finite energy";
    }
  } catch (const std::exception &e) {
    out.ok = false;
    out.reachable = false;
    out.error = e.what();
  }
  return out;
}

BeamProbabilityResult BeamSession::probability(
    const std::vector<BasePair> &target_pairs) const {
  BeamProbabilityResult out;
  out.wmb_closure_policy = config_.wmb_closure_policy;
  out.beam_size = config_.beam_size;
  if (!validate_inputs(&out.error) ||
      !detail::validate_pair_map(
          target_pairs, static_cast<int>(sequence_.size()), &out.error)) {
    return out;
  }
  try {
    std::string sequence = sequence_;
    std::string scaffold = scaffold_structure_;
    weights::RuntimeSnapshot pf(sequence, scaffold);
    replay::WFinalLogWeightContext ctx(pf, config_.wmb_closure_policy);
    const auto n = static_cast<cand_pos_t>(sequence_.size());
    const auto traceback_pairs = detail::to_traceback_pairs(target_pairs);
    out.target_generated_pairs = detail::to_base_pairs(
        engine::generated_traceback_pairs(traceback_pairs, ctx));
    const auto &prepared = prepared_ranking(ranking_cache_, ctx, n, config_);
    const auto result = replay::run_w_final_vp_wmbp_online_target_probability(
        ctx, n, prepared.widths, traceback_pairs, prepared.ranking);
    out.denominator_reachable = result.denominator_reachable;
    out.target_included = result.target_included;
    out.log_partition_weight = result.log_partition;
    out.log_target_weight = result.log_target_weight;
    out.log_probability = result.log_probability;
    out.probability = result.probability;

    if (!out.denominator_reachable) {
      out.ok = true;
      return out;
    }
    if (!out.target_included) {
      out.ok = true;
      return out;
    }
    if (!result.ok) {
      out.error =
          "beam target probability is non-finite or exceeds one";
      return out;
    }
    out.target_free_energy =
        pf.partition_free_energy_from_log_weight(out.log_target_weight);
    out.ok = std::isfinite(out.target_free_energy);
    if (!out.ok) {
      out.error = "beam target mass returned a non-finite free energy";
    }
  } catch (const std::exception &e) {
    out.error = e.what();
  }
  return out;
}

BeamProbabilityResult BeamSession::probability(
    const std::string &target_structure) const {
  BeamProbabilityResult out;
  out.wmb_closure_policy = config_.wmb_closure_policy;
  out.beam_size = config_.beam_size;
  if (target_structure.size() != sequence_.size()) {
    out.error = "target structure length differs from sequence length";
    return out;
  }
  std::vector<BasePair> target_pairs;
  if (!detail::parse_dot_bracket_pairs(
          target_structure, &target_pairs, &out.error)) {
    return out;
  }
  return probability(target_pairs);
}

bool BeamSession::validate_inputs(std::string *error) const {
  if (!detail::validate_session_inputs(
          sequence_, scaffold_structure_, error)) {
    return false;
  }
  return detail::validate_beam_config(config_, error);
}

BeamInsideResult beam_inside(const std::string &sequence,
                             const std::string &scaffold_structure,
                             BeamConfig config) {
  return BeamSession(sequence, scaffold_structure, config).inside();
}

BeamOutsideResult beam_outside(const std::string &sequence,
                               const std::string &scaffold_structure,
                               BeamConfig config) {
  return BeamSession(sequence, scaffold_structure, config).outside();
}

BeamViterbiResult beam_viterbi(const std::string &sequence,
                               const std::string &scaffold_structure,
                               BeamConfig config) {
  return BeamSession(sequence, scaffold_structure, config).viterbi();
}

BeamProbabilityResult beam_probability(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::vector<BasePair> &target_pairs,
    BeamConfig config) {
  return BeamSession(sequence, scaffold_structure, config)
      .probability(target_pairs);
}

BeamProbabilityResult beam_probability(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::string &target_structure,
    BeamConfig config) {
  return BeamSession(sequence, scaffold_structure, config)
      .probability(target_structure);
}

}  // namespace cparty::scfg2::api

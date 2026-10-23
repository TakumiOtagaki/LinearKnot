#ifndef CPARTY_SCFG2_PROFILES_BEAM_PROFILES_HH
#define CPARTY_SCFG2_PROFILES_BEAM_PROFILES_HH

#include "scfg2/beam/run_term_observer.hh"
#include "scfg2/beam/prepare_ranking.hh"
#include "scfg2/beam/retained_outside.hh"
#include "scfg2/profiles/loop_context.hh"
#include "scfg2/profiles/profile_assembly.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/semirings/logsumexp_tape.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <cmath>
#include <string>
#include <vector>

namespace cparty::scfg2::profiles {

namespace detail {

// Rules the beam folds into start-prefix aggregates or the VPR right-run
// prefix.  Their posteriors come from RunTermObserver, not from the
// aggregate candidate.
inline bool beam_aggregated_run_rule(grammar::RuleId rule) {
  switch (rule) {
    case grammar::RuleId::WM_START_V:
    case grammar::RuleId::WM_START_WMB:
    case grammar::RuleId::WIP_BASEPAIR_V:
    case grammar::RuleId::WIP_BASEPAIR_WMB:
    case grammar::RuleId::VM_SPLIT_WMp_BASE:
    case grammar::RuleId::VPL_SPLIT_VP:
    case grammar::RuleId::VPR_SPLIT_VP_BASEPAIR:
      return true;
    default:
      return false;
  }
}

inline LoopContext aggregated_run_context(grammar::RuleId rule) {
  if (!beam_aggregated_run_rule(rule))
    throw std::logic_error("rule has no aggregated unpaired run");
  return fixed_rule_context(rule);
}

// Chain family id for the VPR right-run prefix, distinct from nonterminals.
inline constexpr int kVprRunChainFamily = 1 << 16;
inline constexpr int kLeftRunChainFamily = 1 << 17;

template <class Value>
class RunTermRecorder : public beam::RunTermObserver<Value> {
 public:
  struct Term {
    grammar::RuleId rule;
    cand_pos_t lo;
    cand_pos_t hi;
    Value value;
  };
  struct VprCell {
    cand_pos_t i;
    cand_pos_t t;
    Value total;
    bool carried;
  };

  void on_run_term(grammar::RuleId rule, cand_pos_t lo, cand_pos_t hi,
                   const Value &term) override {
    terms.push_back(Term{rule, lo, hi, semirings::checkpoint_value(term)});
  }
  void on_vpr_run_cell(cand_pos_t i, cand_pos_t t, const Value &total,
                       bool carried) override {
    vpr_cells.push_back(VprCell{i, t, semirings::checkpoint_value(total), carried});
  }
  void on_left_run_unit(grammar::RuleId rule, cand_pos_t j, cand_pos_t base,
                        const Value &unit, const Value *alias,
                        bool linked) override {
    left_units.push_back(LeftUnit{rule, j, base, semirings::checkpoint_value(unit),
                                  alias ? semirings::checkpoint_value(*alias) : Value{}, alias != nullptr,
                                  linked});
  }

  struct LeftUnit {
    grammar::RuleId rule;
    cand_pos_t j;
    cand_pos_t base;
    Value unit;
    Value alias;
    bool has_alias;
    bool linked;
  };
  std::vector<Term> terms;
  std::vector<VprCell> vpr_cells;
  std::vector<LeftUnit> left_units;
};

}  // namespace detail

// Shared assembly for production analytic outside and the taped oracle.
template <class Semiring, class Graph>
ProfileResult beam_profiles_with_graph(const std::string &sequence,
                                   const std::string &scaffold,
                                   const replay::WFinalVpWmbpOnlineWidths &widths,
                                   const std::vector<cand_pos_t> &window_lengths,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure,
    beam::RankingMode ranking_mode = beam::RankingMode::Baseline,
    beam::PilotExponent exponent = beam::PilotExponent::TwoThirds,
    bool by_loop = false) {
  using Value = typename Semiring::value_type;
  std::string seq = sequence;
  std::string g = scaffold;
  weights::RuntimeSnapshot pf(seq, g);
  replay::WFinalLogWeightContext ctx(pf, policy);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto prepared = beam::prepare_ranking(ctx, n, widths, widths.v, ranking_mode, exponent);
  const auto is_unpaired = [&](cand_pos_t pos) { return ctx.is_unpaired(pos); };

  struct Recorded {
    deductions::Deduction deduction;
    UnpairedEmission emission;
    Value contribution;
    weights::LogWeight local{weights::log_zero};
  };
  std::vector<Recorded> recorded;
  detail::RunTermRecorder<Value> recorder;

  Graph graph;
  const auto run = [&] {
    typename Semiring::Scope scope(graph);
    beam::ScopedRunTermObserver<Value> observe(&recorder);
    return replay::run_w_final_vp_wmbp_online_observed<Semiring>(
        ctx, n, prepared.widths,
        [&](const deductions::ItemKey &,
            const deductions::Deduction &deduction,
            Value contribution,
            const auto &,
            auto local) {
          if (detail::beam_aggregated_run_rule(deduction.rule)) return;
          const auto emission = unpaired_emission(deduction, is_unpaired);
          if (emission.count == 0) return;
          recorded.push_back(Recorded{deduction, emission,
                                     semirings::checkpoint_value(contribution), local});
        }, prepared.ranking);
  }();
  const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
  const Value root_value = run.chart.get(root);
  if (!std::isfinite(root_value.log_value)) {
    ProfileResult out;
    out.error = "beam retained no structure";
    return out;
  }
  const auto adjoints = [&] {
    if constexpr (std::is_same_v<Graph, semirings::LogSumExpTape>) {
      return graph.reverse(root_value.node);
    } else {
      return graph.posteriors(root_value.node);
    }
  }();
  const auto adjoint = [&](std::size_t node) {
    return node < adjoints.size() ? adjoints[node] : 0.0;
  };

  ProfileAccumulator acc(n, by_loop);
  for (const auto &entry : recorded) {
    const double probability = adjoint(entry.contribution.node);
    if (!entry.emission.chain_step) {
      acc.add_whole_runs(entry.emission, probability);
      continue;
    }
    const auto &parent = entry.deduction.parent;
    auto &step = acc.chain_step(static_cast<int>(parent.nonterminal), parent.i,
                                parent.j, entry.emission.runs[0].context);
    step.log_entry = entry.contribution.log_value;
    step.log_step = weights::log_weight(entry.local);
    step.has_prev = entry.deduction.child_count > 0;
    step.emit_probability = probability;
    const Value state = run.chart.get(parent);
    step.log_state = state.log_value;
    step.state_adjoint = adjoint(state.node);
  }
  for (const auto &term : recorder.terms) {
    const double probability = adjoint(term.value.node);
    if (probability == 0.0 || term.lo > term.hi) continue;
    UnpairedEmission emission;
    emission.runs[0] = UnpairedRun{term.lo, term.hi,
                                   detail::aggregated_run_context(term.rule)};
    emission.count = 1;
    acc.add_whole_runs(emission, probability);
  }
  const double log_cp1 = ctx.vpr_split_vp_basepair_weight(1);
  for (const auto &cell : recorder.vpr_cells) {
    auto &step = acc.chain_step(detail::kVprRunChainFamily, cell.i, cell.t,
                                fixed_rule_context(grammar::RuleId::VPR_SPLIT_VP_BASEPAIR));
    const double state_adjoint = adjoint(cell.total.node);
    step.log_entry = cell.total.log_value;
    step.log_step = log_cp1;
    step.has_prev = cell.carried;
    step.log_state = cell.total.log_value;
    step.state_adjoint = state_adjoint;
    step.emit_probability = state_adjoint;
  }
  // KeepAll start-prefix columns: leftward chains keyed by (family, j).
  for (const auto &unit : recorder.left_units) {
    double log_step = 0.0;
    switch (unit.rule) {
      case grammar::RuleId::WM_START_V:
        log_step = ctx.extend_unpaired_weight();
        break;
      case grammar::RuleId::VM_SPLIT_WMp_BASE:
        log_step = ctx.vm_base_unpaired_step_weight();
        break;
      case grammar::RuleId::WIP_BASEPAIR_V:
        log_step = ctx.wip_extend_unpaired_weight();
        break;
      default:
        log_step = ctx.vpl_split_vp_weight(1);
        break;
    }
    double unit_adjoint = adjoint(unit.unit.node);
    if (unit.has_alias) unit_adjoint += adjoint(unit.alias.node);
    auto &step = acc.chain_step(
        detail::kLeftRunChainFamily + static_cast<int>(unit.rule), unit.j,
        unit.base, detail::aggregated_run_context(unit.rule), true);
    step.log_entry = unit.unit.log_value;
    step.log_step = log_step;
    step.has_prev = unit.linked;
    step.log_state = unit.unit.log_value;
    step.state_adjoint = unit_adjoint;
    step.emit_probability = unit_adjoint;
  }
  auto out = acc.finish(root_value.log_value, window_lengths);
  if constexpr (std::is_same_v<Graph, semirings::FactorizedInsideGraph>) {
    out.reverse_records = graph.states().size();
    out.reverse_terms = graph.terms().size();
    out.reverse_factors = graph.factors().size();
    out.reverse_storage_bytes = graph.storage_bytes();
  }
  out.ensemble_free_energy =
      pf.partition_free_energy_from_log_weight(root_value.log_value);
  out.rt_kcal_mol = pf.partition_free_energy_from_log_weight(-1.0) -
                    pf.partition_free_energy_from_log_weight(0.0);
  return out;
}

inline ProfileResult beam_profiles(const std::string &sequence,
                                   const std::string &scaffold,
                                   const replay::WFinalVpWmbpOnlineWidths &widths,
                                   const std::vector<cand_pos_t> &window_lengths,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure,
    beam::RankingMode ranking_mode = beam::RankingMode::Baseline,
    beam::PilotExponent exponent = beam::PilotExponent::TwoThirds,
    bool by_loop = false) {
  using S = semirings::LogSumExpSemiring;
  using NT = grammar::NonTerminal;
  std::string seq = sequence, g = scaffold;
  weights::RuntimeSnapshot pf(seq, g);
  replay::WFinalLogWeightContext ctx(pf, policy);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto prepared = beam::prepare_ranking(ctx, n, widths, widths.v, ranking_mode, exponent);
  beam::RetainedRuleAdmission<S> rules;
  const auto run = replay::run_w_final_vp_wmbp_online_observed<S>(ctx, n, prepared.widths, rules, prepared.ranking);
  rules.finish(run.chart);
  const deductions::ItemKey root{NT::W, 1, n, -1, -1};
  const double log_z = run.chart.get(root);
  if (!std::isfinite(log_z)) {
    ProfileResult out;
    out.error = "beam retained no structure";
    return out;
  }
  beam::RetainedOutsideProgram<S, replay::WFinalLogWeightContext> program(
      n, run.chart, ctx, rules, prepared.widths, prepared.ranking);
  engine::AllLocalWeightModel<replay::WFinalLogWeightContext> weight(&ctx);
  ProfileAccumulator acc(n, window_lengths, by_loop);
  const auto reverse = program.outside(root, [&](const auto &d, double beta) {
    const auto emission = unpaired_emission(d, [&](cand_pos_t pos) { return ctx.is_unpaired(pos); });
    if (emission.count == 0) return;
    const auto local = weight.local_weight_for(d);
    double term = S::from_weight(local);
    for (std::size_t c = 0; c < d.child_count; ++c)
      term = S::times(term, run.chart.get(engine::RuntimeDeductionChildItem{}(d, c)));
    const double probability = std::exp(beta + term - log_z);
    if (!emission.chain_step) {
      acc.add_whole_runs(emission, probability);
      return;
    }
    const auto &parent = d.parent;
    auto &step = acc.chain_step(static_cast<int>(parent.nonterminal), parent.i,
                                parent.j, emission.runs[0].context);
    step.log_entry = term;
    step.log_step = weights::log_weight(local);
    step.has_prev = d.child_count > 0;
    step.emit_probability = probability;
    step.log_state = run.chart.get(parent);
    step.state_adjoint = std::exp(beta + step.log_state - log_z);
  });
  const double log_cp = ctx.vpr_split_vp_basepair_weight(1);
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (const auto &[i, cell] : program.vpr_cells(j)) {
      const auto found = reverse.vpr_beta[j].find(i);
      if (found == reverse.vpr_beta[j].end()) continue;
      auto &step = acc.chain_step(detail::kVprRunChainFamily, i, j,
                                  fixed_rule_context(grammar::RuleId::VPR_SPLIT_VP_BASEPAIR));
      step.log_entry = cell.alpha;
      step.log_step = log_cp;
      step.has_prev = cell.carried;
      step.log_state = cell.alpha;
      step.state_adjoint = std::exp(found->second + cell.alpha - log_z);
      step.emit_probability = step.state_adjoint;
    }
  }
  auto out = acc.finish(log_z, window_lengths);
  out.reverse_records = run.chart.values().size() + program.private_cell_count();
  out.ensemble_free_energy = pf.partition_free_energy_from_log_weight(log_z);
  out.rt_kcal_mol = pf.partition_free_energy_from_log_weight(-1.0) -
                    pf.partition_free_energy_from_log_weight(0.0);
  return out;
}

// Keep the original arithmetic graph available for differential validation.
inline ProfileResult beam_profiles_taped(const std::string &sequence,
                                         const std::string &scaffold,
                                         const replay::WFinalVpWmbpOnlineWidths &widths,
                                         const std::vector<cand_pos_t> &window_lengths,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure,
    beam::RankingMode ranking_mode = beam::RankingMode::Baseline,
    beam::PilotExponent exponent = beam::PilotExponent::TwoThirds,
    bool by_loop = false) {
  return beam_profiles_with_graph<semirings::TapeLogSumExpSemiring,
                                 semirings::LogSumExpTape>(
      sequence, scaffold, widths, window_lengths, policy, ranking_mode, exponent, by_loop);
}

}  // namespace cparty::scfg2::profiles

#endif  // CPARTY_SCFG2_PROFILES_BEAM_PROFILES_HH

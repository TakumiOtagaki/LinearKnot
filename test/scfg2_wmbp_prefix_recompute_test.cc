#include "scfg2/beam/wmbp_be_vp_prefix_value.hh"
#include "scfg2/beam/wmbp_be_vp_prefix_selection.hh"
#include "scfg2/engine/execute.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {
using namespace cparty::scfg2;

bool close(double a, double b) {
  return a == b || std::fabs(a - b) <=
      1e-11 * std::max({1.0, std::fabs(a), std::fabs(b)});
}

deductions::Deduction source(grammar::RuleId rule, bool empty_wi = false) {
  using grammar::NonTerminal;
  deductions::Deduction d;
  d.parent = {NonTerminal::WMBP, 1, 20, -1, -1};
  d.rule = rule;
  d.child_count = empty_wi ? 2 : 3;
  d.children[0].item = {NonTerminal::BE, 2, 18, 9, 15};
  auto left_nt = NonTerminal::WI;
  if (rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP) left_nt = NonTerminal::WMBP;
  if (rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP) left_nt = NonTerminal::WMBW;
  if (!empty_wi) d.children[1].item = {left_nt, 1, 8, -1, -1};
  d.children[d.child_count - 1].item = {NonTerminal::VP, 9, 20, -1, -1};
  return d;
}

template <class Semiring>
bool check_values() {
  using grammar::RuleId;
  for (const auto &d : {source(RuleId::WMBP_SPLIT_BE_WMBP_VP),
                        source(RuleId::WMBP_SPLIT_BE_WMBW_VP),
                        source(RuleId::WMBP_SPLIT_BE_WI_VP),
                        source(RuleId::WMBP_SPLIT_BE_WI_VP, true)}) {
    for (const double local : {0.0, 7.0}) {
      engine::Chart<typename Semiring::value_type> chart(Semiring::zero());
      for (std::size_t c = 0; c < d.child_count; ++c) {
        const double seed = 2.0 + c;
        chart.set(d.children[c].item, Semiring::plus(
            Semiring::from_weight(seed), Semiring::from_weight(2 * seed)));
      }
      const auto prefix = beam::evaluate_wmbp_be_vp_prefix<Semiring>(d, chart, local);
      if (!prefix) return false;
      const auto right = beam::wmbp_be_vp_right_child(d);
      const auto actual = Semiring::times(*prefix, chart.get(right));
      const auto expected = engine::evaluate_deduction<Semiring>(
          d, chart, [&](const auto &) { return local; });
      if (!close(Semiring::priority(actual), Semiring::priority(expected)) ||
          !close(Semiring::traceback_priority(actual),
                 Semiring::traceback_priority(expected))) return false;
      const double before = Semiring::priority(*prefix);
      chart.set(right, Semiring::zero());
      const auto after = beam::evaluate_wmbp_be_vp_prefix<Semiring>(d, chart, local);
      if (!close(before, Semiring::priority(*after))) return false;
    }
  }
  deductions::Deduction direct;
  direct.rule = RuleId::WMBP_DIRECT_VP;
  engine::Chart<typename Semiring::value_type> chart(Semiring::zero());
  return !beam::evaluate_wmbp_be_vp_prefix<Semiring>(direct, chart, 1.0);
}

bool check_transient_graph_storage() {
  using S = semirings::FactorizedLogSumExpSemiring;
  semirings::FactorizedInsideGraph graph;
  S::Scope scope(graph);
  const auto d = source(grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP);
  engine::Chart<S::value_type> chart(S::zero());
  for (std::size_t c = 0; c < d.child_count; ++c)
    chart.set(d.children[c].item, S::from_weight(2.0 + c));
  const auto states = graph.states().size();
  const auto terms = graph.terms().size();
  const auto factors = graph.factors().size();
  for (std::size_t repeat = 0; repeat < 20000; ++repeat) {
    const double local = 1.0 + repeat;
    const auto prefix = beam::evaluate_wmbp_be_vp_prefix<S>(d, chart, local);
    if (!prefix || !close(prefix->log_value, std::log(local * 2.0 * 3.0))) return false;
  }
  if (graph.states().size() != states || graph.terms().size() != terms ||
      graph.factors().size() != factors) return false;
  const auto root = semirings::checkpoint_value(engine::evaluate_deduction<S>(
      d, chart, [](const auto &) { return 7.0; }));
  const auto beta = graph.outside_log_beta(root.node);
  for (std::size_t c = 0; c < d.child_count; ++c) {
    const auto child = chart.get(d.children[c].item);
    if (!close(beta[child.node], std::log(7.0 * 2 * 3 * 4 / (2.0 + c)))) return false;
  }
  if (graph.states().size() != states + 1 ||
      graph.terms().size() != terms + 1 || graph.factors().size() != factors + 3)
    return false;
  std::cout << "ranking products=20000 added_states=0 selected_added_states=1\n";
  return true;
}

bool check_parent_selection_lifetime() {
  constexpr std::size_t width = 20, candidates = 80, parents = 2000;
  beam::TopWidthWmbpBeVpPrefixSelection selection(width);
  for (std::size_t p = 0; p < parents; ++p) {
    auto d = source(grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP);
    d.parent.j += static_cast<cand_pos_t>(p);
    std::vector<beam::ScoredWmbpBeVpPrefix> scored;
    for (std::size_t c = 0; c < candidates; ++c) {
      d.children[0].item.i = 2 + static_cast<cand_pos_t>(c);
      scored.push_back({beam::wmbp_be_vp_prefix_frame_key(d), double(c)});
    }
    const auto keys = scored;
    selection.commit_parent(d.parent, std::move(scored));
    if (selection.live_selection_count() != width ||
        selection.live_parent_count() != 1) return false;
    for (std::size_t c = 0; c < candidates; ++c)
      if (selection.allows(keys[c].key) != (c >= candidates - width)) return false;
    selection.release_completed_parents();
    if (selection.live_selection_count() != 0 || selection.live_parent_count() != 0 ||
        selection.candidate_count() != (p + 1) * candidates ||
        selection.selected_count() != (p + 1) * width ||
        selection.pruned_count() != (p + 1) * (candidates - width) ||
        !selection.respects_width() || selection.allows(keys.back().key)) return false;
    // Releasing an already empty selection must not recount diagnostics.
    selection.release_completed_parents();
    if (selection.selected_count() != (p + 1) * width) return false;
  }
  std::cout << "selection parents=2000 width=20 max_live_keys=20 final_live_keys=0\n";
  return true;
}
}  // namespace

int main() {
  using namespace cparty::scfg2::semirings;
  if (!check_values<SumProductSemiring>() || !check_values<MaxProductSemiring>() ||
      !check_values<LogSumExpSemiring>() || !check_values<LogMaxSemiring>() ||
      !check_values<CoupledLogSumMaxSemiring>() ||
      !check_values<TargetedCoupledLogSumMaxSemiring>() ||
      !check_transient_graph_storage() || !check_parent_selection_lifetime())
    return EXIT_FAILURE;
  std::cout << "scfg2_wmbp_prefix_recompute=ok\n";
  return EXIT_SUCCESS;
}

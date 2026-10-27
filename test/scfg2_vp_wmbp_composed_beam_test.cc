#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_item_selection.hh"
#include "scfg2/beam/vp_wmbp_online_schedule.hh"
#include "scfg2/beam/wmbp_be_vp_prefix_events.hh"
#include "scfg2/beam/wmbp_item_selection.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

void check_vpl_padding_interval_shape(const std::string &sequence,
                                      const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(sequence.size());
  for (cand_pos_t j = 1; j <= n; ++j) {
    if (ctx.pair_round(j) >= 0) continue;
    for (cand_pos_t k = 2; k < j; ++k) {
      const cand_pos_t left =
          std::max<cand_pos_t>(1, k - ctx.unpaired_prefix(k - 1));
      bool initialized = false;
      bool expected_border = false;
      for (cand_pos_t i = left; i < k; ++i) {
        const auto cap = std::min(
            static_cast<cand_pos_tu>(ctx.border_b(i, j)),
            static_cast<cand_pos_tu>(ctx.border_Bp(i, j)));
        const bool border = static_cast<cand_pos_tu>(k) < cap;
        if (!initialized) {
          initialized = true;
          expected_border = border;
        } else {
          assert(border == expected_border);
        }
      }
    }
  }
}

void check_vpr_suffix_border_shape(const std::string &sequence,
                                   const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(sequence.size());
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  const auto chart = cparty::scfg2::replay::run_w_final_exact_dp_schedule<
      cparty::scfg2::semirings::SumProductSemiring>(
      schedule, ctx,
      [](const Deduction &) { return true; },
      [](const ItemKey &, const Deduction &, double, const auto &, double) {});
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = 1; i < j; ++i) {
      const ItemKey parent{
          cparty::scfg2::grammar::NonTerminal::VPR, i, j, -1, -1};
      if (!cparty::scfg2::deductions::exact_basic_item_parent_active(
              parent, ctx)) {
        continue;
      }
      const cand_pos_t max_bp =
          std::max(ctx.border_B(i, j), ctx.border_bp(i, j));
      bool initialized = false;
      bool expected_border = false;
      for (cand_pos_t k = i + 1; k < j; ++k) {
        if (!ctx.can_pair_right_span(k, j)) continue;
        if (chart.get(ItemKey{
                cparty::scfg2::grammar::NonTerminal::VP,
                i, k, -1, -1}) == 0.0) continue;
        const bool border = k > max_bp;
        if (!initialized) {
          initialized = true;
          expected_border = border;
        } else {
          assert(border == expected_border);
        }
      }
    }
  }
}

bool same_deduction(const Deduction &lhs, const Deduction &rhs) {
  if (!cparty::scfg2::deductions::ItemKeyEq{}(lhs.parent, rhs.parent) ||
      lhs.rule != rhs.rule || lhs.split.k != rhs.split.k ||
      lhs.split.l != rhs.split.l || lhs.split.p != rhs.split.p ||
      lhs.split.q != rhs.split.q || lhs.child_count != rhs.child_count) {
    return false;
  }
  for (std::size_t idx = 0; idx < lhs.child_count; ++idx) {
    if (!cparty::scfg2::deductions::ItemKeyEq{}(
            lhs.children[idx].item, rhs.children[idx].item)) {
      return false;
    }
  }
  return true;
}

bool is_wmbp_phase_item(const ItemKey &item) {
  const auto nt = item.nonterminal;
  return nt == cparty::scfg2::grammar::NonTerminal::VP ||
         nt == cparty::scfg2::grammar::NonTerminal::BE ||
         nt == cparty::scfg2::grammar::NonTerminal::WMBP;
}

bool is_wmbp_consumer_source(const Deduction &deduction) {
  if (deduction.parent.nonterminal ==
      cparty::scfg2::grammar::NonTerminal::WMBP) {
    return false;
  }
  for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
    if (deduction.children[idx].item.nonterminal ==
        cparty::scfg2::grammar::NonTerminal::WMBP) {
      return true;
    }
  }
  return false;
}

struct ComposedRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t vp_candidates = 0;
  std::size_t vp_selected = 0;
  std::size_t wmbp_candidates = 0;
  std::size_t wmbp_selected = 0;
  std::size_t rejected_vp_sources = 0;
  std::size_t rejected_wmbp_sources = 0;
  std::size_t driven_vp_sources = 0;
  std::size_t driven_wmbp_sources = 0;
  std::size_t wmbp_parent_visits = 0;
  std::size_t wmbp_source_parents = 0;
  std::size_t wmbp_split_source_parents = 0;
  std::size_t wmbp_direct_source_parents = 0;
  std::size_t wmbp_direct_only_parents = 0;
  std::size_t sparse_wmbp_parent_visits = 0;
  std::size_t reverse_wmbp_parent_keys = 0;
  std::size_t reverse_parent_vp_probes = 0;
  std::size_t case3_nonleft_wi_queries = 0;
  std::size_t case3_nonleft_geometry_probes = 0;
  std::size_t case3_nonleft_right_last_probes = 0;
  std::size_t case3_nonleft_middle_last_probes = 0;
  std::size_t case3_nonleft_zero_child_probes = 0;
  std::size_t case3_nonleft_emitted_sources = 0;
  std::size_t case3_nonleft_prefix_candidates = 0;
  std::size_t case3_nonleft_active_prefix_probes = 0;
  std::size_t case3_nonleft_max_active_prefixes = 0;
  std::vector<ItemKey> wmbp_source_parent_keys;
  std::vector<ItemKey> reverse_wmbp_parent_key_items;
  std::size_t sparse_vp_parent_visits = 0;
  std::size_t sparse_be_parent_visits = 0;
  std::size_t keepall_be_parent_visits = 0;
  std::size_t v_parent_visits = 0;
};

struct OnlineSparseRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t vp_parent_visits = 0;
  std::size_t be_parent_visits = 0;
  std::size_t wmbp_parent_visits = 0;
  std::size_t reverse_parent_vp_probes = 0;
  std::size_t case3_nonleft_wi_queries = 0;
  std::size_t case3_nonleft_geometry_probes = 0;
  std::size_t case3_nonleft_emitted_sources = 0;
  std::size_t v_parent_visits = 0;
};

template <class Semiring>
OnlineSparseRun run_online_sparse(const std::string &sequence,
                                  const std::string &scaffold,
                                  std::size_t width) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  auto run = cparty::scfg2::replay::run_w_final_vp_wmbp_online<Semiring>(
      ctx, static_cast<cand_pos_t>(sequence.size()), width);
  return OnlineSparseRun{std::move(run.chart),
                         run.schedule.vp_parent_visits,
                         run.schedule.be_parent_visits,
                         run.schedule.wmbp_parent_visits,
                         run.schedule.reverse_parent_vp_probes,
                         run.vp_case3_nonleft_wi_queries,
                         run.vp_case3_nonleft_geometry_probes,
                         run.vp_case3_nonleft_emitted_sources,
                         run.schedule.family_parent_visits[
                             static_cast<std::size_t>(
                                 cparty::scfg2::grammar::NonTerminal::V)]};
}

template <class Semiring>
cparty::scfg2::replay::WFinalVpWmbpOnlineRun<
    typename Semiring::value_type>
run_online_with_widths(
    const std::string &sequence,
    const std::string &scaffold,
    const cparty::scfg2::replay::WFinalVpWmbpOnlineWidths &widths) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  return cparty::scfg2::replay::run_w_final_vp_wmbp_online<Semiring>(
      ctx, static_cast<cand_pos_t>(sequence.size()), widths);
}

template <class Semiring>
void check_online_family_frontiers(const std::string &sequence,
                                   const std::string &scaffold,
                                   bool require_frontier_pruning) {
  using cparty::scfg2::replay::WFinalVpWmbpOnlineWidths;
  const auto keep_all = run_online_with_widths<Semiring>(
      sequence, scaffold, WFinalVpWmbpOnlineWidths{});

  auto case3_one_widths = WFinalVpWmbpOnlineWidths{};
  case3_one_widths.case3_producer = 1;
  case3_one_widths.case3_tail = 1;
  case3_one_widths.case3_left = 1;
  const auto case3_one = run_online_with_widths<Semiring>(
      sequence, scaffold, case3_one_widths);
  auto be_one_widths = WFinalVpWmbpOnlineWidths{};
  be_one_widths.be_right_tail = 1;
  be_one_widths.be_right_shell = 1;
  const auto be_one = run_online_with_widths<Semiring>(
      sequence, scaffold, be_one_widths);
  auto prefix_one_widths = WFinalVpWmbpOnlineWidths{};
  prefix_one_widths.wmbp_prefix = 1;
  const auto prefix_one = run_online_with_widths<Semiring>(
      sequence, scaffold, prefix_one_widths);
  auto nonleft_one_widths = WFinalVpWmbpOnlineWidths{};
  nonleft_one_widths.case3_nonleft_prefix = 1;
  const auto nonleft_one = run_online_with_widths<Semiring>(
      sequence, scaffold, nonleft_one_widths);
  auto nonleft_two_widths = WFinalVpWmbpOnlineWidths{};
  nonleft_two_widths.case3_nonleft_prefix = 2;
  const auto nonleft_two = run_online_with_widths<Semiring>(
      sequence, scaffold, nonleft_two_widths);
  auto case12_one_widths = WFinalVpWmbpOnlineWidths{};
  case12_one_widths.case12_prefix = 1;
  const auto case12_one = run_online_with_widths<Semiring>(
      sequence, scaffold, case12_one_widths);
  auto case12_two_widths = WFinalVpWmbpOnlineWidths{};
  case12_two_widths.case12_prefix = 2;
  const auto case12_two = run_online_with_widths<Semiring>(
      sequence, scaffold, case12_two_widths);
  const auto uniform_one = run_online_with_widths<Semiring>(
      sequence, scaffold, WFinalVpWmbpOnlineWidths::uniform(1));

  const auto vp_idx = static_cast<std::size_t>(
      cparty::scfg2::grammar::NonTerminal::VP);
  const auto be_idx = static_cast<std::size_t>(
      cparty::scfg2::grammar::NonTerminal::BE);
  const std::size_t full_vp_triangle =
      sequence.size() * (sequence.size() + 1) / 2;
  assert(keep_all.schedule.vp_parent_visits > 0);
  assert(keep_all.schedule.vp_parent_visits < full_vp_triangle);
  assert(keep_all.item_visits[vp_idx] <= keep_all.schedule.vp_parent_visits);
  assert(keep_all.schedule.be_parent_visits > 0);
  assert(keep_all.item_visits[be_idx] <= keep_all.schedule.be_parent_visits);
  assert(uniform_one.schedule.be_parent_visits <=
         keep_all.schedule.be_parent_visits);
  assert(uniform_one.item_visits[vp_idx] <=
         uniform_one.schedule.vp_parent_visits);

  for (const auto &[item, exact] : keep_all.chart.values()) {
    const double case3_value = case3_one.chart.get(item);
    const double be_value = be_one.chart.get(item);
    const double prefix_value = prefix_one.chart.get(item);
    const double nonleft_one_value = nonleft_one.chart.get(item);
    const double nonleft_two_value = nonleft_two.chart.get(item);
    const double case12_one_value = case12_one.chart.get(item);
    const double case12_two_value = case12_two.chart.get(item);
    const double uniform_value = uniform_one.chart.get(item);
    assert(case3_value <= exact || close(case3_value, exact));
    assert(be_value <= exact || close(be_value, exact));
    assert(prefix_value <= exact || close(prefix_value, exact));
    assert(nonleft_one_value <= nonleft_two_value ||
           close(nonleft_one_value, nonleft_two_value));
    assert(nonleft_two_value <= exact || close(nonleft_two_value, exact));
    assert(case12_one_value <= case12_two_value ||
           close(case12_one_value, case12_two_value));
    assert(case12_two_value <= exact || close(case12_two_value, exact));
    assert(uniform_value <= exact || close(uniform_value, exact));
  }

  assert(keep_all.case3_producer_candidates ==
         keep_all.case3_producer_selected);
  assert(keep_all.case3_tail_candidates == keep_all.case3_tail_selected);
  assert(keep_all.case3_left_candidates == keep_all.case3_left_selected);
  assert(keep_all.be_right_tail_candidates ==
         keep_all.be_right_tail_selected);
  assert(keep_all.be_right_shell_candidates ==
         keep_all.be_right_shell_selected);
  assert(keep_all.wmbp_prefix_candidates ==
         keep_all.wmbp_prefix_selected);
  assert(case3_one.case3_producer_selected <=
         case3_one.case3_producer_candidates);
  assert(case3_one.case3_tail_selected <= case3_one.case3_tail_candidates);
  assert(case3_one.case3_left_selected <= case3_one.case3_left_candidates);
  assert(be_one.be_right_tail_selected <= be_one.be_right_tail_candidates);
  assert(be_one.be_right_shell_selected <= be_one.be_right_shell_candidates);
  assert(prefix_one.wmbp_prefix_selected <=
         prefix_one.wmbp_prefix_candidates);
  assert(nonleft_one.vp_case3_nonleft_max_active_prefixes <= 1);
  assert(nonleft_two.vp_case3_nonleft_max_active_prefixes <= 2);
  assert(uniform_one.case3_producer_selected <=
         uniform_one.case3_producer_candidates);
  assert(uniform_one.case3_tail_selected <=
         uniform_one.case3_tail_candidates);
  assert(uniform_one.case3_left_selected <=
         uniform_one.case3_left_candidates);
  assert(keep_all.vp_case3_nonleft_border_mismatches == 0);
  assert(case3_one.vp_case3_nonleft_border_mismatches == 0);
  assert(be_one.vp_case3_nonleft_border_mismatches == 0);
  assert(prefix_one.vp_case3_nonleft_border_mismatches == 0);
  assert(nonleft_one.vp_case3_nonleft_border_mismatches == 0);
  assert(nonleft_two.vp_case3_nonleft_border_mismatches == 0);
  assert(uniform_one.vp_case3_nonleft_border_mismatches == 0);
  assert(keep_all.vp_case12_border_mismatches == 0);
  assert(case12_one.vp_case12_border_mismatches == 0);
  assert(case12_two.vp_case12_border_mismatches == 0);
  assert(uniform_one.vp_case12_border_mismatches == 0);
  assert(keep_all.vp_case3_nonleft_duplicate_sources == 0);
  assert(keep_all.vp_case12_duplicate_sources == 0);
  assert(case12_one.vp_case12_duplicate_sources == 0);
  assert(case12_two.vp_case12_duplicate_sources == 0);
  assert(uniform_one.vp_case12_duplicate_sources == 0);
  assert(keep_all.vp_recursive_split_duplicate_sources == 0);
  assert(case12_one.vp_recursive_split_duplicate_sources == 0);
  assert(case12_two.vp_recursive_split_duplicate_sources == 0);
  assert(uniform_one.vp_recursive_split_duplicate_sources == 0);
  assert(keep_all.be_base_stack_duplicate_sources == 0);
  assert(case3_one.be_base_stack_duplicate_sources == 0);
  assert(be_one.be_base_stack_duplicate_sources == 0);
  assert(prefix_one.be_base_stack_duplicate_sources == 0);
  assert(case12_one.be_base_stack_duplicate_sources == 0);
  assert(case12_two.be_base_stack_duplicate_sources == 0);
  assert(nonleft_one.be_base_stack_duplicate_sources == 0);
  assert(nonleft_two.be_base_stack_duplicate_sources == 0);
  assert(uniform_one.be_base_stack_duplicate_sources == 0);
  assert(keep_all.be_core_left_duplicate_sources == 0);
  assert(case3_one.be_core_left_duplicate_sources == 0);
  assert(be_one.be_core_left_duplicate_sources == 0);
  assert(prefix_one.be_core_left_duplicate_sources == 0);
  assert(case12_one.be_core_left_duplicate_sources == 0);
  assert(case12_two.be_core_left_duplicate_sources == 0);
  assert(nonleft_one.be_core_left_duplicate_sources == 0);
  assert(nonleft_two.be_core_left_duplicate_sources == 0);
  assert(uniform_one.be_core_left_duplicate_sources == 0);
  assert(keep_all.be_right_tail_duplicate_registrations == 0);
  assert(case3_one.be_right_tail_duplicate_registrations == 0);
  assert(be_one.be_right_tail_duplicate_registrations == 0);
  assert(prefix_one.be_right_tail_duplicate_registrations == 0);
  assert(case12_one.be_right_tail_duplicate_registrations == 0);
  assert(case12_two.be_right_tail_duplicate_registrations == 0);
  assert(nonleft_one.be_right_tail_duplicate_registrations == 0);
  assert(nonleft_two.be_right_tail_duplicate_registrations == 0);
  assert(uniform_one.be_right_tail_duplicate_registrations == 0);
  assert(keep_all.be_right_shell_duplicate_sources == 0);
  assert(case3_one.be_right_shell_duplicate_sources == 0);
  assert(be_one.be_right_shell_duplicate_sources == 0);
  assert(prefix_one.be_right_shell_duplicate_sources == 0);
  assert(case12_one.be_right_shell_duplicate_sources == 0);
  assert(case12_two.be_right_shell_duplicate_sources == 0);
  assert(nonleft_one.be_right_shell_duplicate_sources == 0);
  assert(nonleft_two.be_right_shell_duplicate_sources == 0);
  assert(uniform_one.be_right_shell_duplicate_sources == 0);
  assert(case3_one.vp_case3_nonleft_duplicate_sources == 0);
  assert(be_one.vp_case3_nonleft_duplicate_sources == 0);
  assert(prefix_one.vp_case3_nonleft_duplicate_sources == 0);
  assert(nonleft_one.vp_case3_nonleft_duplicate_sources == 0);
  assert(nonleft_two.vp_case3_nonleft_duplicate_sources == 0);
  assert(uniform_one.vp_case3_nonleft_duplicate_sources == 0);
  assert(uniform_one.vp_case3_nonleft_max_active_prefixes <= 1);
  assert(uniform_one.be_right_tail_selected <=
         uniform_one.be_right_tail_candidates);
  assert(uniform_one.be_right_shell_selected <=
         uniform_one.be_right_shell_candidates);
  assert(uniform_one.wmbp_prefix_selected <=
         uniform_one.wmbp_prefix_candidates);
  for (const auto nt : {
           cparty::scfg2::grammar::NonTerminal::VM,
           cparty::scfg2::grammar::NonTerminal::V,
           cparty::scfg2::grammar::NonTerminal::WMv,
           cparty::scfg2::grammar::NonTerminal::WMB,
           cparty::scfg2::grammar::NonTerminal::WMp,
           cparty::scfg2::grammar::NonTerminal::WM,
           cparty::scfg2::grammar::NonTerminal::WI,
           cparty::scfg2::grammar::NonTerminal::WIP,
           cparty::scfg2::grammar::NonTerminal::WMBW,
           cparty::scfg2::grammar::NonTerminal::VPL,
           cparty::scfg2::grammar::NonTerminal::VPR,
           cparty::scfg2::grammar::NonTerminal::W,
           cparty::scfg2::grammar::NonTerminal::BE}) {
    const auto idx = static_cast<std::size_t>(nt);
    assert(keep_all.frontier_candidates[idx] ==
           keep_all.frontier_selected[idx]);
    assert(uniform_one.frontier_selected[idx] <=
           uniform_one.frontier_candidates[idx]);
  }
  const auto wmv_idx = static_cast<std::size_t>(
      cparty::scfg2::grammar::NonTerminal::WMv);
  const auto wmp_idx = static_cast<std::size_t>(
      cparty::scfg2::grammar::NonTerminal::WMp);
  assert(uniform_one.schedule.family_parent_visits[wmv_idx] <=
         2 * sequence.size());
  assert(uniform_one.schedule.family_parent_visits[wmp_idx] <=
         2 * sequence.size());
  assert(uniform_one.vm_binary_join_probes <= 2 * sequence.size());
  assert(uniform_one.wm_binary_join_probes <= 2 * sequence.size());
  assert(uniform_one.wm_extension_probes <= sequence.size());
  assert(uniform_one.wi_binary_join_probes <= 2 * sequence.size());
  assert(uniform_one.wi_direct_probes <= 2 * sequence.size());
  assert(uniform_one.wi_extension_probes <= sequence.size());
  assert(uniform_one.wip_binary_join_probes <= 2 * sequence.size());
  assert(uniform_one.wip_extension_probes <= sequence.size());
  assert(uniform_one.wmbw_binary_join_probes <= sequence.size());
  assert(uniform_one.vpr_binary_join_probes <= sequence.size());
  assert(uniform_one.vpr_base_state_probes <= 2 * sequence.size());
  assert(uniform_one.vpr_base_state_selected <= sequence.size());
  assert(uniform_one.w_binary_join_probes <= 4 * sequence.size());
  assert(uniform_one.w_direct_probes <= 2 * sequence.size());
  assert(uniform_one.w_extension_probes <= 2 * sequence.size());
  assert(uniform_one.w_root_decomposition_probes <= 10 * sequence.size());
  assert(uniform_one.wmb_btail_join_probes <= sequence.size());
  assert(uniform_one.wmb_btail_direct_probes <= sequence.size());
  assert(uniform_one.wmb_btail_selected <= sequence.size());
  assert(uniform_one.wmb_be_btail_join_probes <= sequence.size());
  assert(uniform_one.wmb_retained_direct_probes <= sequence.size());
  assert(uniform_one.wmb_retained_split_selected <= sequence.size());
  assert(uniform_one.wmbp_shared_wmbp_left_probes <= sequence.size());
  assert(uniform_one.wmbp_shared_wmbw_left_probes <= sequence.size());
  assert(uniform_one.wmbp_shared_emitted_sources <= 2 * sequence.size());
  assert(uniform_one.wmbp_wi_history_insertions <= sequence.size());
  assert(uniform_one.wmbp_wi_be_wi_probes <= sequence.size());
  assert(uniform_one.wmbp_wi_empty_probes <= sequence.size());
  assert(uniform_one.wmbp_wi_prefix_candidates <= 2 * sequence.size());
  assert(uniform_one.wmbp_wi_prefix_vp_probes <= sequence.size());
  assert(uniform_one.wmbp_wi_direct_probes <= sequence.size());
  assert(uniform_one.wmbp_wi_emitted_sources <= 2 * sequence.size());
  assert(uniform_one.vp_loop_parent_probes <=
         sequence.size() * (MAXLOOP + 1) * (MAXLOOP + 1));
  assert(uniform_one.vp_stack_event_sources +
             uniform_one.vp_internal_event_sources <=
         2 * uniform_one.vp_loop_parent_probes);
  assert(uniform_one.case3_max_active_segments <= 1);
  assert(uniform_one.case3_tail_materialization_probes <= sequence.size());
  assert(uniform_one.case3_prepared_score_probes <= sequence.size());
  assert(uniform_one.v_hairpin_parent_probes <=
         sequence.size() * (30 + 2));
  assert(uniform_one.v_internal_parent_probes <=
         sequence.size() * (MAXLOOP + 1) * (MAXLOOP + 1));
  assert(uniform_one.v_vm_parent_probes <= sequence.size());
  assert(uniform_one.case3_heap_registered_segments <= 30 * sequence.size());
  assert(uniform_one.case3_heap_expired_segments <=
         uniform_one.case3_heap_registered_segments);
  assert(uniform_one.case3_max_heap_size <=
         uniform_one.case3_heap_registered_segments);
  assert(uniform_one.case3_endpoint_prefix_registrations <=
         30 * sequence.size());
  assert(uniform_one.case3_endpoint_prefix_late_registrations == 0);
  assert(uniform_one.case3_endpoint_prefix_nonmonotone_queries == 0);
  assert(uniform_one.case3_endpoint_prefix_maximum_size <= 1);
  assert(uniform_one.case3_endpoint_prefix_maximum_returned_per_path <= 1);
  for (std::size_t idx = 0;
       idx < cparty::scfg2::replay::w_final_online_nonterminal_count;
       ++idx) {
    assert(keep_all.positive_items[idx] <= keep_all.item_visits[idx]);
  }

  if (require_frontier_pruning) {
    const std::size_t triangular =
        sequence.size() * (sequence.size() + 1) / 2;
    for (const auto nt : {
             cparty::scfg2::grammar::NonTerminal::WI,
             cparty::scfg2::grammar::NonTerminal::V,
             cparty::scfg2::grammar::NonTerminal::VM,
             cparty::scfg2::grammar::NonTerminal::WMv,
             cparty::scfg2::grammar::NonTerminal::WMp,
             cparty::scfg2::grammar::NonTerminal::WM,
             cparty::scfg2::grammar::NonTerminal::WIP,
             cparty::scfg2::grammar::NonTerminal::VP,
             cparty::scfg2::grammar::NonTerminal::WMBW,
             cparty::scfg2::grammar::NonTerminal::WMB}) {
      assert(keep_all.item_visits[static_cast<std::size_t>(nt)] <
             triangular);
    }
    for (const auto nt : {
             cparty::scfg2::grammar::NonTerminal::WMv,
             cparty::scfg2::grammar::NonTerminal::WMp,
             cparty::scfg2::grammar::NonTerminal::WM,
             cparty::scfg2::grammar::NonTerminal::VM}) {
      const auto idx = static_cast<std::size_t>(nt);
      assert(keep_all.item_visits[idx] == keep_all.positive_items[idx]);
    }
    assert(case3_one.case3_producer_selected <
               case3_one.case3_producer_candidates ||
           case3_one.case3_tail_selected < case3_one.case3_tail_candidates ||
           case3_one.case3_left_selected < case3_one.case3_left_candidates);
    assert(be_one.be_right_tail_selected <
               be_one.be_right_tail_candidates ||
           be_one.be_right_shell_selected <
               be_one.be_right_shell_candidates);
    assert(prefix_one.wmbp_prefix_selected <
           prefix_one.wmbp_prefix_candidates);
    bool pruned_wm_block = false;
    for (const auto nt : {
             cparty::scfg2::grammar::NonTerminal::VM,
             cparty::scfg2::grammar::NonTerminal::V,
             cparty::scfg2::grammar::NonTerminal::WMv,
             cparty::scfg2::grammar::NonTerminal::WMB,
             cparty::scfg2::grammar::NonTerminal::WMp,
             cparty::scfg2::grammar::NonTerminal::WM,
             cparty::scfg2::grammar::NonTerminal::WI,
             cparty::scfg2::grammar::NonTerminal::WIP,
             cparty::scfg2::grammar::NonTerminal::WMBW,
             cparty::scfg2::grammar::NonTerminal::VPL,
             cparty::scfg2::grammar::NonTerminal::VPR,
             cparty::scfg2::grammar::NonTerminal::W,
             cparty::scfg2::grammar::NonTerminal::BE}) {
      const auto idx = static_cast<std::size_t>(nt);
      pruned_wm_block |= uniform_one.frontier_selected[idx] <
                         uniform_one.frontier_candidates[idx];
    }
    assert(pruned_wm_block);
    assert(keep_all.vm_base_prefix_cell_visits > 0);
    assert(uniform_one.vm_base_prefix_cell_visits == 0);
    assert(uniform_one.vm_base_prefix_selected <= sequence.size());
    assert(keep_all.wm_start_prefix_cell_visits > 0);
    assert(uniform_one.wm_start_prefix_cell_visits == 0);
    assert(uniform_one.wm_start_prefix_selected <= sequence.size());
    assert(keep_all.wip_start_prefix_cell_visits > 0);
    assert(uniform_one.wip_start_prefix_cell_visits == 0);
    assert(uniform_one.wip_start_prefix_selected <= sequence.size());
    assert(keep_all.vpl_prefix_cell_visits > 0);
    assert(uniform_one.vpl_prefix_cell_visits == 0);
    assert(uniform_one.vpl_prefix_selected <= sequence.size());
  }
}

template <class Semiring, bool Driven>
ComposedRun run_composed(const std::string &sequence,
                         const std::string &scaffold,
                         std::size_t width,
                         const std::vector<ItemKey> *sparse_wmbp_parents =
                             nullptr) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_wmbp_barrier_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  if (sparse_wmbp_parents) {
    const std::unordered_set<ItemKey,
                             cparty::scfg2::deductions::ItemKeyHash,
                             cparty::scfg2::deductions::ItemKeyEq>
        retained(sparse_wmbp_parents->begin(), sparse_wmbp_parents->end());
    schedule.erase(
        std::remove_if(
            schedule.begin(), schedule.end(), [&](const ItemKey &item) {
              return item.nonterminal ==
                         cparty::scfg2::grammar::NonTerminal::WMBP &&
                     retained.find(item) == retained.end();
            }),
        schedule.end());
  }
  cparty::scfg2::beam::TopWidthVpItemSelection vp_selection(width);
  cparty::scfg2::beam::TopWidthWmbpItemSelection wmbp_selection(width);
  cand_pos_t vp_committed_j = -1;
  cand_pos_t wmbp_committed_j = -1;
  std::size_t rejected_vp_sources = 0;
  std::size_t rejected_wmbp_sources = 0;
  std::size_t driven_vp_sources = 0;
  std::size_t driven_wmbp_sources = 0;
  std::size_t wmbp_parent_visits = 0;
  std::size_t wmbp_source_parents = 0;
  std::size_t wmbp_split_source_parents = 0;
  std::size_t wmbp_direct_source_parents = 0;
  std::size_t wmbp_direct_only_parents = 0;
  std::vector<ItemKey> wmbp_source_parent_keys;

  const auto before_item = [&](const ItemKey &item, const auto &) {
    if (item.nonterminal != cparty::scfg2::grammar::NonTerminal::VP &&
        item.j != vp_committed_j) {
      vp_selection.commit_column(item.j);
      vp_committed_j = item.j;
    }
    if (!is_wmbp_phase_item(item) && item.j != wmbp_committed_j) {
      wmbp_selection.commit_column(item.j);
      wmbp_committed_j = item.j;
    }
  };
  const auto after_item = [&](const ItemKey &item, const auto &chart) {
    vp_selection.observe_completed(item, chart.get(item));
    wmbp_selection.observe_completed(item, chart.get(item));
  };

  auto chart = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_provider_item_observers_and_contribution<
          Semiring>(
          schedule,
          ctx,
          [&](const ItemKey &item) {
            auto exact = cparty::scfg2::deductions::
                generate_exact_basic_deductions_for_item(item, ctx);
            if constexpr (!Driven) return exact;

            std::vector<Deduction> driven;
            if (item.nonterminal ==
                cparty::scfg2::grammar::NonTerminal::WMBP) {
              if constexpr (Driven) ++wmbp_parent_visits;
              driven = cparty::scfg2::beam::
                  generate_wmbp_parent_deductions_for_retained_vp_starts(
                      item, vp_selection.selected_starts(item.j), ctx);
              if constexpr (Driven) {
                wmbp_source_parents += driven.empty() ? 0 : 1;
                if (!driven.empty()) wmbp_source_parent_keys.push_back(item);
                bool has_split = false;
                bool has_direct = false;
                for (const auto &deduction : driven) {
                  has_split = has_split ||
                      cparty::scfg2::beam::is_wmbp_be_vp_rule(
                          deduction.rule);
                  has_direct = has_direct ||
                      deduction.rule == cparty::scfg2::grammar::RuleId::
                          WMBP_DIRECT_VP;
                }
                wmbp_split_source_parents += has_split ? 1 : 0;
                wmbp_direct_source_parents += has_direct ? 1 : 0;
                wmbp_direct_only_parents += has_direct && !has_split ? 1 : 0;
              }
              for (const auto &deduction : driven) {
                driven_vp_sources +=
                    cparty::scfg2::beam::is_wmbp_be_vp_rule(deduction.rule)
                    ? 1
                    : 0;
              }
              exact.erase(
                  std::remove_if(
                      exact.begin(), exact.end(), [&](const Deduction &d) {
                        return !vp_selection.allows_wmbp_source(d);
                      }),
                  exact.end());
            } else if (item.nonterminal ==
                           cparty::scfg2::grammar::NonTerminal::WMBW ||
                       item.nonterminal ==
                           cparty::scfg2::grammar::NonTerminal::WMB) {
              driven = cparty::scfg2::deductions::
                  generate_exact_wmbp_consumer_deductions_for_retained_ends(
                      item, ctx, wmbp_selection.selected_ends(item.i));
              for (const auto &deduction : driven) {
                driven_wmbp_sources +=
                    is_wmbp_consumer_source(deduction) ? 1 : 0;
              }
              exact.erase(
                  std::remove_if(
                      exact.begin(), exact.end(), [&](const Deduction &d) {
                        return !wmbp_selection.allows_consumer_source(d);
                      }),
                  exact.end());
            } else {
              return exact;
            }
            assert(driven.size() == exact.size());
            for (std::size_t idx = 0; idx < driven.size(); ++idx) {
              assert(same_deduction(driven[idx], exact[idx]));
            }
            return driven;
          },
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          },
          before_item,
          after_item,
          [&](const Deduction &deduction, const auto &current_chart, double local) {
            if constexpr (!Driven) {
              if (!vp_selection.allows_wmbp_source(deduction)) {
                ++rejected_vp_sources;
                return Semiring::zero();
              }
              if (!wmbp_selection.allows_consumer_source(deduction)) {
                ++rejected_wmbp_sources;
                return Semiring::zero();
              }
            }
            return cparty::scfg2::replay::
                runtime_deduction_contribution<Semiring>(
                    current_chart, deduction, local);
          });
  std::size_t reverse_wmbp_parent_keys = 0;
  std::size_t reverse_parent_vp_probes = 0;
  std::vector<ItemKey> reverse_wmbp_parent_key_items;
  if constexpr (Driven) {
    std::unordered_set<ItemKey,
                       cparty::scfg2::deductions::ItemKeyHash,
                       cparty::scfg2::deductions::ItemKeyEq>
        reverse_keys;
    for (cand_pos_t j = 1; j <= static_cast<cand_pos_t>(sequence.size());
         ++j) {
      const auto discovery = cparty::scfg2::beam::
          discover_wmbp_source_parents_from_retained_vp_starts(
              j, vp_selection.selected_starts(j), ctx);
      reverse_parent_vp_probes += discovery.parent_vp_probes;
      for (const auto &parent : discovery.parents) reverse_keys.insert(parent);
    }
    const std::unordered_set<ItemKey,
                             cparty::scfg2::deductions::ItemKeyHash,
                             cparty::scfg2::deductions::ItemKeyEq>
        scanned_keys(wmbp_source_parent_keys.begin(),
                     wmbp_source_parent_keys.end());
    assert(reverse_keys.size() == scanned_keys.size());
    for (const auto &parent : reverse_keys) {
      assert(scanned_keys.find(parent) != scanned_keys.end());
    }
    reverse_wmbp_parent_keys = reverse_keys.size();
    reverse_wmbp_parent_key_items.assign(reverse_keys.begin(),
                                         reverse_keys.end());
  }
  assert(vp_selection.respects_width());
  assert(wmbp_selection.respects_width());
  return ComposedRun{std::move(chart),
                     vp_selection.candidate_count(),
                     vp_selection.selected_count(),
                     wmbp_selection.candidate_count(),
                     wmbp_selection.selected_count(),
                     rejected_vp_sources,
                     rejected_wmbp_sources,
                     driven_vp_sources,
                     driven_wmbp_sources,
                     wmbp_parent_visits,
                     wmbp_source_parents,
                     wmbp_split_source_parents,
                     wmbp_direct_source_parents,
                     wmbp_direct_only_parents,
                     0,
                     reverse_wmbp_parent_keys,
                     reverse_parent_vp_probes,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     std::move(wmbp_source_parent_keys),
                     std::move(reverse_wmbp_parent_key_items)};
}

template <class Semiring>
ComposedRun check_widths(const std::string &sequence,
                         const std::string &scaffold) {
  const auto unlimited = std::numeric_limits<std::size_t>::max();
  const auto filtered_all = run_composed<Semiring, false>(
      sequence, scaffold, unlimited);
  const auto driven_all = run_composed<Semiring, true>(
      sequence, scaffold, unlimited);
  const auto sparse_all = run_composed<Semiring, true>(
      sequence,
      scaffold,
      unlimited,
      &driven_all.reverse_wmbp_parent_key_items);
  const auto online_all =
      run_online_sparse<Semiring>(sequence, scaffold, unlimited);
  const auto filtered_one = run_composed<Semiring, false>(
      sequence, scaffold, 1);
  const auto driven_one = run_composed<Semiring, true>(
      sequence, scaffold, 1);
  const auto sparse_one = run_composed<Semiring, true>(
      sequence, scaffold, 1, &driven_one.reverse_wmbp_parent_key_items);
  const auto online_one = run_online_sparse<Semiring>(sequence, scaffold, 1);
  const auto online_uniform_one = run_online_with_widths<Semiring>(
      sequence,
      scaffold,
      cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(1));
  for (const auto &[item, exact] : filtered_all.chart.values()) {
    assert(close(driven_all.chart.get(item), exact));
    assert(close(sparse_all.chart.get(item), exact));
    assert(close(online_all.chart.get(item), exact));
    assert(close(filtered_one.chart.get(item), driven_one.chart.get(item)));
    assert(close(sparse_one.chart.get(item), driven_one.chart.get(item)));
    assert(driven_one.chart.get(item) <= exact ||
           close(driven_one.chart.get(item), exact));
    // Online recursive VP/WMBP beams can change which width-limited candidate
    // survives, so two different finite policies need not be cell-wise
    // ordered. Both must remain bounded by their shared KeepAll chart.
    assert(online_one.chart.get(item) <= exact ||
           close(online_one.chart.get(item), exact));
  }
  assert(filtered_one.vp_candidates == driven_one.vp_candidates);
  assert(filtered_one.vp_selected == driven_one.vp_selected);
  assert(filtered_one.wmbp_candidates == driven_one.wmbp_candidates);
  assert(filtered_one.wmbp_selected == driven_one.wmbp_selected);
  assert(sparse_one.vp_selected == driven_one.vp_selected);
  assert(sparse_one.wmbp_selected == driven_one.wmbp_selected);
  assert(sparse_one.wmbp_parent_visits ==
         driven_one.wmbp_source_parent_keys.size());
  assert(sparse_one.wmbp_source_parents ==
         sparse_one.wmbp_parent_visits);
  auto out = filtered_one;
  out.driven_vp_sources = driven_one.driven_vp_sources;
  out.driven_wmbp_sources = driven_one.driven_wmbp_sources;
  out.wmbp_parent_visits = driven_one.wmbp_parent_visits;
  out.wmbp_source_parents = driven_one.wmbp_source_parents;
  out.wmbp_split_source_parents = driven_one.wmbp_split_source_parents;
  out.wmbp_direct_source_parents = driven_one.wmbp_direct_source_parents;
  out.wmbp_direct_only_parents = driven_one.wmbp_direct_only_parents;
  out.sparse_wmbp_parent_visits = sparse_one.wmbp_parent_visits;
  out.reverse_wmbp_parent_keys = driven_one.reverse_wmbp_parent_keys;
  out.reverse_parent_vp_probes = driven_one.reverse_parent_vp_probes;
  out.case3_nonleft_wi_queries =
      online_uniform_one.vp_case3_nonleft_wi_queries;
  out.case3_nonleft_geometry_probes =
      online_uniform_one.vp_case3_nonleft_geometry_probes;
  out.case3_nonleft_right_last_probes =
      online_uniform_one.vp_case3_nonleft_right_last_probes;
  out.case3_nonleft_middle_last_probes =
      online_uniform_one.vp_case3_nonleft_middle_last_probes;
  out.case3_nonleft_zero_child_probes =
      online_uniform_one.vp_case3_nonleft_zero_child_probes;
  out.case3_nonleft_emitted_sources =
      online_uniform_one.vp_case3_nonleft_emitted_sources;
  out.case3_nonleft_prefix_candidates =
      online_uniform_one.vp_case3_nonleft_prefix_candidates;
  out.case3_nonleft_active_prefix_probes =
      online_uniform_one.vp_case3_nonleft_active_prefix_probes;
  out.case3_nonleft_max_active_prefixes =
      online_uniform_one.vp_case3_nonleft_max_active_prefixes;
  out.sparse_vp_parent_visits = online_one.vp_parent_visits;
  out.sparse_be_parent_visits = online_one.be_parent_visits;
  out.keepall_be_parent_visits = online_all.be_parent_visits;
  out.v_parent_visits = online_one.v_parent_visits;
  return out;
}

void check_online_viterbi_choices(const std::string &sequence,
                                  const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto widths =
      cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(1);
  auto run = cparty::scfg2::replay::run_w_final_vp_wmbp_online_viterbi(
      ctx, static_cast<cand_pos_t>(sequence.size()), widths);
  const ItemKey root{cparty::scfg2::grammar::NonTerminal::W,
                     1,
                     static_cast<cand_pos_t>(sequence.size()),
                     -1,
                     -1};
  const auto *root_choice = run.choices.find(root);
  if (run.dp.chart.get(root) == 0.0) {
    assert(!root_choice || !root_choice->has_value);
    return;
  }
  assert(root_choice && root_choice->has_value);
  assert(close(root_choice->contribution, run.dp.chart.get(root)));
  const auto trace = cparty::scfg2::engine::build_traceback_tree(
      root, run.dp.chart, run.choices);
  const auto check_node = [&](const auto &self, const auto &node) -> void {
    if (node.chart_value == 0.0) return;
    assert(node.chosen_deduction.has_value());
    assert(close(node.chosen_contribution, node.chart_value));
    for (const auto &child : node.children) self(self, child);
  };
  check_node(check_node, trace);

  auto traceback = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_viterbi_traceback(
          ctx, static_cast<cand_pos_t>(sequence.size()), widths);
  assert(cparty::scfg2::deductions::ItemKeyEq{}(traceback.root, root));
  assert(close(traceback.run.chart.get(root), run.dp.chart.get(root)));
  assert(traceback.lane_structures.round.size() == sequence.size());
  assert(traceback.lane_structures.square.size() == sequence.size());
  for (const auto &pair : traceback.generated_pairs) {
    assert(!cparty::scfg2::engine::is_fixed_scaffold_pair(pair, ctx));
  }
}

void check_logmax_online_parity(const std::string &sequence,
                                const std::string &scaffold) {
  using cparty::scfg2::semirings::LogMaxSemiring;
  using cparty::scfg2::semirings::MaxProductSemiring;
  using cparty::scfg2::replay::WFinalVpWmbpOnlineWidths;
  const auto widths = WFinalVpWmbpOnlineWidths::uniform(20);
  std::string max_sequence = sequence;
  std::string max_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot max_owner(
      max_sequence, max_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext max_ctx(max_owner);
  auto max_run = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_viterbi_with_semiring<MaxProductSemiring>(
          max_ctx, static_cast<cand_pos_t>(sequence.size()), widths);

  std::string log_sequence = sequence;
  std::string log_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot log_owner(
      log_sequence, log_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext log_ctx(log_owner);
  auto log_run = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_viterbi_with_semiring<LogMaxSemiring>(
          log_ctx, static_cast<cand_pos_t>(sequence.size()), widths);

  const ItemKey root{cparty::scfg2::grammar::NonTerminal::W,
                     1,
                     static_cast<cand_pos_t>(sequence.size()),
                     -1,
                     -1};
  const double max_root = max_run.dp.chart.get(root);
  const double log_root = log_run.dp.chart.get(root);
  assert(max_root > 0.0);
  if (!close(std::log(max_root), log_root)) {
    std::cerr.precision(17);
    std::cerr << "LogMax root mismatch max=" << std::log(max_root)
              << " log=" << log_root
              << " delta=" << std::log(max_root) - log_root << '\n';
    assert(false);
  }
  // Near-tied candidates at a finite-width beam boundary may produce
  // different unreachable internal charts after floating-point rounding.
  // The semantic contract is equality of the root score and its reachable
  // selected path, rather than identity of every discarded side chart.
  const auto max_trace = cparty::scfg2::engine::build_traceback_tree(
      root, max_run.dp.chart, max_run.choices);
  const auto log_trace = cparty::scfg2::engine::build_traceback_tree(
      root, log_run.dp.chart, log_run.choices);
  const auto check_same_selected_path = [&](const auto &self,
                                            const auto &max_node,
                                            const auto &log_node) -> void {
    assert(cparty::scfg2::deductions::ItemKeyEq{}(
        max_node.item, log_node.item));
    assert(max_node.chosen_deduction.has_value() ==
           log_node.chosen_deduction.has_value());
    if (!max_node.chosen_deduction.has_value()) return;
    if (!same_deduction(*max_node.chosen_deduction,
                        *log_node.chosen_deduction)) {
      const auto &max_deduction = *max_node.chosen_deduction;
      const auto &log_deduction = *log_node.chosen_deduction;
      std::cerr << "LogMax path mismatch nt="
                << static_cast<int>(max_node.item.nonterminal)
                << " i=" << max_node.item.i << " j=" << max_node.item.j
                << " max_rule=" << static_cast<int>(max_deduction.rule)
                << " max_split=" << max_deduction.split.k << ','
                << max_deduction.split.l << ',' << max_deduction.split.p
                << ',' << max_deduction.split.q
                << " log_rule=" << static_cast<int>(log_deduction.rule)
                << " log_split=" << log_deduction.split.k << ','
                << log_deduction.split.l << ',' << log_deduction.split.p
                << ',' << log_deduction.split.q
                << " max_log_contribution="
                << std::log(max_node.chosen_contribution)
                << " log_contribution=" << log_node.chosen_contribution
                << '\n';
      assert(false);
    }
    assert(max_node.children.size() == log_node.children.size());
    for (std::size_t idx = 0; idx < max_node.children.size(); ++idx) {
      self(self, max_node.children[idx], log_node.children[idx]);
    }
  };
  check_same_selected_path(check_same_selected_path, max_trace, log_trace);
}

void check_logsumexp_online_parity(const std::string &sequence,
                                   const std::string &scaffold) {
  using cparty::scfg2::replay::WFinalVpWmbpOnlineWidths;
  const auto widths = WFinalVpWmbpOnlineWidths::uniform(20);
  std::string product_sequence = sequence;
  std::string product_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot product_owner(
      product_sequence, product_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext product_ctx(product_owner);
  const auto product_run = cparty::scfg2::replay::run_w_final_vp_wmbp_online<
      cparty::scfg2::semirings::SumProductSemiring>(
      product_ctx, static_cast<cand_pos_t>(sequence.size()), widths);

  std::string log_sequence = sequence;
  std::string log_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot log_owner(
      log_sequence, log_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext log_ctx(log_owner);
  const auto log_run = cparty::scfg2::replay::run_w_final_vp_wmbp_online<
      cparty::scfg2::semirings::LogSumExpSemiring>(
      log_ctx, static_cast<cand_pos_t>(sequence.size()), widths);

  const ItemKey root{cparty::scfg2::grammar::NonTerminal::W,
                     1,
                     static_cast<cand_pos_t>(sequence.size()),
                     -1,
                     -1};
  const double product_root = product_run.chart.get(root);
  const double log_root = log_run.chart.get(root);
  assert(product_root > 0.0);
  if (!close(std::log(product_root), log_root)) {
    std::cerr.precision(17);
    std::cerr << "LogSumExp root mismatch product="
              << std::log(product_root) << " log=" << log_root
              << " delta=" << std::log(product_root) - log_root << '\n';
    assert(false);
  }
}

void check_coupled_online_bound(const std::string &sequence,
                                const std::string &scaffold) {
  using Semiring = cparty::scfg2::semirings::CoupledLogSumMaxSemiring;
  const auto widths =
      cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(20);
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto run = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_viterbi_with_semiring<Semiring>(
          ctx, static_cast<cand_pos_t>(sequence.size()), widths);
  const ItemKey root{cparty::scfg2::grammar::NonTerminal::W,
                     1,
                     static_cast<cand_pos_t>(sequence.size()),
                     -1,
                     -1};
  const auto value = run.dp.chart.get(root);
  assert(value.log_max <= value.log_sum + 1e-12);
  for (const auto &[item, cell] : run.dp.chart.values()) {
    (void)item;
    assert(cell.log_max <= cell.log_sum + 1e-12);
  }
  const auto trace = cparty::scfg2::engine::build_traceback_tree(
      root, run.dp.chart, run.choices);
  const auto check_trace = [&](const auto &self, const auto &node) -> void {
    if (Semiring::is_zero(node.chart_value)) return;
    assert(node.chosen_deduction.has_value());
    assert(close(node.chosen_contribution.log_max,
                 node.chart_value.log_max));
    for (const auto &child : node.children) self(self, child);
  };
  check_trace(check_trace, trace);
  const auto best_target =
      cparty::scfg2::engine::collect_generated_traceback_pairs(trace, ctx);

  std::string inside_sequence = sequence;
  std::string inside_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot inside_owner(
      inside_sequence, inside_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext inside_ctx(inside_owner);
  const auto inside = cparty::scfg2::replay::run_w_final_vp_wmbp_online<
      cparty::scfg2::semirings::LogSumExpSemiring>(
      inside_ctx, static_cast<cand_pos_t>(sequence.size()), widths);
  assert(close(value.log_sum, inside.chart.get(root)));

  std::size_t transformed_deductions = 0;
  const auto identity_transformed = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_observed_transformed<Semiring>(
          ctx,
          static_cast<cand_pos_t>(sequence.size()),
          widths,
          [](const ItemKey &,
             const cparty::scfg2::deductions::Deduction &,
             Semiring::value_type,
             const auto &,
             double) {},
          [&](const cparty::scfg2::deductions::Deduction &,
              Semiring::value_type contribution) {
            ++transformed_deductions;
            return contribution;
          });
  const auto transformed_root = identity_transformed.chart.get(root);
  assert(transformed_deductions > 0);
  assert(close(transformed_root.log_sum, value.log_sum));
  assert(close(transformed_root.log_max, value.log_max));

  using Targeted =
      cparty::scfg2::semirings::TargetedCoupledLogSumMaxSemiring;
  const auto target_rejected = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_observed_transformed<Targeted>(
          ctx,
          static_cast<cand_pos_t>(sequence.size()),
          widths,
          [](const ItemKey &,
             const cparty::scfg2::deductions::Deduction &,
             Targeted::value_type,
             const auto &,
             double) {},
          [](const cparty::scfg2::deductions::Deduction &,
             Targeted::value_type contribution) {
            return Targeted::reject_target(contribution);
          });
  const auto target_rejected_root = target_rejected.chart.get(root);
  assert(close(target_rejected_root.log_sum, value.log_sum));
  assert(close(target_rejected_root.log_max, value.log_max));
  assert(!std::isfinite(target_rejected_root.log_target));

  const auto targeted =
      cparty::scfg2::replay::run_w_final_vp_wmbp_online_targeted(
          ctx,
          static_cast<cand_pos_t>(sequence.size()),
          widths,
          best_target);
  const auto targeted_root = targeted.chart.get(root);
  assert(close(targeted_root.log_sum, value.log_sum));
  assert(close(targeted_root.log_max, value.log_max));
  assert(close(targeted_root.log_target, value.log_max));
  assert(targeted_root.log_target <= targeted_root.log_sum + 1e-12);

  const auto probability = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_target_probability(
          ctx,
          static_cast<cand_pos_t>(sequence.size()),
          widths,
          best_target);
  assert(probability.ok);
  assert(probability.denominator_reachable);
  assert(probability.target_included);
  assert(close(probability.log_partition, value.log_sum));
  assert(close(probability.log_target_weight, value.log_max));
  assert(close(probability.log_probability, value.log_max - value.log_sum));
  assert(probability.probability > 0.0);
  assert(probability.probability <= 1.0);

  const std::vector<cparty::scfg2::engine::TracebackPair> impossible_target{
      {1, 1, cparty::scfg2::engine::TracebackPairLane::Square}};
  const auto absent = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_target_probability(
          ctx,
          static_cast<cand_pos_t>(sequence.size()),
          widths,
          impossible_target);
  assert(absent.ok);
  assert(absent.denominator_reachable);
  assert(!absent.target_included);
  assert(!std::isfinite(absent.log_probability));
  assert(absent.probability == 0.0);
}

void check_target_probability_with_generated_scaffold_interior() {
  using namespace cparty::scfg2;
  std::string sequence =
      "UUUUACCCAGACCUAGUGACCGUCCUGUUUGGGCUAAUAUCGAAAGGUUUAUAACCUACGCUGUAAUGGCAGACCUGGGAGCGUUUGACGAUAUUU";
  std::string scaffold =
      "...................................(((((((..((((.....))))(((((.................)))))....))))))).";
  weights::RuntimeSnapshot owner(sequence, scaffold);
  replay::WFinalLogWeightContext ctx(owner);
  std::vector<engine::TracebackPair> target;
  for (const auto &[i, j] : std::vector<std::pair<int, int>>{
      {6,79},{7,78},{8,77},{9,76},{10,75},{13,72},{14,71},{16,70},
      {17,69},{19,67},{22,34},{23,33},{24,32},{25,31},{36,95},
      {37,94},{38,93},{39,92},{40,91},{41,90},{42,89},{45,57},
      {46,56},{47,55},{48,54},{49,53},{58,84},{59,83},{60,82},
      {61,81},{62,80}}) {
    target.push_back({i, j, engine::TracebackPairLane::Square});
  }
  const ItemKey root{grammar::NonTerminal::W, 1, 96, -1, -1};
  for (bool finite : {false, true}) {
    const auto widths = finite ? replay::WFinalVpWmbpOnlineWidths::uniform(30)
                               : replay::WFinalVpWmbpOnlineWidths{};
    const auto plain = replay::run_w_final_vp_wmbp_online<
        semirings::CoupledLogSumMaxSemiring>(ctx, 96, widths);
    const auto targeted = replay::run_w_final_vp_wmbp_online_targeted(
        ctx, 96, widths, target);
    const auto value = plain.chart.get(root);
    const auto constrained = targeted.chart.get(root);
    assert(close(value.log_sum, constrained.log_sum));
    assert(close(value.log_max, constrained.log_max));
    assert(constrained.log_target <= constrained.log_sum + 1e-10);
    if (!finite) {
      // A VP(19,67) crosses a scaffold component containing the independently
      // generated pair (49,53). Its target mass used to be falsely zero.
      assert(close(constrained.log_sum, 42.537535882221064));
      assert(close(constrained.log_target, 40.63955166800962));
    }
  }
}

}  // namespace

int main(int argc, char **argv) {
  assert(argc == 1 || argc == 3 || argc == 4);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 4 ? argv[3] : "params/rna_Turner04.par",
         1);
  if (argc == 1) check_target_probability_with_generated_scaffold_interior();
  const std::string sequence = argc >= 3
      ? argv[1]
      : "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  const std::string scaffold = argc >= 3
      ? argv[2]
      : "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  assert(sequence.size() == scaffold.size());
  check_vpl_padding_interval_shape(sequence, scaffold);
  check_vpr_suffix_border_shape(sequence, scaffold);
  const auto inside =
      check_widths<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi =
      check_widths<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  check_online_family_frontiers<
      cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold, argc == 1);
  check_online_family_frontiers<
      cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold, argc == 1);
  check_online_viterbi_choices(sequence, scaffold);
  check_logmax_online_parity(sequence, scaffold);
  check_logsumexp_online_parity(sequence, scaffold);
  check_coupled_online_bound(sequence, scaffold);
  (void)viterbi;
  if (argc == 1) {
    assert(inside.rejected_vp_sources > 0);
    assert(inside.rejected_wmbp_sources > 0);
  }
  std::cout << "n\tvp_candidates\tvp_selected"
               "\twmbp_candidates\twmbp_selected"
               "\trejected_vp_sources\trejected_wmbp_sources"
               "\tdriven_vp_sources\tdriven_wmbp_sources"
               "\twmbp_parent_visits\twmbp_source_parents"
               "\twmbp_split_source_parents\twmbp_direct_source_parents"
               "\twmbp_direct_only_parents\tsparse_wmbp_parent_visits"
               "\tsparse_vp_parent_visits"
               "\tsparse_be_parent_visits"
               "\tkeepall_be_parent_visits"
               "\tv_parent_visits"
               "\treverse_wmbp_parent_keys\treverse_parent_vp_probes"
               "\tcase3_nonleft_wi_queries\tcase3_nonleft_geometry_probes"
               "\tcase3_nonleft_right_last_probes"
               "\tcase3_nonleft_middle_last_probes"
               "\tcase3_nonleft_zero_child_probes"
               "\tcase3_nonleft_emitted_sources"
               "\tcase3_nonleft_prefix_candidates"
               "\tcase3_nonleft_active_prefix_probes"
               "\tcase3_nonleft_max_active_prefixes"
               "\tmonotone\n"
            << sequence.size() << '\t' << inside.vp_candidates << '\t'
            << inside.vp_selected << '\t' << inside.wmbp_candidates << '\t'
            << inside.wmbp_selected << '\t' << inside.rejected_vp_sources
            << '\t' << inside.rejected_wmbp_sources << '\t'
            << inside.driven_vp_sources << '\t'
            << inside.driven_wmbp_sources << '\t'
            << inside.wmbp_parent_visits << '\t'
            << inside.wmbp_source_parents << '\t'
            << inside.wmbp_split_source_parents << '\t'
            << inside.wmbp_direct_source_parents << '\t'
            << inside.wmbp_direct_only_parents << '\t'
            << inside.sparse_wmbp_parent_visits << '\t'
            << inside.sparse_vp_parent_visits << '\t'
            << inside.sparse_be_parent_visits << '\t'
            << inside.keepall_be_parent_visits << '\t'
            << inside.v_parent_visits << '\t'
            << inside.reverse_wmbp_parent_keys << '\t'
            << inside.reverse_parent_vp_probes << '\t'
            << inside.case3_nonleft_wi_queries << '\t'
            << inside.case3_nonleft_geometry_probes << '\t'
            << inside.case3_nonleft_right_last_probes << '\t'
            << inside.case3_nonleft_middle_last_probes << '\t'
            << inside.case3_nonleft_zero_child_probes << '\t'
            << inside.case3_nonleft_emitted_sources << '\t'
            << inside.case3_nonleft_prefix_candidates << '\t'
            << inside.case3_nonleft_active_prefix_probes << '\t'
            << inside.case3_nonleft_max_active_prefixes << "\t1\n";
}

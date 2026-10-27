#ifndef CPARTY_SCFG2_REPLAY_W_FINAL_VP_WMBP_ONLINE_HH
#define CPARTY_SCFG2_REPLAY_W_FINAL_VP_WMBP_ONLINE_HH

#include "scfg2/beam/public_item_ranking.hh"
#include "scfg2/beam/run_term_observer.hh"

#include "scfg2/beam/be_left_right_frame_events.hh"
#include "scfg2/beam/be_base_stack_retained_events.hh"
#include "scfg2/beam/be_core_left_prefix_agenda.hh"
#include "scfg2/beam/be_item_selection.hh"
#include "scfg2/beam/be_right_tail_retained_core_index.hh"
#include "scfg2/beam/be_right_shell_source_agenda.hh"
#include "scfg2/beam/family_item_selection.hh"
#include "scfg2/beam/v_retained_agenda.hh"
#include "scfg2/beam/vp_item_selection.hh"
#include "scfg2/beam/vp_loop_retained_events.hh"
#include "scfg2/beam/vp_recursive_split_retained_events.hh"
#include "scfg2/beam/vp_case3_nonleft_last_child_index.hh"
#include "scfg2/beam/vp_case3_nonleft_prefix_interval.hh"
#include "scfg2/beam/vp_case12_prefix_interval.hh"
#include "scfg2/beam/vp_case3_left_selection.hh"
#include "scfg2/beam/vp_case3_tail_selection.hh"
#include "scfg2/beam/vp_case3_wi_frontier.hh"
#include "scfg2/beam/vp_wmbp_online_schedule.hh"
#include "scfg2/beam/vpl_retained_agenda.hh"
#include "scfg2/beam/vpr_retained_agenda.hh"
#include "scfg2/beam/w_retained_join_agenda.hh"
#include "scfg2/beam/w_item_selection.hh"
#include "scfg2/beam/vm_live_parent.hh"
#include "scfg2/beam/vm_retained_join_agenda.hh"
#include "scfg2/beam/wmbp_be_vp_prefix_events.hh"
#include "scfg2/beam/wmbp_item_selection.hh"
#include "scfg2/beam/wmbp_shared_retained_join.hh"
#include "scfg2/beam/wmbp_wi_prefix_store.hh"
#include "scfg2/beam/wmbp_wi_retained_join.hh"
#include "scfg2/beam/wmbw_retained_join_agenda.hh"
#include "scfg2/beam/wmb_btail_store.hh"
#include "scfg2/beam/wmb_item_selection.hh"
#include "scfg2/beam/retained_choices.hh"
#include "scfg2/beam/wmb_retained_agenda.hh"
#include "scfg2/beam/wm_live_parent.hh"
#include "scfg2/beam/wm_retained_join_agenda.hh"
#include "scfg2/beam/wi_retained_join_agenda.hh"
#include "scfg2/beam/wip_retained_join_agenda.hh"
#include "scfg2/beam/wmv_wmp_live_parent.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/engine/choice_map.hh"
#include "scfg2/engine/generated_target_index.hh"
#include "scfg2/engine/traceback.hh"
#include "scfg2/engine/traceback_pairs.hh"

#include <array>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cparty::scfg2::replay {

namespace detail {
template <class Observer, class Item, class Chart>
auto observe_completed_item(Observer &observer, const Item &item,
                            const Chart &chart, int)
    -> decltype(observer.after_item(item, chart), void()) {
  observer.after_item(item, chart);
}
template <class Observer, class Item, class Chart>
void observe_completed_item(Observer &, const Item &, const Chart &, long) {}

template <class Observer, class Items, class Chart, class Store>
auto observe_be_commit(Observer &observer, const Items &items,
                       const Chart &chart, const Store &store, int)
    -> decltype(observer.after_be_commit(items, chart, store), void()) {
  observer.after_be_commit(items, chart, store);
}
template <class Observer, class Items, class Chart, class Store>
void observe_be_commit(Observer &, const Items &, const Chart &, const Store &, long) {}

template <class Observer, class Starts, class Store>
auto observe_wmbp_column(Observer &observer, cand_pos_t j,
                         const Starts &starts, const Store &store, int)
    -> decltype(observer.before_wmbp_column(j, starts, store), void()) {
  observer.before_wmbp_column(j, starts, store);
}
template <class Observer, class Starts, class Store>
void observe_wmbp_column(Observer &, cand_pos_t, const Starts &, const Store &, long) {}
}  // namespace detail

inline constexpr std::size_t w_final_online_nonterminal_count =
    static_cast<std::size_t>(grammar::NonTerminal::BE) + 1;
using WFinalOnlineNonterminalCounts =
    std::array<std::size_t, w_final_online_nonterminal_count>;

// Shared rank hook for numeric, graph, targeted and Viterbi consumers.
// Retained outside reconstructs WMB/BE views with the same fixed hook.
using WFinalPublicItemRanking = beam::PublicItemRanking;

struct WFinalVpWmbpOnlineWidths {
  std::size_t vp = std::numeric_limits<std::size_t>::max();
  std::size_t wmbp = std::numeric_limits<std::size_t>::max();
  std::size_t case3_producer = std::numeric_limits<std::size_t>::max();
  std::size_t case3_tail = std::numeric_limits<std::size_t>::max();
  std::size_t case3_left = std::numeric_limits<std::size_t>::max();
  std::size_t case3_active = std::numeric_limits<std::size_t>::max();
  std::size_t be_right_tail = std::numeric_limits<std::size_t>::max();
  std::size_t be_right_shell = std::numeric_limits<std::size_t>::max();
  std::size_t wmbp_prefix = std::numeric_limits<std::size_t>::max();
  std::size_t vm = std::numeric_limits<std::size_t>::max();
  std::size_t v = std::numeric_limits<std::size_t>::max();
  std::size_t wmv = std::numeric_limits<std::size_t>::max();
  std::size_t wmb = std::numeric_limits<std::size_t>::max();
  std::size_t wmp = std::numeric_limits<std::size_t>::max();
  std::size_t wm = std::numeric_limits<std::size_t>::max();
  std::size_t wi = std::numeric_limits<std::size_t>::max();
  std::size_t wip = std::numeric_limits<std::size_t>::max();
  std::size_t wmbw = std::numeric_limits<std::size_t>::max();
  std::size_t vpl = std::numeric_limits<std::size_t>::max();
  std::size_t vpr = std::numeric_limits<std::size_t>::max();
  std::size_t w = std::numeric_limits<std::size_t>::max();
  std::size_t be = std::numeric_limits<std::size_t>::max();
  std::size_t wmb_btail = std::numeric_limits<std::size_t>::max();
  std::size_t wmbp_wi_history = std::numeric_limits<std::size_t>::max();
  std::size_t wmbp_wi_prefix = std::numeric_limits<std::size_t>::max();
  std::size_t case3_nonleft_prefix = std::numeric_limits<std::size_t>::max();
  std::size_t case12_prefix = std::numeric_limits<std::size_t>::max();
  // Finite-beam seed policy only. KeepAll remains unlimited and therefore
  // exactly matches legacy CParty, whose ViennaRNA energy model extrapolates
  // hairpins longer than 30 rather than rejecting them.
  std::size_t max_hairpin_unpaired =
      std::numeric_limits<std::size_t>::max();
  // Independent finite-geometry horizon. This must not vary with beam width:
  // doing so changes the candidate universe between b and b+1 and can break
  // monotonicity even when each item frontier itself is top-b. It counts only
  // branching scaffold ancestors (at least two pair children), the only
  // possible CASE3 LCAs, so it is exact whenever no endpoint has more.
  std::size_t max_case3_ancestor_depth =
      std::numeric_limits<std::size_t>::max();
  // Items of one scaffold-only derivation (G' empty).  The V, VM, WM, WMv
  // and WMp frontiers always retain them, outside the width, so every
  // scaffold pair stays representable and the root stays reachable at any
  // width.  The runner fills this when it is null and
  // force_scaffold_derivation is set (see scaffold_only_derivation_items).
  std::shared_ptr<const beam::ItemKeySet> forced_items;
  bool force_scaffold_derivation = true;
  beam::MemoryObserver memory_observer = {};
  std::size_t memory_column_interval = 0;
  beam::MemoryPhase memory_phase = beam::MemoryPhase::Main;
  // Numeric lookup acceleration only; the canonical chart and its iteration
  // order remain authoritative for every consumer and retention decision.
  bool numeric_lookup_index = false;

  static WFinalVpWmbpOnlineWidths vp_wmbp(std::size_t width) {
    WFinalVpWmbpOnlineWidths out;
    out.vp = width;
    out.wmbp = width;
    return out;
  }

  static WFinalVpWmbpOnlineWidths uniform(std::size_t width) {
    auto out = WFinalVpWmbpOnlineWidths{
        width, width, width, width, width, width, width, width, width,
        width, width, width, width, width, width, width, width, width, width,
        width, width, width, width, width, width};
    out.case3_nonleft_prefix = width;
    out.case12_prefix = width;
    out.max_hairpin_unpaired = 30;
    out.max_case3_ancestor_depth = 30;
    return out;
  }
};

template <class Value>
struct WFinalVpWmbpOnlineRun {
  engine::Chart<Value> chart;
  beam::VpWmbpOnlineScheduleStats schedule;
  std::size_t vp_candidates = 0;
  std::size_t vp_selected = 0;
  std::size_t wmbp_candidates = 0;
  std::size_t wmbp_selected = 0;
  std::size_t case3_producer_candidates = 0;
  std::size_t case3_producer_selected = 0;
  std::size_t case3_tail_candidates = 0;
  std::size_t case3_tail_selected = 0;
  std::size_t case3_left_candidates = 0;
  std::size_t case3_left_selected = 0;
  std::size_t be_right_tail_candidates = 0;
  std::size_t be_right_tail_selected = 0;
  std::size_t be_right_shell_candidates = 0;
  std::size_t be_right_shell_selected = 0;
  std::size_t wmbp_prefix_candidates = 0;
  std::size_t wmbp_prefix_selected = 0;
  WFinalOnlineNonterminalCounts item_visits{};
  WFinalOnlineNonterminalCounts positive_items{};
  WFinalOnlineNonterminalCounts frontier_candidates{};
  WFinalOnlineNonterminalCounts frontier_selected{};
  std::size_t vm_binary_join_probes = 0;
  std::size_t vm_base_prefix_cell_visits = 0;
  std::size_t vm_base_prefix_selected = 0;
  std::size_t wm_binary_join_probes = 0;
  std::size_t wm_extension_probes = 0;
  std::size_t wm_start_prefix_cell_visits = 0;
  std::size_t wm_start_prefix_selected = 0;
  std::size_t wi_binary_join_probes = 0;
  std::size_t wi_direct_probes = 0;
  std::size_t wi_extension_probes = 0;
  std::size_t wip_binary_join_probes = 0;
  std::size_t wip_extension_probes = 0;
  std::size_t wip_start_prefix_cell_visits = 0;
  std::size_t wip_start_prefix_selected = 0;
  std::size_t wmbw_binary_join_probes = 0;
  std::size_t vpl_prefix_cell_visits = 0;
  std::size_t vpl_prefix_selected = 0;
  std::size_t vpr_binary_join_probes = 0;
  std::size_t vpr_base_state_probes = 0;
  std::size_t vpr_base_state_selected = 0;
  std::size_t w_binary_join_probes = 0;
  std::size_t w_direct_probes = 0;
  std::size_t w_extension_probes = 0;
  std::size_t w_root_decomposition_probes = 0;
  std::size_t wmb_split_positive_sources = 0;
  std::size_t wmb_direct_positive_sources = 0;
  std::size_t wmb_empty_positive_sources = 0;
  std::size_t wmb_btail_join_probes = 0;
  std::size_t wmb_btail_direct_probes = 0;
  std::size_t wmb_btail_selected = 0;
  std::size_t wmb_be_btail_join_probes = 0;
  std::size_t wmb_retained_direct_probes = 0;
  std::size_t wmb_retained_split_selected = 0;
  std::size_t wmbp_shared_wmbp_left_probes = 0;
  std::size_t wmbp_shared_wmbw_left_probes = 0;
  std::size_t wmbp_shared_emitted_sources = 0;
  std::size_t wmbp_wi_history_insertions = 0;
  std::size_t wmbp_wi_be_wi_probes = 0;
  std::size_t wmbp_wi_empty_probes = 0;
  std::size_t wmbp_wi_prefix_candidates = 0;
  std::size_t wmbp_wi_prefix_vp_probes = 0;
  std::size_t wmbp_wi_direct_probes = 0;
  std::size_t wmbp_wi_emitted_sources = 0;
  std::size_t vp_loop_parent_probes = 0;
  std::size_t vp_stack_event_sources = 0;
  std::size_t vp_internal_event_sources = 0;
  std::size_t vp_case3_nonleft_wi_queries = 0;
  std::size_t vp_case3_nonleft_geometry_probes = 0;
  std::size_t vp_case3_nonleft_right_last_probes = 0;
  std::size_t vp_case3_nonleft_middle_last_probes = 0;
  std::size_t vp_case3_nonleft_zero_child_probes = 0;
  std::size_t vp_case3_nonleft_emitted_sources = 0;
  std::size_t vp_case3_nonleft_border_mismatches = 0;
  std::size_t vp_case3_nonleft_duplicate_sources = 0;
  std::size_t vp_case3_nonleft_prefix_candidates = 0;
  std::size_t vp_case3_nonleft_active_prefix_probes = 0;
  std::size_t vp_case3_nonleft_max_active_prefixes = 0;
  std::size_t vp_case12_emitted_sources = 0;
  std::size_t vp_case12_border_mismatches = 0;
  std::size_t vp_case12_duplicate_sources = 0;
  std::array<std::size_t, 4> vp_recursive_split_join_probes{};
  std::array<std::size_t, 4> vp_recursive_split_emitted_sources{};
  std::size_t vp_recursive_split_duplicate_sources = 0;
  std::size_t be_stack_sources = 0;
  std::size_t be_base_stack_duplicate_sources = 0;
  std::size_t be_core_left_sources = 0;
  std::size_t be_core_left_duplicate_sources = 0;
  std::size_t be_right_tail_duplicate_registrations = 0;
  std::size_t be_right_shell_sources = 0;
  std::size_t be_right_shell_duplicate_sources = 0;
  std::size_t case3_active_expiry_probes = 0;
  std::size_t case3_tail_materialization_probes = 0;
  std::size_t case3_prepared_score_probes = 0;
  std::size_t case3_max_active_segments = 0;
  std::size_t v_hairpin_parent_probes = 0;
  std::size_t v_internal_parent_probes = 0;
  std::size_t v_vm_parent_probes = 0;
  std::size_t v_hairpin_sources = 0;
  std::size_t v_internal_sources = 0;
  std::size_t v_vm_sources = 0;
  std::size_t case3_heap_registered_segments = 0;
  std::size_t case3_heap_expired_segments = 0;
  std::size_t case3_max_heap_size = 0;
  std::size_t case3_endpoint_prefix_registrations = 0;
  std::size_t case3_endpoint_prefix_late_registrations = 0;
  std::size_t case3_endpoint_prefix_nonmonotone_queries = 0;
  std::size_t case3_endpoint_prefix_maximum_size = 0;
  std::size_t case3_endpoint_prefix_maximum_returned_per_path = 0;
  // Unpaired endpoints whose branching-ancestor path was cut by
  // max_case3_ancestor_depth. Zero means that cap removed no CASE3 source.
  std::size_t case3_truncated_endpoints = 0;
};

template <class Semiring,
          class Context,
          class CandidateObserverFn,
          class DeductionValueTransformFn>
WFinalVpWmbpOnlineRun<typename Semiring::value_type>
run_w_final_vp_wmbp_online_impl(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths,
    CandidateObserverFn &&candidate_observer,
    DeductionValueTransformFn &&deduction_value_transform,
    const WFinalPublicItemRanking &public_item_ranking = {},
    bool private_scaffold_root_relevant_only = false) {
  if (widths.memory_observer) {
    beam::MemoryObservation start;
    start.phase = widths.memory_phase;
    start.event = "begin";
    start.width = widths.v;
    if (widths.forced_items)
      start.containers.push_back(beam::hash_memory("forced_keys", *widths.forced_items));
    widths.memory_observer(start);
  }
  engine::AllLocalWeightModel<Context> weight(&ctx);
  constexpr bool numeric_lookup = std::is_same_v<Semiring, semirings::LogSumExpSemiring> ||
                                  std::is_same_v<Semiring, semirings::LogMaxSemiring>;
  // Private opt-in used only with the scaffold-only generated-pair transform
  // below. It is invalid for an unrestricted deduction-value transform.
  // Keep every root-relevant family's complete finite competition. The
  // ordinary kernel, custom Context and observed paths retain every family.
  const bool omit_scaffold_dead_families = [&] {
    if constexpr (std::is_same_v<Semiring, semirings::LogMaxSemiring> &&
                  std::is_same_v<std::decay_t<Context>, WFinalLogWeightContext>) {
      return private_scaffold_root_relevant_only && widths.numeric_lookup_index &&
          !widths.memory_observer && !public_item_ranking && n == ctx.sequence_length() &&
          beam::run_term_observer<double>() == nullptr;
    } else return false;
  }();
  engine::Chart<typename Semiring::value_type> chart(
      Semiring::zero(), [](const typename Semiring::value_type &value) {
        return Semiring::is_zero(value);
      }, 0, numeric_lookup && widths.numeric_lookup_index ? n : 0);
  // Seal only the standard, unobserved numeric mirror at the existing
  // complete-family commit boundary. Candidate/weight observations are intact.
  const bool interval_lookup = [&] {
    if constexpr (numeric_lookup &&
                  std::is_same_v<std::decay_t<Context>, WFinalLogWeightContext>)
      return widths.numeric_lookup_index && !widths.memory_observer &&
          n == ctx.sequence_length();
    else return false;
  }();
  const auto seal_primary_lookup_row = [&](grammar::NonTerminal family, cand_pos_t end) {
    if (interval_lookup) chart.seal_primary_lookup_row(family, end);
  };
  beam::TopWidthVpItemSelection vp_selection(widths.vp);
  beam::TopWidthWmbpItemSelection wmbp_selection(widths.wmbp);
  beam::TopWidthBeItemSelection be_selection(widths.be);
  const auto case3_child_index =
      beam::ExactVpCase3BorderEventChildIndex::build(
          n, ctx, widths.max_case3_ancestor_depth);
  const auto case3_nonleft_index =
      beam::VpCase3NonleftChildEventIndex::build(
          n, ctx, widths.max_case3_ancestor_depth);
  beam::VpCase3NonleftPrefixIntervalAgenda case3_nonleft_prefixes(
      widths.case3_nonleft_prefix);
  const auto case12_nonempty_query =
      beam::VpCase12DirectNonemptyLeftQuery::build(n, ctx);
  const auto case12_empty_query =
      beam::VpCase12DirectEmptyLeftQuery::build(
          n, ctx, widths.case12_prefix);
  beam::VpCase12PrefixIntervalAgenda case12_prefixes(widths.case12_prefix);
  using value_type = typename Semiring::value_type;
  beam::BasicPendingVpCase3LeftCandidates<value_type> case3_pending(
      widths.case3_active);
  beam::BasicExactVpCase3TailCache<value_type> case3_tails;
  beam::BasicExactVpCase3LeftCache<value_type> case3_left;
  beam::TopWidthVpCase3TailSelection case3_tail_selection(widths.case3_tail);
  beam::TopWidthVpCase3LeftSelection case3_left_selection(widths.case3_left);
  beam::VpCase3WiProducerFrontier case3_producer(
      widths.case3_producer,
      beam::VpCase3WiFrontierScope::EligibleCase3Producers);
  beam::TopWidthVpCase3PrecedingSiblingProducerFrontier
      case3_endpoint_prefix(widths.case3_active);
  auto be_tail_index =
      beam::BasicBeRightTailRetainedCoreIndex<value_type>::build(n, ctx);
  beam::TopWidthBeRightTailSelection be_tail_selection(
      widths.be_right_tail);
  beam::TopWidthBeRightShellFrameSelection be_shell_selection(
      widths.be_right_shell);
  beam::BeRightShellSourceAgenda be_right_shell_sources;
  beam::TopWidthWmbpBeVpPrefixSelection wmbp_prefix_selection(
      widths.wmbp_prefix);
  beam::TopWidthFamilyItemSelection vm_selection(
      grammar::NonTerminal::VM, widths.vm);
  beam::TopWidthFamilyItemSelection v_selection(
      grammar::NonTerminal::V, widths.v);
  beam::TopWidthFamilyItemSelection wmv_selection(
      grammar::NonTerminal::WMv, widths.wmv);
  beam::TopWidthWmbItemSelection wmb_selection(widths.wmb);
  beam::TopWidthFamilyItemSelection wmp_selection(
      grammar::NonTerminal::WMp, widths.wmp);
  beam::TopWidthFamilyItemSelection wm_selection(
      grammar::NonTerminal::WM, widths.wm);
  beam::TopWidthFamilyItemSelection wi_selection(
      grammar::NonTerminal::WI, widths.wi);
  beam::TopWidthFamilyItemSelection wip_selection(
      grammar::NonTerminal::WIP, widths.wip);
  beam::TopWidthFamilyItemSelection wmbw_selection(
      grammar::NonTerminal::WMBW, widths.wmbw);
  beam::TopWidthFamilyItemSelection vpl_selection(
      grammar::NonTerminal::VPL, widths.vpl);
  beam::TopWidthFamilyItemSelection vpr_selection(
      grammar::NonTerminal::VPR, widths.vpr);
  beam::TopWidthWItemSelection w_selection(widths.w);
  // Private numeric canonical frontiers change no retained rank order. The
  // only pruned callbacks below erase zero-valued chart entries; their order
  // cannot change the remaining chart iteration or any semiring arithmetic.
  // Memory observations retain the original container layout.
  if constexpr (numeric_lookup && std::is_same_v<std::decay_t<Context>, WFinalLogWeightContext>) {
    if (widths.numeric_lookup_index && !widths.memory_observer && n == ctx.sequence_length()) {
      vp_selection.enable_numeric_canonical_storage(n);
      wmbp_selection.enable_numeric_canonical_storage(n);
      vm_selection.enable_numeric_canonical_storage(n);
      v_selection.enable_numeric_canonical_storage(n);
      wmv_selection.enable_numeric_canonical_storage(n);
      wmb_selection.enable_numeric_canonical_storage(n);
      wmp_selection.enable_numeric_canonical_storage(n);
      wm_selection.enable_numeric_canonical_storage(n);
      wi_selection.enable_numeric_canonical_storage(n);
      wip_selection.enable_numeric_canonical_storage(n);
      wmbw_selection.enable_numeric_canonical_storage(n);
      vpl_selection.enable_numeric_canonical_storage(n);
      vpr_selection.enable_numeric_canonical_storage(n);
      w_selection.enable_numeric_canonical_storage(n);
    }
  }
  // Forced WM parents by right boundary, for the WM start-prefix selection.
  std::unordered_map<cand_pos_t, std::vector<cand_pos_t>> forced_wm_starts;
  if (widths.forced_items) {
    for (const auto &item : *widths.forced_items) {
      if (item.nonterminal == grammar::NonTerminal::WM) {
        forced_wm_starts[item.j].push_back(item.i);
      }
    }
    v_selection.set_forced(widths.forced_items);
    vm_selection.set_forced(widths.forced_items);
    wm_selection.set_forced(widths.forced_items);
    wmv_selection.set_forced(widths.forced_items);
    wmp_selection.set_forced(widths.forced_items);
  }
  beam::VmRetainedJoinAgenda<Semiring> vm_agenda;
  beam::VRetainedAgenda<Context> v_agenda;
  if constexpr (numeric_lookup && std::is_same_v<std::decay_t<Context>, WFinalLogWeightContext>) {
    if (widths.numeric_lookup_index && !widths.memory_observer && n == ctx.sequence_length())
      v_agenda.enable_numeric_pending(n);
  }

  beam::WmRetainedJoinAgenda<Semiring> wm_agenda;
  beam::WiRetainedJoinAgenda wi_agenda;
  beam::WipRetainedJoinAgenda<Semiring> wip_agenda;
  beam::WmbwRetainedJoinAgenda wmbw_agenda;
  beam::VplRetainedAgenda<Semiring> vpl_agenda;
  beam::VprRetainedAgenda<Semiring> vpr_agenda;
  beam::WRetainedJoinAgenda w_agenda;
  beam::WmbBtailStore<Semiring> wmb_btails;
  beam::WmbRetainedAgenda<Semiring> wmb_agenda;
  beam::WmbpSharedRetainedJoin<Context> wmbp_shared_join;
  beam::WmbpWiPrefixStore<Semiring> wmbp_wi_prefixes;
  beam::WmbpWiRetainedJoin<Semiring, Context> wmbp_wi_join;
  beam::VpLoopRetainedEvents<Context> vp_loop_events;
  beam::VpRecursiveSplitRetainedEvents vp_recursive_events;
  beam::BeBaseStackRetainedEvents be_base_stack_events;
  auto be_core_left_agenda =
      beam::BeCoreLeftPrefixAgenda::build(n, ctx);
  WFinalOnlineNonterminalCounts item_visits{};
  WFinalOnlineNonterminalCounts positive_items{};
  std::size_t wmb_split_positive_sources = 0;
  std::size_t wmb_direct_positive_sources = 0;
  std::size_t wmb_empty_positive_sources = 0;
  std::size_t wmbp_shared_wmbp_left_probes = 0;
  std::size_t wmbp_shared_wmbw_left_probes = 0;
  std::size_t wmbp_shared_emitted_sources = 0;
  std::size_t wmbp_wi_prefix_vp_probes = 0;
  std::size_t wmbp_wi_direct_probes = 0;
  std::size_t wmbp_wi_emitted_sources = 0;

  const auto observe_memory = [&](const char *event, cand_pos_t column,
                                  bool final = false) {
    if (!widths.memory_observer) return;
    auto observation = beam::chart_memory(widths.memory_phase, event,
        widths.v, column, chart, final);
    vp_loop_events.append_memory(observation.containers,
        "vp_loop_sources", "vp_loop_starts");
    observation.containers.push_back(wmb_btails.memory_snapshot("wmb_btail_history"));
    observation.containers.push_back(v_agenda.memory_snapshot("v_loop_sources"));
    observation.containers.push_back(case12_prefixes.memory_snapshot("case12_prefixes"));
    observation.containers.push_back(case3_nonleft_prefixes.memory_snapshot("case3_nonleft_prefixes"));
    const auto selection_memory = [&](const char *name, const auto &selection) {
      observation.containers.push_back(selection.memory_snapshot(name));
    };
    selection_memory("VM_selection", vm_selection);
    selection_memory("V_selection", v_selection);
    selection_memory("WMv_selection", wmv_selection);
    selection_memory("WMp_selection", wmp_selection);
    selection_memory("WM_selection", wm_selection);
    selection_memory("WI_selection", wi_selection);
    selection_memory("WIP_selection", wip_selection);
    selection_memory("WMBW_selection", wmbw_selection);
    selection_memory("VPL_selection", vpl_selection);
    selection_memory("VPR_selection", vpr_selection);
    widths.memory_observer(observation);
  };
  observe_memory("initialized", 0);

  auto provider = [&](const deductions::ItemKey &item) {
    if (omit_scaffold_dead_families) {
      using N = grammar::NonTerminal;
      switch (item.nonterminal) {
        case N::WI:
        case N::WIP:
        case N::BE:
        case N::VP:
        case N::WMBP:
        case N::WMBW:
        case N::VPL:
        case N::VPR:
        case N::WMB:
          // VP is exactly zero after the private generated-pair transform.
          // All WMBP modes need VP; nonempty WMB needs WMBP and WMB_EMPTY
          // is impossible in the standard Context. WMB is the only edge
          // from these families into W/VM/WMv/WMp/WM/V, so the entire set
          // is root-unconsumed in this pass, including for nonempty G.
          return std::vector<deductions::Deduction>{};
        default:
          break;
      }
    }
    if (item.nonterminal == grammar::NonTerminal::V &&
        v_agenda.column() == item.j) {
      return v_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::VP) {
      std::vector<deductions::Deduction> out;
      auto case12_sources = case12_prefixes.take_for(item);
      out.insert(out.end(),
                 std::make_move_iterator(case12_sources.begin()),
                 std::make_move_iterator(case12_sources.end()));
      auto case3_sources = deductions::generate_exact_basic_deductions(
          item, grammar::RuleId::VP_WI_CASE3, ctx);
      for (auto &source : case3_sources) {
        if (beam::vp_case3_partition(source) ==
            beam::VpCase3Partition::LeftNonempty) {
          out.push_back(std::move(source));
        }
      }
      auto nonleft_sources = case3_nonleft_prefixes.take_for(item);
      out.insert(out.end(),
                 std::make_move_iterator(nonleft_sources.begin()),
                 std::make_move_iterator(nonleft_sources.end()));
      auto loop_sources = vp_loop_events.deductions_for(item);
      std::sort(loop_sources.begin(), loop_sources.end(), [](const auto &lhs,
                                                             const auto &rhs) {
        if (lhs.rule != rhs.rule) {
          return lhs.rule == grammar::RuleId::VP_STACK;
        }
        if (lhs.split.k != rhs.split.k) return lhs.split.k < rhs.split.k;
        return lhs.split.l > rhs.split.l;
      });
      out.insert(out.end(), loop_sources.begin(), loop_sources.end());
      auto recursive_sources = vp_recursive_events.take_for(item);
      out.insert(out.end(),
                 std::make_move_iterator(recursive_sources.begin()),
                 std::make_move_iterator(recursive_sources.end()));
      return out;
    }
    if (item.nonterminal == grammar::NonTerminal::BE) {
      std::vector<deductions::Deduction> out;
      const auto &base_stack = be_base_stack_events.deductions_for(item);
      out.insert(out.end(), base_stack.begin(), base_stack.end());
      const auto &core_left = be_core_left_agenda.deductions_for(item);
      out.insert(out.end(), core_left.begin(), core_left.end());
      const auto &right_shell = be_right_shell_sources.deductions_for(item);
      out.insert(out.end(), right_shell.begin(), right_shell.end());
      return out;
    }
    if (item.nonterminal == grammar::NonTerminal::VM &&
        vm_agenda.column() == item.j) {
      // Each agenda parent has one executor visit. Transfer its sources while
      // leaving aggregate-contribution maps live for the fold below.
      return vm_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::WM &&
        wm_agenda.column() == item.j) {
      return wm_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::WI &&
        wi_agenda.column() == item.j) {
      return wi_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::WIP &&
        wip_agenda.column() == item.j) {
      return wip_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::WMBW &&
        wmbw_agenda.column() == item.j) {
      return wmbw_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::VPL &&
        vpl_agenda.column() == item.j) {
      return vpl_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::VPR &&
        vpr_agenda.column() == item.j) {
      return vpr_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::W &&
        w_agenda.column() == item.j) {
      return w_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::WMB &&
        wmb_agenda.column() == item.j) {
      return wmb_agenda.take_for(item);
    }
    if (item.nonterminal == grammar::NonTerminal::WMBP) {
      std::vector<deductions::Deduction> deductions;
      if (wmbp_shared_join.column() == item.j) {
        const auto &shared = wmbp_shared_join.deductions_for(item);
        deductions.insert(deductions.end(), shared.begin(), shared.end());
      }
      if (wmbp_wi_join.column() == item.j) {
        const auto &remaining = wmbp_wi_join.deductions_for(item);
        deductions.insert(
            deductions.end(), remaining.begin(), remaining.end());
      }
      auto scored =
          beam::score_wmbp_be_vp_prefixes_from_deductions<Semiring>(
          deductions, chart, ctx);
      wmbp_prefix_selection.template commit_parent<Semiring>(
          item, std::move(scored));
      return deductions;
    }
    if (item.nonterminal == grammar::NonTerminal::WMBW ||
        item.nonterminal == grammar::NonTerminal::WMB) {
      return deductions::
          generate_exact_wmbp_consumer_deductions_for_retained_ends(
              item, ctx, wmbp_selection.selected_ends(item.i));
    }
    return deductions::generate_exact_basic_deductions_for_item(item, ctx);
  };
  auto allows = [](const deductions::Deduction &) { return true; };
  auto before_item = [](const deductions::ItemKey &, const auto &) {};
  auto after_item = [&](const deductions::ItemKey &item,
                        const auto &current_chart) {
    const auto value = current_chart.get(item);
    const bool live = !Semiring::is_zero(value);
    if (live) ++positive_items[static_cast<std::size_t>(item.nonterminal)];
    const double priority = Semiring::priority(value);
    const bool external_view = item.nonterminal == grammar::NonTerminal::V ||
                               item.nonterminal == grammar::NonTerminal::WMB;
    double existing_priority = priority;
    if (external_view) {
      const auto prefix = item.i <= 1 ? Semiring::one()
          : current_chart.get(deductions::ItemKey{
                grammar::NonTerminal::W, 1, item.i - 1, -1, -1});
      existing_priority = Semiring::priority(Semiring::times(prefix, value));
    }
    const double ranked_priority = live && public_item_ranking
        ? public_item_ranking(item, priority, existing_priority)
        : existing_priority;
    if (live && public_item_ranking && (std::isnan(ranked_priority) ||
        ranked_priority == std::numeric_limits<double>::infinity()))
      throw std::logic_error("invalid public-item ranking priority");
    vp_selection.observe_ranked(item, live, ranked_priority);
    wmbp_selection.observe_ranked(item, live, ranked_priority);
    be_selection.observe_ranked(item, live, ranked_priority);
    vm_selection.observe_ranked(item, live, ranked_priority);
    v_selection.observe_ranked(item, live, ranked_priority);
    wmv_selection.observe_ranked(item, live, ranked_priority);
    wmb_selection.observe_ranked(
        item, live, priority, ranked_priority);
    wmp_selection.observe_ranked(item, live, ranked_priority);
    wm_selection.observe_ranked(item, live, ranked_priority);
    wi_selection.observe_ranked(item, live, ranked_priority);
    wip_selection.observe_ranked(item, live, ranked_priority);
    wmbw_selection.observe_ranked(item, live, ranked_priority);
    vpl_selection.observe_ranked(item, live, ranked_priority);
    vpr_selection.observe_ranked(item, live, ranked_priority);
    w_selection.observe_ranked(item, live, ranked_priority);
    // Observe candidate items before commit; pruned public items must never
    // be consumed by subsequent deductions.
    detail::observe_completed_item(candidate_observer, item, current_chart, 0);
    if (item.nonterminal == grammar::NonTerminal::WMBP)
      wmbp_prefix_selection.release_completed_parents();
  };
  auto raw_contribution = [&](const deductions::Deduction &deduction,
                              const auto &current_chart,
                              auto local) {
    // Selection and private-aggregate ownership are disjoint rule contracts.
    // Dispatch once instead of asking every unrelated owner about each source.
    // Local-weight evaluation and its zero gate remain in the shared executor.
    using Rule = grammar::RuleId;
    switch (deduction.rule) {
      case Rule::VP_WI_CASE3: {
        if (!case3_left_selection.allows(deduction)) return Semiring::zero();
        const auto reconstructed =
            case3_left.template reconstruct_source<Semiring>(deduction, local);
        if (reconstructed.has_value()) return *reconstructed;
        break;
      }
      case Rule::BE_WIP_WIP:
      case Rule::BE_BASEPAIR_WIP:
        if (!be_tail_selection.allows(deduction)) return Semiring::zero();
        if (!be_shell_selection.allows(deduction)) return Semiring::zero();
        break;
      case Rule::WMBP_SPLIT_BE_WMBP_VP:
      case Rule::WMBP_SPLIT_BE_WMBW_VP:
      case Rule::WMBP_SPLIT_BE_WI_VP:
        if (!wmbp_prefix_selection.allows(deduction)) return Semiring::zero();
        break;
      case Rule::VM_SPLIT_WMp_BASE:
        if (const auto *value = vm_agenda.find_base_contribution(deduction))
          return semirings::checkpoint_value(*value);
        break;
      case Rule::WM_START_V:
      case Rule::WM_START_WMB:
        if (const auto *value = wm_agenda.find_start_contribution(deduction))
          return semirings::checkpoint_value(*value);
        break;
      case Rule::WIP_BASE_V:
      case Rule::WIP_BASE_WMB:
      case Rule::WIP_BASEPAIR_V:
      case Rule::WIP_BASEPAIR_WMB:
        if (const auto *value = wip_agenda.find_start_contribution(deduction))
          return semirings::checkpoint_value(*value);
        break;
      case Rule::VPL_SPLIT_VP:
        if (const auto *value = vpl_agenda.find_contribution(deduction))
          return semirings::checkpoint_value(*value);
        break;
      case Rule::VPR_SPLIT_VP_BASEPAIR:
        if (const auto *value = vpr_agenda.find_base_contribution(deduction))
          return semirings::checkpoint_value(*value);
        break;
      case Rule::WMB_SPLIT_BE_WMBP_WI:
        if (const auto *value = wmb_agenda.find_split_contribution(deduction))
          return semirings::checkpoint_value(*value);
        break;
      default:
        break;
    }
    return runtime_deduction_contribution<Semiring>(
        current_chart, deduction, local);
  };
  const auto contribution = [&](const deductions::Deduction &deduction,
                                const auto &current_chart,
                                auto local) {
    return deduction_value_transform(
        deduction, raw_contribution(deduction, current_chart, local));
  };
  const auto update = [&](const deductions::ItemKey &item) {
    // W_ROOT_WMB_DECOMPOSITION reads a suffix W in this same column.
    // Starts are visited in descending order, so commit all suffixes before
    // folding the root. The root itself is retained separately (b + 1).
    if (item.nonterminal == grammar::NonTerminal::W && item.i == 1)
      w_selection.commit_nonroot_column(item.j, chart, Semiring::zero());
    if (!deductions::exact_basic_item_parent_active(item, ctx)) return;
    if (!beam::wmv_wmp_parent_has_live_child<Semiring>(item, chart, ctx))
      return;
    // The WM / VM live-child checks scan O(span) per parent.  They pay off
    // only when every start is probed (KeepAll widths).  Finite agendas only
    // offer parents with emitted deductions; a parent whose children are all
    // zero then simply evaluates to zero, which no selection treats as live,
    // so skipping the checks changes no value.
    if (widths.wm == std::numeric_limits<std::size_t>::max() &&
        !beam::wm_parent_has_live_child<Semiring>(item, chart, ctx)) {
      return;
    }
    if (widths.vm == std::numeric_limits<std::size_t>::max() &&
        !beam::vm_parent_has_live_child<Semiring>(item, chart, ctx)) {
      return;
    }
    ++item_visits[static_cast<std::size_t>(item.nonterminal)];
    const auto observed_candidate = [&](
        const deductions::ItemKey &observed_item,
        const deductions::Deduction &deduction,
        typename Semiring::value_type candidate,
        const auto &observed_chart,
        auto local) {
      if (observed_item.nonterminal == grammar::NonTerminal::WMB &&
          !Semiring::is_zero(candidate)) {
        switch (deduction.rule) {
          case grammar::RuleId::WMB_SPLIT_BE_WMBP_WI:
            ++wmb_split_positive_sources;
            break;
          case grammar::RuleId::WMB_DIRECT_WMBP:
            ++wmb_direct_positive_sources;
            break;
          case grammar::RuleId::WMB_EMPTY:
            ++wmb_empty_positive_sources;
            break;
          default:
            break;
        }
      }
      candidate_observer(
          observed_item, deduction, candidate, observed_chart, local);
    };
    update_w_final_exact_dp_item_with_provider_and_contribution<Semiring>(
        item,
        chart,
        weight,
        provider,
        allows,
        observed_candidate,
        before_item,
        after_item,
        contribution);
  };

  const std::vector<cand_pos_t> no_retained_starts;
  const auto schedule = beam::execute_vp_wmbp_online_schedule(
      n,
      [&](cand_pos_t j) {
        // Only this column's CASE3 products are consumed by VP.
        // Future columns regenerate from immutable retained WI.
        // Unselected tails must not accumulate into O(nb^2) history.
        case3_pending.release_target_column(j - 1);
        case3_tails = {};
        case3_left = {};
        be_base_stack_events.prepare_column(j, ctx);
        be_core_left_agenda.prepare_column(j, ctx);
        vp_recursive_events.prepare_column(
            j,
            [&](cand_pos_t right) -> const std::vector<cand_pos_t> & {
              return vp_selection.selected_starts(right);
            },
            [&](cand_pos_t right) -> const std::vector<cand_pos_t> & {
              return wip_selection.retained_starts(right);
            },
            [&](cand_pos_t right) -> const std::vector<cand_pos_t> & {
              return vpl_selection.retained_starts(right);
            },
            [&](cand_pos_t right) -> const std::vector<cand_pos_t> & {
              return vpr_selection.retained_starts(right);
            },
            ctx);
        case12_prefixes.template prepare_column<Semiring>(
            j,
            j > 1
                ? wi_selection.retained_starts(j - 1)
                : no_retained_starts,
            case12_empty_query,
            ctx);
        case3_nonleft_prefixes.template prepare_column<Semiring>(
            j, case3_nonleft_index, ctx);
        const auto be_tails = be_tail_index.take_column(j);
        be_tail_selection.template commit_column<Semiring>(j, be_tails);
        const auto be_frames =
            beam::score_be_right_shell_frames_from_tail_events<Semiring>(
                be_tails,
                chart,
                ctx,
                [&](const auto &key) {
                  return be_tail_selection.allows(key);
                });
        be_shell_selection.template commit_column<Semiring>(j, be_frames);
        be_right_shell_sources.prepare_column(
            be_frames, be_shell_selection, ctx);
        case3_producer.template flush_before_registering<Semiring>(
            j,
            chart,
            case3_endpoint_prefix,
            case3_child_index,
            ctx);
        const auto case3_endpoint_matches =
            case3_endpoint_prefix.query_target_endpoint(
                j, case3_child_index, ctx);
        const auto case3_endpoint_query =
            beam::build_vp_case3_sources_from_endpoint_matches(
                j, case3_endpoint_matches, ctx);
        case3_pending.template observe_target_query<Semiring>(
            case3_endpoint_query, chart);
        case3_pending.template prepare_column<Semiring>(
            j, chart, case3_tails);
        case3_tail_selection.template commit_column<Semiring>(
            j, case3_tails);
        const auto scored =
            case3_pending.template score_prepared_column<Semiring>(
            j,
            case3_tails,
            [&](const auto &key) { return case3_tail_selection.allows(key); });
        case3_left_selection.template commit_column<Semiring>(j, scored);
        (void)case3_left.materialize_scored(
            scored,
            [&](const auto &key) { return case3_left_selection.allows(key); });
      },
      [&](cand_pos_t j) {
        std::vector<cand_pos_t> starts;
        const auto append = [&](const std::vector<cand_pos_t> &more) {
          starts.insert(starts.end(), more.begin(), more.end());
        };
        append(case12_prefixes.parent_starts());
        append(case3_nonleft_prefixes.parent_starts());
        append(case3_left.vp_parent_starts(j));
        append(vp_loop_events.parent_starts(j));
        append(vp_recursive_events.parent_starts());
        std::sort(starts.begin(), starts.end(), std::greater<cand_pos_t>());
        starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
        return starts;
      },
      [&](cand_pos_t j) {
        (void)j;
        std::vector<deductions::ItemKey> parents;
        const auto before = [](const auto &lhs, const auto &rhs) {
          return std::tie(lhs.i, lhs.ip, lhs.jp) > std::tie(rhs.i, rhs.ip, rhs.jp);
        };
        const auto merge = [&](const auto &more) {
          std::vector<deductions::ItemKey> merged;
          merged.reserve(parents.size() + more.size());
          std::merge(parents.begin(), parents.end(), more.begin(), more.end(),
                     std::back_inserter(merged), before);
          parents = std::move(merged);
        };
        merge(be_base_stack_events.parents());
        merge(be_core_left_agenda.parents());
        merge(be_right_shell_sources.parents());
        parents.erase(
            std::unique(parents.begin(), parents.end(), [](const auto &lhs,
                                                           const auto &rhs) {
              return deductions::ItemKeyEq{}(lhs, rhs);
            }),
            parents.end());
        return parents;
      },
      update,
      [&](cand_pos_t j) {
        vp_selection.commit_column(j, chart, Semiring::zero());
        seal_primary_lookup_row(grammar::NonTerminal::VP, j);
        vp_loop_events.observe_committed_column(
            j, vp_selection.selected_starts(j), n, ctx);
        vp_loop_events.release_column(j);
      },
      [&](cand_pos_t j) {
        be_selection.commit_column(j, chart, Semiring::zero());
        be_base_stack_events.observe_committed_column(
            be_selection.retained_items(j), n, ctx);
        be_core_left_agenda.observe_committed_be_column(
            be_selection.retained_items(j), ctx);
        be_tail_index.observe_retained_be_column(
            be_selection.retained_items(j), ctx);
        wmbp_wi_prefixes.commit_be_column(
            be_selection.retained_items(j),
            chart,
            widths.wmbp_wi_prefix);
        detail::observe_be_commit(candidate_observer,
            be_selection.retained_items(j), chart, wmbp_wi_prefixes, 0);
      },
      [&](cand_pos_t j) {
        detail::observe_wmbp_column(candidate_observer, j,
            vp_selection.selected_starts(j), wmbp_wi_prefixes, 0);
        wmbp_shared_join.build_column(
            j,
            vp_selection.selected_starts(j),
            [&](grammar::NonTerminal family,
                cand_pos_t right_boundary)
                -> const std::vector<cand_pos_t> & {
              return family == grammar::NonTerminal::WMBP
                  ? wmbp_selection.selected_starts(right_boundary)
                  : wmbw_selection.retained_starts(right_boundary);
            },
            ctx);
        wmbp_shared_wmbp_left_probes +=
            wmbp_shared_join.stats().wmbp_left_probes;
        wmbp_shared_wmbw_left_probes +=
            wmbp_shared_join.stats().wmbw_left_probes;
        wmbp_shared_emitted_sources +=
            wmbp_shared_join.stats().emitted_sources;
        wmbp_wi_join.build_column(
            j, vp_selection.selected_starts(j), wmbp_wi_prefixes, ctx);
        wmbp_wi_prefix_vp_probes +=
            wmbp_wi_join.stats().prefix_vp_probes;
        wmbp_wi_direct_probes += wmbp_wi_join.stats().direct_probes;
        wmbp_wi_emitted_sources += wmbp_wi_join.stats().emitted_sources;
        beam::WmbpSourceParentKeyDiscovery remaining;
        remaining.parents = wmbp_wi_join.parents();
        remaining.parent_vp_probes =
            wmbp_wi_join.stats().prefix_vp_probes +
            wmbp_wi_join.stats().direct_probes;
        std::unordered_set<deductions::ItemKey,
                           deductions::ItemKeyHash,
                           deductions::ItemKeyEq>
            seen(remaining.parents.begin(), remaining.parents.end());
        for (const auto &parent : wmbp_shared_join.parents()) {
          if (seen.insert(parent).second) remaining.parents.push_back(parent);
        }
        std::sort(remaining.parents.begin(),
                  remaining.parents.end(),
                  [](const auto &lhs, const auto &rhs) {
                    return lhs.i > rhs.i;
                  });
        remaining.parent_vp_probes +=
            wmbp_shared_join.stats().wmbp_left_probes +
            wmbp_shared_join.stats().wmbw_left_probes;
        return remaining;
      },
      [&](cand_pos_t j) {
        wmbp_selection.commit_column(j, chart, Semiring::zero());
        seal_primary_lookup_row(grammar::NonTerminal::WMBP, j);
      },
      [&](grammar::NonTerminal nonterminal, cand_pos_t j)
          -> std::vector<cand_pos_t> {
        const auto retained_starts = [&](grammar::NonTerminal family,
                                         cand_pos_t right_boundary)
            -> const std::vector<cand_pos_t> & {
          switch (family) {
            case grammar::NonTerminal::V:
              return v_selection.retained_starts(right_boundary);
            case grammar::NonTerminal::WMv:
              return wmv_selection.retained_starts(right_boundary);
            case grammar::NonTerminal::WMB:
              return wmb_selection.internal_starts(right_boundary);
            case grammar::NonTerminal::WMp:
              return wmp_selection.retained_starts(right_boundary);
            case grammar::NonTerminal::WM:
              return wm_selection.retained_starts(right_boundary);
            case grammar::NonTerminal::WI:
              return wi_selection.retained_starts(right_boundary);
            case grammar::NonTerminal::WIP:
              return wip_selection.retained_starts(right_boundary);
            case grammar::NonTerminal::VPL:
              return vpl_selection.retained_starts(right_boundary);
            case grammar::NonTerminal::VPR:
              return vpr_selection.retained_starts(right_boundary);
            case grammar::NonTerminal::W:
              return w_selection.retained_starts(right_boundary);
            default: {
              static const std::vector<cand_pos_t> empty;
              return empty;
            }
          }
        };
        if (nonterminal == grammar::NonTerminal::VM) {
          vm_agenda.build_column(
              j, chart, ctx, retained_starts, widths.vm);
          return vm_agenda.parent_starts();
        }
        if (nonterminal == grammar::NonTerminal::V) {
          v_agenda.build_column(
              j,
              vm_selection.retained_starts(j),
              widths.max_hairpin_unpaired,
              ctx);
          return v_agenda.parent_starts();
        }
        if (nonterminal == grammar::NonTerminal::WM) {
          const auto forced_wm = forced_wm_starts.find(j);
          wm_agenda.build_column(
              j, chart, ctx, retained_starts, widths.wm,
              forced_wm == forced_wm_starts.end() ? nullptr : &forced_wm->second);
          return wm_agenda.parent_starts();
        }
        if (nonterminal == grammar::NonTerminal::WI) {
          wi_agenda.build_column(j, ctx, retained_starts);
          return wi_agenda.parent_starts();
        }
        if (nonterminal == grammar::NonTerminal::WIP) {
          wip_agenda.build_column(
              j, chart, ctx, retained_starts, widths.wip);
          return wip_agenda.parent_starts();
        }
        if (nonterminal == grammar::NonTerminal::WMBW) {
          wmb_btails.build(
              j + 1,
              chart,
              ctx,
              wi_selection.retained_starts(j),
              [&](cand_pos_t right_boundary)
                  -> const std::vector<cand_pos_t> & {
                return wmbp_selection.selected_starts(right_boundary);
              },
              widths.wmb_btail);
          wmbw_agenda.build_column(
              j,
              ctx,
              [&](cand_pos_t right_boundary)
                  -> const std::vector<cand_pos_t> & {
                return wi_selection.retained_starts(right_boundary);
              },
              [&](cand_pos_t right_boundary)
                  -> const std::vector<cand_pos_t> & {
                return wmbp_selection.selected_starts(right_boundary);
              });
          return wmbw_agenda.parent_starts();
        }
        if (nonterminal == grammar::NonTerminal::WMB) {
          wmb_agenda.build_column(
              j,
              chart,
              ctx,
              wmbp_selection.selected_starts(j),
              be_selection.retained_items(j),
              wmb_btails);
          return wmb_agenda.parent_starts();
        }
        if (nonterminal == grammar::NonTerminal::VPL) {
          vpl_agenda.build_column(
              j, chart, ctx, vp_selection.selected_starts(j), widths.vpl);
          return vpl_agenda.parent_starts();
        }
        if (nonterminal == grammar::NonTerminal::VPR) {
          vpr_agenda.build_column(
              j,
              chart,
              ctx,
              wip_selection.retained_starts(j),
              [&](cand_pos_t right_boundary)
                  -> const std::vector<cand_pos_t> & {
                return vp_selection.selected_starts(right_boundary);
              },
              widths.vpr);
          return vpr_agenda.parent_starts();
        }
        if (nonterminal == grammar::NonTerminal::W) {
          // W is the only root-prefix consumer of WMB.  Every other family
          // reads the context-free internal WMB view above.
          w_agenda.build_column(
              j,
              ctx,
              [&](grammar::NonTerminal family,
                  cand_pos_t right_boundary)
                  -> const std::vector<cand_pos_t> & {
                return family == grammar::NonTerminal::WMB
                    ? wmb_selection.root_suffix_starts(right_boundary)
                    : retained_starts(family, right_boundary);
              },
              [&](cand_pos_t i, cand_pos_t k) {
                return wmb_selection.root_suffix_contains(i, k);
              });
          return w_agenda.parent_starts();
        }
        const std::vector<cand_pos_t> *stem_starts = nullptr;
        const std::vector<cand_pos_t> *extension_starts = nullptr;
        if (nonterminal == grammar::NonTerminal::WMv) {
          stem_starts = &v_selection.retained_starts(j);
          extension_starts = &wmv_selection.retained_starts(j - 1);
        } else if (nonterminal == grammar::NonTerminal::WMp) {
          stem_starts = &wmb_selection.internal_starts(j);
          extension_starts = &wmp_selection.retained_starts(j - 1);
        } else {
          throw std::logic_error("missing sparse SCFG2 family agenda");
        }
        std::vector<cand_pos_t> starts;
        starts.reserve(stem_starts->size() + extension_starts->size());
        starts.insert(starts.end(), stem_starts->begin(), stem_starts->end());
        starts.insert(
            starts.end(), extension_starts->begin(), extension_starts->end());
        std::sort(starts.begin(), starts.end(), std::greater<cand_pos_t>());
        starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
        return starts;
      },
      [&](grammar::NonTerminal nonterminal, cand_pos_t j) {
        switch (nonterminal) {
          case grammar::NonTerminal::VM:
            vm_selection.commit_column(j, chart, Semiring::zero());
            break;
          case grammar::NonTerminal::V:
            v_selection.commit_column(j, chart, Semiring::zero());
            v_agenda.observe_committed_v_column(
                j, v_selection.retained_starts(j), n, ctx);
            v_agenda.release_column();
            break;
          case grammar::NonTerminal::WMv:
            wmv_selection.commit_column(j, chart, Semiring::zero());
            break;
          case grammar::NonTerminal::WMB:
            wmb_selection.commit_column(j, chart, Semiring::zero());
            break;
          case grammar::NonTerminal::WMp:
            wmp_selection.commit_column(j, chart, Semiring::zero());
            break;
          case grammar::NonTerminal::WM:
            wm_selection.commit_column(j, chart, Semiring::zero());
            break;
          case grammar::NonTerminal::WI:
            wi_selection.commit_column(j, chart, Semiring::zero());
            for (const cand_pos_t start :
                 wi_selection.retained_starts(j)) {
              const deductions::ItemKey middle{
                  grammar::NonTerminal::WI, start, j, -1, -1};
              case3_producer.template observe_completed_child<Semiring>(
                  middle, chart, case3_child_index, ctx);
              case12_prefixes.template observe_left_wi<Semiring>(
                  middle,
                  chart.get(middle),
                  case12_nonempty_query,
                  ctx);
              case3_nonleft_prefixes.template observe_middle_wi<Semiring>(
                  middle,
                  chart.get(middle),
                  case3_nonleft_index,
                  ctx);
            }
            wmbp_wi_prefixes.commit_wi_column(
                j,
                wi_selection.retained_starts(j),
                chart,
                widths.wmbp_wi_history);
            break;
          case grammar::NonTerminal::WIP:
            wip_selection.commit_column(j, chart, Semiring::zero());
            for (const cand_pos_t start : wip_selection.retained_starts(j)) {
              be_tail_index.template observe_retained_child<Semiring>(
                  deductions::ItemKey{
                      grammar::NonTerminal::WIP, start, j, -1, -1}, chart);
            }
            be_core_left_agenda.observe_retained_wip_column(
                wip_selection.retained_starts(j), j, ctx);
            break;
          case grammar::NonTerminal::WMBW:
            wmbw_selection.commit_column(j, chart, Semiring::zero());
            break;
          case grammar::NonTerminal::VPL:
            vpl_selection.commit_column(j, chart, Semiring::zero());
            break;
          case grammar::NonTerminal::VPR:
            vpr_selection.commit_column(j, chart, Semiring::zero());
            break;
          case grammar::NonTerminal::W:
            w_selection.commit_column(j, chart, Semiring::zero());
            break;
          default:
            break;
        }
        // W is sealed here after its root update and final commit, never
        // during the earlier nonroot commit used by the root's suffix reads.
        seal_primary_lookup_row(nonterminal, j);
        if (nonterminal == grammar::NonTerminal::W &&
            widths.memory_column_interval &&
            static_cast<std::size_t>(j) % widths.memory_column_interval == 0)
          observe_memory("column_end", j);
      });
  observe_memory("forward_end", n, true);
  WFinalOnlineNonterminalCounts frontier_candidates{};
  WFinalOnlineNonterminalCounts frontier_selected{};
  const auto record_frontier = [&](grammar::NonTerminal nonterminal,
                                   const auto &selection) {
    const auto idx = static_cast<std::size_t>(nonterminal);
    frontier_candidates[idx] = selection.candidate_count();
    frontier_selected[idx] = selection.selected_count();
  };
  record_frontier(grammar::NonTerminal::VM, vm_selection);
  record_frontier(grammar::NonTerminal::V, v_selection);
  record_frontier(grammar::NonTerminal::WMv, wmv_selection);
  record_frontier(grammar::NonTerminal::WMB, wmb_selection);
  record_frontier(grammar::NonTerminal::WMp, wmp_selection);
  record_frontier(grammar::NonTerminal::WM, wm_selection);
  record_frontier(grammar::NonTerminal::WI, wi_selection);
  record_frontier(grammar::NonTerminal::WIP, wip_selection);
  record_frontier(grammar::NonTerminal::WMBW, wmbw_selection);
  record_frontier(grammar::NonTerminal::VPL, vpl_selection);
  record_frontier(grammar::NonTerminal::VPR, vpr_selection);
  record_frontier(grammar::NonTerminal::W, w_selection);
  record_frontier(grammar::NonTerminal::BE, be_selection);
  return WFinalVpWmbpOnlineRun<typename Semiring::value_type>{
      std::move(chart),
      schedule,
      vp_selection.candidate_count(),
      vp_selection.selected_count(),
      wmbp_selection.candidate_count(),
      wmbp_selection.selected_count(),
      case3_producer.candidate_count(),
      case3_producer.selected_count(),
      case3_tail_selection.candidate_count(),
      case3_tail_selection.selected_count(),
      case3_left_selection.candidate_count(),
      case3_left_selection.selected_count(),
      be_tail_selection.candidate_count(),
      be_tail_selection.selected_count(),
      be_shell_selection.candidate_count(),
      be_shell_selection.selected_count(),
      wmbp_prefix_selection.candidate_count(),
      wmbp_prefix_selection.selected_count(),
      item_visits,
      positive_items,
      frontier_candidates,
      frontier_selected,
      vm_agenda.stats().binary_join_probes,
      vm_agenda.stats().base_prefix_cell_visits,
      vm_agenda.stats().base_prefix_selected,
      wm_agenda.stats().binary_join_probes,
      wm_agenda.stats().extension_probes,
      wm_agenda.stats().start_prefix_cell_visits,
      wm_agenda.stats().start_prefix_selected,
      wi_agenda.stats().binary_join_probes,
      wi_agenda.stats().direct_probes,
      wi_agenda.stats().extension_probes,
      wip_agenda.stats().binary_join_probes,
      wip_agenda.stats().extension_probes,
      wip_agenda.stats().start_prefix_cell_visits,
      wip_agenda.stats().start_prefix_selected,
      wmbw_agenda.stats().binary_join_probes,
      vpl_agenda.stats().prefix_cell_visits,
      vpl_agenda.stats().prefix_selected,
      vpr_agenda.stats().binary_join_probes,
      vpr_agenda.stats().base_state_probes,
      vpr_agenda.stats().base_state_selected,
      w_agenda.stats().binary_join_probes,
      w_agenda.stats().direct_probes,
      w_agenda.stats().extension_probes,
      w_agenda.stats().root_decomposition_probes,
      wmb_split_positive_sources,
      wmb_direct_positive_sources,
      wmb_empty_positive_sources,
      wmb_btails.stats().join_probes,
      wmb_btails.stats().direct_probes,
      wmb_btails.stats().selected,
      wmb_agenda.stats().be_btail_join_probes,
      wmb_agenda.stats().direct_probes,
      wmb_agenda.stats().split_selected,
      wmbp_shared_wmbp_left_probes,
      wmbp_shared_wmbw_left_probes,
      wmbp_shared_emitted_sources,
      wmbp_wi_prefixes.stats().wi_history_insertions,
      wmbp_wi_prefixes.stats().be_wi_probes,
      wmbp_wi_prefixes.stats().empty_probes,
      wmbp_wi_prefixes.stats().prefix_candidates,
      wmbp_wi_prefix_vp_probes,
      wmbp_wi_direct_probes,
      wmbp_wi_emitted_sources,
      vp_loop_events.stats().parent_probes,
      vp_loop_events.stats().stack_sources,
      vp_loop_events.stats().internal_sources,
      case3_nonleft_prefixes.stats().middle_queries,
      case3_nonleft_prefixes.stats().boundary_probes +
          case3_nonleft_prefixes.stats().active_prefix_probes,
      0,
      0,
      0,
      case3_nonleft_prefixes.stats().emitted_sources,
      case3_nonleft_prefixes.stats().border_mismatches,
      case3_nonleft_prefixes.stats().duplicate_sources,
      case3_nonleft_prefixes.stats().prefix_candidates,
      case3_nonleft_prefixes.stats().active_prefix_probes,
      case3_nonleft_prefixes.stats().max_active_prefixes,
      case12_prefixes.stats().emitted_sources,
      case12_prefixes.stats().border_mismatches,
      case12_prefixes.stats().duplicate_sources,
      vp_recursive_events.stats().join_probes,
      vp_recursive_events.stats().emitted_sources,
      vp_recursive_events.stats().duplicate_sources,
      be_base_stack_events.stats().stack_sources,
      be_base_stack_events.stats().duplicate_sources,
      be_core_left_agenda.stats().emitted_sources,
      be_core_left_agenda.stats().duplicate_sources,
      be_tail_index.stats().duplicate_registrations,
      be_right_shell_sources.stats().emitted_sources,
      be_right_shell_sources.stats().duplicate_sources,
      case3_pending.active_expiry_probe_count(),
      case3_pending.tail_materialization_probe_count(),
      case3_pending.prepared_score_probe_count(),
      case3_pending.max_active_segment_count(),
      v_agenda.stats().hairpin_parent_probes,
      v_agenda.stats().internal_parent_probes,
      v_agenda.stats().vm_parent_probes,
      v_agenda.stats().hairpin_sources,
      v_agenda.stats().internal_sources,
      v_agenda.stats().vm_sources,
      case3_pending.finite_registered_segment_count(),
      case3_pending.finite_expired_segment_count(),
      case3_pending.finite_max_heap_size(),
      case3_endpoint_prefix.candidate_registration_count(),
      case3_endpoint_prefix.late_registration_count(),
      case3_endpoint_prefix.nonmonotone_query_count(),
      case3_endpoint_prefix.maximum_prefix_size(),
      case3_endpoint_prefix.maximum_returned_per_path(),
      case3_child_index.truncated_endpoint_count()};
}

// Items of one scaffold-only derivation (G' empty) for the forced set of
// WFinalVpWmbpOnlineWidths.  Every structure contains all scaffold pairs, so
// a finite beam that prunes the items realising them in the pseudoknot-free
// way can lose a pair entirely and disconnect the root.  This runs the same
// beam in LogMax with every deduction that generates a G' pair zeroed, so
// only scaffold-only derivations compete, using the same finite beam with
// width at least 1000. The traceback of its root is one such derivation.
// Returns an empty set when none exists.
template <class Context>
std::shared_ptr<const beam::ItemKeySet> scaffold_only_derivation_items(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths) {
  using Semiring = semirings::LogMaxSemiring;
  using Value = typename Semiring::value_type;
  auto wide = WFinalVpWmbpOnlineWidths::uniform(std::max<std::size_t>(
      1000, std::min(widths.v, std::numeric_limits<std::size_t>::max() / 2)));
  wide.max_hairpin_unpaired = widths.max_hairpin_unpaired;
  wide.max_case3_ancestor_depth = widths.max_case3_ancestor_depth;
  wide.force_scaffold_derivation = false;
  wide.memory_observer = widths.memory_observer;
  wide.memory_column_interval = widths.memory_column_interval;
  wide.memory_phase = beam::MemoryPhase::Scaffold;
  wide.numeric_lookup_index = widths.numeric_lookup_index;
  engine::ChoiceMap<Value> choices;
  beam::RetainedChoices<Semiring> retention(&choices);
  auto dp = run_w_final_vp_wmbp_online_impl<Semiring>(
      ctx, n, wide,
      [&](const deductions::ItemKey &item,
          const deductions::Deduction &deduction,
          Value contribution,
          const engine::Chart<Value> &chart,
          auto) {
        if (Semiring::is_zero(contribution)) return;
        // With standard contexts VP is zero in this private G'=empty pass:
        // fixed endpoints are inadmissible and novel generated pairs are zeroed.
        // Every WMBP mode needs VP, so no nonempty WMB connects WI/WIP/BE to W.
        // The private numeric path may omit these root-unconsumed families.
        // Fallbacks keep their numeric states but still omit unreachable choices.
        if constexpr (std::is_same_v<Context, WFinalLogWeightContext> ||
                      std::is_same_v<Context, WFinalExactInsideContext>) {
          if (!widths.memory_observer) {
            using N = grammar::NonTerminal;
            switch (item.nonterminal) {
              case N::WI:
              case N::WIP:
              case N::BE:
                return;
              case N::VM:
                // V_VM is the only active VM consumer, at this same span.
                if (!engine::is_fixed_scaffold_pair(engine::TracebackPair{
                    item.i, item.j, engine::TracebackPairLane::Square}, ctx)) return;
                break;
              default:
                break;
            }
          }
        }
        retention.before_parent(item, chart);
        choices.observe_with_tie_break(
            item, deduction, contribution,
            [](Value lhs, Value rhs) {
              return Semiring::traceback_priority(lhs) >
                     Semiring::traceback_priority(rhs) + 1e-12;
            },
            [](const auto &lhs, const auto &rhs) {
              const auto key = [](const deductions::Deduction &d) {
                return std::make_tuple(static_cast<int>(d.rule), d.split.k,
                                       d.split.l, d.split.p, d.split.q);
              };
              return key(lhs) < key(rhs);
            });
      },
      [&](const deductions::Deduction &deduction, Value contribution) {
        std::vector<engine::TracebackPair> generated;
        engine::append_deduction_generated_pairs(deduction, ctx, &generated);
        return generated.empty() ? contribution : Semiring::zero();
      }, {}, true);
  retention.finish(dp.chart);
  if (widths.memory_observer) {
    auto observation = beam::chart_memory(beam::MemoryPhase::Scaffold,
        "choices_end", wide.v, n, dp.chart);
    observation.containers.push_back(beam::hash_memory("choices", choices.values()));
    widths.memory_observer(observation);
  }
  auto items = std::make_shared<beam::ItemKeySet>();
  const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
  if (Semiring::is_zero(dp.chart.get(root))) return items;
  const auto trace = engine::build_traceback_tree(root, dp.chart, choices);
  std::vector<const engine::DeductionTraceNode<Value> *> stack{&trace};
  while (!stack.empty()) {
    const auto *node = stack.back();
    stack.pop_back();
    items->insert(node->item);
    for (const auto &child : node->children) stack.push_back(&child);
  }
  return items;
}

template <class Semiring,
          class Context,
          class CandidateObserverFn,
          class DeductionValueTransformFn>
WFinalVpWmbpOnlineRun<typename Semiring::value_type>
run_w_final_vp_wmbp_online_observed_transformed(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths,
    CandidateObserverFn &&candidate_observer,
    DeductionValueTransformFn &&deduction_value_transform,
    const WFinalPublicItemRanking &public_item_ranking = {}) {
  constexpr auto keep_all = std::numeric_limits<std::size_t>::max();
  const bool finite = widths.v != keep_all || widths.vm != keep_all ||
                      widths.wm != keep_all || widths.wmv != keep_all ||
                      widths.wmp != keep_all;
  if (widths.force_scaffold_derivation && !widths.forced_items && finite) {
    auto forced = widths;
    forced.forced_items = scaffold_only_derivation_items(ctx, n, widths);
    return run_w_final_vp_wmbp_online_impl<Semiring>(
        ctx, n, forced,
        std::forward<CandidateObserverFn>(candidate_observer),
        std::forward<DeductionValueTransformFn>(deduction_value_transform),
        public_item_ranking);
  }
  return run_w_final_vp_wmbp_online_impl<Semiring>(
      ctx, n, widths,
      std::forward<CandidateObserverFn>(candidate_observer),
      std::forward<DeductionValueTransformFn>(deduction_value_transform),
      public_item_ranking);
}

// Compatibility entry for research observations. Other consumers use the
// same hook through the generic observed runner; regeneration must receive
// that hook too, especially for WMB views and BE prefix FIFO ordering.
template <class Context>
WFinalVpWmbpOnlineRun<semirings::LogSumExpSemiring::value_type>
run_w_final_vp_wmbp_online_ranked_inside(
    Context &ctx, cand_pos_t n, const WFinalVpWmbpOnlineWidths &widths,
    const WFinalPublicItemRanking &ranking) {
  return run_w_final_vp_wmbp_online_observed_transformed<
      semirings::LogSumExpSemiring>(
      ctx, n, widths,
      [](const auto &, const auto &, auto, const auto &, auto) {},
      [](const deductions::Deduction &, auto contribution) { return contribution; },
      ranking);
}

template <class Semiring, class Context, class CandidateObserverFn>
WFinalVpWmbpOnlineRun<typename Semiring::value_type>
run_w_final_vp_wmbp_online_observed(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths,
    CandidateObserverFn &&candidate_observer,
    const WFinalPublicItemRanking &public_item_ranking = {}) {
  return run_w_final_vp_wmbp_online_observed_transformed<Semiring>(
      ctx,
      n,
      widths,
      std::forward<CandidateObserverFn>(candidate_observer),
      [](const deductions::Deduction &,
         typename Semiring::value_type contribution) {
        return contribution;
      }, public_item_ranking);
}

template <class Semiring, class Context>
WFinalVpWmbpOnlineRun<typename Semiring::value_type>
run_w_final_vp_wmbp_online(Context &ctx,
                           cand_pos_t n,
                           const WFinalVpWmbpOnlineWidths &widths,
                           const WFinalPublicItemRanking &public_item_ranking = {}) {
  return run_w_final_vp_wmbp_online_observed<Semiring>(
      ctx,
      n,
      widths,
      [](const deductions::ItemKey &,
         const deductions::Deduction &,
         typename Semiring::value_type,
         const auto &,
         auto) {}, public_item_ranking);
}

template <class Context>
WFinalVpWmbpOnlineRun<
    semirings::TargetedCoupledLogSumMaxSemiring::value_type>
run_w_final_vp_wmbp_online_targeted(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths,
    const std::vector<engine::TracebackPair> &target_pairs,
    const WFinalPublicItemRanking &public_item_ranking = {}) {
  using Semiring = semirings::TargetedCoupledLogSumMaxSemiring;
  const engine::GeneratedTargetIndex target_index(
      engine::generated_traceback_pairs(target_pairs, ctx));
  const bool target_compatible = target_index.compatible_with_scaffold(ctx, n);
  return run_w_final_vp_wmbp_online_observed_transformed<Semiring>(
      ctx,
      n,
      widths,
      [](const deductions::ItemKey &,
         const deductions::Deduction &,
         Semiring::value_type,
         const auto &,
         auto) {},
      [&](const deductions::Deduction &deduction,
          Semiring::value_type contribution) {
        if (target_compatible &&
            target_index.matches_local_emissions(deduction, ctx)) {
          return contribution;
        }
        return Semiring::reject_target(contribution);
      }, public_item_ranking);
}

struct WFinalVpWmbpOnlineTargetProbabilityResult {
  bool ok = false;
  bool denominator_reachable = false;
  bool target_included = false;
  double log_partition = -std::numeric_limits<double>::infinity();
  double log_target_weight = -std::numeric_limits<double>::infinity();
  double log_probability = -std::numeric_limits<double>::infinity();
  double probability = 0.0;
};

template <class Context>
WFinalVpWmbpOnlineTargetProbabilityResult
run_w_final_vp_wmbp_online_target_probability(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths,
    const std::vector<engine::TracebackPair> &target_pairs,
    const WFinalPublicItemRanking &public_item_ranking = {}) {
  const auto run =
      run_w_final_vp_wmbp_online_targeted(ctx, n, widths, target_pairs, public_item_ranking);
  const deductions::ItemKey root{
      grammar::NonTerminal::W, 1, n, -1, -1};
  const auto value = run.chart.get(root);
  WFinalVpWmbpOnlineTargetProbabilityResult out;
  out.log_partition = value.log_sum;
  out.log_target_weight = value.log_target;
  out.denominator_reachable = std::isfinite(value.log_sum);
  out.target_included = std::isfinite(value.log_target);
  if (!out.denominator_reachable) return out;
  if (!out.target_included) {
    out.ok = true;
    return out;
  }
  out.log_probability = value.log_target - value.log_sum;
  if (!std::isfinite(out.log_probability) || out.log_probability > 1e-12) {
    return out;
  }
  out.probability = std::exp(std::min(0.0, out.log_probability));
  out.ok = std::isfinite(out.probability) && out.probability <= 1.0;
  return out;
}

template <class Semiring, class Context>
WFinalVpWmbpOnlineRun<typename Semiring::value_type>
run_w_final_vp_wmbp_online(Context &ctx,
                           cand_pos_t n,
                           std::size_t width) {
  return run_w_final_vp_wmbp_online<Semiring>(
      ctx, n, WFinalVpWmbpOnlineWidths::vp_wmbp(width));
}

template <class Value, class Choices = engine::ChoiceMap<Value>>
struct BasicWFinalVpWmbpOnlineViterbiRun {
  WFinalVpWmbpOnlineRun<Value> dp;
  Choices choices;
};

using WFinalVpWmbpOnlineViterbiRun =
    BasicWFinalVpWmbpOnlineViterbiRun<double>;
using WFinalVpWmbpOnlineCoupledViterbiRun =
    BasicWFinalVpWmbpOnlineViterbiRun<
        semirings::CoupledLogSumMaxSemiring::value_type>;

inline auto online_traceback_deduction_tie_key(
    const deductions::Deduction &deduction) {
  const auto item_key = [](const deductions::ItemKey &item) {
    return std::make_tuple(static_cast<unsigned char>(item.nonterminal),
                           item.i, item.j, item.ip, item.jp);
  };
  return std::make_tuple(static_cast<unsigned short>(deduction.rule),
                         deduction.split.k,
                         deduction.split.l,
                         deduction.split.p,
                         deduction.split.q,
                         deduction.child_count,
                         item_key(deduction.children[0].item),
                         item_key(deduction.children[1].item),
                         item_key(deduction.children[2].item));
}

template <class Semiring, class Context>
BasicWFinalVpWmbpOnlineViterbiRun<typename Semiring::value_type>
run_w_final_vp_wmbp_online_viterbi_with_semiring(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths,
    const WFinalPublicItemRanking &public_item_ranking = {});

template <class Context>
WFinalVpWmbpOnlineViterbiRun run_w_final_vp_wmbp_online_viterbi(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths) {
  return run_w_final_vp_wmbp_online_viterbi_with_semiring<
      semirings::MaxProductSemiring>(ctx, n, widths);
}

template <class Context>
WFinalVpWmbpOnlineCoupledViterbiRun
run_w_final_vp_wmbp_online_coupled_viterbi(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths) {
  return run_w_final_vp_wmbp_online_viterbi_with_semiring<
      semirings::CoupledLogSumMaxSemiring>(ctx, n, widths);
}

template <class Choices>
auto append_choice_fallback_memory(const Choices &choices,
    beam::MemoryObservation &observation, int)
    -> decltype(choices.fallback_values(), void()) {
  observation.containers.push_back(
      beam::hash_memory("choice_fallbacks", choices.fallback_values()));
}
template <class Choices>
void append_choice_fallback_memory(const Choices &, beam::MemoryObservation &, long) {}

template <class Semiring, class Choices, class Context>
BasicWFinalVpWmbpOnlineViterbiRun<typename Semiring::value_type, Choices>
run_w_final_vp_wmbp_online_viterbi_stored(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths,
    const WFinalPublicItemRanking &public_item_ranking) {
  using Value = typename Semiring::value_type;
  Choices choices;
  beam::RetainedChoices<Semiring, Choices> retention(&choices);
  auto dp = run_w_final_vp_wmbp_online_observed<Semiring>(
      ctx,
      n,
      widths,
      [&](const deductions::ItemKey &item,
          const deductions::Deduction &deduction,
          Value contribution,
          const engine::Chart<Value> &chart,
          auto) {
        if (Semiring::is_zero(contribution)) return;
        retention.before_parent(item, chart);
        choices.observe_with_tie_break(
            item,
            deduction,
            contribution,
            [](Value lhs, Value rhs) {
              constexpr double tie_tolerance = 1e-12;
              return Semiring::traceback_priority(lhs) >
                     Semiring::traceback_priority(rhs) + tie_tolerance;
            },
            [](const auto &lhs, const auto &rhs) {
              return online_traceback_deduction_tie_key(lhs) <
                     online_traceback_deduction_tie_key(rhs);
            });
      }, public_item_ranking);
  retention.finish(dp.chart);
  if (widths.memory_observer) {
    auto observation = beam::chart_memory(widths.memory_phase,
        "choices_end", widths.v, n, dp.chart);
    observation.containers.push_back(beam::hash_memory("choices", choices.values()));
    append_choice_fallback_memory(choices, observation, 0);
    widths.memory_observer(observation);
  }
  return BasicWFinalVpWmbpOnlineViterbiRun<Value, Choices>{
      std::move(dp), std::move(choices)};
}

template <class Semiring, class Context>
BasicWFinalVpWmbpOnlineViterbiRun<typename Semiring::value_type>
run_w_final_vp_wmbp_online_viterbi_with_semiring(
    Context &ctx, cand_pos_t n, const WFinalVpWmbpOnlineWidths &widths,
    const WFinalPublicItemRanking &public_item_ranking) {
  return run_w_final_vp_wmbp_online_viterbi_stored<Semiring,
      engine::ChoiceMap<typename Semiring::value_type>>(
          ctx, n, widths, public_item_ranking);
}

template <class Context>
engine::ViterbiTracebackResult<double>
run_w_final_vp_wmbp_online_viterbi_traceback(Context &ctx,
                                              cand_pos_t n,
                                              const WFinalVpWmbpOnlineWidths &
                                                  widths) {
  auto online = run_w_final_vp_wmbp_online_viterbi(ctx, n, widths);
  const deductions::ItemKey root{
      grammar::NonTerminal::W, 1, n, -1, -1};
  return engine::build_viterbi_traceback_result(
      root,
      engine::ViterbiRunResult<double>{
          std::move(online.dp.chart), std::move(online.choices)},
      ctx);
}

template <class Context>
engine::ViterbiTracebackResult<
    semirings::CoupledLogSumMaxSemiring::value_type>
run_w_final_vp_wmbp_online_coupled_viterbi_traceback(
    Context &ctx,
    cand_pos_t n,
    const WFinalVpWmbpOnlineWidths &widths) {
  auto online = run_w_final_vp_wmbp_online_coupled_viterbi(ctx, n, widths);
  const deductions::ItemKey root{
      grammar::NonTerminal::W, 1, n, -1, -1};
  return engine::build_viterbi_traceback_result(
      root,
      engine::ViterbiRunResult<
          semirings::CoupledLogSumMaxSemiring::value_type>{
          std::move(online.dp.chart), std::move(online.choices)},
      ctx);
}

}  // namespace cparty::scfg2::replay

#endif  // CPARTY_SCFG2_REPLAY_W_FINAL_VP_WMBP_ONLINE_HH

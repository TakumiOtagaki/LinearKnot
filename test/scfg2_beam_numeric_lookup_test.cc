#include "scfg2/beam/prepare_ranking.hh"
#include "scfg2/weights/runtime_snapshot.hh"

#include <array>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace sc = cparty::scfg2;
using Key = sc::deductions::ItemKey;
using Deduction = sc::deductions::Deduction;
using NT = sc::grammar::NonTerminal;
using Sum = sc::semirings::LogSumExpSemiring;
using Max = sc::semirings::LogMaxSemiring;
using Chart = sc::engine::Chart<double>;
using Choices = sc::engine::ChoiceMap<double>;
using Run = sc::replay::WFinalVpWmbpOnlineRun<double>;
using Widths = sc::replay::WFinalVpWmbpOnlineWidths;
using Observation = sc::beam::MemoryObservation;
std::size_t cases = 0, compared_items = 0, compared_choices = 0;
std::size_t retained_be = 0, retained_vp = 0, first_winner_ties = 0;
std::size_t guide_cases = 0, guide_links = 0, guide_misses = 0;
std::size_t contribution_queries = 0, compared_candidates = 0;
std::array<std::size_t, 6> aggregate_hits{};
std::string case_description;

void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(case_description + ": " + message);
}
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(a)) == 0; }
bool equal(const Key &a, const Key &b) { return sc::deductions::ItemKeyEq{}(a, b); }
bool equal(const Deduction &a, const Deduction &b) {
  if (!equal(a.parent, b.parent) || a.rule != b.rule || a.child_count != b.child_count ||
      a.split.k != b.split.k || a.split.l != b.split.l ||
      a.split.p != b.split.p || a.split.q != b.split.q) return false;
  // Include unused child slots: serialized/observed deduction identities must
  // remain exactly the same, rather than only describing the same product.
  for (std::size_t c = 0; c < a.children.size(); ++c)
    if (!equal(a.children[c].item, b.children[c].item)) return false;
  return true;
}

struct Candidate { Key key; Deduction deduction; double contribution, local; };

void compare_candidates(const std::vector<Candidate> &a, const std::vector<Candidate> &b) {
  require(a.size() == b.size(), "candidate stream length differs");
  for (std::size_t c = 0; c < a.size(); ++c) {
    require(equal(a[c].key, b[c].key) && equal(a[c].deduction, b[c].deduction) &&
                same(a[c].contribution, b[c].contribution) && same(a[c].local, b[c].local),
            "candidate stream order/full deduction/contribution/local bits differ");
    ++compared_candidates;
  }
}

void compare_charts(const Chart &reference, const Chart &indexed, cand_pos_t n) {
  require(reference.primary_lookup_length() == 0 &&
              indexed.primary_lookup_length() == n,
          "numeric lookup flag did not reach the main chart");
  require(reference.dense_length() == 0 && indexed.dense_length() == 0,
          "finite beam enabled dense chart storage");
  const auto &a = reference.values();
  const auto &b = indexed.values();
  require(a.size() == b.size() && reference.size() == indexed.size(),
          "retained chart sizes differ");
  require(a.bucket_count() == b.bucket_count(), "canonical map bucket growth differs");
  auto right = b.begin();
  for (const auto &[key, value] : a) {
    require(right != b.end() && equal(key, right->first) && same(value, right->second),
            "retained map iteration/key/value bits differ");
    require(indexed.contains(key) && same(value, indexed.get(key)),
            "retained numeric lookup differs from canonical map");
    retained_be += key.nonterminal == NT::BE;
    retained_vp += key.nonterminal == NT::VP;
    ++right;
    ++compared_items;
  }
  require(right == b.end(), "indexed chart added retained keys");
  auto cursor = a.begin();
  indexed.for_each_entry([&](const Key &key, double value) {
    require(cursor != a.end() && equal(key, cursor->first) && same(value, cursor->second),
            "for_each_entry changed canonical iteration");
    ++cursor;
  });
  require(cursor == a.end(), "for_each_entry lost retained keys");
  // Check all primary misses too, including cells that were positive before
  // pruning. This detects stale values after retention-driven erasures.
  for (int family = 0; family < static_cast<int>(NT::BE); ++family)
    for (cand_pos_t j = 1; j <= n; ++j)
      for (cand_pos_t i = 1; i <= j; ++i) {
        const Key key{static_cast<NT>(family), i, j, -1, -1};
        require(reference.contains(key) == indexed.contains(key) &&
                    same(reference.get(key), indexed.get(key)),
                "pruned primary lookup differs from canonical map");
      }
}

void compare_work(const Run &a, const Run &b) {
#define SAME_WORK(field) require(a.field == b.field, "work changed: " #field)
  SAME_WORK(item_visits); SAME_WORK(positive_items);
  SAME_WORK(frontier_candidates); SAME_WORK(frontier_selected);
  SAME_WORK(schedule.vp_parent_visits); SAME_WORK(schedule.be_parent_visits);
  SAME_WORK(schedule.wmbp_parent_visits); SAME_WORK(schedule.reverse_parent_vp_probes);
  SAME_WORK(schedule.family_parent_visits);
  SAME_WORK(vp_candidates); SAME_WORK(vp_selected);
  SAME_WORK(wmbp_candidates); SAME_WORK(wmbp_selected);
  SAME_WORK(case3_producer_candidates); SAME_WORK(case3_producer_selected);
  SAME_WORK(case3_tail_candidates); SAME_WORK(case3_tail_selected);
  SAME_WORK(case3_left_candidates); SAME_WORK(case3_left_selected);
  SAME_WORK(be_right_tail_candidates); SAME_WORK(be_right_tail_selected);
  SAME_WORK(be_right_shell_candidates); SAME_WORK(be_right_shell_selected);
  SAME_WORK(wmbp_prefix_candidates); SAME_WORK(wmbp_prefix_selected);
  SAME_WORK(vm_binary_join_probes); SAME_WORK(vm_base_prefix_cell_visits);
  SAME_WORK(vm_base_prefix_selected);
  SAME_WORK(wm_binary_join_probes); SAME_WORK(wm_extension_probes);
  SAME_WORK(wm_start_prefix_cell_visits); SAME_WORK(wm_start_prefix_selected);
  SAME_WORK(wi_binary_join_probes); SAME_WORK(wi_direct_probes); SAME_WORK(wi_extension_probes);
  SAME_WORK(wip_binary_join_probes); SAME_WORK(wip_extension_probes);
  SAME_WORK(wip_start_prefix_cell_visits); SAME_WORK(wip_start_prefix_selected);
  SAME_WORK(wmbw_binary_join_probes); SAME_WORK(vpl_prefix_cell_visits);
  SAME_WORK(vpl_prefix_selected); SAME_WORK(vpr_binary_join_probes);
  SAME_WORK(vpr_base_state_probes); SAME_WORK(vpr_base_state_selected);
  SAME_WORK(w_binary_join_probes); SAME_WORK(w_direct_probes);
  SAME_WORK(w_extension_probes); SAME_WORK(w_root_decomposition_probes);
  SAME_WORK(wmb_split_positive_sources); SAME_WORK(wmb_direct_positive_sources);
  SAME_WORK(wmb_empty_positive_sources); SAME_WORK(wmb_btail_join_probes);
  SAME_WORK(wmb_btail_direct_probes); SAME_WORK(wmb_btail_selected);
  SAME_WORK(wmb_be_btail_join_probes); SAME_WORK(wmb_retained_direct_probes);
  SAME_WORK(wmb_retained_split_selected);
  SAME_WORK(wmbp_shared_wmbp_left_probes); SAME_WORK(wmbp_shared_wmbw_left_probes);
  SAME_WORK(wmbp_shared_emitted_sources); SAME_WORK(wmbp_wi_history_insertions);
  SAME_WORK(wmbp_wi_be_wi_probes); SAME_WORK(wmbp_wi_empty_probes);
  SAME_WORK(wmbp_wi_prefix_candidates); SAME_WORK(wmbp_wi_prefix_vp_probes);
  SAME_WORK(wmbp_wi_direct_probes); SAME_WORK(wmbp_wi_emitted_sources);
  SAME_WORK(vp_loop_parent_probes); SAME_WORK(vp_stack_event_sources);
  SAME_WORK(vp_internal_event_sources); SAME_WORK(vp_recursive_split_join_probes);
  SAME_WORK(vp_recursive_split_emitted_sources); SAME_WORK(vp_recursive_split_duplicate_sources);
  SAME_WORK(vp_case3_nonleft_wi_queries); SAME_WORK(vp_case3_nonleft_geometry_probes);
  SAME_WORK(vp_case3_nonleft_right_last_probes); SAME_WORK(vp_case3_nonleft_middle_last_probes);
  SAME_WORK(vp_case3_nonleft_zero_child_probes); SAME_WORK(vp_case3_nonleft_emitted_sources);
  SAME_WORK(vp_case3_nonleft_border_mismatches); SAME_WORK(vp_case3_nonleft_duplicate_sources);
  SAME_WORK(vp_case3_nonleft_prefix_candidates); SAME_WORK(vp_case3_nonleft_active_prefix_probes);
  SAME_WORK(vp_case3_nonleft_max_active_prefixes); SAME_WORK(vp_case12_emitted_sources);
  SAME_WORK(vp_case12_border_mismatches); SAME_WORK(vp_case12_duplicate_sources);
  SAME_WORK(be_stack_sources); SAME_WORK(be_base_stack_duplicate_sources);
  SAME_WORK(be_core_left_sources); SAME_WORK(be_core_left_duplicate_sources);
  SAME_WORK(be_right_tail_duplicate_registrations); SAME_WORK(be_right_shell_sources);
  SAME_WORK(be_right_shell_duplicate_sources); SAME_WORK(case3_active_expiry_probes);
  SAME_WORK(case3_tail_materialization_probes); SAME_WORK(case3_prepared_score_probes);
  SAME_WORK(case3_max_active_segments); SAME_WORK(v_hairpin_parent_probes);
  SAME_WORK(v_internal_parent_probes); SAME_WORK(v_vm_parent_probes);
  SAME_WORK(v_hairpin_sources); SAME_WORK(v_internal_sources); SAME_WORK(v_vm_sources);
  SAME_WORK(case3_heap_registered_segments); SAME_WORK(case3_heap_expired_segments);
  SAME_WORK(case3_max_heap_size); SAME_WORK(case3_endpoint_prefix_registrations);
  SAME_WORK(case3_endpoint_prefix_late_registrations);
  SAME_WORK(case3_endpoint_prefix_nonmonotone_queries);
  SAME_WORK(case3_endpoint_prefix_maximum_size);
  SAME_WORK(case3_endpoint_prefix_maximum_returned_per_path);
  SAME_WORK(case3_truncated_endpoints);
#undef SAME_WORK
}

void compare_choices(const Choices &reference, const Choices &indexed) {
  const auto &a = reference.values();
  const auto &b = indexed.values();
  require(!a.empty() && a.size() == b.size() && a.bucket_count() == b.bucket_count(),
          "stored choice count/bucket growth differs");
  auto right = b.begin();
  for (const auto &[key, value] : a) {
    require(right != b.end() && equal(key, right->first) &&
                value.has_value == right->second.has_value &&
                same(value.contribution, right->second.contribution) &&
                equal(value.deduction, right->second.deduction),
            "stored choice iteration/key/contribution/deduction differs");
    ++right;
    ++compared_choices;
  }
  require(right == b.end(), "lookup added stored choices");
}

void compare_memory(const std::vector<Observation> &reference,
                    const std::vector<Observation> &indexed, bool outer) {
  require(reference.size() == indexed.size(), "memory observation event count differs");
  bool scaffold = false, pilot = false, main = false, outer_guide = false;
  for (std::size_t i = 0; i < reference.size(); ++i) {
    const auto &a = reference[i];
    const auto &b = indexed[i];
    require(a.phase == b.phase && std::string(a.event) == b.event && a.width == b.width &&
                a.column == b.column && a.chart_items == b.chart_items,
            "memory observation chart/event metadata differs");
    require(a.families.size() == b.families.size(), "observed family count differs");
    for (std::size_t f = 0; f < a.families.size(); ++f)
      require(a.families[f].nonterminal == b.families[f].nonterminal &&
                  a.families[f].live_items == b.families[f].live_items,
              "observed retained family sizes differ");
    std::size_t cursor = 0, lookup_bytes = 0;
    for (const auto &container : b.containers) {
      if (container.name == "chart_primary_lookup") {
        require(container.vector_capacity_bytes > 0 && container.accounted_bytes() > 0,
                "memory observer did not account lookup capacity");
        lookup_bytes += container.vector_capacity_bytes;
        continue;
      }
      require(cursor < a.containers.size(), "lookup changed canonical memory containers");
      const auto &expected = a.containers[cursor++];
      require(expected.name == container.name && expected.entries == container.entries &&
                  expected.buckets == container.buckets && expected.value_bytes == container.value_bytes &&
                  expected.vector_size_bytes == container.vector_size_bytes &&
                  expected.vector_capacity_bytes == container.vector_capacity_bytes,
              "lookup changed canonical/container memory accounting");
    }
    require(cursor == a.containers.size(), "lookup omitted canonical memory containers");
    if (lookup_bytes) {
      scaffold |= b.phase == sc::beam::MemoryPhase::Scaffold;
      pilot |= b.phase == sc::beam::MemoryPhase::Pilot;
      main |= b.phase == sc::beam::MemoryPhase::Main;
      outer_guide |= b.phase == sc::beam::MemoryPhase::Outer &&
          std::string(b.event) == "guide_end";
    }
  }
  require(scaffold && main && pilot == outer && outer_guide == outer,
          "lookup opt-in/capacity did not propagate through forcing/pilot/guide/main");
}

struct Input { std::string name, sequence, scaffold; };

Input nested_input() {
  Input out{"nested", "GCGCGCGCGCGCGCGCGCGCGCGC", "....((....))............"};
  return out;
}

void check(const Input &input, std::size_t width, sc::WmbClosurePolicy closure,
           sc::beam::RankingMode ranking) {
  const bool outer = ranking == sc::beam::RankingMode::Outer;
  case_description = input.name + " width=" + std::to_string(width) +
      (outer ? " outer" : " baseline") +
      (closure == sc::WmbClosurePolicy::LegacyCompatible ? " legacy" : " weak");
  require(input.sequence.size() == input.scaffold.size(), "fixture length differs");
  const cand_pos_t n = input.sequence.size();
  sc::weights::RuntimeSnapshot owner(input.sequence, input.scaffold);
  sc::replay::WFinalLogWeightContext reference_context(owner, closure), indexed_context(owner, closure);
  std::vector<Observation> observations[2];
  auto plain_widths = Widths::uniform(width), indexed_widths = plain_widths;
  indexed_widths.numeric_lookup_index = true;
  plain_widths.memory_column_interval = indexed_widths.memory_column_interval = n;
  plain_widths.memory_observer = [&](const auto &o) { observations[0].push_back(o); };
  indexed_widths.memory_observer = [&](const auto &o) { observations[1].push_back(o); };
  auto reference = sc::beam::prepare_ranking(reference_context, n, plain_widths, width, ranking);
  auto indexed = sc::beam::prepare_ranking(indexed_context, n, indexed_widths, width, ranking);
  require(!reference.widths.numeric_lookup_index && indexed.widths.numeric_lookup_index,
          "prepare_ranking dropped the numeric lookup setting");
  require(reference.pilot_width == indexed.pilot_width && reference.guide_items == indexed.guide_items,
          "lookup changed pilot width or guide size");
  if (outer) {
    require(reference.widths.forced_items && indexed.widths.forced_items &&
                !reference.widths.forced_items->empty() &&
                reference.widths.forced_items->size() == indexed.widths.forced_items->size(),
            "forcing prepass produced different/vacuous items");
    for (const auto &key : *reference.widths.forced_items)
      require(indexed.widths.forced_items->count(key), "forcing prepass changed a retained key");
  }
  std::vector<Candidate> candidates[2];
  const auto observer = [&](std::size_t side) {
    return [&, side](const Key &key, const Deduction &d, double value, const Chart &, auto local) {
      candidates[side].push_back({key, d, value, sc::weights::log_weight(local)});
    };
  };
  const auto sum_reference = sc::replay::run_w_final_vp_wmbp_online_observed<Sum>(
      reference_context, n, reference.widths, observer(0), reference.ranking);
  const auto sum_indexed = sc::replay::run_w_final_vp_wmbp_online_observed<Sum>(
      indexed_context, n, indexed.widths, observer(1), indexed.ranking);
  compare_candidates(candidates[0], candidates[1]);
  compare_charts(sum_reference.chart, sum_indexed.chart, n);
  compare_work(sum_reference, sum_indexed);
  const Key root{NT::W, 1, n, -1, -1};
  require(std::isfinite(sum_reference.chart.get(root)), "finite-beam root is unreachable");
  const auto max_reference = sc::replay::run_w_final_vp_wmbp_online_viterbi_with_semiring<Max>(
      reference_context, n, reference.widths, reference.ranking);
  const auto max_indexed = sc::replay::run_w_final_vp_wmbp_online_viterbi_with_semiring<Max>(
      indexed_context, n, indexed.widths, indexed.ranking);
  compare_charts(max_reference.dp.chart, max_indexed.dp.chart, n);
  compare_work(max_reference.dp, max_indexed.dp);
  compare_choices(max_reference.choices, max_indexed.choices);
  compare_memory(observations[0], observations[1], outer);
  // The numeric frontier storage is private to unobserved runs. Compare its
  // complete stream, chart and first-winner choices with the observed path,
  // whose original containers remain part of the memory API contract.
  auto numeric_widths = indexed.widths;
  numeric_widths.memory_observer = {};
  candidates[1].clear();
  const auto sum_numeric = sc::replay::run_w_final_vp_wmbp_online_observed<Sum>(
      indexed_context, n, numeric_widths, observer(1), indexed.ranking);
  compare_candidates(candidates[0], candidates[1]);
  compare_charts(sum_reference.chart, sum_numeric.chart, n);
  compare_work(sum_indexed, sum_numeric);
  const auto max_numeric = sc::replay::run_w_final_vp_wmbp_online_viterbi_with_semiring<Max>(
      indexed_context, n, numeric_widths, indexed.ranking);
  compare_charts(max_reference.dp.chart, max_numeric.dp.chart, n);
  compare_work(max_indexed.dp, max_numeric.dp);
  compare_choices(max_indexed.choices, max_numeric.choices);
  ++cases;
}

void compare_guide_work(const sc::beam::outer_detail::Result &a,
                        const sc::beam::outer_detail::Result &b) {
#define SAME_GUIDE(field) require(a.field == b.field, "guide work changed: " #field)
  SAME_GUIDE(outer.visited); SAME_GUIDE(outer.generated);
  SAME_GUIDE(outer.eligible); SAME_GUIDE(outer.selected);
  SAME_GUIDE(outer.max_parent_candidates); SAME_GUIDE(outer.peak_pending);
  SAME_GUIDE(outer.peak_guide_states); SAME_GUIDE(outer.evicted);
  SAME_GUIDE(outer.generated_rules); SAME_GUIDE(outer.eligible_rules);
  SAME_GUIDE(work.input_records); SAME_GUIDE(work.guide_records);
  SAME_GUIDE(work.index_records); SAME_GUIDE(work.comparisons);
  SAME_GUIDE(work.build_moves); SAME_GUIDE(work.undo_moves);
  SAME_GUIDE(work.cursor_steps); SAME_GUIDE(work.query_outputs);
  SAME_GUIDE(work.proposals); SAME_GUIDE(work.max_proposals);
  SAME_GUIDE(work.checks); SAME_GUIDE(work.weight_evaluations);
  SAME_GUIDE(work.child_alpha_lookups); SAME_GUIDE(work.leaf_deductions);
  SAME_GUIDE(work.zero_sibling_deductions);
#undef SAME_GUIDE
}

template <bool TrackLinks>
void direct_guide_lookup(sc::replay::WFinalLogWeightContext &ctx,
                         const Chart &pilot, cand_pos_t n, std::size_t width) {
  const Key root{NT::W, 1, n, -1, -1};
  // Enable the selected/full-generator consistency validator in both runs.
  // Profiling is disabled; elapsed/validation timings are deliberately not
  // part of the differential contract.
  const auto reference = sc::beam::bounded_outer<sc::replay::WFinalLogWeightContext, true, TrackLinks>(
      pilot, root, ctx, width, 4 * width, true, false, false);
  const auto indexed = sc::beam::bounded_outer<sc::replay::WFinalLogWeightContext, true, TrackLinks>(
      pilot, root, ctx, width, 4 * width, true, false, true);
  compare_charts(reference.outer.values, indexed.outer.values, n);
  compare_guide_work(reference, indexed);
  require(reference.work.checks > 0 && reference.outer.values.size() > 1,
          "direct guide lookup witness is vacuous");
  const auto &a = reference.outer.links;
  const auto &b = indexed.outer.links;
  require(a.size() == b.size() && a.bucket_count() == b.bucket_count(),
          "guide link size/bucket growth changed");
  auto cursor = b.begin();
  for (const auto &[key, link] : a) {
    require(cursor != b.end() && equal(key, cursor->first) &&
                equal(link.parent, cursor->second.parent) && same(link.edge, cursor->second.edge),
            "guide link iteration/full key/parent/edge bits changed");
    ++cursor;
    ++guide_links;
  }
  require(cursor == b.end(), "indexed guide added predecessor links");
  if constexpr (TrackLinks)
    require(!a.empty() && a.size() + 1 == reference.outer.values.size(),
            "stored guide predecessor witness is vacuous");
  else require(a.empty(), "TrackLinks=false retained a predecessor map");

  // These complete keys cannot be generated by the guide. Primary keys
  // carrying secondary fields must fall back, and BE misses must retain all
  // four coordinates even when the primary span matches a live guide item.
  std::vector<Key> misses{
      {NT::W, 1, n, 0, -1}, {NT::V, 1, n, -1, 0},
      {NT::BE, 1, n, -1, -1}, {NT::BE, 1, n, n + 1, n + 2},
      {NT::W, 0, n, -1, -1}, {NT::WI, 2, 1, -1, -1},
      {NT::V, 1, n + 1, -1, -1}, {NT::WM, 1, 0, -1, -1},
      {static_cast<NT>(255), 1, n, -1, -1}};
  for (const auto &[key, value] : reference.outer.values.values()) {
    (void)value;
    if (key.nonterminal == NT::BE) {
      misses.push_back({NT::BE, key.i, key.j, key.ip, n + 1});
      misses.push_back({NT::BE, key.i, key.j, n + 1, key.jp});
      break;
    }
  }
  for (const auto &key : misses) {
    require(!reference.outer.values.contains(key) && !indexed.outer.values.contains(key) &&
                same(reference.outer.values.get(key), Sum::zero()) &&
                same(reference.outer.values.get(key), indexed.outer.values.get(key)),
            "full-key guide miss aliased a completed score");
    ++guide_misses;
  }
  ++guide_cases;
}

void direct_guide_contract(const Input &input, std::size_t width, sc::WmbClosurePolicy closure) {
  case_description = input.name + " direct guide width=" + std::to_string(width) +
      (closure == sc::WmbClosurePolicy::LegacyCompatible ? " legacy" : " weak");
  sc::weights::RuntimeSnapshot owner(input.sequence, input.scaffold);
  sc::replay::WFinalLogWeightContext ctx(owner, closure);
  const cand_pos_t n = input.sequence.size();
  // Use one immutable, ordinary pilot chart for all four guide variants so
  // this checks completed-guide indexing independently of pilot indexing.
  const auto pilot = sc::replay::run_w_final_vp_wmbp_online<Sum>(ctx, n, Widths::uniform(width));
  direct_guide_lookup<true>(ctx, pilot.chart, n, width);
  direct_guide_lookup<false>(ctx, pilot.chart, n, width);
}

template <class Reconstructs, class Contribution, class Find>
void compare_aggregate_lookup(NT family, cand_pos_t j, std::size_t slot,
                              Reconstructs reconstructs, Contribution contribution, Find find) {
  const auto check = [&](const Deduction &d) {
    const bool expected = reconstructs(d);
    const auto *actual = find(d);
    require((actual != nullptr) == expected,
            "one-probe aggregate lookup changed rule/full-parent ownership");
    if (actual) {
      require(same(*actual, contribution(d.parent)),
              "one-probe aggregate lookup changed value bits");
      ++aggregate_hits[slot];
    }
    ++contribution_queries;
  };
  // The old reconstructs_* predicate and value getter remain independent
  // references. Test every rule, including rules owned by the other five
  // agendas, and full-key misses even when their primary span is present.
  for (cand_pos_t i = 0; i <= j + 1; ++i) {
    for (auto rule : sc::grammar::all_rule_ids()) {
      Deduction d{};
      d.parent = {family, i, j, -1, -1};
      d.rule = rule;
      check(d);
      d.parent.ip = 0;
      check(d);
      d.parent.ip = -1;
      d.parent.jp = 0;
      check(d);
      d.parent.jp = -1;
      d.parent.nonterminal = NT::BE;
      check(d);
      d.parent.nonterminal = family;
      ++d.parent.j;
      check(d);
    }
  }
}

template <class Semiring>
void aggregate_lookup_contract(const Input &input, sc::WmbClosurePolicy closure) {
  case_description = input.name + " private aggregate lookup";
  sc::weights::RuntimeSnapshot owner(input.sequence, input.scaffold);
  sc::replay::WFinalLogWeightContext ctx(owner, closure);
  const cand_pos_t n = input.sequence.size();
  // Synthetic positive retained cells exercise each aggregate independently
  // of whether another beam happens to prune its supporting scaffold item.
  Chart chart(Semiring::zero());
  std::vector<std::vector<cand_pos_t>> starts_by_end(n + 1);
  std::vector<std::vector<Key>> be_by_end(n + 1);
  const std::vector<cand_pos_t> empty;
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = 1; i <= j; ++i) {
      starts_by_end[j].push_back(i);
      for (int f = 0; f < static_cast<int>(NT::BE); ++f)
        chart.set({static_cast<NT>(f), i, j, -1, -1}, 0.125 * (f + i % 3));
      if (ctx.pair_partner(i) != j) continue;
      for (cand_pos_t ip = i; ip < j; ++ip) {
        const cand_pos_t jp = ctx.pair_partner(ip);
        if (jp <= ip || jp > j) continue;
        const Key be{NT::BE, i, j, ip, jp};
        be_by_end[j].push_back(be);
        chart.set(be, 0.25 * (ip - i));
      }
    }
  }
  const auto starts = [&](NT, cand_pos_t end) -> const std::vector<cand_pos_t> & {
    return end > 0 && end <= n ? starts_by_end[end] : empty;
  };
  const auto vp_starts = [&](cand_pos_t end) -> const std::vector<cand_pos_t> & {
    return starts(NT::VP, end);
  };
  sc::beam::VmRetainedJoinAgenda<Semiring> vm;
  sc::beam::WmRetainedJoinAgenda<Semiring> wm;
  sc::beam::WipRetainedJoinAgenda<Semiring> wip;
  sc::beam::VplRetainedAgenda<Semiring> vpl;
  sc::beam::VprRetainedAgenda<Semiring> vpr;
  sc::beam::WmbRetainedAgenda<Semiring> wmb;
  sc::beam::WmbBtailStore<Semiring> tails;
  for (cand_pos_t j = 1; j <= n; ++j) {
    vm.build_column(j, chart, ctx, starts, 20);
    wm.build_column(j, chart, ctx, starts, 20);
    wip.build_column(j, chart, ctx, starts, 20);
    vpl.build_column(j, chart, ctx, starts_by_end[j], 20);
    vpr.build_column(j, chart, ctx, starts_by_end[j], vp_starts, 20);
    tails.build(j, chart, ctx, starts(NT::WI, j - 1), vp_starts, 20);
    wmb.build_column(j, chart, ctx, starts_by_end[j], be_by_end[j], tails);
#define COMPARE_AGGREGATE(agenda, family, slot, reconstruct, contribution, find) \
    compare_aggregate_lookup(NT::family, j, slot, \
        [&](const auto &d) { return agenda.reconstruct(d); }, \
        [&](const auto &parent) { return agenda.contribution(parent); }, \
        [&](const auto &d) { return agenda.find(d); })
    COMPARE_AGGREGATE(vm, VM, 0, reconstructs_base, base_contribution, find_base_contribution);
    COMPARE_AGGREGATE(wm, WM, 1, reconstructs_start, start_contribution, find_start_contribution);
    COMPARE_AGGREGATE(wip, WIP, 2, reconstructs_start, start_contribution, find_start_contribution);
    COMPARE_AGGREGATE(vpl, VPL, 3, reconstructs, contribution, find_contribution);
    COMPARE_AGGREGATE(vpr, VPR, 4, reconstructs_base, base_contribution, find_base_contribution);
    COMPARE_AGGREGATE(wmb, WMB, 5, reconstructs_split, split_contribution, find_split_contribution);
#undef COMPARE_AGGREGATE
  }
}

struct FirstWinners {
  Choices choices;
  sc::beam::RetainedChoices<Max> retention{&choices};
  std::vector<Candidate> candidates;
  std::size_t ties = 0;
  void observe(const Key &key, const Deduction &deduction, double contribution,
               const Chart &chart, double local) {
    candidates.push_back({key, deduction, contribution, local});
    if (Max::is_zero(contribution)) return;
    retention.before_parent(key, chart);
    const auto *previous = choices.find(key);
    if (previous && previous->has_value && same(previous->contribution, contribution) &&
        !equal(previous->deduction, deduction)) ++ties;
    choices.observe(key, deduction, contribution, Max::better);
  }
};

void first_winner_contract() {
  const auto input = nested_input();
  case_description = "equal-contribution first winner";
  const cand_pos_t n = input.sequence.size();
  sc::weights::RuntimeSnapshot owner(input.sequence, input.scaffold);
  sc::replay::WFinalLogWeightContext ctx(owner);
  auto widths = Widths::uniform(20);
  FirstWinners reference, indexed;
  const auto run = [&](FirstWinners &observer) {
    return sc::replay::run_w_final_vp_wmbp_online_observed_transformed<Max>(
        ctx, n, widths,
        [&](const Key &key, const Deduction &d, double value, const Chart &chart, auto local) {
          observer.observe(key, d, value, chart, sc::weights::log_weight(local));
        }, [](const Deduction &, double value) { return Max::is_zero(value) ? value : 0.0; });
  };
  const auto a = run(reference);
  widths.numeric_lookup_index = true;
  const auto b = run(indexed);
  reference.retention.finish(a.chart);
  indexed.retention.finish(b.chart);
  compare_charts(a.chart, b.chart, n);
  compare_work(a, b);
  compare_choices(reference.choices, indexed.choices);
  compare_candidates(reference.candidates, indexed.candidates);
  require(reference.ties > 0 && reference.ties == indexed.ties &&
              reference.candidates.size() == indexed.candidates.size(),
          "equal-score first-winner witness is vacuous or reordered");
  Choices expected;
  for (std::size_t c = 0; c < reference.candidates.size(); ++c) {
    const auto &x = reference.candidates[c];
    const auto &y = indexed.candidates[c];
    require(equal(x.key, y.key) && equal(x.deduction, y.deduction) &&
                same(x.contribution, y.contribution),
            "equal-score candidate stream changed order or deduction bits");
    if (!Max::is_zero(x.contribution) && a.chart.contains(x.key))
      expected.observe(x.key, x.deduction, x.contribution, Max::better);
  }
  require(expected.values().size() == indexed.choices.values().size(),
          "first-winner choices do not match the retained candidate parents");
  for (const auto &[key, value] : expected.values()) {
    const auto *actual = indexed.choices.find(key);
    require(actual && same(actual->contribution, value.contribution) &&
                equal(actual->deduction, value.deduction),
            "equal-score choice did not preserve the first candidate");
  }
  first_winner_ties += reference.ties;
}
}  // namespace

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    Input branching{"branching",
        "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
        "....(...)((....)..)((.(.((...))))(.(....)((...))))"};
    auto forbidden = branching;
    forbidden.name = "branching-x";
    forbidden.scaffold[0] = 'x';
    forbidden.scaffold[forbidden.scaffold.find('.', forbidden.scaffold.size() / 2)] = 'x';
    forbidden.scaffold[forbidden.scaffold.rfind('.')] = 'x';
    const std::vector<Input> inputs{nested_input(), branching, forbidden};
    for (const auto &input : inputs)
      for (std::size_t width : {1, 3, 20})
        for (auto closure : {sc::WmbClosurePolicy::RequireWeakClosure,
                             sc::WmbClosurePolicy::LegacyCompatible})
          for (auto ranking : {sc::beam::RankingMode::Baseline, sc::beam::RankingMode::Outer})
            check(input, width, closure, ranking);
    check(branching, 100, sc::WmbClosurePolicy::RequireWeakClosure, sc::beam::RankingMode::Outer);
    for (const auto &input : {branching, forbidden})
      for (std::size_t width : {1, 3, 20})
        for (auto closure : {sc::WmbClosurePolicy::RequireWeakClosure,
                             sc::WmbClosurePolicy::LegacyCompatible})
          direct_guide_contract(input, width, closure);
    first_winner_contract();
    for (auto closure : {sc::WmbClosurePolicy::RequireWeakClosure,
                         sc::WmbClosurePolicy::LegacyCompatible}) {
      aggregate_lookup_contract<Sum>(branching, closure);
      aggregate_lookup_contract<Max>(branching, closure);
    }
    for (const auto count : aggregate_hits)
      require(count > 0, "one-probe aggregate ownership witness is vacuous");
    require(retained_be > 0 && retained_vp > 0, "BE fallback/VP primary beam coverage is vacuous");
    require(guide_cases == 24 && guide_links > 0 && guide_misses > 0,
            "completed-guide index coverage is vacuous");
    std::cout << "scfg2_beam_numeric_lookup=ok cases=" << cases
              << " compared_items=" << compared_items << " compared_choices=" << compared_choices
              << " first_winner_ties=" << first_winner_ties
              << " guide_cases=" << guide_cases << " guide_links=" << guide_links
              << " guide_misses=" << guide_misses
              << " compared_candidates=" << compared_candidates
              << " contribution_queries=" << contribution_queries << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_border_event_plan.hh"
#include "scfg2/beam/vp_case3_nonleft_last_child_index.hh"
#include "scfg2/beam/vp_case12_last_child_index.hh"
#include "scfg2/beam/be_right_tail_retained_core_index.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace {

using namespace cparty::scfg2;

std::string pairing_sequence(const std::string &scaffold,
                             std::mt19937 *rng = nullptr) {
  static const char bases[] = {'A', 'C', 'G', 'U'};
  std::string sequence(scaffold.size(), 'A');
  for (std::size_t pos = 0; pos < scaffold.size(); ++pos) {
    if (scaffold[pos] == '(') sequence[pos] = 'G';
    if (scaffold[pos] == ')') sequence[pos] = 'C';
    if (scaffold[pos] == '.' && rng != nullptr) {
      sequence[pos] = bases[(*rng)() % 4];
    }
  }
  return sequence;
}

// Unfiltered endpoint path: every (ancestor, child) step to the root.
std::vector<beam::ScaffoldBranchingPathEntry> raw_path(
    cand_pos_t endpoint,
    const replay::WFinalExactInsideContext &ctx) {
  std::vector<beam::ScaffoldBranchingPathEntry> path;
  if (!ctx.is_unpaired(endpoint)) return path;
  cand_pos_t child = ctx.parent_index(endpoint);
  if (child <= 0 || ctx.pair_partner(child) <= endpoint) return path;
  for (cand_pos_t ancestor = ctx.parent_index(child); ancestor >= 0;
       ancestor = ctx.parent_index(child)) {
    path.push_back({ancestor, child});
    if (ancestor == 0) break;
    child = ancestor;
  }
  return path;
}

std::vector<std::size_t> pair_child_counts(
    cand_pos_t n,
    const replay::WFinalExactInsideContext &ctx) {
  std::vector<std::size_t> counts(static_cast<std::size_t>(n + 1), 0);
  for (cand_pos_t pos = 1; pos <= n; ++pos) {
    if (ctx.pair_partner(pos) > pos) {
      ++counts[static_cast<std::size_t>(ctx.parent_index(pos))];
    }
  }
  return counts;
}

// Sum of LeftNonempty CASE3 sources over every recognized left WI child.
std::size_t left_query_sources(
    cand_pos_t n,
    const beam::ExactVpCase3BorderEventChildIndex &index,
    const replay::WFinalExactInsideContext &ctx) {
  std::size_t total = 0;
  for (cand_pos_t i = 1; i <= n; ++i) {
    if (!ctx.is_unpaired(i) || ctx.parent_index(i) <= 0) continue;
    const cand_pos_t p = ctx.pair_partner(ctx.parent_index(i));
    const deductions::ItemKey left{grammar::NonTerminal::WI, i + 1, p - 1,
                                   -1, -1};
    if (!index.recognizes_completed_wi(left, ctx)) continue;
    const auto query = index.query_completed_wi(left, ctx);
    assert(query.border_mismatch_count == 0);
    total += query.sources.size();
  }
  return total;
}

std::string random_block(std::mt19937 &rng, int level) {
  const int stem = 1 + static_cast<int>(rng() % 8);
  std::string inner;
  const int children = level < 3 ? static_cast<int>(rng() % 3) : 0;
  inner.append(1 + rng() % 3, '.');
  for (int child = 0; child < children; ++child) {
    inner += random_block(rng, level + 1);
    inner.append(rng() % 3, '.');
  }
  if (children == 0) inner.append(2, '.');
  return std::string(stem, '(') + inner + std::string(stem, ')');
}

std::string random_scaffold(std::mt19937 &rng) {
  std::string scaffold(rng() % 3, '.');
  const int blocks = 1 + static_cast<int>(rng() % 3);
  for (int block = 0; block < blocks; ++block) {
    scaffold += random_block(rng, 0);
    scaffold.append(rng() % 3, '.');
  }
  return scaffold;
}

void check_single_stem_has_no_case3_lca() {
  constexpr std::size_t depth = 40;
  std::string scaffold(depth, '(');
  scaffold.append(40, '.');
  scaffold.append(depth, ')');
  std::string sequence = pairing_sequence(scaffold);
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto left = beam::ExactVpCase3BorderEventChildIndex::build(n, ctx);
  const auto nonleft = beam::VpCase3NonleftChildEventIndex::build(n, ctx);
  assert(left.path_entry_count() == 0);
  assert(nonleft.path_entry_count() == 0);
  assert(beam::ExactVpCase3BorderEventPlan::discover(n, ctx).source_count() ==
         0);
}

void check_cap_bounds_deep_branching() {
  std::string scaffold = "(....)";
  for (int level = 0; level < 40; ++level) {
    scaffold = "((....)" + scaffold + ")";
  }
  std::string sequence = pairing_sequence(scaffold);
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto exact_left = beam::ExactVpCase3BorderEventChildIndex::build(n, ctx);
  const auto exact_nonleft = beam::VpCase3NonleftChildEventIndex::build(n, ctx);
  const auto capped_left =
      beam::ExactVpCase3BorderEventChildIndex::build(n, ctx, 3);
  const auto capped_nonleft =
      beam::VpCase3NonleftChildEventIndex::build(n, ctx, 3);
  assert(exact_left.path_entry_count() > capped_left.path_entry_count());
  assert(exact_nonleft.path_entry_count() >
         capped_nonleft.path_entry_count());
  assert(capped_left.path_entry_count() <= 3 * sequence.size());
  assert(capped_nonleft.path_entry_count() <= 3 * sequence.size());
  std::size_t expected_truncated = 0;
  for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
    assert(capped_left.paths_for_endpoint(endpoint).size() <= 3);
    if (exact_left.paths_for_endpoint(endpoint).size() > 3) {
      ++expected_truncated;
    }
  }
  assert(exact_left.truncated_endpoint_count() == 0);
  assert(exact_nonleft.truncated_endpoint_count() == 0);
  assert(expected_truncated > 0);
  assert(capped_left.truncated_endpoint_count() == expected_truncated);
  assert(capped_nonleft.truncated_endpoint_count() == expected_truncated);
  std::cout << "deep_branching exact_left_paths="
            << exact_left.path_entry_count()
            << " capped_left_paths=" << capped_left.path_entry_count() << '\n';
}

// Two kissing hairpins whose stems are longer than the default cap. The
// exterior loop is the only branching ancestor, 35 pairs above each loop.
void check_long_stem_kissing_hairpins_survive_default_cap() {
  constexpr std::size_t stem = 35;
  const std::string hairpin =
      std::string(stem, '(') + "......" + std::string(stem, ')');
  std::string scaffold = hairpin + hairpin;
  std::string sequence = pairing_sequence(scaffold);
  for (std::size_t pos = stem; pos < stem + 6; ++pos) {
    sequence[pos] = 'G';
    sequence[hairpin.size() + pos] = 'C';
  }
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  constexpr std::size_t default_cap = 30;
  assert(replay::WFinalVpWmbpOnlineWidths::uniform(1).max_case3_ancestor_depth ==
         default_cap);

  const auto loop_endpoint = static_cast<cand_pos_t>(stem + 2);
  assert(raw_path(loop_endpoint, ctx).size() == stem);

  const auto exact_left = beam::ExactVpCase3BorderEventChildIndex::build(n, ctx);
  const auto capped_left =
      beam::ExactVpCase3BorderEventChildIndex::build(n, ctx, default_cap);
  assert(capped_left.paths_for_endpoint(loop_endpoint).size() == 1);
  assert(capped_left.truncated_endpoint_count() == 0);
  const std::size_t oracle =
      beam::ExactVpCase3BorderEventPlan::discover(n, ctx).source_count();
  const std::size_t exact_sources = left_query_sources(n, exact_left, ctx);
  const std::size_t capped_sources = left_query_sources(n, capped_left, ctx);
  assert(oracle > 0);
  assert(exact_sources == oracle);
  assert(capped_sources == oracle);

  const auto exact_nonleft = beam::VpCase3NonleftChildEventIndex::build(n, ctx);
  const auto capped_nonleft =
      beam::VpCase3NonleftChildEventIndex::build(n, ctx, default_cap);
  const auto right_open = static_cast<cand_pos_t>(hairpin.size() + 1);
  const auto exact_prefixes =
      exact_nonleft.query_empty_middle_prefixes(right_open, ctx);
  const auto capped_prefixes =
      capped_nonleft.query_empty_middle_prefixes(right_open, ctx);
  assert(!exact_prefixes.empty());
  assert(capped_prefixes.size() == exact_prefixes.size());
  std::cout << "long_stem_kissing left_sources=" << capped_sources
            << " nonleft_prefixes=" << capped_prefixes.size() << '\n';
}

void check_random_scaffolds_match_filtered_raw_paths() {
  std::mt19937 rng(20260924);
  std::size_t checked_sources = 0;
  for (int trial = 0; trial < 40; ++trial) {
    std::string scaffold = random_scaffold(rng);
    std::string sequence = pairing_sequence(scaffold, &rng);
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    replay::WFinalExactInsideContext ctx(owner);
    const auto n = static_cast<cand_pos_t>(sequence.size());
    const auto counts = pair_child_counts(n, ctx);
    const auto left = beam::ExactVpCase3BorderEventChildIndex::build(n, ctx);
    const auto nonleft = beam::VpCase3NonleftChildEventIndex::build(n, ctx);
    const auto all_empty = beam::VpCase12DirectEmptyLeftQuery::build(n, ctx);
    for (std::size_t width : {0, 1, 2, 5}) {
      const auto finite_empty = beam::VpCase12DirectEmptyLeftQuery::build(n, ctx, width);
      for (cand_pos_t j = 1; j <= n; ++j) {
        std::vector<beam::VpCase12PrefixGeometry> expected;
        std::size_t case1 = 0;
        for (const auto &prefix : all_empty.query_activation(j, ctx)) {
          if (width == 0) continue;
          if (prefix.rule == grammar::RuleId::VP_WI_CASE1 && case1++ >= width) continue;
          expected.push_back(prefix);
        }
        assert(finite_empty.query_activation(j, ctx) == expected);
      }
    }
    std::size_t max_branching = 0;
    for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
      std::vector<beam::ScaffoldBranchingPathEntry> filtered;
      for (const auto &step : raw_path(endpoint, ctx)) {
        if (counts[static_cast<std::size_t>(step.ancestor)] >= 2) {
          filtered.push_back(step);
        }
      }
      const auto &compressed = left.paths_for_endpoint(endpoint);
      assert(compressed.size() == filtered.size());
      for (std::size_t k = 0; k < filtered.size(); ++k) {
        assert(compressed[k].ancestor == filtered[k].ancestor);
        assert(compressed[k].child == filtered[k].child);
      }
      max_branching = std::max(max_branching, filtered.size());
    }
    assert(nonleft.path_entry_count() == left.path_entry_count());
    const std::size_t oracle =
        beam::ExactVpCase3BorderEventPlan::discover(n, ctx).source_count();
    assert(left_query_sources(n, left, ctx) == oracle);
    assert(left.truncated_endpoint_count() == 0);
    if (max_branching > 0) {
      const auto tight = beam::ExactVpCase3BorderEventChildIndex::build(
          n, ctx, max_branching);
      assert(tight.truncated_endpoint_count() == 0);
      assert(left_query_sources(n, tight, ctx) == oracle);
      const auto short_by_one = beam::ExactVpCase3BorderEventChildIndex::build(
          n, ctx, max_branching - 1);
      assert(short_by_one.truncated_endpoint_count() > 0);
    }
    checked_sources += oracle;
  }
  assert(checked_sources > 0);
  std::cout << "random_scaffolds oracle_sources=" << checked_sources << '\n';
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  check_single_stem_has_no_case3_lca();
  check_cap_bounds_deep_branching();
  check_long_stem_kissing_hairpins_survive_default_cap();
  check_random_scaffolds_match_filtered_raw_paths();

  std::string depth_scaffold =
      std::string(3, '(') + std::string(10, '.') + std::string(3, ')');
  std::string depth_sequence = pairing_sequence(depth_scaffold);
  cparty::scfg2::weights::RuntimeSnapshot depth_owner(
      depth_sequence, depth_scaffold);
  replay::WFinalExactInsideContext depth_ctx(depth_owner);
  const auto depth_n = static_cast<cand_pos_t>(depth_sequence.size());
  const auto capped_case12_empty =
      beam::VpCase12DirectEmptyLeftQuery::build(depth_n, depth_ctx, 3);
  const auto be_parent_index =
      beam::BeRightTailRetainedCoreIndex::build(depth_n, depth_ctx);
  assert(capped_case12_empty.case1_summary_entry_count() <=
         3 * depth_sequence.size());
  assert(be_parent_index.scaffold_parent_link_count() <= depth_sequence.size());
  const auto width_one = replay::WFinalVpWmbpOnlineWidths::uniform(1);
  const auto width_two = replay::WFinalVpWmbpOnlineWidths::uniform(2);
  assert(width_one.max_case3_ancestor_depth == 30);
  assert(width_two.max_case3_ancestor_depth == 30);

  std::string sibling_scaffold = "(";
  for (int child = 0; child < 20; ++child) {
    sibling_scaffold += "(....)";
  }
  sibling_scaffold += ")";
  std::string sibling_sequence = pairing_sequence(sibling_scaffold);
  cparty::scfg2::weights::RuntimeSnapshot sibling_owner(
      sibling_sequence, sibling_scaffold);
  replay::WFinalExactInsideContext sibling_ctx(sibling_owner);
  const auto sibling_n = static_cast<cand_pos_t>(sibling_sequence.size());
  const auto sibling_index =
      beam::ExactVpCase3BorderEventChildIndex::build(
          sibling_n, sibling_ctx);
  const deductions::ItemKey first_child_wi{
      grammar::NonTerminal::WI, 4, 6, -1, -1};
  assert(sibling_index.recognizes_completed_wi(
      first_child_wi, sibling_ctx));
  const auto unlimited_query = sibling_index.query_completed_wi(
      first_child_wi, sibling_ctx);
  assert(unlimited_query.endpoint_pair_count > 5);
  std::cout << "unlimited_endpoint_probes="
            << unlimited_query.endpoint_pair_count << '\n';
  return EXIT_SUCCESS;
}

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case12_last_child_index.hh"
#include "scfg2/beam/vp_case12_prefix_interval.hh"
#include "scfg2/beam/vp_case3_factor.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

bool close(double lhs, double rhs) {
  const double scale = std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  return std::fabs(lhs - rhs) <= 1e-10 * scale;
}

template <class Semiring>
struct CheckResult {
  std::array<std::size_t, 2> positive{};
  cparty::scfg2::beam::VpCase12LastChildStats stats{};
  cparty::scfg2::beam::VpCase12PrefixPlanStats prefix_stats{};
  cparty::scfg2::beam::VpCase12LeftEventStats left_event_stats{};
};

template <class Semiring>
CheckResult<Semiring> check_index(cparty::scfg2::weights::RuntimeSnapshot &owner) {
  using namespace cparty::scfg2;
  replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const auto schedule =
      engine::build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  const auto chart = replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      [](const Deduction &) { return true; },
      [](const ItemKey &, const Deduction &,
         typename Semiring::value_type, const auto &, double) {});
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);
  const auto exact_index =
      beam::ExactVpCase12LastChildIndex<
          replay::WFinalExactInsideContext>::build(n, ctx);
  const auto scaffold_index =
      beam::ScaffoldVpCase12LastChildIndex<
          replay::WFinalExactInsideContext>::build(n, ctx);
  const auto prefix_plan =
      beam::ExactVpCase12PrefixPlan<
          replay::WFinalExactInsideContext>::build(n, ctx);
  const auto scaffold_prefix_plan =
      beam::ScaffoldVpCase12PrefixPlan<
          replay::WFinalExactInsideContext>::build(n, ctx);
  if (scaffold_prefix_plan.prefixes() != prefix_plan.prefixes()) {
    std::cerr << "CASE1/2 scaffold prefix mismatch exact="
              << prefix_plan.stats().prefixes[0] << ','
              << prefix_plan.stats().prefixes[1] << " scaffold="
              << scaffold_prefix_plan.stats().prefixes[0] << ','
              << scaffold_prefix_plan.stats().prefixes[1] << '\n';
    std::size_t printed = 0;
    for (const auto &prefix : scaffold_prefix_plan.prefixes()) {
      if (std::binary_search(prefix_plan.prefixes().begin(),
                             prefix_plan.prefixes().end(), prefix)) {
        continue;
      }
      std::cerr << "extra prefix rule=" << grammar::rule_name(prefix.rule)
                << " i=" << prefix.i << " p=" << prefix.p
                << " q=" << prefix.q << " target="
                << prefix.target_begin << ',' << prefix.target_end << '\n';
      if (++printed == 12) break;
    }
  }
  assert(scaffold_prefix_plan.prefixes() == prefix_plan.prefixes());
  const auto left_event_index =
      beam::VpCase12LeftEventIndex<
          replay::WFinalExactInsideContext>::build(n, ctx);
  const auto direct_nonempty_query =
      beam::VpCase12DirectNonemptyLeftQuery::build(n, ctx);
  const auto direct_empty_query =
      beam::VpCase12DirectEmptyLeftQuery::build(n, ctx);
  const auto direct_empty_width1 =
      beam::VpCase12DirectEmptyLeftQuery::build(n, ctx, 1);
  const auto &exact_stats = exact_index.stats();
  const auto &scaffold_stats = scaffold_index.stats();
  if (scaffold_stats.sources != exact_stats.sources) {
    std::cerr << "CASE1/2 scaffold source mismatch exact="
              << exact_stats.sources[0] << ',' << exact_stats.sources[1]
              << " scaffold=" << scaffold_stats.sources[0] << ','
              << scaffold_stats.sources[1] << " candidates="
              << scaffold_stats.endpoint_candidates << '\n';
    using SourceId =
        std::tuple<grammar::RuleId, cand_pos_t, cand_pos_t>;
    const auto collect = [&](const auto &source_index) {
      std::set<SourceId> ids;
      for (cand_pos_t right = 1; right <= n; ++right) {
        for (const auto &source : source_index.zero_child_sources(right)) {
          ids.emplace(source.rule, source.parent.i, source.parent.j);
        }
        for (cand_pos_t left = 1; left <= right; ++left) {
          const ItemKey child{
              grammar::NonTerminal::WI, left, right, -1, -1};
          for (const auto &source :
               source_index.sources_for_last_child(child)) {
            ids.emplace(source.rule, source.parent.i, source.parent.j);
          }
        }
      }
      return ids;
    };
    const auto exact_ids = collect(exact_index);
    const auto scaffold_ids = collect(scaffold_index);
    std::size_t printed = 0;
    for (const auto &id : exact_ids) {
      if (scaffold_ids.count(id) != 0) continue;
      std::cerr << "missing rule=" << grammar::rule_name(std::get<0>(id))
                << " parent=" << std::get<1>(id) << ',' << std::get<2>(id)
                << '\n';
      if (++printed == 8) break;
    }
  }
  assert(scaffold_stats.sources == exact_stats.sources);
  assert(scaffold_stats.child_counts == exact_stats.child_counts);
  assert(scaffold_stats.left_empty == exact_stats.left_empty);
  assert(scaffold_stats.right_empty == exact_stats.right_empty);
  assert(scaffold_stats.zero_child == exact_stats.zero_child);
  assert(scaffold_stats.child_indexed == exact_stats.child_indexed);
  assert(scaffold_stats.border_mismatches == 0);
  assert(scaffold_stats.duplicate_sources == 0);

  using Totals = std::unordered_map<ItemKey,
                                    typename Semiring::value_type,
                                    deductions::ItemKeyHash,
                                    deductions::ItemKeyEq>;
  Totals driven;
  Totals prefix_driven;
  using SourceId =
      std::tuple<grammar::RuleId,
                 cand_pos_t,
                 cand_pos_t,
                 cand_pos_t,
                 cand_pos_t>;
  std::set<SourceId> prefix_source_ids;
  std::set<std::pair<cand_pos_t, cand_pos_t>> driven_seen;
  std::array<std::size_t, 2> driven_positive{};
  auto accumulate = [&](const Deduction &source) {
    const double local = weight.weight_for(source);
    const auto term = replay::runtime_deduction_contribution<Semiring>(
        chart, source, local);
    auto [found, inserted] = driven.emplace(source.parent, Semiring::zero());
    found->second = Semiring::plus(found->second, term);
    driven_seen.emplace(source.parent.i, source.parent.j);
    if (term != Semiring::zero()) {
      ++driven_positive[beam::vp_case12_rule_index(source.rule)];
    }
  };

  for (cand_pos_t j = 1; j <= n; ++j) {
    for (const auto &source : scaffold_index.zero_child_sources(j)) {
      accumulate(source);
    }
    for (cand_pos_t i = 1; i <= j; ++i) {
      const ItemKey child{grammar::NonTerminal::WI, i, j, -1, -1};
      if (chart.get(child) == Semiring::zero()) continue;
      for (const auto &source : scaffold_index.sources_for_last_child(child)) {
        accumulate(source);
      }
    }
  }

  std::array<std::size_t, 2> prefix_positive{};
  std::set<beam::VpCase12PrefixGeometry> event_prefixes;
  for (cand_pos_t column = 1; column <= n; ++column) {
    const auto &indexed_empty =
        left_event_index.query_empty_activation(column);
    const auto direct_empty = direct_empty_query.query_activation(column, ctx);
    if (direct_empty != indexed_empty) {
      std::cerr << "CASE1/2 direct empty query mismatch column=" << column
                << " indexed=" << indexed_empty.size()
                << " direct=" << direct_empty.size() << '\n';
    }
    assert(direct_empty == indexed_empty);
    const auto width1_empty =
        direct_empty_width1.query_activation(column, ctx);
    assert(width1_empty.size() <= 2);
    for (const auto &prefix : width1_empty) {
      assert(std::binary_search(
          direct_empty.begin(), direct_empty.end(), prefix));
    }
    for (const auto &prefix : direct_empty) {
      assert(prefix.target_begin == column);
      assert(prefix.i + 1 > prefix.p - 1);
      assert(event_prefixes.insert(prefix).second);
    }
    for (cand_pos_t start = 1; start <= column; ++start) {
      const ItemKey child{
          grammar::NonTerminal::WI, start, column, -1, -1};
      if (chart.get(child) == Semiring::zero()) continue;
      const auto &indexed = left_event_index.query_completed_left(child);
      const auto direct = direct_nonempty_query.query_completed_left(
          child, ctx);
      if (direct != indexed) {
        std::cerr << "CASE1/2 direct nonempty query mismatch child="
                  << start << ',' << column << " indexed="
                  << indexed.size() << " direct=" << direct.size() << '\n';
      }
      assert(direct == indexed);
      const auto width1 = direct_nonempty_query.query_completed_left(
          child, ctx, 1);
      assert(width1.size() <= 2);
      for (const auto &prefix : width1) {
        assert(std::binary_search(direct.begin(), direct.end(), prefix));
      }
      for (const auto &prefix : direct) {
        assert(prefix.i + 1 == start);
        assert(prefix.p - 1 == column);
        assert(event_prefixes.insert(prefix).second);
      }
    }
  }
  assert(event_prefixes ==
         std::set<beam::VpCase12PrefixGeometry>(
             prefix_plan.prefixes().begin(), prefix_plan.prefixes().end()));

  const auto agenda_nonempty_query =
      beam::VpCase12DirectNonemptyLeftQuery::build(n, ctx);
  const auto agenda_empty_query =
      beam::VpCase12DirectEmptyLeftQuery::build(n, ctx);
  beam::VpCase12PrefixIntervalAgenda agenda;
  Totals agenda_driven;
  std::set<SourceId> agenda_source_ids;
  for (cand_pos_t j = 1; j <= n; ++j) {
    std::vector<cand_pos_t> right_starts;
    if (j > 1) {
      for (cand_pos_t start = 1; start <= j - 1; ++start) {
        if (chart.get(ItemKey{
                grammar::NonTerminal::WI, start, j - 1, -1, -1}) !=
            Semiring::zero()) {
          right_starts.push_back(start);
        }
      }
    }
    agenda.prepare_column(j, right_starts, agenda_empty_query, ctx);
    for (cand_pos_t i = 1; i <= j; ++i) {
      const ItemKey parent{grammar::NonTerminal::VP, i, j, -1, -1};
      for (const auto &source : agenda.take_for(parent)) {
        const SourceId id{source.rule,
                          source.parent.i,
                          source.parent.j,
                          source.split.p,
                          source.split.q};
        assert(agenda_source_ids.insert(id).second);
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, source, weight.weight_for(source));
        auto [found, inserted] =
            agenda_driven.emplace(parent, Semiring::zero());
        found->second = Semiring::plus(found->second, term);
      }
    }
    for (cand_pos_t start = 1; start <= j; ++start) {
      const ItemKey left{
          grammar::NonTerminal::WI, start, j, -1, -1};
      const auto value = chart.get(left);
      if (value == Semiring::zero()) continue;
      agenda.observe_left_wi(left, value, agenda_nonempty_query, ctx);
    }
  }
  assert(agenda.stats().border_mismatches == 0);
  assert(agenda.stats().duplicate_sources == 0);
  for (const auto &prefix : prefix_plan.prefixes()) {
    for (cand_pos_t j = prefix.target_begin; j <= prefix.target_end; ++j) {
      const ItemKey parent{
          grammar::NonTerminal::VP, prefix.i, j, -1, -1};
      for (const auto &source :
           deductions::generate_exact_basic_deductions(
               parent, prefix.rule, ctx)) {
        if (!engine::legacy_span_schedule_supports(source)) continue;
        const double local = weight.weight_for(source);
        if (local == 0.0) continue;
        const bool matches = source.split.p == prefix.p &&
            (prefix.rule == grammar::RuleId::VP_WI_CASE2 ||
             source.split.q == prefix.q);
        if (!matches) continue;
        const SourceId id{source.rule,
                          source.parent.i,
                          source.parent.j,
                          source.split.p,
                          source.split.q};
        if (!prefix_source_ids.insert(id).second) {
          std::cerr << "duplicate CASE1/2 prefix source rule="
                    << grammar::rule_name(source.rule) << " parent="
                    << source.parent.i << ',' << source.parent.j
                    << " split=" << source.split.p << ',' << source.split.q
                    << '\n';
          assert(false);
        }
        const double left = beam::case3_wi_span_total(
            chart, source.parent.i + 1, source.split.p - 1);
        const double right = beam::case3_wi_span_total(
            chart, source.split.q + 1, source.parent.j - 1);
        const auto reconstructed = Semiring::times(
            Semiring::times(left, right), local);
        const auto direct = replay::runtime_deduction_contribution<Semiring>(
            chart, source, local);
        if (!close(reconstructed, direct)) {
          std::cerr << "CASE1/2 prefix factor mismatch rule="
                    << grammar::rule_name(source.rule) << " parent="
                    << source.parent.i << ',' << source.parent.j
                    << " split=" << source.split.p << ',' << source.split.q
                    << " reconstructed=" << reconstructed
                    << " direct=" << direct << '\n';
          assert(false);
        }
        auto [found, inserted] =
            prefix_driven.emplace(parent, Semiring::zero());
        found->second = Semiring::plus(found->second, reconstructed);
        if (direct != Semiring::zero()) {
          ++prefix_positive[beam::vp_case12_rule_index(source.rule)];
        }
      }
    }
  }

  std::set<std::pair<cand_pos_t, cand_pos_t>> exact_seen;
  std::set<SourceId> exact_source_ids;
  std::array<std::size_t, 2> exact_positive{};
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = 1; i <= j; ++i) {
      const ItemKey parent{grammar::NonTerminal::VP, i, j, -1, -1};
      auto exact = Semiring::zero();
      for (const auto rule : {grammar::RuleId::VP_WI_CASE1,
                              grammar::RuleId::VP_WI_CASE2}) {
        for (const auto &source :
             deductions::generate_exact_basic_deductions(parent, rule, ctx)) {
          if (!engine::legacy_span_schedule_supports(source)) continue;
          const double local = weight.weight_for(source);
          if (local == 0.0) continue;
          exact_seen.emplace(parent.i, parent.j);
          exact_source_ids.emplace(source.rule,
                                   source.parent.i,
                                   source.parent.j,
                                   source.split.p,
                                   source.split.q);
          const auto term = replay::runtime_deduction_contribution<Semiring>(
              chart, source, local);
          exact = Semiring::plus(exact, term);
          if (term != Semiring::zero()) {
            ++exact_positive[beam::vp_case12_rule_index(rule)];
          }
        }
      }
      const auto found = driven.find(parent);
      const auto indexed = found == driven.end()
          ? Semiring::zero()
          : found->second;
      if (!close(exact, indexed)) {
        std::cerr << "CASE1/2 mismatch parent=" << i << ',' << j
                  << " exact=" << exact << " indexed=" << indexed << '\n';
        assert(false);
      }
      const auto prefix_found = prefix_driven.find(parent);
      const auto prefix_indexed = prefix_found == prefix_driven.end()
          ? Semiring::zero()
          : prefix_found->second;
      if (!close(exact, prefix_indexed)) {
        std::cerr << "CASE1/2 prefix aggregate mismatch parent=" << i << ','
                  << j << " exact=" << exact
                  << " prefix=" << prefix_indexed << '\n';
        assert(false);
      }
      const auto agenda_found = agenda_driven.find(parent);
      const auto agenda_indexed = agenda_found == agenda_driven.end()
          ? Semiring::zero()
          : agenda_found->second;
      if (!close(exact, agenda_indexed)) {
        std::cerr << "CASE1/2 agenda aggregate mismatch parent=" << i << ','
                  << j << " exact=" << exact
                  << " agenda=" << agenda_indexed << '\n';
        assert(false);
      }
    }
  }
  assert(driven_seen == exact_seen);
  assert(driven_positive == exact_positive);
  assert(prefix_source_ids == exact_source_ids);
  assert(agenda_source_ids == exact_source_ids);
  assert(prefix_positive == exact_positive);
  assert(prefix_plan.stats().sources == exact_stats.sources);
  return CheckResult<Semiring>{
      exact_positive,
      exact_stats,
      prefix_plan.stats(),
      left_event_index.stats()};
}

}  // namespace

int main(int argc, char **argv) {
  assert(argc == 1 || argc == 3 || argc == 4);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 4 ? argv[3] : "params/rna_Turner04.par",
         1);
  std::string sequence = argc >= 3
      ? argv[1]
      : "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  std::string scaffold = argc >= 3
      ? argv[2]
      : "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  const auto inside =
      check_index<cparty::scfg2::semirings::SumProductSemiring>(owner);
  const auto viterbi =
      check_index<cparty::scfg2::semirings::MaxProductSemiring>(owner);
  if (argc == 1) {
    assert(inside.positive[0] > 0 && inside.positive[1] > 0);
    assert(viterbi.positive[0] > 0 && viterbi.positive[1] > 0);
  }
  assert(inside.stats.sources == viterbi.stats.sources);
  assert(inside.stats.child_counts == viterbi.stats.child_counts);
  assert(inside.prefix_stats.sources == viterbi.prefix_stats.sources);
  assert(inside.prefix_stats.prefixes == viterbi.prefix_stats.prefixes);
  assert(inside.left_event_stats.nonempty ==
         viterbi.left_event_stats.nonempty);
  assert(inside.left_event_stats.empty == viterbi.left_event_stats.empty);
  std::cout << "case1_sources\tcase2_sources\t"
               "case1_prefixes\tcase2_prefixes\t"
               "case1_left_event_nonempty\tcase1_left_event_empty\t"
               "case2_left_event_nonempty\tcase2_left_event_empty\t"
               "case1_child0\tcase1_child1\tcase1_child2\t"
               "case2_child0\tcase2_child1\tcase2_child2\t"
               "case1_left_empty\tcase1_right_empty\t"
               "case2_left_empty\tcase2_right_empty\t"
               "case1_inside_positive\tcase2_inside_positive\t"
               "case1_viterbi_positive\tcase2_viterbi_positive\n"
            << inside.stats.sources[0] << '\t' << inside.stats.sources[1]
            << '\t' << inside.prefix_stats.prefixes[0] << '\t'
            << inside.prefix_stats.prefixes[1]
            << '\t' << inside.left_event_stats.nonempty[0] << '\t'
            << inside.left_event_stats.empty[0] << '\t'
            << inside.left_event_stats.nonempty[1] << '\t'
            << inside.left_event_stats.empty[1]
            << '\t' << inside.stats.child_counts[0][0] << '\t'
            << inside.stats.child_counts[0][1] << '\t'
            << inside.stats.child_counts[0][2] << '\t'
            << inside.stats.child_counts[1][0] << '\t'
            << inside.stats.child_counts[1][1] << '\t'
            << inside.stats.child_counts[1][2] << '\t'
            << inside.stats.left_empty[0] << '\t'
            << inside.stats.right_empty[0] << '\t'
            << inside.stats.left_empty[1] << '\t'
            << inside.stats.right_empty[1] << '\t'
            << inside.positive[0] << '\t' << inside.positive[1] << '\t'
            << viterbi.positive[0] << '\t' << viterbi.positive[1] << '\n';
  return EXIT_SUCCESS;
}

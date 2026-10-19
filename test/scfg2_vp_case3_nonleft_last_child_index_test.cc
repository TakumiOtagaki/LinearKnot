#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_nonleft_last_child_index.hh"
#include "scfg2/beam/vp_case3_factor.hh"
#include "scfg2/beam/vp_case3_nonleft_prefix_interval.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
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
  std::size_t boundary_sources = 0;
  std::size_t other_sources = 0;
  std::size_t zero_child_sources = 0;
  std::size_t positive_sources = 0;
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
  const auto index =
      beam::ExactVpCase3NonleftLastChildIndex<
          replay::WFinalExactInsideContext>::build(n, ctx);
  const auto scaffold_index =
      beam::ScaffoldVpCase3NonleftLastChildIndex<
          replay::WFinalExactInsideContext>::build(n, ctx);
  const auto &exact_stats = index.stats();
  const auto &stats = scaffold_index.stats();
  assert(stats.boundary_sources == exact_stats.boundary_sources);
  assert(stats.other_sources == exact_stats.other_sources);
  assert(stats.zero_child_sources == exact_stats.zero_child_sources);
  assert(stats.child_indexed_sources == exact_stats.child_indexed_sources);
  assert(stats.border_mismatches == 0);
  assert(stats.duplicate_sources == 0);
  const auto event_index = beam::VpCase3NonleftChildEventIndex::build(n, ctx);

  using Totals = std::unordered_map<ItemKey,
                                    typename Semiring::value_type,
                                    deductions::ItemKeyHash,
                                    deductions::ItemKeyEq>;
  Totals driven;
  Totals prefix_driven;
  Totals interval_driven;
  std::set<std::pair<cand_pos_t, cand_pos_t>> event_seen;
  std::set<std::pair<cand_pos_t, cand_pos_t>> prefix_seen;
  std::set<std::pair<cand_pos_t, cand_pos_t>> interval_seen;
  std::size_t driven_positive = 0;
  auto accumulate = [&](const Deduction &source) {
    const auto term = replay::runtime_deduction_contribution<Semiring>(
        chart, source, weight.weight_for(source));
    auto [found, inserted] = driven.emplace(source.parent, Semiring::zero());
    found->second = Semiring::plus(found->second, term);
    driven_positive += term != Semiring::zero() ? 1 : 0;
  };

  for (cand_pos_t j = 1; j <= n; ++j) {
    const auto zero_query = event_index.query_zero_child_column(j, ctx);
    assert(zero_query.border_mismatches == 0);
    assert(zero_query.duplicate_sources == 0);
    for (const auto &source : zero_query.sources) {
      assert(event_seen.emplace(source.parent.i, source.parent.j).second);
      accumulate(source);
    }
    for (cand_pos_t i = 1; i <= j; ++i) {
      const ItemKey child{grammar::NonTerminal::WI, i, j, -1, -1};
      if (chart.get(child) == Semiring::zero()) continue;
      const auto query = event_index.query_completed_wi(child, ctx);
      assert(query.border_mismatches == 0);
      assert(query.duplicate_sources == 0);
      for (const auto &source : query.sources) {
        assert(event_seen.emplace(source.parent.i, source.parent.j).second);
        accumulate(source);
      }
    }
  }

  beam::VpCase3NonleftPrefixIntervalAgenda interval_agenda;
  for (cand_pos_t j = 1; j <= n; ++j) {
    if (j > 1) {
      for (cand_pos_t start = 1; start <= j - 1; ++start) {
        const ItemKey middle{
            grammar::NonTerminal::WI, start, j - 1, -1, -1};
        interval_agenda.observe_middle_wi(
            middle, chart.get(middle), event_index, ctx);
      }
    }
    interval_agenda.prepare_column(j, event_index, ctx);
    for (cand_pos_t i = 1; i <= j; ++i) {
      const ItemKey parent{
          grammar::NonTerminal::VP, i, j, -1, -1};
      for (const auto &source : interval_agenda.take_for(parent)) {
        assert(interval_seen.emplace(source.parent.i, source.parent.j).second);
        const double local = weight.weight_for(source);
        const double term = replay::runtime_deduction_contribution<Semiring>(
            chart, source, local);
        auto [found, inserted] =
            interval_driven.emplace(source.parent, Semiring::zero());
        found->second = Semiring::plus(found->second, term);
      }
    }
  }

  std::set<beam::VpCase3NonleftPrefixGeometry> prefixes;
  for (cand_pos_t k = 1; k <= n; ++k) {
    for (const auto &prefix :
         event_index.query_empty_middle_prefixes(k, ctx)) {
      prefixes.insert(prefix);
    }
    if (k <= 1) continue;
    for (cand_pos_t start = 1; start <= k - 1; ++start) {
      const ItemKey middle{
          grammar::NonTerminal::WI, start, k - 1, -1, -1};
      for (const auto &prefix :
           event_index.query_middle_prefixes(middle, ctx)) {
        prefixes.insert(prefix);
      }
    }
  }
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (const auto &prefix : prefixes) {
      if (j < prefix.k || j > prefix.target_close) continue;
      const ItemKey parent{
          grammar::NonTerminal::VP, prefix.i, j, -1, -1};
      for (const auto &source :
           deductions::generate_exact_vp_case3_deductions_for_item(
               parent, ctx)) {
        const auto partition = beam::vp_case3_partition(source);
        if (!partition.has_value() ||
            *partition == beam::VpCase3Partition::LeftNonempty ||
            source.split.p != prefix.p || source.split.q != prefix.q ||
            source.split.k != prefix.k ||
            !engine::legacy_span_schedule_supports(source)) {
          continue;
        }
        const double local = weight.weight_for(source);
        if (local == 0.0) continue;
        assert(prefix_seen.emplace(source.parent.i, source.parent.j).second);
        const double middle = beam::case3_wi_span_total(
            chart, prefix.q + 1, prefix.k - 1);
        const double right = beam::case3_wi_span_total(
            chart, source.split.l + 1, source.parent.j - 1);
        const double reconstructed = local * middle * right;
        const double direct = replay::runtime_deduction_contribution<Semiring>(
            chart, source, local);
        assert(close(reconstructed, direct));
        auto [found, inserted] =
            prefix_driven.emplace(source.parent, Semiring::zero());
        found->second = Semiring::plus(found->second, reconstructed);
      }
    }
  }

  std::size_t exact_positive = 0;
  std::set<std::pair<cand_pos_t, cand_pos_t>> exact_seen;
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = 1; i <= j; ++i) {
      const ItemKey parent{grammar::NonTerminal::VP, i, j, -1, -1};
      auto exact = Semiring::zero();
      for (const auto &source :
           deductions::generate_exact_vp_case3_deductions_for_item(
               parent, ctx)) {
        const auto partition = beam::vp_case3_partition(source);
        if (!partition.has_value() ||
            *partition == beam::VpCase3Partition::LeftNonempty ||
            !engine::legacy_span_schedule_supports(source)) {
          continue;
        }
        const double local = weight.weight_for(source);
        if (local == 0.0) continue;
        exact_seen.emplace(source.parent.i, source.parent.j);
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, source, local);
        exact = Semiring::plus(exact, term);
        exact_positive += term != Semiring::zero() ? 1 : 0;
      }
      const auto found = driven.find(parent);
      const auto indexed = found == driven.end()
          ? Semiring::zero()
          : found->second;
      const auto prefix_found = prefix_driven.find(parent);
      const auto prefix_value = prefix_found == prefix_driven.end()
          ? Semiring::zero()
          : prefix_found->second;
      const auto interval_found = interval_driven.find(parent);
      const auto interval_value = interval_found == interval_driven.end()
          ? Semiring::zero()
          : interval_found->second;
      if (!close(exact, indexed)) {
        std::cerr << "CASE3 non-left mismatch parent=" << i << ',' << j
                  << " exact=" << exact << " indexed=" << indexed << '\n';
        for (const auto &source :
             deductions::generate_exact_vp_case3_deductions_for_item(
                 parent, ctx)) {
          const auto partition = beam::vp_case3_partition(source);
          if (!partition.has_value() ||
              *partition == beam::VpCase3Partition::LeftNonempty) continue;
          std::cerr << "split p,q,k,l=" << source.split.p << ','
                    << source.split.q << ',' << source.split.k << ','
                    << source.split.l << " children=" << source.child_count;
          for (std::size_t child = 0; child < source.child_count; ++child) {
            const auto &item = source.children[child].item;
            std::cerr << " [" << item.i << ',' << item.j << ']';
          }
          std::cerr << '\n';
        }
        assert(false);
      }
      if (!close(exact, prefix_value)) {
        std::cerr << "CASE3 non-left prefix mismatch parent=" << i << ',' << j
                  << " exact=" << exact << " prefix=" << prefix_value
                  << '\n';
        assert(false);
      }
      if (!close(exact, interval_value)) {
        std::cerr << "CASE3 non-left interval mismatch parent=" << i << ','
                  << j << " exact=" << exact
                  << " interval=" << interval_value << '\n';
        assert(false);
      }
    }
  }
  if (event_seen != exact_seen) {
    std::cerr << "CASE3 non-left identity mismatch exact="
              << exact_seen.size() << " event=" << event_seen.size() << '\n';
    for (const auto &parent : exact_seen) {
      if (event_seen.count(parent) == 0) {
        std::cerr << "missing parent=" << parent.first << ','
                  << parent.second << '\n';
      }
    }
    for (const auto &parent : event_seen) {
      if (exact_seen.count(parent) == 0) {
        std::cerr << "extra parent=" << parent.first << ','
                  << parent.second << '\n';
      }
    }
    assert(false);
  }
  assert(prefix_seen == exact_seen);
  assert(interval_seen == exact_seen);
  assert(interval_agenda.stats().border_mismatches == 0);
  assert(interval_agenda.stats().duplicate_prefixes == 0);
  assert(interval_agenda.stats().duplicate_sources == 0);
  assert(driven_positive == exact_positive);
  assert(stats.zero_child_sources + stats.child_indexed_sources ==
         stats.boundary_sources + stats.other_sources);
  return CheckResult<Semiring>{stats.boundary_sources,
                               stats.other_sources,
                               stats.zero_child_sources,
                               exact_positive};
}

// A long left hairpin has Theta(n) ineligible endpoints before its one
// eligible prefix. Repeated finite queries must not rescan those endpoints.
void check_bounded_prefix_query() {
  using namespace cparty::scfg2;
  struct CountingContext {
    replay::WFinalExactInsideContext &ctx;
    mutable std::size_t reads = 0;
    cand_pos_t parent_index(cand_pos_t pos) const { ++reads; return ctx.parent_index(pos); }
    cand_pos_t pair_partner(cand_pos_t pos) const { ++reads; return ctx.pair_partner(pos); }
    bool is_unpaired(cand_pos_t pos) const { return ctx.is_unpaired(pos); }
  };
  for (const cand_pos_t n : {100, 1000, 4000}) {
    std::string sequence(n, 'A');
    sequence[n - 6] = sequence[n - 1] = 'U';
    std::string scaffold = "(" + std::string(n - 7, '.') + ")(...)";
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    replay::WFinalExactInsideContext ctx(owner);
    CountingContext observed{ctx};
    const auto index = beam::VpCase3NonleftChildEventIndex::build(n, observed, 30);
    observed.reads = 0;
    for (int q = 0; q < 64; ++q) {
      const auto prefixes = index.query_empty_middle_prefixes(n - 4, observed, 1);
      assert(prefixes.size() == 1);
      assert(prefixes[0].i == n - 6 && prefixes[0].p == n - 5);
      assert(prefixes[0].k == n - 4 && prefixes[0].target_close == n);
    }
    assert(observed.reads <= 64 * 8);
  }
}

}  // namespace

int main(int argc, char **argv) {
  assert(argc == 1 || argc == 3 || argc == 4);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 4 ? argv[3] : "params/rna_Turner04.par",
         1);
  check_bounded_prefix_query();
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
  assert(inside.boundary_sources == viterbi.boundary_sources);
  assert(inside.other_sources == viterbi.other_sources);
  assert(inside.zero_child_sources == viterbi.zero_child_sources);
  if (argc == 1) {
    assert(inside.boundary_sources > 0);
    assert(inside.other_sources > 0);
    assert(inside.positive_sources > 0);
    assert(viterbi.positive_sources > 0);
  }
  std::cout << "boundary\tother\tzero_child\tinside_positive\t"
               "viterbi_positive\n"
            << inside.boundary_sources << '\t' << inside.other_sources << '\t'
            << inside.zero_child_sources << '\t' << inside.positive_sources
            << '\t' << viterbi.positive_sources << '\n';
  return EXIT_SUCCESS;
}

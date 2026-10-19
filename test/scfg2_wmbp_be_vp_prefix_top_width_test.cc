#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/wmbp_be_vp_prefix_events.hh"
#include "scfg2/beam/wmbp_be_vp_split_child_index.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

void check_typed_cache_weights() {
  using namespace cparty::scfg2;
  using S = semirings::LogSumExpSemiring;
  engine::Chart<double> chart(S::zero());
  Deduction d;
  d.rule = grammar::RuleId::WMBP_SPLIT_BE_WI_VP;
  d.children[0].item = {grammar::NonTerminal::BE, 1, 6, 2, 5};
  d.children[1].item = {grammar::NonTerminal::WI, 7, 8, -1, -1};
  d.children[2].item = {grammar::NonTerminal::VP, 9, 12, -1, -1};
  d.child_count = 3;
  chart.set(d.children[0].item, std::log(2.0));
  chart.set(d.children[1].item, std::log(3.0));
  chart.set(d.children[2].item, std::log(4.0));
  const auto check = [&](auto local, double log_local) {
    beam::ExactWmbpBeVpPrefixCache cache;
    const auto prefix = cache.materialize_prefix<S>(d, chart, local);
    const auto source = cache.reconstruct_source<S>(d, chart, local);
    assert(prefix && source && std::isfinite(*prefix) && std::isfinite(*source));
    assert(close(*prefix, log_local + std::log(6.0)));
    assert(close(*source, log_local + std::log(24.0)));
  };
  check(0.5, std::log(0.5));
  check(weights::AnyLocalWeight{0.5}, std::log(0.5));
  // Underflowing raw projection must stay finite in the log semiring.
  check(weights::LogWeight{-1000.0}, -1000.0);
  check(weights::AnyLocalWeight{weights::LogWeight{-1000.0}}, -1000.0);
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

struct SelectionRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t candidates = 0;
  std::size_t selected = 0;
  std::size_t pruned = 0;
  std::size_t rejected_sources = 0;
  std::size_t cache_prefixes = 0;
  std::size_t cache_reuses = 0;
  std::size_t inconsistent_reuses = 0;
  std::size_t combined_l_visits = 0;
  std::size_t separate_l_visits = 0;
  std::size_t indexed_l_visits = 0;
  std::size_t index_bucket_visits = 0;
};

template <class Semiring>
SelectionRun run_width(const std::string &sequence,
                       const std::string &scaffold,
                       std::size_t width) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto plan = cparty::scfg2::beam::
      ExactWmbpBeVpPrefixParentPlan::build(schedule, ctx);
  const auto split_index = cparty::scfg2::beam::
      ExactWmbpBeVpSplitChildIndex::build(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  cparty::scfg2::beam::ExactWmbpBeVpPrefixCache cache;
  cparty::scfg2::beam::ExactWmbpBeVpPrefixCache oracle_cache;
  cparty::scfg2::beam::TopWidthWmbpBeVpPrefixSelection selection(width);
  std::size_t rejected_sources = 0;
  std::size_t combined_l_visits = 0;
  std::size_t separate_l_visits = 0;
  std::size_t indexed_l_visits = 0;
  std::size_t index_bucket_visits = 0;
  std::optional<ItemKey> pending_parent;
  std::vector<Deduction> pending_deductions;
  auto chart = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_provider_item_observers_and_contribution<
          Semiring>(
          schedule,
          ctx,
          [&](const ItemKey &item) {
            if (item.nonterminal !=
                cparty::scfg2::grammar::NonTerminal::WMBP) {
              return cparty::scfg2::deductions::
                  generate_exact_basic_deductions_for_item(item, ctx);
            }
            assert(pending_parent.has_value());
            assert(cparty::scfg2::deductions::ItemKeyEq{}(
                *pending_parent, item));
            pending_parent.reset();
            return std::exchange(pending_deductions, {});
          },
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          },
          [&](const ItemKey &item, const auto &current_chart) {
            if (item.nonterminal !=
                cparty::scfg2::grammar::NonTerminal::WMBP) {
              return;
            }
            const auto planned =
                plan.template materialize_parent<Semiring>(
                    item, current_chart, oracle_cache);
            auto discovery = cparty::scfg2::beam::
                discover_wmbp_be_vp_prefixes_for_parent<Semiring>(
                    item, current_chart, ctx, cache);
            const auto recomputed = cparty::scfg2::beam::
                score_wmbp_be_vp_prefixes_from_deductions<Semiring>(
                    discovery.deductions, current_chart, ctx);
            assert(cparty::scfg2::beam::same_scored_wmbp_be_vp_prefixes(
                recomputed, discovery.scored));
            const bool active = cparty::scfg2::deductions::
                exact_wmbp_item_active(item, ctx);
            if (!active) {
              assert(discovery.deductions.empty());
              assert(discovery.scored.empty());
            }
            auto shared_combined = cparty::scfg2::deductions::
                generate_exact_wmbp_shared_be_vp_deductions(
                    item, ctx, true, true);
            auto shared_wmbp = cparty::scfg2::deductions::
                generate_exact_wmbp_shared_be_vp_deductions(
                    item, ctx, true, false);
            auto shared_wmbw = cparty::scfg2::deductions::
                generate_exact_wmbp_shared_be_vp_deductions(
                    item, ctx, false, true);
            shared_wmbp.insert(
                shared_wmbp.end(),
                std::make_move_iterator(shared_wmbw.begin()),
                std::make_move_iterator(shared_wmbw.end()));
            assert(shared_combined.size() == shared_wmbp.size());
            for (std::size_t idx = 0; idx < shared_wmbp.size(); ++idx) {
              assert(same_deduction(shared_combined[idx], shared_wmbp[idx]));
            }
            const auto family_combined = cparty::scfg2::deductions::
                generate_exact_wmbp_be_vp_family_deductions(
                    item, ctx, true, true, true);
            const auto family_shared = cparty::scfg2::deductions::
                generate_exact_wmbp_be_vp_family_deductions(
                    item, ctx, true, true, false);
            const auto family_wi = cparty::scfg2::deductions::
                generate_exact_wmbp_be_vp_family_deductions(
                    item, ctx, false, false, true);
            const auto same_list = [&](const auto &lhs, const auto &rhs) {
              assert(lhs.size() == rhs.size());
              for (std::size_t idx = 0; idx < lhs.size(); ++idx) {
                assert(same_deduction(lhs[idx], rhs[idx]));
              }
            };
            same_list(family_combined.wmbp, family_shared.wmbp);
            same_list(family_combined.wmbw, family_shared.wmbw);
            same_list(family_combined.wi, family_wi.wi);
            assert(family_combined.l_visits <=
                   family_shared.l_visits + family_wi.l_visits);
            const auto coordinate_scan = cparty::scfg2::deductions::
                scan_exact_wmbp_be_vp_split_geometries(
                    item, ctx, true, true);
            const auto indexed_scan = cparty::scfg2::beam::
                scan_exact_wmbp_be_vp_split_geometries_from_index(
                    item, ctx, split_index, true, true);
            assert(cparty::scfg2::beam::
                same_wmbp_be_vp_split_geometry_scans(
                    coordinate_scan, indexed_scan.geometry));
            assert(indexed_scan.geometry.l_visits <=
                   coordinate_scan.l_visits);
            if (active) {
              combined_l_visits += family_combined.l_visits;
              separate_l_visits +=
                  family_shared.l_visits + family_wi.l_visits;
              indexed_l_visits += indexed_scan.geometry.l_visits;
              index_bucket_visits += indexed_scan.parent_bucket_visits;
            }
            assert(cparty::scfg2::beam::same_scored_wmbp_be_vp_prefixes(
                planned, discovery.scored));
            selection.template commit_parent<Semiring>(
                item, std::move(discovery.scored));
            pending_parent = item;
            pending_deductions = std::move(discovery.deductions);
          },
          [](const ItemKey &, const auto &) {},
          [&](const Deduction &deduction,
              const auto &current_chart,
              double local) {
            if (!selection.allows(deduction)) {
              ++rejected_sources;
              return Semiring::zero();
            }
            const auto reconstructed =
                cache.template reconstruct_source<Semiring>(
                    deduction, current_chart, local);
            if (reconstructed.has_value()) return *reconstructed;
            return cparty::scfg2::replay::
                runtime_deduction_contribution<Semiring>(
                    current_chart, deduction, local);
          });
  assert(!pending_parent.has_value());
  assert(selection.respects_width());
  assert(selection.candidate_count() == plan.prefix_frame_count());
  return SelectionRun{std::move(chart),
                      selection.candidate_count(),
                      selection.selected_count(),
                      selection.pruned_count(),
                      rejected_sources,
                      cache.prefix_count(),
                      cache.reuse_count(),
                      cache.inconsistent_reuse_count(),
                      combined_l_visits,
                      separate_l_visits,
                      indexed_l_visits,
                      index_bucket_visits};
}

template <class Semiring>
cparty::scfg2::engine::Chart<double> run_direct(
    const std::string &sequence,
    const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  return cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      [](const Deduction &) { return true; },
      [](const ItemKey &, const Deduction &, double, const auto &, double) {});
}

struct WidthChecks {
  SelectionRun keep_all;
  SelectionRun width_one;
};

template <class Semiring>
WidthChecks check_widths(const std::string &sequence,
                         const std::string &scaffold) {
  const auto keep_all = run_width<Semiring>(
      sequence, scaffold, std::numeric_limits<std::size_t>::max());
  const auto width_one = run_width<Semiring>(sequence, scaffold, 1);
  const auto direct = run_direct<Semiring>(sequence, scaffold);
  assert(keep_all.candidates == keep_all.selected);
  assert(keep_all.pruned == 0);
  assert(keep_all.rejected_sources == 0);
  assert(keep_all.inconsistent_reuses == 0);
  assert(width_one.inconsistent_reuses == 0);
  for (const auto &[item, exact] : direct.values()) {
    assert(close(keep_all.chart.get(item), exact));
    const double approximate = width_one.chart.get(item);
    assert(approximate <= exact || close(approximate, exact));
  }
  return {keep_all, width_one};
}

}  // namespace

int main(int argc, char **argv) {
  check_typed_cache_weights();
  assert(argc == 1 || argc == 3 || argc == 4);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 4 ? argv[3] : "params/rna_Turner04.par",
         1);
  const std::string sequence = argc >= 3
      ? argv[1]
      : "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  const std::string scaffold = argc >= 3
      ? argv[2]
      : "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  assert(sequence.size() == scaffold.size());

  const auto inside =
      check_widths<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi =
      check_widths<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  assert(inside.keep_all.candidates == viterbi.keep_all.candidates);
  assert(inside.width_one.selected == viterbi.width_one.selected);
  assert(inside.width_one.pruned == viterbi.width_one.pruned);
  assert(inside.keep_all.combined_l_visits ==
         viterbi.keep_all.combined_l_visits);
  assert(inside.keep_all.separate_l_visits ==
         viterbi.keep_all.separate_l_visits);
  assert(inside.keep_all.combined_l_visits <=
         inside.keep_all.separate_l_visits);
  assert(inside.keep_all.indexed_l_visits ==
         viterbi.keep_all.indexed_l_visits);
  assert(inside.keep_all.indexed_l_visits <=
         inside.keep_all.combined_l_visits);
  if (argc == 1) assert(inside.width_one.pruned > 0);

  std::cout << "n\tcandidates\tselected\tpruned\trejected_sources"
               "\tcache_prefixes\tcache_reuses\tcombined_l_visits"
               "\tseparate_l_visits\tindexed_l_visits"
               "\tindex_bucket_visits\tmonotone\n"
            << sequence.size() << '\t' << inside.width_one.candidates << '\t'
            << inside.width_one.selected << '\t'
            << inside.width_one.pruned << '\t'
            << inside.width_one.rejected_sources << '\t'
            << inside.width_one.cache_prefixes << '\t'
            << inside.width_one.cache_reuses << '\t'
            << inside.width_one.combined_l_visits << '\t'
            << inside.width_one.separate_l_visits << '\t'
            << inside.width_one.indexed_l_visits << '\t'
            << inside.width_one.index_bucket_visits << "\t1\n";
}

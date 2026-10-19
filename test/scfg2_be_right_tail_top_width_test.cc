#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/be_right_tail_selection.hh"
#include "scfg2/beam/be_right_tail_child_index.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

struct TailRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t candidates = 0;
  std::size_t selected = 0;
  std::size_t pruned = 0;
  std::size_t planned_sources = 0;
  std::size_t rejected_sources = 0;
  std::size_t scaffold_pairs = 0;
  std::size_t core_outer_visits = 0;
  std::size_t anchor_visits = 0;
  std::size_t coordinate_core_visits = 0;
  std::size_t coordinate_anchor_visits = 0;
};

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

template <class Semiring>
TailRun run_width(const std::string &sequence,
                  const std::string &scaffold,
                  std::size_t width,
                  bool online) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto plan =
      cparty::scfg2::beam::ExactBeRightTailColumnPlan::build(schedule, ctx);
  cparty::scfg2::beam::ExactBeRightTailChildIndex child_index(plan);
  auto scaffold_index =
      cparty::scfg2::beam::ExactBeRightTailChildIndex::build_from_scaffold(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto coordinate_oracle = cparty::scfg2::beam::
      ExactBeRightTailChildIndex::build_from_scaffold_coordinate_oracle(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto pair_tree = cparty::scfg2::beam::ExactBeRightTailChildIndex::
      build_from_scaffold_pair_tree(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  assert(scaffold_index.has_same_keys(child_index));
  assert(scaffold_index.has_same_keys(coordinate_oracle));
  assert(scaffold_index.has_same_keys(pair_tree));
  cparty::scfg2::beam::ExactBeFrameCache cache;
  cparty::scfg2::beam::TopWidthBeRightTailSelection selection(width);
  cand_pos_t committed_j = -1;
  std::size_t rejected_sources = 0;
  auto chart = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_item_observers_and_contribution<
          Semiring>(
          schedule,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          },
          [&](const ItemKey &item, const auto &current_chart) {
            if (item.j == committed_j) return;
            committed_j = item.j;
            auto candidates = online
                ? scaffold_index.take_column(item.j)
                : plan.template materialize_column<Semiring>(
                      item.j, current_chart, cache);
            selection.template commit_column<Semiring>(
                item.j, std::move(candidates));
          },
          [&](const ItemKey &item, const auto &current_chart) {
            if (online) {
              scaffold_index.template observe_completed_child<Semiring>(
                  item, current_chart, cache);
            }
          },
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
  assert(selection.respects_width());
  assert(selection.candidate_count() == plan.tail_count());
  if (online) {
    assert(scaffold_index.record_count() == plan.tail_count());
    assert(scaffold_index.materialized_count() == plan.tail_count());
    assert(scaffold_index.committed_count() == plan.tail_count());
    assert(scaffold_index.pending_empty());
  }
  return TailRun{std::move(chart),
                 selection.candidate_count(),
                 selection.selected_count(),
                 selection.pruned_count(),
                 plan.source_count(),
                 rejected_sources,
                 pair_tree.scaffold_pair_count(),
                 pair_tree.core_outer_visit_count(),
                 pair_tree.anchor_visit_count(),
                 coordinate_oracle.coordinate_core_visit_count(),
                 coordinate_oracle.coordinate_anchor_visit_count()};
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

template <class Semiring>
TailRun check_width(const std::string &sequence,
                    const std::string &scaffold) {
  const auto keep_all = run_width<Semiring>(
      sequence,
      scaffold,
      std::numeric_limits<std::size_t>::max(),
      false);
  const auto width_one = run_width<Semiring>(sequence, scaffold, 1, false);
  const auto online_keep_all = run_width<Semiring>(
      sequence,
      scaffold,
      std::numeric_limits<std::size_t>::max(),
      true);
  const auto online_width_one =
      run_width<Semiring>(sequence, scaffold, 1, true);
  const auto direct = run_direct<Semiring>(sequence, scaffold);
  assert(keep_all.candidates == keep_all.selected);
  assert(keep_all.pruned == 0);
  assert(keep_all.rejected_sources == 0);
  for (const auto &[item, exact] : direct.values()) {
    assert(close(keep_all.chart.get(item), exact));
  }
  assert(width_one.candidates == keep_all.candidates);
  assert(width_one.selected <= width_one.candidates);
  if (width_one.pruned > 0) assert(width_one.rejected_sources > 0);
  for (const auto &[item, exact] : keep_all.chart.values()) {
    const double approximate = width_one.chart.get(item);
    assert(approximate <= exact || close(approximate, exact));
    assert(close(online_keep_all.chart.get(item), exact));
    assert(close(online_width_one.chart.get(item), approximate));
  }
  assert(online_width_one.candidates == width_one.candidates);
  assert(online_width_one.selected == width_one.selected);
  assert(online_width_one.pruned == width_one.pruned);
  assert(online_width_one.rejected_sources == width_one.rejected_sources);
  return online_width_one;
}

}  // namespace

int main(int argc, char **argv) {
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
      check_width<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi =
      check_width<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  assert(inside.candidates == viterbi.candidates);
  assert(inside.selected == viterbi.selected);
  assert(inside.pruned == viterbi.pruned);
  assert(inside.planned_sources == viterbi.planned_sources);
  assert(inside.scaffold_pairs == viterbi.scaffold_pairs);
  assert(inside.core_outer_visits == viterbi.core_outer_visits);
  assert(inside.anchor_visits == viterbi.anchor_visits);
  assert(inside.coordinate_core_visits == viterbi.coordinate_core_visits);
  assert(inside.coordinate_anchor_visits ==
         viterbi.coordinate_anchor_visits);
  if (argc == 1) assert(inside.pruned > 0);
  std::cout << "n\tcandidates\tselected\tpruned\tplanned_sources"
               "\trejected_sources\tscaffold_pairs\tcore_outer_visits"
               "\tanchor_visits\tcoordinate_core_visits"
               "\tcoordinate_anchor_visits\tmonotone\n"
            << sequence.size() << '\t' << inside.candidates << '\t'
            << inside.selected << '\t' << inside.pruned << '\t'
            << inside.planned_sources << '\t' << inside.rejected_sources
            << '\t' << inside.scaffold_pairs << '\t'
            << inside.core_outer_visits << '\t' << inside.anchor_visits
            << '\t' << inside.coordinate_core_visits << '\t'
            << inside.coordinate_anchor_visits << "\t1\n";
  return EXIT_SUCCESS;
}

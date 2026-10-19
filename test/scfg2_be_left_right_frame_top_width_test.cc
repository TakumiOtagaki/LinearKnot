#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/be_left_right_frame_selection.hh"
#include "scfg2/beam/be_left_right_frame_events.hh"
#include "scfg2/beam/be_right_tail_child_index.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

struct FrameRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t tail_candidates = 0;
  std::size_t tail_selected = 0;
  std::size_t tail_pruned = 0;
  std::size_t frame_candidates = 0;
  std::size_t frame_selected = 0;
  std::size_t frame_pruned = 0;
  std::size_t tail_rejected_sources = 0;
  std::size_t frame_rejected_sources = 0;
};

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

template <class Semiring>
FrameRun run_widths(const std::string &sequence,
                    const std::string &scaffold,
                    std::size_t tail_width,
                    std::size_t frame_width) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  auto tail_index =
      cparty::scfg2::beam::ExactBeRightTailChildIndex::build_from_scaffold(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto frame_plan = cparty::scfg2::beam::
      ExactBeRightShellFrameColumnPlan::build(schedule, ctx);
  cparty::scfg2::beam::ExactBeFrameCache cache;
  cparty::scfg2::beam::TopWidthBeRightTailSelection tail_selection(
      tail_width);
  cparty::scfg2::beam::TopWidthBeRightShellFrameSelection frame_selection(
      frame_width);
  cand_pos_t committed_j = -1;
  std::size_t tail_rejected_sources = 0;
  std::size_t frame_rejected_sources = 0;
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
            const auto tail_candidates = tail_index.take_column(item.j);
            tail_selection.template commit_column<Semiring>(
                item.j, tail_candidates);
            const auto planned = frame_plan.template score_column<Semiring>(
                item.j,
                current_chart,
                cache,
                [&](const auto &key) { return tail_selection.allows(key); });
            const auto event_scored = cparty::scfg2::beam::
                score_be_right_shell_frames_from_tail_events<Semiring>(
                    tail_candidates,
                    current_chart,
                    ctx,
                    [&](const auto &key) {
                      return tail_selection.allows(key);
                    });
            assert(cparty::scfg2::beam::same_scored_be_right_shell_frames(
                planned, event_scored));
            frame_selection.template commit_column<Semiring>(
                item.j, event_scored);
          },
          [&](const ItemKey &item, const auto &current_chart) {
            tail_index.template observe_completed_child<Semiring>(
                item, current_chart, cache);
          },
          [&](const Deduction &deduction,
              const auto &current_chart,
              double local) {
            if (!tail_selection.allows(deduction)) {
              ++tail_rejected_sources;
              return Semiring::zero();
            }
            if (!frame_selection.allows(deduction)) {
              ++frame_rejected_sources;
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
  assert(tail_selection.respects_width());
  assert(frame_selection.respects_width());
  assert(tail_index.pending_empty());
  return FrameRun{std::move(chart),
                  tail_selection.candidate_count(),
                  tail_selection.selected_count(),
                  tail_selection.pruned_count(),
                  frame_selection.candidate_count(),
                  frame_selection.selected_count(),
                  frame_selection.pruned_count(),
                  tail_rejected_sources,
                  frame_rejected_sources};
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
  FrameRun frame_one;
  FrameRun tail_one;
  FrameRun combined_one;
};

template <class Semiring>
WidthChecks check_widths(const std::string &sequence,
                         const std::string &scaffold) {
  constexpr std::size_t all = std::numeric_limits<std::size_t>::max();
  const auto keep_all = run_widths<Semiring>(sequence, scaffold, all, all);
  const auto frame_one = run_widths<Semiring>(sequence, scaffold, all, 1);
  const auto tail_one = run_widths<Semiring>(sequence, scaffold, 1, all);
  const auto combined_one = run_widths<Semiring>(sequence, scaffold, 1, 1);
  const auto direct = run_direct<Semiring>(sequence, scaffold);
  assert(keep_all.tail_candidates == keep_all.tail_selected);
  assert(keep_all.frame_candidates == keep_all.frame_selected);
  assert(keep_all.tail_rejected_sources == 0);
  assert(keep_all.frame_rejected_sources == 0);
  assert(frame_one.frame_candidates == keep_all.frame_candidates);
  assert(tail_one.frame_candidates <= keep_all.frame_candidates);
  assert(combined_one.frame_candidates == tail_one.frame_candidates);
  for (const auto &[item, exact] : direct.values()) {
    assert(close(keep_all.chart.get(item), exact));
    for (const auto *run : {&frame_one, &tail_one, &combined_one}) {
      const double approximate = run->chart.get(item);
      assert(approximate <= exact || close(approximate, exact));
    }
  }
  return WidthChecks{frame_one, tail_one, combined_one};
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
      check_widths<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi =
      check_widths<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  assert(inside.frame_one.frame_candidates ==
         viterbi.frame_one.frame_candidates);
  assert(inside.frame_one.frame_selected ==
         viterbi.frame_one.frame_selected);
  assert(inside.tail_one.tail_selected == viterbi.tail_one.tail_selected);
  assert(inside.tail_one.frame_candidates ==
         viterbi.tail_one.frame_candidates);
  assert(inside.combined_one.frame_candidates ==
         viterbi.combined_one.frame_candidates);
  assert(inside.combined_one.frame_selected ==
         viterbi.combined_one.frame_selected);
  if (argc == 1) {
    assert(inside.frame_one.frame_pruned > 0);
    assert(inside.tail_one.tail_pruned > 0);
  }
  std::cout
      << "n\tframe_candidates\tframe_selected\tframe_pruned"
         "\tframe_rejected_sources\ttail_candidates\ttail_selected"
         "\ttail_pruned\ttail_result_frames\ttail_rejected_sources"
         "\tcombined_frames\tcombined_selected\tcombined_pruned"
         "\tcombined_tail_rejected\tcombined_frame_rejected\tmonotone\n"
      << sequence.size() << '\t' << inside.frame_one.frame_candidates << '\t'
      << inside.frame_one.frame_selected << '\t'
      << inside.frame_one.frame_pruned << '\t'
      << inside.frame_one.frame_rejected_sources << '\t'
      << inside.tail_one.tail_candidates << '\t'
      << inside.tail_one.tail_selected << '\t'
      << inside.tail_one.tail_pruned << '\t'
      << inside.tail_one.frame_candidates << '\t'
      << inside.tail_one.tail_rejected_sources << '\t'
      << inside.combined_one.frame_candidates << '\t'
      << inside.combined_one.frame_selected << '\t'
      << inside.combined_one.frame_pruned << '\t'
      << inside.combined_one.tail_rejected_sources << '\t'
      << inside.combined_one.frame_rejected_sources << "\t1\n";
  return EXIT_SUCCESS;
}

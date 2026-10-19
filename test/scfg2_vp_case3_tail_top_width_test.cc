#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_left_selection.hh"
#include "scfg2/beam/vp_case3_middle_reuse_profile.hh"
#include "scfg2/beam/vp_case3_right_reuse_profile.hh"
#include "scfg2/beam/vp_case3_tail_selection.hh"
#include "scfg2/beam/vp_case3_wi_frontier.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {

struct TailRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t tail_candidates = 0;
  std::size_t tail_selected = 0;
  std::size_t tail_pruned = 0;
  std::size_t left_candidates = 0;
  std::size_t rejected_sources = 0;
  cparty::scfg2::beam::VpCase3RightReuseProfile right_profile{};
  cparty::scfg2::beam::VpCase3MiddleReuseProfile middle_profile{};
};

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

template <class Semiring>
TailRun run_tail_width(const std::string &sequence,
                       const std::string &scaffold,
                       std::size_t tail_width) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(sequence.size());
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  const auto child_index =
      cparty::scfg2::beam::ExactVpCase3BorderEventChildIndex::build(n, ctx);
  cparty::scfg2::beam::PendingVpCase3LeftCandidates pending;
  cparty::scfg2::beam::ExactVpCase3TailCache tails;
  cparty::scfg2::beam::ExactVpCase3LeftCache left;
  cparty::scfg2::beam::TopWidthVpCase3TailSelection tail_selection(
      tail_width);
  cparty::scfg2::beam::TopWidthVpCase3LeftSelection left_selection(
      std::numeric_limits<std::size_t>::max());
  auto producer =
      cparty::scfg2::beam::VpCase3WiProducerFrontier::keep_all();
  cand_pos_t materialized_j = -1;
  std::size_t rejected_sources = 0;
  auto chart = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_item_observers_and_contribution<
          Semiring>(
          schedule,
          ctx,
          [](const auto &) { return true; },
          [](const auto &, const auto &, double, const auto &, double) {},
          [&](const auto &item, const auto &current_chart) {
            if (item.j == materialized_j) return;
            materialized_j = item.j;
            producer.flush_before(
                item.j, current_chart, pending, child_index, ctx);
            pending.prepare_column(item.j, current_chart, tails);
            tail_selection.commit_column(item.j, tails);
            const auto scored = pending.score_prepared_column(
                item.j,
                tails,
                [&](const auto &key) { return tail_selection.allows(key); });
            left_selection.commit_column(item.j, scored);
            (void)left.materialize_scored(
                scored,
                [&](const auto &key) { return left_selection.allows(key); });
          },
          [&](const auto &item, const auto &current_chart) {
            producer.observe_completed_child(
                item, current_chart, child_index, ctx);
          },
          [&](const auto &deduction,
              const auto &current_chart,
              double local) {
            if (!left_selection.allows(deduction)) {
              ++rejected_sources;
              return Semiring::zero();
            }
            const auto reconstructed =
                left.reconstruct_source(deduction, local);
            if (reconstructed.has_value()) return *reconstructed;
            return cparty::scfg2::replay::
                runtime_deduction_contribution<Semiring>(
                    current_chart, deduction, local);
          });
  assert(tail_selection.respects_width());
  assert(tail_selection.record_probe_count() ==
         tail_selection.candidate_count());
  assert(left_selection.selected_count() ==
         left_selection.candidate_count());
  const auto right_profile =
      cparty::scfg2::beam::profile_vp_case3_right_reuse(tails);
  const auto middle_profile =
      cparty::scfg2::beam::profile_vp_case3_middle_reuse(tails);
  assert(tails.right_factor_count() == right_profile.right_view_count);
  assert(tails.right_factor_materialization_count() ==
         right_profile.right_view_count);
  assert(tails.right_factor_reuse_count() + tails.right_factor_count() ==
         tails.tail_count());
  assert(tails.middle_factor_count() == middle_profile.middle_view_count);
  assert(tails.middle_factor_materialization_count() ==
         middle_profile.middle_view_count);
  assert(tails.middle_factor_reuse_count() + tails.middle_factor_count() ==
         tails.tail_count());
  return TailRun{std::move(chart),
                 tail_selection.candidate_count(),
                 tail_selection.selected_count(),
                 tail_selection.pruned_count(),
                 left_selection.candidate_count(),
                 rejected_sources,
                 right_profile,
                 middle_profile};
}

template <class Semiring>
TailRun check_tail_width(const std::string &sequence,
                         const std::string &scaffold) {
  const auto keep_all = run_tail_width<Semiring>(
      sequence, scaffold, std::numeric_limits<std::size_t>::max());
  const auto width_one = run_tail_width<Semiring>(sequence, scaffold, 1);
  assert(keep_all.tail_candidates == keep_all.tail_selected);
  assert(keep_all.tail_pruned == 0);
  assert(keep_all.rejected_sources == 0);
  assert(width_one.tail_candidates == keep_all.tail_candidates);
  assert(width_one.tail_selected <= width_one.tail_candidates);
  assert(width_one.left_candidates <= keep_all.left_candidates);
  if (width_one.tail_pruned > 0) {
    assert(width_one.rejected_sources > 0);
  }
  for (const auto &[item, exact] : keep_all.chart.values()) {
    const double approximate = width_one.chart.get(item);
    assert(approximate <= exact || close(approximate, exact));
  }
  return width_one;
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
      check_tail_width<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi =
      check_tail_width<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  assert(inside.tail_candidates == viterbi.tail_candidates);
  assert(inside.tail_selected == viterbi.tail_selected);
  assert(inside.tail_pruned == viterbi.tail_pruned);
  assert(inside.left_candidates == viterbi.left_candidates);
  assert(inside.right_profile.right_view_count ==
         viterbi.right_profile.right_view_count);
  assert(inside.right_profile.reused_right_view_count ==
         viterbi.right_profile.reused_right_view_count);
  assert(inside.right_profile.maximum_tails_per_right_view ==
         viterbi.right_profile.maximum_tails_per_right_view);
  assert(inside.right_profile.multi_right_view_column_count == 0);
  assert(inside.right_profile.maximum_right_views_per_column <= 1);
  assert(inside.middle_profile.middle_view_count ==
         viterbi.middle_profile.middle_view_count);
  assert(inside.middle_profile.reused_middle_view_count ==
         viterbi.middle_profile.reused_middle_view_count);
  assert(inside.middle_profile.maximum_tails_per_middle_view ==
         viterbi.middle_profile.maximum_tails_per_middle_view);
  if (argc == 1) assert(inside.tail_pruned > 0);
  std::cout << "n\ttail_candidates\ttail_selected\ttail_pruned"
               "\tleft_candidates\trejected_sources\tright_views"
               "\tnonempty_right_views\tempty_right_views\treused_right_views"
               "\tmax_tails_per_right\tright_columns\tmulti_right_columns"
               "\tmax_right_views_per_column\tmiddle_views"
               "\tnonempty_middle_views\tempty_middle_views"
               "\treused_middle_views\tmax_tails_per_middle"
               "\tmiddle_columns\tmulti_middle_columns"
               "\tmax_middle_views_per_column\tmonotone\n"
            << sequence.size() << '\t' << inside.tail_candidates << '\t'
            << inside.tail_selected << '\t' << inside.tail_pruned << '\t'
            << inside.left_candidates << '\t' << inside.rejected_sources
            << '\t' << inside.right_profile.right_view_count << '\t'
            << inside.right_profile.nonempty_right_view_count << '\t'
            << inside.right_profile.empty_right_view_count << '\t'
            << inside.right_profile.reused_right_view_count << '\t'
            << inside.right_profile.maximum_tails_per_right_view << '\t'
            << inside.right_profile.column_count << '\t'
            << inside.right_profile.multi_right_view_column_count << '\t'
            << inside.right_profile.maximum_right_views_per_column << '\t'
            << inside.middle_profile.middle_view_count << '\t'
            << inside.middle_profile.nonempty_middle_view_count << '\t'
            << inside.middle_profile.empty_middle_view_count << '\t'
            << inside.middle_profile.reused_middle_view_count << '\t'
            << inside.middle_profile.maximum_tails_per_middle_view << '\t'
            << inside.middle_profile.column_count << '\t'
            << inside.middle_profile.multi_middle_view_column_count << '\t'
            << inside.middle_profile.maximum_middle_views_per_column
            << "\t1\n";
  return EXIT_SUCCESS;
}

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_left_selection.hh"
#include "scfg2/beam/vp_case3_border_event_plan.hh"
#include "scfg2/beam/vp_case3_left_child_index.hh"
#include "scfg2/beam/vp_case3_tail_selection.hh"
#include "scfg2/beam/vp_case3_wi_frontier.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

struct FiniteRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t candidates = 0;
  std::size_t selected = 0;
  std::size_t materialized = 0;
  std::size_t pruned = 0;
  std::size_t rejected_sources = 0;
  std::size_t producer_candidates = 0;
  std::size_t producer_selected = 0;
  std::size_t producer_pruned = 0;
  std::size_t queried_endpoint_pairs = 0;
  std::size_t queried_sources = 0;
  std::size_t materialized_tails = 0;
  std::size_t tail_candidates = 0;
  std::size_t tail_selected = 0;
  std::size_t tail_pruned = 0;
};

struct WidthChecks {
  FiniteRun left_one;
  FiniteRun tail_one;
  FiniteRun producer_one;
  FiniteRun combined_one;
  FiniteRun all_wi_one;
};

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

template <class Semiring>
FiniteRun run_width(const std::string &sequence,
                    const std::string &scaffold,
                    std::size_t left_width,
                    std::size_t producer_width,
                    std::size_t tail_width,
                    cparty::scfg2::beam::VpCase3WiFrontierScope scope =
                        cparty::scfg2::beam::VpCase3WiFrontierScope::
                            EligibleCase3Producers) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto tail_plan =
      cparty::scfg2::beam::build_vp_case3_tail_column_plan_from_border_events(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto segments =
      cparty::scfg2::beam::build_vp_case3_tail_attachment_segments(tail_plan);
  const cparty::scfg2::beam::ExactVpCase3LeftColumnPlan left_plan(segments);
  const cparty::scfg2::beam::ExactVpCase3LeftChildIndex child_index(segments);
  const auto child_event_index =
      cparty::scfg2::beam::ExactVpCase3BorderEventChildIndex::build(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  cparty::scfg2::beam::PendingVpCase3LeftCandidates pending(left_width);
  cparty::scfg2::beam::ExactVpCase3TailCache tails;
  cparty::scfg2::beam::ExactVpCase3LeftCache left;
  cparty::scfg2::beam::TopWidthVpCase3LeftSelection selection(left_width);
  cparty::scfg2::beam::TopWidthVpCase3TailSelection tail_selection(
      tail_width);
  cparty::scfg2::beam::VpCase3WiProducerFrontier producer(
      producer_width, scope);
  cand_pos_t materialized_j = -1;
  std::size_t rejected_sources = 0;
  auto chart = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_item_observers_and_contribution<
          Semiring>(
          schedule,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {},
          [&](const ItemKey &item, const auto &current_chart) {
            if (item.j == materialized_j) return;
            materialized_j = item.j;
            producer.flush_before(
                item.j, current_chart, pending, child_event_index, ctx);
            pending.prepare_column(item.j, current_chart, tails);
            tail_selection.commit_column(item.j, tails);
            const auto scored = pending.score_prepared_column(
                item.j,
                tails,
                [&](const auto &key) { return tail_selection.allows(key); });
            selection.commit_column(item.j, scored);
            (void)left.materialize_scored(
                scored,
                [&](const cparty::scfg2::beam::VpCase3LeftKey &key) {
                  return selection.allows(key);
                });
          },
          [&](const ItemKey &item, const auto &current_chart) {
            producer.observe_completed_child(
                item, current_chart, child_event_index, ctx);
          },
          [&](const Deduction &deduction,
              const auto &current_chart,
              double local) {
            if (!selection.allows(deduction)) {
              ++rejected_sources;
              return Semiring::zero();
            }
            const auto reconstructed = left.reconstruct_source(deduction, local);
            if (reconstructed.has_value()) return *reconstructed;
            assert(cparty::scfg2::beam::vp_case3_partition(deduction) !=
                   cparty::scfg2::beam::VpCase3Partition::LeftNonempty);
            return cparty::scfg2::replay::
                runtime_deduction_contribution<Semiring>(
                    current_chart, deduction, local);
          });
  assert(selection.respects_width());
  assert(tail_selection.respects_width());
  assert(producer.respects_width());
  assert(pending.border_mismatch_count() == 0);
  assert(pending.duplicate_source_count() == 0);
  if (left_width != std::numeric_limits<std::size_t>::max()) {
    assert(pending.max_active_segment_count() <= left_width);
    assert(pending.tail_materialization_probe_count() <=
           sequence.size() * left_width);
    assert(pending.prepared_score_probe_count() <=
           sequence.size() * left_width);
  }
  if (producer.pruned_count() == 0) {
    assert(pending.observed_attachment_count() ==
           child_index.attachment_count());
    assert(pending.registered_segment_count() ==
           child_index.future_segment_count());
    assert(pending.queried_source_count() == tail_plan.source_count());
    assert(child_index.attachment_count() == left_plan.geometry_count());
  }
  assert(left.candidate_count() == selection.selected_count());
  return FiniteRun{std::move(chart),
                   selection.candidate_count(),
                   selection.selected_count(),
                   left.candidate_count(),
                   selection.pruned_count(),
                   rejected_sources,
                   producer.candidate_count(),
                   producer.selected_count(),
                   producer.pruned_count(),
                   pending.queried_endpoint_pair_count(),
                   pending.queried_source_count(),
                   tails.tail_count(),
                   tail_selection.candidate_count(),
                   tail_selection.selected_count(),
                   tail_selection.pruned_count()};
}

template <class Semiring>
WidthChecks check_width_contract(const std::string &sequence,
                                 const std::string &scaffold) {
  constexpr std::size_t keep_all_width =
      std::numeric_limits<std::size_t>::max();
  const auto keep_all = run_width<Semiring>(
      sequence, scaffold, keep_all_width, keep_all_width, keep_all_width);
  const auto width_one =
      run_width<Semiring>(sequence, scaffold, 1, keep_all_width, keep_all_width);
  const auto tail_one =
      run_width<Semiring>(sequence, scaffold, keep_all_width, keep_all_width, 1);
  const auto producer_one =
      run_width<Semiring>(sequence, scaffold, keep_all_width, 1, keep_all_width);
  const auto combined_one = run_width<Semiring>(sequence, scaffold, 1, 1, 1);
  const auto all_wi_keep_all = run_width<Semiring>(
      sequence,
      scaffold,
      keep_all_width,
      keep_all_width,
      keep_all_width,
      cparty::scfg2::beam::VpCase3WiFrontierScope::AllWiItems);
  const auto all_wi_one = run_width<Semiring>(
      sequence,
      scaffold,
      keep_all_width,
      1,
      keep_all_width,
      cparty::scfg2::beam::VpCase3WiFrontierScope::AllWiItems);
  assert(keep_all.candidates == keep_all.selected);
  assert(keep_all.materialized == keep_all.selected);
  assert(keep_all.pruned == 0);
  assert(keep_all.rejected_sources == 0);
  assert(keep_all.tail_candidates == keep_all.tail_selected);
  assert(keep_all.tail_pruned == 0);
  assert(width_one.candidates <= keep_all.candidates);
  assert(width_one.materialized == width_one.selected);
  if (width_one.pruned > 0) {
    assert(width_one.selected < width_one.candidates);
    assert(width_one.rejected_sources > 0);
  }
  assert(tail_one.tail_candidates == keep_all.tail_candidates);
  if (tail_one.tail_pruned > 0) {
    assert(tail_one.tail_selected < tail_one.tail_candidates);
    assert(tail_one.candidates <= keep_all.candidates);
    assert(tail_one.rejected_sources > 0);
  }
  for (const auto &[item, exact] : keep_all.chart.values()) {
    const double approximate = width_one.chart.get(item);
    assert(approximate <= exact || close(approximate, exact));
    const double producer_approximate = producer_one.chart.get(item);
    assert(producer_approximate <= exact ||
           close(producer_approximate, exact));
    const double tail_approximate = tail_one.chart.get(item);
    assert(tail_approximate <= exact || close(tail_approximate, exact));
    const double combined_approximate = combined_one.chart.get(item);
    assert(combined_approximate <= exact || close(combined_approximate, exact));
    assert(close(all_wi_keep_all.chart.get(item), exact));
    const double all_wi_approximate = all_wi_one.chart.get(item);
    assert(all_wi_approximate <= exact || close(all_wi_approximate, exact));
  }
  assert(keep_all.producer_pruned == 0);
  if (producer_one.producer_pruned > 0) {
    assert(producer_one.candidates <= keep_all.candidates);
  }
  assert(producer_one.queried_endpoint_pairs <=
         keep_all.queried_endpoint_pairs);
  assert(producer_one.queried_sources <= keep_all.queried_sources);
  assert(producer_one.materialized_tails <= keep_all.materialized_tails);
  assert(combined_one.queried_endpoint_pairs <=
         producer_one.queried_endpoint_pairs);
  assert(combined_one.queried_sources <= producer_one.queried_sources);
  assert(combined_one.materialized_tails <= producer_one.materialized_tails);
  assert(combined_one.tail_candidates <= tail_one.tail_candidates);
  assert(combined_one.candidates <= producer_one.candidates);
  assert(combined_one.selected <= combined_one.candidates);
  assert(all_wi_keep_all.rejected_sources == 0);
  assert(all_wi_one.queried_endpoint_pairs <=
         all_wi_keep_all.queried_endpoint_pairs);
  assert(all_wi_one.queried_sources <= all_wi_keep_all.queried_sources);
  return WidthChecks{
      width_one, tail_one, producer_one, combined_one, all_wi_one};
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
      check_width_contract<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi =
      check_width_contract<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  assert(inside.left_one.candidates == viterbi.left_one.candidates);
  assert(inside.left_one.selected == viterbi.left_one.selected);
  assert(inside.left_one.materialized == viterbi.left_one.materialized);
  assert(inside.left_one.pruned == viterbi.left_one.pruned);
  assert(inside.tail_one.tail_candidates ==
         viterbi.tail_one.tail_candidates);
  assert(inside.tail_one.tail_selected == viterbi.tail_one.tail_selected);
  assert(inside.tail_one.tail_pruned == viterbi.tail_one.tail_pruned);
  assert(inside.producer_one.producer_candidates ==
         viterbi.producer_one.producer_candidates);
  assert(inside.producer_one.producer_selected ==
         viterbi.producer_one.producer_selected);
  assert(inside.producer_one.producer_pruned ==
         viterbi.producer_one.producer_pruned);
  assert(inside.producer_one.queried_endpoint_pairs ==
         viterbi.producer_one.queried_endpoint_pairs);
  assert(inside.producer_one.queried_sources ==
         viterbi.producer_one.queried_sources);
  assert(inside.producer_one.materialized_tails ==
         viterbi.producer_one.materialized_tails);
  assert(inside.producer_one.candidates == viterbi.producer_one.candidates);
  assert(inside.combined_one.producer_candidates ==
         viterbi.combined_one.producer_candidates);
  assert(inside.combined_one.producer_selected ==
         viterbi.combined_one.producer_selected);
  assert(inside.combined_one.tail_candidates ==
         viterbi.combined_one.tail_candidates);
  assert(inside.combined_one.tail_selected ==
         viterbi.combined_one.tail_selected);
  assert(inside.combined_one.candidates == viterbi.combined_one.candidates);
  assert(inside.combined_one.selected == viterbi.combined_one.selected);
  assert(inside.all_wi_one.producer_candidates ==
         viterbi.all_wi_one.producer_candidates);
  assert(inside.all_wi_one.producer_selected ==
         viterbi.all_wi_one.producer_selected);
  assert(inside.all_wi_one.queried_sources ==
         viterbi.all_wi_one.queried_sources);
  assert(inside.all_wi_one.candidates == viterbi.all_wi_one.candidates);
  if (argc == 1) {
    assert(inside.left_one.rejected_sources > 0);
    assert(inside.tail_one.tail_pruned > 0);
    assert(inside.producer_one.producer_pruned > 0);
    assert(inside.combined_one.producer_pruned > 0);
    assert(inside.combined_one.rejected_sources > 0);
    assert(inside.all_wi_one.producer_pruned > 0);
  }
  std::cout << "n\tcandidates\tselected\tmaterialized\tpruned\trejected_sources"
               "\tproducer_candidates\tproducer_selected\tproducer_pruned"
               "\tproducer_rejected_sources\tproducer_endpoint_pairs"
               "\tproducer_sources\tproducer_tails"
               "\tproducer_left_candidates\tall_wi_candidates"
               "\tall_wi_selected\tall_wi_sources\tall_wi_left_candidates"
               "\ttail_candidates\ttail_selected\ttail_left_candidates"
               "\tcombined_producer_candidates\tcombined_producer_selected"
               "\tcombined_producer_pruned\tcombined_endpoint_pairs"
               "\tcombined_sources\tcombined_materialized_tails"
               "\tcombined_tail_candidates\tcombined_tail_selected"
               "\tcombined_tail_pruned\tcombined_left_candidates"
               "\tcombined_left_selected\tcombined_left_pruned"
               "\tcombined_rejected_sources"
               "\tmonotone\n"
            << sequence.size() << '\t' << inside.left_one.candidates << '\t'
            << inside.left_one.selected << '\t'
            << inside.left_one.materialized << '\t'
            << inside.left_one.pruned << '\t'
            << inside.left_one.rejected_sources << '\t'
            << inside.producer_one.producer_candidates << '\t'
            << inside.producer_one.producer_selected << '\t'
            << inside.producer_one.producer_pruned << '\t'
            << inside.producer_one.rejected_sources << '\t'
            << inside.producer_one.queried_endpoint_pairs << '\t'
            << inside.producer_one.queried_sources << '\t'
            << inside.producer_one.materialized_tails << '\t'
            << inside.producer_one.candidates << '\t'
            << inside.all_wi_one.producer_candidates << '\t'
            << inside.all_wi_one.producer_selected << '\t'
            << inside.all_wi_one.queried_sources << '\t'
            << inside.all_wi_one.candidates << '\t'
            << inside.tail_one.tail_candidates << '\t'
            << inside.tail_one.tail_selected << '\t'
            << inside.tail_one.candidates << '\t'
            << inside.combined_one.producer_candidates << '\t'
            << inside.combined_one.producer_selected << '\t'
            << inside.combined_one.producer_pruned << '\t'
            << inside.combined_one.queried_endpoint_pairs << '\t'
            << inside.combined_one.queried_sources << '\t'
            << inside.combined_one.materialized_tails << '\t'
            << inside.combined_one.tail_candidates << '\t'
            << inside.combined_one.tail_selected << '\t'
            << inside.combined_one.tail_pruned << '\t'
            << inside.combined_one.candidates << '\t'
            << inside.combined_one.selected << '\t'
            << inside.combined_one.pruned << '\t'
            << inside.combined_one.rejected_sources << "\t1\n";
  return EXIT_SUCCESS;
}

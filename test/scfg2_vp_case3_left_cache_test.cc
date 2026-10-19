#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_left_cache.hh"
#include "scfg2/beam/vp_case3_border_event_plan.hh"
#include "scfg2/beam/vp_case3_left_child_index.hh"
#include "scfg2/beam/vp_case3_left_selection.hh"
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

struct Result {
  std::size_t tails = 0;
  std::size_t left_candidates = 0;
  std::size_t sources = 0;
  std::size_t segments = 0;
  std::size_t wi_children = 0;
  std::size_t max_child_fanout = 0;
  cparty::scfg2::beam::VpCase3LeftChildIndexProfile child_profile{};
};

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

void check_column_parent_start_index() {
  using namespace cparty::scfg2;
  beam::ExactVpCase3LeftCache cache;
  std::vector<beam::VpCase3LeftRecord> records;
  for (cand_pos_t j = 10; j <= 30; j += 10) {
    for (cand_pos_t i = 1; i <= 2; ++i) {
      records.push_back(beam::VpCase3LeftRecord{
          beam::VpCase3LeftKey{
              beam::VpCase3TailKey{j, 3, 4, 5}, i, 6},
          1.0,
          beam::vp_case3_carrier_bit(grammar::NonTerminal::VP)});
    }
  }
  assert(cache.materialize_scored(
             records, [](const beam::VpCase3LeftKey &) { return true; }) ==
         records.size());
  const auto public_starts = cache.vp_parent_starts(20);
  const auto repeated_starts = cache.parent_starts(
      20, grammar::NonTerminal::VP);
  assert((public_starts == std::vector<cand_pos_t>{2, 1}));
  assert(repeated_starts == public_starts);
  assert(cache.parent_start_probe_count() == 4);
}

template <class Semiring>
Result check_left_factorization(const std::string &sequence,
                                const std::string &scaffold) {
  std::string baseline_sequence = sequence;
  std::string baseline_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot baseline_owner(
      baseline_sequence, baseline_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext baseline_ctx(baseline_owner);
  const auto baseline_schedule =
      cparty::scfg2::engine::build_pkfree_basic_exact_schedule(
          static_cast<cand_pos_t>(sequence.size()), baseline_ctx);
  const auto baseline =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          baseline_schedule,
          baseline_ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {});

  std::string factorized_sequence = sequence;
  std::string factorized_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot factorized_owner(
      factorized_sequence, factorized_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext factorized_ctx(
      factorized_owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), factorized_ctx);
  const auto tail_plan =
      cparty::scfg2::beam::build_vp_case3_tail_column_plan_from_border_events(
          static_cast<cand_pos_t>(sequence.size()), factorized_ctx);
  const auto segments =
      cparty::scfg2::beam::build_vp_case3_tail_attachment_segments(tail_plan);
  const cparty::scfg2::beam::ExactVpCase3LeftColumnPlan left_plan(segments);
  const cparty::scfg2::beam::ExactVpCase3LeftChildIndex child_index(segments);
  const auto child_event_index =
      cparty::scfg2::beam::ExactVpCase3BorderEventChildIndex::build(
          static_cast<cand_pos_t>(sequence.size()), factorized_ctx);
  cparty::scfg2::beam::PendingVpCase3LeftCandidates pending;
  cparty::scfg2::beam::ExactVpCase3TailCache tails;
  cparty::scfg2::beam::ExactVpCase3LeftCache left;
  cparty::scfg2::beam::KeepAllVpCase3LeftSelection selection;
  cand_pos_t materialized_j = -1;
  const auto factorized = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_item_observers_and_contribution<
          Semiring>(
          schedule,
          factorized_ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {},
          [&](const ItemKey &item, const auto &chart) {
            if (item.j == materialized_j) return;
            materialized_j = item.j;
            const auto scored =
                pending.template finalize_column<Semiring>(
                    item.j, chart, tails);
            (void)left.materialize_scored(
                scored,
                [](const cparty::scfg2::beam::VpCase3LeftKey &) {
                  return true;
                });
            selection.commit_column(item.j, left.records());
          },
          [&](const ItemKey &item, const auto &chart) {
            pending.template observe_completed_child<Semiring>(
                item, chart, child_event_index, factorized_ctx);
          },
          [&](const Deduction &deduction, const auto &chart, double local) {
            assert(selection.allows(deduction));
            const auto reconstructed =
                left.template reconstruct_source<Semiring>(deduction, local);
            if (reconstructed.has_value()) return *reconstructed;
            assert(cparty::scfg2::beam::vp_case3_partition(deduction) !=
                   cparty::scfg2::beam::VpCase3Partition::LeftNonempty);
            return cparty::scfg2::replay::
                runtime_deduction_contribution<Semiring>(
                    chart, deduction, local);
          });

  assert(baseline_schedule.size() == schedule.size());
  for (const auto &item : baseline_schedule) {
    assert(factorized.contains(item));
    assert(close(baseline.get(item), factorized.get(item)));
  }
  assert(left.candidate_count() == left_plan.geometry_count());
  assert(child_index.attachment_count() == left_plan.geometry_count());
  assert(pending.observed_attachment_count() ==
         child_index.attachment_count());
  assert(pending.registered_segment_count() ==
         child_index.future_segment_count());
  assert(pending.queried_source_count() == tail_plan.source_count());
  assert(pending.border_mismatch_count() == 0);
  assert(pending.duplicate_source_count() == 0);
  assert(left.materialization_count() == left.candidate_count());
  assert(left.attachment_count() == tail_plan.source_count());
  assert(selection.selected_count() == left.candidate_count());
  assert(selection.column_count() == left_plan.column_count());
  assert(selection.columns_are_ranked());
  assert(tails.tail_count() == tail_plan.tail_count());
  const auto child_profile = child_index.profile();
  assert(child_profile.signature_segment_count ==
         child_index.future_segment_count());
  return Result{tails.tail_count(),
                left.candidate_count(),
                left.attachment_count(),
                segments.size(),
                child_index.child_count(),
                child_index.max_attachments_per_child(),
                child_profile};
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
  check_column_parent_start_index();
  const auto inside =
      check_left_factorization<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi =
      check_left_factorization<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  assert(inside.tails == viterbi.tails);
  assert(inside.left_candidates == viterbi.left_candidates);
  assert(inside.sources == viterbi.sources);
  assert(inside.segments == viterbi.segments);
  assert(inside.wi_children == viterbi.wi_children);
  assert(inside.max_child_fanout == viterbi.max_child_fanout);
  assert(inside.child_profile.child_target_group_count ==
         viterbi.child_profile.child_target_group_count);
  assert(inside.child_profile.target_run_count ==
         viterbi.child_profile.target_run_count);
  assert(inside.child_profile.signature_segment_count ==
         viterbi.child_profile.signature_segment_count);
  std::cout << "n\ttails\tleft_candidates\tsegments\tsources\twi_children"
               "\tchild_target_groups\ttarget_runs\tsignature_segments"
               "\tmax_child_fanout\tmax_same_target\tmax_target_runs"
               "\tmax_signature_segments\tparity\n"
            << sequence.size() << '\t' << inside.tails << '\t'
            << inside.left_candidates << '\t' << inside.segments << '\t'
            << inside.sources << '\t' << inside.wi_children << '\t'
            << inside.child_profile.child_target_group_count << '\t'
            << inside.child_profile.target_run_count << '\t'
            << inside.child_profile.signature_segment_count << '\t'
            << inside.max_child_fanout << '\t'
            << inside.child_profile.maximum_same_target_count << '\t'
            << inside.child_profile.maximum_target_runs_per_child << '\t'
            << inside.child_profile.maximum_signature_segments_per_child
            << "\t1\n";
  return EXIT_SUCCESS;
}

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/source_selection.hh"
#include "scfg2/beam/virtual_candidate.hh"
#include "scfg2/beam/virtual_profile.hh"
#include "scfg2/beam/vp_case3_factor.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using cparty::scfg2::beam::VirtualBeamFamily;
using cparty::scfg2::beam::VpCase3Partition;
using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;
using cparty::scfg2::grammar::NonTerminal;
using cparty::scfg2::grammar::RuleId;
using cparty::scfg2::grammar::RuleSplit;

Deduction case3(cand_pos_t i,
                cand_pos_t j,
                cand_pos_t p,
                cand_pos_t q,
                cand_pos_t k,
                cand_pos_t l) {
  Deduction out{};
  out.parent = ItemKey{NonTerminal::VP, i, j, -1, -1};
  out.rule = RuleId::VP_WI_CASE3;
  out.split = RuleSplit{k, l, p, q};
  return out;
}

void check_synthetic_partition() {
  const Deduction boundary = case3(4, 30, 5, 10, 11, 24);
  const Deduction left = case3(4, 30, 8, 10, 14, 24);
  const Deduction other = case3(4, 30, 5, 10, 14, 24);

  assert(cparty::scfg2::beam::vp_case3_partition(boundary) ==
         VpCase3Partition::Boundary);
  assert(cparty::scfg2::beam::vp_case3_partition(left) ==
         VpCase3Partition::LeftNonempty);
  assert(cparty::scfg2::beam::vp_case3_partition(other) ==
         VpCase3Partition::Other);

  const auto boundary_candidates =
      cparty::scfg2::beam::virtual_beam_candidates_for_deduction(boundary);
  assert(boundary_candidates.size == 1);
  assert(boundary_candidates.candidates[0].full_identity.family ==
         VirtualBeamFamily::VpCase3Boundary);
  assert(boundary_candidates.candidates[0].parent_facing);

  const auto left_candidates =
      cparty::scfg2::beam::virtual_beam_candidates_for_deduction(left);
  assert(left_candidates.size == 2);
  assert(left_candidates.candidates[0].full_identity.family ==
         VirtualBeamFamily::VpCase3LeftTail);
  assert(!left_candidates.candidates[0].parent_facing);
  assert(left_candidates.candidates[1].full_identity.family ==
         VirtualBeamFamily::VpCase3Left);
  assert(left_candidates.candidates[1].parent_facing);
  assert(left_candidates.candidates[1].bucket.right_boundary == left.parent.j);
  assert(left_candidates.candidates[1].full_identity.source_split.k ==
         left.split.k);
  assert(left_candidates.candidates[1].full_identity.source_split.l ==
         left.split.l);
  assert(left_candidates.candidates[1].full_identity.source_split.p ==
         left.split.p);
  assert(left_candidates.candidates[1].full_identity.source_split.q ==
         left.split.q);

  const auto other_candidates =
      cparty::scfg2::beam::virtual_beam_candidates_for_deduction(other);
  assert(other_candidates.size == 1);
  assert(other_candidates.candidates[0].full_identity.family ==
         VirtualBeamFamily::VpCase3Other);

  Deduction not_case3 = left;
  not_case3.rule = RuleId::VP_STACK;
  assert(!cparty::scfg2::beam::vp_case3_partition(not_case3).has_value());
  assert(cparty::scfg2::beam::virtual_beam_candidates_for_deduction(not_case3)
             .size == 0);

  Deduction be{};
  be.parent = ItemKey{NonTerminal::BE, 3, 40, 12, 25};
  be.rule = RuleId::BE_WIP_WIP;
  be.split = RuleSplit{8, 34, -1, -1};
  const auto be_candidates =
      cparty::scfg2::beam::virtual_beam_candidates_for_deduction(be);
  assert(be_candidates.size == 3);
  assert(be_candidates.candidates[0].full_identity.family ==
         VirtualBeamFamily::BeCore);
  assert(be_candidates.candidates[0].bucket.right_boundary == be.split.l);
  assert(!be_candidates.candidates[0].parent_facing);
  assert(be_candidates.candidates[1].full_identity.family ==
         VirtualBeamFamily::BeRightTail);
  assert(be_candidates.candidates[1].bucket.right_boundary == be.parent.j);
  assert(!be_candidates.candidates[1].parent_facing);
  assert(be_candidates.candidates[2].full_identity.family ==
         VirtualBeamFamily::BeLeftRightTail);
  assert(be_candidates.candidates[2].parent_facing);
  assert(be_candidates.candidates[2].full_identity.source_parent.ip == 12);
  assert(be_candidates.candidates[2].full_identity.source_parent.jp == 25);

  cparty::scfg2::beam::VirtualBeamProfile profile;
  const Deduction left_tie = case3(4, 30, 9, 10, 14, 24);
  profile.observe(left, 2.0, false);
  profile.observe(left_tie, 2.0, true);
  profile.observe(case3(4, 30, 10, 10, 14, 24), 0.0, true);
  const auto *left_summary = profile.find({VirtualBeamFamily::VpCase3Left, 30});
  assert(left_summary);
  assert(left_summary->candidate_count == 2);
  assert(left_summary->selected_source_count == 1);
  assert(left_summary->selected_source_count_in_top(1) == 0);
  assert(left_summary->selected_source_count_in_top(2) == 1);
  assert(left_summary->proposed_unique_identity_count() == 2);
  const auto *tail_summary =
      profile.find({VirtualBeamFamily::VpCase3LeftTail, 30});
  assert(tail_summary);
  assert(tail_summary->candidate_count == 2);
  assert(tail_summary->proposed_unique_identity_count() == 1);
  assert(tail_summary->max_sources_per_proposed_identity() == 2);
  assert(tail_summary->selected_source_count_in_top_unique(1) == 1);
  assert(tail_summary->top_unique_source_contribution_sum(1) == 4.0);

  cparty::scfg2::engine::Chart<double> factor_chart(0.0);
  factor_chart.set(ItemKey{NonTerminal::WI, 5, 7, -1, -1}, 2.0);
  factor_chart.set(ItemKey{NonTerminal::WI, 11, 13, -1, -1}, 3.0);
  factor_chart.set(ItemKey{NonTerminal::WI, 25, 29, -1, -1}, 5.0);
  const auto factors =
      cparty::scfg2::beam::vp_case3_left_factor_values(left, factor_chart);
  assert(factors.available);
  assert(factors.left_child == 2.0);
  assert(factors.private_tail == 15.0);
  assert(factors.parent_facing_left == 30.0);
  assert(!cparty::scfg2::beam::vp_case3_left_factor_values(
              boundary, factor_chart)
              .available);
}

void check_exact_witness_classification() {
  std::string sequence =
      "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  std::string scaffold =
      "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  std::string viterbi_sequence = sequence;
  std::string viterbi_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot viterbi_owner(
      viterbi_sequence, viterbi_scaffold);
  const auto viterbi =
      cparty::scfg2::replay::run_w_final_exact_viterbi_traceback(
          viterbi_owner);
  std::vector<Deduction> selected_deductions;
  cparty::scfg2::beam::collect_selected_source_deductions(
      viterbi.trace, &selected_deductions);
  assert(!selected_deductions.empty());

  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::build_pkfree_basic_exact_schedule(
      static_cast<cand_pos_t>(sequence.size()), ctx);

  std::size_t direct_count = 0;
  std::size_t primary_count = 0;
  std::size_t left_tail_count = 0;
  std::size_t be_wip_wip_count = 0;
  for (const auto &item : schedule) {
    const auto deductions =
        cparty::scfg2::deductions::generate_exact_basic_deductions_for_item(
            item, ctx);
    for (const auto &deduction : deductions) {
      if (deduction.rule == RuleId::BE_WIP_WIP) {
        ++be_wip_wip_count;
        const auto candidates =
            cparty::scfg2::beam::virtual_beam_candidates_for_deduction(
                deduction);
        assert(candidates.size == 3);
        assert(candidates.candidates[0].bucket.right_boundary ==
               deduction.split.l);
        assert(candidates.candidates[1].bucket.right_boundary ==
               deduction.parent.j);
        assert(candidates.candidates[2].bucket.right_boundary ==
               deduction.parent.j);
      }
      if (deduction.rule != RuleId::VP_WI_CASE3) continue;
      ++direct_count;
      const auto partition = cparty::scfg2::beam::vp_case3_partition(deduction);
      assert(partition.has_value());
      const auto candidates =
          cparty::scfg2::beam::virtual_beam_candidates_for_deduction(deduction);
      std::size_t local_primary_count = 0;
      for (std::size_t idx = 0; idx < candidates.size; ++idx) {
        const auto &candidate = candidates.candidates[idx];
        assert(candidate.full_identity.source_parent.nonterminal ==
               deduction.parent.nonterminal);
        assert(candidate.bucket.right_boundary == deduction.parent.j);
        if (candidate.parent_facing) {
          ++local_primary_count;
          ++primary_count;
        } else {
          assert(candidate.full_identity.family ==
                 VirtualBeamFamily::VpCase3LeftTail);
          ++left_tail_count;
        }
      }
      assert(local_primary_count == 1);
      assert(candidates.size ==
             (*partition == VpCase3Partition::LeftNonempty ? 2u : 1u));
    }
  }

  assert(direct_count > 0);
  assert(be_wip_wip_count > 0);
  assert(primary_count == direct_count);

  const auto baseline =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<
          cparty::scfg2::semirings::SumProductSemiring>(
          schedule,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &,
             const Deduction &,
             double,
             const cparty::scfg2::engine::Chart<double> &,
             double) {});
  cparty::scfg2::beam::VirtualBeamProfile profile;
  const auto observed =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<
          cparty::scfg2::semirings::SumProductSemiring>(
          schedule,
          ctx,
          [](const Deduction &) { return true; },
          [&](const ItemKey &,
              const Deduction &deduction,
              double contribution,
              const cparty::scfg2::engine::Chart<double> &chart,
              double local) {
            const bool selected_source =
                cparty::scfg2::beam::source_deduction_is_selected(
                    deduction, selected_deductions);
            const auto factors =
                cparty::scfg2::beam::vp_case3_left_factor_values(
                    deduction, chart);
            cparty::scfg2::beam::VirtualBeamPriorityOverrides priorities;
            if (factors.available) {
              priorities.vp_case3_left_tail = factors.private_tail;
              priorities.vp_case3_left = factors.parent_facing_left;
              assert(std::fabs(
                         contribution -
                         local * factors.parent_facing_left) <=
                     std::max(1e-15, std::fabs(contribution) * 1e-12));
            }
            profile.observe(
                deduction, contribution, selected_source, priorities);
          });
  for (const auto &item : schedule) {
    assert(baseline.get(item) == observed.get(item));
  }
  assert(profile.bucket_count() > 0);
  assert(profile.total_candidate_observations() > 0);
  std::size_t selected_virtual_observations = 0;
  for (cand_pos_t j = 1; j <= static_cast<cand_pos_t>(sequence.size()); ++j) {
    for (const auto family : {
             VirtualBeamFamily::VpCase3Boundary,
             VirtualBeamFamily::VpCase3LeftTail,
             VirtualBeamFamily::VpCase3Left,
             VirtualBeamFamily::VpCase3Other,
             VirtualBeamFamily::BeCore,
             VirtualBeamFamily::BeRightTail,
             VirtualBeamFamily::BeLeftRightTail}) {
      const auto *summary = profile.find({family, j});
      if (!summary) continue;
      assert(summary->candidate_count ==
             summary->source_observations.size());
      assert(summary->full_source_identities_are_unique());
      assert(summary->proposed_identities_have_consistent_ranking_values());
      assert(summary->source_contribution_sum >= 0.0);
      assert(std::isfinite(summary->source_contribution_sum));
      assert(summary->top_source_contribution_sum(10) <=
             summary->source_contribution_sum + 1e-12);
      assert(summary->selected_source_count_in_top(10) <=
             summary->selected_source_count);
      if (family == VirtualBeamFamily::VpCase3LeftTail ||
          family == VirtualBeamFamily::VpCase3Left) {
        assert(summary->uses_factor_priority);
      }
      selected_virtual_observations += summary->selected_source_count;
      if (family == VirtualBeamFamily::VpCase3LeftTail ||
          family == VirtualBeamFamily::BeCore ||
          family == VirtualBeamFamily::BeRightTail) {
        assert(summary->parent_facing_count == 0);
      } else {
        assert(summary->parent_facing_count == summary->candidate_count);
      }
    }
  }
  assert(selected_virtual_observations > 0);
  std::cout << "scfg2_virtual_beam_candidate=ok direct_case3="
            << direct_count << " left_tail=" << left_tail_count
            << " be_wip_wip=" << be_wip_wip_count
            << " profiled=" << profile.total_candidate_observations()
            << " selected_virtual=" << selected_virtual_observations << '\n';
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  check_synthetic_partition();
  check_exact_witness_classification();
  return EXIT_SUCCESS;
}

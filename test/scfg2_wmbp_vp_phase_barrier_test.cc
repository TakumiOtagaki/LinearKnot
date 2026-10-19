#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/wmbp_vp_phase_barrier_profile.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

void check_wmbp_phase_order(const std::vector<ItemKey> &schedule) {
  cand_pos_t last_j = -1;
  int phase = 0;
  for (const auto &item : schedule) {
    if (item.j != last_j) {
      last_j = item.j;
      phase = 0;
    }
    const auto nt = item.nonterminal;
    const int item_phase = nt == cparty::scfg2::grammar::NonTerminal::VP
        ? 0
        : (nt == cparty::scfg2::grammar::NonTerminal::BE ||
                   nt == cparty::scfg2::grammar::NonTerminal::WMBP
               ? 1
               : 2);
    assert(item_phase >= phase);
    phase = item_phase;
  }
}

template <class Semiring>
void check_schedule_parity(const std::string &sequence,
                           const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto exact = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto barrier = cparty::scfg2::beam::
      build_wmbp_vp_phase_barrier_schedule(exact, ctx);
  const auto prerequisite = cparty::scfg2::beam::
      build_vp_prerequisite_barrier_schedule(exact, ctx);
  const auto static_barrier = cparty::scfg2::engine::
      build_pkfree_basic_exact_vp_barrier_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto wmbp_phase = cparty::scfg2::beam::
      build_nonterminal_prerequisite_barrier_schedule(
          static_barrier,
          ctx,
          cparty::scfg2::grammar::NonTerminal::WMBP);
  const auto be_wmbp_target_mask =
      (std::uint64_t{1} << static_cast<unsigned>(
           cparty::scfg2::grammar::NonTerminal::BE)) |
      (std::uint64_t{1} << static_cast<unsigned>(
           cparty::scfg2::grammar::NonTerminal::WMBP));
  const auto be_wmbp_phase = cparty::scfg2::beam::
      build_nonterminal_set_prerequisite_barrier_schedule(
          static_barrier, ctx, be_wmbp_target_mask);
  const auto static_wmbp_phase = cparty::scfg2::engine::
      build_pkfree_basic_exact_wmbp_barrier_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  assert(barrier.size() == exact.size());
  assert(prerequisite.feasible(exact.size()));
  assert(static_barrier.size() == exact.size());
  assert(static_barrier.size() == prerequisite.schedule.size());
  assert(wmbp_phase.feasible(exact.size()));
  assert(wmbp_phase.target_item_count > 0);
  assert(be_wmbp_phase.feasible(exact.size()));
  assert(static_wmbp_phase.size() == exact.size());
  check_wmbp_phase_order(static_wmbp_phase);
  for (std::size_t idx = 0; idx < static_barrier.size(); ++idx) {
    assert(cparty::scfg2::deductions::ItemKeyEq{}(
        static_barrier[idx], prerequisite.schedule[idx]));
  }
  const auto exact_chart =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          exact,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          });
  const auto barrier_chart =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          barrier,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          });
  const auto prerequisite_chart =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          prerequisite.schedule,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          });
  const auto static_chart =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          static_barrier,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          });
  const auto wmbp_phase_chart =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          wmbp_phase.schedule,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          });
  const auto be_wmbp_phase_chart =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          be_wmbp_phase.schedule,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          });
  const auto static_wmbp_phase_chart =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          static_wmbp_phase,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          });
  for (const auto &item : exact) {
    assert(close(exact_chart.get(item), barrier_chart.get(item)));
    assert(close(exact_chart.get(item), prerequisite_chart.get(item)));
    assert(close(exact_chart.get(item), static_chart.get(item)));
    assert(close(exact_chart.get(item), wmbp_phase_chart.get(item)));
    assert(close(exact_chart.get(item), be_wmbp_phase_chart.get(item)));
    assert(close(exact_chart.get(item), static_wmbp_phase_chart.get(item)));
  }
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
  assert(sequence.size() == scaffold.size());
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto profile = cparty::scfg2::beam::
      profile_wmbp_vp_phase_barrier(schedule, ctx);
  assert(profile.column_count == sequence.size());
  assert(profile.item_count == schedule.size());
  assert(profile.same_column_edge_count > 0);
  assert(profile.vp_barrier_edge_count > 0);
  assert(profile.barrier_wmbp_edge_count > 0);
  assert(profile.feasible());
  const auto prerequisite = cparty::scfg2::beam::
      build_vp_prerequisite_barrier_schedule(schedule, ctx);
  const auto static_barrier = cparty::scfg2::engine::
      build_pkfree_basic_exact_vp_barrier_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto wmbp_phase = cparty::scfg2::beam::
      build_nonterminal_prerequisite_barrier_schedule(
          static_barrier,
          ctx,
          cparty::scfg2::grammar::NonTerminal::WMBP);
  const auto be_wmbp_target_mask =
      (std::uint64_t{1} << static_cast<unsigned>(
           cparty::scfg2::grammar::NonTerminal::BE)) |
      (std::uint64_t{1} << static_cast<unsigned>(
           cparty::scfg2::grammar::NonTerminal::WMBP));
  const auto be_wmbp_phase = cparty::scfg2::beam::
      build_nonterminal_set_prerequisite_barrier_schedule(
          static_barrier, ctx, be_wmbp_target_mask);
  assert(prerequisite.feasible(schedule.size()));
  assert(prerequisite.prerequisite_nonterminal_mask ==
         (std::uint64_t{1} << static_cast<unsigned>(
              cparty::scfg2::grammar::NonTerminal::VP)));
  assert(wmbp_phase.feasible(schedule.size()));
  assert(wmbp_phase.target_item_count > 0);
  assert(be_wmbp_phase.feasible(schedule.size()));
  const auto wmbp_phase_family_mask =
      (std::uint64_t{1} << static_cast<unsigned>(
           cparty::scfg2::grammar::NonTerminal::VP)) |
      be_wmbp_target_mask;
  assert((be_wmbp_phase.prerequisite_nonterminal_mask &
          ~wmbp_phase_family_mask) == 0);
  check_schedule_parity<cparty::scfg2::semirings::SumProductSemiring>(
      sequence, scaffold);
  check_schedule_parity<cparty::scfg2::semirings::MaxProductSemiring>(
      sequence, scaffold);
  std::cout << "n\tcolumns\titems\tsame_column_edges\tvp_barrier_edges"
               "\tbarrier_wmbp_edges\tcyclic_columns\tmissing_children"
               "\tbefore_mask\tafter_mask\tsplit_mask"
               "\tvp_prerequisites\twmbp_prerequisites"
               "\tprerequisite_mask"
               "\twmbp_phase_prerequisites\twmbp_phase_targets"
               "\twmbp_phase_vp\twmbp_phase_be\twmbp_phase_wmbp"
               "\twmbp_phase_missing\twmbp_phase_mask"
               "\tbe_wmbp_phase_prerequisites\tbe_wmbp_phase_targets"
               "\tbe_wmbp_phase_vp\tbe_wmbp_phase_be"
               "\tbe_wmbp_phase_wmbp\tbe_wmbp_phase_missing"
               "\tbe_wmbp_phase_mask\tfeasible\n"
            << sequence.size() << '\t' << profile.column_count << '\t'
            << profile.item_count << '\t'
            << profile.same_column_edge_count << '\t'
            << profile.vp_barrier_edge_count << '\t'
            << profile.barrier_wmbp_edge_count << '\t'
            << profile.cyclic_column_count << '\t'
            << profile.missing_child_count << '\t'
            << profile.nonterminal_before_barrier_mask << '\t'
            << profile.nonterminal_after_barrier_mask << '\t'
            << profile.split_nonterminal_mask() << '\t'
            << prerequisite.prerequisite_item_count << '\t'
            << prerequisite.wmbp_prerequisite_count << '\t'
            << prerequisite.prerequisite_nonterminal_mask << '\t'
            << wmbp_phase.prerequisite_item_count << '\t'
            << wmbp_phase.target_item_count << '\t'
            << wmbp_phase.prerequisite_count(
                   cparty::scfg2::grammar::NonTerminal::VP) << '\t'
            << wmbp_phase.prerequisite_count(
                   cparty::scfg2::grammar::NonTerminal::BE) << '\t'
            << wmbp_phase.prerequisite_count(
                   cparty::scfg2::grammar::NonTerminal::WMBP) << '\t'
            << wmbp_phase.missing_child_count << '\t'
            << wmbp_phase.prerequisite_nonterminal_mask << '\t'
            << be_wmbp_phase.prerequisite_item_count << '\t'
            << be_wmbp_phase.target_item_count << '\t'
            << be_wmbp_phase.prerequisite_count(
                   cparty::scfg2::grammar::NonTerminal::VP) << '\t'
            << be_wmbp_phase.prerequisite_count(
                   cparty::scfg2::grammar::NonTerminal::BE) << '\t'
            << be_wmbp_phase.prerequisite_count(
                   cparty::scfg2::grammar::NonTerminal::WMBP) << '\t'
            << be_wmbp_phase.missing_child_count << '\t'
            << be_wmbp_phase.prerequisite_nonterminal_mask << '\t'
            << static_cast<int>(profile.feasible()) << '\n';
}

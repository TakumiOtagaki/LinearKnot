#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_item_selection.hh"
#include "scfg2/beam/wmbp_be_vp_prefix_events.hh"
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

struct BeamRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t candidates = 0;
  std::size_t selected = 0;
  std::size_t pruned = 0;
  std::size_t rejected_sources = 0;
  std::size_t generated_sources = 0;
  std::size_t evaluated_sources = 0;
};

template <class Semiring>
BeamRun run_width(const std::string &sequence,
                  const std::string &scaffold,
                  std::size_t width) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_vp_barrier_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  cparty::scfg2::beam::TopWidthVpItemSelection selection(width);
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
          [&](const ItemKey &item, const auto &) {
            if (item.nonterminal == cparty::scfg2::grammar::NonTerminal::VP ||
                item.j == committed_j) {
              return;
            }
            selection.commit_column(item.j);
            committed_j = item.j;
          },
          [&](const ItemKey &item, const auto &current_chart) {
            selection.observe_completed(item, current_chart.get(item));
          },
          [&](const Deduction &deduction,
              const auto &current_chart,
              double local) {
            if (!selection.allows_wmbp_source(deduction)) {
              ++rejected_sources;
              return Semiring::zero();
            }
            return cparty::scfg2::replay::
                runtime_deduction_contribution<Semiring>(
                    current_chart, deduction, local);
          });
  assert(selection.respects_width());
  return BeamRun{std::move(chart),
                 selection.candidate_count(),
                 selection.selected_count(),
                 selection.pruned_count(),
                 rejected_sources,
                 0,
                 0};
}

template <class Semiring>
BeamRun run_driven_width(const std::string &sequence,
                         const std::string &scaffold,
                         std::size_t width) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_vp_barrier_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  cparty::scfg2::beam::TopWidthVpItemSelection selection(width);
  cand_pos_t committed_j = -1;
  std::optional<ItemKey> pending_parent;
  std::vector<Deduction> pending_deductions;
  std::size_t generated_sources = 0;
  std::size_t evaluated_sources = 0;
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
          [&](const ItemKey &,
              const Deduction &deduction,
              double,
              const auto &,
              double) {
            evaluated_sources +=
                cparty::scfg2::beam::is_wmbp_be_vp_rule(deduction.rule)
                ? 1
                : 0;
          },
          [&](const ItemKey &item, const auto &) {
            if (item.nonterminal != cparty::scfg2::grammar::NonTerminal::VP &&
                item.j != committed_j) {
              selection.commit_column(item.j);
              committed_j = item.j;
            }
            if (item.nonterminal ==
                cparty::scfg2::grammar::NonTerminal::WMBP) {
              pending_parent = item;
              pending_deductions = cparty::scfg2::beam::
                  generate_wmbp_parent_deductions_for_retained_vp_starts(
                      item, selection.selected_starts(item.j), ctx);
              for (const auto &deduction : pending_deductions) {
                generated_sources +=
                    cparty::scfg2::beam::is_wmbp_be_vp_rule(deduction.rule)
                    ? 1
                    : 0;
              }
            }
          },
          [&](const ItemKey &item, const auto &current_chart) {
            selection.observe_completed(item, current_chart.get(item));
          },
          [&](const Deduction &deduction,
              const auto &current_chart,
              double local) {
            return cparty::scfg2::replay::
                runtime_deduction_contribution<Semiring>(
                    current_chart, deduction, local);
          });
  assert(!pending_parent.has_value());
  assert(selection.respects_width());
  return BeamRun{std::move(chart),
                 selection.candidate_count(),
                 selection.selected_count(),
                 selection.pruned_count(),
                 0,
                 generated_sources,
                 evaluated_sources};
}

template <class Semiring>
auto run_direct(const std::string &sequence, const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_vp_barrier_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  return cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      [](const Deduction &) { return true; },
      [](const ItemKey &, const Deduction &, double, const auto &, double) {});
}

template <class Semiring>
BeamRun check_widths(const std::string &sequence,
                     const std::string &scaffold) {
  const auto keep_all = run_width<Semiring>(
      sequence, scaffold, std::numeric_limits<std::size_t>::max());
  const auto width_one = run_width<Semiring>(sequence, scaffold, 1);
  const auto driven_keep_all = run_driven_width<Semiring>(
      sequence, scaffold, std::numeric_limits<std::size_t>::max());
  const auto driven_width_one =
      run_driven_width<Semiring>(sequence, scaffold, 1);
  const auto direct = run_direct<Semiring>(sequence, scaffold);
  assert(keep_all.candidates == keep_all.selected);
  assert(keep_all.pruned == 0);
  for (const auto &[item, exact] : direct.values()) {
    assert(close(keep_all.chart.get(item), exact));
    assert(close(driven_keep_all.chart.get(item), exact));
    const double approximate = width_one.chart.get(item);
    assert(close(driven_width_one.chart.get(item), approximate));
    assert(approximate <= exact || close(approximate, exact));
  }
  assert(width_one.candidates == driven_width_one.candidates);
  assert(width_one.selected == driven_width_one.selected);
  auto out = width_one;
  out.generated_sources = driven_width_one.generated_sources;
  out.evaluated_sources = driven_width_one.evaluated_sources;
  return out;
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
  assert(inside.candidates == viterbi.candidates);
  if (argc == 1) {
    assert(inside.pruned > 0);
    assert(inside.rejected_sources > 0);
  }
  std::cout << "n\tvp_candidates\tvp_selected\tvp_pruned"
               "\trejected_wmbp_sources\tdriven_wmbp_sources"
               "\tdriven_admitted_sources\tmonotone\n"
            << sequence.size() << '\t' << inside.candidates << '\t'
            << inside.selected << '\t' << inside.pruned << '\t'
            << inside.rejected_sources << '\t'
            << inside.generated_sources << '\t'
            << inside.evaluated_sources << "\t1\n";
}

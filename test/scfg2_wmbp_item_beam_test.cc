#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/wmbp_item_selection.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
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

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
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

bool is_wmbp_phase_item(const ItemKey &item) {
  const auto nt = item.nonterminal;
  return nt == cparty::scfg2::grammar::NonTerminal::VP ||
         nt == cparty::scfg2::grammar::NonTerminal::BE ||
         nt == cparty::scfg2::grammar::NonTerminal::WMBP;
}

struct BeamRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t candidates = 0;
  std::size_t selected = 0;
  std::size_t pruned = 0;
  std::size_t rejected_consumer_sources = 0;
  std::size_t driven_consumer_sources = 0;
};

bool is_wmbp_consumer_source(const Deduction &deduction) {
  if (deduction.parent.nonterminal ==
      cparty::scfg2::grammar::NonTerminal::WMBP) {
    return false;
  }
  for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
    if (deduction.children[idx].item.nonterminal ==
        cparty::scfg2::grammar::NonTerminal::WMBP) {
      return true;
    }
  }
  return false;
}

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
      build_pkfree_basic_exact_wmbp_barrier_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  cparty::scfg2::beam::TopWidthWmbpItemSelection selection(width);
  cand_pos_t committed_j = -1;
  std::size_t rejected_consumer_sources = 0;
  auto chart = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_item_observers_and_contribution<
          Semiring>(
          schedule,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          },
          [&](const ItemKey &item, const auto &) {
            if (is_wmbp_phase_item(item) || item.j == committed_j) return;
            selection.commit_column(item.j);
            committed_j = item.j;
          },
          [&](const ItemKey &item, const auto &current_chart) {
            selection.observe_completed(item, current_chart.get(item));
          },
          [&](const Deduction &deduction,
              const auto &current_chart,
              double local) {
            if (!selection.allows_consumer_source(deduction)) {
              ++rejected_consumer_sources;
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
                 rejected_consumer_sources,
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
      build_pkfree_basic_exact_wmbp_barrier_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  cparty::scfg2::beam::TopWidthWmbpItemSelection selection(width);
  cand_pos_t committed_j = -1;
  std::size_t driven_consumer_sources = 0;
  auto chart = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_provider_item_observers_and_contribution<
          Semiring>(
          schedule,
          ctx,
          [&](const ItemKey &item) {
            if (item.nonterminal ==
                    cparty::scfg2::grammar::NonTerminal::WMBW ||
                item.nonterminal ==
                    cparty::scfg2::grammar::NonTerminal::WMB) {
              auto deductions = cparty::scfg2::deductions::
                  generate_exact_wmbp_consumer_deductions_for_retained_ends(
                      item, ctx, selection.selected_ends(item.i));
              auto filtered_oracle = cparty::scfg2::deductions::
                  generate_exact_basic_deductions_for_item(item, ctx);
              filtered_oracle.erase(
                  std::remove_if(
                      filtered_oracle.begin(),
                      filtered_oracle.end(),
                      [&](const Deduction &deduction) {
                        return !selection.allows_consumer_source(deduction);
                      }),
                  filtered_oracle.end());
              assert(deductions.size() == filtered_oracle.size());
              for (std::size_t idx = 0; idx < deductions.size(); ++idx) {
                assert(same_deduction(deductions[idx], filtered_oracle[idx]));
              }
              for (const auto &deduction : deductions) {
                driven_consumer_sources +=
                    is_wmbp_consumer_source(deduction) ? 1 : 0;
              }
              return deductions;
            }
            return cparty::scfg2::deductions::
                generate_exact_basic_deductions_for_item(item, ctx);
          },
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          },
          [&](const ItemKey &item, const auto &) {
            if (is_wmbp_phase_item(item) || item.j == committed_j) return;
            selection.commit_column(item.j);
            committed_j = item.j;
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
  assert(selection.respects_width());
  return BeamRun{std::move(chart),
                 selection.candidate_count(),
                 selection.selected_count(),
                 selection.pruned_count(),
                 0,
                 driven_consumer_sources};
}

template <class Semiring>
auto run_direct(const std::string &sequence, const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_wmbp_barrier_schedule(
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
  out.driven_consumer_sources = driven_width_one.driven_consumer_sources;
  return out;
}

}  // namespace

int main(int argc, char **argv) {
  {
    // Both increasing runtime commits and the older helper's out-of-order
    // commits must expose the same sorted endpoint index to consumers.
    using namespace cparty::scfg2;
    auto selection = beam::TopWidthWmbpItemSelection::keep_all();
    for (cand_pos_t j : {3, 5, 4, 2}) {
      selection.observe_completed({grammar::NonTerminal::WMBP, 1, j, -1, -1}, 1.);
      selection.commit_column(j);
    }
    assert((selection.selected_ends(1) == std::vector<cand_pos_t>{2, 3, 4, 5}));
  }
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
  assert(inside.selected == viterbi.selected);
  if (argc == 1) {
    assert(inside.pruned > 0);
    assert(inside.rejected_consumer_sources > 0);
  }
  std::cout << "n\twmbp_candidates\twmbp_selected\twmbp_pruned"
               "\trejected_consumer_sources\tdriven_consumer_sources"
               "\tmonotone\n"
            << sequence.size() << '\t' << inside.candidates << '\t'
            << inside.selected << '\t' << inside.pruned << '\t'
            << inside.rejected_consumer_sources << '\t'
            << inside.driven_consumer_sources << "\t1\n";
}

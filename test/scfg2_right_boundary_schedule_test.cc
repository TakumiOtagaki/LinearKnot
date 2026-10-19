#include "scfg2/weights/runtime_snapshot.hh"
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

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

struct ScheduleParityResult {
  std::size_t item_count = 0;
  std::size_t legacy_excluded_deduction_count = 0;
};

template <class Semiring>
ScheduleParityResult check_schedule_parity(const std::string &sequence,
                                           const std::string &scaffold) {
  std::string span_sequence = sequence;
  std::string span_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot span_owner(
      span_sequence, span_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext span_ctx(span_owner);
  const auto span_schedule =
      cparty::scfg2::engine::build_pkfree_basic_exact_schedule(
          static_cast<cand_pos_t>(sequence.size()), span_ctx);
  std::size_t legacy_excluded_deduction_count = 0;
  for (const auto &item : span_schedule) {
    const auto deductions = cparty::scfg2::deductions::
        generate_exact_basic_deductions_for_item(item, span_ctx);
    for (const auto &deduction : deductions) {
      if (!cparty::scfg2::engine::legacy_span_schedule_supports(deduction)) {
        ++legacy_excluded_deduction_count;
      }
    }
  }
  const auto span_chart =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          span_schedule,
          span_ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &,
             const Deduction &,
             double,
             const cparty::scfg2::engine::Chart<double> &,
             double) {});

  std::string right_sequence = sequence;
  std::string right_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot right_owner(
      right_sequence, right_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext right_ctx(right_owner);
  const auto right_schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), right_ctx);
  const auto right_chart =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          right_schedule,
          right_ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &,
             const Deduction &,
             double,
             const cparty::scfg2::engine::Chart<double> &,
             double) {});

  assert(span_schedule.size() == right_schedule.size());
  for (const auto &item : span_schedule) {
    assert(right_chart.contains(item));
    if (!close(span_chart.get(item), right_chart.get(item))) {
      std::cerr << "schedule mismatch nt="
                << cparty::scfg2::grammar::nonterminal_name(item.nonterminal)
                << " i=" << item.i << " j=" << item.j
                << " ip=" << item.ip << " jp=" << item.jp
                << " span=" << span_chart.get(item)
                << " right=" << right_chart.get(item) << '\n';
      assert(false);
    }
  }
  for (const auto &item : right_schedule) {
    assert(span_chart.contains(item));
  }
  return ScheduleParityResult{
      span_schedule.size(), legacy_excluded_deduction_count};
}

}  // namespace

int main(int argc, char **argv) {
  assert(argc == 1 || argc == 3 || argc == 4);
  if (argc == 4) {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[3], 1);
  } else {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  }
  const std::string sequence = argc >= 3
      ? argv[1]
      : "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  const std::string scaffold = argc >= 3
      ? argv[2]
      : "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  assert(sequence.size() == scaffold.size());
  const auto inside_items =
      check_schedule_parity<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi_items =
      check_schedule_parity<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  assert(inside_items.item_count == viterbi_items.item_count);
  assert(inside_items.legacy_excluded_deduction_count ==
         viterbi_items.legacy_excluded_deduction_count);
  std::cout << "n\titems\tlegacy_excluded\tinside_parity\tviterbi_parity\n"
            << sequence.size() << '\t' << inside_items.item_count << '\t'
            << inside_items.legacy_excluded_deduction_count << "\t1\t1\n";
  return EXIT_SUCCESS;
}

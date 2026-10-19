#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/be_right_tail_retained_core_index.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
      std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

template <class Semiring>
std::array<std::size_t, 4> check_index(cparty::scfg2::weights::RuntimeSnapshot &owner) {
  using namespace cparty::scfg2;
  replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const auto schedule =
      engine::build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  const auto chart = replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      [](const Deduction &) { return true; },
      [](const ItemKey &, const Deduction &,
         typename Semiring::value_type, const auto &, double) {});
  const auto plan = beam::ExactBeRightTailColumnPlan::build(schedule, ctx);
  auto dynamic = beam::BeRightTailRetainedCoreIndex::build(n, ctx);
  beam::ExactBeFrameCache planned_cache;
  beam::ExactBeFrameCache dynamic_cache;

  std::size_t positive_planned = 0;
  for (cand_pos_t j = 1; j <= n; ++j) {
    auto expected = plan.materialize_column(j, chart, planned_cache);
    expected.erase(
        std::remove_if(expected.begin(), expected.end(), [](const auto &tail) {
          return tail.factor == 0.0;
        }),
        expected.end());
    auto actual = dynamic.take_column(j);
    const auto rank = [](const auto &lhs, const auto &rhs) {
      return beam::be_right_tail_tie_key(lhs.key) <
          beam::be_right_tail_tie_key(rhs.key);
    };
    std::sort(expected.begin(), expected.end(), rank);
    std::sort(actual.begin(), actual.end(), rank);
    assert(expected.size() == actual.size());
    positive_planned += expected.size();
    for (std::size_t idx = 0; idx < expected.size(); ++idx) {
      assert(expected[idx].key == actual[idx].key);
      assert(close(expected[idx].factor, actual[idx].factor));
    }

    std::vector<ItemKey> retained_be;
    for (const auto &item : schedule) {
      if (item.nonterminal == grammar::NonTerminal::BE && item.j == j &&
          chart.get(item) != Semiring::zero()) {
        retained_be.push_back(item);
      }
    }
    dynamic.observe_retained_be_column(retained_be, ctx);
    for (const auto &item : schedule) {
      if (item.nonterminal == grammar::NonTerminal::WIP && item.j == j &&
          chart.get(item) != Semiring::zero()) {
        dynamic.observe_completed_child(item, chart, dynamic_cache);
      }
    }
  }
  assert(dynamic.pending_empty());
  assert(dynamic.stats().duplicate_registrations == 0);
  assert(dynamic.stats().materialized_tails == positive_planned);
  assert(dynamic.stats().committed_tails == positive_planned);
  return {positive_planned,
          dynamic.stats().retained_core_events,
          dynamic.stats().immediate_parent_probes,
          dynamic.stats().registered_tails};
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
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  const auto inside =
      check_index<cparty::scfg2::semirings::SumProductSemiring>(owner);
  const auto viterbi =
      check_index<cparty::scfg2::semirings::MaxProductSemiring>(owner);
  assert(inside == viterbi);
  std::cout << "positive_tails\tretained_cores\tparent_probes"
               "\tregistered_tails\n"
            << inside[0] << '\t' << inside[1] << '\t' << inside[2] << '\t'
            << inside[3] << '\n';
  return EXIT_SUCCESS;
}

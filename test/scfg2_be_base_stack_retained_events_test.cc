#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/be_base_stack_retained_events.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

bool close(double lhs, double rhs) {
  const double scale = std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  return std::fabs(lhs - rhs) <= 1e-10 * scale;
}

template <class Semiring>
std::array<std::size_t, 2> check_events(cparty::scfg2::weights::RuntimeSnapshot &owner) {
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
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);

  using Totals = std::unordered_map<ItemKey,
                                    typename Semiring::value_type,
                                    deductions::ItemKeyHash,
                                    deductions::ItemKeyEq>;
  using SourceId = std::tuple<grammar::RuleId,
                              cand_pos_t,
                              cand_pos_t,
                              cand_pos_t,
                              cand_pos_t>;
  Totals driven;
  std::set<SourceId> driven_ids;
  beam::BeBaseStackRetainedEvents events;
  for (cand_pos_t j = 1; j <= n; ++j) {
    events.prepare_column(j, ctx);
    for (const auto &parent : events.parents()) {
      for (const auto &source : events.deductions_for(parent)) {
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, source, weight.weight_for(source));
        assert(term != Semiring::zero());
        assert(driven_ids.emplace(source.rule,
                                  parent.i,
                                  parent.j,
                                  parent.ip,
                                  parent.jp)
                   .second);
        auto [found, inserted] = driven.emplace(parent, Semiring::zero());
        found->second = Semiring::plus(found->second, term);
      }
    }
    std::vector<ItemKey> retained;
    for (const auto &item : schedule) {
      if (item.nonterminal == grammar::NonTerminal::BE && item.j == j &&
          chart.get(item) != Semiring::zero()) {
        retained.push_back(item);
      }
    }
    events.observe_committed_column(retained, n, ctx);
  }
  assert(events.stats().duplicate_sources == 0);

  const std::array<grammar::RuleId, 2> rules{
      grammar::RuleId::BE_BASE_SAMEPAIR,
      grammar::RuleId::BE_STACK};
  std::set<SourceId> exact_ids;
  std::array<std::size_t, 2> counts{};
  for (const auto &parent : schedule) {
    if (parent.nonterminal != grammar::NonTerminal::BE) continue;
    auto aggregate = Semiring::zero();
    for (std::size_t idx = 0; idx < rules.size(); ++idx) {
      for (const auto &source : deductions::generate_exact_basic_deductions(
               parent, rules[idx], ctx)) {
        if (!engine::legacy_span_schedule_supports(source)) continue;
        const double local = weight.weight_for(source);
        if (local == 0.0) continue;
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, source, local);
        if (term == Semiring::zero()) continue;
        aggregate = Semiring::plus(aggregate, term);
        exact_ids.emplace(
            source.rule, parent.i, parent.j, parent.ip, parent.jp);
        ++counts[idx];
      }
    }
    const auto found = driven.find(parent);
    const auto actual = found == driven.end()
        ? Semiring::zero()
        : found->second;
    if (!close(aggregate, actual)) {
      std::cerr << "BE base/stack mismatch parent=" << parent.i << ','
                << parent.j << ',' << parent.ip << ',' << parent.jp
                << " exact=" << aggregate << " driven=" << actual << '\n';
      assert(false);
    }
  }
  assert(driven_ids == exact_ids);
  assert(events.stats().base_events == counts[0]);
  assert(events.stats().stack_sources == counts[1]);
  return counts;
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
      check_events<cparty::scfg2::semirings::SumProductSemiring>(owner);
  const auto viterbi =
      check_events<cparty::scfg2::semirings::MaxProductSemiring>(owner);
  assert(inside == viterbi);
  if (argc == 1) assert(inside[0] > 0 && inside[1] > 0);
  std::cout << "be_base\tbe_stack\n"
            << inside[0] << '\t' << inside[1] << '\n';
  return EXIT_SUCCESS;
}

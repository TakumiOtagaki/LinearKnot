#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

bool close(double lhs, double rhs) {
  const double scale = std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  return std::fabs(lhs - rhs) <= 1e-10 * scale;
}

template <class Semiring>
void check_factor(cparty::scfg2::weights::RuntimeSnapshot &owner) {
  using namespace cparty::scfg2;
  replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const auto schedule =
      engine::build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  const auto chart = replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule, ctx,
      [](const deductions::Deduction &) { return true; },
      [](const deductions::ItemKey &, const deductions::Deduction &,
         typename Semiring::value_type, const auto &, double) {});
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);
  std::size_t checked = 0;
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = 1; i < j; ++i) {
      const deductions::ItemKey parent{
          grammar::NonTerminal::WMB, i, j, -1, -1};
      auto exact = Semiring::zero();
      auto factored = Semiring::zero();
      std::unordered_map<cand_pos_t, typename Semiring::value_type> tails;
      std::unordered_map<cand_pos_t, deductions::Deduction> representatives;
      for (const auto &deduction : deductions::generate_exact_basic_deductions(
               parent, grammar::RuleId::WMB_SPLIT_BE_WMBP_WI, ctx)) {
        const auto direct_term =
            replay::runtime_deduction_contribution<Semiring>(
                chart, deduction, weight.weight_for(deduction));
        exact = Semiring::plus(exact, direct_term);
        const cand_pos_t t = deduction.split.l;
        auto tail = chart.get(deductions::ItemKey{
            grammar::NonTerminal::WMBP,
            i,
            deduction.split.k,
            -1,
            -1});
        if (deduction.split.k + 1 <= t - 1) {
          tail = Semiring::times(
              tail,
              chart.get(deductions::ItemKey{
                  grammar::NonTerminal::WI,
                  deduction.split.k + 1,
                  t - 1,
                  -1,
                  -1}));
        }
        tails[t] = Semiring::plus(tails[t], tail);
        representatives.emplace(t, deduction);
      }
      for (const auto &[t, tail] : tails) {
        const auto &deduction = representatives.at(t);
        const auto be = chart.get(deductions::ItemKey{
            grammar::NonTerminal::BE,
            deduction.split.p,
            j,
            deduction.split.q,
            t});
        const auto term = Semiring::times(
            Semiring::times(be, tail), weight.weight_for(deduction));
        factored = Semiring::plus(factored, term);
      }
      assert(close(factored, exact));
      if (exact != Semiring::zero()) ++checked;
    }
  }
  assert(checked > 0);
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  std::string sequence =
      "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  std::string scaffold =
      "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  check_factor<cparty::scfg2::semirings::SumProductSemiring>(owner);
  check_factor<cparty::scfg2::semirings::MaxProductSemiring>(owner);
  std::cout << "scfg2_wmb_btail_factor=ok\n";
  return EXIT_SUCCESS;
}

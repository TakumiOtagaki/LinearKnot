#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vpr_base_prefix.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

bool close(double lhs, double rhs) {
  const double scale = std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  return std::fabs(lhs - rhs) <= 1e-10 * scale;
}

template <class Semiring>
void check_prefix(cparty::scfg2::weights::RuntimeSnapshot &owner) {
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
  beam::VprBasePrefixColumns<Semiring> prefix;
  std::size_t checked = 0;
  for (cand_pos_t j = 1; j <= n; ++j) {
    prefix.build_column(j, chart, ctx);
    for (cand_pos_t i = 1; i < j; ++i) {
      const deductions::ItemKey parent{
          grammar::NonTerminal::VPR, i, j, -1, -1};
      auto exact = Semiring::zero();
      auto best = Semiring::zero();
      cand_pos_t best_k = -1;
      for (const auto &deduction : deductions::generate_exact_basic_deductions(
               parent, grammar::RuleId::VPR_SPLIT_VP_BASEPAIR, ctx)) {
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, deduction, weight.weight_for(deduction));
        exact = Semiring::plus(exact, term);
        if (Semiring::better(term, best)) {
          best = term;
          best_k = deduction.split.k;
        }
      }
      if (!deductions::exact_basic_item_parent_active(parent, ctx)) {
        exact = Semiring::zero();
      }
      const cand_pos_t max_bp =
          std::max(ctx.border_B(i, j), ctx.border_bp(i, j));
      const auto &cell = prefix.cell(i);
      const auto factored =
          deductions::exact_basic_item_parent_active(parent, ctx) &&
                  cell.best_k > max_bp
              ? cell.total
              : Semiring::zero();
      assert(close(factored, exact));
      if (exact != Semiring::zero()) {
        assert(cell.best_k == best_k);
        ++checked;
      }
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
  check_prefix<cparty::scfg2::semirings::SumProductSemiring>(owner);
  check_prefix<cparty::scfg2::semirings::MaxProductSemiring>(owner);
  std::cout << "scfg2_vpr_base_prefix=ok\n";
  return EXIT_SUCCESS;
}

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vpl_start_prefix.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
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
#include <vector>

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
  std::size_t checked = 0;
  for (cand_pos_t j = 1; j <= n; ++j) {
    beam::VplStartPrefixColumn<Semiring> prefix;
    prefix.build(j, chart, ctx);
    std::vector<cand_pos_t> vp_starts;
    std::vector<std::pair<double, cand_pos_t>> ranked;
    for (cand_pos_t k = 1; k < j; ++k) {
      if (chart.get(deductions::ItemKey{
              grammar::NonTerminal::VP, k, j, -1, -1}) != Semiring::zero()) {
        vp_starts.push_back(k);
      }
    }
    for (cand_pos_t i = 1; i < j; ++i) {
      const deductions::ItemKey parent{
          grammar::NonTerminal::VPL, i, j, -1, -1};
      auto exact = Semiring::zero();
      auto best = Semiring::zero();
      cand_pos_t best_k = -1;
      for (const auto &deduction : deductions::generate_exact_basic_deductions(
               parent, grammar::RuleId::VPL_SPLIT_VP, ctx)) {
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, deduction, weight.weight_for(deduction));
        exact = Semiring::plus(exact, term);
        if (Semiring::better(term, best)) {
          best = term;
          best_k = deduction.split.k;
        }
      }
      const bool active = deductions::exact_basic_item_parent_active(parent, ctx);
      if (!active) exact = Semiring::zero();
      const auto factored = active ? prefix.contribution(i) : Semiring::zero();
      assert(close(factored, exact));
      if (exact != Semiring::zero()) {
        ranked.emplace_back(exact, i);
        assert(prefix.cell(i).best_k == best_k);
        ++checked;
      }
    }
    std::sort(ranked.begin(), ranked.end(), [](const auto &lhs,
                                                const auto &rhs) {
      if (lhs.first != rhs.first) return lhs.first > rhs.first;
      return lhs.second < rhs.second;
    });
    // Very large finite widths and KeepAll must produce the full feasible
    // list without attempting to reserve one allocation slot per width.
    for (const std::size_t width : {std::size_t{1}, std::size_t{3},
                                    std::size_t{10},
                                    std::numeric_limits<std::size_t>::max() / 2,
                                    std::numeric_limits<std::size_t>::max()}) {
      const auto selected = beam::select_vpl_start_prefix_top_width<Semiring>(
          j, chart, ctx, vp_starts, width);
      assert(selected.size() == std::min(width, ranked.size()));
      for (std::size_t rank = 0; rank < selected.size(); ++rank) {
        assert(selected[rank].i == ranked[rank].second);
        assert(close(selected[rank].total, ranked[rank].first));
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
  std::cout << "scfg2_vpl_start_prefix=ok\n";
  return EXIT_SUCCESS;
}

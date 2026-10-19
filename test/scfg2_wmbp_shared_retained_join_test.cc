#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/wmbp_shared_retained_join.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

bool close(double lhs, double rhs) {
  const double scale = std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  return std::fabs(lhs - rhs) <= 1e-10 * scale;
}

template <class Semiring>
std::size_t check_join(cparty::scfg2::weights::RuntimeSnapshot &owner) {
  using namespace cparty::scfg2;
  replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const auto schedule =
      engine::build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  const auto chart = replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      [](const deductions::Deduction &) { return true; },
      [](const deductions::ItemKey &, const deductions::Deduction &,
         typename Semiring::value_type, const auto &, double) {});
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);

  std::vector<std::vector<cand_pos_t>> vp_starts(n + 1);
  std::vector<std::vector<cand_pos_t>> wmbp_starts(n + 1);
  std::vector<std::vector<cand_pos_t>> wmbw_starts(n + 1);
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = 1; i <= j; ++i) {
      if (chart.get(deductions::ItemKey{
              grammar::NonTerminal::VP, i, j, -1, -1}) != Semiring::zero()) {
        vp_starts[j].push_back(i);
      }
      if (chart.get(deductions::ItemKey{
              grammar::NonTerminal::WMBP, i, j, -1, -1}) !=
          Semiring::zero()) {
        wmbp_starts[j].push_back(i);
      }
      if (chart.get(deductions::ItemKey{
              grammar::NonTerminal::WMBW, i, j, -1, -1}) !=
          Semiring::zero()) {
        wmbw_starts[j].push_back(i);
      }
    }
  }

  std::size_t exact_positive_sources = 0;
  std::size_t joined_positive_sources = 0;
  std::size_t total_probes = 0;
  for (cand_pos_t j = 1; j <= n; ++j) {
    beam::WmbpSharedRetainedJoin<replay::WFinalExactInsideContext> join;
    join.build_column(
        j,
        vp_starts[j],
        [&](grammar::NonTerminal family,
            cand_pos_t right) -> const std::vector<cand_pos_t> & {
          assert(right >= 0 && right <= n);
          return family == grammar::NonTerminal::WMBP
              ? wmbp_starts[right]
              : wmbw_starts[right];
        },
        ctx);
    total_probes += join.stats().wmbp_left_probes;
    total_probes += join.stats().wmbw_left_probes;

    std::unordered_map<deductions::ItemKey,
                       typename Semiring::value_type,
                       deductions::ItemKeyHash,
                       deductions::ItemKeyEq>
        joined;
    for (const auto &parent : join.parents()) {
      auto &total = joined[parent];
      for (const auto &deduction : join.deductions_for(parent)) {
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, deduction, weight.weight_for(deduction));
        total = Semiring::plus(total, term);
        if (term != Semiring::zero()) ++joined_positive_sources;
      }
    }

    for (cand_pos_t i = 1; i <= j; ++i) {
      const deductions::ItemKey parent{
          grammar::NonTerminal::WMBP, i, j, -1, -1};
      auto exact = Semiring::zero();
      for (const auto &deduction :
           deductions::generate_exact_basic_deductions_for_item(parent, ctx)) {
        if (deduction.rule !=
                grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP &&
            deduction.rule !=
                grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP) {
          continue;
        }
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, deduction, weight.weight_for(deduction));
        exact = Semiring::plus(exact, term);
        if (term != Semiring::zero()) ++exact_positive_sources;
      }
      const auto found = joined.find(parent);
      const auto factored = found == joined.end() ? Semiring::zero()
                                                  : found->second;
      assert(close(exact, factored));
    }
  }
  assert(joined_positive_sources == exact_positive_sources);
  if (exact_positive_sources > 0) assert(total_probes > 0);
  return exact_positive_sources;
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
  const auto inside_sources =
      check_join<cparty::scfg2::semirings::SumProductSemiring>(owner);
  const auto viterbi_sources =
      check_join<cparty::scfg2::semirings::MaxProductSemiring>(owner);
  if (argc == 1) {
    assert(inside_sources > 0);
    assert(viterbi_sources > 0);
  }
  std::cout << "scfg2_wmbp_shared_retained_join=ok\n";
  return EXIT_SUCCESS;
}

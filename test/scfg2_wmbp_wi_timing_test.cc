#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/wmbp_be_vp_prefix_profile.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

struct TimingCounts {
  std::size_t before = 0;
  std::size_t same = 0;
  std::size_t after = 0;
  std::size_t empty = 0;
  std::size_t nonempty = 0;
};

TimingCounts observe(cparty::scfg2::weights::RuntimeSnapshot &owner) {
  using namespace cparty::scfg2;
  replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = engine::build_pkfree_basic_exact_right_boundary_schedule(
      static_cast<cand_pos_t>(owner.length()), ctx);
  TimingCounts out;
  (void)replay::run_w_final_exact_dp_schedule<
      semirings::SumProductSemiring>(
      schedule,
      ctx,
      [](const deductions::Deduction &) { return true; },
      [&](const deductions::ItemKey &,
          const deductions::Deduction &deduction,
          double candidate,
          const auto &,
          double) {
        if (deduction.rule != grammar::RuleId::WMBP_SPLIT_BE_WI_VP ||
            candidate == 0.0) {
          return;
        }
        const auto &be = deduction.children[0].item;
        const auto &vp = beam::wmbp_be_vp_right_child(deduction);
        if (be.j < vp.i) {
          ++out.before;
        } else if (be.j == vp.i) {
          ++out.same;
        } else {
          ++out.after;
        }
        if (deduction.child_count == 2) {
          ++out.empty;
        } else {
          ++out.nonempty;
        }
      });
  assert(out.before + out.same + out.after == out.empty + out.nonempty);
  return out;
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
  const auto counts = observe(owner);
  assert(counts.before == 0);
  assert(counts.same == 0);
  if (argc == 1) assert(counts.after > 0);
  std::cout << "before\tsame\tafter\tempty\tnonempty\n"
            << counts.before << '\t' << counts.same << '\t' << counts.after
            << '\t' << counts.empty << '\t' << counts.nonempty << '\n';
  return EXIT_SUCCESS;
}

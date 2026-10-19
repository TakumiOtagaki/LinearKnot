#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <cassert>
#include <cstdlib>
#include <string>

namespace {

template <class Semiring>
void check_widths(cparty::scfg2::replay::WFinalExactInsideContext &ctx,
                  cand_pos_t n) {
  using namespace cparty::scfg2;
  const auto run = [&](std::size_t width) {
    const auto result = replay::run_w_final_vp_wmbp_online<Semiring>(
        ctx, n, replay::WFinalVpWmbpOnlineWidths::uniform(width));
    return result.chart.get(deductions::ItemKey{
        grammar::NonTerminal::W, 1, n, -1, -1});
  };
  const double width20 = run(20);
  const double width40 = run(40);
  assert(width20 > 0.0);
  assert(width40 >= width20);
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  std::string sequence =
      "GGGUCUGAACGCAAGAACCGUUUUAAUAUGAAAUAAUGCGUUAACUGACUAUGGCGAAGGAUCAUCUGUCUCCGGAGAUCUUGAUUCUAUAGUGUUGUCCUAAGCAUCAUGAGGCUUCGCAGGUGCCAUCGGUUCGACAGCGUUCAAGGUUUGGACACCAAGGCCUGGGGGGGAGAAUGAUUCCAAUCUGGAUGGCGGUU";
  std::string scaffold =
      ".(...)(..(....)(....))(.(...)(.(....).))(..(.(....)))(....((.....))(...)(((...))(....(....)((...)(.....)))(.....))(.(....))(...)(...(...)((...)(...))((.....)((.(....))..)))(..(....)(.....(..(....)))))";
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  check_widths<cparty::scfg2::semirings::SumProductSemiring>(ctx, n);
  check_widths<cparty::scfg2::semirings::MaxProductSemiring>(ctx, n);
  return EXIT_SUCCESS;
}

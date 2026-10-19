#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/beam.hh"
#include "scfg2/api/exact.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

bool round_only_balanced(const std::string &structure) {
  int depth = 0;
  for (char symbol : structure) {
    if (symbol == '(') {
      ++depth;
    } else if (symbol == ')') {
      if (--depth < 0) return false;
    } else if (symbol != '.') {
      return false;
    }
  }
  return depth == 0;
}

}  // namespace

int main() {
  using namespace cparty::scfg2;
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);

  const std::string sequence = "UUAUAUGGGUAUUAUUCUCA";
  const std::string empty_scaffold(sequence.size(), '.');
  std::string owner_sequence = sequence;
  std::string owner_scaffold = empty_scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      owner_sequence, owner_scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());

  // With no fixed scaffold pair, every interval is weakly closed. Therefore
  // the public/open VP family that creates Square crossing carriers is absent.
  for (cand_pos_t i = 1; i <= n; ++i) {
    for (cand_pos_t j = i + 1; j <= n; ++j) {
      if (!ctx.weakly_closed(i, j) ||
          deductions::exact_round_open_family_active(ctx, i, j) ||
          deductions::exact_vp_public_item_active(
              deductions::ItemKey{
                  grammar::NonTerminal::VP, i, j, -1, -1},
              ctx)) {
        std::cerr << "empty-scaffold VP gate failed at "
                  << i << ',' << j << '\n';
        return EXIT_FAILURE;
      }
    }
  }

  api::ExactSession exact(sequence, empty_scaffold);
  const auto exact_best = exact.viterbi();
  api::BeamSession beam(
      sequence, empty_scaffold, api::BeamConfig{160, 30, 30});
  const auto beam_best = beam.viterbi();
  if (!exact_best.ok || !beam_best.ok || !beam_best.reachable ||
      !round_only_balanced(exact_best.generated_structure) ||
      !round_only_balanced(exact_best.union_structure) ||
      !round_only_balanced(beam_best.generated_structure) ||
      !round_only_balanced(beam_best.union_structure) ||
      beam_best.work.vp_candidates != 0 ||
      beam_best.work.wmbp_candidates != 0 ||
      beam_best.work.case3_candidates != 0 ||
      beam_best.work.vp_loop_parent_probes != 0) {
    std::cerr << "empty-scaffold traceback was not pseudoknot-free: exact="
              << exact_best.error << " beam=" << beam_best.error << '\n';
    return EXIT_FAILURE;
  }

  std::cout << "scfg2_empty_scaffold_pkfree_contract=ok\n";
  return EXIT_SUCCESS;
}

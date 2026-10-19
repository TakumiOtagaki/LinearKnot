#include "scfg2/api/exact.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct GoldenInsideCase {
  const char *sequence;
  const char *scaffold;
  double ensemble_free_energy;
  cparty::scfg2::WmbClosurePolicy policy =
      cparty::scfg2::WmbClosurePolicy::LegacyCompatible;
};

bool check_case(const GoldenInsideCase &golden) {
  cparty::scfg2::api::ExactSession session(golden.sequence, golden.scaffold,
      golden.policy);
  const auto inside = session.inside();
  const auto best = session.viterbi();
  if (!inside.ok || !best.ok || !(inside.partition_weight > 0.0) ||
      !(best.weight > 0.0)) {
    std::cerr << "exact run failed for " << golden.sequence
              << ": inside=" << inside.error << " viterbi=" << best.error
              << "\n";
    return false;
  }
  if (std::fabs(inside.ensemble_free_energy - golden.ensemble_free_energy) >
      1e-10) {
    std::cerr << "exact-inside golden mismatch for " << golden.sequence
              << ": got=" << inside.ensemble_free_energy
              << " expected=" << golden.ensemble_free_energy << "\n";
    return false;
  }
  if (best.weight > inside.partition_weight * (1.0 + 1e-12)) {
    std::cerr << "ViterbiWeight > InsideWeight for " << golden.sequence
              << "\n";
    return false;
  }

  const auto evaluated = session.evaluate_energy(best.generated_structure);
  const auto probability = session.probability(best.generated_structure);
  if (!evaluated.ok || !evaluated.reachable || !probability.ok ||
      !probability.reachable ||
      std::fabs(evaluated.energy - best.energy) > 1e-10 ||
      probability.probability < 0.0 ||
      probability.probability > 1.0 + 1e-12) {
    std::cerr << "cross-mode consistency failed for " << golden.sequence
              << ": evaluator=" << evaluated.error
              << " probability=" << probability.error << "\n";
    return false;
  }
  return true;
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);

  // These values were frozen from the broadly legacy-inside-aligned scfg2
  // checkpoint immediately before retiring src/scfg, including the
  // paper-faithful VP internal-loop border correction: CParty SI Case (5)
  // permits l == max(bp, B) + 1, while the legacy implementation accidentally
  // rejected it by combining a +1 shift with a strict comparison.
  // Values affected by multiloop closing neighbors were remeasured from
  // patched legacy including 0009. Weak-closure values use the independent
  // external legacy closure gate, not production's own inside result.
  const std::vector<GoldenInsideCase> cases = {
      {"AAAAAA", "......", 0.0},
      {"GAGAAACAC", ".........", -0.0073609393941053676},
      // `x` is pair-forbidden (-1), not an ordinary free endpoint (-2).
      // These values come from the frozen legacy W_final_pf::hfold_pf oracle;
      // the public legacy CLI cannot serve them because its preceding MFE
      // path crashes on `x` before inside is reached.
      {"GAGAAACAC", "x........", -0.0015876059375017433},
      // Patched legacy accepts uppercase X with exactly the same sentinel.
      {"GAGAAACAC", "X........", -0.0015876059375017433},
      {"GAGAAACAAAAAAGAAACC", "...................", -0.01555920819919634},
      // Independent corrected legacy: free and fully scaffolded boundary
      // multiloops. The latter has no pairable free bases, so Z is one target.
      {"GAGGGAAAACCCAGGGAAAACCCAC", ".........................", -7.2678525734899075},
      {"GAGGGAAAACCCAGGGAAAACCCAC", "(.(((....))).(((....))).)", -1.7000000000000013},
      // Patched legacy inside, full precision. VP(2,25) contains one maximal
      // G pair (11,17), but its WIP*VPR split k=10 has positive weight.
      // Contracting that sum to k=11 loses this term (the original checkpoint's
      // E shift was 4.87e-9, before the multiloop correction).
      {"GGGAAAAACGGAAAAACAAAACCAC", "(.........(.....)....)...",
       7.559778863810862},
      // Independent patched legacy including WIP 0008 and multiloop 0009.
      {"GAGAAACAAGAAGAAAAAAAAACCAAAAC", "(.(...)..(............)).....",
       11.897487329970824},
      {"UUAUAUGGGUAUUAUUCUCA", "((....))((...))(...)", 8.7625014673487396},
      {"UUAUAUGGGUAUUAUUCUCA", "((x...))((...))(...)", 8.9979393490826496},
      {"UUAUAUGGGUAUUAUUCUCA", "((X...))((...))(...)", 8.9979393490826496},
      {"UUAUAUGGGUAUUAUUCUCA", "((....))((.x.))(...)", 9.2060219269771366},
      {"UUAUAUGGGUAUUAUUCUCA", "((....))((.X.))(...)", 9.2060219269771366},
      // Long scaffold BE loops: current patched legacy inside, including the
      // explicit logarithmic PF extrapolation (0007), WIP boundary (0008)
      // and multiloop closing (0009).
      // The old out-of-bounds weights gave 14.408023917346283 instead.
      {"UCUACGUCAGUAGAUCUAACCGUGGAUAGCUUCUACAGGACAUGGGAUCCCAUUGGGUAAGACUU",
       "(.(X..(.X.X.X)..X).X.X..X..XX..X.XX.X.....X......X.X..).X(.XX)...",
       14.327075804769787},
      // Independent current patched legacy probe with only completed WMB on
      // non-weakly-closed intervals zeroed. Seed 20261008, rand50_0024 and
      // rand100_0098, remeasured after 0009; historical values are preserved
      // in results/wmb_closure_20261010/cases.tsv. See the WIP correction note.
      {"CCAUGAAAGGACCACUGUUUUGUCGUACUGGGCCCGUAGGGUUCCACGGU",
       "(...)(...(..(...)(..(...)(...))(.((....)).))(...))",
       37.8242027529809,
       cparty::scfg2::WmbClosurePolicy::RequireWeakClosure},
      {"CUCCCGAACUUGAAUGCUAGUUAUUAGCGUCGAGUGGGCUUUGGCGGUCAGAUAGUAUUCGCAAGCAAUAUUUGCAAGGAGGCUCGAACCGAGUUCAGUG",
       "....(..(.(...).(....)).((...)(.(....(...)))(.(.(...))(...)(...(.....)))(.....((...).)))((...)(...)))",
       96.64100479928429,
       cparty::scfg2::WmbClosurePolicy::RequireWeakClosure},
  };
  for (const auto &golden : cases) {
    if (!check_case(golden)) return EXIT_FAILURE;
  }

  std::cout << "scfg2_exact_inside_oracle=ok\n";
  return EXIT_SUCCESS;
}

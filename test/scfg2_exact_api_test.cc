#include "scfg2/api/exact.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);

  const std::string sequence = "UUAUAUGGGUAUUAUUCUCA";
  const std::string scaffold = "((....))((...))(...)";
  cparty::scfg2::api::ExactSession session(sequence, scaffold);

  const auto partition = session.inside();
  const auto outside = session.outside();
  const auto best = session.viterbi();
  if (!partition.ok || !outside.ok || !outside.reachable || !best.ok ||
      std::fabs(outside.log_partition_weight -
                std::log(partition.partition_weight)) > 1e-10 ||
      outside.unpaired_probabilities.size() != sequence.size() + 1) {
    std::cerr << "scfg2 exact API setup failed: inside=" << partition.error
              << " outside=" << outside.error
              << " viterbi=" << best.error << "\n";
    return EXIT_FAILURE;
  }
  std::vector<double> paired(sequence.size() + 1, 0.0);
  for (const auto &pair : outside.pair_probabilities) {
    paired[static_cast<std::size_t>(pair.i)] += pair.union_probability;
    paired[static_cast<std::size_t>(pair.j)] += pair.union_probability;
  }
  for (std::size_t i = 1; i <= sequence.size(); ++i) {
    if (std::fabs(
            paired[i] + outside.unpaired_probabilities[i] - 1.0) > 1e-9) {
      std::cerr << "scfg2 exact outside residue normalization failed\n";
      return EXIT_FAILURE;
    }
  }

  const auto evaluated = session.evaluate_energy(best.generated_structure);
  const auto probability = session.probability(best.generated_structure);
  if (!evaluated.ok || !evaluated.reachable ||
      !probability.ok || !probability.reachable ||
      std::fabs(evaluated.energy - best.energy) > 1e-10 ||
      probability.probability > 1.0 + 1e-12 ||
      evaluated.scaffold_structure != best.scaffold_structure ||
      evaluated.generated_structure != best.generated_structure ||
      evaluated.union_structure != best.union_structure ||
      probability.union_structure != best.union_structure) {
    std::cerr << "scfg2 exact API consistency failed: evaluator="
              << evaluated.error << " probability=" << probability.error
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "scfg2_exact_api=ok\n";
  return EXIT_SUCCESS;
}

#include "scfg2/api/beam.hh"
#include "scfg2/api/exact.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

bool close(double lhs, double rhs, double tolerance = 1e-10) {
  const double scale = std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  return std::fabs(lhs - rhs) <= tolerance * scale;
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);

  const std::string sequence = "UUAUAUGGGUAUUAUUCUCA";
  const std::string scaffold = "((....))((...))(...)";
  const cparty::scfg2::api::BeamConfig config{160, 30, 30};
  cparty::scfg2::api::BeamSession beam(sequence, scaffold, config);

  const auto inside = beam.inside();
  const auto cached_inside = beam.inside();
  const auto outside = beam.outside();
  const auto best = beam.viterbi();
  if (!inside.ok || !inside.reachable ||
      !cached_inside.ok || !cached_inside.reachable ||
      !outside.ok || !outside.reachable ||
      !best.ok || !best.reachable ||
      !std::isfinite(inside.log_partition_weight) ||
      !std::isfinite(inside.ensemble_free_energy) ||
      inside.elapsed_seconds <= 0.0 ||
      inside.work.total_item_visits == 0 ||
      inside.work.families.empty() ||
      !std::isfinite(best.log_weight) ||
      !std::isfinite(best.energy) ||
      inside.beam_size != config.beam_size ||
      outside.beam_size != config.beam_size ||
      best.beam_size != config.beam_size ||
      inside.log_partition_weight != cached_inside.log_partition_weight ||
      !close(inside.log_partition_weight, outside.log_partition_weight) ||
      outside.unpaired_probabilities.size() != sequence.size() + 1 ||
      outside.work.total_item_visits != inside.work.total_item_visits ||
      best.scaffold_structure != scaffold ||
      best.generated_structure.size() != sequence.size() ||
      best.union_structure.size() != sequence.size() ||
      best.work.families.empty() ||
      best.work.total_item_visits == 0 ||
      best.work.total_frontier_candidates <
          best.work.total_frontier_selected) {
    std::cerr << "beam API setup failed: inside=" << inside.error
              << " outside=" << outside.error
              << " viterbi=" << best.error << '\n';
    return EXIT_FAILURE;
  }

  cparty::scfg2::api::ExactSession exact(sequence, scaffold);
  const auto exact_inside = exact.inside();
  const auto exact_outside = exact.outside();
  const auto exact_best = exact.viterbi();
  if (!exact_inside.ok || !exact_outside.ok || !exact_best.ok ||
      !close(inside.log_partition_weight,
             std::log(exact_inside.partition_weight)) ||
      outside.pair_probabilities.size() !=
          exact_outside.pair_probabilities.size() ||
      !close(best.energy, exact_best.energy) ||
      best.generated_structure != exact_best.generated_structure ||
      best.union_structure != exact_best.union_structure) {
    std::cerr << "high-beam/exact API parity failed\n";
    return EXIT_FAILURE;
  }
  for (std::size_t index = 0;
       index < outside.pair_probabilities.size();
       ++index) {
    const auto &beam_pair = outside.pair_probabilities[index];
    const auto &exact_pair = exact_outside.pair_probabilities[index];
    if (beam_pair.i != exact_pair.i || beam_pair.j != exact_pair.j ||
        !close(beam_pair.round_probability,
               exact_pair.round_probability) ||
        !close(beam_pair.square_probability,
               exact_pair.square_probability) ||
        !close(beam_pair.union_probability,
               exact_pair.union_probability)) {
      std::cerr << "high-beam/exact outside marginal parity failed\n";
      return EXIT_FAILURE;
    }
  }

  const auto probability = beam.probability(best.generated_structure);
  const auto pair_probability = beam.probability(best.generated_pairs);
  const auto union_probability = beam.probability(best.union_structure);
  if (!probability.ok || !probability.denominator_reachable ||
      !probability.target_included ||
      !pair_probability.ok || !pair_probability.denominator_reachable ||
      !pair_probability.target_included ||
      !union_probability.ok || !union_probability.denominator_reachable ||
      !union_probability.target_included ||
      probability.probability < 0.0 ||
      probability.probability > 1.0 + 1e-12 ||
      probability.log_probability > 1e-12 ||
      !close(probability.log_probability,
             pair_probability.log_probability) ||
      !close(probability.log_probability,
             union_probability.log_probability) ||
      probability.target_generated_pairs.size() !=
          best.generated_pairs.size()) {
    std::cerr << "beam API probability failed: string="
              << probability.error << " pairs=" << pair_probability.error
              << '\n';
    return EXIT_FAILURE;
  }

  const auto absent = beam.probability(
      std::vector<cparty::scfg2::api::BasePair>{{1, 20}});
  if (!absent.ok || !absent.denominator_reachable ||
      absent.target_included || absent.probability != 0.0) {
    std::cerr << "beam API fail-closed target contract failed: "
              << absent.error << '\n';
    return EXIT_FAILURE;
  }

  cparty::scfg2::api::BeamSession invalid(sequence, scaffold, 0);
  const auto invalid_inside = invalid.inside();
  if (invalid_inside.ok ||
      invalid_inside.error != "beam size must be positive") {
    std::cerr << "beam API invalid-width contract failed\n";
    return EXIT_FAILURE;
  }

  const auto free_best =
      cparty::scfg2::api::beam_viterbi(sequence, scaffold, config);
  if (!free_best.ok || !free_best.reachable ||
      free_best.generated_structure != best.generated_structure) {
    std::cerr << "beam API free-function wrapper failed\n";
    return EXIT_FAILURE;
  }

  if (inside.work.case3_truncated_endpoints != 0 ||
      best.work.case3_truncated_endpoints != 0) {
    std::cerr << "beam API default depth cap unexpectedly truncated\n";
    return EXIT_FAILURE;
  }

  // Two branching levels: each loop inside the inner multiloop has two
  // branching ancestors, so a cap of one cuts all 8 of those positions.
  std::string nested_scaffold = "((....)((....)(....)))";
  std::string nested_sequence(nested_scaffold.size(), 'A');
  for (std::size_t pos = 0; pos < nested_scaffold.size(); ++pos) {
    if (nested_scaffold[pos] == '(') nested_sequence[pos] = 'G';
    if (nested_scaffold[pos] == ')') nested_sequence[pos] = 'C';
  }
  const auto capped_inside =
      cparty::scfg2::api::BeamSession(
          nested_sequence, nested_scaffold,
          cparty::scfg2::api::BeamConfig{20, 30, 1})
          .inside();
  const auto uncapped_inside =
      cparty::scfg2::api::BeamSession(
          nested_sequence, nested_scaffold,
          cparty::scfg2::api::BeamConfig{20, 30, 30})
          .inside();
  if (!capped_inside.ok || !uncapped_inside.ok ||
      capped_inside.work.case3_truncated_endpoints != 8 ||
      uncapped_inside.work.case3_truncated_endpoints != 0) {
    std::cerr << "beam API depth-cap truncation count failed: capped="
              << capped_inside.work.case3_truncated_endpoints
              << " uncapped=" << uncapped_inside.work.case3_truncated_endpoints
              << '\n';
    return EXIT_FAILURE;
  }

  std::cout << "scfg2_beam_api=ok\n";
  return EXIT_SUCCESS;
}

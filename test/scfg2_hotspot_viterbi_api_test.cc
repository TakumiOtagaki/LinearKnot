#include "scfg2/api/hotspot_viterbi.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using cparty::scfg2::api::BasePair;
using cparty::scfg2::api::ScaffoldViterbiRecord;

ScaffoldViterbiRecord record(std::string scaffold,
                             std::string union_structure,
                             double energy,
                             std::size_t pass,
                             std::uint64_t origin) {
  ScaffoldViterbiRecord out;
  out.proposal.pass = pass;
  out.proposal.origin = origin;
  out.scaffold_structure = std::move(scaffold);
  out.viterbi.ok = true;
  out.viterbi.reachable = true;
  out.viterbi.energy = energy;
  out.viterbi.union_structure = std::move(union_structure);
  return out;
}

bool same_pair(const BasePair &pair, int i, int j) {
  return pair.i == i && pair.j == j;
}

}  // namespace

int main() {
  std::vector<ScaffoldViterbiRecord> synthetic;
  synthetic.push_back(record("(......)", "((....))", -4.0, 2, 7));
  synthetic.push_back(record(".(....).", "((....))", -5.0, 1, 8));
  synthetic.push_back(record("..(....)", "..(....)", -5.0, 0, 9));
  const auto grouped =
      cparty::scfg2::api::select_minimum_energy_records(
          std::move(synthetic));
  if (!grouped.ok || !grouped.reachable ||
      grouped.records.size() != 3 || grouped.structures.size() != 2) {
    std::cerr << "synthetic grouping setup failed: " << grouped.error << '\n';
    return EXIT_FAILURE;
  }
  const auto &selected = grouped.structures[grouped.selected_structure];
  if (selected.completed_pairs.size() != 2 ||
      !same_pair(selected.completed_pairs[0], 1, 8) ||
      !same_pair(selected.completed_pairs[1], 2, 7) ||
      selected.record_indices.size() != 2 ||
      selected.selected_record != 1 ||
      selected.minimum_energy != -5.0 ||
      selected.maximum_energy != -4.0 ||
      selected.energy_spread != 1.0) {
    std::cerr << "equal-R grouping or deterministic selection failed\n";
    return EXIT_FAILURE;
  }

  auto failed_record = record("........", "........", 0.0, 0, 0);
  failed_record.viterbi.ok = false;
  failed_record.viterbi.error = "injected";
  const auto failed =
      cparty::scfg2::api::select_minimum_energy_records({failed_record});
  if (failed.ok || failed.error.find("injected") == std::string::npos) {
    std::cerr << "record failure did not fail closed\n";
    return EXIT_FAILURE;
  }

  auto unreachable_record = record("........", "........", 0.0, 0, 0);
  unreachable_record.viterbi.reachable = false;
  const auto unreachable =
      cparty::scfg2::api::select_minimum_energy_records(
          {unreachable_record});
  if (!unreachable.ok || unreachable.reachable ||
      !unreachable.structures.empty()) {
    std::cerr << "unreachable-record contract failed\n";
    return EXIT_FAILURE;
  }

  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  cparty::scfg2::proposal::WavefrontHotspotConfig proposal_config;
  proposal_config.anchor_cap = 2;
  proposal_config.beam_width = 4;
  proposal_config.passes = 1;
  proposal_config.max_results = 2;
  proposal_config.diversity_mode =
      cparty::scfg2::proposal::HotspotDiversityMode::None;
  cparty::scfg2::api::BeamConfig beam_config{64, 30, 30};
  const std::string sequence = "UUAUAUGGGUAUUAUUCUCA";
  const auto predicted = cparty::scfg2::api::predict_hotspot_viterbi(
      sequence, proposal_config, beam_config);
  if (!predicted.ok || !predicted.reachable ||
      predicted.records.empty() || predicted.records.size() > 2 ||
      predicted.structures.empty() ||
      predicted.proposal_elapsed_seconds < 0.0 ||
      predicted.viterbi_elapsed_seconds <= 0.0 ||
      predicted.total_elapsed_seconds <
          predicted.viterbi_elapsed_seconds) {
    std::cerr << "hotspot/Viterbi integration failed: "
              << predicted.error << '\n';
    return EXIT_FAILURE;
  }

  bool saw_empty_scaffold = false;
  double minimum = std::numeric_limits<double>::infinity();
  for (const auto &candidate : predicted.records) {
    if (candidate.scaffold_structure == std::string(sequence.size(), '.')) {
      saw_empty_scaffold = true;
    }
    if (candidate.viterbi.beam_size != beam_config.beam_size ||
        !candidate.viterbi.ok) {
      std::cerr << "candidate did not use the common BeamConfig\n";
      return EXIT_FAILURE;
    }
    if (candidate.viterbi.reachable) {
      minimum = std::min(minimum, candidate.viterbi.energy);
    }
  }
  const auto &winner =
      predicted.structures[predicted.selected_structure];
  if (!saw_empty_scaffold || !std::isfinite(minimum) ||
      winner.minimum_energy != minimum) {
    std::cerr << "minimum-energy winner contract failed\n";
    return EXIT_FAILURE;
  }

  std::cout << "scfg2_hotspot_viterbi_api=ok"
            << " records=" << predicted.records.size()
            << " groups=" << predicted.structures.size()
            << " energy=" << winner.minimum_energy << '\n';
  return EXIT_SUCCESS;
}

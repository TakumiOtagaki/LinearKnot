#ifndef CPARTY_SCFG2_API_HOTSPOT_VITERBI_HH
#define CPARTY_SCFG2_API_HOTSPOT_VITERBI_HH

#include "scfg2/api/beam.hh"
#include "scfg2/proposal/hotspot_wavefront.hh"

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace cparty::scfg2::api {

struct ScaffoldViterbiRecord {
  proposal::HotspotProposal proposal;
  std::string scaffold_structure;
  BeamViterbiResult viterbi;
  double elapsed_seconds = 0.0;
};

struct CompletedStructureGroup {
  std::vector<BasePair> completed_pairs;
  std::string representative_union_structure;
  std::vector<std::size_t> record_indices;
  double minimum_energy = std::numeric_limits<double>::infinity();
  double maximum_energy = -std::numeric_limits<double>::infinity();
  double energy_spread = std::numeric_limits<double>::quiet_NaN();
  std::size_t selected_record = 0;
};

struct HotspotViterbiPrediction {
  bool ok = false;
  bool reachable = false;
  std::string error;
  proposal::WavefrontHotspotStats proposal_stats;
  double proposal_elapsed_seconds = 0.0;
  double viterbi_elapsed_seconds = 0.0;
  double total_elapsed_seconds = 0.0;
  std::vector<ScaffoldViterbiRecord> records;
  std::vector<CompletedStructureGroup> structures;
  std::size_t selected_structure = 0;
};

// Groups precomputed records by completed pair map R and selects the minimum
// finite CParty energy. This pure aggregation boundary is also useful to
// consumers that supply their own scaffold proposal source.
HotspotViterbiPrediction select_minimum_energy_records(
    std::vector<ScaffoldViterbiRecord> records);

// Generates bounded pseudoknot-free scaffold candidates, runs exactly one
// BeamSession::viterbi() with the same beam configuration for every candidate,
// groups equal completed pair maps R, and selects the observed minimum energy.
HotspotViterbiPrediction predict_hotspot_viterbi(
    const std::string &sequence,
    proposal::WavefrontHotspotConfig proposal_config = {},
    BeamConfig beam_config = {});

}  // namespace cparty::scfg2::api

#endif  // CPARTY_SCFG2_API_HOTSPOT_VITERBI_HH

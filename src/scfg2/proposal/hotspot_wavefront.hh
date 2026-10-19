#ifndef CPARTY_SCFG2_PROPOSAL_HOTSPOT_WAVEFRONT_HH
#define CPARTY_SCFG2_PROPOSAL_HOTSPOT_WAVEFRONT_HH

#include "scfg2/api/types.hh"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cparty::scfg2::proposal {

enum class HotspotDiversityMode {
  None,
  DeterministicPassDiversification,
};

struct WavefrontHotspotConfig {
  std::size_t anchor_cap = 8;
  std::size_t beam_width = 16;
  std::size_t passes = 8;
  std::size_t max_results = 24;
  std::size_t seed_pairs = 3;
  std::size_t minimum_pairs = 3;
  std::size_t minimum_loop_unpaired = 3;
  bool allow_gu = true;
  std::uint64_t seed = 1;

  // The initial scorer is intentionally inexpensive. Larger is better.
  double gc_pair_score = 3.0;
  double au_pair_score = 2.0;
  double gu_pair_score = 1.0;
  double stack_bonus = 1.0;
  double bulge_penalty = 1.0;
  double internal_loop_penalty = 2.0;
  double minimum_score = 0.0;

  // Diversity changes bounded anchor representatives and beam-retention
  // priority between passes. It never changes heuristic_score or the CParty
  // energy used by the final repeated-Viterbi selector.
  HotspotDiversityMode diversity_mode =
      HotspotDiversityMode::DeterministicPassDiversification;
  double diversity_priority_jitter_amplitude = 0.25;
};

struct HotspotProposal {
  double heuristic_score = 0.0;
  std::vector<api::BasePair> scaffold_pairs;
  std::size_t pass = 0;
  std::uint64_t origin = 0;
};

struct WavefrontHotspotStats {
  std::size_t boundaries = 0;
  std::size_t anchors_retained = 0;
  std::size_t successor_candidates = 0;
  std::size_t states_retained = 0;
  std::size_t completed_candidates = 0;
};

struct WavefrontHotspotResult {
  bool ok = false;
  std::string error;
  std::vector<HotspotProposal> proposals;
  WavefrontHotspotStats stats;
};

// HotKnots-inspired sparse stem proposal. The production contract is bounded
// by anchor_cap and beam_width; it never enumerates all complementary pairs.
WavefrontHotspotResult propose_hotspots_wavefront(
    const std::string &sequence,
    const WavefrontHotspotConfig &config = {});

}  // namespace cparty::scfg2::proposal

#endif  // CPARTY_SCFG2_PROPOSAL_HOTSPOT_WAVEFRONT_HH

#include "scfg2/proposal/hotspot_wavefront.hh"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

std::size_t parse_size(const char *text, const char *name) {
  char *end = nullptr;
  const auto value = std::strtoull(text, &end, 10);
  if (!text[0] || !end || *end != '\0' || value == 0) {
    std::cerr << name << " must be a positive integer\n";
    std::exit(EXIT_FAILURE);
  }
  return static_cast<std::size_t>(value);
}

std::string render_scaffold(
    std::size_t n,
    const std::vector<cparty::scfg2::api::BasePair> &pairs) {
  std::string out(n, '.');
  for (const auto &pair : pairs) {
    if (pair.i < 1 || pair.j <= pair.i ||
        static_cast<std::size_t>(pair.j) > n ||
        out[static_cast<std::size_t>(pair.i - 1)] != '.' ||
        out[static_cast<std::size_t>(pair.j - 1)] != '.') {
      return {};
    }
    out[static_cast<std::size_t>(pair.i - 1)] = '(';
    out[static_cast<std::size_t>(pair.j - 1)] = ')';
  }
  return out;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 7) {
    std::cerr << "usage: scfg2_hotspot_proposal_tool SEQUENCE"
              << " ANCHOR_CAP HOTSPOT_WIDTH PASSES MAX_RESULTS"
              << " DIVERSITY:none|deterministic\n";
    return EXIT_FAILURE;
  }
  const std::string sequence = argv[1];
  cparty::scfg2::proposal::WavefrontHotspotConfig config;
  config.anchor_cap = parse_size(argv[2], "anchor cap");
  config.beam_width = parse_size(argv[3], "hotspot width");
  config.passes = parse_size(argv[4], "passes");
  config.max_results = parse_size(argv[5], "max results");
  const std::string diversity = argv[6];
  if (diversity == "none") {
    config.diversity_mode =
        cparty::scfg2::proposal::HotspotDiversityMode::None;
  } else if (diversity == "deterministic") {
    config.diversity_mode =
        cparty::scfg2::proposal::HotspotDiversityMode::
            DeterministicPassDiversification;
  } else {
    std::cerr << "diversity must be none or deterministic\n";
    return EXIT_FAILURE;
  }

  const auto begin = std::chrono::steady_clock::now();
  const auto proposed =
      cparty::scfg2::proposal::propose_hotspots_wavefront(sequence, config);
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - begin).count();
  std::cout << std::setprecision(12)
            << "SUMMARY\tok=" << proposed.ok
            << "\trecords=" << proposed.proposals.size()
            << "\tproposal_seconds=" << seconds
            << "\terror=" << proposed.error << '\n'
            << "PROPOSAL_HEADER\tindex\tscaffold_pairs\theuristic_score"
            << "\tpass\torigin\tscaffold_G\n";
  for (std::size_t index = 0; index < proposed.proposals.size(); ++index) {
    const auto &proposal = proposed.proposals[index];
    const auto scaffold =
        render_scaffold(sequence.size(), proposal.scaffold_pairs);
    if (scaffold.size() != sequence.size()) {
      std::cerr << "proposal contains an invalid scaffold pair map\n";
      return EXIT_FAILURE;
    }
    std::cout << "PROPOSAL\t" << index
              << '\t' << proposal.scaffold_pairs.size()
              << '\t' << proposal.heuristic_score
              << '\t' << proposal.pass
              << '\t' << proposal.origin
              << '\t' << scaffold << '\n';
  }
  return proposed.ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

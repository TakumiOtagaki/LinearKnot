#include "scfg2/api/beam.hh"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
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

}  // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: scfg2_fixed_scaffold_viterbi_tool"
              << " SEQUENCE PARAM_FILE BEAM_WIDTH SCAFFOLD\n";
    return EXIT_FAILURE;
  }
  const std::string sequence = argv[1];
  const std::string parameter_file = argv[2];
  const std::size_t beam_width = parse_size(argv[3], "beam width");
  const std::string scaffold = argv[4];
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", parameter_file.c_str(), 1);

  const auto begin = std::chrono::steady_clock::now();
  cparty::scfg2::api::BeamSession session(
      sequence, scaffold,
      cparty::scfg2::api::BeamConfig{beam_width, 30, 30});
  const auto result = session.viterbi();
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - begin).count();
  std::cout << std::setprecision(12)
            << "SUMMARY\tok=" << result.ok
            << "\treachable=" << result.reachable
            << "\tenergy="
            << (result.reachable
                    ? result.energy
                    : std::numeric_limits<double>::infinity())
            << "\tseconds=" << seconds
            << "\tgenerated_Gprime=" << result.generated_structure
            << "\tunion_R=" << result.union_structure
            << "\tcase3_truncated_endpoints="
            << result.work.case3_truncated_endpoints
            << "\terror=" << result.error << '\n';
  return result.ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

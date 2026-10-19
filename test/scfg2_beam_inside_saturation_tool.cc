#include "scfg2/api/beam.hh"
#include "scfg2/api/exact.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::size_t> parse_widths(const std::string &text) {
  std::vector<std::size_t> out;
  std::stringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (!token.empty()) out.push_back(std::stoull(token));
  }
  return out;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4 || argc > 6) {
    std::cerr << "usage: scfg2_beam_inside_saturation_tool"
              << " SEQUENCE SCAFFOLD|empty WIDTHS_CSV"
              << " [PARAM_FILE] [exact|noexact]\n";
    return EXIT_FAILURE;
  }
  const std::string sequence = argv[1];
  std::string scaffold = argv[2];
  if (scaffold == "empty") scaffold = std::string(sequence.size(), '.');
  const auto widths = parse_widths(argv[3]);
  if (sequence.empty() || scaffold.size() != sequence.size() ||
      widths.empty() || widths.front() == 0 ||
      !std::is_sorted(widths.begin(), widths.end()) ||
      std::adjacent_find(widths.begin(), widths.end()) != widths.end()) {
    std::cerr << "invalid sequence, scaffold, or increasing width list\n";
    return EXIT_FAILURE;
  }
  setenv(
      "CPARTY_FIXED_ENERGY_PARAM_FILE",
      argc >= 5 ? argv[4] : "params/rna_Turner04.par",
      1);
  const bool run_exact = argc == 6 && std::string(argv[5]) == "exact";
  if (argc == 6 && !run_exact && std::string(argv[5]) != "noexact") {
    std::cerr << "last argument must be exact or noexact\n";
    return EXIT_FAILURE;
  }

  double exact_log_z = std::numeric_limits<double>::quiet_NaN();
  double exact_free_energy = std::numeric_limits<double>::quiet_NaN();
  if (run_exact) {
    cparty::scfg2::api::ExactSession exact(sequence, scaffold);
    const auto result = exact.inside();
    if (!result.ok || !result.reachable) {
      std::cerr << "exact inside failed: " << result.error << '\n';
      return EXIT_FAILURE;
    }
    exact_log_z = result.log_partition_weight;
    exact_free_energy = result.ensemble_free_energy;
  }

  std::cout << std::setprecision(17)
            << "n\twidth\treachable\tlog_z\tfree_energy"
               "\tdelta_log_z\tdelta_free_energy\tabs_delta_free_energy"
               "\tmonotone_from_previous\texact_log_z\texact_free_energy"
               "\tlog_z_gap_to_exact\tseconds\titem_visits"
               "\tfrontier_candidates\tfrontier_selected\tvp_candidates"
               "\twmbp_candidates\tcase3_candidates\tbinary_join_probes"
               "\textension_probes\tv_parent_probes"
               "\tvp_loop_parent_probes\tcase3_truncated_endpoints\n";

  double previous_log_z = -std::numeric_limits<double>::infinity();
  double previous_free_energy = std::numeric_limits<double>::quiet_NaN();
  for (std::size_t width : widths) {
    std::cerr << "inside width=" << width << " start\n";
    cparty::scfg2::api::BeamSession session(
        sequence,
        scaffold,
        cparty::scfg2::api::BeamConfig{width, 30, 30});
    const auto result = session.inside();
    if (!result.ok) {
      std::cerr << "beam inside failed at width " << width
                << ": " << result.error << '\n';
      return EXIT_FAILURE;
    }
    std::cerr << "inside width=" << width
              << " done seconds=" << result.elapsed_seconds << '\n';
    const bool has_previous = std::isfinite(previous_log_z);
    const double delta_log_z = has_previous
        ? result.log_partition_weight - previous_log_z
        : std::numeric_limits<double>::quiet_NaN();
    const double delta_free_energy = has_previous
        ? result.ensemble_free_energy - previous_free_energy
        : std::numeric_limits<double>::quiet_NaN();
    const bool monotone =
        !has_previous ||
        result.log_partition_weight + 1e-12 >= previous_log_z;
    const auto &work = result.work;
    std::cout << sequence.size() << '\t' << width << '\t'
              << result.reachable << '\t'
              << result.log_partition_weight << '\t'
              << result.ensemble_free_energy << '\t'
              << delta_log_z << '\t' << delta_free_energy << '\t'
              << std::fabs(delta_free_energy) << '\t'
              << monotone << '\t' << exact_log_z << '\t'
              << exact_free_energy << '\t'
              << (run_exact
                      ? exact_log_z - result.log_partition_weight
                      : std::numeric_limits<double>::quiet_NaN())
              << '\t' << result.elapsed_seconds
              << '\t' << work.total_item_visits
              << '\t' << work.total_frontier_candidates
              << '\t' << work.total_frontier_selected
              << '\t' << work.vp_candidates
              << '\t' << work.wmbp_candidates
              << '\t' << work.case3_candidates
              << '\t' << work.binary_join_probes
              << '\t' << work.extension_probes
              << '\t' << work.v_parent_probes
              << '\t' << work.vp_loop_parent_probes
              << '\t' << work.case3_truncated_endpoints << '\n';
    previous_log_z = result.log_partition_weight;
    previous_free_energy = result.ensemble_free_energy;
  }
  return EXIT_SUCCESS;
}

#include "scfg2/api/beam.hh"

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

struct Difference {
  double pair_max = 0.0;
  double pair_mean = 0.0;
  std::size_t pair_max_i = 0;
  std::size_t pair_max_j = 0;
  double unpaired_max = 0.0;
  double unpaired_mean = 0.0;
  std::size_t unpaired_max_position = 0;
};

Difference compare(
    const cparty::scfg2::api::BeamOutsideResult &lhs,
    const cparty::scfg2::api::BeamOutsideResult &rhs,
    std::size_t n) {
  std::vector<double> lhs_pairs((n + 1) * (n + 1), 0.0);
  std::vector<double> rhs_pairs((n + 1) * (n + 1), 0.0);
  const auto fill = [n](const auto &pairs, auto *matrix) {
    for (const auto &pair : pairs) {
      (*matrix)[static_cast<std::size_t>(pair.i) * (n + 1) +
                static_cast<std::size_t>(pair.j)] =
          pair.union_probability;
    }
  };
  fill(lhs.pair_probabilities, &lhs_pairs);
  fill(rhs.pair_probabilities, &rhs_pairs);

  Difference out;
  std::size_t pair_count = 0;
  for (std::size_t i = 1; i <= n; ++i) {
    for (std::size_t j = i + 1; j <= n; ++j) {
      const auto index = i * (n + 1) + j;
      const double delta = std::fabs(lhs_pairs[index] - rhs_pairs[index]);
      if (delta > out.pair_max) {
        out.pair_max = delta;
        out.pair_max_i = i;
        out.pair_max_j = j;
      }
      out.pair_mean += delta;
      ++pair_count;
    }
    const double delta = std::fabs(
        lhs.unpaired_probabilities[i] - rhs.unpaired_probabilities[i]);
    if (delta > out.unpaired_max) {
      out.unpaired_max = delta;
      out.unpaired_max_position = i;
    }
    out.unpaired_mean += delta;
  }
  if (pair_count > 0) {
    out.pair_mean /= static_cast<double>(pair_count);
  }
  if (n > 0) out.unpaired_mean /= static_cast<double>(n);
  return out;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4 || argc > 5) {
    std::cerr << "usage: scfg2_beam_outside_saturation_tool"
              << " SEQUENCE SCAFFOLD|empty WIDTHS_CSV [PARAM_FILE]\n";
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
      argc == 5 ? argv[4] : "params/rna_Turner04.par",
      1);

  std::vector<cparty::scfg2::api::BeamOutsideResult> results;
  results.reserve(widths.size());
  for (const auto width : widths) {
    std::cerr << "outside width=" << width << " start\n";
    cparty::scfg2::api::BeamSession session(
        sequence,
        scaffold,
        cparty::scfg2::api::BeamConfig{width, 30, 30});
    auto result = session.outside();
    if (!result.ok || !result.reachable) {
      std::cerr << "beam outside failed at width " << width
                << ": " << result.error << '\n';
      return EXIT_FAILURE;
    }
    std::cerr << "outside width=" << width
              << " done seconds=" << result.elapsed_seconds
              << " reverse_records=" << result.reverse_records << '\n';
    results.push_back(std::move(result));
  }

  std::cout << std::setprecision(17)
            << "n\twidth\tlog_z\tfree_energy\tseconds\treverse_records"
               "\tpair_count\tpair_max_vs_previous"
               "\tpair_max_i_vs_previous\tpair_max_j_vs_previous"
               "\tpair_mean_vs_previous\tunpaired_max_vs_previous"
               "\tunpaired_max_position_vs_previous"
               "\tunpaired_mean_vs_previous\tpair_max_vs_largest"
               "\tpair_max_i_vs_largest\tpair_max_j_vs_largest"
               "\tpair_mean_vs_largest\tunpaired_max_vs_largest"
               "\tunpaired_max_position_vs_largest"
               "\tunpaired_mean_vs_largest\titem_visits"
               "\tfrontier_candidates\tfrontier_selected\n";
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  const auto &largest = results.back();
  for (std::size_t index = 0; index < results.size(); ++index) {
    const auto &result = results[index];
    const Difference previous =
        index == 0 ? Difference{nan, nan, 0, 0, nan, nan, 0}
                   : compare(result, results[index - 1], sequence.size());
    const Difference reference = compare(result, largest, sequence.size());
    std::cout << sequence.size() << '\t' << widths[index] << '\t'
              << result.log_partition_weight << '\t'
              << result.ensemble_free_energy << '\t'
              << result.elapsed_seconds << '\t'
              << result.reverse_records << '\t'
              << result.pair_probabilities.size() << '\t'
              << previous.pair_max << '\t'
              << previous.pair_max_i << '\t'
              << previous.pair_max_j << '\t'
              << previous.pair_mean << '\t'
              << previous.unpaired_max << '\t'
              << previous.unpaired_max_position << '\t'
              << previous.unpaired_mean << '\t'
              << reference.pair_max << '\t'
              << reference.pair_max_i << '\t'
              << reference.pair_max_j << '\t'
              << reference.pair_mean << '\t'
              << reference.unpaired_max << '\t'
              << reference.unpaired_max_position << '\t'
              << reference.unpaired_mean << '\t'
              << result.work.total_item_visits << '\t'
              << result.work.total_frontier_candidates << '\t'
              << result.work.total_frontier_selected << '\n';
  }
  return EXIT_SUCCESS;
}

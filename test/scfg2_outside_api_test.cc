#include "scfg2/api/beam.hh"
#include "scfg2/api/exact.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

bool close(double lhs, double rhs, double tolerance = 1e-9) {
  return std::fabs(lhs - rhs) <=
         tolerance * std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
}

const cparty::scfg2::api::PairProbability *find_pair(
    const std::vector<cparty::scfg2::api::PairProbability> &pairs,
    int i,
    int j) {
  for (const auto &pair : pairs) {
    if (pair.i == i && pair.j == j) return &pair;
  }
  return nullptr;
}

bool check_normalization(const cparty::scfg2::api::OutsideResult &outside,
                         std::size_t n) {
  std::vector<double> paired(n + 1, 0.0);
  for (const auto &pair : outside.pair_probabilities) {
    paired[static_cast<std::size_t>(pair.i)] += pair.union_probability;
    paired[static_cast<std::size_t>(pair.j)] += pair.union_probability;
  }
  for (std::size_t i = 1; i <= n; ++i) {
    if (!close(paired[i] + outside.unpaired_probabilities[i], 1.0)) {
      return false;
    }
  }
  return true;
}

bool check_case(const std::string &sequence,
                const std::string &scaffold,
                bool expect_square) {
  cparty::scfg2::api::ExactSession exact(sequence, scaffold);
  const auto exact_outside = exact.outside();
  cparty::scfg2::api::BeamSession beam(
      sequence, scaffold, cparty::scfg2::api::BeamConfig{160, 30, 30});
  const auto beam_outside = beam.outside();
  if (!exact_outside.ok || !exact_outside.reachable ||
      !beam_outside.ok || !beam_outside.reachable ||
      !(beam_outside.elapsed_seconds > 0.0) ||
      beam_outside.reverse_records == 0 || beam_outside.tape_nodes != 0 ||
      beam_outside.reverse_terms != 0 || beam_outside.reverse_factors != 0 ||
      !close(exact_outside.log_partition_weight,
             beam_outside.log_partition_weight) ||
      exact_outside.pair_probabilities.size() !=
          beam_outside.pair_probabilities.size() ||
      !check_normalization(exact_outside, sequence.size())) {
    std::cerr << "outside setup/parity failed sequence=" << sequence
              << " exact=" << exact_outside.error
              << " beam=" << beam_outside.error << '\n';
    return false;
  }
  bool has_square = false;
  for (std::size_t index = 0;
       index < exact_outside.pair_probabilities.size();
       ++index) {
    const auto &lhs = exact_outside.pair_probabilities[index];
    const auto &rhs = beam_outside.pair_probabilities[index];
    has_square = has_square || lhs.square_probability > 1e-12;
    if (lhs.i != rhs.i || lhs.j != rhs.j ||
        !close(lhs.round_probability, rhs.round_probability) ||
        !close(lhs.square_probability, rhs.square_probability) ||
        !close(lhs.union_probability, rhs.union_probability)) {
      std::cerr << "outside exact/beam pair mismatch\n";
      return false;
    }
  }
  if (has_square != expect_square) {
    std::cerr << "outside Square-lane expectation failed sequence="
              << sequence << " expected=" << expect_square
              << " observed=" << has_square << '\n';
    return false;
  }
  return true;
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);

  if (!check_case("GGGAAACCCU", "..........", false)) {
    return EXIT_FAILURE;
  }

  const std::string round_sequence = "UGUGGCAAUUUGAAAUCUAA";
  const std::string round_scaffold = "...(..(...).....)...";
  if (!check_case(round_sequence, round_scaffold, true)) {
    return EXIT_FAILURE;
  }
  cparty::scfg2::api::ExactSession round_exact(
      round_sequence, round_scaffold);
  const auto round_outside = round_exact.outside();
  for (const auto coordinate :
       {std::pair<int, int>{5, 16}, std::pair<int, int>{6, 12}}) {
    const auto *pair =
        find_pair(round_outside.pair_probabilities,
                  coordinate.first,
                  coordinate.second);
    if (!pair || !(pair->round_probability > 0.0) ||
        pair->square_probability != 0.0) {
      std::cerr << "generated V pair was not Round-only\n";
      return EXIT_FAILURE;
    }
  }
  const auto *fixed = find_pair(round_outside.pair_probabilities, 4, 17);
  if (!fixed || fixed->round_probability != 1.0 ||
      fixed->square_probability != 0.0 ||
      fixed->union_probability != 1.0) {
    std::cerr << "fixed scaffold outside probability is not one\n";
    return EXIT_FAILURE;
  }

  if (!check_case("GCGCGCGCGCGC", "..(....)....", true)) {
    return EXIT_FAILURE;
  }

  std::cout << "scfg2_outside_api=ok\n";
  return EXIT_SUCCESS;
}

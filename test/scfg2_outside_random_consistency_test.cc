#include "scfg2/api/beam.hh"
#include "scfg2/api/exact.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

bool close(double lhs, double rhs, double tolerance = 1e-9) {
  return std::fabs(lhs - rhs) <=
         tolerance * std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
}

std::string random_scaffold(std::mt19937 &rng, int n) {
  std::string structure(static_cast<std::size_t>(n), '.');
  std::uniform_real_distribution<double> probability(0.0, 1.0);
  std::function<void(int, int)> fill = [&](int left, int right) {
    if (right - left <= 3) return;
    if (probability(rng) < 0.38) {
      std::uniform_int_distribution<int> partner(left + 4, right);
      const int mate = partner(rng);
      structure[static_cast<std::size_t>(left - 1)] = '(';
      structure[static_cast<std::size_t>(mate - 1)] = ')';
      fill(left + 1, mate - 1);
      fill(mate + 1, right);
    } else {
      fill(left + 1, right);
    }
  };
  fill(1, n);
  return structure;
}

std::string compatible_sequence(std::mt19937 &rng,
                                const std::string &structure) {
  static constexpr std::pair<char, char> pairs[] = {
      {'A', 'U'}, {'U', 'A'}, {'C', 'G'},
      {'G', 'C'}, {'G', 'U'}, {'U', 'G'}};
  static constexpr char bases[] = {'A', 'C', 'G', 'U'};
  std::uniform_int_distribution<int> pair_choice(0, 5);
  std::uniform_int_distribution<int> base_choice(0, 3);
  std::vector<int> stack;
  std::string sequence(structure.size(), 'A');
  for (std::size_t pos = 0; pos < structure.size(); ++pos) {
    if (structure[pos] == '(') {
      stack.push_back(static_cast<int>(pos));
    } else if (structure[pos] == ')') {
      const int left = stack.back();
      stack.pop_back();
      const auto pair = pairs[static_cast<std::size_t>(pair_choice(rng))];
      sequence[static_cast<std::size_t>(left)] = pair.first;
      sequence[pos] = pair.second;
    } else {
      sequence[pos] = bases[static_cast<std::size_t>(base_choice(rng))];
    }
  }
  return sequence;
}

std::vector<std::pair<int, int>> scaffold_pairs(
    const std::string &scaffold) {
  std::vector<int> stack;
  std::vector<std::pair<int, int>> out;
  for (std::size_t pos = 0; pos < scaffold.size(); ++pos) {
    if (scaffold[pos] == '(') {
      stack.push_back(static_cast<int>(pos + 1));
    } else if (scaffold[pos] == ')') {
      out.emplace_back(stack.back(), static_cast<int>(pos + 1));
      stack.pop_back();
    }
  }
  return out;
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

struct Stats {
  int checked = 0;
  int nonempty_scaffolds = 0;
  int round_positive_cases = 0;
  int square_positive_cases = 0;
};

bool check_case(const std::string &sequence,
                const std::string &scaffold,
                Stats *stats) {
  cparty::scfg2::api::ExactSession exact(sequence, scaffold);
  const auto exact_inside = exact.inside();
  const auto exact_outside = exact.outside();
  cparty::scfg2::api::BeamSession beam(
      sequence, scaffold, cparty::scfg2::api::BeamConfig{160, 30, 30});
  const auto beam_outside = beam.outside();
  if (!exact_inside.ok || !exact_outside.ok || !exact_outside.reachable ||
      !beam_outside.ok || !beam_outside.reachable ||
      !close(std::log(exact_inside.partition_weight),
             exact_outside.log_partition_weight) ||
      !close(exact_outside.log_partition_weight,
             beam_outside.log_partition_weight) ||
      exact_outside.pair_probabilities.size() !=
          beam_outside.pair_probabilities.size()) {
    std::cerr << "outside random setup/root parity failed sequence="
              << sequence << " scaffold=" << scaffold << '\n';
    return false;
  }
  std::vector<double> paired(sequence.size() + 1, 0.0);
  bool round_positive = false;
  bool square_positive = false;
  for (std::size_t index = 0;
       index < exact_outside.pair_probabilities.size();
       ++index) {
    const auto &lhs = exact_outside.pair_probabilities[index];
    const auto &rhs = beam_outside.pair_probabilities[index];
    if (lhs.i != rhs.i || lhs.j != rhs.j ||
        !close(lhs.round_probability, rhs.round_probability) ||
        !close(lhs.square_probability, rhs.square_probability) ||
        !close(lhs.union_probability, rhs.union_probability)) {
      std::cerr << "outside random exact/beam marginal mismatch\n";
      return false;
    }
    paired[static_cast<std::size_t>(lhs.i)] += lhs.union_probability;
    paired[static_cast<std::size_t>(lhs.j)] += lhs.union_probability;
    round_positive = round_positive || lhs.round_probability > 1e-12;
    square_positive = square_positive || lhs.square_probability > 1e-12;
  }
  for (std::size_t i = 1; i <= sequence.size(); ++i) {
    if (!close(
            paired[i] + exact_outside.unpaired_probabilities[i], 1.0)) {
      std::cerr << "outside random residue normalization failed\n";
      return false;
    }
  }
  const auto fixed = scaffold_pairs(scaffold);
  for (const auto &[i, j] : fixed) {
    const auto *pair = find_pair(exact_outside.pair_probabilities, i, j);
    if (!pair || pair->round_probability != 1.0 ||
        pair->square_probability != 0.0 ||
        pair->union_probability != 1.0) {
      std::cerr << "outside random fixed-pair contract failed\n";
      return false;
    }
  }
  if (fixed.empty() && square_positive) {
    std::cerr << "empty scaffold produced positive Square marginal\n";
    return false;
  }
  ++stats->checked;
  stats->nonempty_scaffolds += fixed.empty() ? 0 : 1;
  stats->round_positive_cases += round_positive ? 1 : 0;
  stats->square_positive_cases += square_positive ? 1 : 0;
  return true;
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  std::mt19937 rng(20260728);
  Stats stats;
  for (const int length : {10, 12, 16, 20}) {
    for (int sample = 0; sample < 8; ++sample) {
      std::string scaffold =
          sample < 2 ? std::string(static_cast<std::size_t>(length), '.')
                     : random_scaffold(rng, length);
      const std::string sequence = compatible_sequence(rng, scaffold);
      if (!check_case(sequence, scaffold, &stats)) return EXIT_FAILURE;
    }
  }
  if (stats.nonempty_scaffolds == 0 ||
      stats.round_positive_cases == 0 ||
      stats.square_positive_cases == 0) {
    std::cerr << "outside random sweep lacked required semantic coverage\n";
    return EXIT_FAILURE;
  }
  std::cout << "scfg2_outside_random_consistency=ok"
            << " checked=" << stats.checked
            << " nonempty_scaffolds=" << stats.nonempty_scaffolds
            << " round_positive_cases=" << stats.round_positive_cases
            << " square_positive_cases=" << stats.square_positive_cases
            << '\n';
  return EXIT_SUCCESS;
}

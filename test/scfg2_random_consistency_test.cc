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

using Pair = std::pair<int, int>;

std::vector<Pair> pairs_for_lane(const std::string &structure,
                                 char open,
                                 char close) {
  std::vector<int> stack;
  std::vector<Pair> out;
  for (std::size_t pos = 0; pos < structure.size(); ++pos) {
    if (structure[pos] == open) {
      stack.push_back(static_cast<int>(pos + 1));
    } else if (structure[pos] == close) {
      if (stack.empty()) return {};
      out.emplace_back(stack.back(), static_cast<int>(pos + 1));
      stack.pop_back();
    }
  }
  return stack.empty() ? out : std::vector<Pair>{};
}

bool crosses(const Pair &lhs, const Pair &rhs) {
  return (lhs.first < rhs.first && rhs.first < lhs.second &&
          lhs.second < rhs.second) ||
         (rhs.first < lhs.first && lhs.first < rhs.second &&
          rhs.second < lhs.second);
}

bool generated_rendering_is_semantic(
    const std::string &scaffold,
    const cparty::scfg2::api::ViterbiResult &best) {
  const auto scaffold_round = pairs_for_lane(scaffold, '(', ')');
  const auto generated_round =
      pairs_for_lane(best.generated_structure, '(', ')');
  const auto generated_square =
      pairs_for_lane(best.generated_structure, '[', ']');

  if (generated_round.size() + generated_square.size() !=
      best.generated_pairs.size()) {
    return false;
  }
  for (const auto &pair : best.generated_pairs) {
    const Pair needle{pair.i, pair.j};
    const bool rendered =
        std::find(generated_round.begin(), generated_round.end(), needle) !=
            generated_round.end() ||
        std::find(generated_square.begin(), generated_square.end(), needle) !=
            generated_square.end();
    if (!rendered) return false;
  }
  for (const auto &square : generated_square) {
    bool has_crossing_round = false;
    for (const auto &round : scaffold_round) {
      if (crosses(square, round)) has_crossing_round = true;
    }
    for (const auto &round : generated_round) {
      if (crosses(square, round)) has_crossing_round = true;
    }
    if (!has_crossing_round) return false;
  }
  return true;
}

std::string compose_expected_union(const std::string &scaffold,
                                   const std::string &generated) {
  std::string out = scaffold;
  for (std::size_t pos = 0; pos < generated.size(); ++pos) {
    if (generated[pos] != '.') out[pos] = generated[pos];
  }
  return out;
}

struct SweepStats {
  int round_output_cases = 0;
  int square_output_cases = 0;
  int mixed_output_cases = 0;
};

std::string random_scaffold(std::mt19937 &rng, int n) {
  std::string structure(static_cast<std::size_t>(n), '.');
  std::uniform_real_distribution<double> probability(0.0, 1.0);
  std::function<void(int, int)> fill = [&](int left, int right) {
    if (right - left <= 3) return;
    if (probability(rng) < 0.32) {
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

bool check_case(const std::string &sequence,
                const std::string &scaffold,
                SweepStats *stats = nullptr) {
  cparty::scfg2::api::ExactSession session(sequence, scaffold);
  const auto inside = session.inside();
  const auto best = session.viterbi();
  if (!inside.ok || !best.ok || !(inside.partition_weight > 0.0) ||
      !(best.weight > 0.0) ||
      best.weight > inside.partition_weight * (1.0 + 1e-12)) {
    std::cerr << "inside_ok=" << inside.ok << " inside_weight="
              << inside.partition_weight << " inside_error=" << inside.error
              << " best_ok=" << best.ok << " best_weight=" << best.weight
              << " best_error=" << best.error << "\n";
    return false;
  }
  if (!generated_rendering_is_semantic(scaffold, best)) {
    std::cerr << "generated Round/Square rendering contract failed\n";
    return false;
  }
  if (stats) {
    const bool has_round =
        !pairs_for_lane(best.generated_structure, '(', ')').empty();
    const bool has_square =
        !pairs_for_lane(best.generated_structure, '[', ']').empty();
    if (has_round) ++stats->round_output_cases;
    if (has_square) ++stats->square_output_cases;
    if (has_round && has_square) ++stats->mixed_output_cases;
  }
  const auto evaluated = session.evaluate_energy(best.generated_structure);
  const auto probability = session.probability(best.generated_structure);
  const bool ok = evaluated.ok && evaluated.reachable && probability.ok &&
                  probability.reachable &&
         std::fabs(evaluated.energy - best.energy) <= 1e-10 &&
         probability.probability >= 0.0 &&
                  probability.probability <= 1.0 + 1e-12 &&
                  best.union_structure == compose_expected_union(
                      best.scaffold_structure, best.generated_structure) &&
                  evaluated.union_structure == best.union_structure &&
                  probability.union_structure == best.union_structure;
  if (!ok) {
    std::cerr << "best_energy=" << best.energy
              << " target=" << best.generated_structure
              << " generated_pairs=";
    for (const auto &pair : best.generated_pairs) {
      std::cerr << pair.i << '-' << pair.j << ',';
    }
    std::cerr
              << " eval_ok=" << evaluated.ok
              << " eval_reachable=" << evaluated.reachable
              << " eval_energy=" << evaluated.energy
              << " eval_error=" << evaluated.error
              << " probability_ok=" << probability.ok
              << " probability_reachable=" << probability.reachable
              << " probability=" << probability.probability
              << " probability_error=" << probability.error << "\n";
  }
  return ok;
}

bool has_pair(const std::vector<cparty::scfg2::api::BasePair> &pairs,
              int i,
              int j) {
  for (const auto &pair : pairs) {
    if (pair.i == i && pair.j == j) return true;
  }
  return false;
}

int sweep_count() {
  const char *raw = std::getenv("SCFG2_RANDOM_CONSISTENCY_COUNT");
  if (!raw || !*raw) return 12;
  const int parsed = std::atoi(raw);
  return parsed > 0 ? parsed : 12;
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);

  // Regression for the case that exposed the distinction between grammar
  // lane and fixed-vs-generated ownership. V(5,16) and V(6,12) are Round
  // grammar pairs, but they are not in scaffold G and therefore must appear
  // in the public generated pair map and its rendered structure.
  const std::string regression_sequence = "UGUGGCAAUUUGAAAUCUAA";
  const std::string regression_scaffold = "...(..(...).....)...";
  cparty::scfg2::api::ExactSession regression(
      regression_sequence, regression_scaffold);
  const auto regression_best = regression.viterbi();
  if (!regression_best.ok ||
      !has_pair(regression_best.generated_pairs, 5, 16) ||
      !has_pair(regression_best.generated_pairs, 6, 12) ||
      regression_best.generated_structure != "....((.....)...)...." ||
      regression_best.union_structure != "...((((...))...))...") {
    std::cerr << "generated Round-V regression was not exported\n";
    return EXIT_FAILURE;
  }
  if (!check_case(regression_sequence, regression_scaffold)) {
    std::cerr << "generated Round-V regression failed cross-mode replay\n";
    return EXIT_FAILURE;
  }

  std::mt19937 rng(20260717);
  int checked = 0;
  SweepStats stats;
  const int count = sweep_count();
  for (const int length : {20, 30, 50, 80, 100}) {
    for (int sample = 0; sample < count; ++sample) {
      const std::string scaffold = random_scaffold(rng, length);
      const std::string sequence = compatible_sequence(rng, scaffold);
      if (!check_case(sequence, scaffold, &stats)) {
        std::cerr << "random consistency failure: length=" << length
                  << " sample=" << sample << " sequence=" << sequence
                  << " scaffold=" << scaffold << "\n";
        return EXIT_FAILURE;
      }
      ++checked;
    }
  }
  std::cout << "scfg2_random_consistency=ok checked=" << checked
            << " round_output_cases=" << stats.round_output_cases
            << " square_output_cases=" << stats.square_output_cases
            << " mixed_output_cases=" << stats.mixed_output_cases << "\n";
  return EXIT_SUCCESS;
}

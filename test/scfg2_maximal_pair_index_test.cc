// MaximalPairIndex::geometry(i, j) must agree with the scan
// exact_direct_child_geometry(ctx, i, j) on every span: same "exactly one"
// verdict, and the same child and gaps when there is one.

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/deductions/maximal_pair_index.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

using namespace cparty::scfg2;
std::size_t exhaustive_spans = 0;

struct ScaffoldContext : deductions::PairabilityContext {
  std::vector<cand_pos_t> round, square;
  explicit ScaffoldContext(const std::string &g, bool mixed = false)
      : round(g.size() + 1, -1), square(g.size() + 1, -1) {
    std::vector<cand_pos_t> stack;
    for (cand_pos_t pos = 1; pos <= sequence_length(); ++pos) {
      if (g[pos - 1] == '(') stack.push_back(pos);
      if (g[pos - 1] == ')') {
        const cand_pos_t left = stack.back();
        stack.pop_back();
        auto &lane = mixed && left % 2 ? square : round;
        lane[left] = pos;
        lane[pos] = left;
      }
    }
  }
  cand_pos_t turn() const override { return 0; }
  cand_pos_t sequence_length() const override { return static_cast<cand_pos_t>(round.size() - 1); }
  bool can_pair(cand_pos_t, cand_pos_t) const override { return false; }
  cand_pos_t unpaired_prefix(cand_pos_t) const override { return 0; }
  cand_pos_t pair_round(cand_pos_t pos) const override { return round[pos]; }
  cand_pos_t pair_square(cand_pos_t pos) const override { return square[pos]; }
};

bool agrees(const deductions::ExactDirectChildGeometry &want,
            const deductions::ExactDirectChildGeometry &got) {
  return got.count == std::min(want.count, 2) &&
      want.child_i == got.child_i && want.child_j == got.child_j &&
      want.left_gap == got.left_gap && want.right_gap == got.right_gap;
}

void check_scaffold(const std::string &g) {
  for (bool mixed : {false, true}) {
    ScaffoldContext ctx(g, mixed);
    deductions::MaximalPairIndex index(ctx);
    if (!index.usable()) std::exit(EXIT_FAILURE);
    for (cand_pos_t i = 1; i < ctx.sequence_length(); ++i) {
      for (cand_pos_t j = i + 1; j <= ctx.sequence_length(); ++j) {
        if (!agrees(deductions::exact_direct_child_geometry(ctx, i, j), index.geometry(i, j))) {
          std::cerr << "exhaustive maximal-pair mismatch G=" << g
                    << " mixed=" << mixed << " i=" << i << " j=" << j << '\n';
          std::exit(EXIT_FAILURE);
        }
        ++exhaustive_spans;
      }
    }
  }
}

void enumerate(std::string &g, int balance, int remaining) {
  if (!remaining) {
    if (!balance) check_scaffold(g);
    return;
  }
  for (char ch : {'.', '(', ')'}) {
    const int next = balance + (ch == '(' ? 1 : ch == ')' ? -1 : 0);
    if (next < 0 || next > remaining - 1) continue;
    g.push_back(ch);
    enumerate(g, next, remaining - 1);
    g.pop_back();
  }
}

void check_long() {
  for (int n : {4096, 16384, 65536, 262144}) {
    for (int shape = 0; shape < 3; ++shape) {
      std::string g(n, '.');
      if (shape == 1) for (int x = 2; x < n / 2; ++x) {
        g[x] = '(';
        g[n - x - 1] = ')';
      }
      if (shape == 2) for (int x = 2; x + 3 < n; x += 3) {
        g[x] = '(';
        g[x + 2] = ')';
      }
      ScaffoldContext ctx(g);
      deductions::MaximalPairIndex index(ctx);
      if (!index.usable() || index.storage_bytes() > 512 * static_cast<std::size_t>(n + 1))
        std::exit(EXIT_FAILURE);
      const auto got = index.geometry(1, n);
      if (got.count != shape || (shape == 1 && (got.child_i != 3 || got.child_j != n - 2))) {
        std::cerr << "long maximal-pair mismatch n=" << n << " shape=" << shape << '\n';
        std::exit(EXIT_FAILURE);
      }
      // Short spans sample every part of the long input, using the original
      // scan as oracle while avoiding quadratic test work in n.
      for (int i = 1; i + 12 <= n; i += 997) {
        if (!agrees(deductions::exact_direct_child_geometry(ctx, i, i + 12), index.geometry(i, i + 12)))
          std::exit(EXIT_FAILURE);
      }
      std::cout << "n=" << n << " shape=" << shape << " storage_bytes=" << index.storage_bytes() << '\n';
    }
  }
  ScaffoldContext crossing("......");
  crossing.round[2] = 4; crossing.round[4] = 2;
  crossing.square[3] = 5; crossing.square[5] = 3;
  if (deductions::MaximalPairIndex(crossing).usable() ||
      deductions::exact_direct_child_geometry(crossing, 1, 6).count != 2)
    std::exit(EXIT_FAILURE);
  ScaffoldContext dual("....");
  dual.round[1] = dual.square[1] = 4;
  dual.round[4] = dual.square[4] = 1;
  if (deductions::MaximalPairIndex(dual).usable()) std::exit(EXIT_FAILURE);
  ScaffoldContext asymmetric("....");
  asymmetric.round[1] = 4;
  if (deductions::MaximalPairIndex(asymmetric).usable()) std::exit(EXIT_FAILURE);
}

char complement(char c) {
  switch (c) {
    case 'A': return 'U';
    case 'U': return 'A';
    case 'G': return 'C';
    default: return 'G';
  }
}

// Random pseudoknot-free scaffold: recursive stems with random helix lengths,
// branching, and unpaired gaps.
void fill(std::mt19937 &rng, std::string &g, std::string &seq, int lo, int hi,
          int depth) {
  std::uniform_int_distribution<int> coin(0, 99);
  int pos = lo;
  while (pos + 6 <= hi) {
    if (coin(rng) < 35) {
      ++pos;
      continue;
    }
    const int max_span = hi - pos + 1;
    const int span = std::uniform_int_distribution<int>(6, max_span)(rng);
    const int a = pos, b = pos + span - 1;
    const int helix = std::uniform_int_distribution<int>(1, std::max(1, (span - 4) / 2))(rng);
    int t = 0;
    for (; t < helix && b - a - 2 * t >= 4; ++t) {
      g[a + t] = '(';
      g[b - t] = ')';
      seq[b - t] = complement(seq[a + t]);
    }
    if (depth < 4) fill(rng, g, seq, a + t, b - t, depth + 1);
    pos = b + 1;
  }
}

}  // namespace

int main(int argc, char **argv) {
  using namespace cparty::scfg2;
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  const int cases = argc > 1 ? std::atoi(argv[1]) : 40;
  std::mt19937 rng(20261008);
  std::size_t spans = 0, singles = 0;
  for (int n = 1; n <= 10; ++n) {
    std::string g;
    enumerate(g, 0, n);
  }
  check_long();
  for (int c = 0; c < cases; ++c) {
    const int n = std::uniform_int_distribution<int>(10, 120)(rng);
    std::string seq(n, 'A');
    for (auto &ch : seq) ch = "ACGU"[std::uniform_int_distribution<int>(0, 3)(rng)];
    std::string g(n, '.');
    fill(rng, g, seq, 0, n - 1, 0);
    cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
    replay::WFinalExactInsideContext ctx(owner);
    const deductions::MaximalPairIndex index(ctx);
    if (!index.usable()) {
      std::cerr << "index unusable on a round-only scaffold\n";
      return EXIT_FAILURE;
    }
    for (cand_pos_t i = 1; i <= n; ++i) {
      for (cand_pos_t j = i + 1; j <= n; ++j) {
        const auto want = deductions::exact_direct_child_geometry(ctx, i, j);
        const auto got = index.geometry(i, j);
        ++spans;
        const bool want_single = want.count == 1;
        const bool got_single = got.count == 1;
        if (want_single != got_single ||
            (want.count == 0) != (got.count == 0) ||
            (want_single && (want.child_i != got.child_i || want.child_j != got.child_j ||
                             want.left_gap != got.left_gap ||
                             want.right_gap != got.right_gap))) {
          std::cerr << "mismatch n=" << n << " G=" << g << " i=" << i << " j=" << j
                    << " want count=" << want.count << " (" << want.child_i << ","
                    << want.child_j << ") got count=" << got.count << " ("
                    << got.child_i << "," << got.child_j << ")\n";
          return EXIT_FAILURE;
        }
        if (want_single) ++singles;
      }
    }
  }
  std::cout << "spans=" << spans << " single_child_spans=" << singles << "\n";
  std::cout << "exhaustive_spans=" << exhaustive_spans << '\n';
  std::cout << "scfg2_maximal_pair_index=ok\n";
  return EXIT_SUCCESS;
}

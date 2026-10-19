// A finite beam must never lose a scaffold pair.  Every structure contains
// all pairs of G, so the beam keeps one scaffold-only derivation (G' empty)
// outside the width.  For random scaffolds, including hairpin loops longer
// than the beam hairpin cap and wide multiloops, and for very small widths,
// check that the root is reachable and that every scaffold pair has a
// nonzero V item.

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <utility>

namespace {

char complement(char c) {
  switch (c) {
    case 'A': return 'U';
    case 'U': return 'A';
    case 'G': return 'C';
    default: return 'G';
  }
}

void stem(std::string &g, std::string &seq, int a, int b, int len) {
  for (int t = 0; t < len && b - a - 2 * t >= 4; ++t) {
    g[a + t] = '(';
    g[b - t] = ')';
    seq[b - t] = complement(seq[a + t]);
  }
}

// An outer stem closing a multiloop of `branches` hairpins; one hairpin has
// a loop of `long_loop` bases (beyond the beam hairpin cap of 30).
std::string scaffold(std::mt19937 &rng, std::string &seq, int branches, int long_loop) {
  const int n = static_cast<int>(seq.size());
  std::string g(n, '.');
  stem(g, seq, 0, n - 1, 3);
  int pos = 5;
  for (int k = 0; k < branches; ++k) {
    const int gap = std::uniform_int_distribution<int>(2, 9)(rng);
    const int loop = k == 0 ? long_loop : std::uniform_int_distribution<int>(3, 12)(rng);
    const int len = std::uniform_int_distribution<int>(2, 5)(rng);
    const int a = pos + gap, b = a + 2 * len + loop - 1;
    if (b >= n - 6) break;
    stem(g, seq, a, b, len);
    pos = b + 1;
  }
  return g;
}

}  // namespace

int main(int argc, char **argv) {
  using namespace cparty::scfg2;
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  const int cases = argc > 1 ? std::atoi(argv[1]) : 12;
  // --no-force: run without the scaffold-only derivation, to show the test
  // detects the failure it guards against.
  const bool force = !(argc > 2 && std::string(argv[2]) == "--no-force");
  if (force) {
    // V stays KeepAll while WM is finite. The scaffold forcing probe uses
    // a huge finite width; reserving that many prefix slots used to throw.
    std::string seq = "GGGAAACCC", g = "(((...)))";
    cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
    replay::WFinalExactInsideContext ctx(owner);
    auto mixed = replay::WFinalVpWmbpOnlineWidths{};
    mixed.wm = 20;
    const auto reference = replay::run_w_final_vp_wmbp_online<semirings::LogSumExpSemiring>(
        ctx, seq.size(), replay::WFinalVpWmbpOnlineWidths{});
    const auto actual = replay::run_w_final_vp_wmbp_online<semirings::LogSumExpSemiring>(
        ctx, seq.size(), mixed);
    const deductions::ItemKey root{grammar::NonTerminal::W, 1,
        static_cast<cand_pos_t>(seq.size()), -1, -1};
    if (!std::isfinite(actual.chart.get(root)) ||
        std::fabs(actual.chart.get(root) - reference.chart.get(root)) > 1e-12) {
      std::cerr << "mixed KeepAll/finite scaffold forcing changed the saturated root\n";
      return EXIT_FAILURE;
    }
  }
  std::mt19937 rng(20261009);
  int checked = 0;
  double total_seconds = 0.0;
  int failures = 0;
  std::size_t maximum_forced_per_family_column = 0;
  for (int c = 0; c < cases; ++c) {
    const int n = std::uniform_int_distribution<int>(120, 320)(rng);
    std::string seq(n, 'A');
    for (auto &ch : seq) ch = "ACGU"[std::uniform_int_distribution<int>(0, 3)(rng)];
    std::string g = scaffold(rng, seq, std::uniform_int_distribution<int>(2, 6)(rng),
                             std::uniform_int_distribution<int>(31, 45)(rng));
    cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
    replay::WFinalExactInsideContext ctx(owner);
    if (force) {
      // The standard forced set is a single G'=empty traceback. Its pure
      // pk-free branches partition the sequence, so no family occurs twice
      // at the same right boundary. This is stronger than an arbitrary
      // caller-supplied forced set and bounds each forced family's width by b+1.
      const auto forced_items = replay::scaffold_only_derivation_items(
          ctx, n, replay::WFinalVpWmbpOnlineWidths::uniform(20));
      std::map<std::pair<grammar::NonTerminal, cand_pos_t>, std::size_t> counts;
      for (const auto &item : *forced_items) {
        using NT = grammar::NonTerminal;
        if (item.nonterminal != NT::W && item.nonterminal != NT::V &&
            item.nonterminal != NT::VM && item.nonterminal != NT::WM &&
            item.nonterminal != NT::WMv) {
          std::cerr << "scaffold-only traceback contains a PK family\n";
          return EXIT_FAILURE;
        }
        const auto count = ++counts[{item.nonterminal, item.j}];
        maximum_forced_per_family_column = std::max(maximum_forced_per_family_column, count);
        if (count > 1) {
          std::cerr << "multiple scaffold-only forced items in one family/column G=" << g << '\n';
          return EXIT_FAILURE;
        }
      }
    }
    for (const std::size_t b : {1, 2, 5, 20}) {
      auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(b);
      widths.force_scaffold_derivation = force;
      const auto started = std::chrono::steady_clock::now();
      const auto run = replay::run_w_final_vp_wmbp_online<semirings::LogSumExpSemiring>(
          ctx, n, widths);
      total_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - started).count();
      const auto value = [&](grammar::NonTerminal nt, cand_pos_t i, cand_pos_t j) {
        return run.chart.get(deductions::ItemKey{nt, i, j, -1, -1});
      };
      if (!std::isfinite(value(grammar::NonTerminal::W, 1, n))) {
        ++failures;
        if (!force) continue;
        std::cerr << "root unreachable b=" << b << " G=" << g << "\n";
        return EXIT_FAILURE;
      }
      for (cand_pos_t i = 1; i <= n; ++i) {
        const cand_pos_t j = ctx.pair_partner(i);
        if (j > i && !std::isfinite(value(grammar::NonTerminal::V, i, j))) {
          ++failures;
          if (!force) break;
          std::cerr << "scaffold pair (" << i << "," << j << ") lost at b=" << b
                    << " G=" << g << "\n";
          return EXIT_FAILURE;
        }
      }
      ++checked;
    }
  }
  std::cout << "checked_runs=" << checked << " failures=" << failures
            << " beam_seconds=" << total_seconds
            << " max_forced_per_family_column=" << maximum_forced_per_family_column
            << "\nscfg2_beam_scaffold_forcing=ok\n";
  return EXIT_SUCCESS;
}

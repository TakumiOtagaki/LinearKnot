// Diagnostic: where does a finite beam disconnect the root?
//
// Prints the first column j at which the exterior prefix W(1,j) becomes zero
// and, for the scaffold pairs around that column, whether their V item
// survived.  Usage: SEQUENCE SCAFFOLD WIDTH [PARAM_FILE]

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

int main(int argc, char **argv) {
  using namespace cparty::scfg2;
  if (argc < 4) {
    std::cerr << "usage: " << argv[0] << " SEQUENCE SCAFFOLD WIDTH [PARAM_FILE]\n";
    return EXIT_FAILURE;
  }
  std::string sequence = argv[1];
  std::string scaffold = argv[2];
  const std::size_t width = std::stoull(argv[3]);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argc > 4 ? argv[4] : "params/rna_Turner04.par", 1);
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto run = replay::run_w_final_vp_wmbp_online<semirings::LogSumExpSemiring>(
      ctx, n, replay::WFinalVpWmbpOnlineWidths::uniform(width));
  const auto value = [&](grammar::NonTerminal nt, cand_pos_t i, cand_pos_t j) {
    return run.chart.get(deductions::ItemKey{nt, i, j, -1, -1});
  };
  {
    const auto forced = replay::scaffold_only_derivation_items(
        ctx, n, replay::WFinalVpWmbpOnlineWidths::uniform(width));
    std::size_t forced_v = 0;
    for (const auto &item : *forced) forced_v += item.nonterminal == grammar::NonTerminal::V;
    std::cout << "scaffold_only_derivation items=" << forced->size()
              << " V items=" << forced_v << "\n";
  }
  {
    const auto forced = replay::scaffold_only_derivation_items(
        ctx, n, replay::WFinalVpWmbpOnlineWidths::uniform(width));
    std::vector<deductions::ItemKey> missing;
    for (const auto &item : *forced) {
      if (!std::isfinite(run.chart.get(item))) missing.push_back(item);
    }
    std::sort(missing.begin(), missing.end(), [](const auto &a, const auto &b) {
      return a.j - a.i < b.j - b.i;
    });
    std::cout << "forced items with zero value in this run: " << missing.size() << "\n";
    for (std::size_t k = 0; k < std::min<std::size_t>(missing.size(), 12); ++k) {
      const auto &m = missing[k];
      std::cout << "  " << grammar::nonterminal_name(m.nonterminal) << "(" << m.i << ","
                << m.j << ")\n";
    }
  }
  std::cout << "n=" << n << " width=" << width
            << " log W(1,n)=" << value(grammar::NonTerminal::W, 1, n) << "\n";
  // First weakly closed prefix [1, j] whose W(1, j) is zero.
  for (cand_pos_t j = 1; j <= n; ++j) {
    if (ctx.weakly_closed(1, j) &&
        !std::isfinite(value(grammar::NonTerminal::W, 1, j))) {
      std::cout << "first_zero_weakly_closed_prefix=" << j << "\n";
      break;
    }
  }
  // Scaffold pairs whose V item is zero, innermost (shortest) first.  A pair
  // crossed by G' pairs may legitimately live only inside a band (BE).
  std::vector<std::pair<cand_pos_t, cand_pos_t>> zero_pairs;
  for (cand_pos_t i = 1; i <= n; ++i) {
    const cand_pos_t j = ctx.pair_partner(i);
    if (j > i && !std::isfinite(value(grammar::NonTerminal::V, i, j))) {
      zero_pairs.push_back({i, j});
    }
  }
  std::sort(zero_pairs.begin(), zero_pairs.end(), [](const auto &a, const auto &b) {
    return a.second - a.first < b.second - b.first;
  });
  std::cout << "scaffold pairs with zero V: " << zero_pairs.size() << "\n";
  for (std::size_t k = 0; k < std::min<std::size_t>(zero_pairs.size(), 8); ++k) {
    const auto [i, j] = zero_pairs[k];
    cand_pos_t inner_pairs = 0;
    for (cand_pos_t x = i + 1; x < j; ++x) inner_pairs += ctx.pair_partner(x) > x;
    std::cout << "  (" << i << "," << j << ") span=" << j - i + 1
              << " inner_scaffold_pairs=" << inner_pairs << " G=" << scaffold.substr(i - 1, std::min<cand_pos_t>(j - i + 1, 80))
              << "\n";
  }
  return EXIT_SUCCESS;
}

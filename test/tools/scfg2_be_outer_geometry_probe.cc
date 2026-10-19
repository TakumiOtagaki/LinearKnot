// Observation only: compare a single scaffold-child candidate with the shared
// exact BE generator. No inside, outer, weights or pruning policy is changed.
#include <cassert>
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/geometry/static_scaffold_geometry.hh"

#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace sc = cparty::scfg2;
using Rule = sc::grammar::RuleId;

// Geometry queries match WFinal's static-scaffold boundary; this probe has
// no sequence-dependent energy or legacy-oracle interpretation.
struct Context final : sc::deductions::PairabilityContext {
  sc::geometry::StaticScaffoldGeometry tree;
  explicit Context(const std::string &g) : tree(g, g.size()) {}
  cand_pos_t turn() const override { return 3; }
  cand_pos_t sequence_length() const override { return tree.size(); }
  bool can_pair(cand_pos_t, cand_pos_t) const override { return false; }
  cand_pos_t pair_partner(cand_pos_t p) const override { return tree.partner(p); }
  cand_pos_t unpaired_prefix(cand_pos_t p) const override {
    return tree.unpaired_prefix(p);
  }
  bool is_unpaired(cand_pos_t p) const override {
    return p >= 1 && p <= tree.size() && tree.partner(p) < 0;
  }
  bool is_empty_region(cand_pos_t i, cand_pos_t j) const override {
    return j > i && (j == i + 1 || unpaired_prefix(j - 1) >= j - i - 1);
  }
  bool weakly_closed(cand_pos_t i, cand_pos_t j) const override {
    return tree.weakly_closed(i, j);
  }
};

struct Counts {
  std::uint64_t cases = 0, parents = 0, full_positions = 0;
  std::uint64_t direct_probes = 0, deductions = 0;
};

void check(const std::string &g, Counts &counts) {
  Context ctx(g);
  ++counts.cases;
  constexpr std::array rules{Rule::BE_INTERNAL_LOOP, Rule::BE_WIP_WIP,
                             Rule::BE_WIP_BASEPAIR, Rule::BE_BASEPAIR_WIP};
  for (cand_pos_t i = 1; i <= ctx.tree.size(); ++i) {
    const auto j = ctx.tree.partner(i);
    if (j <= i) continue;
    for (cand_pos_t ip = i; ip < j; ++ip) {
      const auto jp = ctx.tree.partner(ip);
      if (jp <= ip || jp > j) continue;
      ++counts.parents;
      const sc::deductions::ItemKey parent{sc::grammar::NonTerminal::BE,
                                           i, j, ip, jp};
      const auto k = ip == i ? -1 : ctx.tree.child_below(i, ip);
      const auto l = ctx.tree.partner(k);
      for (const auto rule : rules) {
        counts.full_positions += ip - i;
        const auto full = sc::deductions::generate_exact_basic_deductions(
            parent, rule, ctx);
        bool direct = false;
        if (k > i) {
          ++counts.direct_probes;
          sc::grammar::RuleSplit split{};
          split.k = k;
          split.l = l;
          direct = sc::deductions::exact_be_split_allows(parent, rule, k, l, ctx)
              && sc::deductions::is_valid(
                  sc::deductions::expand_basic_rule(parent, rule, split));
        }
        if (full.size() != static_cast<std::size_t>(direct) ||
            (!full.empty() && (full[0].split.k != k || full[0].split.l != l)))
          throw std::runtime_error("BE scaffold-child mismatch: " + g);
        counts.deductions += full.size();
      }
    }
  }
}

void random_region(std::string &g, int lo, int hi, std::mt19937 &rng) {
  if (lo >= hi) return;
  if (hi - lo < 2 || rng() % 4 == 0) {
    if (rng() % 5 == 0) g[lo] = 'x';
    random_region(g, lo + 1, hi, rng);
    return;
  }
  const int right = lo + 1 + rng() % (hi - lo - 1);
  g[lo] = '(';
  g[right] = ')';
  random_region(g, lo + 1, right, rng);
  random_region(g, right + 1, hi, rng);
}

int main() {
  try {
    Counts counts;
    std::vector<std::vector<std::string>> balanced(11);
    balanced[0].push_back("");
    for (std::size_t n = 1; n < balanced.size(); ++n) {
      for (const auto &rest : balanced[n - 1]) balanced[n].push_back("." + rest);
      for (std::size_t inner = 0; inner + 2 <= n; ++inner)
        for (const auto &a : balanced[inner])
          for (const auto &b : balanced[n - inner - 2])
            balanced[n].push_back("(" + a + ")" + b);
      for (const auto &g : balanced[n]) check(g, counts);
    }
    std::mt19937 rng(20261013);
    for (int sample = 0; sample < 200; ++sample) {
      std::string g(20 + rng() % 141, '.');
      random_region(g, 0, g.size(), rng);
      check(g, counts);
    }
    for (int depth : {8, 32, 64, 128})
      check(std::string(depth, '(') + "...." + std::string(depth, ')'), counts);
    std::cout << "cases\tparents\tfull_split_positions\tdirect_probes\tdeductions\n"
              << counts.cases << '\t' << counts.parents << '\t'
              << counts.full_positions << '\t' << counts.direct_probes << '\t'
              << counts.deductions << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

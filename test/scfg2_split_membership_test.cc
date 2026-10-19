// The single-split membership tests used by the beam event builders must
// agree exactly with the generators they replace:
//   exact_v_internal_split_allows        vs V_INTERNAL
//   exact_vp_internal_loop_split_allows  vs VP_INTERNAL_LOOP
//   exact_be_split_allows                vs BE_INTERNAL_LOOP / BE_WIP_WIP /
//                                           BE_WIP_BASEPAIR / BE_BASEPAIR_WIP
// For every parent and every candidate split, the predicate is true iff the
// generator emits that split.

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/engine/basic_schedule.hh"

#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using namespace cparty::scfg2;
using SplitSet = std::set<std::pair<cand_pos_t, cand_pos_t>>;

char complement(char c) {
  switch (c) {
    case 'A': return 'U';
    case 'U': return 'A';
    case 'G': return 'C';
    default: return 'G';
  }
}

void fill(std::mt19937 &rng, std::string &g, std::string &seq, int lo, int hi,
          int depth) {
  std::uniform_int_distribution<int> coin(0, 99);
  int pos = lo;
  while (pos + 6 <= hi) {
    if (coin(rng) < 40) {
      ++pos;
      continue;
    }
    const int span = std::uniform_int_distribution<int>(6, hi - pos + 1)(rng);
    const int a = pos, b = pos + span - 1;
    const int helix =
        std::uniform_int_distribution<int>(1, std::max(1, (span - 4) / 2))(rng);
    int t = 0;
    for (; t < helix && b - a - 2 * t >= 4; ++t) {
      // occasional bulge inside a helix
      if (t > 0 && coin(rng) < 15 && b - a - 2 * t >= 6) {
        ++t;
        continue;
      }
      g[a + t] = '(';
      g[b - t] = ')';
      seq[b - t] = complement(seq[a + t]);
    }
    if (depth < 3) fill(rng, g, seq, a + t, b - t, depth + 1);
    pos = b + 1;
  }
}

SplitSet generated(const deductions::ItemKey &parent, grammar::RuleId rule,
                   const deductions::PairabilityContext &ctx) {
  SplitSet out;
  for (const auto &d : deductions::generate_exact_basic_deductions(parent, rule, ctx)) {
    out.insert({d.split.k, d.split.l});
  }
  return out;
}

bool report(const std::string &what, const std::string &g,
            const deductions::ItemKey &parent, cand_pos_t k, cand_pos_t l,
            bool want) {
  std::cerr << what << " mismatch G=" << g << " parent=(" << parent.i << ","
            << parent.j << "," << parent.ip << "," << parent.jp << ") split=("
            << k << "," << l << ") generator=" << want << "\n";
  return false;
}

bool check_selected(const deductions::ItemKey &parent,
                    const deductions::PairabilityContext &ctx, int n) {
  const auto full = deductions::generate_exact_basic_deductions_for_item(parent, ctx);
  using Id = std::tuple<grammar::RuleId, cand_pos_t, cand_pos_t, cand_pos_t, cand_pos_t>;
  const auto id = [](const auto &d) { return Id{d.rule, d.split.k, d.split.l, d.split.p, d.split.q}; };
  std::vector<std::vector<cand_pos_t>> selections{{}, {-1, 1, n / 2, n, n + 1}, {}};
  for (int k = 1; k <= n; ++k) if (k % 3 == 1) selections[2].push_back(k);
  for (auto &positions : selections) {
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
    std::vector<Id> want, got;
    for (const auto &d : full) {
      const auto coordinate = deductions::exact_linear_split_coordinate(d);
      if (!coordinate || std::binary_search(positions.begin(), positions.end(), *coordinate))
        want.push_back(id(d));
    }
    for (const auto &d : deductions::generate_exact_basic_deductions_for_item(parent, ctx, &positions))
      got.push_back(id(d));
    if (want != got) {
      std::cerr << "selected-split mismatch family=" << grammar::nonterminal_name(parent.nonterminal)
                << " parent=" << parent.i << ',' << parent.j << '\n';
      return false;
    }
  }
  return true;
}

}  // namespace

int main(int argc, char **argv) {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  const int cases = argc > 1 ? std::atoi(argv[1]) : 25;
  std::mt19937 rng(20261009);
  std::size_t checked = 0, positives = 0;
  for (int c = 0; c < cases; ++c) {
    const int n = std::uniform_int_distribution<int>(12, 60)(rng);
    std::string seq(n, 'A');
    for (auto &ch : seq) ch = "ACGU"[std::uniform_int_distribution<int>(0, 3)(rng)];
    std::string g(n, '.');
    fill(rng, g, seq, 0, n - 1, 0);
    cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
    replay::WFinalExactInsideContext ctx(owner);
    const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
    constexpr auto family_count = engine::basic_exact_span_family_order().size();
    std::array<int, family_count> phase{};
    int ordinal = 0;
    for (auto family : engine::basic_exact_span_family_order()) {
      if (engine::legacy_span_phase(family) != ordinal) return EXIT_FAILURE;
      phase.at(static_cast<std::size_t>(family)) = ordinal++;
    }
    int previous = -1;
    for (const auto &item : schedule) {
      const int bucket = (item.j - item.i) * family_count + phase.at(static_cast<std::size_t>(item.nonterminal));
      if (bucket < previous) return EXIT_FAILURE;
      previous = bucket;
      // Sample coordinates, but cover every family and each anchored BE.
      if (item.nonterminal == grammar::NonTerminal::BE ||
          (item.i * 7 + item.j * 11 + c) % 29 == 0)
        if (!check_selected(item, ctx, n)) return EXIT_FAILURE;
    }

    for (cand_pos_t i = 1; i <= n; ++i) {
      for (cand_pos_t j = i + 1; j <= n; ++j) {
        for (const auto nt : {grammar::NonTerminal::V, grammar::NonTerminal::VP}) {
          const deductions::ItemKey parent{nt, i, j, -1, -1};
          const auto rule = nt == grammar::NonTerminal::V
                                ? grammar::RuleId::V_INTERNAL
                                : grammar::RuleId::VP_INTERNAL_LOOP;
          const auto want = generated(parent, rule, ctx);
          for (cand_pos_t k = i + 1; k < j; ++k) {
            for (cand_pos_t l = k + 1; l < j; ++l) {
              const bool got =
                  nt == grammar::NonTerminal::V
                      ? deductions::exact_v_internal_split_allows(parent, k, l, ctx)
                      : deductions::exact_vp_internal_loop_split_allows(parent, k, l, ctx);
              const bool expected = want.count({k, l}) > 0;
              ++checked;
              positives += expected;
              if (got != expected) {
                report(nt == grammar::NonTerminal::V ? "V_INTERNAL" : "VP_INTERNAL_LOOP",
                       g, parent, k, l, expected);
                return EXIT_FAILURE;
              }
            }
          }
        }
      }
    }
    // BE parents: outer G pair (i, j) and an inner G pair (ip, jp) of its band.
    for (cand_pos_t i = 1; i <= n; ++i) {
      const cand_pos_t j = ctx.pair_partner(i);
      if (j <= i) continue;
      for (cand_pos_t ip = i; ip < j; ++ip) {
        const cand_pos_t jp = ctx.pair_partner(ip);
        if (jp <= ip || jp > j) continue;
        const deductions::ItemKey parent{grammar::NonTerminal::BE, i, j, ip, jp};
        for (const auto rule :
             {grammar::RuleId::BE_INTERNAL_LOOP, grammar::RuleId::BE_WIP_WIP,
              grammar::RuleId::BE_WIP_BASEPAIR, grammar::RuleId::BE_BASEPAIR_WIP}) {
          const auto want = generated(parent, rule, ctx);
          for (cand_pos_t k = i + 1; k < j; ++k) {
            const cand_pos_t l = ctx.pair_partner(k);
            if (l <= k) continue;
            const bool got = deductions::exact_be_split_allows(parent, rule, k, l, ctx);
            const bool expected = want.count({k, l}) > 0;
            ++checked;
            positives += expected;
            if (got != expected) {
              report("BE", g, parent, k, l, expected);
              return EXIT_FAILURE;
            }
          }
        }
      }
    }
  }
  std::cout << "checked=" << checked << " emitted=" << positives << "\n";
  std::cout << "scfg2_split_membership=ok\n";
  return EXIT_SUCCESS;
}

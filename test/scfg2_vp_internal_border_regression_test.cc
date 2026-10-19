#include "scfg2/api/exact.hh"
#include "scfg2/deductions/generate_exact_basic.hh"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

struct RegressionCase {
  const char *id;
  const char *sequence;
  const char *scaffold;
  const char *target;
};

// An observation context lets the support contract be checked independently
// of energy tables, the scaffold index, and the optimized interval builder.
struct LoopContext : cparty::scfg2::deductions::PairabilityContext {
  int Bp = -1, b = -1, B = -1, bp = -1;
  std::array<bool, 65> blocked{};
  cand_pos_t turn() const override { return 3; }
  cand_pos_t sequence_length() const override { return 64; }
  cand_pos_t pair_round(cand_pos_t) const override {
    return -2;
  }
  bool is_unpaired(cand_pos_t pos) const override {
    return pos > 0 && pos <= 64 && !blocked[pos];
  }
  bool can_pair(cand_pos_t i, cand_pos_t j) const override {
    return j - i > turn() && is_unpaired(i) && is_unpaired(j);
  }
  cand_pos_t unpaired_prefix(cand_pos_t pos) const override {
    return pos;
  }
  cand_pos_t border_Bp(cand_pos_t, cand_pos_t) const override { return Bp; }
  cand_pos_t border_b(cand_pos_t, cand_pos_t) const override { return b; }
  cand_pos_t border_B(cand_pos_t, cand_pos_t) const override { return B; }
  cand_pos_t border_bp(cand_pos_t, cand_pos_t) const override { return bp; }
};

bool check_loop_domain() {
  using namespace cparty::scfg2;
  const deductions::ItemKey parent{grammar::NonTerminal::VP, 2, 48, -1, -1};
  std::size_t cases = 0;
  for (const int Bp : {-2, -1, 0, 8, 16})
    for (const int b : {-2, 0, 12})
      for (const int B : {-2, -1, 0, 16, 30})
        for (const int bp : {-1, 0, 24})
          for (const int blocked : {0, 5, 45}) {
            LoopContext ctx;
            ctx.Bp = Bp; ctx.b = b; ctx.B = B; ctx.bp = bp;
            if (blocked) ctx.blocked[blocked] = true;
            std::vector<std::pair<int, int>> expected, actual;
            // Enumerate the entire rectangle, then apply inequalities. No
            // clipped bounds or edge-coordinate formula from production.
            for (int k = parent.i + 1; k < parent.j; ++k) {
              for (int l = parent.j - 1; l > k; --l) {
                const int gaps = (k - parent.i - 1) + (parent.j - l - 1);
                if (gaps == 0 || gaps >= MAXLOOP) continue;
                // The current context contract clips this left edge only
                // when both border values are positive; retain its sentinels.
                if (Bp > 0 && b > 0 && (k >= Bp || k >= b)) continue;
                if (l <= B || l <= bp || !ctx.can_pair(k, l) ||
                    !ctx.is_empty_region(parent.i, k) ||
                    !ctx.is_empty_region(l, parent.j)) continue;
                expected.emplace_back(k, l);
              }
            }
            const auto ds = deductions::generate_exact_basic_deductions(
                parent, grammar::RuleId::VP_INTERNAL_LOOP, ctx);
            for (const auto &d : ds) {
              if (d.child_count != 1 || d.children[0].item.i != d.split.k ||
                  d.children[0].item.j != d.split.l) return false;
              actual.emplace_back(d.split.k, d.split.l);
            }
            if (actual != expected) {
              std::cerr << "VP loop domain/order mismatch: Bp=" << Bp
                        << " b=" << b << " B=" << B << " bp=" << bp
                        << " blocked=" << blocked << '\n';
              return false;
            }
            ++cases;
          }
  std::cout << "vp_loop_domain_cases=" << cases << '\n';
  return true;
}

bool check_case(const RegressionCase &test_case) {
  cparty::scfg2::api::ExactSession session(
      test_case.sequence, test_case.scaffold);

  const auto evaluated = session.evaluate_energy(test_case.target);
  if (!evaluated.ok || !evaluated.reachable ||
      !std::isfinite(evaluated.energy)) {
    std::cerr << test_case.id
              << " target energy is not reachable: ok=" << evaluated.ok
              << " reachable=" << evaluated.reachable
              << " energy=" << evaluated.energy
              << " error=" << evaluated.error << '\n';
    return false;
  }

  const auto best = session.viterbi();
  if (!best.ok || !std::isfinite(best.energy) ||
      best.energy > evaluated.energy + 1e-10) {
    std::cerr << test_case.id
              << " Viterbi/target ordering is invalid: ok=" << best.ok
              << " viterbi_energy=" << best.energy
              << " target_energy=" << evaluated.energy
              << " error=" << best.error << '\n';
    return false;
  }

  const auto probability = session.probability(test_case.target);
  if (!probability.ok || !probability.reachable ||
      !std::isfinite(probability.log_probability) ||
      probability.probability < 0.0 ||
      probability.probability > 1.0 + 1e-12) {
    std::cerr << test_case.id
              << " target probability is invalid: ok=" << probability.ok
              << " reachable=" << probability.reachable
              << " probability=" << probability.probability
              << " log_probability=" << probability.log_probability
              << " error=" << probability.error << '\n';
    return false;
  }

  return true;
}

}  // namespace

int main() {
  if (!check_loop_domain()) return EXIT_FAILURE;
  setenv(
      "CPARTY_FIXED_ENERGY_PARAM_FILE",
      "params/rna_DirksPierce09.par",
      1);

  const RegressionCase cases[] = {
      {
          "PKB00048",
          "AAAACACAUAACUUCCGACAACCGCCGGGACGGAGGAAAAAAAAAAAAAAAAAACCGGCGG",
          ".....................(((((((..........................)))))))",
          "...........[[[[[[.[.........].]]]]]].........................",
      },
      {
          "PKB436",
          "CGUCGACUACCGCAGGGCGGAUGGUAGAGGUAGUGGAAAAGACCGCUGUCCGACGAAAAAAGAACGCCCA",
          "..............(((((.............................................))))).",
          "[[[[[.[[[[[[[......].]]]]]].[[[[[[[[......]]]]]]]]]]]]]...............",
      },
  };

  for (const auto &test_case : cases) {
    if (!check_case(test_case)) return EXIT_FAILURE;
  }

  std::cout << "scfg2_vp_internal_border_regression=ok\n";
  return EXIT_SUCCESS;
}

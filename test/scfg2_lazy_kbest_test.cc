#include "scfg2/replay/lazy_kbest.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace {
using namespace cparty::scfg2;
using S = semirings::FactorizedLogSumExpSemiring;
void require(bool ok, const char *message) { if (!ok) throw std::runtime_error(message); }

// Independent exhaustive Cartesian expansion, used only on tiny forests.
std::vector<double> exhaustive(const semirings::FactorizedInsideGraph &g, std::uint32_t node) {
  std::vector<double> all;
  const auto &state = g.states()[node];
  for (std::size_t t = state.first_term; t < std::size_t(state.first_term) + state.term_count; ++t) {
    const auto &term = g.terms()[t];
    std::vector<double> partial{term.local_log_weight};
    for (std::size_t f = 0; f < term.factor_count; ++f) {
      const auto child = exhaustive(g, g.factors()[term.first_factor + f]);
      std::vector<double> joined;
      for (const auto a : partial) for (const auto b : child) joined.push_back(a + b);
      partial = std::move(joined);
    }
    all.insert(all.end(), partial.begin(), partial.end());
  }
  std::sort(all.begin(), all.end(), std::greater<double>());
  return all;
}

void check_cartesian() {
  semirings::FactorizedInsideGraph g;
  S::Scope scope(g);
  const auto x = S::plus(S::from_log_weight(0.9), S::plus(
      S::from_log_weight(-0.4), S::from_log_weight(-1.8))).checkpointed();
  const auto y = S::plus(S::times(x, S::from_log_weight(2.0)),
                         S::from_log_weight(-3.0)).checkpointed();
  // Three factors, including two independent occurrences of the same node;
  // shared child ranks generate many duplicate candidate paths in the grid.
  const auto product = S::times(S::times(x, y), x).checkpointed();
  const auto root = S::plus(S::times(product, S::from_log_weight(-0.6)),
                            S::times(y, S::from_log_weight(0.3))).checkpointed();
  const auto want = exhaustive(g, root.node);
  replay::LazyKBest lazy(g);
  for (std::size_t k = 0; k < want.size(); ++k) {
    require(lazy.ensure(root.node, k), "lazy list ended early");
    require(std::abs(lazy.derivation(root.node, k).log_weight - want[k]) < 1e-12,
            "lazy rank differs from exhaustive sort");
  }
  require(!lazy.ensure(root.node, want.size()), "lazy list repeated a candidate");
}

void check_local_precision_and_depth() {
  semirings::FactorizedInsideGraph g;
  S::Scope scope(g);
  const auto huge = S::from_log_weight(1e16).checkpointed();
  const auto shifted = S::times(huge, S::from_log_weight(0.25)).checkpointed();
  require(g.terms()[g.states()[shifted.node].first_term].local_log_weight == 0.25,
          "local weight was reconstructed by cancellation");
  auto root = S::plus(S::from_log_weight(0.0), S::from_log_weight(-1.0)).checkpointed();
  for (int i = 0; i < 30000; ++i) root = S::times(root, S::one()).checkpointed();
  replay::LazyKBest lazy(g);
  require(lazy.ensure(root.node, 1), "unary chain lost its second derivation");
  require(lazy.derivation(root.node, 1).log_weight == -1.0, "wrong unary chain score");
  std::size_t visited = 0;
  lazy.visit_derivation(root.node, 1, [&](auto) { ++visited; });
  require(visited >= 30000, "traceback did not traverse the chain");
  require(!lazy.ensure(root.node, 2), "unary chain did not exhaust");
}
}

int main() {
  try { check_cartesian(); check_local_precision_and_depth(); }
  catch (const std::exception &e) { std::cerr << e.what() << '\n'; return EXIT_FAILURE; }
  std::cout << "scfg2_lazy_kbest=ok\n";
}

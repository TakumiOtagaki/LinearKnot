#include "scfg2/engine/chart.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {
using namespace cparty::scfg2;
using S = semirings::FactorizedLogSumExpSemiring;

bool close(double a, double b) { return std::fabs(a - b) < 1e-10; }

bool check_graph() {
  semirings::FactorizedInsideGraph graph;
  S::Scope scope(graph);
  const auto x = S::from_weight(2.0).checkpointed();
  const auto y = S::from_weight(3.0).checkpointed();
  const auto product = S::times(x, y).checkpointed();
  const auto root = S::plus(product, S::times(product, x)).checkpointed();
  const auto posterior = graph.posteriors(root.node);
  // Z = xy + x^2 y = 18. The shared product is visited once on both
  // alternatives; x is visited once or twice, including the repeated factor.
  if (!close(root.log_value, std::log(18.0)) ||
      !close(posterior[root.node], 1.0) ||
      !close(posterior[product.node], 1.0) ||
      !close(posterior[x.node], 5.0 / 3.0) ||
      !close(posterior[y.node], 1.0)) return false;
  const auto beta = graph.outside_log_beta(root.node);
  for (auto node : {x.node, y.node, product.node, root.node}) {
    if (!close(std::exp(graph.states()[node].log_value + beta[node] - root.log_value),
               posterior[node])) return false;
    constexpr double eps = 1e-5;
    const double derivative = (graph.evaluate_with_log_shift(root.node, node, eps) -
        graph.evaluate_with_log_shift(root.node, node, -eps)) / (2.0 * eps);
    if (std::fabs(derivative - posterior[node]) > 1e-8) return false;
  }
  // Frozen prefixes survive extending a sum, and are counted as a shared
  // state if they have another reader. A zero contribution records no edge.
  auto total = S::plus(x, y);
  const auto prefix = total;
  total = S::plus(total, x);
  const auto both = S::plus(prefix, total).checkpointed();
  if (!close(both.log_value, std::log(12.0))) return false;
  const auto p = graph.posteriors(both.node);
  if (!close(p[x.node], 0.5) || !close(p[y.node], 0.5)) return false;
  const auto count = graph.states().size();
  if (!S::is_zero(S::times(x, S::zero())) ||
      S::plus(S::zero(), both).node != both.node ||
      graph.states().size() != count) return false;
  return true;
}

bool check_prune_capture() {
  semirings::FactorizedInsideGraph graph;
  S::Scope scope(graph);
  engine::Chart<S::value_type> chart(S::zero());
  const deductions::ItemKey wip{grammar::NonTerminal::WIP, 2, 8, -1, -1};
  const auto seed = S::from_weight(3.0).checkpointed();
  chart.set(wip, S::times(seed, S::from_weight(2.0)));
  const auto pre = chart.get(wip);
  const auto tail = S::times(pre, S::from_weight(5.0)).checkpointed();
  chart.set(wip, S::zero());
  if (chart.contains(wip) || !S::is_zero(chart.get(wip))) return false;
  const auto posterior = graph.posteriors(tail.node);
  return close(posterior[pre.node], 1.0) && close(posterior[seed.node], 1.0);
}

bool check_nested_sum_after_shared_product() {
  semirings::FactorizedInsideGraph graph;
  S::Scope scope(graph);
  const auto x = S::from_weight(2.0).checkpointed();
  const auto y = S::from_weight(3.0).checkpointed();
  const auto product = S::times(x, y);
  const auto nested = S::plus(product, x);
  const auto root = S::plus(product, nested).checkpointed();
  const auto nested_node = nested.checkpointed().node;
  // Z = xy + (xy + x). The unsealed product has two readers; sealing the
  // nested sum must finish before either of the root's terms is recorded.
  if (!close(root.log_value, std::log(14.0)) ||
      graph.states().size() != 4 || graph.terms().size() != 6 ||
      graph.factors().size() != 6) return false;
  const auto posterior = graph.posteriors(root.node);
  if (!close(posterior[x.node], 1.0) ||
      !close(posterior[y.node], 6.0 / 7.0)) return false;
  const auto &state = graph.states()[root.node];
  if (state.term_count != 2) return false;
  const auto &first = graph.terms()[state.first_term];
  const auto &second = graph.terms()[state.first_term + 1];
  return first.factor_count == 2 && second.factor_count == 1 &&
      graph.factors()[first.first_factor] == x.node &&
      graph.factors()[first.first_factor + 1] == y.node &&
      graph.factors()[second.first_factor] == nested_node;
}
}  // namespace

int main() {
  if (!check_graph() || !check_prune_capture() ||
      !check_nested_sum_after_shared_product()) {
    std::cerr << "factorized analytic outside failed\n";
    return EXIT_FAILURE;
  }
  std::cout << "scfg2_factorized_outside=ok\n";
}

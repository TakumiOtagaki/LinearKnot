#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/engine/execute.hh"
#include "scfg2/engine/outside.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

namespace {
using namespace cparty::scfg2;
using deductions::Deduction;
using deductions::ItemKey;
using grammar::NonTerminal;
using grammar::RuleId;

bool close(double a, double b) {
  return a == b || std::fabs(a - b) <=
      1e-11 * std::max({1.0, std::fabs(a), std::fabs(b)});
}

// Actual shared BE geometry/local weights, with independent leaf values. The
// frozen graph is only a test oracle; the reverse deduction stores no edges.
bool check_be_sources() {
  std::string sequence = "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  std::string scaffold = "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  engine::AllLocalWeightModel<decltype(ctx)> weights(&ctx);
  const auto schedule = engine::build_pkfree_basic_exact_schedule(
      static_cast<cand_pos_t>(sequence.size()), ctx);
  std::set<RuleId> checked;
  bool diagonal_stack = false;
  for (const auto &parent : schedule) {
    if (parent.nonterminal != NonTerminal::BE) continue;
    for (const auto &ded : deductions::generate_exact_basic_deductions_for_item(
             parent, ctx)) {
      if (!engine::legacy_span_schedule_supports(ded)) continue;
      const double local = weights.weight_for(ded);
      if (!(local > 0.0)) continue;
      if (ded.rule == RuleId::BE_STACK &&
          parent.i == parent.ip && parent.j == parent.jp) {
        const auto &child = ded.children[0].item;
        if (child.ip != parent.i + 1 || child.jp != parent.j - 1) return false;
        diagonal_stack = true;
      }
      if (!checked.insert(ded.rule).second) continue;
      using Plain = semirings::LogSumExpSemiring;
      using Frozen = semirings::FactorizedLogSumExpSemiring;
      engine::Chart<double> alpha(Plain::zero());
      engine::Chart<double> beta(Plain::zero());
      engine::Chart<Frozen::value_type> frozen(Frozen::zero());
      semirings::FactorizedInsideGraph graph;
      Frozen::value_type root;
      {
        Frozen::Scope scope(graph);
        for (std::size_t c = 0; c < ded.child_count; ++c) {
          const double value = 1.25 + 0.75 * c;
          alpha.set(ded.children[c].item, std::log(value));
          frozen.set(ded.children[c].item, Frozen::from_weight(value));
        }
        root = semirings::checkpoint_value(engine::evaluate_deduction<Frozen>(
            ded, frozen, [&](const Deduction &d) { return weights.weight_for(d); }));
      }
      const double parent_beta = std::log(2.3);
      engine::propagate_outside_deduction<Plain>(
          ded, alpha, parent_beta, &beta,
          [&](const Deduction &d) { return weights.weight_for(d); });
      const auto expected = graph.outside_log_beta(root.node);
      for (std::size_t c = 0; c < ded.child_count; ++c) {
        const auto &child = ded.children[c].item;
        if (!close(beta.get(child), expected[frozen.get(child).node] + parent_beta)) {
          std::cerr << "BE reverse mismatch rule=" << static_cast<int>(ded.rule)
                    << " child=" << c << '\n';
          return false;
        }
      }
      if (ded.child_count == 0 && !beta.values().empty()) return false;
    }
  }
  std::cout << "BE regenerated rules=" << checked.size()
            << " diagonal_stack=" << diagonal_stack << '\n';
  return checked.size() == 6 && diagonal_stack;
}

bool check_zero_and_repeated_occurrences() {
  using Plain = semirings::SumProductSemiring;
  const ItemKey x{NonTerminal::WIP, 2, 4, -1, -1};
  const ItemKey y{NonTerminal::BE, 5, 12, 7, 10};
  Deduction ded;
  ded.parent = {NonTerminal::BE, 1, 15, 7, 10};
  ded.rule = RuleId::BE_WIP_WIP;
  ded.child_count = 3;
  ded.children[0].item = x;
  ded.children[1].item = y;
  // An engine-level multiplicity witness, not a valid BE geometry witness.
  ded.children[2].item = x;
  for (const auto values : {std::array<double, 2>{3, 5},
                            std::array<double, 2>{0, 5},
                            std::array<double, 2>{3, 0},
                            std::array<double, 2>{0, 0}}) {
    engine::Chart<double> alpha(0), beta(0);
    alpha.set(x, values[0]);
    alpha.set(y, values[1]);
    engine::propagate_outside_deduction<Plain>(
        ded, alpha, 7.0, &beta, [](const Deduction &) { return 2.0; });
    // d(14*x*x*y)/dx and d(14*x*x*y)/dy, including the zero boundary.
    if (!close(beta.get(x), 28 * values[0] * values[1]) ||
        !close(beta.get(y), 14 * values[0] * values[0])) return false;
    engine::Chart<double> log_alpha(semirings::LogSumExpSemiring::zero());
    engine::Chart<double> log_beta(semirings::LogSumExpSemiring::zero());
    log_alpha.set(x, std::log(values[0]));
    log_alpha.set(y, std::log(values[1]));
    engine::propagate_outside_deduction<semirings::LogSumExpSemiring>(
        ded, log_alpha, std::log(7.0), &log_beta,
        [](const Deduction &) { return 2.0; });
    if (!close(std::exp(log_beta.get(x)), beta.get(x)) ||
        !close(std::exp(log_beta.get(y)), beta.get(y))) return false;
  }
  return true;
}

// An artificial engine DAG isolates retention semantics. The raw B=2*A branch
// is pruned before consumers run. Both passes use one filtered provider, so no
// old B value or branch contributes to Z or to outside.
bool check_retained_provider() {
  using Plain = semirings::SumProductSemiring;
  using Frozen = semirings::FactorizedLogSumExpSemiring;
  const ItemKey a{NonTerminal::W, 1, 1, -1, -1};
  const ItemKey b{NonTerminal::W, 1, 2, -1, -1};
  const ItemKey c{NonTerminal::W, 1, 3, -1, -1};
  const ItemKey root{NonTerminal::W, 1, 4, -1, -1};
  const std::vector<ItemKey> retained{a, c, root};
  const auto live = [&](const ItemKey &key) {
    return std::any_of(retained.begin(), retained.end(), [&](const auto &item) {
      return deductions::ItemKeyEq{}(key, item);
    });
  };
  const auto branch = [](ItemKey parent, int weight, std::vector<ItemKey> children) {
    Deduction ded;
    ded.parent = parent;
    ded.split.k = weight;  // Test-local generic weight, not a grammar weight.
    ded.child_count = children.size();
    for (std::size_t i = 0; i < children.size(); ++i)
      ded.children[i].item = children[i];
    return ded;
  };
  const std::vector<Deduction> raw{
      branch(a, 2, {}), branch(b, 2, {a}), branch(c, 3, {a}),
      branch(root, 7, {b}), branch(root, 1, {c})};
  const auto generate = [&](const ItemKey &parent) {
    std::vector<Deduction> out;
    for (const auto &ded : raw) {
      if (!deductions::ItemKeyEq{}(ded.parent, parent) || !live(parent)) continue;
      bool admitted = true;
      for (std::size_t i = 0; i < ded.child_count; ++i)
        admitted = admitted && live(ded.children[i].item);
      if (admitted) out.push_back(ded);
    }
    return out;
  };
  const auto weight = [](const Deduction &ded) { return double(ded.split.k); };
  engine::Chart<double> alpha(0);
  engine::Chart<Frozen::value_type> frozen(Frozen::zero());
  semirings::FactorizedInsideGraph graph;
  {
    Frozen::Scope scope(graph);
    for (const auto &parent : retained) {
      engine::run_cell_update<Plain>(parent, generate(parent), &alpha, weight);
      engine::run_cell_update<Frozen>(parent, generate(parent), &frozen, weight);
    }
  }
  std::size_t observed = 0;
  double observed_leaf_beta = 0.0;
  const auto beta = engine::run_regenerated_outside<Plain>(
      retained, alpha, root, generate, weight, engine::DeductionChildItem{},
      [&](const Deduction &ded, double parent_beta) {
        ++observed;
        if (ded.child_count == 0) observed_leaf_beta = parent_beta;
      });
  const auto reference = graph.outside_log_beta(frozen.get(root).node);
  for (const auto &parent : retained)
    if (!close(beta.get(parent), std::exp(reference[frozen.get(parent).node])))
      return false;
  return alpha.get(root) == 6.0 && beta.get(a) == 3.0 && beta.get(b) == 0.0 &&
         generate(b).empty() && generate(root).size() == 1 &&
         observed == 3 && observed_leaf_beta == 3.0;
}
}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  if (!check_be_sources() || !check_zero_and_repeated_occurrences() ||
      !check_retained_provider())
    return EXIT_FAILURE;
  std::cout << "scfg2_regenerated_outside=ok\n";
  return EXIT_SUCCESS;
}

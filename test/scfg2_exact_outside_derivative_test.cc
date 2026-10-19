#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/engine/chart.hh"
#include "scfg2/engine/outside.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;
using cparty::scfg2::grammar::NonTerminal;
using cparty::scfg2::grammar::RuleId;
using cparty::scfg2::semirings::LogSumExpSemiring;

ItemKey effective_child(const Deduction &deduction, std::size_t child_index) {
  const auto child = deduction.children[child_index].item;
  if (child.nonterminal == NonTerminal::VP &&
      deduction.rule == RuleId::WMBP_DIRECT_VP) {
    return ItemKey{NonTerminal::VP, child.i, child.j, -1, -1};
  }
  return child;
}

template <class Context>
cparty::scfg2::engine::Chart<double> direct_outside(
    const std::vector<ItemKey> &schedule,
    Context &ctx,
    const cparty::scfg2::engine::Chart<double> &inside,
    const ItemKey &root) {
  cparty::scfg2::engine::AllLocalWeightModel<Context> weights(&ctx);
  return cparty::scfg2::engine::run_regenerated_outside<LogSumExpSemiring>(
      schedule, inside, root,
      [&](const ItemKey &parent) {
        auto candidates = cparty::scfg2::deductions::
            generate_exact_basic_deductions_for_item(parent, ctx);
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
            [](const Deduction &deduction) {
              return !cparty::scfg2::engine::legacy_span_schedule_supports(
                  deduction);
            }), candidates.end());
        return candidates;
      },
      [&](const Deduction &deduction) { return weights.weight_for(deduction); },
      effective_child);
}

bool close(double lhs, double rhs, double tolerance) {
  return std::fabs(lhs - rhs) <=
         tolerance * std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
}

bool check_case(const std::string &sequence, const std::string &scaffold) {
  std::string mutable_sequence = sequence;
  std::string mutable_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      mutable_sequence, mutable_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto schedule =
      cparty::scfg2::engine::build_pkfree_basic_exact_schedule(n, ctx);
  const ItemKey root{NonTerminal::W, 1, n, -1, -1};

  const auto ordinary =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<
          LogSumExpSemiring>(
          schedule,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          });

  cparty::scfg2::semirings::LogSumExpTape tape;
  cparty::scfg2::engine::Chart<
      cparty::scfg2::semirings::LogSumExpTapeValue>
      taped(cparty::scfg2::semirings::TapeLogSumExpSemiring::zero());
  {
    cparty::scfg2::semirings::TapeLogSumExpSemiring::Scope scope(tape);
    taped = cparty::scfg2::replay::run_w_final_exact_dp_schedule<
        cparty::scfg2::semirings::TapeLogSumExpSemiring>(
        schedule,
        ctx,
        [](const Deduction &) { return true; },
        [](const ItemKey &, const Deduction &, auto, const auto &, double) {});
  }

  const auto taped_root = taped.get(root);
  if (!close(ordinary.get(root), taped_root.log_value, 1e-12)) {
    std::cerr << "root mismatch sequence=" << sequence
              << " ordinary=" << ordinary.get(root)
              << " taped=" << taped_root.log_value << '\n';
    return false;
  }
  const auto adjoints = tape.reverse(taped_root.node);
  const auto outside = direct_outside(schedule, ctx, ordinary, root);
  std::size_t compared = 0;
  for (const auto &[item, taped_value] : taped.values()) {
    if (taped_value.node ==
        cparty::scfg2::semirings::invalid_tape_node) {
      continue;
    }
    const double item_inside = ordinary.get(item);
    const double item_outside = outside.get(item);
    const double expected =
        LogSumExpSemiring::is_zero(item_inside) ||
                LogSumExpSemiring::is_zero(item_outside)
            ? 0.0
            : std::exp(item_inside + item_outside - ordinary.get(root));
    const double observed = adjoints[taped_value.node];
    if (!close(expected, observed, 1e-9)) {
      std::cerr << "item adjoint mismatch sequence=" << sequence
                << " nt="
                << cparty::scfg2::grammar::nonterminal_name(item.nonterminal)
                << " item=(" << item.i << ',' << item.j << ',' << item.ip
                << ',' << item.jp << ") expected=" << expected
                << " observed=" << observed << '\n';
      return false;
    }
    ++compared;
  }
  if (compared == 0) {
    std::cerr << "no taped items compared\n";
    return false;
  }
  std::cout << "case\tn=" << n << "\tcompared=" << compared
            << "\tlog_z=" << ordinary.get(root) << '\n';
  return true;
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  if (!check_case("GAGAAACAC", ".........")) return EXIT_FAILURE;
  if (!check_case("GCGCGCGCGCGC", "..(....)....")) return EXIT_FAILURE;
  if (!check_case("GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
                  "....(...)((....)..)((.(.((...))))(.(....)((...))))"))
    return EXIT_FAILURE;
  std::cout << "scfg2_exact_outside_derivative=ok\n";
  return EXIT_SUCCESS;
}

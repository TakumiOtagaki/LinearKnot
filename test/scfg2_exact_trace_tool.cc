#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/exact.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <array>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using TraceNode = cparty::scfg2::engine::DeductionTraceNode<double>;

std::vector<cparty::scfg2::engine::TracebackPair> parse_target(
    const std::string &structure) {
  std::array<std::vector<int>, 2> stacks;
  std::vector<cparty::scfg2::engine::TracebackPair> out;
  for (std::size_t pos = 0; pos < structure.size(); ++pos) {
    const char symbol = structure[pos];
    for (std::size_t lane = 0; lane < stacks.size(); ++lane) {
      const char open = lane == 0 ? '(' : '[';
      const char close = lane == 0 ? ')' : ']';
      if (symbol == open) {
        stacks[lane].push_back(static_cast<int>(pos + 1));
      } else if (symbol == close) {
        if (stacks[lane].empty()) throw std::runtime_error("unbalanced target");
        out.push_back({stacks[lane].back(), static_cast<cand_pos_t>(pos + 1),
                       cparty::scfg2::engine::TracebackPairLane::Square});
        stacks[lane].pop_back();
      }
    }
  }
  for (const auto &stack : stacks) {
    if (!stack.empty()) throw std::runtime_error("unbalanced target");
  }
  return out;
}

void print_trace(const TraceNode &node,
                 const cparty::scfg2::replay::WFinalExactInsideContext &ctx,
                 int depth = 0) {
  const auto &item = node.item;
  std::cout << std::string(static_cast<std::size_t>(depth * 2), ' ')
            << cparty::scfg2::grammar::nonterminal_name(item.nonterminal)
            << '(' << item.i << ',' << item.j;
  if (item.ip >= 0 || item.jp >= 0) std::cout << ';' << item.ip << ',' << item.jp;
  std::cout << ") value=" << node.chart_value;
  if (node.chosen_deduction.has_value()) {
    const auto &deduction = *node.chosen_deduction;
    std::cout << " rule=" << cparty::scfg2::grammar::rule_name(deduction.rule)
              << " split=" << deduction.split.k << ',' << deduction.split.l
              << ',' << deduction.split.p << ',' << deduction.split.q
              << " contribution=" << node.chosen_contribution;
    const auto output =
        cparty::scfg2::engine::deduction_output_pairs(deduction, ctx);
    std::cout << " output=";
    for (const auto &pair : output) {
      std::cout << pair.i << '-' << pair.j
                << (pair.lane == cparty::scfg2::engine::TracebackPairLane::Round
                        ? 'R'
                        : 'S')
                << ',';
    }
  }
  std::cout << '\n';
  for (const auto &child : node.children) print_trace(child, ctx, depth + 1);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "usage: scfg2_exact_trace_tool SEQUENCE SCAFFOLD [TARGET] [PARAM]\n";
    return 2;
  }
  if (argc == 5) setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[4], 1);
  std::string sequence = argv[1];
  std::string scaffold = argv[2];
  cparty::scfg2::api::ExactSession session(sequence, scaffold);
  if (argc >= 4) {
    const auto summary = session.evaluate_energy(argv[3]);
    if (!summary.ok || !summary.reachable) {
      std::cerr << "constrained evaluation failed: " << summary.error << '\n';
      return 1;
    }
    cparty::scfg2::weights::RuntimeSnapshot pf(
      sequence, scaffold);
    cparty::scfg2::replay::WFinalExactInsideContext ctx(pf);
    const auto result =
        cparty::scfg2::replay::run_w_final_exact_output_constrained_viterbi_traceback(
            pf, parse_target(argv[3]));
    std::cout << "mode=constrained energy=" << summary.energy
              << " G=" << summary.scaffold_structure
              << " G'=" << summary.generated_structure
              << " R=" << summary.union_structure
              << '\n';
    print_trace(result.trace, ctx);
  } else {
    const auto summary = session.viterbi();
    if (!summary.ok) {
      std::cerr << "Viterbi failed: " << summary.error << '\n';
      return 1;
    }
    cparty::scfg2::weights::RuntimeSnapshot pf(
      sequence, scaffold);
    cparty::scfg2::replay::WFinalExactInsideContext ctx(pf);
    const auto result =
        cparty::scfg2::replay::run_w_final_exact_viterbi_traceback(pf);
    std::cout << "mode=viterbi energy=" << summary.energy
              << " G=" << summary.scaffold_structure
              << " G'=" << summary.generated_structure
              << " R=" << summary.union_structure
              << " grammar_square=" << result.lane_structures.square
              << '\n';
    print_trace(result.trace, ctx);
  }
  return 0;
}

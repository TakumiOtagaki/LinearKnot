// Diagnostic for beam outside memory: how much of the logsumexp tape and the
// chart is reachable from the root / nonzero.  Usage: SEQUENCE SCAFFOLD WIDTH

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  using namespace cparty::scfg2;
  if (argc < 4) {
    std::cerr << "usage: " << argv[0] << " SEQUENCE SCAFFOLD WIDTH\n";
    return EXIT_FAILURE;
  }
  std::string sequence = argv[1];
  std::string scaffold = argv[2];
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  semirings::LogSumExpTape tape;
  const auto run = [&] {
    semirings::TapeLogSumExpSemiring::Scope scope(tape);
    return replay::run_w_final_vp_wmbp_online<semirings::TapeLogSumExpSemiring>(
        ctx, n, replay::WFinalVpWmbpOnlineWidths::uniform(std::stoull(argv[3])));
  }();
  const auto &nodes = tape.nodes();
  const auto root = run.chart.get(deductions::ItemKey{grammar::NonTerminal::W, 1, n, -1, -1});
  std::vector<char> reachable(nodes.size(), 0);
  std::vector<std::size_t> stack{root.node};
  std::size_t count = 0, constants = 0;
  while (!stack.empty()) {
    const auto x = stack.back();
    stack.pop_back();
    if (x >= nodes.size() || reachable[x]) continue;
    reachable[x] = 1;
    ++count;
    if (nodes[x].operation == semirings::LogSumExpTapeNode::Operation::Constant) {
      ++constants;
      continue;
    }
    stack.push_back(nodes[x].lhs);
    stack.push_back(nodes[x].rhs);
  }
  std::size_t by_op[3] = {0, 0, 0}, reach_by_op[3] = {0, 0, 0};
  std::size_t product_with_constant = 0, reach_product_with_constant = 0;
  for (std::size_t x = 0; x < nodes.size(); ++x) {
    const auto op = static_cast<int>(nodes[x].operation);
    ++by_op[op];
    if (reachable[x]) ++reach_by_op[op];
    if (nodes[x].operation == semirings::LogSumExpTapeNode::Operation::Product &&
        (nodes[nodes[x].lhs].operation == semirings::LogSumExpTapeNode::Operation::Constant ||
         nodes[nodes[x].rhs].operation == semirings::LogSumExpTapeNode::Operation::Constant)) {
      ++product_with_constant;
      if (reachable[x]) ++reach_product_with_constant;
    }
  }
  std::cout << "constant=" << by_op[0] << " (reachable " << reach_by_op[0] << ") product="
            << by_op[1] << " (reachable " << reach_by_op[1] << ", with a constant side "
            << product_with_constant << "/" << reach_product_with_constant << ") logadd="
            << by_op[2] << " (reachable " << reach_by_op[2] << ")\n";
  std::size_t chart_total = 0, chart_nonzero = 0;
  for (const auto &[item, value] : run.chart.values()) {
    ++chart_total;
    chart_nonzero += std::isfinite(value.log_value);
  }
  std::cout << "tape_nodes=" << nodes.size() << " reachable=" << count << " ("
            << 100.0 * count / nodes.size() << "%) reachable_constants=" << constants
            << "\nchart_entries=" << chart_total << " nonzero=" << chart_nonzero << " ("
            << 100.0 * chart_nonzero / chart_total << "%)\n";
  return EXIT_SUCCESS;
}

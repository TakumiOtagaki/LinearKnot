#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::size_t> parse_widths(const std::string &text) {
  std::vector<std::size_t> out;
  std::stringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (!token.empty()) out.push_back(std::stoull(token));
  }
  return out;
}

bool no_greater_than(double lhs, double rhs) {
  if (!std::isfinite(lhs)) return true;
  if (!std::isfinite(rhs)) return false;
  return lhs <= rhs + 1e-12;
}

bool same_log_value(double lhs, double rhs) {
  if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
    return !std::isfinite(lhs) && !std::isfinite(rhs) &&
           std::signbit(lhs) == std::signbit(rhs);
  }
  return std::fabs(lhs - rhs) <= 1e-12;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4 || argc > 5) {
    std::cerr << "usage: " << argv[0]
              << " SEQUENCE SCAFFOLD WIDTHS_CSV [PARAM_FILE]\n";
    return EXIT_FAILURE;
  }
  const std::string sequence = argv[1];
  const std::string scaffold = argv[2];
  const auto widths = parse_widths(argv[3]);
  if (sequence.size() != scaffold.size() || widths.empty() ||
      widths.front() == 0 || !std::is_sorted(widths.begin(), widths.end()) ||
      std::adjacent_find(widths.begin(), widths.end()) != widths.end()) {
    std::cerr << "invalid sequence/scaffold length or widths\n";
    return EXIT_FAILURE;
  }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 5 ? argv[4] : "params/rna_Turner04.par",
         1);

  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const cparty::scfg2::deductions::ItemKey root{
      cparty::scfg2::grammar::NonTerminal::W, 1, n, -1, -1};

  std::cout << std::setprecision(17)
            << "n\twidth\tlog_inside_root\tinside_finite"
               "\tinside_monotone\tlog_viterbi_root\tviterbi_finite"
               "\tviterbi_monotone\tviterbi_le_inside"
               "\tcoupled_log_inside_root\tcoupled_log_viterbi_root"
               "\tcoupled_viterbi_le_inside\tcoupled_inside_matches"
               "\telapsed_ms\tinside_item_visits\tviterbi_item_visits\n";
  double previous_inside = -std::numeric_limits<double>::infinity();
  double previous_viterbi = -std::numeric_limits<double>::infinity();
  for (const std::size_t width : widths) {
    const auto selected_widths =
        cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(width);
    const auto started = std::chrono::steady_clock::now();
    const auto inside_run =
        cparty::scfg2::replay::run_w_final_vp_wmbp_online<
            cparty::scfg2::semirings::LogSumExpSemiring>(
            ctx, n, selected_widths);
    const auto viterbi_run = cparty::scfg2::replay::run_w_final_vp_wmbp_online<
        cparty::scfg2::semirings::LogMaxSemiring>(
        ctx, n, selected_widths);
    const auto coupled_run =
        cparty::scfg2::replay::run_w_final_vp_wmbp_online<
            cparty::scfg2::semirings::CoupledLogSumMaxSemiring>(
            ctx, n, selected_widths);
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - started)
            .count();
    const double inside = inside_run.chart.get(root);
    const double viterbi = viterbi_run.chart.get(root);
    const auto coupled = coupled_run.chart.get(root);
    std::size_t inside_item_visits = 0;
    for (const std::size_t visits : inside_run.item_visits)
      inside_item_visits += visits;
    std::size_t viterbi_item_visits = 0;
    for (const std::size_t visits : viterbi_run.item_visits)
      viterbi_item_visits += visits;
    std::cout << n << '\t' << width << '\t' << inside << '\t'
              << (std::isfinite(inside) ? 1 : 0) << '\t'
              << (no_greater_than(previous_inside, inside) ? 1 : 0) << '\t'
              << viterbi << '\t' << (std::isfinite(viterbi) ? 1 : 0) << '\t'
              << (no_greater_than(previous_viterbi, viterbi) ? 1 : 0) << '\t'
              << (no_greater_than(viterbi, inside) ? 1 : 0) << '\t'
              << coupled.log_sum << '\t' << coupled.log_max << '\t'
              << (no_greater_than(coupled.log_max, coupled.log_sum) ? 1 : 0)
              << '\t'
              << (same_log_value(coupled.log_sum, inside) ? 1 : 0)
              << '\t'
              << elapsed_ms << '\t' << inside_item_visits << '\t'
              << viterbi_item_visits << '\n';
    previous_inside = inside;
    previous_viterbi = viterbi;
  }
  return EXIT_SUCCESS;
}

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

namespace {

using cparty::scfg2::deductions::ItemKey;
using cparty::scfg2::grammar::NonTerminal;
using cparty::scfg2::replay::WFinalOnlineNonterminalCounts;
using cparty::scfg2::replay::WFinalVpWmbpOnlineRun;

bool close(double lhs, double rhs, double tolerance = 1e-11) {
  if (lhs == rhs) return true;
  return std::fabs(lhs - rhs) <=
         tolerance * std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
}

bool same_counts(const WFinalOnlineNonterminalCounts &lhs,
                 const WFinalOnlineNonterminalCounts &rhs) {
  return lhs == rhs;
}

template <class TapedRun>
bool same_selection(
    const WFinalVpWmbpOnlineRun<double> &ordinary,
    const TapedRun &taped) {
  return ordinary.vp_candidates == taped.vp_candidates &&
         ordinary.vp_selected == taped.vp_selected &&
         ordinary.wmbp_candidates == taped.wmbp_candidates &&
         ordinary.wmbp_selected == taped.wmbp_selected &&
         ordinary.case3_producer_candidates ==
             taped.case3_producer_candidates &&
         ordinary.case3_producer_selected ==
             taped.case3_producer_selected &&
         ordinary.case3_tail_candidates == taped.case3_tail_candidates &&
         ordinary.case3_tail_selected == taped.case3_tail_selected &&
         ordinary.case3_left_candidates == taped.case3_left_candidates &&
         ordinary.case3_left_selected == taped.case3_left_selected &&
         same_counts(ordinary.item_visits, taped.item_visits) &&
         same_counts(ordinary.positive_items, taped.positive_items) &&
         same_counts(
             ordinary.frontier_candidates, taped.frontier_candidates) &&
         same_counts(ordinary.frontier_selected, taped.frontier_selected);
}

bool check_case(const std::string &sequence,
                const std::string &scaffold,
                std::size_t width,
                bool force_scaffold = true) {
  std::string mutable_sequence = sequence;
  std::string mutable_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      mutable_sequence, mutable_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  auto widths =
      cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(width);
  widths.force_scaffold_derivation = force_scaffold;
  const ItemKey root{NonTerminal::W, 1, n, -1, -1};

  const auto ordinary =
      cparty::scfg2::replay::run_w_final_vp_wmbp_online<
          cparty::scfg2::semirings::LogSumExpSemiring>(
          ctx, n, widths);

  cparty::scfg2::semirings::LogSumExpTape tape;
  const auto taped = [&] {
    cparty::scfg2::semirings::TapeLogSumExpSemiring::Scope scope(tape);
    return cparty::scfg2::replay::run_w_final_vp_wmbp_online<
        cparty::scfg2::semirings::TapeLogSumExpSemiring>(
        ctx, n, widths);
  }();

  const double ordinary_root = ordinary.chart.get(root);
  const auto taped_root = taped.chart.get(root);
  if (!close(ordinary_root, taped_root.log_value)) {
    std::cerr << "beam tape root mismatch n=" << n << " width=" << width
              << " ordinary=" << ordinary_root
              << " taped=" << taped_root.log_value << '\n';
    return false;
  }
  if (!same_selection(ordinary, taped)) {
    std::cerr << "beam tape selection mismatch n=" << n
              << " width=" << width << '\n';
    return false;
  }
  if (ordinary.chart.values().size() != taped.chart.values().size()) {
    std::cerr << "beam tape chart size mismatch n=" << n
              << " width=" << width << " ordinary="
              << ordinary.chart.values().size() << " taped="
              << taped.chart.values().size() << '\n';
    return false;
  }
  std::size_t finite_items = 0;
  std::map<std::pair<NonTerminal, cand_pos_t>, std::size_t> retained_counts;
  for (const auto &[item, ordinary_value] : ordinary.chart.values()) {
    if (!std::isfinite(ordinary_value)) {
      std::cerr << "pruned/zero item occupies numeric inside chart\n";
      return false;
    }
    ++retained_counts[{item.nonterminal, item.j}];
    if (!taped.chart.contains(item)) {
      std::cerr << "beam tape missing retained item\n";
      return false;
    }
    const auto taped_value = taped.chart.get(item);
    if (!close(ordinary_value, taped_value.log_value)) {
      std::cerr << "beam tape item value mismatch nt="
                << cparty::scfg2::grammar::nonterminal_name(item.nonterminal)
                << " item=(" << item.i << ',' << item.j << ',' << item.ip
                << ',' << item.jp << ") ordinary=" << ordinary_value
                << " taped=" << taped_value.log_value << '\n';
      return false;
    }
    if (std::isfinite(ordinary_value)) ++finite_items;
  }
  if (!force_scaffold) {
    for (const auto &[family_column, count] : retained_counts) {
      const auto family = family_column.first;
      const auto limit = family == NonTerminal::WMB ? 2 * width
                       : family == NonTerminal::W ? width + 1 : width;
      if (count > limit) {
        std::cerr << "retained chart exceeds family budget nt="
                  << int(family) << " j=" << family_column.second
                  << " count=" << count << " limit=" << limit << '\n';
        return false;
      }
    }
  }
  if (taped_root.node ==
      cparty::scfg2::semirings::invalid_tape_node) {
    if (!force_scaffold && !std::isfinite(ordinary_root)) return true;
    std::cerr << "beam tape root has no arithmetic node\n";
    return false;
  }
  const auto adjoints = tape.reverse(taped_root.node);
  if (!close(adjoints[taped_root.node], 1.0, 1e-12)) {
    std::cerr << "beam tape root adjoint is not one\n";
    return false;
  }
  for (double adjoint : adjoints) {
    if (!std::isfinite(adjoint) || adjoint < 0.0) {
      std::cerr << "beam tape has invalid adjoint\n";
      return false;
    }
  }
  std::size_t derivative_checks = 0;
  constexpr double epsilon = 1e-5;
  for (const auto &[item, value] : taped.chart.values()) {
    (void)item;
    if (value.node ==
            cparty::scfg2::semirings::invalid_tape_node ||
        adjoints[value.node] <= 1e-12) {
      continue;
    }
    const double finite_difference =
        (tape.evaluate_with_log_shift(
             taped_root.node, value.node, epsilon) -
         tape.evaluate_with_log_shift(
             taped_root.node, value.node, -epsilon)) /
        (2.0 * epsilon);
    if (!close(adjoints[value.node], finite_difference, 1e-7)) {
      std::cerr << "beam frozen-DAG derivative mismatch adjoint="
                << adjoints[value.node]
                << " finite_difference=" << finite_difference << '\n';
      return false;
    }
    if (++derivative_checks == 24) break;
  }
  if (derivative_checks == 0) {
    std::cerr << "beam tape checked no frozen-DAG derivatives\n";
    return false;
  }
  std::cout << "case\tn=" << n << "\twidth=" << width
            << "\tfinite_items=" << finite_items
            << "\tderivatives=" << derivative_checks
            << "\ttape_nodes=" << tape.nodes().size()
            << "\tlog_z=" << ordinary_root << '\n';
  return true;
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  if (!check_case("GAGAAACAAAAAAGAAACC", "...................", 5)) {
    return EXIT_FAILURE;
  }
  if (!check_case(
          "GCGCGCGCGCGC", "..(....)....", 5)) {
    return EXIT_FAILURE;
  }
  if (!check_case(
          "GCGCGCGCGCGC", "..(....)....", 20)) {
    return EXIT_FAILURE;
  }
  for (std::size_t width : {1, 3, 20}) {
    if (!check_case(
            "AGGUACCCUACAAUUACGGCAUCAGUCUAAAAACUCUGGGAUGGGCA",
            "(.(....).....)..(..(....)........(....)).(....)",
            width, false)) return EXIT_FAILURE;
  }
  std::cout << "scfg2_beam_outside_retained_dag=ok\n";
  return EXIT_SUCCESS;
}

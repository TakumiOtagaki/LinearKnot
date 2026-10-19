#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKeyEq;
using cparty::scfg2::engine::DeductionTraceNode;

bool same_deduction(const Deduction &lhs, const Deduction &rhs) {
  if (!ItemKeyEq{}(lhs.parent, rhs.parent) || lhs.rule != rhs.rule ||
      lhs.split.k != rhs.split.k || lhs.split.l != rhs.split.l ||
      lhs.split.p != rhs.split.p || lhs.split.q != rhs.split.q ||
      lhs.child_count != rhs.child_count) {
    return false;
  }
  for (std::size_t idx = 0; idx < lhs.child_count; ++idx) {
    if (!ItemKeyEq{}(lhs.children[idx].item, rhs.children[idx].item)) {
      return false;
    }
  }
  return true;
}

struct TraceCoverage {
  std::size_t items = 0;
  std::size_t retained_items = 0;
  std::size_t choices = 0;
  std::size_t matching_choices = 0;
};

void observe_coverage(const DeductionTraceNode<double> &exact,
                      const cparty::scfg2::engine::Chart<double> &finite_chart,
                      const cparty::scfg2::engine::ChoiceMap<double> &finite,
                      TraceCoverage *out) {
  ++out->items;
  if (finite_chart.get(exact.item) != 0.0) ++out->retained_items;
  if (exact.chosen_deduction.has_value()) {
    ++out->choices;
    const auto *choice = finite.find(exact.item);
    if (choice && choice->has_value &&
        same_deduction(*exact.chosen_deduction, choice->deduction)) {
      ++out->matching_choices;
    }
  }
  for (const auto &child : exact.children) {
    observe_coverage(child, finite_chart, finite, out);
  }
}

double ratio(double numerator, double denominator) {
  return denominator == 0.0 ? 0.0 : numerator / denominator;
}

bool no_greater_than(double value, double upper) {
  if (std::isnan(value) || std::isnan(upper)) return false;
  if (value <= upper) return true;
  const auto scale = std::max(std::abs(value), std::abs(upper));
  return value - upper <=
      64.0 * std::numeric_limits<double>::epsilon() * scale;
}

std::vector<std::size_t> parse_widths(const std::string &text) {
  std::vector<std::size_t> out;
  std::stringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (!token.empty()) out.push_back(std::stoull(token));
  }
  return out;
}

bool apply_geometry_mode(
    const std::string &mode,
    cparty::scfg2::replay::WFinalVpWmbpOnlineWidths *widths) {
  const auto unlimited = std::numeric_limits<std::size_t>::max();
  if (mode == "finite") return true;
  if (mode == "beam-only" || mode == "ancestor-only") {
    widths->max_hairpin_unpaired = unlimited;
  }
  if (mode == "beam-only" || mode == "hairpin-only") {
    widths->max_case3_ancestor_depth = unlimited;
  }
  return mode == "beam-only" || mode == "hairpin-only" ||
      mode == "ancestor-only";
}

std::size_t pair_intersection(
    const std::vector<cparty::scfg2::engine::TracebackPair> &lhs,
    const std::vector<cparty::scfg2::engine::TracebackPair> &rhs) {
  std::size_t count = 0;
  for (const auto &pair : lhs) {
    if (cparty::scfg2::engine::contains_traceback_pair(rhs, pair)) ++count;
  }
  return count;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4 || argc > 6) {
    std::cerr << "usage: " << argv[0]
              << " SEQUENCE SCAFFOLD WIDTHS_CSV [PARAM_FILE]"
                 " [finite|beam-only|hairpin-only|ancestor-only]\n";
    return EXIT_FAILURE;
  }
  const std::string sequence = argv[1];
  const std::string scaffold = argv[2];
  const auto widths = parse_widths(argv[3]);
  if (sequence.size() != scaffold.size() || widths.empty() ||
      widths.front() == 0 ||
      !std::is_sorted(widths.begin(), widths.end()) ||
      std::adjacent_find(widths.begin(), widths.end()) != widths.end()) {
    std::cerr << "invalid sequence/scaffold length or widths\n";
    return EXIT_FAILURE;
  }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc >= 5 ? argv[4] : "params/rna_Turner04.par",
         1);
  const std::string geometry_mode = argc >= 6 ? argv[5] : "finite";

  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const cparty::scfg2::replay::WFinalVpWmbpOnlineWidths keep_all_widths;
  const auto root = cparty::scfg2::deductions::ItemKey{
      cparty::scfg2::grammar::NonTerminal::W, 1, n, -1, -1};

  const auto exact_inside = cparty::scfg2::replay::run_w_final_vp_wmbp_online<
      cparty::scfg2::semirings::SumProductSemiring>(
      ctx, n, keep_all_widths);
  const auto exact_inside_root = exact_inside.chart.get(root);
  auto exact_viterbi =
      cparty::scfg2::replay::run_w_final_vp_wmbp_online_viterbi(
          ctx, n, keep_all_widths);
  const auto exact_viterbi_root = exact_viterbi.dp.chart.get(root);
  const auto exact_trace = cparty::scfg2::engine::build_traceback_tree(
      root, exact_viterbi.dp.chart, exact_viterbi.choices);
  const auto exact_pairs =
      cparty::scfg2::engine::collect_generated_traceback_pairs(
          exact_trace, ctx);

  std::cout << std::setprecision(17)
            << "n\tgeometry\twidth\texact_inside\tfinite_inside\tinside_ratio"
               "\texact_viterbi\tfinite_viterbi\tviterbi_ratio"
               "\ttrace_items\tretained_trace_items\titem_coverage"
               "\ttrace_choices\tmatching_trace_choices\tchoice_coverage"
               "\tfull_trace_match\texact_generated_pairs"
               "\tfinite_generated_pairs\tpair_intersection"
               "\tpair_precision\tpair_recall"
               "\tinside_monotone\tviterbi_monotone\n";
  double previous_inside_root = 0.0;
  double previous_viterbi_root = 0.0;
  for (const auto width : widths) {
    auto finite_widths =
        cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(width);
    if (!apply_geometry_mode(geometry_mode, &finite_widths)) {
      std::cerr << "invalid geometry mode: " << geometry_mode << '\n';
      return EXIT_FAILURE;
    }
    const auto finite_inside =
        cparty::scfg2::replay::run_w_final_vp_wmbp_online<
            cparty::scfg2::semirings::SumProductSemiring>(
            ctx, n, finite_widths);
    auto finite_viterbi =
        cparty::scfg2::replay::run_w_final_vp_wmbp_online_viterbi(
            ctx, n, finite_widths);
    const auto finite_trace = cparty::scfg2::engine::build_traceback_tree(
        root, finite_viterbi.dp.chart, finite_viterbi.choices);
    const auto finite_pairs =
        cparty::scfg2::engine::collect_generated_traceback_pairs(
            finite_trace, ctx);
    TraceCoverage coverage;
    observe_coverage(exact_trace,
                     finite_viterbi.dp.chart,
                     finite_viterbi.choices,
                     &coverage);
    const auto intersection = pair_intersection(exact_pairs, finite_pairs);
    const bool full_trace = coverage.retained_items == coverage.items &&
        coverage.matching_choices == coverage.choices;
    const auto finite_inside_root = finite_inside.chart.get(root);
    const auto finite_viterbi_root = finite_viterbi.dp.chart.get(root);
    const bool inside_monotone =
        no_greater_than(previous_inside_root, finite_inside_root);
    const bool viterbi_monotone =
        no_greater_than(previous_viterbi_root, finite_viterbi_root);
    if (!no_greater_than(finite_inside_root, exact_inside_root) ||
        !no_greater_than(finite_viterbi_root, exact_viterbi_root)) {
      std::cerr << "finite quality invariant failed at width " << width
                << ": inside(previous/finite/exact)="
                << previous_inside_root << '/' << finite_inside_root << '/'
                << exact_inside_root << ", viterbi(previous/finite/exact)="
                << previous_viterbi_root << '/' << finite_viterbi_root << '/'
                << exact_viterbi_root
                << ", endpoint-prefix(reg/late/nonmonotone/max/max-returned)="
                << finite_viterbi.dp.case3_endpoint_prefix_registrations << '/'
                << finite_viterbi.dp.case3_endpoint_prefix_late_registrations
                << '/'
                << finite_viterbi.dp.case3_endpoint_prefix_nonmonotone_queries
                << '/' << finite_viterbi.dp.case3_endpoint_prefix_maximum_size
                << '/'
                << finite_viterbi.dp
                       .case3_endpoint_prefix_maximum_returned_per_path
                << '\n';
      return EXIT_FAILURE;
    }
    std::cout << n << '\t' << geometry_mode << '\t' << width << '\t'
              << exact_inside_root << '\t'
              << finite_inside_root << '\t'
              << ratio(finite_inside_root, exact_inside_root) << '\t'
              << exact_viterbi_root << '\t' << finite_viterbi_root << '\t'
              << ratio(finite_viterbi_root, exact_viterbi_root) << '\t'
              << coverage.items << '\t' << coverage.retained_items << '\t'
              << ratio(coverage.retained_items, coverage.items) << '\t'
              << coverage.choices << '\t' << coverage.matching_choices << '\t'
              << ratio(coverage.matching_choices, coverage.choices) << '\t'
              << (full_trace ? 1 : 0) << '\t' << exact_pairs.size() << '\t'
              << finite_pairs.size() << '\t' << intersection << '\t'
              << ratio(intersection, finite_pairs.size()) << '\t'
              << ratio(intersection, exact_pairs.size()) << '\t'
              << (inside_monotone ? 1 : 0) << '\t'
              << (viterbi_monotone ? 1 : 0) << '\n';
    previous_inside_root = finite_inside_root;
    previous_viterbi_root = finite_viterbi_root;
  }
  return EXIT_SUCCESS;
}

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;
using cparty::scfg2::deductions::ItemKeyEq;

bool same_deduction(const Deduction &lhs, const Deduction &rhs) {
  if (!ItemKeyEq{}(lhs.parent, rhs.parent) || lhs.rule != rhs.rule ||
      lhs.split.k != rhs.split.k || lhs.split.l != rhs.split.l ||
      lhs.split.p != rhs.split.p || lhs.split.q != rhs.split.q ||
      lhs.child_count != rhs.child_count) {
    return false;
  }
  for (std::size_t idx = 0; idx < lhs.child_count; ++idx) {
    if (!ItemKeyEq{}(lhs.children[idx].item, rhs.children[idx].item))
      return false;
  }
  return true;
}

struct Candidate {
  Deduction deduction{};
  double contribution = 0.0;
};

using CandidateMap = std::unordered_map<
    ItemKey,
    std::vector<Candidate>,
    cparty::scfg2::deductions::ItemKeyHash,
    ItemKeyEq>;

const Candidate *find_candidate(const CandidateMap &candidates,
                                const ItemKey &item,
                                const Deduction &deduction) {
  const auto found = candidates.find(item);
  if (found == candidates.end()) return nullptr;
  for (const auto &candidate : found->second) {
    if (same_deduction(candidate.deduction, deduction)) return &candidate;
  }
  return nullptr;
}

std::pair<std::size_t, std::size_t> logsum_family_rank(
    const CandidateMap &candidates,
    const ItemKey &target) {
  std::vector<std::pair<ItemKey, double>> ranked;
  for (const auto &[item, item_candidates] : candidates) {
    if (item.nonterminal != target.nonterminal || item.j != target.j)
      continue;
    double total =
        cparty::scfg2::semirings::LogSumExpSemiring::zero();
    for (const auto &candidate : item_candidates) {
      total = cparty::scfg2::semirings::LogSumExpSemiring::plus(
          total, candidate.contribution);
    }
    ranked.push_back({item, total});
  }
  std::sort(ranked.begin(), ranked.end(), [](const auto &lhs,
                                              const auto &rhs) {
    if (lhs.second != rhs.second) return lhs.second > rhs.second;
    return cparty::scfg2::beam::right_boundary_item_tie_key(lhs.first) <
        cparty::scfg2::beam::right_boundary_item_tie_key(rhs.first);
  });
  for (std::size_t idx = 0; idx < ranked.size(); ++idx) {
    if (ItemKeyEq{}(ranked[idx].first, target))
      return {idx + 1, ranked.size()};
  }
  return {0, ranked.size()};
}

void print_item(const ItemKey &item) {
  std::cout << cparty::scfg2::grammar::nonterminal_name(item.nonterminal)
            << '(' << item.i << ',' << item.j;
  if (item.ip >= 0 || item.jp >= 0)
    std::cout << ';' << item.ip << ',' << item.jp;
  std::cout << ')';
}

void diagnose(
    const cparty::scfg2::engine::DeductionTraceNode<double> &node,
    const cparty::scfg2::engine::Chart<double> &inside_chart,
    const CandidateMap &inside_candidates,
    std::vector<ItemKey> *path,
    std::size_t *missing) {
  path->push_back(node.item);
  if (!node.chosen_deduction.has_value() && std::isfinite(node.chart_value)) {
    std::cout << "NO_MAX_CHOICE depth=" << path->size() - 1 << " path=";
    for (const auto &item : *path) {
      print_item(item);
      std::cout << '/';
    }
    std::cout << " logmax_item=" << node.chart_value
              << " logsum_item=" << inside_chart.get(node.item) << '\n';
  }
  if (std::isfinite(node.chart_value) &&
      !std::isfinite(inside_chart.get(node.item))) {
    const auto rank = logsum_family_rank(inside_candidates, node.item);
    std::cout << "PRUNED_IN_LOGSUM depth=" << path->size() - 1
              << " path=";
    for (const auto &item : *path) {
      print_item(item);
      std::cout << '/';
    }
    std::cout << " logmax_item=" << node.chart_value
              << " logsum_family_rank=" << rank.first << '/' << rank.second
              << '\n';
  }
  if (node.chosen_deduction.has_value() &&
      !find_candidate(
          inside_candidates, node.item, *node.chosen_deduction)) {
    ++*missing;
    std::cout << "MISSING depth=" << path->size() - 1 << " path=";
    for (const auto &item : *path) {
      print_item(item);
      std::cout << '/';
    }
    std::cout << " rule=" << cparty::scfg2::grammar::rule_name(
        node.chosen_deduction->rule)
              << " logmax_item=" << node.chart_value
              << " logsum_item=" << inside_chart.get(node.item)
              << " logsum_parent_candidate_count=";
    const auto found = inside_candidates.find(node.item);
    std::cout << (found == inside_candidates.end() ? 0 : found->second.size())
              << " children=";
    for (std::size_t idx = 0;
         idx < node.chosen_deduction->child_count;
         ++idx) {
      const auto &child = node.chosen_deduction->children[idx].item;
      print_item(child);
      std::cout << ':' << inside_chart.get(child) << ',';
    }
    std::cout << '\n';
  }
  for (const auto &child : node.children)
    diagnose(child, inside_chart, inside_candidates, path, missing);
  path->pop_back();
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4 || argc > 5) {
    std::cerr << "usage: " << argv[0]
              << " SEQUENCE SCAFFOLD WIDTH [PARAM_FILE]\n";
    return EXIT_FAILURE;
  }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 5 ? argv[4] : "params/rna_Turner04.par",
         1);
  std::string sequence = argv[1];
  std::string scaffold = argv[2];
  if (sequence.size() != scaffold.size()) return EXIT_FAILURE;
  const std::size_t width = std::stoull(argv[3]);
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(sequence.size());
  const auto widths =
      cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(width);

  auto max_run = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_viterbi_with_semiring<
          cparty::scfg2::semirings::LogMaxSemiring>(ctx, n, widths);
  CandidateMap inside_candidates;
  auto inside_run = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_observed<
          cparty::scfg2::semirings::LogSumExpSemiring>(
          ctx,
          n,
          widths,
          [&](const ItemKey &item,
              const Deduction &deduction,
              double contribution,
              const auto &,
              double) {
            if (!std::isfinite(contribution)) return;
            inside_candidates[item].push_back(Candidate{
                deduction,
                contribution});
          });
  const ItemKey root{
      cparty::scfg2::grammar::NonTerminal::W, 1, n, -1, -1};
  const auto trace = cparty::scfg2::engine::build_traceback_tree(
      root, max_run.dp.chart, max_run.choices);
  std::cout << std::setprecision(17)
            << "n=" << n << " width=" << width
            << " logsum_root=" << inside_run.chart.get(root)
            << " logmax_root=" << max_run.dp.chart.get(root)
            << " delta=" << max_run.dp.chart.get(root) -
                                  inside_run.chart.get(root)
            << '\n';
  std::vector<ItemKey> path;
  std::size_t missing = 0;
  diagnose(trace, inside_run.chart, inside_candidates, &path, &missing);
  std::cout << "missing_selected_deductions=" << missing << '\n';
  return EXIT_SUCCESS;
}

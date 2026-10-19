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
    if (!ItemKeyEq{}(lhs.children[idx].item, rhs.children[idx].item)) {
      return false;
    }
  }
  return true;
}

struct Candidate {
  Deduction deduction{};
  double contribution = 0.0;
};

using CandidateMap = std::unordered_map<ItemKey,
                                        std::vector<Candidate>,
                                        cparty::scfg2::deductions::ItemKeyHash,
                                        ItemKeyEq>;

struct ObservedRun {
  cparty::scfg2::replay::WFinalVpWmbpOnlineRun<double> dp;
  cparty::scfg2::engine::ChoiceMap<double> choices;
  CandidateMap candidates;
};

template <class Context>
ObservedRun run_observed(
    Context &ctx,
    cand_pos_t n,
    std::size_t width) {
  cparty::scfg2::engine::ChoiceMap<double> choices;
  CandidateMap candidates;
  auto dp = cparty::scfg2::replay::run_w_final_vp_wmbp_online_observed<
      cparty::scfg2::semirings::LogMaxSemiring>(
      ctx,
      n,
      cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(width),
      [&](const ItemKey &item,
          const Deduction &deduction,
          double contribution,
          const auto &,
          double) {
        if (!std::isfinite(contribution)) return;
        candidates[item].push_back(Candidate{deduction, contribution});
        choices.observe(
            item,
            deduction,
            contribution,
            [](double lhs, double rhs) { return lhs > rhs; });
      });
  return ObservedRun{
      std::move(dp), std::move(choices), std::move(candidates)};
}

void print_item(const ItemKey &item) {
  std::cout << cparty::scfg2::grammar::nonterminal_name(item.nonterminal)
            << '(' << item.i << ',' << item.j;
  if (item.ip >= 0 || item.jp >= 0) {
    std::cout << ';' << item.ip << ',' << item.jp;
  }
  std::cout << ')';
}

const Candidate *find_candidate(
    const CandidateMap &candidates,
    const ItemKey &item,
    const Deduction &deduction) {
  const auto found = candidates.find(item);
  if (found == candidates.end()) return nullptr;
  for (const auto &candidate : found->second) {
    if (same_deduction(candidate.deduction, deduction)) return &candidate;
  }
  return nullptr;
}

std::pair<std::size_t, std::size_t> family_column_rank(
    const ObservedRun &run,
    const ItemKey &target) {
  const auto selection_priority = [&](const ItemKey &item,
                                      double value) {
    if (item.nonterminal !=
        cparty::scfg2::grammar::NonTerminal::V) {
      return value;
    }
    if (item.i <= 1) return value;
    return cparty::scfg2::semirings::LogMaxSemiring::times(
        value,
        run.dp.chart.get(ItemKey{
            cparty::scfg2::grammar::NonTerminal::W,
            1,
            item.i - 1,
            -1,
            -1}));
  };
  std::vector<std::pair<ItemKey, double>> ranked;
  for (const auto &[item, choice] : run.choices.values()) {
    if (!choice.has_value || item.nonterminal != target.nonterminal ||
        item.j != target.j) {
      continue;
    }
    ranked.push_back(
        {item, selection_priority(item, choice.contribution)});
  }
  std::sort(ranked.begin(), ranked.end(), [](const auto &lhs,
                                              const auto &rhs) {
    if (lhs.second != rhs.second) return lhs.second > rhs.second;
    return cparty::scfg2::beam::right_boundary_item_tie_key(lhs.first) <
        cparty::scfg2::beam::right_boundary_item_tie_key(rhs.first);
  });
  for (std::size_t idx = 0; idx < ranked.size(); ++idx) {
    if (ItemKeyEq{}(ranked[idx].first, target)) {
      return {idx + 1, ranked.size()};
    }
  }
  return {0, ranked.size()};
}

void diagnose_trace(
    const cparty::scfg2::engine::DeductionTraceNode<double> &node,
    const ObservedRun &narrow,
    const ObservedRun &wide,
    std::vector<ItemKey> *path,
    std::size_t *loss_count,
    std::size_t *choice_difference_count) {
  path->push_back(node.item);
  const auto wide_value = wide.dp.chart.get(node.item);
  const auto *wide_choice = wide.choices.find(node.item);
  const Candidate *retained_candidate = nullptr;
  if (node.chosen_deduction.has_value()) {
    retained_candidate = find_candidate(
        wide.candidates, node.item, *node.chosen_deduction);
  }
  const bool choice_same =
      node.chosen_deduction.has_value() && wide_choice &&
      wide_choice->has_value &&
      same_deduction(*node.chosen_deduction, wide_choice->deduction);
  if (!choice_same) ++*choice_difference_count;
  if (!std::isfinite(wide_value) ||
      (node.chosen_deduction.has_value() && !retained_candidate)) {
    ++*loss_count;
    std::cout << "LOSS depth=" << path->size() - 1 << " path=";
    for (const auto &item : *path) {
      print_item(item);
      std::cout << '/';
    }
    std::cout << " narrow_value=" << node.chart_value
              << " wide_value=" << wide_value;
    if (node.chosen_deduction.has_value()) {
      std::cout << " narrow_rule="
                << cparty::scfg2::grammar::rule_name(
                       node.chosen_deduction->rule);
    }
    if (wide_choice && wide_choice->has_value) {
      std::cout << " wide_rule="
                << cparty::scfg2::grammar::rule_name(
                       wide_choice->deduction.rule)
                << " wide_choice=" << wide_choice->contribution;
    }
    const auto candidates = wide.candidates.find(node.item);
    std::cout << " wide_candidate_count="
              << (candidates == wide.candidates.end()
                      ? 0
                      : candidates->second.size())
              ;
    if (!std::isfinite(wide_value) &&
        wide_choice && wide_choice->has_value) {
      const auto narrow_rank = family_column_rank(narrow, node.item);
      const auto wide_rank = family_column_rank(wide, node.item);
      std::cout << " narrow_family_rank=" << narrow_rank.first << '/'
                << narrow_rank.second
                << " wide_family_rank=" << wide_rank.first << '/'
                << wide_rank.second;
    }
    std::cout << '\n';
  } else if (!choice_same && node.chosen_deduction.has_value()) {
    std::cout << "RERANK depth=" << path->size() - 1 << ' ';
    print_item(node.item);
    std::cout << " narrow_value=" << node.chart_value
              << " wide_value=" << wide_value
              << " narrow_rule="
              << cparty::scfg2::grammar::rule_name(
                     node.chosen_deduction->rule)
              << " narrow_rule_at_wide=" << retained_candidate->contribution;
    if (wide_choice && wide_choice->has_value) {
      std::cout << " wide_rule="
                << cparty::scfg2::grammar::rule_name(
                       wide_choice->deduction.rule)
                << " wide_choice=" << wide_choice->contribution;
    }
    std::cout << '\n';
  }
  for (const auto &child : node.children) {
    diagnose_trace(
        child, narrow, wide, path, loss_count, choice_difference_count);
  }
  path->pop_back();
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 5 || argc > 6) {
    std::cerr << "usage: " << argv[0]
              << " SEQUENCE SCAFFOLD NARROW_WIDTH WIDE_WIDTH [PARAM_FILE]\n";
    return EXIT_FAILURE;
  }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 6 ? argv[5] : "params/rna_Turner04.par",
         1);
  std::string sequence = argv[1];
  std::string scaffold = argv[2];
  if (scaffold == "empty") scaffold = std::string(sequence.size(), '.');
  if (sequence.size() != scaffold.size()) return EXIT_FAILURE;
  const auto narrow_width = std::stoull(argv[3]);
  const auto wide_width = std::stoull(argv[4]);
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  auto narrow = run_observed(ctx, n, narrow_width);
  auto wide = run_observed(ctx, n, wide_width);
  const ItemKey root{
      cparty::scfg2::grammar::NonTerminal::W, 1, n, -1, -1};
  const auto trace = cparty::scfg2::engine::build_traceback_tree(
      root, narrow.dp.chart, narrow.choices);
  std::cout << std::setprecision(17)
            << "n=" << n << " narrow_width=" << narrow_width
            << " wide_width=" << wide_width
            << " narrow_root=" << narrow.dp.chart.get(root)
            << " wide_root=" << wide.dp.chart.get(root) << '\n';
  std::vector<ItemKey> path;
  std::size_t loss_count = 0;
  std::size_t choice_difference_count = 0;
  diagnose_trace(
      trace, narrow, wide, &path, &loss_count, &choice_difference_count);
  std::cout << "loss_count=" << loss_count
            << " choice_difference_count=" << choice_difference_count
            << '\n';
  return EXIT_SUCCESS;
}

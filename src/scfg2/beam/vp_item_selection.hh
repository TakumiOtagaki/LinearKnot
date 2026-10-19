#ifndef CPARTY_SCFG2_BEAM_VP_ITEM_SELECTION_HH
#define CPARTY_SCFG2_BEAM_VP_ITEM_SELECTION_HH

#include "scfg2/beam/right_boundary_item_frontier.hh"
#include "scfg2/beam/wmbp_be_vp_prefix_profile.hh"

#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

class TopWidthVpItemSelection {
 public:
  explicit TopWidthVpItemSelection(std::size_t width) : frontier_(width) {}

  void observe_completed(const deductions::ItemKey &item, double value) {
    observe_ranked(item, value != 0.0, value);
  }

  void observe_ranked(const deductions::ItemKey &item,
                      bool live,
                      double priority) {
    if (item.nonterminal != grammar::NonTerminal::VP || !live) return;
    frontier_.observe(item, priority);
  }

  void commit_column(cand_pos_t j) {
    commit_column_partition(j, [](const deductions::ItemKey &) {});
  }

  template <class Chart, class Zero>
  void commit_column(cand_pos_t j, Chart &chart, Zero zero) {
    commit_column_partition(j, [&](const deductions::ItemKey &item) {
      chart.set(item, zero);
    });
  }

  template <class PrunedFn>
  void commit_column_partition(cand_pos_t j, PrunedFn &&on_pruned) {
    auto &starts = selected_starts_by_column_[j];
    frontier_.flush_column_partition(j, [&](const deductions::ItemKey &item) {
      selected_.insert(item);
      starts.push_back(item.i);
    }, on_pruned);
    std::sort(starts.begin(), starts.end());
  }

  bool allows(const deductions::ItemKey &vp) const {
    return selected_.find(vp) != selected_.end();
  }

  bool allows_wmbp_source(const deductions::Deduction &deduction) const {
    if (deduction.rule == grammar::RuleId::WMBP_DIRECT_VP) {
      return allows(deductions::ItemKey{grammar::NonTerminal::VP,
                                        deduction.parent.i,
                                        deduction.parent.j,
                                        -1,
                                        -1});
    }
    if (!is_wmbp_be_vp_rule(deduction.rule)) return true;
    return allows(wmbp_be_vp_right_child(deduction));
  }

  const std::vector<cand_pos_t> &selected_starts(cand_pos_t j) const {
    const auto found = selected_starts_by_column_.find(j);
    if (found == selected_starts_by_column_.end()) {
      static const std::vector<cand_pos_t> empty;
      return empty;
    }
    return found->second;
  }

  std::size_t candidate_count() const { return frontier_.candidate_count(); }
  std::size_t selected_count() const { return frontier_.selected_count(); }
  std::size_t pruned_count() const { return frontier_.pruned_count(); }
  bool respects_width() const { return frontier_.respects_width(); }

 private:
  RightBoundaryTopWidthItemFrontier frontier_;
  std::unordered_set<deductions::ItemKey,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      selected_;
  std::unordered_map<cand_pos_t, std::vector<cand_pos_t>>
      selected_starts_by_column_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_ITEM_SELECTION_HH

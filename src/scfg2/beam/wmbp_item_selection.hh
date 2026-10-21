#ifndef CPARTY_SCFG2_BEAM_WMBP_ITEM_SELECTION_HH
#define CPARTY_SCFG2_BEAM_WMBP_ITEM_SELECTION_HH

#include "scfg2/beam/right_boundary_item_frontier.hh"
#include "scfg2/deductions/deduction.hh"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

class TopWidthWmbpItemSelection {
 public:
  static TopWidthWmbpItemSelection keep_all() {
    return TopWidthWmbpItemSelection(std::numeric_limits<std::size_t>::max());
  }

  explicit TopWidthWmbpItemSelection(std::size_t width) : frontier_(width) {}
  void enable_numeric_canonical_storage(cand_pos_t n) {
    frontier_.enable_numeric_canonical_storage(grammar::NonTerminal::WMBP, n);
  }

  void observe_completed(const deductions::ItemKey &item, double value) {
    observe_ranked(item, value != 0.0, value);
  }

  void observe_ranked(const deductions::ItemKey &item,
                      bool live,
                      double priority) {
    if (item.nonterminal != grammar::NonTerminal::WMBP || !live) {
      return;
    }
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
    frontier_.flush_column_partition(j, [&](const deductions::ItemKey &item) {
      selected_.insert(item);
      if (selected_starts_by_right_.size() <=
          static_cast<std::size_t>(item.j)) {
        selected_starts_by_right_.resize(
            static_cast<std::size_t>(item.j) + 1);
      }
      selected_starts_by_right_[static_cast<std::size_t>(item.j)].push_back(
          item.i);
      auto &ends = selected_ends_by_left_[item.i];
      // Production commits increasing columns. Keep the older helper's
      // out-of-order contract as a fallback without searching growing
      // histories on the monotone path.
      if (ends.empty() || ends.back() <= item.j) ends.push_back(item.j);
      else ends.insert(std::lower_bound(ends.begin(), ends.end(), item.j), item.j);
    }, on_pruned);
  }

  const std::vector<cand_pos_t> &selected_starts(cand_pos_t j) const {
    if (j < 0 || static_cast<std::size_t>(j) >=
                     selected_starts_by_right_.size()) {
      static const std::vector<cand_pos_t> empty;
      return empty;
    }
    return selected_starts_by_right_[static_cast<std::size_t>(j)];
  }

  bool allows(const deductions::ItemKey &wmbp) const {
    return selected_.find(wmbp) != selected_.end();
  }

  const std::vector<cand_pos_t> &selected_ends(cand_pos_t i) const {
    const auto found = selected_ends_by_left_.find(i);
    if (found == selected_ends_by_left_.end()) {
      static const std::vector<cand_pos_t> empty;
      return empty;
    }
    return found->second;
  }

  bool allows_consumer_source(
      const deductions::Deduction &deduction) const {
    // This predicate remains the compatibility contract for the older
    // consumer-filter path.  The online retained-join executor separately
    // uses selected_starts() to beam recursive shared-BE WMBP sources.
    if (deduction.parent.nonterminal == grammar::NonTerminal::WMBP) {
      return true;
    }
    for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
      const auto &child = deduction.children[idx].item;
      if (child.nonterminal == grammar::NonTerminal::WMBP && !allows(child)) {
        return false;
      }
    }
    return true;
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
      selected_ends_by_left_;
  std::vector<std::vector<cand_pos_t>> selected_starts_by_right_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMBP_ITEM_SELECTION_HH

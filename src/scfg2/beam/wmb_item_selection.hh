#ifndef CPARTY_SCFG2_BEAM_WMB_ITEM_SELECTION_HH
#define CPARTY_SCFG2_BEAM_WMB_ITEM_SELECTION_HH

#include "scfg2/beam/right_boundary_item_frontier.hh"

#include <cstddef>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

// WMB has two genuinely different consumers.  WMp uses WMB as an internal
// child and needs a local-value view, while W uses it as a suffix and benefits
// from a completed-root-prefix view.  Retain the union of both top-width
// views in the shared chart, then expose each view to its matching consumer.
class TopWidthWmbItemSelection {
 public:
  explicit TopWidthWmbItemSelection(std::size_t width)
      : internal_(width), root_suffix_(width) {}
  void enable_numeric_canonical_storage(cand_pos_t n) {
    internal_.enable_numeric_canonical_storage(grammar::NonTerminal::WMB, n);
    root_suffix_.enable_numeric_canonical_storage(grammar::NonTerminal::WMB, n);
  }

  void observe_completed(const deductions::ItemKey &item,
                         double value,
                         double root_suffix_priority) {
    observe_ranked(
        item, value != 0.0, value, root_suffix_priority);
  }

  void observe_ranked(const deductions::ItemKey &item,
                      bool live,
                      double internal_priority,
                      double root_suffix_priority) {
    if (item.nonterminal != grammar::NonTerminal::WMB || !live) return;
    internal_.observe(item, internal_priority);
    root_suffix_.observe(item, root_suffix_priority);
  }

  template <class Chart, class Zero>
  void commit_column(cand_pos_t j, Chart &chart, Zero zero) {
    ensure_column(j);
    auto &internal_starts = internal_starts_[static_cast<std::size_t>(j)];
    auto &root_starts = root_starts_[static_cast<std::size_t>(j)];
    std::unordered_set<deductions::ItemKey,
                       deductions::ItemKeyHash,
                       deductions::ItemKeyEq>
        retained;
    std::unordered_set<deductions::ItemKey,
                       deductions::ItemKeyHash,
                       deductions::ItemKeyEq>
        pruned;
    internal_.flush_column_partition(
        j,
        [&](const deductions::ItemKey &item) {
          internal_starts.push_back(item.i);
          retained.insert(item);
        },
        [&](const deductions::ItemKey &item) { pruned.insert(item); });
    root_suffix_.flush_column_partition(
        j,
        [&](const deductions::ItemKey &item) {
          root_starts.push_back(item.i);
          retained.insert(item);
        },
        [&](const deductions::ItemKey &item) { pruned.insert(item); });
    for (const auto &item : pruned) {
      if (retained.find(item) == retained.end()) chart.set(item, zero);
    }
    selected_count_ += retained.size();
    auto &membership = root_membership_[static_cast<std::size_t>(j)];
    membership = root_starts;
    std::sort(membership.begin(), membership.end());
  }

  const std::vector<cand_pos_t> &internal_starts(cand_pos_t j) const {
    return starts_at(internal_starts_, j);
  }

  const std::vector<cand_pos_t> &root_suffix_starts(cand_pos_t j) const {
    return starts_at(root_starts_, j);
  }

  bool root_suffix_contains(cand_pos_t i, cand_pos_t j) const {
    const auto &starts = starts_at(root_membership_, j);
    return std::binary_search(starts.begin(), starts.end(), i);
  }

  std::size_t candidate_count() const { return internal_.candidate_count(); }
  std::size_t selected_count() const { return selected_count_; }

 private:
  void ensure_column(cand_pos_t j) {
    const auto size = static_cast<std::size_t>(j) + 1;
    if (internal_starts_.size() < size) internal_starts_.resize(size);
    if (root_starts_.size() < size) root_starts_.resize(size);
    if (root_membership_.size() < size) root_membership_.resize(size);
  }

  static const std::vector<cand_pos_t> &starts_at(
      const std::vector<std::vector<cand_pos_t>> &columns,
      cand_pos_t j) {
    if (j < 0 || static_cast<std::size_t>(j) >= columns.size()) {
      static const std::vector<cand_pos_t> empty;
      return empty;
    }
    return columns[static_cast<std::size_t>(j)];
  }

  RightBoundaryTopWidthItemFrontier internal_;
  RightBoundaryTopWidthItemFrontier root_suffix_;
  std::size_t selected_count_ = 0;
  std::vector<std::vector<cand_pos_t>> internal_starts_;
  std::vector<std::vector<cand_pos_t>> root_starts_;
  std::vector<std::vector<cand_pos_t>> root_membership_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMB_ITEM_SELECTION_HH

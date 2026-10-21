#ifndef CPARTY_SCFG2_BEAM_FAMILY_ITEM_SELECTION_HH
#define CPARTY_SCFG2_BEAM_FAMILY_ITEM_SELECTION_HH

#include "scfg2/beam/right_boundary_item_frontier.hh"

#include <cstddef>
#include <memory>
#include <vector>

namespace cparty::scfg2::beam {

// One semantic nonterminal view per right boundary.  This is the Phase-5
// retained-support boundary: candidate enumeration may still be exact/broad,
// but pruned items are removed from the chart before later family phases read
// them.  Active-child agenda construction is deliberately a separate step.
class TopWidthFamilyItemSelection {
 public:
  TopWidthFamilyItemSelection(grammar::NonTerminal family,
                              std::size_t width)
      : family_(family), frontier_(width) {}

  void set_forced(std::shared_ptr<const ItemKeySet> forced) {
    frontier_.set_forced(std::move(forced));
  }
  void enable_numeric_canonical_storage(cand_pos_t n) {
    frontier_.enable_numeric_canonical_storage(family_, n);
  }

  void observe_ranked(const deductions::ItemKey &item,
                      bool live,
                      double priority) {
    if (item.nonterminal != family_ || !live) return;
    frontier_.observe(item, priority);
  }

  template <class Chart, class Zero>
  void commit_column(cand_pos_t j, Chart &chart, Zero zero) {
    if (retained_starts_by_column_.size() <= static_cast<std::size_t>(j)) {
      retained_starts_by_column_.resize(static_cast<std::size_t>(j) + 1);
    }
    auto &starts = retained_starts_by_column_[static_cast<std::size_t>(j)];
    frontier_.flush_column_partition(
        j,
        [&](const deductions::ItemKey &item) { starts.push_back(item.i); },
        [&](const deductions::ItemKey &item) { chart.set(item, zero); });
  }

  const std::vector<cand_pos_t> &retained_starts(cand_pos_t j) const {
    const auto idx = static_cast<std::size_t>(j);
    if (idx >= retained_starts_by_column_.size()) {
      static const std::vector<cand_pos_t> empty;
      return empty;
    }
    return retained_starts_by_column_[idx];
  }

  std::size_t candidate_count() const { return frontier_.candidate_count(); }
  std::size_t selected_count() const { return frontier_.selected_count(); }

  ContainerMemory memory_snapshot(const char *name) const {
    auto out = frontier_.memory_snapshot(name);
    out.add(nested_vector_memory("", retained_starts_by_column_));
    return out;
  }

 private:
  grammar::NonTerminal family_;
  RightBoundaryTopWidthItemFrontier frontier_;
  std::vector<std::vector<cand_pos_t>> retained_starts_by_column_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_FAMILY_ITEM_SELECTION_HH

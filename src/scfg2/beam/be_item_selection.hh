#ifndef CPARTY_SCFG2_BEAM_BE_ITEM_SELECTION_HH
#define CPARTY_SCFG2_BEAM_BE_ITEM_SELECTION_HH

#include "scfg2/beam/right_boundary_item_frontier.hh"

#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

class TopWidthBeItemSelection {
 public:
  explicit TopWidthBeItemSelection(std::size_t width) : frontier_(width) {}

  void observe_ranked(const deductions::ItemKey &item,
                      bool live,
                      double priority) {
    if (item.nonterminal == grammar::NonTerminal::BE && live)
      frontier_.observe(item, priority);
  }

  template <class Chart, class Zero>
  void commit_column(cand_pos_t j, Chart &chart, Zero zero) {
    if (retained_.size() <= static_cast<std::size_t>(j))
      retained_.resize(static_cast<std::size_t>(j) + 1);
    auto &items = retained_[static_cast<std::size_t>(j)];
    frontier_.flush_column_partition(
        j,
        [&](const deductions::ItemKey &item) { items.push_back(item); },
        [&](const deductions::ItemKey &item) { chart.set(item, zero); });
  }

  const std::vector<deductions::ItemKey> &retained_items(cand_pos_t j) const {
    if (j < 0 || static_cast<std::size_t>(j) >= retained_.size()) {
      static const std::vector<deductions::ItemKey> empty;
      return empty;
    }
    return retained_[static_cast<std::size_t>(j)];
  }
  std::size_t candidate_count() const { return frontier_.candidate_count(); }
  std::size_t selected_count() const { return frontier_.selected_count(); }

 private:
  RightBoundaryTopWidthItemFrontier frontier_;
  std::vector<std::vector<deductions::ItemKey>> retained_;
};

}  // namespace cparty::scfg2::beam

#endif

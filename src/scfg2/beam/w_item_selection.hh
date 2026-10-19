#ifndef CPARTY_SCFG2_BEAM_W_ITEM_SELECTION_HH
#define CPARTY_SCFG2_BEAM_W_ITEM_SELECTION_HH

#include "scfg2/beam/right_boundary_item_frontier.hh"

#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

class TopWidthWItemSelection {
 public:
  explicit TopWidthWItemSelection(std::size_t nonroot_width)
      : nonroot_(nonroot_width) {}

  void observe_ranked(const deductions::ItemKey &item,
                      bool live,
                      double priority) {
    if (item.nonterminal != grammar::NonTerminal::W || !live) return;
    ++candidates_;
    if (item.i == 1) {
      root_ = item;
      root_live_ = true;
      return;
    }
    nonroot_.observe(item, priority);
  }

  template <class Chart, class Zero>
  void commit_nonroot_column(cand_pos_t j, Chart &chart, Zero zero) {
    if (retained_.size() <= static_cast<std::size_t>(j)) {
      retained_.resize(static_cast<std::size_t>(j) + 1);
    }
    auto &starts = retained_[static_cast<std::size_t>(j)];
    nonroot_.flush_column_partition(
        j,
        [&](const deductions::ItemKey &item) {
          starts.push_back(item.i);
          ++selected_;
        },
        [&](const deductions::ItemKey &item) { chart.set(item, zero); });
  }

  template <class Chart, class Zero>
  void commit_column(cand_pos_t j, Chart &chart, Zero zero) {
    commit_nonroot_column(j, chart, zero);
    auto &starts = retained_[static_cast<std::size_t>(j)];
    if (root_.j == j && root_live_) {
      starts.push_back(1);
      ++selected_;
    }
    root_ = deductions::ItemKey{};
    root_live_ = false;
  }

  const std::vector<cand_pos_t> &retained_starts(cand_pos_t j) const {
    if (j < 0 || static_cast<std::size_t>(j) >= retained_.size()) {
      static const std::vector<cand_pos_t> empty;
      return empty;
    }
    return retained_[static_cast<std::size_t>(j)];
  }

  std::size_t candidate_count() const { return candidates_; }
  std::size_t selected_count() const { return selected_; }

 private:
  RightBoundaryTopWidthItemFrontier nonroot_;
  deductions::ItemKey root_{};
  bool root_live_ = false;
  std::size_t candidates_ = 0;
  std::size_t selected_ = 0;
  std::vector<std::vector<cand_pos_t>> retained_;
};

}  // namespace cparty::scfg2::beam

#endif

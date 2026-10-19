#ifndef CPARTY_SCFG2_ENGINE_GENERATED_TARGET_INDEX_HH
#define CPARTY_SCFG2_ENGINE_GENERATED_TARGET_INDEX_HH

#include "scfg2/engine/traceback_pairs.hh"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

namespace cparty::scfg2::engine {

// Observation-only index of the fixed generated output. A merge-sort tree
// counts target pairs wholly inside [lo,hi] in O(log^2(P+2)), after
// O(P log(P+2)) construction and storage. Deduction checks use a constant
// number of intervals, rather than scanning the target for every item.
class GeneratedTargetIndex {
 public:
  explicit GeneratedTargetIndex(std::vector<TracebackPair> target)
      : target_(std::move(target)) {
    std::sort(target_.begin(), target_.end(), before);
    target_.erase(std::unique(target_.begin(), target_.end(), traceback_pair_equal),
                  target_.end());
    tree_.resize(4 * target_.size());
    if (!target_.empty()) build(0, 0, target_.size());
  }

  struct Stats {
    std::size_t interval_queries = 0;
    std::size_t range_node_visits = 0;
  };
  const Stats &stats() const { return stats_; }

  template <class Context>
  bool exactly_partitions(const deductions::Deduction &deduction,
                          const Context &ctx) const {
    const auto required = regions(deduction.parent);
    std::vector<Interval> covered;
    std::vector<Interval> clipped;
    for (std::size_t c = 0; c < deduction.child_count; ++c) {
      for (const auto &child : regions(deduction.children[c].item)) {
        covered.push_back(child);
        for (const auto &parent : required)
          clipped.emplace_back(std::max(child.first, parent.first),
                               std::min(child.second, parent.second));
      }
    }
    const auto covered_count = count_union(covered);
    if (covered_count != count_union(clipped)) return false;

    std::vector<TracebackPair> emitted;
    append_deduction_generated_pairs(deduction, ctx, &emitted);
    std::size_t new_emissions = 0;
    for (const auto &pair : emitted) {
      if (!std::binary_search(target_.begin(), target_.end(), pair, before) ||
          !owns(required, pair)) return false;
      if (!owns(covered, pair)) ++new_emissions;
    }
    return count_union(required) == covered_count + new_emissions;
  }

 private:
  using Interval = std::pair<cand_pos_t, cand_pos_t>;

  static bool before(const TracebackPair &a, const TracebackPair &b) {
    return traceback_pair_identity(a) < traceback_pair_identity(b);
  }

  static std::vector<Interval> regions(const deductions::ItemKey &item) {
    if (item.nonterminal == grammar::NonTerminal::VM)
      return {{item.i + 1, item.j - 1}};
    if (item.nonterminal == grammar::NonTerminal::BE)
      return {{item.i + 1, item.ip - 1}, {item.jp + 1, item.j - 1}};
    return {{item.i, item.j}};
  }

  static bool owns(const std::vector<Interval> &regions, const TracebackPair &pair) {
    for (const auto &[lo, hi] : regions)
      if (lo <= pair.i && pair.j <= hi) return true;
    return false;
  }

  // Pair ownership is a union of "both endpoints inside one interval".
  // Adjacent intervals must not be merged: that would admit spanning pairs.
  // In opening-order, overlap with all preceding intervals is determined by
  // their largest closing boundary, so only one intersection is subtracted.
  std::size_t count_union(std::vector<Interval> intervals) const {
    std::sort(intervals.begin(), intervals.end());
    cand_pos_t previous_hi = std::numeric_limits<cand_pos_t>::min();
    std::size_t total = 0;
    for (const auto &[lo, hi] : intervals) {
      if (lo > hi) continue;
      total += count_interval(lo, hi);
      total -= count_interval(lo, std::min(hi, previous_hi));
      previous_hi = std::max(previous_hi, hi);
    }
    return total;
  }

  std::size_t count_interval(cand_pos_t lo, cand_pos_t hi) const {
    ++stats_.interval_queries;
    if (lo > hi || target_.empty()) return 0;
    const auto first = std::lower_bound(target_.begin(), target_.end(), lo,
        [](const TracebackPair &p, cand_pos_t x) { return p.i < x; });
    const auto last = std::upper_bound(first, target_.end(), hi,
        [](cand_pos_t x, const TracebackPair &p) { return x < p.i; });
    if (first == last) return 0;
    return count_range(0, 0, target_.size(), first - target_.begin(),
                       last - target_.begin(), hi);
  }

  void build(std::size_t node, std::size_t lo, std::size_t hi) {
    if (hi - lo == 1) {
      tree_[node].push_back(target_[lo].j);
      return;
    }
    const auto middle = lo + (hi - lo) / 2;
    const auto left = 2 * node + 1, right = left + 1;
    build(left, lo, middle);
    build(right, middle, hi);
    tree_[node].resize(hi - lo);
    std::merge(tree_[left].begin(), tree_[left].end(),
               tree_[right].begin(), tree_[right].end(), tree_[node].begin());
  }

  std::size_t count_range(std::size_t node, std::size_t lo, std::size_t hi,
                          std::size_t first, std::size_t last,
                          cand_pos_t closing) const {
    ++stats_.range_node_visits;
    if (first >= hi || last <= lo) return 0;
    if (first <= lo && hi <= last)
      return std::upper_bound(tree_[node].begin(), tree_[node].end(), closing)
             - tree_[node].begin();
    const auto middle = lo + (hi - lo) / 2;
    return count_range(2 * node + 1, lo, middle, first, last, closing) +
           count_range(2 * node + 2, middle, hi, first, last, closing);
  }

  std::vector<TracebackPair> target_;
  std::vector<std::vector<cand_pos_t>> tree_;
  mutable Stats stats_;
};

}  // namespace cparty::scfg2::engine

#endif

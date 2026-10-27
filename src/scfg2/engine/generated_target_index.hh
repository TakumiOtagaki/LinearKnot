#ifndef CPARTY_SCFG2_ENGINE_GENERATED_TARGET_INDEX_HH
#define CPARTY_SCFG2_ENGINE_GENERATED_TARGET_INDEX_HH

#include "scfg2/engine/traceback_pairs.hh"
#include "scfg2/profiles/loop_context.hh"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

namespace cparty::scfg2::engine {

// Observation-only index of the fixed generated output. Pair membership and
// emitted-unpaired endpoint queries support complete-output conditioning.
// A merge-sort tree also retains the narrower interval-partition observation,
// counting pairs inside [lo,hi] in O(log^2(P+2)), after O(P log(P+2)) setup.
class GeneratedTargetIndex {
 public:
  explicit GeneratedTargetIndex(std::vector<TracebackPair> target)
      : target_(std::move(target)) {
    std::sort(target_.begin(), target_.end(), before);
    target_.erase(std::unique(target_.begin(), target_.end(), traceback_pair_equal),
                  target_.end());
    for (const auto &pair : target_) {
      endpoints_.push_back(pair.i);
      endpoints_.push_back(pair.j);
    }
    std::sort(endpoints_.begin(), endpoints_.end());
    tree_.resize(4 * target_.size());
    if (!target_.empty()) build(0, 0, target_.size());
  }

  struct Stats {
    std::size_t interval_queries = 0;
    std::size_t range_node_visits = 0;
  };
  const Stats &stats() const { return stats_; }

  template <class Context>
  bool compatible_with_scaffold(const Context &ctx, cand_pos_t n) const {
    if (std::adjacent_find(endpoints_.begin(), endpoints_.end()) != endpoints_.end())
      return false;
    for (const auto &pair : target_) {
      if (pair.i < 1 || pair.j > n || pair.i >= pair.j ||
          !ctx.is_unpaired(pair.i) || !ctx.is_unpaired(pair.j)) return false;
    }
    return true;
  }

  // Complete-output conditioning on the shared paired/unpaired emission
  // surface. Every generated pair must be requested, and no requested pair
  // endpoint may be emitted unpaired. For valid complete derivations, the
  // per-base emission partition therefore forces the entire target output.
  // Unlike a bounding-span partition, this also handles VP crossing a
  // scaffold component whose ordinary interior is consumed elsewhere.
  template <class Context>
  bool matches_local_emissions(const deductions::Deduction &deduction,
                               const Context &ctx) const {
    std::vector<TracebackPair> emitted;
    append_deduction_generated_pairs(deduction, ctx, &emitted);
    for (const auto &pair : emitted) {
      if (!std::binary_search(target_.begin(), target_.end(), pair, before))
        return false;
    }
    const auto unpaired = profiles::unpaired_emission(
        deduction, [&](cand_pos_t position) { return ctx.is_unpaired(position); });
    for (std::size_t k = 0; k < unpaired.count; ++k) {
      const auto &run = unpaired.runs[k];
      const auto first = std::lower_bound(endpoints_.begin(), endpoints_.end(), run.lo);
      if (first != endpoints_.end() && *first <= run.hi) return false;
    }
    return true;
  }

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
  std::vector<cand_pos_t> endpoints_;
  std::vector<std::vector<cand_pos_t>> tree_;
  mutable Stats stats_;
};

}  // namespace cparty::scfg2::engine

#endif

#ifndef CPARTY_SCFG2_BEAM_BOUNDED_RANKED_BUCKET_HH
#define CPARTY_SCFG2_BEAM_BOUNDED_RANKED_BUCKET_HH

#include "scfg2/beam/memory_observation.hh"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct BoundedRankedBucketStats {
  std::size_t update_comparisons = 0;
  std::size_t view_comparisons = 0;
  std::size_t view_builds = 0;
  std::size_t view_cells = 0;
};

// Rank(a,b) means a is preferred to b. The heap root is the worst retained
// entry, so each fixed-width update takes O(log(b+2)), including rejection.
// Equivalent primary ranks are resolved by arrival order (earlier first).
// A ranked vector is built only on demand, O(b log(b+2)), then cached until
// the retained set changes. Views must not be used across mutations.
// Changing the width cannot restore earlier discarded entries; shrinking
// may evict multiple old entries. KeepAll has no finite-b complexity bound.
template <class Cell, class Rank>
class BoundedRankedBucket {
 public:
  void insert(Cell cell, std::size_t width) {
    auto compare = [&](const Entry &a, const Entry &b) {
      ++stats_.update_comparisons;
      return before(a, b);
    };
    while (heap_.size() > width) {
      std::pop_heap(heap_.begin(), heap_.end(), compare);
      heap_.pop_back();
      view_valid_ = false;
    }
    Entry entry{std::move(cell), next_arrival_++};
    if (!width) return;
    if (heap_.size() < width) {
      heap_.push_back(std::move(entry));
      std::push_heap(heap_.begin(), heap_.end(), compare);
      view_valid_ = false;
    } else if (compare(entry, heap_.front())) {
      std::pop_heap(heap_.begin(), heap_.end(), compare);
      heap_.back() = std::move(entry);
      std::push_heap(heap_.begin(), heap_.end(), compare);
      view_valid_ = false;
    }
  }

  const std::vector<Cell> &ranked() const {
    if (view_valid_) return view_;
    std::vector<std::size_t> order(heap_.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      ++stats_.view_comparisons;
      return before(heap_[a], heap_[b]);
    });
    view_.clear();
    view_.reserve(heap_.size());
    for (std::size_t idx : order) view_.push_back(heap_[idx].cell);
    ++stats_.view_builds;
    stats_.view_cells += heap_.size();
    view_valid_ = true;
    return view_;
  }

  std::size_t size() const { return heap_.size(); }
  // Observation of the exact cutoff, including the heap's FIFO tie break.
  // The pointer is valid only until the next mutation of this bucket.
  const Cell *worst() const {
    return heap_.empty() ? nullptr : &heap_.front().cell;
  }
  const BoundedRankedBucketStats &stats() const { return stats_; }

  ContainerMemory memory_snapshot(const char *name) const {
    auto out = vector_memory(name, heap_);
    out.add(vector_memory("", view_)); // Never build/refresh the cached view.
    return out;
  }

 private:
  struct Entry {
    Cell cell;
    std::size_t arrival;
  };
  bool before(const Entry &a, const Entry &b) const {
    if (rank_(a.cell, b.cell)) return true;
    if (rank_(b.cell, a.cell)) return false;
    return a.arrival < b.arrival;
  }

  Rank rank_{};
  std::vector<Entry> heap_;
  std::size_t next_arrival_ = 0;
  mutable std::vector<Cell> view_;
  mutable bool view_valid_ = false;
  mutable BoundedRankedBucketStats stats_{};
};

}  // namespace cparty::scfg2::beam
#endif

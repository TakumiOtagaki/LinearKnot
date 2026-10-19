#ifndef CPARTY_SCFG2_BEAM_BOUNDED_PREFIX_PENDING_HH
#define CPARTY_SCFG2_BEAM_BOUNDED_PREFIX_PENDING_HH

#include "scfg2/beam/bounded_ranked_bucket.hh"
#include "scfg2/beam/prefix_agenda_storage.hh"
#include "scfg2/types.hh"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

// A coordinate-indexed agenda for immutable, uniquely registered candidates.
// Each activation group keeps its original top-width order. Bucket presence
// includes arrivals discarded at width zero so take() consumes their column.
template <class Record, class Rank>
class BoundedPrefixPending {
 public:
  struct Taken {
    std::vector<Record> records;
  };

  void insert(cand_pos_t column, Record record, std::size_t width) {
    if (static_cast<std::size_t>(column) >= columns_.size()) columns_.resize(column + 1);
    auto &slot = columns_[column];
    ++slot.candidates;
    const auto before = slot.heap.stats().update_comparisons;
    slot.heap.insert(std::move(record), width);
    update_comparisons_ += slot.heap.stats().update_comparisons - before;
  }

  Taken take(cand_pos_t column) {
    if (column < 0 || static_cast<std::size_t>(column) >= columns_.size()) return {};
    auto &slot = columns_[column];
    if (!slot.candidates) return {};
    Taken out{slot.heap.ranked()};
    view_comparisons_ += slot.heap.stats().view_comparisons;
    slot = Slot{}; // release both heap and its cached view
    return out;
  }

  void append_storage(PrefixAgendaStorage &out) const {
    for (const auto &slot : columns_) {
      out.pending_records += slot.heap.size();
      out.largest_pending_bucket = std::max(out.largest_pending_bucket, slot.heap.size());
    }
  }
  std::size_t comparisons() const { return update_comparisons_ + view_comparisons_; }

  ContainerMemory memory_snapshot(const char *name) const {
    auto out = vector_memory(name, columns_);
    for (const auto &slot : columns_) out.add(slot.heap.memory_snapshot(""));
    return out;
  }

 private:
  struct Slot {
    BoundedRankedBucket<Record, Rank> heap;
    std::size_t candidates = 0;
  };
  std::vector<Slot> columns_;
  std::size_t update_comparisons_ = 0, view_comparisons_ = 0;
};
} // namespace cparty::scfg2::beam
#endif

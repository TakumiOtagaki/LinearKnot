#ifndef CPARTY_SCFG2_BEAM_RIGHT_BOUNDARY_ITEM_FRONTIER_HH
#define CPARTY_SCFG2_BEAM_RIGHT_BOUNDARY_ITEM_FRONTIER_HH

#include "scfg2/deductions/item_key.hh"
#include "scfg2/beam/memory_observation.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

using ItemKeySet = std::unordered_set<deductions::ItemKey,
                                      deductions::ItemKeyHash,
                                      deductions::ItemKeyEq>;

struct RankedRightBoundaryItem {
  deductions::ItemKey item{};
  double priority = 0.0;
};

inline auto right_boundary_item_tie_key(const deductions::ItemKey &item) {
  return std::make_tuple(static_cast<int>(item.nonterminal),
                         item.i,
                         item.j,
                         item.ip,
                         item.jp);
}

inline bool right_boundary_item_rank_before(
    const RankedRightBoundaryItem &lhs,
    const RankedRightBoundaryItem &rhs) {
  if (lhs.priority != rhs.priority) return lhs.priority > rhs.priority;
  return right_boundary_item_tie_key(lhs.item) <
         right_boundary_item_tie_key(rhs.item);
}

// Semantics-free right-boundary bucket. Grammar-family wrappers decide which
// completed items enter one instance and what retaining an item triggers.
// Keeping separate instances is the context-specific view/beam contract.
class RightBoundaryTopWidthItemFrontier {
 public:
  static RightBoundaryTopWidthItemFrontier keep_all() {
    return RightBoundaryTopWidthItemFrontier(
        std::numeric_limits<std::size_t>::max());
  }

  explicit RightBoundaryTopWidthItemFrontier(std::size_t width)
      : width_(width) {}
  RightBoundaryTopWidthItemFrontier(const RightBoundaryTopWidthItemFrontier &) = default;
  RightBoundaryTopWidthItemFrontier(RightBoundaryTopWidthItemFrontier &&) noexcept = default;
  RightBoundaryTopWidthItemFrontier &operator=(RightBoundaryTopWidthItemFrontier &&) noexcept = default;
  RightBoundaryTopWidthItemFrontier &operator=(const RightBoundaryTopWidthItemFrontier &other) {
    if (this != &other) {
      RightBoundaryTopWidthItemFrontier copy(other);
      *this = std::move(copy);
    }
    return *this;
  }

  // Explicit private numeric opt-in. One primary family and one pending
  // column use start-indexed authoritative storage. The kept callback order
  // remains the total rank order; the pruned tail's order is unspecified in
  // this mode. Callbacks may observe counts or erase numeric chart entries,
  // but must not submit/flush this frontier recursively. Default callers keep
  // the original map and callback ordering, including unusual keys/priorities.
  // Every unsupported key, simultaneous second column or NaN restores that
  // map by replaying first insertions, then permanently leaves direct mode.
  // Infinities retain a total rank order and are supported directly.
  bool enable_numeric_canonical_storage(grammar::NonTerminal family, cand_pos_t n) {
    if (!pending_by_column_.empty() || !direct_candidates_.empty() || n < 1 ||
        static_cast<std::size_t>(family) >= static_cast<std::size_t>(grammar::NonTerminal::BE))
      return false;
    try {
      std::vector<std::size_t> positions(static_cast<std::size_t>(n) + 1, 0);
      std::vector<RankedRightBoundaryItem> candidates;
      candidates.reserve(static_cast<std::size_t>(n));
      direct_positions_.swap(positions);
      direct_candidates_.swap(candidates);
      direct_family_ = family;
      direct_length_ = n;
      direct_column_ = -1;
      direct_enabled_ = true;
      return true;
    } catch (const std::bad_alloc &) {
      return false;
    } catch (const std::length_error &) {
      return false;
    }
  }
  bool numeric_canonical_storage_enabled() const { return direct_enabled_; }

  // Items that are always retained when live, in addition to the top
  // `width_` of the others (the scaffold-only derivation; see
  // scaffold_only_derivation_items).
  void set_forced(std::shared_ptr<const ItemKeySet> forced) {
    forced_ = std::move(forced);
  }

  void observe(const deductions::ItemKey &item, double priority) {
    if (direct_enabled_) {
      if (!direct_key(item) || std::isnan(priority) ||
          (!direct_candidates_.empty() && direct_column_ != item.j)) {
        restore_map();
      } else {
        if (direct_candidates_.empty()) direct_column_ = item.j;
        auto &position = direct_positions_[static_cast<std::size_t>(item.i)];
        if (position) {
          auto &old = direct_candidates_[position - 1];
          if (priority > old.priority) old.priority = priority;
          return;
        }
        // Normal runs reserve n entries once. Copies may lose that spare
        // capacity; a failed optional growth falls back without losing data.
        try {
          direct_candidates_.push_back({item, priority});
        } catch (const std::bad_alloc &) {
          restore_map();
          observe(item, priority);
          return;
        } catch (const std::length_error &) {
          restore_map();
          observe(item, priority);
          return;
        }
        position = direct_candidates_.size();
        ++candidate_count_;
        return;
      }
    }
    auto &column = pending_by_column_[item.j];
    const auto [found, inserted] = column.try_emplace(
        item, RankedRightBoundaryItem{item, priority});
    if (inserted) {
      ++candidate_count_;
    } else if (priority > found->second.priority) {
      found->second.priority = priority;
    }
  }

  template <class RetainedFn>
  void flush_before(cand_pos_t target_j, RetainedFn &&on_retained) {
    if (direct_enabled_ && !direct_candidates_.empty() && direct_column_ < target_j)
      flush_column_partition(direct_column_, on_retained, [](const deductions::ItemKey &) {});
    auto bucket = pending_by_column_.begin();
    while (bucket != pending_by_column_.end() && bucket->first < target_j) {
      auto ranked = ranked_candidates(bucket->second);
      const std::size_t kept = select_top(ranked);
      for (std::size_t idx = 0; idx < kept; ++idx) {
        on_retained(ranked[idx].item);
      }
      selected_count_ += kept;
      maximum_selected_per_column_ =
          std::max(maximum_selected_per_column_, kept - last_forced_);
      bucket = pending_by_column_.erase(bucket);
    }
  }

  template <class RetainedFn>
  void flush_column(cand_pos_t j, RetainedFn &&on_retained) {
    flush_column_partition(
        j,
        std::forward<RetainedFn>(on_retained),
        [](const deductions::ItemKey &) {});
  }

  template <class RetainedFn, class PrunedFn>
  void flush_column_partition(cand_pos_t j,
                              RetainedFn &&on_retained,
                              PrunedFn &&on_pruned) {
    const bool direct = direct_enabled_ && !direct_candidates_.empty() && direct_column_ == j;
    const auto found = pending_by_column_.find(j);
    if (!direct && found == pending_by_column_.end()) return;
    auto ranked = direct ? direct_candidates_ : ranked_candidates(found->second);
    const std::size_t kept = select_top(ranked);
    for (std::size_t idx = 0; idx < kept; ++idx) {
      on_retained(ranked[idx].item);
    }
    for (std::size_t idx = kept; idx < ranked.size(); ++idx) {
      on_pruned(ranked[idx].item);
    }
    selected_count_ += kept;
    maximum_selected_per_column_ =
        std::max(maximum_selected_per_column_, kept - last_forced_);
    if (direct) clear_direct_column();
    else pending_by_column_.erase(found);
  }

  std::size_t candidate_count() const { return candidate_count_; }
  std::size_t selected_count() const { return selected_count_; }
  std::size_t pruned_count() const {
    return candidate_count_ - selected_count_ - pending_candidate_count();
  }
  bool respects_width() const {
    return maximum_selected_per_column_ <= width_;
  }
  std::size_t selection_comparisons() const { return selection_comparisons_; }

  ContainerMemory memory_snapshot(const char *name) const {
    ContainerMemory out;
    out.name = name;
    out.value_bytes = pending_by_column_.size() * sizeof(decltype(pending_by_column_)::value_type);
    for (const auto &[column, candidates] : pending_by_column_) {
      (void)column;
      out.add(hash_memory("", candidates));
    }
    out.add(vector_memory("", direct_positions_));
    out.add(vector_memory("", direct_candidates_));
    return out;
  }

 private:
  using ColumnCandidates =
      std::unordered_map<deductions::ItemKey,
                         RankedRightBoundaryItem,
                         deductions::ItemKeyHash,
                         deductions::ItemKeyEq>;

  // Moves the top `width_` candidates to the front in rank order and returns
  // how many were kept.  The rank is a strict total order (priority, then
  // the item key). Heap selection retains the same prefix as a full sort
  // with an O(C log(b+2)) bound even when C is much larger than b^2.
  // The pruned tail is left unordered.
  std::size_t select_top(std::vector<RankedRightBoundaryItem> &ranked) {
    const auto before = [&](const auto &a, const auto &b) {
      ++selection_comparisons_;
      return right_boundary_item_rank_before(a, b);
    };
    // Forced items go first and do not use a width slot.
    std::size_t forced = 0;
    if (forced_ && !forced_->empty()) {
      const auto split = std::partition(
          ranked.begin(), ranked.end(),
          [&](const RankedRightBoundaryItem &c) { return forced_->count(c.item) > 0; });
      forced = static_cast<std::size_t>(split - ranked.begin());
    }
    const auto rest_begin = ranked.begin() + static_cast<std::ptrdiff_t>(forced);
    const std::size_t rest = ranked.size() - forced;
    const std::size_t kept = forced + std::min(width_, rest);
    last_forced_ = forced;
    if (kept < ranked.size() && kept > forced) {
      std::partial_sort(rest_begin,
                        ranked.begin() + static_cast<std::ptrdiff_t>(kept),
                        ranked.end(), before);
    }
    std::sort(ranked.begin(),
              ranked.begin() + static_cast<std::ptrdiff_t>(kept), before);
    return kept;
  }

  static std::vector<RankedRightBoundaryItem> ranked_candidates(
      const ColumnCandidates &column) {
    std::vector<RankedRightBoundaryItem> ranked;
    ranked.reserve(column.size());
    for (const auto &[item, candidate] : column) {
      (void)item;
      ranked.push_back(candidate);
    }
    return ranked;
  }

  std::size_t pending_candidate_count() const {
    std::size_t count = direct_candidates_.size();
    for (const auto &[j, candidates] : pending_by_column_) {
      (void)j;
      count += candidates.size();
    }
    return count;
  }

  bool direct_key(const deductions::ItemKey &item) const {
    return !direct_positions_.empty() && item.nonterminal == direct_family_ &&
        item.ip == -1 && item.jp == -1 && item.i >= 1 &&
        item.i <= item.j && item.j <= direct_length_;
  }
  void clear_direct_column() {
    for (const auto &candidate : direct_candidates_)
      direct_positions_[static_cast<std::size_t>(candidate.item.i)] = 0;
    direct_candidates_.clear();
    direct_column_ = -1;
  }
  void restore_map() {
    // Nothing has been erased from this pending column. Replaying the first
    // arrival order recreates its original bucket/iteration layout; duplicate
    // priority improvements never reordered old map nodes. Commit only after
    // complete replay, so a failed mandatory map allocation preserves state.
    std::map<cand_pos_t, ColumnCandidates> restored;
    if (!direct_candidates_.empty()) {
      auto &column = restored[direct_column_];
      for (const auto &candidate : direct_candidates_)
        column.try_emplace(candidate.item, candidate);
    }
    pending_by_column_.swap(restored);
    direct_enabled_ = false;
    direct_column_ = -1;
    std::vector<std::size_t>{}.swap(direct_positions_);
    std::vector<RankedRightBoundaryItem>{}.swap(direct_candidates_);
  }

  std::size_t width_ = 0;
  std::shared_ptr<const ItemKeySet> forced_;
  std::size_t last_forced_ = 0;
  std::size_t selection_comparisons_ = 0;
  std::map<cand_pos_t, ColumnCandidates> pending_by_column_;
  bool direct_enabled_ = false;
  grammar::NonTerminal direct_family_ = grammar::NonTerminal::W;
  cand_pos_t direct_length_ = 0, direct_column_ = -1;
  std::vector<std::size_t> direct_positions_;
  std::vector<RankedRightBoundaryItem> direct_candidates_;
  std::size_t candidate_count_ = 0;
  std::size_t selected_count_ = 0;
  std::size_t maximum_selected_per_column_ = 0;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_RIGHT_BOUNDARY_ITEM_FRONTIER_HH

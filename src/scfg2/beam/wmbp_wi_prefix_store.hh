#ifndef CPARTY_SCFG2_BEAM_WMBP_WI_PREFIX_STORE_HH
#define CPARTY_SCFG2_BEAM_WMBP_WI_PREFIX_STORE_HH

#include "scfg2/deductions/item_key.hh"
#include "scfg2/beam/bounded_ranked_bucket.hh"

#include <cstddef>
#include <optional>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct WmbpWiHistoryBoundary {
  double priority;
  cand_pos_t right;
};

template <class Value>
struct WmbpWiPrefixCell {
  deductions::ItemKey be{};
  deductions::ItemKey wi{};
  bool has_wi = false;
  cand_pos_t vp_i = -1;
  Value factor{};
};

struct WmbpWiPrefixStoreStats {
  std::size_t wi_history_insertions = 0;
  std::size_t be_wi_probes = 0;
  std::size_t empty_probes = 0;
  std::size_t prefix_candidates = 0;
  std::size_t bucket_update_comparisons = 0;
  std::size_t bucket_view_comparisons = 0;
  std::size_t bucket_view_builds = 0;
  std::size_t bucket_view_cells = 0;
};

template <class Semiring>
class WmbpWiPrefixStore {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart>
  void commit_wi_column(cand_pos_t j,
                        const std::vector<cand_pos_t> &starts,
                        const Chart &chart,
                        std::size_t width) {
    for (const cand_pos_t i : starts) {
      const deductions::ItemKey item{
          grammar::NonTerminal::WI, i, j, -1, -1};
      const Value value = chart.get(item);
      if (Semiring::is_zero(value)) continue;
      auto &cells = wi_by_left_[i];
      ++stats_.wi_history_insertions;
      insert_bucket(cells, WiCell{item, value}, width);
    }
  }

  template <class Chart>
  void commit_be_column(const std::vector<deductions::ItemKey> &be_items,
                        const Chart &chart,
                        std::size_t width) {
    for (const auto &be : be_items) {
      const Value be_value = chart.get(be);
      if (Semiring::is_zero(be_value)) continue;
      ++stats_.empty_probes;
      add_prefix(WmbpWiPrefixCell<Value>{
          be, {}, false, be.ip + 1, be_value}, width);
      const auto found = wi_by_left_.find(be.ip + 1);
      if (found == wi_by_left_.end()) continue;
      for (const auto &wi : ranked_view(found->second)) {
        ++stats_.be_wi_probes;
        // Positive WI-mode sources complete BE strictly after the right-VP
        // start. The targeted timing oracle guards this online-order contract.
        if (wi.item.j + 1 >= be.j) continue;
        add_prefix(WmbpWiPrefixCell<Value>{
            be,
            wi.item,
            true,
            wi.item.j + 1,
            Semiring::times(be_value, wi.value)},
            width);
      }
    }
  }

  const std::vector<WmbpWiPrefixCell<Value>> &retained(cand_pos_t vp_i) const {
    const auto found = prefixes_by_vp_start_.find(vp_i);
    if (found == prefixes_by_vp_start_.end()) {
      static const std::vector<WmbpWiPrefixCell<Value>> empty;
      return empty;
    }
    return ranked_view(found->second);
  }

  const WmbpWiPrefixStoreStats &stats() const { return stats_; }

  // Constant-time observations; these neither build ranked views nor save
  // numerical expressions/branch lists for outside.
  std::optional<WmbpWiHistoryBoundary> wi_history_boundary(cand_pos_t i) const {
    const auto found = wi_by_left_.find(i);
    if (found == wi_by_left_.end()) return std::nullopt;
    const auto *last = found->second.worst();
    if (!last) return std::nullopt;
    return WmbpWiHistoryBoundary{Semiring::priority(last->value), last->item.j};
  }
  const WmbpWiPrefixCell<Value> *prefix_boundary(cand_pos_t vp_i) const {
    const auto found = prefixes_by_vp_start_.find(vp_i);
    return found == prefixes_by_vp_start_.end() ? nullptr : found->second.worst();
  }

 private:
  struct WiCell {
    deductions::ItemKey item{};
    Value value{};
  };

  struct WiRank {
    bool operator()(const WiCell &lhs, const WiCell &rhs) const {
      const double lhs_priority = Semiring::priority(lhs.value);
      const double rhs_priority = Semiring::priority(rhs.value);
      if (lhs_priority != rhs_priority) return lhs_priority > rhs_priority;
      return lhs.item.j < rhs.item.j;
    }
  };
  struct PrefixRank {
    bool operator()(const WmbpWiPrefixCell<Value> &lhs,
                    const WmbpWiPrefixCell<Value> &rhs) const {
      const double lhs_priority = Semiring::priority(lhs.factor);
      const double rhs_priority = Semiring::priority(rhs.factor);
      if (lhs_priority != rhs_priority) return lhs_priority > rhs_priority;
      if (lhs.be.i != rhs.be.i) return lhs.be.i < rhs.be.i;
      return lhs.be.j < rhs.be.j;
    }
  };

  template <class Bucket, class Cell>
  void insert_bucket(Bucket &bucket, Cell cell, std::size_t width) {
    const std::size_t before = bucket.stats().update_comparisons;
    bucket.insert(std::move(cell), width);
    stats_.bucket_update_comparisons += bucket.stats().update_comparisons - before;
  }

  template <class Bucket>
  const auto &ranked_view(const Bucket &bucket) const {
    const auto before = bucket.stats();
    const auto &cells = bucket.ranked();
    const auto &after = bucket.stats();
    stats_.bucket_view_comparisons += after.view_comparisons - before.view_comparisons;
    stats_.bucket_view_builds += after.view_builds - before.view_builds;
    stats_.bucket_view_cells += after.view_cells - before.view_cells;
    return cells;
  }

  void add_prefix(WmbpWiPrefixCell<Value> cell, std::size_t width) {
    auto &cells = prefixes_by_vp_start_[cell.vp_i];
    insert_bucket(cells, std::move(cell), width);
    ++stats_.prefix_candidates;
  }

  std::unordered_map<cand_pos_t, BoundedRankedBucket<WiCell, WiRank>> wi_by_left_;
  std::unordered_map<cand_pos_t, BoundedRankedBucket<WmbpWiPrefixCell<Value>, PrefixRank>>
      prefixes_by_vp_start_;
  mutable WmbpWiPrefixStoreStats stats_{};
};

}  // namespace cparty::scfg2::beam

#endif

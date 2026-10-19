#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_LEFT_CACHE_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_LEFT_CACHE_HH

#include "scfg2/beam/vp_case3_tail_attachment_segments.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase3LeftKey {
  VpCase3TailKey tail{};
  cand_pos_t i = -1;
  cand_pos_t p = -1;
};

inline constexpr bool operator==(const VpCase3LeftKey &lhs,
                                 const VpCase3LeftKey &rhs) {
  return lhs.tail == rhs.tail && lhs.i == rhs.i && lhs.p == rhs.p;
}

struct VpCase3LeftKeyHash {
  std::size_t operator()(const VpCase3LeftKey &key) const noexcept {
    std::size_t out = VpCase3TailKeyHash{}(key.tail);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.i + 13);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.p + 17);
    return out;
  }
};

template <class Value>
struct BasicVpCase3LeftRecord {
  VpCase3LeftKey key{};
  Value factor{};
  unsigned carrier_mask = 0;
};

using VpCase3LeftRecord = BasicVpCase3LeftRecord<double>;

// Independent segment/column counts for the retained-child reference checks.
class ExactVpCase3LeftColumnPlan {
 public:
  explicit ExactVpCase3LeftColumnPlan(
      const std::vector<VpCase3TailAttachmentSegment> &segments) {
    std::unordered_set<cand_pos_t> columns;
    for (const auto &segment : segments) {
      columns.insert(segment.tail.j);
      geometry_count_ +=
          static_cast<std::size_t>(segment.i_end - segment.i_begin + 1);
    }
    column_count_ = columns.size();
  }

  std::size_t geometry_count() const { return geometry_count_; }
  std::size_t column_count() const { return column_count_; }

 private:
  std::size_t column_count_ = 0;
  std::size_t geometry_count_ = 0;
};

template <class Value>
class BasicExactVpCase3LeftCache {
 public:
  template <class AllowsKeyFn>
  std::size_t materialize_scored(
      const std::vector<BasicVpCase3LeftRecord<Value>> &scored,
      AllowsKeyFn &&allows_key) {
    std::size_t inserted = 0;
    for (const auto &record : scored) {
      if (!allows_key(record.key)) continue;
      auto sealed = record;
      sealed.factor = semirings::checkpoint_value(record.factor);
      const auto [found, did_insert] = records_.emplace(record.key, std::move(sealed));
      (void)found;
      if (!did_insert) continue;
      parent_starts_by_column_[record.key.tail.j].push_back(
          ParentStartRecord{
              record.key.i,
              record.carrier_mask});
      ++inserted;
    }
    return inserted;
  }

  template <class Semiring, class LocalWeight>
  std::optional<Value> reconstruct_source(
      const deductions::Deduction &deduction,
      LocalWeight local_weight) {
    if (vp_case3_partition(deduction) != VpCase3Partition::LeftNonempty) {
      return std::nullopt;
    }
    const VpCase3LeftKey key{
        VpCase3TailKey{deduction.parent.j,
                       deduction.split.q,
                       deduction.split.k,
                       deduction.split.l},
        deduction.parent.i,
        deduction.split.p};
    const auto record = records_.find(key);
    if (record == records_.end()) return std::nullopt;
    const unsigned carrier =
        vp_case3_carrier_bit(deduction.parent.nonterminal);
    if ((record->second.carrier_mask & carrier) == 0) return std::nullopt;
    ++attachment_count_;
    return Semiring::times(
        Semiring::from_weight(local_weight), record->second.factor);
  }

  std::optional<Value> reconstruct_source(
      const deductions::Deduction &deduction,
      double local_weight) {
    return reconstruct_source<semirings::SumProductSemiring>(
        deduction, local_weight);
  }

  std::size_t candidate_count() const { return records_.size(); }
  std::size_t materialization_count() const { return records_.size(); }
  std::size_t attachment_count() const { return attachment_count_; }
  const auto &records() const { return records_; }

  std::vector<cand_pos_t> parent_starts(
      cand_pos_t j, grammar::NonTerminal carrier) const {
    std::vector<cand_pos_t> out;
    const unsigned carrier_bit = vp_case3_carrier_bit(carrier);
    const auto column = parent_starts_by_column_.find(j);
    if (column == parent_starts_by_column_.end()) return out;
    for (const auto &record : column->second) {
      ++parent_start_probe_count_;
      if ((record.carrier_mask & carrier_bit) != 0) out.push_back(record.i);
    }
    std::sort(out.begin(), out.end(), std::greater<cand_pos_t>());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  std::vector<cand_pos_t> vp_parent_starts(cand_pos_t j) const {
    return parent_starts(j, grammar::NonTerminal::VP);
  }

  std::size_t parent_start_probe_count() const {
    return parent_start_probe_count_;
  }

 private:
  struct ParentStartRecord {
    cand_pos_t i = -1;
    unsigned carrier_mask = 0;
  };

  std::unordered_map<VpCase3LeftKey,
                     BasicVpCase3LeftRecord<Value>,
                     VpCase3LeftKeyHash>
      records_;
  std::unordered_map<cand_pos_t, std::vector<ParentStartRecord>>
      parent_starts_by_column_;
  std::size_t attachment_count_ = 0;
  mutable std::size_t parent_start_probe_count_ = 0;
};

using ExactVpCase3LeftCache = BasicExactVpCase3LeftCache<double>;

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_LEFT_CACHE_HH

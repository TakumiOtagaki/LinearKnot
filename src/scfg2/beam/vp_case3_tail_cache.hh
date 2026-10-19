#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_CACHE_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_CACHE_HH

#include "scfg2/beam/vp_case3_factor.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase3TailKey {
  cand_pos_t j = -1;
  cand_pos_t q = -1;
  cand_pos_t k = -1;
  cand_pos_t l = -1;
};

inline constexpr bool operator==(const VpCase3TailKey &lhs,
                                 const VpCase3TailKey &rhs) {
  return lhs.j == rhs.j && lhs.q == rhs.q && lhs.k == rhs.k &&
         lhs.l == rhs.l;
}

struct VpCase3TailKeyHash {
  std::size_t operator()(const VpCase3TailKey &key) const noexcept {
    std::size_t out = static_cast<std::size_t>(key.j + 3);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.q + 5);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.k + 7);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.l + 11);
    return out;
  }
};

struct VpCase3RightFactorKey {
  cand_pos_t j = -1;
  cand_pos_t l = -1;
};

inline constexpr bool operator==(const VpCase3RightFactorKey &lhs,
                                 const VpCase3RightFactorKey &rhs) {
  return lhs.j == rhs.j && lhs.l == rhs.l;
}

struct VpCase3RightFactorKeyHash {
  std::size_t operator()(const VpCase3RightFactorKey &key) const noexcept {
    return (static_cast<std::size_t>(key.j + 3) * 1315423911u) ^
           static_cast<std::size_t>(key.l + 5);
  }
};

struct VpCase3MiddleFactorKey {
  cand_pos_t q = -1;
  cand_pos_t k = -1;
};

inline constexpr bool operator==(const VpCase3MiddleFactorKey &lhs,
                                 const VpCase3MiddleFactorKey &rhs) {
  return lhs.q == rhs.q && lhs.k == rhs.k;
}

struct VpCase3MiddleFactorKeyHash {
  std::size_t operator()(const VpCase3MiddleFactorKey &key) const noexcept {
    return (static_cast<std::size_t>(key.q + 3) * 1315423911u) ^
           static_cast<std::size_t>(key.k + 5);
  }
};

// Only source identity is observed after reconstruction; numerical factors
// belong to the chart and the shared tail, not to each attachment record.
struct VpCase3SourceAttachment {
  deductions::ItemKey source_parent{};
  grammar::RuleSplit source_split{};
};

inline bool same_attachment_source(
    const VpCase3SourceAttachment &lhs,
    const VpCase3SourceAttachment &rhs) {
  return deductions::ItemKeyEq{}(lhs.source_parent, rhs.source_parent) &&
         lhs.source_split.k == rhs.source_split.k &&
         lhs.source_split.l == rhs.source_split.l &&
         lhs.source_split.p == rhs.source_split.p &&
         lhs.source_split.q == rhs.source_split.q;
}

// Exact/private cache contract. The cache shares only the algebraic
// middle-right product. Standalone reconstruction also records each complete
// source identity for VP_WI_CASE3 reference checks, separately from the values.
template <class Value>
class BasicExactVpCase3TailCache {
 public:
  template <class Semiring, class Chart>
  bool materialize_tail(const VpCase3TailKey &key, const Chart &chart) {
    if (tail_factors_.find(key) != tail_factors_.end()) return false;
    const VpCase3MiddleFactorKey middle_key{key.q, key.k};
    auto middle = middle_factors_.find(middle_key);
    if (middle == middle_factors_.end()) {
      middle = middle_factors_
                   .emplace(middle_key,
                            case3_wi_span_total<Semiring>(
                                chart, key.q + 1, key.k - 1))
                   .first;
    } else {
      ++middle_factor_reuse_count_;
    }
    const VpCase3RightFactorKey right_key{key.j, key.l};
    auto right = right_factors_.find(right_key);
    if (right == right_factors_.end()) {
      right = right_factors_
                  .emplace(right_key,
                           case3_wi_span_total<Semiring>(
                               chart, key.l + 1, key.j - 1))
                  .first;
    } else {
      ++right_factor_reuse_count_;
    }
    tail_factors_.emplace(
        key,
        semirings::checkpoint_value(
            Semiring::times(middle->second, right->second)));
    keys_by_column_[key.j].push_back(key);
    return true;
  }

  template <class Semiring, class Chart>
  std::optional<Value> reconstruct_materialized_left_source(
      const deductions::Deduction &deduction,
      const Chart &chart,
      double local_weight) {
    if (vp_case3_partition(deduction) != VpCase3Partition::LeftNonempty) {
      return std::nullopt;
    }
    const VpCase3TailKey key{deduction.parent.j,
                             deduction.split.q,
                             deduction.split.k,
                             deduction.split.l};
    const auto it = tail_factors_.find(key);
    if (it == tail_factors_.end()) return std::nullopt;
    ++reuse_count_;
    const Value left = case3_wi_span_total<Semiring>(
        chart, deduction.parent.i + 1, deduction.split.p - 1);
    const Value reconstructed = Semiring::times(
        Semiring::from_weight(local_weight),
        Semiring::times(left, it->second));
    source_attachments_[key].push_back(VpCase3SourceAttachment{
        deduction.parent,
        deduction.split});
    ++attachment_count_;
    return reconstructed;
  }

  template <class Semiring, class Chart>
  std::optional<Value> observe_left_source(
      const deductions::Deduction &deduction,
      const Chart &chart,
      double local_weight) {
    if (vp_case3_partition(deduction) != VpCase3Partition::LeftNonempty) {
      return std::nullopt;
    }

    const VpCase3TailKey key{deduction.parent.j,
                             deduction.split.q,
                             deduction.split.k,
                             deduction.split.l};
    const bool inserted = materialize_tail<Semiring>(key, chart);
    const auto reconstructed =
        reconstruct_materialized_left_source<Semiring>(
            deduction, chart, local_weight);
    if (inserted) --reuse_count_;
    return reconstructed;
  }

  template <class Chart>
  std::optional<Value> observe_left_source(
      const deductions::Deduction &deduction,
      const Chart &chart,
      double local_weight) {
    return observe_left_source<semirings::SumProductSemiring>(
        deduction, chart, local_weight);
  }

  std::size_t tail_count() const { return tail_factors_.size(); }
  std::size_t attachment_count() const { return attachment_count_; }
  std::size_t materialization_count() const { return tail_factors_.size(); }
  std::size_t reuse_count() const { return reuse_count_; }
  std::size_t right_factor_count() const { return right_factors_.size(); }
  std::size_t right_factor_materialization_count() const {
    return right_factors_.size();
  }
  std::size_t right_factor_reuse_count() const {
    return right_factor_reuse_count_;
  }
  std::size_t middle_factor_count() const { return middle_factors_.size(); }
  std::size_t middle_factor_materialization_count() const {
    return middle_factors_.size();
  }
  std::size_t middle_factor_reuse_count() const {
    return middle_factor_reuse_count_;
  }

  std::optional<Value> tail_factor(const VpCase3TailKey &key) const {
    const auto it = tail_factors_.find(key);
    if (it == tail_factors_.end()) return std::nullopt;
    return it->second;
  }

  // Reference observations: reading metadata never materializes a tail or
  // reconstructs a source. Production consumes only factors and column keys.
  // Reacquire returned views after any cache mutation.
  const std::vector<VpCase3SourceAttachment> &source_attachments(
      const VpCase3TailKey &key) const {
    const auto found = source_attachments_.find(key);
    if (found == source_attachments_.end()) {
      static const std::vector<VpCase3SourceAttachment> empty;
      return empty;
    }
    return found->second;
  }

  template <class Fn>
  void for_each_tail_key(Fn &&fn) const {
    for (const auto &[key, factor] : tail_factors_) {
      (void)factor;
      fn(key);
    }
  }

  std::size_t max_attachments_per_tail() const {
    std::size_t maximum = 0;
    for (const auto &[key, attachments] : source_attachments_) {
      (void)key;
      maximum = std::max(maximum, attachments.size());
    }
    return maximum;
  }

  bool full_source_attachments_are_unique() const {
    std::vector<const VpCase3SourceAttachment *> seen;
    for (const auto &[key, attachments] : source_attachments_) {
      (void)key;
      for (const auto &attachment : attachments) {
        for (const auto *prior : seen) {
          if (same_attachment_source(attachment, *prior)) return false;
        }
        seen.push_back(&attachment);
      }
    }
    return true;
  }

  const std::vector<VpCase3TailKey> &keys_for_column(cand_pos_t j) const {
    const auto found = keys_by_column_.find(j);
    if (found == keys_by_column_.end()) {
      static const std::vector<VpCase3TailKey> empty;
      return empty;
    }
    return found->second;
  }

 private:
  std::unordered_map<VpCase3TailKey, Value, VpCase3TailKeyHash> tail_factors_;
  // Populated only by standalone reconstruction/reference consumers.
  // Online inside uses the numeric map without allocating attachment records.
  std::unordered_map<VpCase3TailKey,
                     std::vector<VpCase3SourceAttachment>,
                     VpCase3TailKeyHash> source_attachments_;
  std::unordered_map<VpCase3RightFactorKey,
                     Value,
                     VpCase3RightFactorKeyHash>
      right_factors_;
  std::unordered_map<VpCase3MiddleFactorKey,
                     Value,
                     VpCase3MiddleFactorKeyHash>
      middle_factors_;
  std::unordered_map<cand_pos_t, std::vector<VpCase3TailKey>>
      keys_by_column_;
  std::size_t attachment_count_ = 0;
  std::size_t reuse_count_ = 0;
  std::size_t right_factor_reuse_count_ = 0;
  std::size_t middle_factor_reuse_count_ = 0;
};

using ExactVpCase3TailCache = BasicExactVpCase3TailCache<double>;

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_CACHE_HH

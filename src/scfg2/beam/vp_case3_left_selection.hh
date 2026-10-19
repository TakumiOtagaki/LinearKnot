#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_LEFT_SELECTION_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_LEFT_SELECTION_HH

#include "scfg2/beam/vp_case3_left_cache.hh"

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

struct RankedVpCase3LeftCandidate {
  VpCase3LeftKey key{};
  double priority = 0.0;
};

inline auto vp_case3_left_tie_key(const VpCase3LeftKey &key) {
  return std::make_tuple(key.i,
                         key.p,
                         key.tail.q,
                         key.tail.k,
                         key.tail.l,
                         key.tail.j);
}

inline bool vp_case3_left_rank_before(
    const RankedVpCase3LeftCandidate &lhs,
    const RankedVpCase3LeftCandidate &rhs) {
  if (lhs.priority != rhs.priority) return lhs.priority > rhs.priority;
  return vp_case3_left_tie_key(lhs.key) < vp_case3_left_tie_key(rhs.key);
}

inline VpCase3LeftKey vp_case3_left_source_key(
    const deductions::Deduction &deduction) {
  return VpCase3LeftKey{
      VpCase3TailKey{deduction.parent.j,
                     deduction.split.q,
                     deduction.split.k,
                     deduction.split.l},
      deduction.parent.i,
      deduction.split.p};
}

template <class Semiring, class Value>
inline std::vector<RankedVpCase3LeftCandidate>
rank_vp_case3_left_column(
    cand_pos_t j,
    const std::unordered_map<VpCase3LeftKey,
                             BasicVpCase3LeftRecord<Value>,
                             VpCase3LeftKeyHash> &records) {
  std::vector<RankedVpCase3LeftCandidate> ranked;
  for (const auto &[key, record] : records) {
    if (key.tail.j == j) {
      ranked.push_back({key, Semiring::priority(record.factor)});
    }
  }
  std::sort(ranked.begin(), ranked.end(), vp_case3_left_rank_before);
  return ranked;
}

inline std::vector<RankedVpCase3LeftCandidate> rank_vp_case3_left_column(
    cand_pos_t j,
    const std::unordered_map<VpCase3LeftKey,
                             VpCase3LeftRecord,
                             VpCase3LeftKeyHash> &records) {
  return rank_vp_case3_left_column<semirings::SumProductSemiring>(j, records);
}

template <class Semiring, class Value>
inline std::vector<RankedVpCase3LeftCandidate>
rank_vp_case3_left_column(
    const std::vector<BasicVpCase3LeftRecord<Value>> &records) {
  std::vector<RankedVpCase3LeftCandidate> ranked;
  ranked.reserve(records.size());
  for (const auto &record : records) {
    ranked.push_back({record.key, Semiring::priority(record.factor)});
  }
  std::sort(ranked.begin(), ranked.end(), vp_case3_left_rank_before);
  return ranked;
}

class KeepAllVpCase3LeftSelection {
 public:
  void commit_column(
      cand_pos_t j,
      const std::unordered_map<VpCase3LeftKey,
                               VpCase3LeftRecord,
                               VpCase3LeftKeyHash> &records) {
    auto ranked = rank_vp_case3_left_column(j, records);
    if (ranked.empty()) return;
    for (const auto &candidate : ranked) selected_.insert(candidate.key);
    ranked_columns_[j] = std::move(ranked);
  }

  bool allows(const deductions::Deduction &deduction) const {
    if (vp_case3_partition(deduction) != VpCase3Partition::LeftNonempty) {
      return true;
    }
    const VpCase3LeftKey key = vp_case3_left_source_key(deduction);
    return selected_.find(key) != selected_.end();
  }

  std::size_t selected_count() const { return selected_.size(); }
  std::size_t column_count() const { return ranked_columns_.size(); }

  bool columns_are_ranked() const {
    for (const auto &[j, ranked] : ranked_columns_) {
      (void)j;
      if (!std::is_sorted(
              ranked.begin(), ranked.end(), vp_case3_left_rank_before)) {
        return false;
      }
    }
    return true;
  }

 private:
  std::unordered_set<VpCase3LeftKey, VpCase3LeftKeyHash> selected_;
  std::unordered_map<cand_pos_t, std::vector<RankedVpCase3LeftCandidate>>
      ranked_columns_;
};

class TopWidthVpCase3LeftSelection {
 public:
  explicit TopWidthVpCase3LeftSelection(std::size_t width) : width_(width) {}

  template <class Semiring, class Value>
  void commit_column(cand_pos_t j,
                     const std::vector<BasicVpCase3LeftRecord<Value>> &records) {
    auto ranked = rank_vp_case3_left_column<Semiring>(records);
    if (ranked.empty()) return;
    const std::size_t kept = std::min(width_, ranked.size());
    for (std::size_t idx = 0; idx < kept; ++idx) {
      selected_.insert(ranked[idx].key);
    }
    selected_per_column_[j] = kept;
    candidate_per_column_[j] = ranked.size();
  }

  void commit_column(cand_pos_t j,
                     const std::vector<VpCase3LeftRecord> &records) {
    commit_column<semirings::SumProductSemiring>(j, records);
  }

  bool allows(const VpCase3LeftKey &key) const {
    return selected_.find(key) != selected_.end();
  }

  bool allows(const deductions::Deduction &deduction) const {
    if (vp_case3_partition(deduction) != VpCase3Partition::LeftNonempty) {
      return true;
    }
    return allows(vp_case3_left_source_key(deduction));
  }

  std::size_t selected_count() const { return selected_.size(); }
  std::size_t candidate_count() const {
    std::size_t count = 0;
    for (const auto &[j, value] : candidate_per_column_) {
      (void)j;
      count += value;
    }
    return count;
  }
  std::size_t pruned_count() const { return candidate_count() - selected_count(); }
  bool respects_width() const {
    for (const auto &[j, count] : selected_per_column_) {
      (void)j;
      if (count > width_) return false;
    }
    return true;
  }

 private:
  std::size_t width_ = 0;
  std::unordered_set<VpCase3LeftKey, VpCase3LeftKeyHash> selected_;
  std::unordered_map<cand_pos_t, std::size_t> selected_per_column_;
  std::unordered_map<cand_pos_t, std::size_t> candidate_per_column_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_LEFT_SELECTION_HH

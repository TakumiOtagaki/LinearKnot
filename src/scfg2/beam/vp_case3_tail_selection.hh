#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_SELECTION_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_SELECTION_HH

#include "scfg2/beam/vp_case3_tail_cache.hh"

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

struct RankedVpCase3TailCandidate {
  VpCase3TailKey key{};
  double priority = 0.0;
};

inline auto vp_case3_tail_tie_key(const VpCase3TailKey &key) {
  return std::make_tuple(key.q, key.k, key.l, key.j);
}

inline bool vp_case3_tail_rank_before(
    const RankedVpCase3TailCandidate &lhs,
    const RankedVpCase3TailCandidate &rhs) {
  if (lhs.priority != rhs.priority) return lhs.priority > rhs.priority;
  return vp_case3_tail_tie_key(lhs.key) < vp_case3_tail_tie_key(rhs.key);
}

class TopWidthVpCase3TailSelection {
 public:
  explicit TopWidthVpCase3TailSelection(std::size_t width) : width_(width) {}

  template <class Semiring, class Value>
  void commit_column(cand_pos_t j,
                     const BasicExactVpCase3TailCache<Value> &cache) {
    std::vector<RankedVpCase3TailCandidate> ranked;
    for (const auto &key : cache.keys_for_column(j)) {
      ++record_probe_count_;
      const auto factor = cache.tail_factor(key);
      if (factor.has_value()) {
        ranked.push_back({key, Semiring::priority(*factor)});
      }
    }
    std::sort(ranked.begin(), ranked.end(), vp_case3_tail_rank_before);
    const std::size_t kept = std::min(width_, ranked.size());
    for (std::size_t idx = 0; idx < kept; ++idx) {
      selected_.insert(ranked[idx].key);
    }
    candidate_per_column_[j] = ranked.size();
    selected_per_column_[j] = kept;
  }

  void commit_column(cand_pos_t j, const ExactVpCase3TailCache &cache) {
    commit_column<semirings::SumProductSemiring>(j, cache);
  }

  bool allows(const VpCase3TailKey &key) const {
    return selected_.find(key) != selected_.end();
  }

  std::size_t selected_count() const { return selected_.size(); }
  std::size_t candidate_count() const {
    std::size_t count = 0;
    for (const auto &[j, candidates] : candidate_per_column_) {
      (void)j;
      count += candidates;
    }
    return count;
  }
  std::size_t pruned_count() const {
    return candidate_count() - selected_count();
  }
  bool respects_width() const {
    for (const auto &[j, selected] : selected_per_column_) {
      (void)j;
      if (selected > width_) return false;
    }
    return true;
  }
  std::size_t record_probe_count() const { return record_probe_count_; }

 private:
  std::size_t width_ = 0;
  std::unordered_set<VpCase3TailKey, VpCase3TailKeyHash> selected_;
  std::unordered_map<cand_pos_t, std::size_t> candidate_per_column_;
  std::unordered_map<cand_pos_t, std::size_t> selected_per_column_;
  std::size_t record_probe_count_ = 0;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_SELECTION_HH

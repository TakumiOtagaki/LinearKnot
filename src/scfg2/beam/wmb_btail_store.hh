#ifndef CPARTY_SCFG2_BEAM_WMB_BTAIL_STORE_HH
#define CPARTY_SCFG2_BEAM_WMB_BTAIL_STORE_HH

#include "scfg2/deductions/item_key.hh"
#include "scfg2/beam/memory_observation.hh"
#include "scfg2/beam/wmb_btail_seeds.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

template <class Value>
struct WmbBtailCell {
  cand_pos_t i = -1;
  Value total{};
  cand_pos_t best_k = -1;
};

struct WmbBtailStoreStats {
  std::size_t join_probes = 0;
  std::size_t direct_probes = 0;
  std::size_t selected = 0;
};

template <class Semiring>
class WmbBtailStore {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context, class WmbpStartsFn>
  void build(cand_pos_t t,
             const Chart &chart,
             const Context &ctx,
             const std::vector<cand_pos_t> &wi_starts,
             WmbpStartsFn &&wmbp_starts,
             std::size_t width) {
    if (t < 2 || ctx.pair_partner(t) < 1 || ctx.pair_partner(t) >= t) return;
    const cand_pos_t right = t - 1;
    // Only the aggregate and representative survive this build. Keep the
    // comparison value alive during selection, then release its expression.
    struct Candidate {
      WmbBtailCell<Value> cell;
      Value best_term;
    };
    std::unordered_map<cand_pos_t, Candidate> candidates;
    const auto add_seed = [&](cand_pos_t i, cand_pos_t k, Value wi) {
      const Value wmbp = chart.get(deductions::ItemKey{
          grammar::NonTerminal::WMBP, i, k, -1, -1});
      const Value term = Semiring::times(wmbp, wi);
      auto found = candidates.try_emplace(
          i,
          Candidate{{i, Semiring::zero(), -1}, Semiring::zero()}).first;
      auto &candidate = found->second;
      auto &cell = candidate.cell;
      cell.total = Semiring::plus(cell.total, term);
      constexpr double tie_tolerance = 1e-12;
      const double term_priority = Semiring::traceback_priority(term);
      const double best_priority =
          Semiring::traceback_priority(candidate.best_term);
      const bool better = term_priority > best_priority + tie_tolerance;
      const bool tied = cell.best_k >= 0 &&
          std::fabs(term_priority - best_priority) <= tie_tolerance;
      if (cell.best_k < 0 || better || (tied && k < cell.best_k)) {
        candidate.best_term = term;
        cell.best_k = k;
      }
    };
    const auto counts = for_each_wmb_btail_seed(t, ctx, wi_starts, wmbp_starts,
        [&](cand_pos_t i, cand_pos_t k, bool has_wi) {
          add_seed(i, k, has_wi ? chart.get(deductions::ItemKey{
              grammar::NonTerminal::WI, k + 1, right, -1, -1}) : Semiring::one());
        });
    stats_.direct_probes += counts.direct;
    stats_.join_probes += counts.joined;
    std::vector<WmbBtailCell<Value>> ranked;
    for (const auto &[i, candidate] : candidates)
      if (!Semiring::is_zero(candidate.cell.total))
        ranked.push_back(candidate.cell);
    std::sort(ranked.begin(), ranked.end(), [](const auto &lhs,
                                                const auto &rhs) {
      const double lhs_priority = Semiring::priority(lhs.total);
      const double rhs_priority = Semiring::priority(rhs.total);
      if (lhs_priority != rhs_priority) return lhs_priority > rhs_priority;
      return lhs.i < rhs.i;
    });
    if (width != std::numeric_limits<std::size_t>::max() &&
        ranked.size() > width)
      ranked.resize(width);
    stats_.selected += ranked.size();
    for (auto &cell : ranked) {
      cell.total = semirings::checkpoint_value(cell.total);
    }
    // resize(width) retains the O(b^2) candidate vector's allocation. Persist
    // only selected cells: otherwise all closing columns retain O(nb^2)
    // capacity even though their logical cell count is O(nb).
    if (width != std::numeric_limits<std::size_t>::max() &&
        ranked.capacity() > ranked.size()) {
      std::vector<WmbBtailCell<Value>>(ranked.begin(), ranked.end()).swap(ranked);
    }
    by_close_[t] = std::move(ranked);
  }

  const std::vector<WmbBtailCell<Value>> &retained(cand_pos_t t) const {
    const auto found = by_close_.find(t);
    if (found == by_close_.end()) {
      static const std::vector<WmbBtailCell<Value>> empty;
      return empty;
    }
    return found->second;
  }
  const WmbBtailStoreStats &stats() const { return stats_; }

  ContainerMemory memory_snapshot(const char *name) const {
    return hash_vector_memory(name, by_close_);
  }

 private:
  WmbBtailStoreStats stats_{};
  std::unordered_map<cand_pos_t, std::vector<WmbBtailCell<Value>>> by_close_;
};

}  // namespace cparty::scfg2::beam

#endif

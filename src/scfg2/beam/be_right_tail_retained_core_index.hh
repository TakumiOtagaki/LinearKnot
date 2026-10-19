#ifndef CPARTY_SCFG2_BEAM_BE_RIGHT_TAIL_RETAINED_CORE_INDEX_HH
#define CPARTY_SCFG2_BEAM_BE_RIGHT_TAIL_RETAINED_CORE_INDEX_HH

#include "scfg2/beam/be_right_tail_selection.hh"

#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

struct BeRightTailRetainedCoreStats {
  std::size_t retained_core_events = 0;
  std::size_t immediate_parent_probes = 0;
  std::size_t registered_tails = 0;
  std::size_t completed_wip_events = 0;
  std::size_t materialized_tails = 0;
  std::size_t committed_tails = 0;
  std::size_t duplicate_registrations = 0;
};

// Dynamic inverse of the right-WIP dependency. A retained anchored BE core
// registers at most one future BE_RTAIL, using its immediate scaffold parent.
// Completion of that tail's WIP child materializes the factor for the target
// outer-pair column.
template <class Value>
class BasicBeRightTailRetainedCoreIndex {
 public:
  template <class Context>
  static BasicBeRightTailRetainedCoreIndex build(cand_pos_t n,
                                                  const Context &ctx) {
    BasicBeRightTailRetainedCoreIndex out;
    std::vector<cand_pos_t> open;
    for (cand_pos_t pos = 1; pos <= n; ++pos) {
      const cand_pos_t partner = ctx.pair_partner(pos);
      if (partner > pos) {
        if (!open.empty()) out.parent_by_open_[pos] = open.back();
        open.push_back(pos);
      } else if (partner >= 1 && partner < pos && !open.empty() &&
                 open.back() == partner) {
        open.pop_back();
      }
    }
    return out;
  }

  template <class Context>
  void observe_retained_be_column(
      const std::vector<deductions::ItemKey> &retained_be,
      const Context &ctx) {
    for (const auto &core : retained_be) {
      if (core.nonterminal != grammar::NonTerminal::BE) continue;
      ++stats_.retained_core_events;
      const auto parent = parent_by_open_.find(core.i);
      if (parent == parent_by_open_.end()) continue;
      ++stats_.immediate_parent_probes;
      const cand_pos_t i = parent->second;
      const cand_pos_t j = ctx.pair_partner(i);
      if (!(j > core.j && ctx.weakly_closed(core.j + 1, j - 1) &&
            (ctx.is_empty_region(i, core.i) ||
             ctx.weakly_closed(i + 1, core.i - 1)))) {
        continue;
      }
      const BeRightTailKey key{
          j, BeCoreKey{core.i, core.j, core.ip, core.jp}};
      if (!keys_.insert(key).second) {
        ++stats_.duplicate_registrations;
        continue;
      }
      const deductions::ItemKey right_child{
          grammar::NonTerminal::WIP,
          core.j + 1,
          j - 1,
          -1,
          -1};
      by_right_child_[right_child].push_back(key);
      ++stats_.registered_tails;
    }
  }

  template <class Semiring, class Chart>
  void observe_retained_child(const deductions::ItemKey &item,
                              const Chart &chart) {
    if (item.nonterminal != grammar::NonTerminal::WIP) return;
    const auto found = by_right_child_.find(item);
    if (found == by_right_child_.end()) return;
    ++stats_.completed_wip_events;
    for (const auto &key : found->second) {
      const Value factor = Semiring::times(
          chart.get(deductions::ItemKey{grammar::NonTerminal::BE,
                                       key.core.k, key.core.l,
                                       key.core.ip, key.core.jp}),
          chart.get(item));
      if (Semiring::is_zero(factor)) continue;
      pending_[key.j].push_back({key, factor});
      ++stats_.materialized_tails;
    }
    // No future column can complete this same right child again.
    by_right_child_.erase(found);
  }

  template <class Semiring, class Chart>
  void observe_completed_child(const deductions::ItemKey &item,
                               const Chart &chart,
                               BasicExactBeFrameCache<Value> &cache) {
    if (item.nonterminal != grammar::NonTerminal::WIP) return;
    const auto found = by_right_child_.find(item);
    if (found == by_right_child_.end()) return;
    ++stats_.completed_wip_events;
    for (const auto &key : found->second) {
      const Value factor =
          cache.template materialize_right_tail<Semiring>(key, chart);
      if (Semiring::is_zero(factor)) continue;
      pending_[key.j].push_back({key, factor});
      ++stats_.materialized_tails;
    }
  }

  template <class Chart>
  void observe_completed_child(const deductions::ItemKey &item,
                               const Chart &chart,
                               ExactBeFrameCache &cache) {
    observe_completed_child<semirings::SumProductSemiring>(
        item, chart, cache);
  }

  std::vector<BasicScoredBeRightTail<Value>> take_column(cand_pos_t j) {
    const auto found = pending_.find(j);
    if (found == pending_.end()) return {};
    auto out = std::move(found->second);
    pending_.erase(found);
    stats_.committed_tails += out.size();
    return out;
  }

  const BeRightTailRetainedCoreStats &stats() const { return stats_; }

  // Observation of registration geometry, before completion removes a child.
  // A reverse consumer can retain only these IDs and recompute its products.
  template <class Fn>
  void for_each_pending_child_tail(Fn &&fn) const {
    for (const auto &[child, tails] : by_right_child_)
      for (const auto &tail : tails) fn(child, tail);
  }
  std::size_t scaffold_parent_link_count() const {
    return parent_by_open_.size();
  }
  bool pending_empty() const { return pending_.empty(); }

 private:
  std::unordered_map<cand_pos_t, cand_pos_t> parent_by_open_;
  std::unordered_map<deductions::ItemKey,
                     std::vector<BeRightTailKey>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      by_right_child_;
  std::unordered_set<BeRightTailKey, BeRightTailKeyHash> keys_;
  std::unordered_map<cand_pos_t,
                     std::vector<BasicScoredBeRightTail<Value>>>
      pending_;
  BeRightTailRetainedCoreStats stats_{};
};

using BeRightTailRetainedCoreIndex =
    BasicBeRightTailRetainedCoreIndex<double>;

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_RIGHT_TAIL_RETAINED_CORE_INDEX_HH

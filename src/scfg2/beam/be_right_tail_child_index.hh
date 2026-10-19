#ifndef CPARTY_SCFG2_BEAM_BE_RIGHT_TAIL_CHILD_INDEX_HH
#define CPARTY_SCFG2_BEAM_BE_RIGHT_TAIL_CHILD_INDEX_HH

#include "scfg2/beam/be_right_tail_selection.hh"

#include <cstddef>
#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

// Derived discovery index. The right WIP is the last dependency to complete
// for BE_RTAIL^R; its completion materializes every future tail using the
// already-complete anchored BE core. The exact column plan remains the oracle
// for geometry only and is not consulted by the execution callbacks.
class ExactBeRightTailChildIndex {
 public:
  explicit ExactBeRightTailChildIndex(
      const ExactBeRightTailColumnPlan &plan) {
    plan.for_each_record([&](const BeRightTailPlanRecord &record) {
      const auto &deduction = record.representative;
      const auto &right_child =
          deduction.children[deduction.child_count - 1].item;
      insert(right_child, record.key);
    });
  }

  template <class Context>
  static ExactBeRightTailChildIndex build_from_scaffold_coordinate_oracle(
      cand_pos_t n,
      const Context &ctx) {
    ExactBeRightTailChildIndex out;
    for (cand_pos_t j = 1; j <= n; ++j) {
      const cand_pos_t i = ctx.pair_partner(j);
      if (!(i >= 1 && i < j)) continue;
      for (cand_pos_t k = i + 1; k < j; ++k) {
        ++out.coordinate_core_visit_count_;
        const cand_pos_t l = ctx.pair_partner(k);
        if (!(l > k && l < j) ||
            !ctx.weakly_closed(l + 1, j - 1)) {
          continue;
        }
        if (!(ctx.is_empty_region(i, k) ||
              ctx.weakly_closed(i + 1, k - 1))) {
          continue;
        }
        for (cand_pos_t ip = k; ip < l; ++ip) {
          ++out.coordinate_anchor_visit_count_;
          const cand_pos_t jp = ctx.pair_partner(ip);
          if (!(jp > ip && jp <= l)) continue;
          const BeRightTailKey key{
              j, BeCoreKey{k, l, ip, jp}};
          out.insert(deductions::ItemKey{grammar::NonTerminal::WIP,
                                         l + 1,
                                         j - 1,
                                         -1,
                                         -1},
                     key);
        }
      }
    }
    return out;
  }

  template <class Context>
  static ExactBeRightTailChildIndex build_from_scaffold(cand_pos_t n,
                                                        const Context &ctx) {
    return build_from_scaffold_coordinate_oracle(n, ctx);
  }

  template <class Context>
  static ExactBeRightTailChildIndex build_from_scaffold_pair_tree(
      cand_pos_t n,
      const Context &ctx) {
    struct PairNode {
      cand_pos_t i = -1;
      cand_pos_t j = -1;
      std::size_t parent = no_node();
      std::vector<std::size_t> children;

      static constexpr std::size_t no_node() {
        return static_cast<std::size_t>(-1);
      }
    };

    ExactBeRightTailChildIndex out;
    std::vector<PairNode> pairs;
    std::vector<std::size_t> open;
    for (cand_pos_t pos = 1; pos <= n; ++pos) {
      const cand_pos_t partner = ctx.pair_partner(pos);
      if (partner > pos) {
        const std::size_t parent =
            open.empty() ? PairNode::no_node() : open.back();
        const std::size_t node = pairs.size();
        pairs.push_back(PairNode{pos, partner, parent, {}});
        if (parent != PairNode::no_node()) {
          pairs[parent].children.push_back(node);
        }
        open.push_back(node);
      } else if (partner >= 1 && partner < pos) {
        if (!open.empty() && pairs[open.back()].j == pos) open.pop_back();
      }
    }
    out.scaffold_pair_count_ = pairs.size();

    for (std::size_t core_idx = 0; core_idx < pairs.size(); ++core_idx) {
      const auto &core = pairs[core_idx];
      for (std::size_t outer_idx = core.parent;
           outer_idx != PairNode::no_node();
           outer_idx = pairs[outer_idx].parent) {
        ++out.core_outer_visit_count_;
        const auto &outer = pairs[outer_idx];
        if (!ctx.weakly_closed(core.j + 1, outer.j - 1) ||
            !(ctx.is_empty_region(outer.i, core.i) ||
              ctx.weakly_closed(outer.i + 1, core.i - 1))) {
          continue;
        }
        std::vector<std::size_t> pending{core_idx};
        while (!pending.empty()) {
          const std::size_t anchor_idx = pending.back();
          pending.pop_back();
          const auto &anchor = pairs[anchor_idx];
          const BeRightTailKey key{
              outer.j,
              BeCoreKey{core.i, core.j, anchor.i, anchor.j}};
          out.insert(deductions::ItemKey{grammar::NonTerminal::WIP,
                                         core.j + 1,
                                         outer.j - 1,
                                         -1,
                                         -1},
                     key);
          ++out.anchor_visit_count_;
          pending.insert(
              pending.end(), anchor.children.begin(), anchor.children.end());
        }
      }
    }
    return out;
  }

  template <class Semiring, class Chart>
  void observe_completed_child(const deductions::ItemKey &item,
                               const Chart &chart,
                               ExactBeFrameCache &cache) {
    if (item.nonterminal != grammar::NonTerminal::WIP) return;
    const auto found = by_right_child_.find(item);
    if (found == by_right_child_.end()) return;
    for (const auto &key : found->second) {
      const double factor =
          cache.template materialize_right_tail<Semiring>(key, chart);
      pending_[key.j].push_back({key, factor});
      ++materialized_count_;
    }
  }

  template <class Chart>
  void observe_completed_child(const deductions::ItemKey &item,
                               const Chart &chart,
                               ExactBeFrameCache &cache) {
    observe_completed_child<semirings::SumProductSemiring>(
        item, chart, cache);
  }

  std::vector<ScoredBeRightTail> take_column(cand_pos_t j) {
    const auto found = pending_.find(j);
    if (found == pending_.end()) return {};
    auto out = std::move(found->second);
    pending_.erase(found);
    committed_count_ += out.size();
    return out;
  }

  std::size_t record_count() const { return record_count_; }
  std::size_t materialized_count() const { return materialized_count_; }
  std::size_t committed_count() const { return committed_count_; }
  bool pending_empty() const { return pending_.empty(); }
  std::size_t scaffold_pair_count() const { return scaffold_pair_count_; }
  std::size_t core_outer_visit_count() const {
    return core_outer_visit_count_;
  }
  std::size_t anchor_visit_count() const { return anchor_visit_count_; }
  std::size_t coordinate_core_visit_count() const {
    return coordinate_core_visit_count_;
  }
  std::size_t coordinate_anchor_visit_count() const {
    return coordinate_anchor_visit_count_;
  }
  bool has_same_keys(const ExactBeRightTailChildIndex &other) const {
    return keys_ == other.keys_;
  }

 private:
  ExactBeRightTailChildIndex() = default;

  void insert(const deductions::ItemKey &right_child,
              const BeRightTailKey &key) {
    if (!keys_.insert(key).second) return;
    by_right_child_[right_child].push_back(key);
    ++record_count_;
  }

  std::unordered_map<deductions::ItemKey,
                     std::vector<BeRightTailKey>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      by_right_child_;
  std::unordered_set<BeRightTailKey, BeRightTailKeyHash> keys_;
  std::unordered_map<cand_pos_t, std::vector<ScoredBeRightTail>> pending_;
  std::size_t record_count_ = 0;
  std::size_t materialized_count_ = 0;
  std::size_t committed_count_ = 0;
  std::size_t scaffold_pair_count_ = 0;
  std::size_t core_outer_visit_count_ = 0;
  std::size_t anchor_visit_count_ = 0;
  std::size_t coordinate_core_visit_count_ = 0;
  std::size_t coordinate_anchor_visit_count_ = 0;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_RIGHT_TAIL_CHILD_INDEX_HH

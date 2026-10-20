#ifndef CPARTY_SCFG2_BEAM_VPR_RETAINED_AGENDA_HH
#define CPARTY_SCFG2_BEAM_VPR_RETAINED_AGENDA_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/run_term_observer.hh"
#include "scfg2/beam/vpr_base_prefix.hh"
#include "scfg2/deductions/expand_basic.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct VprRetainedAgendaStats {
  std::size_t binary_join_probes = 0;
  std::size_t base_state_probes = 0;
  std::size_t base_state_selected = 0;
};

template <class Semiring>
class VprRetainedAgenda {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context, class VpStartsFn>
  void build_column(cand_pos_t j,
                    const Chart &chart,
                    const Context &ctx,
                    const std::vector<cand_pos_t> &wip_starts,
                    VpStartsFn &&vp_starts,
                    std::size_t width) {
    const bool follows_previous_column = column_ == j - 1;
    column_ = j;
    deductions_by_parent_.clear();
    base_contribution_by_parent_.clear();
    parents_.clear();

    std::unordered_map<cand_pos_t, VprBasePrefixCell<Value>> candidates;
    const bool extends = j > 1 && ctx.can_pair_right_span(j - 1, j);
    std::unordered_map<cand_pos_t, bool> carried_starts;
    if (extends && follows_previous_column) {
      for (const auto &[i, cell] : base_states_) carried_starts[i] = true;
    }
    if (extends) {
      // Every term VP(i,k) * cp^(j-k) gains one unpaired base per column.
      const Value one_unpaired =
          Semiring::from_weight(weights::local_weight(ctx, ctx.vpr_split_vp_basepair_weight(1)));
      for (const auto &[i, cell] : base_states_) {
        ++stats_.base_state_probes;
        candidates.emplace(
            i, VprBasePrefixCell<Value>{
                   Semiring::times(cell.total, one_unpaired),
                   Semiring::times(cell.best_term, one_unpaired),
                   cell.best_k});
      }
      const cand_pos_t k = j - 1;
      for (const cand_pos_t i : vp_starts(k)) {
        ++stats_.base_state_probes;
        if (i >= k) continue;
        auto [candidate, inserted] = candidates.try_emplace(
            i,
            VprBasePrefixCell<Value>{
                Semiring::zero(), Semiring::zero(), -1});
        (void)inserted;
        auto &cell = candidate->second;
        const Value term = Semiring::times(
            chart.get(deductions::ItemKey{
                grammar::NonTerminal::VP, i, k, -1, -1}),
            one_unpaired);
        cell.total = Semiring::plus(cell.total, term);
        if (Semiring::better(term, cell.best_term)) {
          cell.best_term = term;
          cell.best_k = k;
        }
      }
    }
    std::vector<std::pair<cand_pos_t, VprBasePrefixCell<Value>>> ranked(
        candidates.begin(), candidates.end());
    ranked.erase(std::remove_if(ranked.begin(), ranked.end(), [&](const auto &x) {
                   return Semiring::is_zero(x.second.total);
                 }), ranked.end());
    std::sort(ranked.begin(), ranked.end(), [](const auto &lhs,
                                                const auto &rhs) {
      const double lhs_priority = Semiring::priority(lhs.second.total);
      const double rhs_priority = Semiring::priority(rhs.second.total);
      if (lhs_priority != rhs_priority)
        return lhs_priority > rhs_priority;
      return lhs.first < rhs.first;
    });
    if (width != std::numeric_limits<std::size_t>::max() &&
        ranked.size() > width) {
      ranked.resize(width);
    }
    base_states_.clear();
    for (const auto &[i, cell] : ranked) {
      auto sealed = cell;
      sealed.total = semirings::checkpoint_value(cell.total);
      base_states_[i] = sealed;
      ++stats_.base_state_selected;
      if (auto *observer = run_term_observer<Value>()) {
        observer->on_vpr_run_cell(i, j, cell.total, carried_starts.count(i) > 0);
      }
    }

    const auto append = [&](const deductions::ItemKey &parent,
                            grammar::RuleId rule,
                            cand_pos_t k) {
      grammar::RuleSplit split{};
      split.k = k;
      const auto deduction = deductions::expand_basic_rule(parent, rule, split);
      if (!deductions::is_valid(deduction)) return;
      auto [found, inserted] = deductions_by_parent_.try_emplace(parent);
      if (inserted) parents_.push_back(parent);
      found->second.push_back(deduction);
    };

    for (const auto &[i, cell] : base_states_) {
      const deductions::ItemKey parent{
          grammar::NonTerminal::VPR, i, j, -1, -1};
      if (!deductions::exact_basic_item_parent_active(parent, ctx)) continue;
      const cand_pos_t max_bp =
          std::max(ctx.border_B(i, j), ctx.border_bp(i, j));
      if (cell.best_k <= max_bp) continue;
      append(parent, grammar::RuleId::VPR_SPLIT_VP_BASEPAIR, cell.best_k);
      base_contribution_by_parent_[parent] = cell.total;
    }

    for (const cand_pos_t wip_i : wip_starts) {
      const cand_pos_t k = wip_i - 1;
      if (k < 1) continue;
      for (const cand_pos_t i : vp_starts(k)) {
        ++stats_.binary_join_probes;
        const deductions::ItemKey parent{
            grammar::NonTerminal::VPR, i, j, -1, -1};
        if (!deductions::exact_basic_item_parent_active(parent, ctx)) continue;
        const cand_pos_t max_bp =
            std::max(ctx.border_B(i, j), ctx.border_bp(i, j));
        if (k <= i || k <= max_bp) continue;
        append(parent, grammar::RuleId::VPR_SPLIT_VP_WIP, k);
      }
    }

    std::sort(parents_.begin(), parents_.end(), [](const auto &lhs,
                                                   const auto &rhs) {
      return lhs.i > rhs.i;
    });
  }

  std::vector<cand_pos_t> parent_starts() const {
    std::vector<cand_pos_t> starts;
    for (const auto &parent : parents_) starts.push_back(parent.i);
    return starts;
  }

  const std::vector<deductions::Deduction> &deductions_for(
      const deductions::ItemKey &parent) const {
    const auto found = deductions_by_parent_.find(parent);
    if (found == deductions_by_parent_.end()) {
      static const std::vector<deductions::Deduction> empty;
      return empty;
    }
    return found->second;
  }

  // Consume this completed parent's source buffer once. Aggregate values and
  // parent lists stay available until the next column; source views expire.
  std::vector<deductions::Deduction> take_for(
      const deductions::ItemKey &parent) {
    const auto found = deductions_by_parent_.find(parent);
    if (found == deductions_by_parent_.end()) {
      return {};
    }
    return std::move(found->second);
  }

  // Observe a private aggregate with one map probe. The source remains
  // owned by this column; the pointer expires at the next build_column().
  const Value *find_base_contribution(const deductions::Deduction &deduction) const {
    if (deduction.rule != grammar::RuleId::VPR_SPLIT_VP_BASEPAIR) return nullptr;
    const auto found = base_contribution_by_parent_.find(deduction.parent);
    return found == base_contribution_by_parent_.end() ? nullptr : &found->second;
  }

  bool reconstructs_base(const deductions::Deduction &deduction) const {
    return deduction.rule == grammar::RuleId::VPR_SPLIT_VP_BASEPAIR &&
           base_contribution_by_parent_.find(deduction.parent) !=
               base_contribution_by_parent_.end();
  }
  Value base_contribution(const deductions::ItemKey &parent) const {
    const auto found = base_contribution_by_parent_.find(parent);
    return found == base_contribution_by_parent_.end() ? Semiring::zero()
                                                       : found->second;
  }
  cand_pos_t column() const { return column_; }
  const VprRetainedAgendaStats &stats() const { return stats_; }

 private:
  cand_pos_t column_ = -1;
  VprRetainedAgendaStats stats_{};
  std::unordered_map<cand_pos_t, VprBasePrefixCell<Value>> base_states_;
  std::unordered_map<deductions::ItemKey, std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash, deductions::ItemKeyEq>
      deductions_by_parent_;
  std::unordered_map<deductions::ItemKey, Value,
                     deductions::ItemKeyHash, deductions::ItemKeyEq>
      base_contribution_by_parent_;
  std::vector<deductions::ItemKey> parents_;
};

}  // namespace cparty::scfg2::beam

#endif

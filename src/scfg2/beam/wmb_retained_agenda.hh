#ifndef CPARTY_SCFG2_BEAM_WMB_RETAINED_AGENDA_HH
#define CPARTY_SCFG2_BEAM_WMB_RETAINED_AGENDA_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/wmb_btail_store.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/weights/wmb_transition.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct WmbRetainedAgendaStats {
  std::size_t be_btail_join_probes = 0;
  std::size_t direct_probes = 0;
  std::size_t split_selected = 0;
};

// Shared expansion of a selected BTAIL factor against a retained BE. The
// representative split owns exactly the same guards and local weight as in
// forward; reverse consumes this view rather than inventing a second gate.
template <class Context, class Semiring, class EmitFn>
std::size_t for_each_wmb_be_btail_term(
    cand_pos_t j, const Context &ctx,
    const std::vector<deductions::ItemKey> &be_items,
    const WmbBtailStore<Semiring> &btails, EmitFn &&emit) {
  std::size_t probes = 0;
  for (const auto &be : be_items) {
    if (be.nonterminal != grammar::NonTerminal::BE || be.j != j) continue;
    for (const auto &tail : btails.retained(be.jp)) {
      ++probes;
      if (tail.i >= be.i) continue;
      const deductions::ItemKey parent{grammar::NonTerminal::WMB, tail.i, j, -1, -1};
      std::vector<deductions::Deduction> generated;
      deductions::append_exact_wmb_split_for_wmbp_end(generated, parent, ctx, tail.best_k);
      if (generated.empty()) continue;
      const auto &source = generated.front();
      if (source.split.p == be.i && source.split.q == be.ip && source.split.l == be.jp)
        emit(source, be, tail, weights::local_weight(ctx, weights::transition_weight_wmb(source, ctx)));
    }
  }
  return probes;
}

template <class Semiring>
class WmbRetainedAgenda {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context>
  void build_column(cand_pos_t j,
                    const Chart &chart,
                    const Context &ctx,
                    const std::vector<cand_pos_t> &direct_wmbp_starts,
                    const std::vector<deductions::ItemKey> &be_items,
                    const WmbBtailStore<Semiring> &btails) {
    column_ = j;
    deductions_by_parent_.clear();
    split_contribution_by_parent_.clear();
    parents_.clear();
    const auto append = [&](const deductions::Deduction &deduction) {
      if (!deductions::is_valid(deduction)) return;
      auto [found, inserted] =
          deductions_by_parent_.try_emplace(deduction.parent);
      if (inserted) parents_.push_back(deduction.parent);
      found->second.push_back(deduction);
    };
    for (const cand_pos_t i : direct_wmbp_starts) {
      ++stats_.direct_probes;
      const deductions::ItemKey parent{
          grammar::NonTerminal::WMB, i, j, -1, -1};
      if (!deductions::exact_basic_item_parent_active(parent, ctx)) continue;
      append(deductions::expand_basic_rule(
          parent, grammar::RuleId::WMB_DIRECT_WMBP, grammar::RuleSplit{}));
    }

    struct Aggregate {
      Value total = Semiring::zero();
      Value best = Semiring::zero();
      deductions::Deduction representative{};
      bool has_representative = false;
    };
    std::unordered_map<deductions::ItemKey, Aggregate,
                       deductions::ItemKeyHash, deductions::ItemKeyEq>
        aggregates;
    stats_.be_btail_join_probes += for_each_wmb_be_btail_term(
        j, ctx, be_items, btails,
        [&](const auto &deduction, const auto &be_item, const auto &tail, auto local) {
          const auto &parent = deduction.parent;
          const Value be = chart.get(be_item);
          const Value term = Semiring::times(
              Semiring::times(be, tail.total),
              Semiring::from_weight(local));
          auto &aggregate = aggregates[parent];
          aggregate.total = Semiring::plus(aggregate.total, term);
          constexpr double tie_tolerance = 1e-12;
          const double term_priority = Semiring::traceback_priority(term);
          const double best_priority =
              Semiring::traceback_priority(aggregate.best);
          const bool better = term_priority > best_priority + tie_tolerance;
          const bool tied = aggregate.has_representative &&
              std::fabs(term_priority - best_priority) <= tie_tolerance;
          const auto tie_key = [](const deductions::Deduction &source) {
            return std::make_tuple(static_cast<unsigned short>(source.rule),
                                   source.split.k,
                                   source.split.l,
                                   source.split.p,
                                   source.split.q);
          };
          if (!aggregate.has_representative || better ||
              (tied && tie_key(deduction) <
                           tie_key(aggregate.representative))) {
            aggregate.best = term;
            aggregate.representative = deduction;
            aggregate.has_representative = true;
          }
        });
    for (const auto &[parent, aggregate] : aggregates) {
      if (!aggregate.has_representative || Semiring::is_zero(aggregate.total))
        continue;
      append(aggregate.representative);
      split_contribution_by_parent_[parent] = aggregate.total;
      ++stats_.split_selected;
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
  bool reconstructs_split(const deductions::Deduction &deduction) const {
    return deduction.rule == grammar::RuleId::WMB_SPLIT_BE_WMBP_WI &&
           split_contribution_by_parent_.find(deduction.parent) !=
               split_contribution_by_parent_.end();
  }
  Value split_contribution(const deductions::ItemKey &parent) const {
    const auto found = split_contribution_by_parent_.find(parent);
    return found == split_contribution_by_parent_.end() ? Semiring::zero()
                                                        : found->second;
  }
  cand_pos_t column() const { return column_; }
  const WmbRetainedAgendaStats &stats() const { return stats_; }

 private:
  cand_pos_t column_ = -1;
  WmbRetainedAgendaStats stats_{};
  std::unordered_map<deductions::ItemKey, std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash, deductions::ItemKeyEq>
      deductions_by_parent_;
  std::unordered_map<deductions::ItemKey, Value,
                     deductions::ItemKeyHash, deductions::ItemKeyEq>
      split_contribution_by_parent_;
  std::vector<deductions::ItemKey> parents_;
};

}  // namespace cparty::scfg2::beam

#endif

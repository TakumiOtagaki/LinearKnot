#ifndef CPARTY_SCFG2_BEAM_VPL_RETAINED_AGENDA_HH
#define CPARTY_SCFG2_BEAM_VPL_RETAINED_AGENDA_HH

#include "scfg2/beam/column_store.hh"

#include "scfg2/beam/vpl_start_prefix.hh"
#include "scfg2/deductions/expand_basic.hh"
#include "scfg2/deductions/generate_exact_basic.hh"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct VplRetainedAgendaStats {
  std::size_t prefix_cell_visits = 0;
  std::size_t prefix_selected = 0;
};

template <class Semiring>
class VplRetainedAgenda {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context>
  void build_column(cand_pos_t j,
                    const Chart &chart,
                    const Context &ctx,
                    const std::vector<cand_pos_t> &vp_starts,
                    std::size_t width) {
    column_ = j;
    deductions_by_parent_.reset_column(j);
    contribution_by_parent_.reset_column(j);
    parents_.clear();
    const auto append = [&](cand_pos_t i, Value aggregate, cand_pos_t best_k) {
      const deductions::ItemKey parent{
          grammar::NonTerminal::VPL, i, j, -1, -1};
      if (!deductions::exact_basic_item_parent_active(parent, ctx)) return;
      grammar::RuleSplit split{};
      split.k = best_k;
      const auto deduction = deductions::expand_basic_rule(
          parent, grammar::RuleId::VPL_SPLIT_VP, split);
      if (!deductions::is_valid(deduction)) return;
      deductions_by_parent_[parent].push_back(deduction);
      contribution_by_parent_[parent] = aggregate;
      parents_.push_back(parent);
      ++stats_.prefix_selected;
    };
    if (width == std::numeric_limits<std::size_t>::max()) {
      prefix_.build(j, chart, ctx);
      for (cand_pos_t i = j - 1; i >= 1; --i) {
        ++stats_.prefix_cell_visits;
        const Value aggregate = prefix_.contribution(i);
        if (Semiring::is_zero(aggregate)) continue;
        append(i, aggregate, prefix_.cell(i).best_k);
      }
    } else {
      const auto selected = select_vpl_start_prefix_top_width<Semiring>(
          j, chart, ctx, vp_starts, width);
      for (const auto &cell : selected) {
        append(cell.i, cell.total, cell.best_k);
      }
    }
    std::sort(parents_.begin(), parents_.end(), [](const auto &lhs,
                                                   const auto &rhs) {
      return lhs.i > rhs.i;
    });
  }

  std::vector<cand_pos_t> parent_starts() const {
    std::vector<cand_pos_t> starts;
    starts.reserve(parents_.size());
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
  const Value *find_contribution(const deductions::Deduction &deduction) const {
    if (deduction.rule != grammar::RuleId::VPL_SPLIT_VP) return nullptr;
    const auto found = contribution_by_parent_.find(deduction.parent);
    return found == contribution_by_parent_.end() ? nullptr : &found->second;
  }

  bool reconstructs(const deductions::Deduction &deduction) const {
    return deduction.rule == grammar::RuleId::VPL_SPLIT_VP &&
           contribution_by_parent_.find(deduction.parent) !=
               contribution_by_parent_.end();
  }

  Value contribution(const deductions::ItemKey &parent) const {
    const auto found = contribution_by_parent_.find(parent);
    return found == contribution_by_parent_.end() ? Semiring::zero()
                                                  : found->second;
  }

  cand_pos_t column() const { return column_; }
  const VplRetainedAgendaStats &stats() const { return stats_; }

 private:
  cand_pos_t column_ = -1;
  VplRetainedAgendaStats stats_{};
  ColumnStore<std::vector<deductions::Deduction>, grammar::NonTerminal::VPL, 8> deductions_by_parent_;
  ColumnStore<Value, grammar::NonTerminal::VPL> contribution_by_parent_;
  std::vector<deductions::ItemKey> parents_;
  VplStartPrefixColumn<Semiring> prefix_;
};

}  // namespace cparty::scfg2::beam

#endif

#ifndef CPARTY_SCFG2_BEAM_VM_RETAINED_JOIN_AGENDA_HH
#define CPARTY_SCFG2_BEAM_VM_RETAINED_JOIN_AGENDA_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/expand_basic.hh"
#include "scfg2/beam/vm_wmp_base_prefix.hh"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct VmRetainedJoinAgendaStats {
  std::size_t binary_join_probes = 0;
  std::size_t base_prefix_cell_visits = 0;
  std::size_t base_prefix_selected = 0;
};

template <class Semiring>
class VmRetainedJoinAgenda {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context, class RetainedStartsFn>
  void build_column(cand_pos_t j,
                    const Chart &chart,
                    const Context &ctx,
                    RetainedStartsFn &&retained_starts,
                    std::size_t base_prefix_width) {
    column_ = j;
    deductions_by_parent_.clear();
    base_contribution_by_parent_.clear();
    parents_.clear();

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

    const auto append_binary_family = [&](grammar::NonTerminal right_family,
                                          grammar::RuleId rule) {
      for (const cand_pos_t k : retained_starts(right_family, j - 1)) {
        for (const cand_pos_t left_start :
             retained_starts(grammar::NonTerminal::WM, k - 1)) {
          ++stats_.binary_join_probes;
          const cand_pos_t parent_i = left_start - 1;
          if (parent_i < 1 || k < parent_i + 2 ||
              k > j - ctx.turn() - 1) {
            continue;
          }
          append(deductions::ItemKey{
                     grammar::NonTerminal::VM, parent_i, j, -1, -1},
                 rule,
                 k);
        }
      }
    };

    append_binary_family(grammar::NonTerminal::WMv,
                         grammar::RuleId::VM_SPLIT_WM_WMv);
    append_binary_family(grammar::NonTerminal::WMp,
                         grammar::RuleId::VM_SPLIT_WM_WMp);

    const auto append_base = [&](cand_pos_t parent_i,
                                 Value tail,
                                 cand_pos_t best_k) {
      const deductions::ItemKey parent{
          grammar::NonTerminal::VM, parent_i, j, -1, -1};
      append(parent, grammar::RuleId::VM_SPLIT_WMp_BASE, best_k);
      base_contribution_by_parent_[parent] = Semiring::times(
          Semiring::from_weight(weights::local_weight(ctx, ctx.split_wm_wmp_weight(parent_i, j))),
          tail);
      ++stats_.base_prefix_selected;
    };
    if (base_prefix_width == std::numeric_limits<std::size_t>::max()) {
      base_prefix_.build(j, chart, ctx);
      const cand_pos_t max_k = j - ctx.turn() - 1;
      for (cand_pos_t parent_i = max_k - 1; parent_i >= 1; --parent_i) {
        ++stats_.base_prefix_cell_visits;
        const auto &cell = base_prefix_.cell(parent_i);
        if (Semiring::is_zero(cell.tail)) continue;
        append_base(parent_i, cell.tail, cell.best_k);
      }
    } else {
      const auto selected = select_vm_wmp_base_prefix_top_width<Semiring>(
          j,
          chart,
          ctx,
          retained_starts(grammar::NonTerminal::WMp, j - 1),
          base_prefix_width);
      for (const auto &cell : selected) {
        append_base(cell.i, cell.tail, cell.best_k);
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
  const Value *find_base_contribution(const deductions::Deduction &deduction) const {
    if (deduction.rule != grammar::RuleId::VM_SPLIT_WMp_BASE) return nullptr;
    const auto found = base_contribution_by_parent_.find(deduction.parent);
    return found == base_contribution_by_parent_.end() ? nullptr : &found->second;
  }

  bool reconstructs_base(const deductions::Deduction &deduction) const {
    return deduction.rule == grammar::RuleId::VM_SPLIT_WMp_BASE &&
           base_contribution_by_parent_.find(deduction.parent) !=
               base_contribution_by_parent_.end();
  }

  Value base_contribution(const deductions::ItemKey &parent) const {
    const auto found = base_contribution_by_parent_.find(parent);
    return found == base_contribution_by_parent_.end() ? Semiring::zero()
                                                       : found->second;
  }

  cand_pos_t column() const { return column_; }
  const VmRetainedJoinAgendaStats &stats() const { return stats_; }

 private:
  cand_pos_t column_ = -1;
  VmRetainedJoinAgendaStats stats_{};
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      deductions_by_parent_;
  std::unordered_map<deductions::ItemKey,
                     Value,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      base_contribution_by_parent_;
  std::vector<deductions::ItemKey> parents_;
  VmWmpBasePrefixColumn<Semiring> base_prefix_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VM_RETAINED_JOIN_AGENDA_HH

#ifndef CPARTY_SCFG2_BEAM_WM_RETAINED_JOIN_AGENDA_HH
#define CPARTY_SCFG2_BEAM_WM_RETAINED_JOIN_AGENDA_HH

#include "scfg2/beam/wm_start_prefix.hh"
#include "scfg2/deductions/expand_basic.hh"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct WmRetainedJoinAgendaStats {
  std::size_t binary_join_probes = 0;
  std::size_t extension_probes = 0;
  std::size_t start_prefix_cell_visits = 0;
  std::size_t start_prefix_selected = 0;
};

template <class Semiring>
class WmRetainedJoinAgenda {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context, class RetainedStartsFn>
  void build_column(cand_pos_t j,
                    const Chart &chart,
                    const Context &ctx,
                    RetainedStartsFn &&retained_starts,
                    std::size_t start_prefix_width,
                    const std::vector<cand_pos_t> *forced_start_parents = nullptr) {
    column_ = j;
    deductions_by_parent_.clear();
    start_contribution_by_parent_.clear();
    parents_.clear();

    const auto append = [&](const deductions::ItemKey &parent,
                            grammar::RuleId rule,
                            cand_pos_t k = -1) {
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
      for (const cand_pos_t k : retained_starts(right_family, j)) {
        for (const cand_pos_t parent_i :
             retained_starts(grammar::NonTerminal::WM, k - 1)) {
          ++stats_.binary_join_probes;
          if (parent_i < 1 || k <= parent_i || j - parent_i + 1 < 4) {
            continue;
          }
          append(deductions::ItemKey{
                     grammar::NonTerminal::WM, parent_i, j, -1, -1},
                 rule,
                 k);
        }
      }
    };
    append_binary_family(grammar::NonTerminal::V,
                         grammar::RuleId::WM_SPLIT_V);
    append_binary_family(grammar::NonTerminal::WMB,
                         grammar::RuleId::WM_SPLIT_WMB);

    for (const cand_pos_t parent_i :
         retained_starts(grammar::NonTerminal::WM, j - 1)) {
      ++stats_.extension_probes;
      if (!ctx.is_unpaired(j) || j - parent_i + 1 < 4) continue;
      append(deductions::ItemKey{
                 grammar::NonTerminal::WM, parent_i, j, -1, -1},
             grammar::RuleId::WM_EXTEND_UNPAIRED);
    }

    const auto append_start = [&](cand_pos_t parent_i,
                                  Value aggregate,
                                  grammar::RuleId best_rule,
                                  cand_pos_t best_k) {
      const deductions::ItemKey parent{
          grammar::NonTerminal::WM, parent_i, j, -1, -1};
      append(parent, best_rule, best_k);
      start_contribution_by_parent_[parent] = aggregate;
      ++stats_.start_prefix_selected;
    };
    if (start_prefix_width == std::numeric_limits<std::size_t>::max()) {
      std::vector<bool> seed_positions(static_cast<std::size_t>(j) + 1, false);
      for (const auto k : retained_starts(grammar::NonTerminal::V, j))
        seed_positions[k] = true;
      for (const auto k : retained_starts(grammar::NonTerminal::WMB, j))
        seed_positions[k] = true;
      start_prefix_.build(j, chart, ctx,
          [&](cand_pos_t i) { return seed_positions[i]; });
      for (cand_pos_t parent_i = j - 3; parent_i >= 1; --parent_i) {
        ++stats_.start_prefix_cell_visits;
        const Value aggregate = start_prefix_.contribution(parent_i);
        if (Semiring::is_zero(aggregate)) continue;
        const auto &cell = start_prefix_.cell(parent_i);
        append_start(parent_i, aggregate, cell.best_rule, cell.best_k);
      }
    } else {
      const auto selected = select_wm_start_prefix_top_width<Semiring>(
          j,
          chart,
          ctx,
          retained_starts(grammar::NonTerminal::V, j),
          retained_starts(grammar::NonTerminal::WMB, j),
          start_prefix_width,
          forced_start_parents);
      for (const auto &cell : selected) {
        append_start(
            cell.i, cell.total, cell.best_rule, cell.best_k);
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
  const Value *find_start_contribution(const deductions::Deduction &deduction) const {
    if (deduction.rule != grammar::RuleId::WM_START_V &&
        deduction.rule != grammar::RuleId::WM_START_WMB) return nullptr;
    const auto found = start_contribution_by_parent_.find(deduction.parent);
    return found == start_contribution_by_parent_.end() ? nullptr : &found->second;
  }

  bool reconstructs_start(const deductions::Deduction &deduction) const {
    return (deduction.rule == grammar::RuleId::WM_START_V ||
            deduction.rule == grammar::RuleId::WM_START_WMB) &&
           start_contribution_by_parent_.find(deduction.parent) !=
               start_contribution_by_parent_.end();
  }

  Value start_contribution(const deductions::ItemKey &parent) const {
    const auto found = start_contribution_by_parent_.find(parent);
    return found == start_contribution_by_parent_.end() ? Semiring::zero()
                                                        : found->second;
  }

  cand_pos_t column() const { return column_; }
  const WmRetainedJoinAgendaStats &stats() const { return stats_; }

 private:
  cand_pos_t column_ = -1;
  WmRetainedJoinAgendaStats stats_{};
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      deductions_by_parent_;
  std::unordered_map<deductions::ItemKey,
                     Value,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      start_contribution_by_parent_;
  std::vector<deductions::ItemKey> parents_;
  WmStartPrefixColumn<Semiring> start_prefix_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WM_RETAINED_JOIN_AGENDA_HH

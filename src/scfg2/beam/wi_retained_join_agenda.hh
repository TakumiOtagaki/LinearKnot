#ifndef CPARTY_SCFG2_BEAM_WI_RETAINED_JOIN_AGENDA_HH
#define CPARTY_SCFG2_BEAM_WI_RETAINED_JOIN_AGENDA_HH

#include "scfg2/deductions/expand_basic.hh"
#include "scfg2/beam/column_store.hh"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct WiRetainedJoinAgendaStats {
  std::size_t binary_join_probes = 0;
  std::size_t direct_probes = 0;
  std::size_t extension_probes = 0;
};

class WiRetainedJoinAgenda {
 public:
  template <class Context, class RetainedStartsFn>
  void build_column(cand_pos_t j,
                    const Context &ctx,
                    RetainedStartsFn &&retained_starts) {
    column_ = j;
    deductions_by_parent_.reset_column(j);
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

    if (ctx.is_unpaired(j)) {
      append(deductions::ItemKey{
                 grammar::NonTerminal::WI, j, j, -1, -1},
             grammar::RuleId::WI_BASE_SINGLE);
    }

    const auto append_branch_family = [&](grammar::NonTerminal right_family,
                                          grammar::RuleId rule) {
      for (const cand_pos_t k : retained_starts(right_family, j)) {
        ++stats_.direct_probes;
        append(deductions::ItemKey{
                   grammar::NonTerminal::WI, k, j, -1, -1},
               rule,
               k);
        for (const cand_pos_t parent_i :
             retained_starts(grammar::NonTerminal::WI, k - 1)) {
          ++stats_.binary_join_probes;
          if (parent_i >= k) continue;
          append(deductions::ItemKey{
                     grammar::NonTerminal::WI, parent_i, j, -1, -1},
                 rule,
                 k);
        }
      }
    };
    append_branch_family(grammar::NonTerminal::V,
                         grammar::RuleId::WI_SPLIT_V);
    append_branch_family(grammar::NonTerminal::WMB,
                         grammar::RuleId::WI_SPLIT_WMB);

    for (const cand_pos_t parent_i :
         retained_starts(grammar::NonTerminal::WI, j - 1)) {
      ++stats_.extension_probes;
      if (!ctx.is_unpaired(j)) continue;
      append(deductions::ItemKey{
                 grammar::NonTerminal::WI, parent_i, j, -1, -1},
             grammar::RuleId::WI_EXTEND_UNPAIRED);
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

  cand_pos_t column() const { return column_; }
  const WiRetainedJoinAgendaStats &stats() const { return stats_; }

 private:
  cand_pos_t column_ = -1;
  WiRetainedJoinAgendaStats stats_{};
  // Only canonical primary parents of this one column are authoritative.
  // Source arrival order remains within each buffer; parents_ owns traversal.
  ColumnStore<std::vector<deductions::Deduction>, grammar::NonTerminal::WI, 8>
      deductions_by_parent_;
  std::vector<deductions::ItemKey> parents_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WI_RETAINED_JOIN_AGENDA_HH

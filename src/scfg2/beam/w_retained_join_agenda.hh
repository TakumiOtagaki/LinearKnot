#ifndef CPARTY_SCFG2_BEAM_W_RETAINED_JOIN_AGENDA_HH
#define CPARTY_SCFG2_BEAM_W_RETAINED_JOIN_AGENDA_HH

#include "scfg2/beam/column_store.hh"

#include "scfg2/deductions/generate_exact_basic.hh"

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct WRetainedJoinAgendaStats {
  std::size_t binary_join_probes = 0;
  std::size_t direct_probes = 0;
  std::size_t extension_probes = 0;
  std::size_t root_decomposition_probes = 0;
};

class WRetainedJoinAgenda {
 public:
  template <class Context, class RetainedStartsFn>
  void build_column(cand_pos_t j,
                    const Context &ctx,
                    RetainedStartsFn &&retained_starts) {
    build_column(j, ctx, retained_starts, [&](cand_pos_t i, cand_pos_t k) {
      const auto &starts = retained_starts(grammar::NonTerminal::WMB, k);
      return std::find(starts.begin(), starts.end(), i) != starts.end();
    });
  }

  template <class Context, class RetainedStartsFn, class WmbContainsFn>
  void build_column(cand_pos_t j, const Context &ctx,
                    RetainedStartsFn &&retained_starts, WmbContainsFn &&wmb_contains) {
    column_ = j;
    deductions_by_parent_.reset_column(j);
    parents_.clear();
    const auto append_deduction = [&](const deductions::Deduction &deduction) {
      if (!deductions::is_valid(deduction)) return;
      auto [found, inserted] =
          deductions_by_parent_.try_emplace(deduction.parent);
      if (inserted) parents_.push_back(deduction.parent);
      found->second.push_back(deduction);
    };
    const auto append = [&](const deductions::ItemKey &parent,
                            auto rule,
                            cand_pos_t k = -1,
                            cand_pos_t p = -1) {
      grammar::RuleSplit split{};
      split.k = k;
      split.p = p;
      // Fixed binary tags select the existing shared schema. Extension and
      // root decomposition keep the original runtime expansion route.
      const auto deduction = [&] {
        using RuleSelector = decltype(rule);
        if constexpr (std::is_same_v<RuleSelector, grammar::RuleId>)
          return deductions::expand_basic_rule(parent, rule, split);
        else
          return deductions::expand_basic_rule<RuleSelector::value>(parent, split);
      }();
      append_deduction(deduction);
    };

    append_deduction(deductions::make_zero_child_deduction(
        deductions::ItemKey{grammar::NonTerminal::W, j, j, -1, -1},
        grammar::RuleId::W_EXTEND_UNPAIRED,
        grammar::RuleSplit{}));

    for (const cand_pos_t i :
         retained_starts(grammar::NonTerminal::W, j - 1)) {
      ++stats_.extension_probes;
      if (!ctx.is_unpaired(j) || i >= j) continue;
      append(deductions::ItemKey{
                 grammar::NonTerminal::W, i, j, -1, -1},
             grammar::RuleId::W_EXTEND_UNPAIRED);
    }

    const auto append_branch = [&](grammar::NonTerminal family,
                                   auto rule) {
      for (const cand_pos_t k : retained_starts(family, j)) {
        ++stats_.direct_probes;
        const deductions::ItemKey direct{
            grammar::NonTerminal::W, k, j, -1, -1};
        const bool common = ctx.weakly_closed(1, j) &&
                            ctx.weakly_closed(1, k - 1);
        const bool direct_allowed =
            common && (rule != grammar::RuleId::W_SPLIT_V ||
                       ctx.can_pair(k, j));
        if (direct_allowed) append(direct, rule, k);
        for (const cand_pos_t i :
             retained_starts(grammar::NonTerminal::W, k - 1)) {
          ++stats_.binary_join_probes;
          if (i >= k || !common) continue;
          if (rule == grammar::RuleId::W_SPLIT_V) {
            if (!ctx.can_pair(k, j)) continue;
          } else if (!ctx.weakly_closed(k, j)) {
            continue;
          }
          append(deductions::ItemKey{
                     grammar::NonTerminal::W, i, j, -1, -1},
                 rule,
                 k);
        }
      }
    };
    append_branch(grammar::NonTerminal::V,
                  std::integral_constant<grammar::RuleId,
                                         grammar::RuleId::W_SPLIT_V>{});
    append_branch(grammar::NonTerminal::WMB,
                  std::integral_constant<grammar::RuleId,
                                         grammar::RuleId::W_SPLIT_WMB>{});

    const auto ordinary_parents = parents_;
    for (const auto &right_parent : ordinary_parents) {
      const cand_pos_t k = right_parent.i - 1;
      if (k < 1 || k > j) continue;
      ++stats_.root_decomposition_probes;
      const cand_pos_t p = deductions::root_wmb_open(ctx, 1, k);
      if (!(p > 0 && p < k) || !ctx.weakly_closed(1, p - 1) ||
          !(k == j || ctx.weakly_closed(k + 1, j))) {
        continue;
      }
      if (!wmb_contains(p, k)) {
        continue;
      }
      append(deductions::ItemKey{
                 grammar::NonTerminal::W, 1, j, -1, -1},
             grammar::RuleId::W_ROOT_WMB_DECOMPOSITION,
             k,
             p);
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
  cand_pos_t column() const { return column_; }
  const WRetainedJoinAgendaStats &stats() const { return stats_; }

 private:
  cand_pos_t column_ = -1;
  WRetainedJoinAgendaStats stats_{};
  ColumnStore<std::vector<deductions::Deduction>, grammar::NonTerminal::W, 8>
      deductions_by_parent_;
  std::vector<deductions::ItemKey> parents_;
};

}  // namespace cparty::scfg2::beam

#endif

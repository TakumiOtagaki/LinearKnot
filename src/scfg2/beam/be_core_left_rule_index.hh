#ifndef CPARTY_SCFG2_BEAM_BE_CORE_LEFT_RULE_INDEX_HH
#define CPARTY_SCFG2_BEAM_BE_CORE_LEFT_RULE_INDEX_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

// Stores retained child identities and one event per left WIP, rather than
// the Cartesian product of WIP prefixes and anchored BE cores. Both forward
// evaluation and reverse outside can enumerate any column from this index.
// With <= b retained BE/WIP states per column, the index occupies O(nb), and
// visiting every column costs O(nb^2). No numerical snapshots are stored.
class BeCoreLeftRuleIndex {
 public:
  // Outside can rebuild this view once from the final chart, without saving
  // a forward agenda or replaying its history. Each committed child identity
  // must be registered once (the online executor already guarantees this).
  template <class Semiring, class Chart, class Context>
  static BeCoreLeftRuleIndex from_retained_chart(
      cand_pos_t n, const Chart &chart, const Context &ctx) {
    auto out = build(n, ctx);
    std::vector<std::vector<cand_pos_t>> wip_starts(
        static_cast<std::size_t>(n) + 1);
    for (const auto &[item, value] : chart.values()) {
      if (Semiring::is_zero(value)) continue;
      if (item.nonterminal == grammar::NonTerminal::BE)
        out.observe_committed_be_column({item});
      else if (item.nonterminal == grammar::NonTerminal::WIP)
        wip_starts[item.j].push_back(item.i);
    }
    for (cand_pos_t j = 1; j <= n; ++j)
      out.observe_retained_wip_column(wip_starts[j], j, ctx);
    return out;
  }

  template <class Context>
  static BeCoreLeftRuleIndex build(cand_pos_t n, const Context &ctx) {
    BeCoreLeftRuleIndex out;
    out.events_.resize(static_cast<std::size_t>(n) + 1);
    out.cores_by_open_.resize(static_cast<std::size_t>(n) + 1);
    out.left_counts_by_core_.resize(static_cast<std::size_t>(n) + 1);
    std::vector<cand_pos_t> open;
    for (cand_pos_t pos = 1; pos <= n; ++pos) {
      const cand_pos_t partner = ctx.pair_partner(pos);
      if (partner > pos) {
        if (!open.empty() && ctx.is_empty_region(open.back(), pos)) {
          const auto close = ctx.pair_partner(open.back());
          out.events_[close].push_back(
              Event{open.back(), pos, grammar::RuleId::BE_INTERNAL_LOOP});
          ++out.event_count_;
          out.last_event_column_ = std::max(out.last_event_column_, close);
        }
        open.push_back(pos);
      } else if (partner >= 1 && partner < pos && !open.empty() &&
                 open.back() == partner) {
        open.pop_back();
      }
    }
    return out;
  }

  template <class Context>
  void observe_retained_wip_column(
      const std::vector<cand_pos_t> &starts, cand_pos_t right,
      const Context &ctx) {
    const cand_pos_t k = right + 1;
    const cand_pos_t l = ctx.pair_partner(k);
    if (l <= k) return;
    for (const cand_pos_t a : starts) {
      const cand_pos_t i = a - 1;
      if (!(i >= 1 && i < k && ctx.pair_partner(i) > l &&
            ctx.weakly_closed(i + 1, k - 1))) continue;
      const auto close = ctx.pair_partner(i);
      events_[close].push_back(
          Event{i, k, grammar::RuleId::BE_WIP_BASEPAIR});
      last_event_column_ = std::max(last_event_column_, close);
      ++event_count_;
      const auto count = ++left_counts_by_core_[k];
      max_left_prefixes_ = std::max(max_left_prefixes_, count);
    }
  }

  void observe_committed_be_column(
      const std::vector<deductions::ItemKey> &retained) {
    for (const auto &core : retained) {
      if (core.nonterminal != grammar::NonTerminal::BE) continue;
      cores_by_open_[core.i].push_back(core);
      ++core_count_;
    }
  }

  template <class Context, class EmitFn, class ProbeFn>
  void for_each_column_source(cand_pos_t j, Context &ctx,
                              EmitFn &&emit, ProbeFn &&probe) const {
    if (j < 1 || static_cast<std::size_t>(j) >= events_.size()) return;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    for (const auto &event : events_[j]) {
      for (const auto &core : cores_by_open_[event.k]) {
        probe(event.rule);
        const deductions::ItemKey parent{
            grammar::NonTerminal::BE, event.i, j, core.ip, core.jp};
        if (!deductions::exact_be_split_allows(
                parent, event.rule, core.i, core.j, ctx)) continue;
        grammar::RuleSplit split{};
        split.k = core.i;
        split.l = core.j;
        const auto source =
            deductions::expand_basic_rule(parent, event.rule, split);
        if (!deductions::is_valid(source)) continue;
        engine::require_legacy_span_schedule_support(source);
        if (weights::is_zero(weight.local_weight_for(source))) continue;
        const std::size_t core_slot =
            event.rule == grammar::RuleId::BE_WIP_BASEPAIR ? 1 : 0;
        if (!deductions::ItemKeyEq{}(source.children[core_slot].item, core))
          continue;
        emit(source);
      }
    }
  }

  template <class Context, class EmitFn>
  void for_each_column_source(cand_pos_t j, Context &ctx, EmitFn &&emit) const {
    for_each_column_source(j, ctx, emit, [](grammar::RuleId) {});
  }

  std::size_t core_count() const { return core_count_; }
  std::size_t event_count() const { return event_count_; }
  std::size_t max_left_prefixes_per_core() const { return max_left_prefixes_; }
  bool has_events_after(cand_pos_t j) const {
    return last_event_column_ > j;
  }

 private:
  struct Event {
    cand_pos_t i;
    cand_pos_t k;
    grammar::RuleId rule;
  };
  std::vector<std::vector<Event>> events_;
  std::vector<std::vector<deductions::ItemKey>> cores_by_open_;
  std::vector<std::size_t> left_counts_by_core_;
  cand_pos_t last_event_column_ = -1;
  std::size_t core_count_ = 0;
  std::size_t event_count_ = 0;
  std::size_t max_left_prefixes_ = 0;
};

}  // namespace cparty::scfg2::beam

#endif

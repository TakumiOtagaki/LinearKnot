#ifndef CPARTY_SCFG2_BEAM_WIP_START_PREFIX_HH
#define CPARTY_SCFG2_BEAM_WIP_START_PREFIX_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/run_term_observer.hh"
#include "scfg2/deductions/item_key.hh"
#include "scfg2/grammar/rules.hh"

#include <algorithm>
#include <cstddef>
#include <queue>
#include <vector>

namespace cparty::scfg2::beam {

template <class Value>
struct WipStartPrefixCell {
  Value total{};
  Value best_term{};
  cand_pos_t best_k = -1;
  grammar::RuleId best_rule = grammar::RuleId::WIP_BASE_V;
};

template <class Semiring>
class WipStartPrefixColumn {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context>
  void build(cand_pos_t j, const Chart &chart, const Context &ctx) {
    build(j, chart, ctx, [](cand_pos_t) { return true; });
  }

  // The agenda supplies the same merged V/internal-WMB seed positions as
  // its finite selector. A root-only WMB may still enter at a V seed.
  template <class Chart, class Context, class SeedAllowsFn>
  void build(cand_pos_t j, const Chart &chart, const Context &ctx,
             SeedAllowsFn &&seed_allows) {
    cells_.assign(static_cast<std::size_t>(j + 2),
                  WipStartPrefixCell<Value>{Semiring::zero(),
                                            Semiring::zero(),
                                            -1,
                                            grammar::RuleId::WIP_BASE_V});
    Value tail = Semiring::zero();
    Value tail_best = Semiring::zero();
    cand_pos_t tail_best_k = -1;
    grammar::RuleId tail_best_base_rule = grammar::RuleId::WIP_BASE_V;
    cand_pos_t last_unit_base = -1;
    for (cand_pos_t i = j; i >= 1; --i) {
      auto &cell = cells_[static_cast<std::size_t>(i)];
      const bool seed = seed_allows(i);
      const Value v_seed = seed ? Semiring::times(
          chart.get(deductions::ItemKey{
              grammar::NonTerminal::V, i, j, -1, -1}),
          Semiring::from_weight(weights::local_weight(ctx, ctx.wip_base_v_weight()))) : Semiring::zero();
      const Value wmb_seed = seed ? Semiring::times(
          chart.get(deductions::ItemKey{
              grammar::NonTerminal::WMB, i, j, -1, -1}),
          Semiring::from_weight(weights::local_weight(ctx, ctx.wip_base_wmb_weight()))) : Semiring::zero();
      cell.total = Semiring::plus(v_seed, wmb_seed);
      if (Semiring::better(wmb_seed, v_seed)) {
        cell.best_term = wmb_seed;
        cell.best_rule = grammar::RuleId::WIP_BASE_WMB;
      } else {
        cell.best_term = v_seed;
        cell.best_rule = grammar::RuleId::WIP_BASE_V;
      }
      if (!Semiring::is_zero(cell.best_term)) cell.best_k = i;

      Value extended_unit = Semiring::zero();
      if (ctx.can_pair_left_span(i, i + 1)) {
        const Value extended =
            Semiring::times(
                Semiring::from_weight(weights::local_weight(ctx, ctx.wip_extend_unpaired_weight())),
                tail);
        cell.total = Semiring::plus(cell.total, extended);
        extended_unit = extended;
        const Value extended_best = Semiring::times(
            Semiring::from_weight(weights::local_weight(ctx, ctx.wip_extend_unpaired_weight())),
            tail_best);
        if (Semiring::better(extended_best, cell.best_term)) {
          cell.best_term = extended_best;
          cell.best_k = tail_best_k;
          cell.best_rule =
              tail_best_base_rule == grammar::RuleId::WIP_BASE_V
                  ? grammar::RuleId::WIP_BASEPAIR_V
                  : grammar::RuleId::WIP_BASEPAIR_WMB;
        }
      }

      const bool extendable_seed = i <= j - ctx.turn() - 1;
      Value next_tail = Semiring::zero();
      Value next_tail_best = Semiring::zero();
      cand_pos_t next_tail_best_k = -1;
      grammar::RuleId next_tail_rule = grammar::RuleId::WIP_BASE_V;
      if (extendable_seed) {
        next_tail = Semiring::plus(v_seed, wmb_seed);
        if (Semiring::better(wmb_seed, v_seed)) {
          next_tail_best = wmb_seed;
          next_tail_rule = grammar::RuleId::WIP_BASE_WMB;
        } else {
          next_tail_best = v_seed;
          next_tail_rule = grammar::RuleId::WIP_BASE_V;
        }
        if (!Semiring::is_zero(next_tail_best)) next_tail_best_k = i;
      }
      if (ctx.can_pair_left_span(i, i + 1)) {
        const Value carried =
            Semiring::times(
                Semiring::from_weight(weights::local_weight(ctx, ctx.wip_extend_unpaired_weight())),
                tail);
        next_tail = Semiring::plus(next_tail, carried);
        if (!Semiring::is_zero(carried)) {
          notify_left_run_unit(grammar::RuleId::WIP_BASEPAIR_V, j, i,
                               extended_unit, &carried,
                               last_unit_base == i + 1);
          last_unit_base = i;
        }
        const Value carried_best = Semiring::times(
            Semiring::from_weight(weights::local_weight(ctx, ctx.wip_extend_unpaired_weight())),
            tail_best);
        if (Semiring::better(carried_best, next_tail_best)) {
          next_tail_best = carried_best;
          next_tail_best_k = tail_best_k;
          next_tail_rule = tail_best_base_rule;
        }
      }
      cell.total = semirings::checkpoint_value(cell.total);
      tail = semirings::checkpoint_value(next_tail);
      tail_best = next_tail_best;
      tail_best_k = next_tail_best_k;
      tail_best_base_rule = next_tail_rule;
    }
  }

  Value contribution(cand_pos_t i) const {
    return i < 1 || static_cast<std::size_t>(i) >= cells_.size()
               ? Semiring::zero()
               : cells_[static_cast<std::size_t>(i)].total;
  }

  const WipStartPrefixCell<Value> &cell(cand_pos_t i) const {
    return cells_[static_cast<std::size_t>(i)];
  }

 private:
  std::vector<WipStartPrefixCell<Value>> cells_;
};

template <class Value>
struct WipStartPrefixSelectedCell : WipStartPrefixCell<Value> {
  cand_pos_t i = -1;
};

template <class Semiring, class Chart, class Context>
std::vector<WipStartPrefixSelectedCell<typename Semiring::value_type>>
select_wip_start_prefix_top_width(
    cand_pos_t j,
    const Chart &chart,
    const Context &ctx,
    const std::vector<cand_pos_t> &v_starts,
    const std::vector<cand_pos_t> &wmb_starts,
    std::size_t width) {
  using Value = typename Semiring::value_type;
  struct Seed {
    cand_pos_t k = -1;
    cand_pos_t left = -1;
    Value v;
    Value wmb;
  };
  std::vector<cand_pos_t> starts = v_starts;
  starts.insert(starts.end(), wmb_starts.begin(), wmb_starts.end());
  std::sort(starts.begin(), starts.end());
  starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
  std::vector<Seed> seeds;
  std::vector<cand_pos_t> cuts;
  for (const cand_pos_t k : starts) {
    const Value v = chart.get(deductions::ItemKey{
        grammar::NonTerminal::V, k, j, -1, -1});
    const Value wmb = chart.get(deductions::ItemKey{
        grammar::NonTerminal::WMB, k, j, -1, -1});
    if (Semiring::is_zero(v) && Semiring::is_zero(wmb)) continue;
    const cand_pos_t left = k <= j - ctx.turn() - 1
                                ? std::max<cand_pos_t>(
                                      1, k - ctx.unpaired_prefix(k - 1))
                                : k;
    seeds.push_back(Seed{k, left, v, wmb});
    cuts.push_back(left);
    cuts.push_back(k + 1);
  }
  if (seeds.empty() || width == 0) return {};
  // WIP exists only on the round weakly-closed family (plus the diagonal).
  // Split the heap intervals at that family boundary so an inactive endpoint
  // cannot hide a later active cell in the same monotone interval.
  const cand_pos_t unpaired_suffix = ctx.unpaired_prefix(j);
  if (unpaired_suffix > 0) {
    cuts.push_back(std::max<cand_pos_t>(1, j - unpaired_suffix + 1));
  }
  cuts.push_back(j);
  cuts.push_back(j + 1);
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

  const auto evaluate = [&](cand_pos_t i) {
    WipStartPrefixSelectedCell<Value> out;
    out.i = i;
    out.total = Semiring::zero();
    out.best_term = Semiring::zero();
    if (i != j && !ctx.weakly_closed(i, j)) return out;
    for (const auto &seed : seeds) {
      if (i < seed.left || i > seed.k ||
          !ctx.can_pair_left_span(i, seed.k)) {
        continue;
      }
      const bool direct = i == seed.k;
      const Value v_term = Semiring::times(
          seed.v,
          Semiring::from_weight(weights::local_weight(ctx,
              direct ? ctx.wip_base_v_weight()
                     : ctx.wip_basepair_v_weight(seed.k - i))));
      const Value wmb_term = Semiring::times(
          seed.wmb,
          Semiring::from_weight(weights::local_weight(ctx,
              direct ? ctx.wip_base_wmb_weight()
                     : ctx.wip_basepair_wmb_weight(seed.k - i))));
      out.total = Semiring::plus(out.total, v_term);
      out.total = Semiring::plus(out.total, wmb_term);
      if (!direct) {
        notify_run_term(grammar::RuleId::WIP_BASEPAIR_V, i, seed.k - 1, v_term);
        notify_run_term(grammar::RuleId::WIP_BASEPAIR_WMB, i, seed.k - 1, wmb_term);
      }
      if (Semiring::better(v_term, out.best_term)) {
        out.best_term = v_term;
        out.best_k = seed.k;
        out.best_rule = direct ? grammar::RuleId::WIP_BASE_V
                               : grammar::RuleId::WIP_BASEPAIR_V;
      }
      if (Semiring::better(wmb_term, out.best_term)) {
        out.best_term = wmb_term;
        out.best_k = seed.k;
        out.best_rule = direct ? grammar::RuleId::WIP_BASE_WMB
                               : grammar::RuleId::WIP_BASEPAIR_WMB;
      }
    }
    return out;
  };

  struct Cursor {
    cand_pos_t i = -1;
    cand_pos_t end = -1;
    cand_pos_t step = 0;
    WipStartPrefixSelectedCell<Value> cell{};
  };
  struct Worse {
    bool operator()(const Cursor &lhs, const Cursor &rhs) const {
      const double lhs_priority = Semiring::priority(lhs.cell.total);
      const double rhs_priority = Semiring::priority(rhs.cell.total);
      if (lhs_priority != rhs_priority) return lhs_priority < rhs_priority;
      return lhs.i > rhs.i;
    }
  };
  std::priority_queue<Cursor, std::vector<Cursor>, Worse> heap;
  const bool toward_left =
      Semiring::priority(
          Semiring::from_weight(weights::local_weight(ctx, ctx.wip_extend_unpaired_weight()))) <
      Semiring::priority(Semiring::one());
  for (std::size_t cut = 0; cut + 1 < cuts.size(); ++cut) {
    const cand_pos_t lo = cuts[cut];
    const cand_pos_t hi = cuts[cut + 1] - 1;
    if (lo > hi) continue;
    const cand_pos_t first = toward_left ? hi : lo;
    const cand_pos_t end = toward_left ? lo : hi;
    const cand_pos_t step = toward_left ? -1 : 1;
    auto cell = evaluate(first);
    if (!Semiring::is_zero(cell.total))
      heap.push(Cursor{first, end, step, cell});
  }
  std::vector<WipStartPrefixSelectedCell<Value>> selected;
  // Width is a retention limit, not a number of feasible parent starts.
  // In particular a mixed KeepAll/finite forcing pass can use SIZE_MAX/2.
  selected.reserve(std::min(width, static_cast<std::size_t>(
      std::max<cand_pos_t>(0, j))));
  while (!heap.empty() && selected.size() < width) {
    auto cursor = heap.top();
    heap.pop();
    selected.push_back(cursor.cell);
    if (cursor.i == cursor.end) continue;
    cursor.i += cursor.step;
    cursor.cell = evaluate(cursor.i);
    if (!Semiring::is_zero(cursor.cell.total)) heap.push(cursor);
  }
  return selected;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WIP_START_PREFIX_HH

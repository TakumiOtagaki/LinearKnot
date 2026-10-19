#ifndef CPARTY_SCFG2_BEAM_VPL_START_PREFIX_HH
#define CPARTY_SCFG2_BEAM_VPL_START_PREFIX_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/run_term_observer.hh"
#include "scfg2/deductions/item_key.hh"

#include <algorithm>
#include <cstddef>
#include <queue>
#include <vector>

namespace cparty::scfg2::beam {

template <class Value>
struct VplStartPrefixCell {
  Value total{};
  Value best_term{};
  cand_pos_t best_k = -1;
};

template <class Context>
bool vpl_seed_border_allows(cand_pos_t k, cand_pos_t j, const Context &ctx) {
  if (k <= 1) return false;
  const cand_pos_t parent_i = k - 1;
  const cand_pos_tu cap = std::min(
      static_cast<cand_pos_tu>(ctx.border_b(parent_i, j)),
      static_cast<cand_pos_tu>(ctx.border_Bp(parent_i, j)));
  return static_cast<cand_pos_tu>(k) < cap;
}

template <class Semiring>
class VplStartPrefixColumn {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context>
  void build(cand_pos_t j, const Chart &chart, const Context &ctx) {
    cells_.assign(static_cast<std::size_t>(j + 1),
                  VplStartPrefixCell<Value>{Semiring::zero(),
                                            Semiring::zero(), -1});
    Value carried = Semiring::zero();
    Value carried_best = Semiring::zero();
    cand_pos_t carried_best_k = -1;
    cand_pos_t last_unit_base = -1;
    for (cand_pos_t i = j - 1; i >= 1; --i) {
      const cand_pos_t k = i + 1;
      Value seed = Semiring::zero();
      if (k < j && vpl_seed_border_allows(k, j, ctx) &&
          ctx.can_pair_left_span(i, k)) {
        seed = Semiring::times(
            chart.get(deductions::ItemKey{
                grammar::NonTerminal::VP, k, j, -1, -1}),
            Semiring::from_weight(weights::local_weight(ctx, ctx.vpl_split_vp_weight(1))));
      }
      auto &cell = cells_[static_cast<std::size_t>(i)];
      cell.total = seed;
      cell.best_term = seed;
      if (!Semiring::is_zero(seed)) cell.best_k = k;
      const bool extended_from_next = ctx.can_pair_left_span(i, i + 1);
      if (extended_from_next) {
        const Value extended = Semiring::times(
            Semiring::from_weight(weights::local_weight(ctx, ctx.vpl_split_vp_weight(1))), carried);
        cell.total = Semiring::plus(cell.total, extended);
        const Value extended_best = Semiring::times(
            Semiring::from_weight(weights::local_weight(ctx, ctx.vpl_split_vp_weight(1))), carried_best);
        if (Semiring::better(extended_best, cell.best_term)) {
          cell.best_term = extended_best;
          cell.best_k = carried_best_k;
        }
      }
      if (!Semiring::is_zero(cell.total)) {
        notify_left_run_unit(grammar::RuleId::VPL_SPLIT_VP, j, i, cell.total,
                             static_cast<const Value *>(nullptr),
                             extended_from_next && last_unit_base == i + 1);
        last_unit_base = i;
      }
      cell.total = semirings::checkpoint_value(cell.total);
      carried = cell.total;
      carried_best = cell.best_term;
      carried_best_k = cell.best_k;
    }
  }

  Value contribution(cand_pos_t i) const {
    return i < 1 || static_cast<std::size_t>(i) >= cells_.size()
               ? Semiring::zero()
               : cells_[static_cast<std::size_t>(i)].total;
  }

  const VplStartPrefixCell<Value> &cell(cand_pos_t i) const {
    return cells_[static_cast<std::size_t>(i)];
  }

 private:
  std::vector<VplStartPrefixCell<Value>> cells_;
};

template <class Value>
struct VplStartPrefixSelectedCell : VplStartPrefixCell<Value> {
  cand_pos_t i = -1;
};

template <class Semiring, class Chart, class Context>
std::vector<VplStartPrefixSelectedCell<typename Semiring::value_type>>
select_vpl_start_prefix_top_width(cand_pos_t j,
                                  const Chart &chart,
                                  const Context &ctx,
                                  const std::vector<cand_pos_t> &vp_starts,
                                  std::size_t width) {
  using Value = typename Semiring::value_type;
  struct Seed { cand_pos_t k; cand_pos_t left; };
  std::vector<Seed> seeds;
  std::vector<cand_pos_t> cuts;
  for (const cand_pos_t k : vp_starts) {
    if (k >= j || !vpl_seed_border_allows(k, j, ctx)) continue;
    if (Semiring::is_zero(chart.get(deductions::ItemKey{
            grammar::NonTerminal::VP, k, j, -1, -1}))) {
      continue;
    }
    const cand_pos_t left =
        std::max<cand_pos_t>(1, k - ctx.unpaired_prefix(k - 1));
    if (left >= k) continue;
    seeds.push_back(Seed{k, left});
    cuts.push_back(left);
    cuts.push_back(k);
  }
  if (seeds.empty() || width == 0) return {};
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

  const auto evaluate = [&](cand_pos_t i) {
    VplStartPrefixSelectedCell<Value> out;
    out.i = i;
    out.total = Semiring::zero();
    out.best_term = Semiring::zero();
    for (const auto &seed : seeds) {
      if (i < seed.left || i >= seed.k) continue;
      const Value term = Semiring::times(
          chart.get(deductions::ItemKey{
              grammar::NonTerminal::VP, seed.k, j, -1, -1}),
          Semiring::from_weight(weights::local_weight(ctx, ctx.vpl_split_vp_weight(seed.k - i))));
      out.total = Semiring::plus(out.total, term);
      notify_run_term(grammar::RuleId::VPL_SPLIT_VP, i, seed.k - 1, term);
      if (Semiring::better(term, out.best_term)) {
        out.best_term = term;
        out.best_k = seed.k;
      }
    }
    return out;
  };
  struct Cursor {
    cand_pos_t i;
    cand_pos_t end;
    cand_pos_t step;
    VplStartPrefixSelectedCell<Value> cell;
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
          Semiring::from_weight(weights::local_weight(ctx, ctx.vpl_split_vp_weight(1)))) <
      Semiring::priority(Semiring::one());
  for (std::size_t cut = 0; cut + 1 < cuts.size(); ++cut) {
    const cand_pos_t lo = cuts[cut];
    const cand_pos_t hi = cuts[cut + 1] - 1;
    if (lo > hi) continue;
    const cand_pos_t first = toward_left ? hi : lo;
    const cand_pos_t end = toward_left ? lo : hi;
    const cand_pos_t step = toward_left ? -1 : 1;
    const auto cell = evaluate(first);
    if (!Semiring::is_zero(cell.total))
      heap.push(Cursor{first, end, step, cell});
  }
  std::vector<VplStartPrefixSelectedCell<Value>> selected;
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

#endif

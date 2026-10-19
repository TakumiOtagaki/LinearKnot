#ifndef CPARTY_SCFG2_BEAM_WM_START_PREFIX_HH
#define CPARTY_SCFG2_BEAM_WM_START_PREFIX_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/run_term_observer.hh"
#include "scfg2/deductions/item_key.hh"
#include "scfg2/grammar/rules.hh"

#include <cstddef>
#include <algorithm>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

template <class Value>
struct WmStartPrefixCell {
  Value total{};
  Value best_term{};
  cand_pos_t best_k = -1;
  grammar::RuleId best_rule = grammar::RuleId::WM_START_V;
};

// Exact factorization of the two empty-left WM start families. The seed at i
// is a V/WMB branch beginning at i; an unpaired i extends every start to its
// right by one ML-base step.
template <class Semiring>
class WmStartPrefixColumn {
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
                  WmStartPrefixCell<Value>{Semiring::zero(),
                                           Semiring::zero(),
                                           -1,
                                           grammar::RuleId::WM_START_V});
    cand_pos_t last_unit_base = -1;
    for (cand_pos_t i = j; i >= 1; --i) {
      auto &cell = cells_[static_cast<std::size_t>(i)];
      const bool seed = seed_allows(i);
      const Value v_term = seed ? Semiring::times(
          chart.get(deductions::ItemKey{
              grammar::NonTerminal::V, i, j, -1, -1}),
          Semiring::from_weight(weights::local_weight(ctx, ctx.start_v_weight(i, j, i)))) : Semiring::zero();
      const Value wmb_term = seed ? Semiring::times(
          chart.get(deductions::ItemKey{
              grammar::NonTerminal::WMB, i, j, -1, -1}),
          Semiring::from_weight(weights::local_weight(ctx, ctx.start_wmb_weight(i, j, i)))) : Semiring::zero();
      cell.total = Semiring::plus(v_term, wmb_term);
      if (Semiring::better(wmb_term, v_term)) {
        cell.best_term = wmb_term;
        cell.best_rule = grammar::RuleId::WM_START_WMB;
      } else {
        cell.best_term = v_term;
        cell.best_rule = grammar::RuleId::WM_START_V;
      }
      if (!Semiring::is_zero(cell.best_term)) cell.best_k = i;

      if (i == j || !ctx.can_pair_left_span(i, i + 1)) {
        cell.total = semirings::checkpoint_value(cell.total);
        continue;
      }
      const auto &next = cells_[static_cast<std::size_t>(i + 1)];
      const Value extended =
          Semiring::times(
              Semiring::from_weight(weights::local_weight(ctx, ctx.extend_unpaired_weight())),
              next.total);
      cell.total = Semiring::plus(cell.total, extended);
      if (!Semiring::is_zero(extended)) {
        notify_left_run_unit(grammar::RuleId::WM_START_V, j, i, extended,
                             static_cast<const Value *>(nullptr),
                             last_unit_base == i + 1);
        last_unit_base = i;
      }
      const Value extended_best =
          Semiring::times(
              Semiring::from_weight(weights::local_weight(ctx, ctx.extend_unpaired_weight())),
              next.best_term);
      if (Semiring::better(extended_best, cell.best_term)) {
        cell.best_term = extended_best;
        cell.best_rule = next.best_rule;
        cell.best_k = next.best_k;
      }
      cell.total = semirings::checkpoint_value(cell.total);
    }
  }

  Value contribution(cand_pos_t i) const {
    if (i < 1 || static_cast<std::size_t>(i) >= cells_.size()) {
      return Semiring::zero();
    }
    return cells_[static_cast<std::size_t>(i)].total;
  }

  const WmStartPrefixCell<Value> &cell(cand_pos_t i) const {
    return cells_[static_cast<std::size_t>(i)];
  }

 private:
  std::vector<WmStartPrefixCell<Value>> cells_;
};

template <class Value>
struct WmStartPrefixSelectedCell {
  cand_pos_t i = -1;
  Value total{};
  Value best_term{};
  cand_pos_t best_k = -1;
  grammar::RuleId best_rule = grammar::RuleId::WM_START_V;
};

// Finite-width exact top selection over the retained V/WMB seeds. Active
// seeds are constant between their left-span endpoints; within one such
// segment every term changes by the same ML-base multiplier, hence the
// aggregate is monotone and can be enumerated from one endpoint.
template <class Semiring, class Chart, class Context>
std::vector<WmStartPrefixSelectedCell<typename Semiring::value_type>>
select_wm_start_prefix_top_width(
    cand_pos_t j,
    const Chart &chart,
    const Context &ctx,
    const std::vector<cand_pos_t> &v_starts,
    const std::vector<cand_pos_t> &wmb_starts,
    std::size_t width,
    const std::vector<cand_pos_t> *forced_starts = nullptr) {
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
    const cand_pos_t left =
        std::max<cand_pos_t>(1, k - ctx.unpaired_prefix(k - 1));
    seeds.push_back(Seed{k, left, v, wmb});
    cuts.push_back(left);
    cuts.push_back(k + 1);
  }
  if (seeds.empty() || width == 0) return {};
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

  const auto evaluate = [&](cand_pos_t i) {
    WmStartPrefixSelectedCell<Value> out;
    out.i = i;
    out.total = Semiring::zero();
    out.best_term = Semiring::zero();
    for (const auto &seed : seeds) {
      if (i < seed.left || i > seed.k ||
          !ctx.can_pair_left_span(i, seed.k)) {
        continue;
      }
      const Value v_term = Semiring::times(
          seed.v,
          Semiring::from_weight(weights::local_weight(ctx, ctx.start_v_weight(i, j, seed.k))));
      const Value wmb_term = Semiring::times(
          seed.wmb,
          Semiring::from_weight(weights::local_weight(ctx, ctx.start_wmb_weight(i, j, seed.k))));
      out.total = Semiring::plus(out.total, v_term);
      out.total = Semiring::plus(out.total, wmb_term);
      notify_run_term(grammar::RuleId::WM_START_V, i, seed.k - 1, v_term);
      notify_run_term(grammar::RuleId::WM_START_WMB, i, seed.k - 1, wmb_term);
      if (Semiring::better(v_term, out.best_term)) {
        out.best_term = v_term;
        out.best_k = seed.k;
        out.best_rule = grammar::RuleId::WM_START_V;
      }
      if (Semiring::better(wmb_term, out.best_term)) {
        out.best_term = wmb_term;
        out.best_k = seed.k;
        out.best_rule = grammar::RuleId::WM_START_WMB;
      }
    }
    return out;
  };

  struct Cursor {
    cand_pos_t i = -1;
    cand_pos_t end = -1;
    cand_pos_t step = 0;
    WmStartPrefixSelectedCell<Value> cell{};
  };
  struct CursorWorse {
    bool operator()(const Cursor &lhs, const Cursor &rhs) const {
      const double lhs_priority = Semiring::priority(lhs.cell.total);
      const double rhs_priority = Semiring::priority(rhs.cell.total);
      if (lhs_priority != rhs_priority) return lhs_priority < rhs_priority;
      return lhs.i > rhs.i;
    }
  };
  std::priority_queue<Cursor, std::vector<Cursor>, CursorWorse> heap;
  const bool toward_left =
      Semiring::priority(
          Semiring::from_weight(weights::local_weight(ctx, ctx.extend_unpaired_weight()))) <
      Semiring::priority(Semiring::one());
  for (std::size_t cut = 0; cut + 1 < cuts.size(); ++cut) {
    cand_pos_t lo = cuts[cut];
    cand_pos_t hi = cuts[cut + 1] - 1;
    hi = std::min<cand_pos_t>(hi, j - 3);
    if (lo > hi) continue;
    const cand_pos_t first = toward_left ? hi : lo;
    const cand_pos_t end = toward_left ? lo : hi;
    const cand_pos_t step = toward_left ? -1 : 1;
    auto cell = evaluate(first);
    if (!Semiring::is_zero(cell.total)) {
      heap.push(Cursor{first, end, step, cell});
    }
  }

  std::vector<WmStartPrefixSelectedCell<Value>> selected;
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
  // Starts of the scaffold-only derivation are kept outside the width.
  if (forced_starts != nullptr) {
    for (const cand_pos_t i : *forced_starts) {
      const bool present = std::any_of(selected.begin(), selected.end(),
                                       [&](const auto &cell) { return cell.i == i; });
      if (present) continue;
      auto cell = evaluate(i);
      if (!Semiring::is_zero(cell.total)) selected.push_back(cell);
    }
  }
  return selected;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WM_START_PREFIX_HH

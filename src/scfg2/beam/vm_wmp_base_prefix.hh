#ifndef CPARTY_SCFG2_BEAM_VM_WMP_BASE_PREFIX_HH
#define CPARTY_SCFG2_BEAM_VM_WMP_BASE_PREFIX_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/run_term_observer.hh"
#include "scfg2/deductions/item_key.hh"

#include <cstddef>
#include <algorithm>
#include <queue>
#include <vector>

namespace cparty::scfg2::beam {

template <class Value>
struct VmWmpBasePrefixCell {
  Value tail{};
  cand_pos_t best_k = -1;
  Value best_term{};
};

// Exact private factor for the VM no-left-WM family. The recurrence shares
// the sum/max over all admissible WMp starts; the parent-dependent closing
// weight remains outside and is applied exactly once by contribution().
template <class Semiring>
class VmWmpBasePrefixColumn {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context>
  void build(cand_pos_t j, const Chart &chart, const Context &ctx) {
    j_ = j;
    cells_.assign(static_cast<std::size_t>(j + 2),
                  VmWmpBasePrefixCell<Value>{
                      Semiring::zero(), -1, Semiring::zero()});
    const cand_pos_t max_k = j - ctx.turn() - 1;
    cand_pos_t last_unit_base = -1;
    for (cand_pos_t i = max_k - 1; i >= 1; --i) {
      auto &cell = cells_[static_cast<std::size_t>(i)];
      const cand_pos_t seed_k = i + 1;
      const Value seed = chart.get(deductions::ItemKey{
          grammar::NonTerminal::WMp, seed_k, j - 1, -1, -1});
      cell.tail = seed;
      cell.best_term = seed;
      if (!Semiring::is_zero(seed)) cell.best_k = seed_k;

      if (!ctx.can_pair_left_span(i + 1, i + 2)) {
        cell.tail = semirings::checkpoint_value(cell.tail);
        continue;
      }
      const auto &next = cells_[static_cast<std::size_t>(i + 1)];
      const Value extended = Semiring::times(
          Semiring::from_weight(weights::local_weight(ctx, ctx.vm_base_unpaired_step_weight())),
          next.tail);
      cell.tail = Semiring::plus(cell.tail, extended);
      if (!Semiring::is_zero(extended)) {
        // Covers base i + 1 of the run [i+1, k-1].
        notify_left_run_unit(grammar::RuleId::VM_SPLIT_WMp_BASE, j, i + 1,
                             extended, static_cast<const Value *>(nullptr),
                             last_unit_base == i + 2);
        last_unit_base = i + 1;
      }
      const Value extended_best = Semiring::times(
          Semiring::from_weight(weights::local_weight(ctx, ctx.vm_base_unpaired_step_weight())),
          next.best_term);
      if (Semiring::better(extended_best, cell.best_term)) {
        cell.best_term = extended_best;
        cell.best_k = next.best_k;
      }
      cell.tail = semirings::checkpoint_value(cell.tail);
    }
  }

  template <class Context>
  Value contribution(cand_pos_t i, const Context &ctx) const {
    if (i < 1 || static_cast<std::size_t>(i) >= cells_.size()) {
      return Semiring::zero();
    }
    return Semiring::times(
        Semiring::from_weight(weights::local_weight(ctx, ctx.split_wm_wmp_weight(i, j_))),
        cells_[static_cast<std::size_t>(i)].tail);
  }

  const VmWmpBasePrefixCell<Value> &cell(cand_pos_t i) const {
    return cells_[static_cast<std::size_t>(i)];
  }

 private:
  cand_pos_t j_ = -1;
  std::vector<VmWmpBasePrefixCell<Value>> cells_;
};

template <class Value>
struct VmWmpBasePrefixSelectedCell {
  cand_pos_t i = -1;
  Value tail{};
  Value best_term{};
  cand_pos_t best_k = -1;
};

// Finite beam over the private tail before the parent-dependent VM closing
// weight. Retained WMp seeds define O(b) active intervals, and each interval
// is geometric, so exact top-width tail cells require no all-i scan.
template <class Semiring, class Chart, class Context>
std::vector<VmWmpBasePrefixSelectedCell<typename Semiring::value_type>>
select_vm_wmp_base_prefix_top_width(
    cand_pos_t j,
    const Chart &chart,
    const Context &ctx,
    const std::vector<cand_pos_t> &wmp_starts,
    std::size_t width) {
  using Value = typename Semiring::value_type;
  struct Seed {
    cand_pos_t k = -1;
    cand_pos_t left_parent = -1;
    Value child{};
  };
  std::vector<Seed> seeds;
  std::vector<cand_pos_t> cuts;
  for (const cand_pos_t k : wmp_starts) {
    const Value child = chart.get(deductions::ItemKey{
        grammar::NonTerminal::WMp, k, j - 1, -1, -1});
    if (Semiring::is_zero(child) || k > j - ctx.turn() - 1) continue;
    const cand_pos_t left_parent = std::max<cand_pos_t>(
        1, k - ctx.unpaired_prefix(k - 1) - 1);
    seeds.push_back(Seed{k, left_parent, child});
    cuts.push_back(left_parent);
    cuts.push_back(k);
  }
  if (seeds.empty() || width == 0) return {};
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

  const auto evaluate = [&](cand_pos_t i) {
    VmWmpBasePrefixSelectedCell<Value> out;
    out.i = i;
    out.tail = Semiring::zero();
    out.best_term = Semiring::zero();
    for (const auto &seed : seeds) {
      if (i < seed.left_parent || i >= seed.k ||
          !ctx.can_pair_left_span(i + 1, seed.k)) {
        continue;
      }
      const Value term = Semiring::times(
          seed.child,
          Semiring::from_weight(weights::local_weight(ctx, ctx.vm_base_gap_weight(seed.k - i - 1))));
      out.tail = Semiring::plus(out.tail, term);
      notify_run_term(grammar::RuleId::VM_SPLIT_WMp_BASE, i + 1, seed.k - 1, term);
      if (Semiring::better(term, out.best_term)) {
        out.best_term = term;
        out.best_k = seed.k;
      }
    }
    return out;
  };

  struct Cursor {
    cand_pos_t i = -1;
    cand_pos_t end = -1;
    cand_pos_t step = 0;
    VmWmpBasePrefixSelectedCell<Value> cell{};
  };
  struct CursorWorse {
    bool operator()(const Cursor &lhs, const Cursor &rhs) const {
      const double lhs_priority = Semiring::priority(lhs.cell.tail);
      const double rhs_priority = Semiring::priority(rhs.cell.tail);
      if (lhs_priority != rhs_priority) return lhs_priority < rhs_priority;
      return lhs.i > rhs.i;
    }
  };
  std::priority_queue<Cursor, std::vector<Cursor>, CursorWorse> heap;
  const bool toward_left =
      Semiring::priority(
          Semiring::from_weight(weights::local_weight(ctx, ctx.vm_base_unpaired_step_weight()))) <
      Semiring::priority(Semiring::one());
  for (std::size_t cut = 0; cut + 1 < cuts.size(); ++cut) {
    const cand_pos_t lo = cuts[cut];
    const cand_pos_t hi = cuts[cut + 1] - 1;
    if (lo > hi) continue;
    const cand_pos_t first = toward_left ? hi : lo;
    const cand_pos_t end = toward_left ? lo : hi;
    const cand_pos_t step = toward_left ? -1 : 1;
    auto cell = evaluate(first);
    if (!Semiring::is_zero(cell.tail)) {
      heap.push(Cursor{first, end, step, cell});
    }
  }

  std::vector<VmWmpBasePrefixSelectedCell<Value>> selected;
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
    if (!Semiring::is_zero(cursor.cell.tail)) heap.push(cursor);
  }
  return selected;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VM_WMP_BASE_PREFIX_HH

#ifndef CPARTY_SCFG2_ENGINE_OUTSIDE_HH
#define CPARTY_SCFG2_ENGINE_OUTSIDE_HH

#include "scfg2/engine/chart.hh"
#include "scfg2/deductions/deduction.hh"

#include <array>
#include <vector>

namespace cparty::scfg2::engine {

struct DeductionChildItem {
  deductions::ItemKey operator()(const deductions::Deduction &deduction,
                                 std::size_t index) const {
    return deduction.children[index].item;
  }
};

struct IgnoreOutsideDeduction {
  template <class Value>
  void operator()(const deductions::Deduction &, Value) const {}
};

// Reverse one admitted sum/product deduction, using the same local weight as
// inside. Each child OCCURRENCE receives beta(parent) * local * its siblings.
// Prefix/suffix products avoid division, including when a child has value zero.
// Admission/retention and any boundary child aliases belong to the caller.
template <class Semiring, class InsideChart, class LocalWeightFn,
          class ChildItemFn = DeductionChildItem>
void propagate_outside_deduction(
    const deductions::Deduction &deduction,
    const InsideChart &inside,
    typename Semiring::value_type parent_outside,
    Chart<typename Semiring::value_type> *outside,
    LocalWeightFn &&local_weight,
    ChildItemFn child_item = {}) {
  if (!outside || Semiring::is_zero(parent_outside) ||
      deduction.child_count == 0) return;
  using Value = typename Semiring::value_type;
  constexpr std::size_t max_children =
      std::tuple_size<decltype(deduction.children)>::value;
  std::array<Value, max_children + 1> prefix, suffix;
  std::array<Value, max_children> alpha;
  std::array<deductions::ItemKey, max_children> children;
  prefix[0] = Semiring::times(
      parent_outside, Semiring::from_weight(local_weight(deduction)));
  if (Semiring::is_zero(prefix[0])) return;
  for (std::size_t c = 0; c < deduction.child_count; ++c) {
    children[c] = child_item(deduction, c);
    alpha[c] = inside.get(children[c]);
    prefix[c + 1] = Semiring::times(prefix[c], alpha[c]);
  }
  suffix[deduction.child_count] = Semiring::one();
  for (std::size_t c = deduction.child_count; c > 0; --c) {
    suffix[c - 1] = Semiring::times(alpha[c - 1], suffix[c]);
  }
  for (std::size_t c = 0; c < deduction.child_count; ++c) {
    const auto contribution = Semiring::times(prefix[c], suffix[c + 1]);
    if (Semiring::is_zero(contribution)) continue;
    outside->set(children[c],
                 Semiring::plus(outside->get(children[c]), contribution));
  }
}

// A derived consumer of the shared deduction surface, with no saved edges.
// Preconditions: one immutable inside value per scheduled state; a topological
// forward schedule; regeneration returns precisely the branches admitted by
// inside, with the same child identities/values and local weights. In particular
// this does not reconstruct the old online runtime's pre-prune captured values.
template <class Semiring, class InsideChart, class DeductionGenerator,
          class LocalWeightFn, class ChildItemFn = DeductionChildItem,
          class DeductionObserver = IgnoreOutsideDeduction>
Chart<typename Semiring::value_type> run_regenerated_outside(
    const std::vector<deductions::ItemKey> &schedule,
    const InsideChart &inside,
    const deductions::ItemKey &root,
    DeductionGenerator &&generate,
    LocalWeightFn &&local_weight,
    ChildItemFn child_item = {},
    DeductionObserver observe = {},
    cand_pos_t dense_length = 0) {
  Chart<typename Semiring::value_type> outside(Semiring::zero(), nullptr, dense_length);
  outside.set(root, Semiring::one());
  for (auto cursor = schedule.rbegin(); cursor != schedule.rend(); ++cursor) {
    const auto beta = outside.get(*cursor);
    if (Semiring::is_zero(beta)) continue;
    for (const auto &deduction : generate(*cursor)) {
      // beta(parent) is final in reverse topological order. Observe leaves
      // too: they may emit whole unpaired runs even though they have no child
      // adjoints. Admission and zero local weights remain caller contracts.
      observe(deduction, beta);
      propagate_outside_deduction<Semiring>(
          deduction, inside, beta, &outside, local_weight, child_item);
    }
  }
  return outside;
}

}  // namespace cparty::scfg2::engine

#endif

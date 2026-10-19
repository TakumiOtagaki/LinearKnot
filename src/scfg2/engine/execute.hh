#ifndef CPARTY_SCFG2_ENGINE_EXECUTE_HH
#define CPARTY_SCFG2_ENGINE_EXECUTE_HH

#include "scfg2/deductions/deduction.hh"
#include "scfg2/engine/chart.hh"
#include "scfg2/engine/choice_map.hh"

#include <utility>
#include <vector>

namespace cparty::scfg2::engine {

template <class Semiring, class LocalWeightFn>
typename Semiring::value_type evaluate_deduction(
    const deductions::Deduction &deduction,
    const Chart<typename Semiring::value_type> &chart,
    LocalWeightFn &&local_weight) {
  using value_type = typename Semiring::value_type;
  value_type total = Semiring::from_weight(local_weight(deduction));
  for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
    total = Semiring::times(total, chart.get(deduction.children[idx].item));
  }
  return total;
}

template <class Semiring, class LocalWeightFn>
typename Semiring::value_type fold_deductions(
    const std::vector<deductions::Deduction> &deductions,
    const Chart<typename Semiring::value_type> &chart,
    LocalWeightFn &&local_weight) {
  using value_type = typename Semiring::value_type;
  value_type total = Semiring::zero();
  for (const auto &deduction : deductions) {
    total = Semiring::plus(
        total,
        evaluate_deduction<Semiring>(deduction, chart, std::forward<LocalWeightFn>(local_weight)));
  }
  return total;
}

template <class Semiring, class LocalWeightFn>
void run_cell_update(
    const deductions::ItemKey &parent,
    const std::vector<deductions::Deduction> &deductions,
    Chart<typename Semiring::value_type> *chart,
    LocalWeightFn &&local_weight) {
  if (!chart) {
    return;
  }
  chart->set(parent, fold_deductions<Semiring>(deductions, *chart, std::forward<LocalWeightFn>(local_weight)));
}

template <class Semiring, class LocalWeightFn>
void run_cell_update_with_choice(
    const deductions::ItemKey &parent,
    const std::vector<deductions::Deduction> &deductions,
    Chart<typename Semiring::value_type> *chart,
    ChoiceMap<typename Semiring::value_type> *choice_map,
    LocalWeightFn &&local_weight) {
  using value_type = typename Semiring::value_type;
  if (!chart) {
    return;
  }
  value_type total = Semiring::zero();
  for (const auto &deduction : deductions) {
    const value_type contribution =
        evaluate_deduction<Semiring>(deduction, *chart, local_weight);
    total = Semiring::plus(total, contribution);
    if (choice_map && !Semiring::is_zero(contribution)) {
      choice_map->observe(
          parent,
          deduction,
          contribution,
          [](value_type lhs, value_type rhs) { return Semiring::better(lhs, rhs); });
    }
  }
  chart->set(parent, total);
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_EXECUTE_HH

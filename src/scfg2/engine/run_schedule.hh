#ifndef CPARTY_SCFG2_ENGINE_RUN_SCHEDULE_HH
#define CPARTY_SCFG2_ENGINE_RUN_SCHEDULE_HH

#include "scfg2/engine/execute.hh"
#include "scfg2/engine/viterbi_result.hh"

#include <vector>

namespace cparty::scfg2::engine {

template <class Semiring, class DeductionProviderFn, class LocalWeightFn>
Chart<typename Semiring::value_type> run_schedule(
    const std::vector<deductions::ItemKey> &schedule,
    DeductionProviderFn &&deduction_provider,
    LocalWeightFn &&local_weight) {
  using value_type = typename Semiring::value_type;
  Chart<value_type> chart(Semiring::zero());
  for (const auto &item : schedule) {
    const auto deductions = deduction_provider(item);
    run_cell_update<Semiring>(item, deductions, &chart, local_weight);
  }
  return chart;
}

template <class Semiring, class DeductionProviderFn, class LocalWeightFn>
Chart<typename Semiring::value_type> run_schedule(
    const std::vector<deductions::ItemKey> &schedule,
    Chart<typename Semiring::value_type> chart,
    DeductionProviderFn &&deduction_provider,
    LocalWeightFn &&local_weight) {
  for (const auto &item : schedule) {
    const auto deductions = deduction_provider(item);
    run_cell_update<Semiring>(item, deductions, &chart, local_weight);
  }
  return chart;
}

template <class Semiring, class DeductionProviderFn, class LocalWeightFn>
ViterbiRunResult<typename Semiring::value_type> run_schedule_with_choices(
    const std::vector<deductions::ItemKey> &schedule,
    DeductionProviderFn &&deduction_provider,
    LocalWeightFn &&local_weight) {
  using value_type = typename Semiring::value_type;
  Chart<value_type> chart(Semiring::zero());
  ChoiceMap<value_type> choices;
  for (const auto &item : schedule) {
    const auto deductions = deduction_provider(item);
    run_cell_update_with_choice<Semiring>(item, deductions, &chart, &choices, local_weight);
  }
  return {std::move(chart), std::move(choices)};
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_RUN_SCHEDULE_HH

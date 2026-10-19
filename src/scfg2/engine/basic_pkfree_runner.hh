#ifndef CPARTY_SCFG2_ENGINE_BASIC_PKFREE_RUNNER_HH
#define CPARTY_SCFG2_ENGINE_BASIC_PKFREE_RUNNER_HH

#include "scfg2/deductions/basic_split_enumerator.hh"
#include "scfg2/deductions/generate_basic.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/engine/local_weight.hh"
#include "scfg2/engine/local_weight_model.hh"
#include "scfg2/engine/run_schedule.hh"
#include "scfg2/semirings/semiring.hh"

#include <vector>

namespace cparty::scfg2::engine {

inline std::vector<deductions::Deduction> generate_pkfree_basic_inside_deductions(
    const deductions::ItemKey &item) {
  using grammar::NonTerminal;
  using grammar::RuleId;

  if (item.nonterminal == NonTerminal::V) {
    return deductions::generate_basic_deductions(
        item, RuleId::V_HAIRPIN, deductions::enumerate_basic_splits(item, RuleId::V_HAIRPIN));
  }

  if (item.nonterminal == NonTerminal::W) {
    std::vector<deductions::Deduction> out;
    const auto unpaired = deductions::generate_basic_deductions(
        item,
        RuleId::W_EXTEND_UNPAIRED,
        deductions::enumerate_basic_splits(item, RuleId::W_EXTEND_UNPAIRED));
    out.insert(out.end(), unpaired.begin(), unpaired.end());
    const auto split_v = deductions::generate_basic_deductions(
        item,
        RuleId::W_SPLIT_V,
        deductions::enumerate_basic_splits(item, RuleId::W_SPLIT_V));
    out.insert(out.end(), split_v.begin(), split_v.end());
    return out;
  }

  return {};
}

template <class Semiring>
Chart<typename Semiring::value_type> run_pkfree_basic_inside(
    cand_pos_t n,
    const LocalWeightModel &weight_model) {
  const auto schedule = build_pkfree_basic_schedule(n);
  return run_schedule<Semiring>(
      schedule,
      [](const deductions::ItemKey &item) {
        return generate_pkfree_basic_inside_deductions(item);
      },
      [&](const deductions::Deduction &deduction) {
        return weight_model.weight_for(deduction);
      });
}

template <class Semiring>
Chart<typename Semiring::value_type> run_pkfree_basic_inside(cand_pos_t n) {
  ToyLocalWeightModel model;
  return run_pkfree_basic_inside<Semiring>(n, model);
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_BASIC_PKFREE_RUNNER_HH

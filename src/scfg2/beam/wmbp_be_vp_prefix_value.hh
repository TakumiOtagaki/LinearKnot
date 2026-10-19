#ifndef CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_VALUE_HH
#define CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_VALUE_HH

#include "scfg2/beam/wmbp_be_vp_prefix_profile.hh"

#include <optional>

namespace cparty::scfg2::beam {

// A transient product for ranking a generated WMBP source before its right VP
// attaches. No cache or checkpoint: discarded ranking candidates own no graph
// state. The ordinary deduction evaluator recomputes the selected full source.
template <class Semiring, class Chart, class LocalWeight>
std::optional<typename Semiring::value_type> evaluate_wmbp_be_vp_prefix(
    const deductions::Deduction &deduction,
    const Chart &chart,
    LocalWeight local_weight) {
  if (!is_wmbp_be_vp_rule(deduction.rule)) return std::nullopt;
  const auto key = wmbp_be_vp_prefix_key(deduction);
  auto value = Semiring::times(
      Semiring::from_weight(local_weight), chart.get(key.be));
  if (key.has_left) value = Semiring::times(value, chart.get(key.left));
  return value;
}

}  // namespace cparty::scfg2::beam

#endif

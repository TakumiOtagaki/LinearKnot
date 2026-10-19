#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_FACTOR_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_FACTOR_HH

#include "scfg2/beam/virtual_candidate.hh"
#include "scfg2/deductions/item_key.hh"
#include "scfg2/semirings/semiring.hh"

namespace cparty::scfg2::beam {

struct VpCase3LeftFactorValues {
  bool available = false;
  double left_child = 0.0;
  double private_tail = 0.0;
  double parent_facing_left = 0.0;
};

template <class Semiring, class Chart>
typename Semiring::value_type case3_wi_span_total(
    const Chart &chart, cand_pos_t left, cand_pos_t right) {
  if (left > right) return Semiring::one();
  return chart.get(deductions::ItemKey{
      grammar::NonTerminal::WI, left, right, -1, -1});
}

template <class Chart>
double case3_wi_span_total(const Chart &chart,
                           cand_pos_t left,
                           cand_pos_t right) {
  return case3_wi_span_total<semirings::SumProductSemiring>(
      chart, left, right);
}

// Algebraic factors for the disjoint LEFT_NONEMPTY family:
//
//   private_tail       = WI(middle) * WI(right)
//   parent_facing_left = WI(left) * private_tail
//
// The original VP_WI_CASE3 local weight remains excluded: it belongs to the
// final VP wrapper once the complete source split is known.
template <class Semiring, class Chart>
VpCase3LeftFactorValues vp_case3_left_factor_values(
    const deductions::Deduction &deduction,
    const Chart &chart) {
  if (vp_case3_partition(deduction) != VpCase3Partition::LeftNonempty) {
    return {};
  }
  const double left = case3_wi_span_total<Semiring>(
      chart, deduction.parent.i + 1, deduction.split.p - 1);
  const double middle = case3_wi_span_total<Semiring>(
      chart, deduction.split.q + 1, deduction.split.k - 1);
  const double right = case3_wi_span_total<Semiring>(
      chart, deduction.split.l + 1, deduction.parent.j - 1);
  const double tail = Semiring::times(middle, right);
  return VpCase3LeftFactorValues{
      true, left, tail, Semiring::times(left, tail)};
}

template <class Chart>
VpCase3LeftFactorValues vp_case3_left_factor_values(
    const deductions::Deduction &deduction,
    const Chart &chart) {
  return vp_case3_left_factor_values<semirings::SumProductSemiring>(
      deduction, chart);
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_FACTOR_HH

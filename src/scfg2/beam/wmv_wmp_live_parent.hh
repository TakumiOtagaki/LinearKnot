#ifndef CPARTY_SCFG2_BEAM_WMV_WMP_LIVE_PARENT_HH
#define CPARTY_SCFG2_BEAM_WMV_WMP_LIVE_PARENT_HH

#include "scfg2/deductions/generate_exact_basic.hh"

namespace cparty::scfg2::beam {

// Exact child-live omission for the two unary stem/extension families.
// False means every admissible deduction has a zero child contribution.
template <class Semiring, class Chart, class Context>
bool wmv_wmp_parent_has_live_child(
    const deductions::ItemKey &parent,
    const Chart &chart,
    const Context &ctx) {
  grammar::NonTerminal stem_child = grammar::NonTerminal::W;
  grammar::NonTerminal extension_child = grammar::NonTerminal::W;
  if (parent.nonterminal == grammar::NonTerminal::WMv) {
    stem_child = grammar::NonTerminal::V;
    extension_child = grammar::NonTerminal::WMv;
  } else if (parent.nonterminal == grammar::NonTerminal::WMp) {
    stem_child = grammar::NonTerminal::WMB;
    extension_child = grammar::NonTerminal::WMp;
  } else {
    return true;
  }

  if (deductions::exact_wmv_wmp_stem_parent_allows(parent, ctx) &&
      !Semiring::is_zero(chart.get(deductions::ItemKey{
          stem_child, parent.i, parent.j, -1, -1}))) {
    return true;
  }
  return deductions::exact_wmv_wmp_extend_parent_allows(parent, ctx) &&
         !Semiring::is_zero(chart.get(deductions::ItemKey{
             extension_child, parent.i, parent.j - 1, -1, -1}));
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMV_WMP_LIVE_PARENT_HH

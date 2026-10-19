#ifndef CPARTY_SCFG2_BEAM_WM_LIVE_PARENT_HH
#define CPARTY_SCFG2_BEAM_WM_LIVE_PARENT_HH

#include "scfg2/deductions/generate_exact_basic.hh"

namespace cparty::scfg2::beam {

// Exact child-live omission for WM. Start owns the empty-left case; recursive
// split is considered only for k > i and requires a live left WM prefix.
template <class Semiring, class Chart, class Context>
bool wm_parent_has_live_child(const deductions::ItemKey &parent,
                              const Chart &chart,
                              const Context &ctx) {
  if (parent.nonterminal != grammar::NonTerminal::WM) return true;
  if (deductions::exact_wm_extend_parent_allows(parent, ctx) &&
      !Semiring::is_zero(
          chart.get(deductions::ItemKey{grammar::NonTerminal::WM,
                                        parent.i,
                                        parent.j - 1,
                                        -1,
                                        -1}))) {
    return true;
  }
  if (!deductions::exact_wm_split_parent_allows(parent)) return false;

  for (cand_pos_t k = parent.i; k <= parent.j; ++k) {
    const bool right_live = !Semiring::is_zero(chart.get(
                                deductions::ItemKey{
                                    grammar::NonTerminal::V,
                                    k, parent.j, -1, -1})) ||
                            !Semiring::is_zero(chart.get(
                                deductions::ItemKey{
                                    grammar::NonTerminal::WMB,
                                    k, parent.j, -1, -1}));
    if (!right_live) continue;
    if (deductions::exact_wm_start_split_allows(parent, k, ctx)) {
      return true;
    }
    if (k > parent.i &&
        !Semiring::is_zero(
            chart.get(deductions::ItemKey{grammar::NonTerminal::WM,
                                          parent.i,
                                          k - 1,
                                          -1,
                                          -1}))) {
      return true;
    }
  }
  return false;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WM_LIVE_PARENT_HH

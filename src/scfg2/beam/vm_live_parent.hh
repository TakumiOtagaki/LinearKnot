#ifndef CPARTY_SCFG2_BEAM_VM_LIVE_PARENT_HH
#define CPARTY_SCFG2_BEAM_VM_LIVE_PARENT_HH

#include "scfg2/deductions/generate_exact_basic.hh"

namespace cparty::scfg2::beam {

// Exact child-live omission for VM. The no-left-WM case belongs exclusively
// to VM_SPLIT_WMp_BASE; both binary families begin at k=i+2.
template <class Semiring, class Chart, class Context>
bool vm_parent_has_live_child(const deductions::ItemKey &parent,
                              const Chart &chart,
                              const Context &ctx) {
  if (parent.nonterminal != grammar::NonTerminal::VM) return true;
  const cand_pos_t end_k = deductions::exact_vm_end_k(parent, ctx);

  for (cand_pos_t k = deductions::exact_vm_base_begin_k(parent);
       k <= end_k;
       ++k) {
    if (deductions::exact_vm_base_split_allows(parent, k, ctx) &&
        !Semiring::is_zero(
            chart.get(deductions::ItemKey{grammar::NonTerminal::WMp,
                                          k,
                                          parent.j - 1,
                                          -1,
                                          -1}))) {
      return true;
    }
  }
  for (cand_pos_t k = deductions::exact_vm_binary_begin_k(parent);
       k <= end_k;
       ++k) {
    if (Semiring::is_zero(
            chart.get(deductions::ItemKey{grammar::NonTerminal::WM,
                                          parent.i + 1,
                                          k - 1,
                                          -1,
                                          -1}))) {
      continue;
    }
    if (!Semiring::is_zero(
            chart.get(deductions::ItemKey{grammar::NonTerminal::WMv,
                                          k,
                                          parent.j - 1,
                                          -1,
                                          -1})) ||
        !Semiring::is_zero(
            chart.get(deductions::ItemKey{grammar::NonTerminal::WMp,
                                          k,
                                          parent.j - 1,
                                          -1,
                                          -1}))) {
      return true;
    }
  }
  return false;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VM_LIVE_PARENT_HH

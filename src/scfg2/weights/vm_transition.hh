#ifndef CPARTY_SCFG2_WEIGHTS_VM_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_VM_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_vm.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_vm(const deductions::Deduction &deduction,
                                   const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::VM_SPLIT_WM_WMv:
      return ctx.split_wm_wmv_weight(deduction.parent.i, deduction.parent.j);
    case grammar::RuleId::VM_SPLIT_WM_WMp:
      return ctx.split_wm_wmp_weight(deduction.parent.i, deduction.parent.j);
    case grammar::RuleId::VM_SPLIT_WMp_BASE:
      return ctx.split_wmp_base_weight(
          deduction.parent.i, deduction.parent.j, deduction.split.k - deduction.parent.i - 1);
    case grammar::RuleId::VM_SCALE2:
      return ctx.vm_wrapper_weight();
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_VM_TRANSITION_HH

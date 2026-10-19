#ifndef CPARTY_SCFG2_WEIGHTS_V_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_V_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_v.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_v(const deductions::Deduction &deduction,
                                  const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::V_HAIRPIN:
      return ctx.hairpin_weight(deduction.parent.i, deduction.parent.j);
    case grammar::RuleId::V_INTERNAL:
      if (deduction.child_count != 1 ||
          deduction.children[0].item.nonterminal != grammar::NonTerminal::V) {
        return impossible_weight<Context>();
      }
      return ctx.internal_weight(deduction.parent.i,
                                 deduction.parent.j,
                                 deduction.children[0].item.i,
                                 deduction.children[0].item.j);
    case grammar::RuleId::V_VM:
      return ctx.vm_closing_weight(deduction.parent.i, deduction.parent.j);
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_V_TRANSITION_HH

#ifndef CPARTY_SCFG2_WEIGHTS_BE_TRANSITION_HH
#define CPARTY_SCFG2_WEIGHTS_BE_TRANSITION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/deduction.hh"
#include "scfg2/weights/context_be.hh"

namespace cparty::scfg2::weights {

template <class Context>
inline double transition_weight_be(const deductions::Deduction &deduction,
                                   const Context &ctx) {
  switch (deduction.rule) {
    case grammar::RuleId::BE_BASE_SAMEPAIR:
      return ctx.be_base_samepair_weight();
    case grammar::RuleId::BE_STACK:
      return ctx.be_stack_weight(deduction.parent.i, deduction.parent.j);
    case grammar::RuleId::BE_INTERNAL_LOOP:
      if (deduction.child_count != 1) {
        return impossible_weight<Context>();
      }
      return ctx.be_internal_loop_weight(deduction.parent.i,
                                         deduction.parent.j,
                                         deduction.children[0].item.i,
                                         deduction.children[0].item.j);
    case grammar::RuleId::BE_WIP_WIP:
      return ctx.be_split_wip_wip_weight();
    case grammar::RuleId::BE_WIP_BASEPAIR:
      return ctx.be_split_wip_basepair_weight(deduction.parent.j - deduction.split.l - 1);
    case grammar::RuleId::BE_BASEPAIR_WIP:
      return ctx.be_split_basepair_wip_weight(deduction.split.k - deduction.parent.i - 1);
    default:
      return impossible_weight<Context>();
  }
}

}  // namespace cparty::scfg2::weights

#endif  // CPARTY_SCFG2_WEIGHTS_BE_TRANSITION_HH

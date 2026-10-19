#include "scfg2/deductions/expand_basic.hh"
#include "scfg2/weights/context_v.hh"
#include "scfg2/weights/v_transition.hh"

#include <cassert>
#include <cmath>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;
using cparty::scfg2::deductions::expand_basic_rule;
using cparty::scfg2::grammar::NonTerminal;
using cparty::scfg2::grammar::RuleId;
using cparty::scfg2::weights::VTransitionContext;

class FakeVTransitionContext final : public VTransitionContext {
 public:
  double hairpin_weight(cand_pos_t i, cand_pos_t j) const override {
    return static_cast<double>(i * 100 + j);
  }

  double internal_weight(cand_pos_t i,
                         cand_pos_t j,
                         cand_pos_t k,
                         cand_pos_t l) const override {
    return static_cast<double>(i * 1000 + j * 100 + k * 10 + l);
  }

  double vm_closing_weight(cand_pos_t i, cand_pos_t j) const override {
    return static_cast<double>(i * 10 + j);
  }
};

void test_hairpin_transition_weight() {
  FakeVTransitionContext ctx;
  const ItemKey parent{NonTerminal::V, 2, 9, -1, -1};
  const auto deduction = expand_basic_rule(parent, RuleId::V_HAIRPIN, {});
  assert(std::fabs(cparty::scfg2::weights::transition_weight_v(deduction, ctx) - 209.0) < 1e-12);
}

void test_internal_transition_weight() {
  FakeVTransitionContext ctx;
  Deduction deduction{};
  deduction.parent = ItemKey{NonTerminal::V, 2, 10, -1, -1};
  deduction.rule = RuleId::V_INTERNAL;
  deduction.child_count = 1;
  deduction.children[0].item = ItemKey{NonTerminal::V, 4, 8, -1, -1};
  assert(std::fabs(cparty::scfg2::weights::transition_weight_v(deduction, ctx) - 3048.0) < 1e-12);
}

void test_vm_transition_weight() {
  FakeVTransitionContext ctx;
  const ItemKey parent{NonTerminal::V, 3, 11, -1, -1};
  const auto deduction = expand_basic_rule(parent, RuleId::V_VM, {});
  assert(std::fabs(cparty::scfg2::weights::transition_weight_v(deduction, ctx) - 41.0) < 1e-12);
}

void test_invalid_internal_shape_returns_zero() {
  FakeVTransitionContext ctx;
  Deduction deduction{};
  deduction.parent = ItemKey{NonTerminal::V, 2, 10, -1, -1};
  deduction.rule = RuleId::V_INTERNAL;
  deduction.child_count = 1;
  deduction.children[0].item = ItemKey{NonTerminal::VM, 4, 8, -1, -1};
  assert(std::fabs(cparty::scfg2::weights::transition_weight_v(deduction, ctx)) < 1e-12);
}

}  // namespace

int main() {
  test_hairpin_transition_weight();
  test_internal_transition_weight();
  test_vm_transition_weight();
  test_invalid_internal_shape_returns_zero();
  return 0;
}

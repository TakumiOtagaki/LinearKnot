#include "scfg2/deductions/generate_basic.hh"
#include "scfg2/engine/basic_exact_weight.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

struct StubContext final : cparty::scfg2::engine::BasicExactWeightContext {
  double exterior_unpaired_weight() const override { return 2.0; }
  double exp_extloop(cand_pos_t i, cand_pos_t j) const override {
    return static_cast<double>(10 * i + j);
  }
  double hairpin_weight(cand_pos_t, cand_pos_t) const override { return 5.0; }
  double internal_weight(cand_pos_t, cand_pos_t, cand_pos_t, cand_pos_t) const override {
    return 7.0;
  }
  double vm_closing_weight(cand_pos_t, cand_pos_t) const override { return 11.0; }
  double split_wm_wmv_weight(cand_pos_t, cand_pos_t) const override { return 13.0; }
  double split_wm_wmp_weight(cand_pos_t, cand_pos_t) const override { return 17.0; }
  double split_wmp_base_weight(cand_pos_t, cand_pos_t, cand_pos_t) const override {
    return 19.0;
  }
  double vm_base_unpaired_step_weight() const override { return 1.5; }
  double vm_base_gap_weight(cand_pos_t left_unpaired) const override {
    return left_unpaired == 0 ? 1.0 : 1.5;
  }
  double vm_wrapper_weight() const override { return 23.0; }
  double wmv_stem_weight(cand_pos_t, cand_pos_t) const override { return 29.0; }
  double wmp_stem_weight(cand_pos_t, cand_pos_t) const override { return 31.0; }
  double wmv_extend_unpaired_weight() const override { return 37.0; }
  double wmp_extend_unpaired_weight() const override { return 41.0; }
  double start_v_weight(cand_pos_t, cand_pos_t, cand_pos_t) const override { return 43.0; }
  double start_wmb_weight(cand_pos_t, cand_pos_t, cand_pos_t) const override { return 47.0; }
  double split_v_weight(cand_pos_t, cand_pos_t) const override { return 53.0; }
  double split_wmb_weight(cand_pos_t, cand_pos_t) const override { return 59.0; }
  double extend_unpaired_weight() const override { return 61.0; }
};

}  // namespace

int main() {
  using namespace cparty::scfg2;
  StubContext ctx;
  engine::BasicExactWeightModel model(&ctx);

  const deductions::ItemKey w{grammar::NonTerminal::W, 1, 8, -1, -1};
  const auto unpaired = deductions::generate_basic_deductions_single(
      w, grammar::RuleId::W_EXTEND_UNPAIRED);
  assert(unpaired.size() == 1);
  assert(std::fabs(model.weight_for(unpaired.front()) - 2.0) < 1e-12);

  const deductions::ItemKey v{grammar::NonTerminal::V, 2, 7, -1, -1};
  const auto hairpin = deductions::generate_basic_deductions_single(
      v, grammar::RuleId::V_HAIRPIN);
  assert(hairpin.size() == 1);
  assert(std::fabs(model.weight_for(hairpin.front()) - 5.0) < 1e-12);

  const deductions::ItemKey vm{grammar::NonTerminal::VM, 2, 7, -1, -1};
  const auto scale2 = deductions::generate_basic_deductions_single(
      vm, grammar::RuleId::VM_SCALE2);
  assert(scale2.size() == 1);
  assert(std::fabs(model.weight_for(scale2.front()) - 23.0) < 1e-12);

  std::cout << "scfg2_rules_core_exact_weight=ok\n";
  return EXIT_SUCCESS;
}

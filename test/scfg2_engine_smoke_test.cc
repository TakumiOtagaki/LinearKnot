#include "scfg2/deductions/basic_split_enumerator.hh"
#include "scfg2/deductions/generate_basic.hh"
#include "scfg2/deductions/exact_context.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/basic_exact_weight.hh"
#include "scfg2/engine/basic_exact_inside_runner.hh"
#include "scfg2/engine/basic_pkfree_runner.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/engine/chart.hh"
#include "scfg2/engine/execute.hh"
#include "scfg2/engine/local_weight.hh"
#include "scfg2/engine/generated_target_index.hh"
#include "scfg2/engine/run_schedule.hh"
#include "scfg2/beam/wmv_wmp_live_parent.hh"
#include "scfg2/semirings/semiring.hh"
#include "scfg2/weights/vm_transition.hh"

#include <cassert>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

namespace {

using cparty::scfg2::deductions::ItemKey;
using cparty::scfg2::engine::Chart;
using cparty::scfg2::grammar::NonTerminal;
using cparty::scfg2::grammar::RuleId;
using cparty::scfg2::grammar::RuleSplit;

struct TestLogMaxSemiring {
  using value_type = double;
  static value_type zero() {
    return -std::numeric_limits<value_type>::infinity();
  }
  static value_type one() { return 0.0; }
  static value_type from_weight(double weight) { return std::log(weight); }
  static bool is_zero(value_type value) { return value == zero(); }
  static double priority(value_type value) { return value; }
  static value_type plus(value_type lhs, value_type rhs) {
    return std::max(lhs, rhs);
  }
  static value_type times(value_type lhs, value_type rhs) {
    return lhs + rhs;
  }
  static bool better(value_type lhs, value_type rhs) { return lhs > rhs; }
};

struct FakeBasicExactWeightContext : cparty::scfg2::engine::BasicExactWeightContext {
  double exterior_unpaired_value = 2.0;
  double hairpin_value = 5.0;
  double internal_factor_value = 17.0;
  double internal_scale_value = 3.0;
  double vm_value = 13.0;
  double vm_split_value = 19.0;
  double vm_wrapper_value = 17.0;
  double wmv_stem_value = 23.0;
  double wmp_stem_value = 24.0;
  double wmv_extend_value = 29.0;
  double wmp_extend_value = 30.0;
  double wm_start_value = 31.0;
  double wm_start_wmb_value = 32.0;
  double wm_split_value = 37.0;
  double wm_split_wmb_value = 38.0;
  double wm_extend_value = 41.0;

  double exterior_unpaired_weight() const override { return exterior_unpaired_value; }
  double exp_extloop(cand_pos_t i, cand_pos_t j) const override {
    return static_cast<double>(i * 10 + j);
  }
  double hairpin_weight(cand_pos_t, cand_pos_t) const override { return hairpin_value; }
  double internal_weight(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const override {
    (void)i;
    (void)j;
    (void)k;
    (void)l;
    return internal_factor_value * internal_scale_value;
  }
  double vm_closing_weight(cand_pos_t, cand_pos_t) const override { return vm_value; }
  double split_wm_wmv_weight(cand_pos_t, cand_pos_t) const override { return vm_split_value; }
  double split_wm_wmp_weight(cand_pos_t, cand_pos_t) const override { return vm_split_value + 1.0; }
  double split_wmp_base_weight(cand_pos_t, cand_pos_t, cand_pos_t left_unpaired) const override {
    return 50.0 + static_cast<double>(left_unpaired);
  }
  double vm_base_unpaired_step_weight() const override { return 1.0; }
  double vm_base_gap_weight(cand_pos_t) const override { return 1.0; }
  double vm_wrapper_weight() const override { return vm_wrapper_value; }
  double wmv_stem_weight(cand_pos_t, cand_pos_t) const override { return wmv_stem_value; }
  double wmp_stem_weight(cand_pos_t, cand_pos_t) const override { return wmp_stem_value; }
  double wmv_extend_unpaired_weight() const override { return wmv_extend_value; }
  double wmp_extend_unpaired_weight() const override { return wmp_extend_value; }
  double start_v_weight(cand_pos_t, cand_pos_t, cand_pos_t) const override { return wm_start_value; }
  double start_wmb_weight(cand_pos_t, cand_pos_t, cand_pos_t) const override { return wm_start_wmb_value; }
  double split_v_weight(cand_pos_t, cand_pos_t) const override { return wm_split_value; }
  double split_wmb_weight(cand_pos_t, cand_pos_t) const override { return wm_split_wmb_value; }
  double extend_unpaired_weight() const override { return wm_extend_value; }
};

struct FakeAllExactWeightContext : FakeBasicExactWeightContext,
                                   cparty::scfg2::weights::WTransitionContext,
                                   cparty::scfg2::weights::WITransitionContext,
                                   cparty::scfg2::weights::WIPTransitionContext,
                                   cparty::scfg2::weights::VPLTransitionContext,
                                   cparty::scfg2::weights::VPRTransitionContext,
                                   cparty::scfg2::weights::VPTransitionContext,
                                   cparty::scfg2::weights::WMBWTransitionContext,
                                   cparty::scfg2::weights::WMBPTransitionContext,
                                   cparty::scfg2::weights::WMBTransitionContext,
                                   cparty::scfg2::weights::BETransitionContext,
                                   cparty::scfg2::deductions::PairabilityContext {
  double w_extend_unpaired_weight() const override { return 2.0; }
  double w_split_v_weight(cand_pos_t, cand_pos_t) const override { return 3.0; }
  double w_split_wmb_weight() const override { return 0.0; }
  double w_root_wmb_decomposition_weight() const override { return 0.0; }

  double wi_base_single_weight() const override { return 19.0; }
  double wi_split_v_weight() const override { return 23.0; }
  double wi_split_wmb_weight() const override { return 0.0; }
  double wi_extend_unpaired_weight() const override { return 29.0; }

  double wip_base_v_weight() const override { return 31.0; }
  double wip_base_wmb_weight() const override { return 0.0; }
  double wip_split_v_weight() const override { return 37.0; }
  double wip_split_wmb_weight() const override { return 0.0; }
  double wip_basepair_v_weight(cand_pos_t) const override { return 41.0; }
  double wip_basepair_wmb_weight(cand_pos_t) const override { return 0.0; }
  double wip_extend_unpaired_weight() const override { return 43.0; }

  double vpl_split_vp_weight(cand_pos_t) const override { return 0.0; }
  double vpr_split_vp_wip_weight() const override { return 0.0; }
  double vpr_split_vp_basepair_weight(cand_pos_t) const override { return 0.0; }
  double vp_wi_case_weight() const override { return 0.0; }
  double vp_stack_weight(cand_pos_t, cand_pos_t) const override { return 0.0; }
  double vp_internal_loop_weight(cand_pos_t, cand_pos_t, cand_pos_t, cand_pos_t) const override {
    return 0.0;
  }
  double vp_split_weight() const override { return 0.0; }
  double wmbw_split_wmbp_wi_weight() const override { return 0.0; }
  double wmbp_split_be_wmbp_vp_weight(const cparty::scfg2::weights::WmbpSharedBeVpAnchors &) const override {
    return 0.0;
  }
  double wmbp_split_be_wmbw_vp_weight(const cparty::scfg2::weights::WmbpSharedBeVpAnchors &) const override {
    return 0.0;
  }
  double wmbp_direct_vp_weight(cand_pos_t, cand_pos_t) const override { return 0.0; }
  double wmbp_split_be_wi_vp_weight(const cparty::scfg2::deductions::WmbpWiVpAnchors &) const override {
    return 0.0;
  }
  double wmb_split_be_wmbp_wi_weight(
      const cparty::scfg2::weights::WmbSplitBeWmbpWiWeightAnchors &) const override {
    return 0.0;
  }
  double wmb_direct_wmbp_weight() const override { return 0.0; }
  double wmb_empty_weight() const override { return 0.0; }
  double be_base_samepair_weight() const override { return 0.0; }
  double be_stack_weight(cand_pos_t, cand_pos_t) const override { return 0.0; }
  double be_internal_loop_weight(cand_pos_t, cand_pos_t, cand_pos_t, cand_pos_t) const override {
    return 0.0;
  }
  double be_split_wip_wip_weight() const override { return 0.0; }
  double be_split_wip_basepair_weight(cand_pos_t) const override { return 0.0; }
  double be_split_basepair_wip_weight(cand_pos_t) const override { return 0.0; }

  cand_pos_t turn() const override { return 0; }
  cand_pos_t sequence_length() const override { return 3; }
  bool can_pair(cand_pos_t i, cand_pos_t j) const override { return i < j; }
  cand_pos_t unpaired_prefix(cand_pos_t pos) const override { return pos; }
  cand_pos_t pair_partner(cand_pos_t pos) const override {
    switch (pos) {
      case 1: return 3;
      case 3: return 1;
      default: return -1;
    }
  }
  cand_pos_t parent_index(cand_pos_t) const override { return -1; }
  cand_pos_t border_b(cand_pos_t, cand_pos_t) const override { return -1; }
  cand_pos_t border_bp(cand_pos_t, cand_pos_t) const override { return -1; }
  cand_pos_t border_B(cand_pos_t, cand_pos_t) const override { return -1; }
  cand_pos_t border_Bp(cand_pos_t, cand_pos_t) const override { return -1; }
  bool is_unpaired(cand_pos_t) const override { return true; }
};

void test_chart_default_and_set() {
  Chart<double> chart(0.0);
  const ItemKey key{NonTerminal::V, 1, 4, -1, -1};
  assert(chart.get(key) == 0.0);
  chart.set(key, 2.5);
  assert(chart.contains(key));
  assert(chart.get(key) == 2.5);
}

void test_sumproduct_update() {
  Chart<double> chart(0.0);
  const ItemKey w_prefix{NonTerminal::W, 1, 4, -1, -1};
  const ItemKey v{NonTerminal::V, 5, 12, -1, -1};
  chart.set(w_prefix, 2.0);
  chart.set(v, 4.0);

  const ItemKey parent{NonTerminal::W, 1, 12, -1, -1};
  RuleSplit split{};
  split.k = 5;
  const auto deductions =
      cparty::scfg2::deductions::generate_basic_deductions(parent, RuleId::W_SPLIT_V, {split});
  assert(deductions.size() == 1);
  assert(deductions[0].rule == RuleId::W_SPLIT_V);
  assert(deductions[0].child_count == 2);
  assert(deductions[0].children[0].item.i == 1);
  assert(deductions[0].children[0].item.j == 4);
  assert(deductions[0].children[1].item.i == 5);
  assert(deductions[0].children[1].item.j == 12);
  const double direct = cparty::scfg2::engine::evaluate_deduction<
      cparty::scfg2::semirings::SumProductSemiring>(deductions[0], chart, cparty::scfg2::engine::toy_local_weight);
  assert(std::fabs(direct - 24.0) < 1e-12);
  cparty::scfg2::engine::run_cell_update<cparty::scfg2::semirings::SumProductSemiring>(
      parent, deductions, &chart, cparty::scfg2::engine::toy_local_weight);
  assert(std::fabs(chart.get(parent) - 24.0) < 1e-12);
}

void test_maxproduct_fold() {
  Chart<double> chart(0.0);
  chart.set(ItemKey{NonTerminal::W, 1, 4, -1, -1}, 2.0);
  chart.set(ItemKey{NonTerminal::V, 5, 12, -1, -1}, 4.0);
  chart.set(ItemKey{NonTerminal::W, 1, 6, -1, -1}, 10.0);
  chart.set(ItemKey{NonTerminal::V, 7, 12, -1, -1}, 1.0);

  const ItemKey parent{NonTerminal::W, 1, 12, -1, -1};
  RuleSplit a{};
  a.k = 5;
  RuleSplit b{};
  b.k = 7;
  const auto deductions =
      cparty::scfg2::deductions::generate_basic_deductions(parent, RuleId::W_SPLIT_V, {a, b});
  const double total = cparty::scfg2::engine::fold_deductions<
      cparty::scfg2::semirings::MaxProductSemiring>(
      deductions, chart, cparty::scfg2::engine::toy_local_weight);
  assert(std::fabs(total - 30.0) < 1e-12);
}

void test_local_weight_enters_semiring_domain() {
  Chart<double> chart(TestLogMaxSemiring::zero());
  chart.set(ItemKey{NonTerminal::W, 1, 4, -1, -1}, std::log(2.0));
  chart.set(ItemKey{NonTerminal::V, 5, 12, -1, -1}, std::log(4.0));
  const ItemKey parent{NonTerminal::W, 1, 12, -1, -1};
  RuleSplit split{};
  split.k = 5;
  const auto deductions = cparty::scfg2::deductions::generate_basic_deductions(
      parent, RuleId::W_SPLIT_V, {split});
  const double contribution =
      cparty::scfg2::engine::evaluate_deduction<TestLogMaxSemiring>(
          deductions.front(), chart, [](const auto &) { return 3.0; });
  assert(std::fabs(contribution - std::log(24.0)) < 1e-12);
}

void test_live_parent_uses_semiring_zero() {
  FakeAllExactWeightContext ctx;
  Chart<double> chart(TestLogMaxSemiring::zero());
  const ItemKey parent{NonTerminal::WMv, 2, 6, -1, -1};
  // Log-domain one is numeric zero and must still count as a live child.
  chart.set(ItemKey{NonTerminal::V, 2, 6, -1, -1},
            TestLogMaxSemiring::one());
  assert(cparty::scfg2::beam::wmv_wmp_parent_has_live_child<
      TestLogMaxSemiring>(parent, chart, ctx));

  using Coupled =
      cparty::scfg2::semirings::CoupledLogSumMaxSemiring;
  Chart<Coupled::value_type> coupled_chart(Coupled::zero());
  coupled_chart.set(ItemKey{NonTerminal::V, 2, 6, -1, -1},
                    Coupled::one());
  assert(cparty::scfg2::beam::wmv_wmp_parent_has_live_child<Coupled>(
      parent, coupled_chart, ctx));
}

void test_targeted_coupled_semiring() {
  using Semiring =
      cparty::scfg2::semirings::TargetedCoupledLogSumMaxSemiring;
  const auto two = Semiring::from_weight(2.0);
  const auto three = Semiring::from_weight(3.0);
  const auto product = Semiring::times(two, three);
  assert(std::fabs(product.log_sum - std::log(6.0)) < 1e-12);
  assert(std::fabs(product.log_max - std::log(6.0)) < 1e-12);
  assert(std::fabs(product.log_target - std::log(6.0)) < 1e-12);

  const auto rejected = Semiring::reject_target(product);
  assert(std::fabs(rejected.log_sum - product.log_sum) < 1e-12);
  assert(std::fabs(rejected.log_max - product.log_max) < 1e-12);
  assert(!std::isfinite(rejected.log_target));
  const auto sum = Semiring::plus(rejected, two);
  assert(std::fabs(sum.log_sum - std::log(8.0)) < 1e-12);
  assert(std::fabs(sum.log_max - std::log(6.0)) < 1e-12);
  assert(std::fabs(sum.log_target - std::log(2.0)) < 1e-12);
}

void test_zero_child_cell_update() {
  Chart<double> chart(0.0);
  const ItemKey parent{NonTerminal::V, 4, 9, -1, -1};
  const auto deductions =
      cparty::scfg2::deductions::generate_basic_deductions_single(parent, RuleId::V_HAIRPIN);
  cparty::scfg2::engine::run_cell_update<cparty::scfg2::semirings::SumProductSemiring>(
      parent, deductions, &chart, cparty::scfg2::engine::toy_local_weight);
  assert(std::fabs(chart.get(parent) - 5.0) < 1e-12);
}

void test_toy_local_weight_surface() {
  const ItemKey parent{NonTerminal::WM, 2, 11, -1, -1};
  const auto deductions =
      cparty::scfg2::deductions::generate_basic_deductions_single(parent, RuleId::WM_EXTEND_UNPAIRED);
  assert(deductions.size() == 1);
  assert(std::fabs(cparty::scfg2::engine::toy_local_weight(deductions[0]) - 83.0) < 1e-12);

  const ItemKey pk_parent{NonTerminal::WMBP, 2, 11, -1, -1};
  const auto pk_deductions =
      cparty::scfg2::deductions::generate_basic_deductions_single(pk_parent, RuleId::WMBP_DIRECT_VP);
  assert(pk_deductions.size() == 1);
  assert(std::fabs(cparty::scfg2::engine::toy_local_weight(pk_deductions[0]) - 211.0) < 1e-12);
}

void test_toy_local_weight_explicit_rule_coverage() {
  using cparty::scfg2::engine::has_explicit_toy_local_weight;
  for (const auto rule : cparty::scfg2::grammar::all_rule_ids()) {
    assert(has_explicit_toy_local_weight(rule));
  }
}

void test_run_schedule_sum_and_max() {
  const ItemKey w14{NonTerminal::W, 1, 4, -1, -1};
  const ItemKey v512{NonTerminal::V, 5, 12, -1, -1};
  const ItemKey w112{NonTerminal::W, 1, 12, -1, -1};

  Chart<double> seeded_sum(0.0);
  seeded_sum.set(w14, 2.0);
  seeded_sum.set(v512, 4.0);

  const std::vector<ItemKey> schedule{w112};
  const auto provider = [&](const ItemKey &item) {
    if (item.nonterminal == NonTerminal::W && item.i == 1 && item.j == 12) {
      RuleSplit a{};
      a.k = 5;
      return cparty::scfg2::deductions::generate_basic_deductions(item, RuleId::W_SPLIT_V, {a});
    }
    return std::vector<cparty::scfg2::deductions::Deduction>{};
  };

  const auto sum_chart = cparty::scfg2::engine::run_schedule<
      cparty::scfg2::semirings::SumProductSemiring>(
      schedule, seeded_sum, provider, cparty::scfg2::engine::toy_local_weight);
  assert(std::fabs(sum_chart.get(w14) - 2.0) < 1e-12);
  assert(std::fabs(sum_chart.get(v512) - 4.0) < 1e-12);
  assert(std::fabs(sum_chart.get(w112) - 24.0) < 1e-12);

  const ItemKey w16{NonTerminal::W, 1, 6, -1, -1};
  const ItemKey v712{NonTerminal::V, 7, 12, -1, -1};
  Chart<double> seeded_max(0.0);
  seeded_max.set(w14, 2.0);
  seeded_max.set(v512, 4.0);
  seeded_max.set(w16, 10.0);
  seeded_max.set(v712, 1.0);

  const std::vector<ItemKey> schedule2{w112};
  const auto provider2 = [&](const ItemKey &item) {
    RuleSplit a{};
    a.k = 5;
    RuleSplit b{};
    b.k = 7;
    return cparty::scfg2::deductions::generate_basic_deductions(item, RuleId::W_SPLIT_V, {a, b});
  };
  auto seeded_max_chart = cparty::scfg2::engine::run_schedule<
      cparty::scfg2::semirings::MaxProductSemiring>(
      schedule2, seeded_max, provider2, cparty::scfg2::engine::toy_local_weight);
  assert(std::fabs(seeded_max_chart.get(w14) - 2.0) < 1e-12);
  assert(std::fabs(seeded_max_chart.get(w112) - 30.0) < 1e-12);
}

void test_basic_schedule_builder() {
  const auto schedule = cparty::scfg2::engine::build_pkfree_basic_schedule(3);
  assert(schedule.size() == 84);
  const auto find_index = [&](NonTerminal nt, cand_pos_t i, cand_pos_t j) {
    for (size_t idx = 0; idx < schedule.size(); ++idx) {
      if (schedule[idx].nonterminal == nt && schedule[idx].i == i && schedule[idx].j == j) {
        return idx;
      }
    }
    assert(false);
    return size_t{0};
  };

  assert(schedule[0].nonterminal == NonTerminal::VM);
  assert(schedule[0].i == 1);
  assert(schedule[0].j == 1);
  assert(find_index(NonTerminal::V, 1, 3) < find_index(NonTerminal::WMv, 1, 3));
  assert(find_index(NonTerminal::VP, 1, 3) < find_index(NonTerminal::WMBP, 1, 3));
  assert(find_index(NonTerminal::VP, 1, 3) < find_index(NonTerminal::WMBP, 1, 3));
  assert(find_index(NonTerminal::WMBP, 1, 3) < find_index(NonTerminal::WMBW, 1, 3));
  assert(find_index(NonTerminal::WMBP, 1, 3) < find_index(NonTerminal::WMB, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::WMp, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::WM, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::WIP, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::WI, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::W, 1, 3));
  assert(schedule.back().nonterminal == NonTerminal::W);
  assert(schedule.back().i == 1);
  assert(schedule.back().j == 3);
}

void test_exact_schedule_builder_dependency_order() {
  FakeAllExactWeightContext ctx{};
  const auto schedule = cparty::scfg2::engine::build_pkfree_basic_exact_schedule(3, ctx);
  const auto find_index = [&](NonTerminal nt,
                              cand_pos_t i,
                              cand_pos_t j,
                              cand_pos_t ip = -1,
                              cand_pos_t jp = -1) {
    for (size_t idx = 0; idx < schedule.size(); ++idx) {
      if (schedule[idx].nonterminal == nt &&
          schedule[idx].i == i &&
          schedule[idx].j == j &&
          schedule[idx].ip == ip &&
          schedule[idx].jp == jp) {
        return idx;
      }
    }
    assert(false);
    return size_t{0};
  };

  assert(find_index(NonTerminal::V, 1, 3) < find_index(NonTerminal::WMv, 1, 3));
  assert(find_index(NonTerminal::VP, 1, 3) < find_index(NonTerminal::BE, 1, 3, 1, 3));
  assert(find_index(NonTerminal::BE, 1, 3, 1, 3) < find_index(NonTerminal::WMBP, 1, 3));
  assert(find_index(NonTerminal::WMBP, 1, 3) < find_index(NonTerminal::WMBW, 1, 3));
  assert(find_index(NonTerminal::WMBP, 1, 3) < find_index(NonTerminal::WMB, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::WMp, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::WM, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::WIP, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::WI, 1, 3));
  assert(find_index(NonTerminal::WMB, 1, 3) < find_index(NonTerminal::W, 1, 3));
}

void test_basic_split_enumerator() {
  const ItemKey w{NonTerminal::W, 2, 4, -1, -1};
  const auto splits = cparty::scfg2::deductions::enumerate_basic_splits(w, RuleId::W_SPLIT_V);
  assert(splits.size() == 3);
  assert(splits[0].k == 2);
  assert(splits[2].k == 4);
}

void test_pkfree_basic_inside_runner() {
  const auto sum_chart =
      cparty::scfg2::engine::run_pkfree_basic_inside<cparty::scfg2::semirings::SumProductSemiring>(3);
  const ItemKey v11{NonTerminal::V, 1, 1, -1, -1};
  const ItemKey wi11{NonTerminal::WI, 1, 1, -1, -1};
  const ItemKey w11{NonTerminal::W, 1, 1, -1, -1};
  const ItemKey w13{NonTerminal::W, 1, 3, -1, -1};
  const ItemKey wip13{NonTerminal::WIP, 1, 3, -1, -1};
  const ItemKey vpl13{NonTerminal::VPL, 1, 3, -1, -1};
  const ItemKey vpr13{NonTerminal::VPR, 1, 3, -1, -1};
  const ItemKey vp13{NonTerminal::VP, 1, 3, -1, -1};
  const ItemKey wmp13{NonTerminal::WMp, 1, 3, -1, -1};
  const ItemKey wmbp13{NonTerminal::WMBP, 1, 3, -1, -1};
  const ItemKey wmb13{NonTerminal::WMB, 1, 3, -1, -1};
  assert(std::fabs(sum_chart.get(v11) - 5.0) < 1e-12);
  assert(sum_chart.get(wi11) >= 0.0);
  assert(sum_chart.get(w11) >= 0.0);
  assert(sum_chart.get(w13) >= 0.0);
  assert(sum_chart.get(wip13) >= 0.0);
  assert(sum_chart.get(vpl13) >= 0.0);
  assert(sum_chart.get(vpr13) >= 0.0);
  assert(sum_chart.get(vp13) >= 0.0);
  assert(sum_chart.get(wmp13) >= 0.0);
  assert(sum_chart.get(wmbp13) >= 0.0);
  assert(sum_chart.get(wmb13) >= 0.0);

  const auto max_chart =
      cparty::scfg2::engine::run_pkfree_basic_inside<cparty::scfg2::semirings::MaxProductSemiring>(3);
  assert(std::fabs(max_chart.get(v11) - 5.0) < 1e-12);
  assert(max_chart.get(wi11) >= 0.0);
  assert(max_chart.get(w11) >= 0.0);
  assert(max_chart.get(w13) >= 0.0);
  assert(max_chart.get(wip13) >= 0.0);
  assert(max_chart.get(vpl13) >= 0.0);
  assert(max_chart.get(vpr13) >= 0.0);
  assert(max_chart.get(vp13) >= 0.0);
  assert(max_chart.get(wmp13) >= 0.0);
  assert(max_chart.get(wmbp13) >= 0.0);
  assert(max_chart.get(wmb13) >= 0.0);
  assert(sum_chart.get(w13) >= max_chart.get(w13));
}

void test_exact_runner_with_all_local_weight_model() {
  FakeAllExactWeightContext ctx{};
  cparty::scfg2::engine::AllLocalWeightModel<FakeAllExactWeightContext> model(&ctx);
  const auto sum_chart =
      cparty::scfg2::engine::run_basic_exact_inside<cparty::scfg2::semirings::SumProductSemiring>(3, model, ctx);
  const ItemKey wi11{NonTerminal::WI, 1, 1, -1, -1};
  const ItemKey vm11{NonTerminal::VM, 1, 1, -1, -1};
  const ItemKey wip13{NonTerminal::WIP, 1, 3, -1, -1};
  const ItemKey vp13{NonTerminal::VP, 1, 3, -1, -1};
  const ItemKey be13{NonTerminal::BE, 1, 3, 1, 3};
  const ItemKey wmp13{NonTerminal::WMp, 1, 3, -1, -1};
  const ItemKey wmbp13{NonTerminal::WMBP, 1, 3, -1, -1};
  const ItemKey wmb13{NonTerminal::WMB, 1, 3, -1, -1};
  assert(std::fabs(sum_chart.get(wi11) - 19.0) < 1e-12);
  assert(std::fabs(sum_chart.get(vm11) - 0.0) < 1e-12);
  assert(sum_chart.get(wip13) > 0.0);
  assert(sum_chart.get(vp13) >= 0.0);
  assert(sum_chart.get(be13) >= 0.0);
  assert(sum_chart.get(wmp13) >= 0.0);
  assert(sum_chart.get(wmbp13) >= 0.0);
  assert(sum_chart.get(wmb13) >= 0.0);
}

void test_exact_viterbi_records_root_choice() {
  FakeAllExactWeightContext ctx{};
  cparty::scfg2::engine::AllLocalWeightModel<FakeAllExactWeightContext> model(&ctx);
  const auto result = cparty::scfg2::engine::run_basic_exact_viterbi(3, model, ctx);
  const ItemKey root{NonTerminal::W, 1, 3, -1, -1};
  const auto max_chart =
      cparty::scfg2::engine::run_basic_exact_inside<
          cparty::scfg2::semirings::MaxProductSemiring>(3, model, ctx);

  assert(std::fabs(result.chart.get(root) - max_chart.get(root)) < 1e-12);
  const auto *choice = result.choices.find(root);
  assert(choice);
  assert(choice->has_value);
  assert(std::fabs(choice->contribution - result.chart.get(root)) < 1e-12);
  assert(choice->deduction.parent.nonterminal == NonTerminal::W);
  assert(choice->deduction.parent.i == 1);
  assert(choice->deduction.parent.j == 3);
}

void test_deduction_output_pair_observation() {
  FakeAllExactWeightContext ctx{};
  const auto v = cparty::scfg2::deductions::expand_basic_rule(
      ItemKey{NonTerminal::V, 2, 8, -1, -1}, RuleId::V_HAIRPIN, {});
  const auto v_pairs = cparty::scfg2::engine::deduction_output_pairs(v, ctx);
  assert(v_pairs.size() == 1);
  assert(v_pairs[0].i == 2 && v_pairs[0].j == 8);
  assert(v_pairs[0].lane == cparty::scfg2::engine::TracebackPairLane::Round);

  const auto vp = cparty::scfg2::deductions::expand_basic_rule(
      ItemKey{NonTerminal::VP, 3, 10, -1, -1}, RuleId::VP_STACK, {});
  const auto vp_pairs = cparty::scfg2::engine::deduction_output_pairs(vp, ctx);
  assert(vp_pairs.size() == 1);
  assert(vp_pairs[0].i == 3 && vp_pairs[0].j == 10);
  assert(vp_pairs[0].lane == cparty::scfg2::engine::TracebackPairLane::Square);

  const auto be = cparty::scfg2::deductions::expand_basic_rule(
      ItemKey{NonTerminal::BE, 1, 12, 4, 9}, RuleId::BE_BASE_SAMEPAIR, {});
  const auto be_pairs = cparty::scfg2::engine::deduction_output_pairs(be, ctx);
  assert(be_pairs.size() == 2);
  assert(be_pairs[0].lane == cparty::scfg2::engine::TracebackPairLane::Round);
  assert(be_pairs[1].lane == cparty::scfg2::engine::TracebackPairLane::Round);
}

void test_exact_generated_target_partition_observation() {
  FakeAllExactWeightContext ctx{};
  cparty::scfg2::deductions::Deduction terminal_vp;
  terminal_vp.parent = ItemKey{NonTerminal::VP, 2, 8, -1, -1};
  terminal_vp.rule = RuleId::VP_WI_CASE1;
  const std::vector<cparty::scfg2::engine::TracebackPair> exact = {
      {2, 8, cparty::scfg2::engine::TracebackPairLane::Square}};
  const std::vector<cparty::scfg2::engine::TracebackPair> missing;
  assert(cparty::scfg2::engine::deduction_exactly_partitions_generated_target(
      terminal_vp, ctx, exact));
  assert(!cparty::scfg2::engine::deduction_exactly_partitions_generated_target(
      terminal_vp, ctx, missing));
  const std::vector<cparty::scfg2::engine::TracebackPair> strict_superset = {
      {2, 8, cparty::scfg2::engine::TracebackPairLane::Square},
      {3, 7, cparty::scfg2::engine::TracebackPairLane::Square},
  };
  assert(cparty::scfg2::engine::deduction_output_is_subset_of(
      terminal_vp, ctx, strict_superset));
  assert(!cparty::scfg2::engine::deduction_exactly_partitions_generated_target(
      terminal_vp, ctx, strict_superset));

  cparty::scfg2::deductions::Deduction split;
  split.parent = ItemKey{NonTerminal::W, 1, 10, -1, -1};
  split.rule = RuleId::W_SPLIT_V;
  split.child_count = 2;
  split.children[0].item = ItemKey{NonTerminal::W, 1, 1, -1, -1};
  split.children[1].item = ItemKey{NonTerminal::V, 2, 8, -1, -1};
  assert(cparty::scfg2::engine::deduction_exactly_partitions_generated_target(
      split, ctx, exact));
}

void test_exact_item_generator_vm_surface() {
  FakeAllExactWeightContext ctx{};
  const ItemKey vm{NonTerminal::VM, 1, 1, -1, -1};
  const auto deductions = cparty::scfg2::deductions::generate_exact_basic_deductions_for_item(vm, ctx);
  bool saw_scale2 = false;
  for (const auto &deduction : deductions) {
    if (deduction.rule == RuleId::VM_SCALE2) {
      saw_scale2 = true;
      break;
    }
  }
  assert(!saw_scale2);
}

void test_basic_exact_weight_model_surface() {
  FakeBasicExactWeightContext ctx{};
  cparty::scfg2::engine::BasicExactWeightModel model(&ctx);

  const ItemKey w{NonTerminal::W, 1, 12, -1, -1};
  RuleSplit w_split{};
  w_split.k = 5;
  const auto w_deductions =
      cparty::scfg2::deductions::generate_basic_deductions(w, RuleId::W_SPLIT_V, {w_split});
  assert(w_deductions.size() == 1);
  assert(std::fabs(model.weight_for(w_deductions[0]) - 62.0) < 1e-12);

  const ItemKey w_unpaired{NonTerminal::W, 1, 3, -1, -1};
  const auto w_unpaired_deductions =
      cparty::scfg2::deductions::generate_basic_deductions_single(w_unpaired, RuleId::W_EXTEND_UNPAIRED);
  assert(w_unpaired_deductions.size() == 1);
  assert(std::fabs(model.weight_for(w_unpaired_deductions[0]) - 2.0) < 1e-12);

  const ItemKey v{NonTerminal::V, 2, 8, -1, -1};
  const auto v_deductions =
      cparty::scfg2::deductions::generate_basic_deductions_single(v, RuleId::V_HAIRPIN);
  assert(v_deductions.size() == 1);
  assert(std::fabs(model.weight_for(v_deductions[0]) - 5.0) < 1e-12);

  RuleSplit v_internal_split{};
  v_internal_split.k = 4;
  v_internal_split.l = 6;
  const auto v_internal_deductions =
      cparty::scfg2::deductions::generate_basic_deductions(v, RuleId::V_INTERNAL, {v_internal_split});
  assert(v_internal_deductions.size() == 1);
  assert(std::fabs(model.weight_for(v_internal_deductions[0]) - 51.0) < 1e-12);

  const auto v_vm_deductions =
      cparty::scfg2::deductions::generate_basic_deductions_single(v, RuleId::V_VM);
  assert(v_vm_deductions.size() == 1);
  assert(std::fabs(model.weight_for(v_vm_deductions[0]) - 13.0) < 1e-12);

  const ItemKey vm{NonTerminal::VM, 2, 7, -1, -1};
  RuleSplit vm_split{};
  vm_split.k = 4;
  const auto vm_deductions =
      cparty::scfg2::deductions::generate_basic_deductions(vm, RuleId::VM_SPLIT_WM_WMv, {vm_split});
  assert(vm_deductions.size() == 1);
  assert(std::fabs(model.weight_for(vm_deductions[0]) - 19.0) < 1e-12);

  const auto vm_wmp =
      cparty::scfg2::deductions::generate_basic_deductions(vm, RuleId::VM_SPLIT_WM_WMp, {vm_split});
  assert(vm_wmp.size() == 1);
  assert(std::fabs(cparty::scfg2::weights::transition_weight_vm(vm_wmp[0], ctx) - 20.0) < 1e-12);
  assert(std::fabs(model.weight_for(vm_wmp[0]) - 20.0) < 1e-12);

  const auto vm_wmp_base =
      cparty::scfg2::deductions::generate_basic_deductions(vm, RuleId::VM_SPLIT_WMp_BASE, {vm_split});
  assert(vm_wmp_base.size() == 1);
  assert(std::fabs(model.weight_for(vm_wmp_base[0]) - 51.0) < 1e-12);

  const auto vm_scale2 =
      cparty::scfg2::deductions::generate_basic_deductions_single(vm, RuleId::VM_SCALE2);
  assert(vm_scale2.size() == 1);
  assert(std::fabs(model.weight_for(vm_scale2[0]) - 17.0) < 1e-12);

  const ItemKey wmv{NonTerminal::WMv, 3, 7, -1, -1};
  const auto wmv_stem =
      cparty::scfg2::deductions::generate_basic_deductions_single(wmv, RuleId::WMv_STEM_V);
  assert(wmv_stem.size() == 1);
  assert(std::fabs(model.weight_for(wmv_stem[0]) - 23.0) < 1e-12);

  const auto wmv_extend =
      cparty::scfg2::deductions::generate_basic_deductions_single(wmv, RuleId::WMv_EXTEND_UNPAIRED);
  assert(wmv_extend.size() == 1);
  assert(std::fabs(model.weight_for(wmv_extend[0]) - 29.0) < 1e-12);

  const ItemKey wm{NonTerminal::WM, 2, 7, -1, -1};
  RuleSplit wm_split{};
  wm_split.k = 4;
  const auto wm_start =
      cparty::scfg2::deductions::generate_basic_deductions(wm, RuleId::WM_START_V, {wm_split});
  assert(wm_start.size() == 1);
  assert(std::fabs(model.weight_for(wm_start[0]) - 31.0) < 1e-12);

  const auto wm_mid =
      cparty::scfg2::deductions::generate_basic_deductions(wm, RuleId::WM_SPLIT_V, {wm_split});
  assert(wm_mid.size() == 1);
  assert(std::fabs(model.weight_for(wm_mid[0]) - 37.0) < 1e-12);

  const auto wm_extend =
      cparty::scfg2::deductions::generate_basic_deductions_single(wm, RuleId::WM_EXTEND_UNPAIRED);
  assert(wm_extend.size() == 1);
  assert(std::fabs(model.weight_for(wm_extend[0]) - 41.0) < 1e-12);
}

void test_pkfree_basic_inside_runner_with_injected_weight_model() {
  FakeBasicExactWeightContext ctx{};
  cparty::scfg2::engine::BasicExactWeightModel model(&ctx);

  const auto sum_chart =
      cparty::scfg2::engine::run_pkfree_basic_inside<cparty::scfg2::semirings::SumProductSemiring>(3, model);
  const ItemKey v11{NonTerminal::V, 1, 1, -1, -1};
  const ItemKey w11{NonTerminal::W, 1, 1, -1, -1};
  const ItemKey w13{NonTerminal::W, 1, 3, -1, -1};
  assert(std::fabs(sum_chart.get(v11) - 5.0) < 1e-12);
  assert(sum_chart.get(w11) >= 0.0);
  assert(sum_chart.get(w13) >= 0.0);

  const auto max_chart =
      cparty::scfg2::engine::run_pkfree_basic_inside<cparty::scfg2::semirings::MaxProductSemiring>(3, model);
  assert(std::fabs(max_chart.get(v11) - 5.0) < 1e-12);
  assert(max_chart.get(w11) >= 0.0);
  assert(max_chart.get(w13) >= 0.0);
  assert(sum_chart.get(w13) >= max_chart.get(w13));
}

void test_bulk_pair_uniqueness_preserves_order_and_lanes() {
  using namespace cparty::scfg2::engine;
  std::mt19937 rng(20261008);
  for (int c = 0; c < 100; ++c) {
    std::vector<TracebackPair> occurrences;
    std::vector<TracebackPair> reference;
    for (int t = 0; t < 80; ++t) {
      const cand_pos_t i = 1 + rng() % 8;
      const cand_pos_t j = i + 1 + rng() % 8;
      const auto lane = rng() % 2 ? TracebackPairLane::Round : TracebackPairLane::Square;
      const TracebackPair pair{i, j, lane};
      occurrences.push_back(pair);
      append_unique_traceback_pair(&reference, pair);
    }
    stable_unique_traceback_pairs(&occurrences);
    assert(occurrences.size() == reference.size());
    for (std::size_t i = 0; i < reference.size(); ++i)
      assert(traceback_pair_equal(occurrences[i], reference[i]));
  }
  std::vector<TracebackPair> many;
  for (cand_pos_t i = 4096; i > 0; --i) {
    const TracebackPair pair{i, 8193 - i, TracebackPairLane::Square};
    many.push_back(pair);
    many.push_back(pair);
  }
  stable_unique_traceback_pairs(&many);
  assert(many.size() == 4096);
  for (std::size_t t = 0; t < many.size(); ++t)
    assert(many[t].i == cand_pos_t(4096 - t));
}

void test_indexed_target_partition_matches_set_oracle() {
  using namespace cparty::scfg2;
  FakeAllExactWeightContext ctx;
  std::mt19937 rng(20261008);
  const NonTerminal families[] = {NonTerminal::W, NonTerminal::V,
      NonTerminal::VP, NonTerminal::VM, NonTerminal::BE};
  const auto item = [&]() {
    const cand_pos_t i = 1 + rng() % 31;
    const cand_pos_t j = i + 1 + rng() % (32 - i);
    const cand_pos_t ip = i + rng() % (j - i + 1);
    const cand_pos_t jp = ip + rng() % (j - ip + 1);
    return ItemKey{families[rng() % 5], i, j, ip, jp};
  };
  for (int c = 0; c < 10000; ++c) {
    std::vector<engine::TracebackPair> raw;
    const int count = c % 5 == 0 ? 0 : rng() % 16;
    for (int p = 0; p < count; ++p) {
      const auto pair = item();
      raw.push_back({pair.i, pair.j, engine::TracebackPairLane::Square});
    }
    const auto target = engine::generated_traceback_pairs(raw, ctx);
    const engine::GeneratedTargetIndex index(target);
    deductions::Deduction deduction;
    deduction.parent = item();
    deduction.child_count = rng() % 4;
    for (std::size_t child = 0; child < deduction.child_count; ++child)
      deduction.children[child].item = item();
    assert(index.exactly_partitions(deduction, ctx) ==
        engine::deduction_exactly_partitions_generated_target(deduction, ctx, target));
  }
  // Every pair is covered by one of two disjoint children. A partition
  // query on a large target must visit logarithmically many tree nodes.
  for (std::size_t pairs : {32u, 256u, 2048u}) {
    std::vector<engine::TracebackPair> target;
    for (std::size_t p = 0; p < pairs; ++p)
      target.push_back({cand_pos_t(2*p+1), cand_pos_t(2*p+2), engine::TracebackPairLane::Square});
    const engine::GeneratedTargetIndex index(target);
    for (std::size_t c = 0; c < 200; ++c) {
      const cand_pos_t split = 2 * (c % pairs);
      deductions::Deduction d;
      d.parent = {NonTerminal::W, 1, cand_pos_t(2*pairs), -1, -1};
      d.child_count = 2;
      d.children[0].item = {NonTerminal::W, 1, split, -1, -1};
      d.children[1].item = {NonTerminal::W, split+1, cand_pos_t(2*pairs), -1, -1};
      assert(index.exactly_partitions(d, ctx));
    }
    std::size_t levels = 0;
    for (std::size_t size = pairs; size; size >>= 1) ++levels;
    assert(index.stats().interval_queries <= 200 * 12);
    assert(index.stats().range_node_visits <=
           index.stats().interval_queries * 8 * (levels + 1));
  }
}

void test_target_emissions_do_not_own_crossed_scaffold_interiors() {
  using namespace cparty::scfg2;
  FakeAllExactWeightContext ctx;
  deductions::Deduction d;
  d.parent = {NonTerminal::VP, 19, 67, -1, -1};
  d.rule = RuleId::VP_WI_CASE2;
  d.child_count = 2;
  d.children[0].item = {NonTerminal::WI, 20, 35, -1, -1};
  d.children[1].item = {NonTerminal::WI, 63, 66, -1, -1};
  const std::vector<engine::TracebackPair> target{
      {19, 67, engine::TracebackPairLane::Square},
      {49, 53, engine::TracebackPairLane::Square}};
  const engine::GeneratedTargetIndex index(target);
  assert(index.compatible_with_scaffold(ctx, 96));
  // The second target pair is emitted by the crossed scaffold's ordinary
  // interior, not by VP or either of its WI children. A bounding interval
  // overconstrains this otherwise valid local transition.
  assert(!index.exactly_partitions(d, ctx));
  assert(index.matches_local_emissions(d, ctx));
  assert(engine::deduction_locally_matches_generated_target(d, ctx, target));
  d.parent.i = 18;
  assert(!index.matches_local_emissions(d, ctx));
  d.parent = {NonTerminal::WI, 45, 49, -1, -1};
  d.rule = RuleId::WI_EXTEND_UNPAIRED;
  d.child_count = 1;
  assert(!index.matches_local_emissions(d, ctx));
  d.parent.j = 50;
  assert(index.matches_local_emissions(d, ctx));
  const engine::GeneratedTargetIndex invalid({
      {2, 8, engine::TracebackPairLane::Square},
      {2, 9, engine::TracebackPairLane::Square}});
  assert(!invalid.compatible_with_scaffold(ctx, 96));
}

}  // namespace

int main() {
  test_bulk_pair_uniqueness_preserves_order_and_lanes();
  test_indexed_target_partition_matches_set_oracle();
  test_target_emissions_do_not_own_crossed_scaffold_interiors();
  test_chart_default_and_set();
  test_sumproduct_update();
  test_maxproduct_fold();
  test_local_weight_enters_semiring_domain();
  test_live_parent_uses_semiring_zero();
  test_targeted_coupled_semiring();
  test_zero_child_cell_update();
  test_toy_local_weight_surface();
  test_toy_local_weight_explicit_rule_coverage();
  test_run_schedule_sum_and_max();
  test_basic_schedule_builder();
  test_exact_schedule_builder_dependency_order();
  test_basic_split_enumerator();
  test_pkfree_basic_inside_runner();
  test_basic_exact_weight_model_surface();
  test_pkfree_basic_inside_runner_with_injected_weight_model();
  test_exact_runner_with_all_local_weight_model();
  test_exact_viterbi_records_root_choice();
  test_deduction_output_pair_observation();
  test_exact_generated_target_partition_observation();
  test_exact_item_generator_vm_surface();
  return 0;
}

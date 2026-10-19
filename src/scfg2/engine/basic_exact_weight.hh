#ifndef CPARTY_SCFG2_ENGINE_BASIC_EXACT_WEIGHT_HH
#define CPARTY_SCFG2_ENGINE_BASIC_EXACT_WEIGHT_HH

#include "scfg2/engine/local_weight_model.hh"
#include "scfg2/weights/context_v.hh"
#include "scfg2/weights/context_vm.hh"
#include "scfg2/weights/context_wm.hh"
#include "scfg2/weights/context_wmv_wmp.hh"
#include "scfg2/weights/v_transition.hh"
#include "scfg2/weights/vm_transition.hh"
#include "scfg2/weights/wm_transition.hh"
#include "scfg2/weights/wmv_wmp_transition.hh"

namespace cparty::scfg2::engine {

class BasicExactWeightContext : public weights::VTransitionContext,
                                public weights::VMTransitionContext,
                                public weights::WMTransitionContext,
                                public weights::WMvWMpTransitionContext {
 public:
  virtual ~BasicExactWeightContext() = default;

  virtual double exterior_unpaired_weight() const = 0;
  virtual double exp_extloop(cand_pos_t i, cand_pos_t j) const = 0;
};

class BasicExactWeightModel : public LocalWeightModel {
 public:
  explicit BasicExactWeightModel(const BasicExactWeightContext *ctx) : ctx_(ctx) {}

  double weight_for(const deductions::Deduction &deduction) const override {
    if (!ctx_) {
      return 0.0;
    }
    switch (deduction.rule) {
      case grammar::RuleId::W_EXTEND_UNPAIRED:
        return ctx_->exterior_unpaired_weight();
      case grammar::RuleId::W_SPLIT_V:
        return ctx_->exp_extloop(deduction.split.k, deduction.parent.j);
      case grammar::RuleId::V_HAIRPIN:
      case grammar::RuleId::V_INTERNAL:
      case grammar::RuleId::V_VM:
        return weights::transition_weight_v(deduction, *ctx_);
      case grammar::RuleId::VM_SPLIT_WM_WMv:
      case grammar::RuleId::VM_SPLIT_WM_WMp:
      case grammar::RuleId::VM_SPLIT_WMp_BASE:
      case grammar::RuleId::VM_SCALE2:
        return weights::transition_weight_vm(deduction, *ctx_);
      case grammar::RuleId::WMv_STEM_V:
      case grammar::RuleId::WMp_STEM_WMB:
      case grammar::RuleId::WMv_EXTEND_UNPAIRED:
      case grammar::RuleId::WMp_EXTEND_UNPAIRED:
        return weights::transition_weight_wmv_wmp(deduction, *ctx_);
      case grammar::RuleId::WM_START_V:
      case grammar::RuleId::WM_START_WMB:
      case grammar::RuleId::WM_SPLIT_V:
      case grammar::RuleId::WM_SPLIT_WMB:
      case grammar::RuleId::WM_EXTEND_UNPAIRED:
        return weights::transition_weight_wm(deduction, *ctx_);
      default:
        return 0.0;
    }
  }

 private:
  const BasicExactWeightContext *ctx_ = nullptr;
};

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_BASIC_EXACT_WEIGHT_HH

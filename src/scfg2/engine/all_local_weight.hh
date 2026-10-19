#ifndef CPARTY_SCFG2_ENGINE_ALL_LOCAL_WEIGHT_HH
#define CPARTY_SCFG2_ENGINE_ALL_LOCAL_WEIGHT_HH

#include "scfg2/engine/local_weight_model.hh"
#include "scfg2/weights/be_transition.hh"
#include "scfg2/weights/v_transition.hh"
#include "scfg2/weights/vm_transition.hh"
#include "scfg2/weights/vp_transition.hh"
#include "scfg2/weights/vpl_transition.hh"
#include "scfg2/weights/vpr_transition.hh"
#include "scfg2/weights/w_transition.hh"
#include "scfg2/weights/wi_transition.hh"
#include "scfg2/weights/wip_transition.hh"
#include "scfg2/weights/wm_transition.hh"
#include "scfg2/weights/wmb_transition.hh"
#include "scfg2/weights/wmbp_transition.hh"
#include "scfg2/weights/wmbw_transition.hh"
#include "scfg2/weights/wmv_wmp_transition.hh"

namespace cparty::scfg2::engine {

template <class Context>
class AllLocalWeightModel : public LocalWeightModel {
 public:
  explicit AllLocalWeightModel(const Context *ctx) : ctx_(ctx) {}

  auto local_weight_for(const deductions::Deduction &deduction) const {
    if constexpr (weights::LogContext<Context>::value)
      return weights::LogWeight{context_weight_for(deduction)};
    else return context_weight_for(deduction);
  }

  double weight_for(const deductions::Deduction &deduction) const override {
    return weights::raw_weight(local_weight_for(deduction));
  }

  double log_weight_for(const deductions::Deduction &deduction) const {
    return weights::log_weight(local_weight_for(deduction));
  }

  double context_weight_for(const deductions::Deduction &deduction) const {
    if (!ctx_) {
      return weights::impossible_weight<Context>();
    }
    switch (deduction.rule) {
      case grammar::RuleId::W_EXTEND_UNPAIRED:
      case grammar::RuleId::W_SPLIT_V:
      case grammar::RuleId::W_SPLIT_WMB:
      case grammar::RuleId::W_ROOT_WMB_DECOMPOSITION:
        return weights::transition_weight_w(deduction, *ctx_);
      case grammar::RuleId::V_HAIRPIN:
      case grammar::RuleId::V_INTERNAL:
      case grammar::RuleId::V_VM:
        return weights::transition_weight_v(deduction, *ctx_);
      case grammar::RuleId::WI_BASE_SINGLE:
      case grammar::RuleId::WI_SPLIT_V:
      case grammar::RuleId::WI_SPLIT_WMB:
      case grammar::RuleId::WI_EXTEND_UNPAIRED:
        return weights::transition_weight_wi(deduction, *ctx_);
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
      case grammar::RuleId::WIP_BASE_V:
      case grammar::RuleId::WIP_BASE_WMB:
      case grammar::RuleId::WIP_SPLIT_V:
      case grammar::RuleId::WIP_SPLIT_WMB:
      case grammar::RuleId::WIP_BASEPAIR_V:
      case grammar::RuleId::WIP_BASEPAIR_WMB:
      case grammar::RuleId::WIP_EXTEND_UNPAIRED:
        return weights::transition_weight_wip(deduction, *ctx_);
      case grammar::RuleId::VPL_SPLIT_VP:
        return weights::transition_weight_vpl(deduction, *ctx_);
      case grammar::RuleId::VPR_SPLIT_VP_WIP:
      case grammar::RuleId::VPR_SPLIT_VP_BASEPAIR:
        return weights::transition_weight_vpr(deduction, *ctx_);
      case grammar::RuleId::VP_WI_CASE1:
      case grammar::RuleId::VP_WI_CASE2:
      case grammar::RuleId::VP_WI_CASE3:
      case grammar::RuleId::VP_STACK:
      case grammar::RuleId::VP_INTERNAL_LOOP:
      case grammar::RuleId::VP_WIP_VP_LEFT:
      case grammar::RuleId::VP_VP_WIP_RIGHT:
      case grammar::RuleId::VP_WIP_VPR:
      case grammar::RuleId::VP_VPL_WIP:
        return weights::transition_weight_vp(deduction, *ctx_);
      case grammar::RuleId::WMBW_SPLIT_WMBP_WI:
        return weights::transition_weight_wmbw(deduction, *ctx_);
      case grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP:
      case grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP:
      case grammar::RuleId::WMBP_DIRECT_VP:
      case grammar::RuleId::WMBP_SPLIT_BE_WI_VP:
        return weights::transition_weight_wmbp(deduction, *ctx_);
      case grammar::RuleId::WMB_SPLIT_BE_WMBP_WI:
      case grammar::RuleId::WMB_DIRECT_WMBP:
      case grammar::RuleId::WMB_EMPTY:
        return weights::transition_weight_wmb(deduction, *ctx_);
      case grammar::RuleId::BE_BASE_SAMEPAIR:
      case grammar::RuleId::BE_STACK:
      case grammar::RuleId::BE_INTERNAL_LOOP:
      case grammar::RuleId::BE_WIP_WIP:
      case grammar::RuleId::BE_WIP_BASEPAIR:
      case grammar::RuleId::BE_BASEPAIR_WIP:
        return weights::transition_weight_be(deduction, *ctx_);
      default:
        return weights::impossible_weight<Context>();
    }
  }

 private:
  const Context *ctx_ = nullptr;
};

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_ALL_LOCAL_WEIGHT_HH

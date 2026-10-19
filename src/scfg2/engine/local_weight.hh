#ifndef CPARTY_SCFG2_ENGINE_LOCAL_WEIGHT_HH
#define CPARTY_SCFG2_ENGINE_LOCAL_WEIGHT_HH

#include "scfg2/engine/local_weight_model.hh"

#include "scfg2/deductions/deduction.hh"

namespace cparty::scfg2::engine {

inline bool has_explicit_toy_local_weight(grammar::RuleId rule) {
  switch (rule) {
    case grammar::RuleId::V_HAIRPIN:
    case grammar::RuleId::V_INTERNAL:
    case grammar::RuleId::V_VM:
    case grammar::RuleId::W_EXTEND_UNPAIRED:
    case grammar::RuleId::W_SPLIT_V:
    case grammar::RuleId::W_SPLIT_WMB:
    case grammar::RuleId::WI_BASE_SINGLE:
    case grammar::RuleId::WI_SPLIT_V:
    case grammar::RuleId::WI_SPLIT_WMB:
    case grammar::RuleId::WI_EXTEND_UNPAIRED:
    case grammar::RuleId::VM_SPLIT_WM_WMv:
    case grammar::RuleId::VM_SPLIT_WM_WMp:
    case grammar::RuleId::VM_SPLIT_WMp_BASE:
    case grammar::RuleId::VM_SCALE2:
    case grammar::RuleId::WMv_STEM_V:
    case grammar::RuleId::WMp_STEM_WMB:
    case grammar::RuleId::WMv_EXTEND_UNPAIRED:
    case grammar::RuleId::WMp_EXTEND_UNPAIRED:
    case grammar::RuleId::WM_START_V:
    case grammar::RuleId::WM_START_WMB:
    case grammar::RuleId::WM_SPLIT_V:
    case grammar::RuleId::WM_SPLIT_WMB:
    case grammar::RuleId::WM_EXTEND_UNPAIRED:
    case grammar::RuleId::W_ROOT_WMB_DECOMPOSITION:
    case grammar::RuleId::WIP_BASE_V:
    case grammar::RuleId::WIP_BASE_WMB:
    case grammar::RuleId::WIP_SPLIT_V:
    case grammar::RuleId::WIP_SPLIT_WMB:
    case grammar::RuleId::WIP_BASEPAIR_V:
    case grammar::RuleId::WIP_BASEPAIR_WMB:
    case grammar::RuleId::WIP_EXTEND_UNPAIRED:
    case grammar::RuleId::VPL_SPLIT_VP:
    case grammar::RuleId::VPR_SPLIT_VP_WIP:
    case grammar::RuleId::VPR_SPLIT_VP_BASEPAIR:
    case grammar::RuleId::VP_WI_CASE1:
    case grammar::RuleId::VP_WI_CASE2:
    case grammar::RuleId::VP_WI_CASE3:
    case grammar::RuleId::VP_STACK:
    case grammar::RuleId::VP_INTERNAL_LOOP:
    case grammar::RuleId::VP_WIP_VP_LEFT:
    case grammar::RuleId::VP_VP_WIP_RIGHT:
    case grammar::RuleId::VP_WIP_VPR:
    case grammar::RuleId::VP_VPL_WIP:
    case grammar::RuleId::WMBW_SPLIT_WMBP_WI:
    case grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP:
    case grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP:
    case grammar::RuleId::WMBP_DIRECT_VP:
    case grammar::RuleId::WMBP_SPLIT_BE_WI_VP:
    case grammar::RuleId::WMB_SPLIT_BE_WMBP_WI:
    case grammar::RuleId::WMB_DIRECT_WMBP:
    case grammar::RuleId::WMB_EMPTY:
    case grammar::RuleId::BE_BASE_SAMEPAIR:
    case grammar::RuleId::BE_STACK:
    case grammar::RuleId::BE_INTERNAL_LOOP:
    case grammar::RuleId::BE_WIP_WIP:
    case grammar::RuleId::BE_WIP_BASEPAIR:
    case grammar::RuleId::BE_BASEPAIR_WIP:
      return true;
    default:
      return false;
  }
}

struct ToyLocalWeightModel : LocalWeightModel {
  double weight_for(const deductions::Deduction &deduction) const {
    switch (deduction.rule) {
      case grammar::RuleId::V_HAIRPIN:
        return 5.0;
      case grammar::RuleId::V_INTERNAL:
        return 11.0;
      case grammar::RuleId::V_VM:
        return 13.0;
      case grammar::RuleId::W_EXTEND_UNPAIRED:
        return 2.0;
      case grammar::RuleId::W_SPLIT_V:
        return 3.0;
      case grammar::RuleId::W_SPLIT_WMB:
        return 17.0;
      case grammar::RuleId::WI_BASE_SINGLE:
        return 19.0;
      case grammar::RuleId::WI_SPLIT_V:
        return 23.0;
      case grammar::RuleId::WI_SPLIT_WMB:
        return 27.0;
      case grammar::RuleId::WI_EXTEND_UNPAIRED:
        return 29.0;
      case grammar::RuleId::VM_SPLIT_WM_WMv:
        return 31.0;
      case grammar::RuleId::VM_SPLIT_WM_WMp:
        return 37.0;
      case grammar::RuleId::VM_SPLIT_WMp_BASE:
        return 41.0;
      case grammar::RuleId::VM_SCALE2:
        return 43.0;
      case grammar::RuleId::WMv_STEM_V:
        return 47.0;
      case grammar::RuleId::WMp_STEM_WMB:
        return 53.0;
      case grammar::RuleId::WMv_EXTEND_UNPAIRED:
        return 59.0;
      case grammar::RuleId::WMp_EXTEND_UNPAIRED:
        return 61.0;
      case grammar::RuleId::WM_START_V:
        return 67.0;
      case grammar::RuleId::WM_START_WMB:
        return 71.0;
      case grammar::RuleId::WM_SPLIT_V:
        return 73.0;
      case grammar::RuleId::WM_SPLIT_WMB:
        return 79.0;
      case grammar::RuleId::WM_EXTEND_UNPAIRED:
        return 83.0;
      case grammar::RuleId::W_ROOT_WMB_DECOMPOSITION:
        return 89.0;
      case grammar::RuleId::WIP_BASE_V:
        return 97.0;
      case grammar::RuleId::WIP_BASE_WMB:
        return 101.0;
      case grammar::RuleId::WIP_SPLIT_V:
        return 103.0;
      case grammar::RuleId::WIP_SPLIT_WMB:
        return 107.0;
      case grammar::RuleId::WIP_BASEPAIR_V:
        return 109.0;
      case grammar::RuleId::WIP_BASEPAIR_WMB:
        return 113.0;
      case grammar::RuleId::WIP_EXTEND_UNPAIRED:
        return 127.0;
      case grammar::RuleId::VPL_SPLIT_VP:
        return 131.0;
      case grammar::RuleId::VPR_SPLIT_VP_WIP:
        return 137.0;
      case grammar::RuleId::VPR_SPLIT_VP_BASEPAIR:
        return 139.0;
      case grammar::RuleId::VP_WI_CASE1:
        return 149.0;
      case grammar::RuleId::VP_WI_CASE2:
        return 151.0;
      case grammar::RuleId::VP_WI_CASE3:
        return 157.0;
      case grammar::RuleId::VP_STACK:
        return 163.0;
      case grammar::RuleId::VP_INTERNAL_LOOP:
        return 167.0;
      case grammar::RuleId::VP_WIP_VP_LEFT:
        return 173.0;
      case grammar::RuleId::VP_VP_WIP_RIGHT:
        return 179.0;
      case grammar::RuleId::VP_WIP_VPR:
        return 181.0;
      case grammar::RuleId::VP_VPL_WIP:
        return 191.0;
      case grammar::RuleId::WMBW_SPLIT_WMBP_WI:
        return 193.0;
      case grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP:
        return 197.0;
      case grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP:
        return 199.0;
      case grammar::RuleId::WMBP_DIRECT_VP:
        return 211.0;
      case grammar::RuleId::WMBP_SPLIT_BE_WI_VP:
        return 223.0;
      case grammar::RuleId::WMB_SPLIT_BE_WMBP_WI:
        return 227.0;
      case grammar::RuleId::WMB_DIRECT_WMBP:
        return 229.0;
      case grammar::RuleId::WMB_EMPTY:
        return 233.0;
      case grammar::RuleId::BE_BASE_SAMEPAIR:
        return 239.0;
      case grammar::RuleId::BE_STACK:
        return 241.0;
      case grammar::RuleId::BE_INTERNAL_LOOP:
        return 251.0;
      case grammar::RuleId::BE_WIP_WIP:
        return 257.0;
      case grammar::RuleId::BE_WIP_BASEPAIR:
        return 263.0;
      case grammar::RuleId::BE_BASEPAIR_WIP:
        return 269.0;
      default:
        return 1.0;
    }
  }
};

inline double toy_local_weight(const deductions::Deduction &deduction) {
  return ToyLocalWeightModel{}.weight_for(deduction);
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_LOCAL_WEIGHT_HH

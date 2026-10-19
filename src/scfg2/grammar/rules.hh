#ifndef CPARTY_SCFG2_GRAMMAR_RULES_HH
#define CPARTY_SCFG2_GRAMMAR_RULES_HH

#include "scfg2/grammar/nonterminal.hh"
#include "scfg2/grammar/surface_kind.hh"

#include <array>
#include <string_view>

namespace cparty::scfg2::grammar {

struct RuleSplit {
  int k = -1;
  int l = -1;
  int p = -1;
  int q = -1;
};

struct RuleChild {
  NonTerminal nonterminal;
  bool uses_outer_span = false;
  bool uses_secondary_span = false;
};

enum class RuleId : unsigned short {
  V_HAIRPIN,
  V_INTERNAL,
  V_VM,
  WI_BASE_SINGLE,
  WI_SPLIT_V,
  WI_SPLIT_WMB,
  WI_EXTEND_UNPAIRED,
  W_EXTEND_UNPAIRED,
  W_SPLIT_V,
  W_SPLIT_WMB,
  W_ROOT_WMB_DECOMPOSITION,
  VM_SPLIT_WM_WMv,
  VM_SPLIT_WM_WMp,
  VM_SPLIT_WMp_BASE,
  VM_SCALE2,
  WMv_STEM_V,
  WMp_STEM_WMB,
  WMv_EXTEND_UNPAIRED,
  WMp_EXTEND_UNPAIRED,
  WM_START_V,
  WM_START_WMB,
  WM_SPLIT_V,
  WM_SPLIT_WMB,
  WM_EXTEND_UNPAIRED,
  WIP_BASE_V,
  WIP_BASE_WMB,
  WIP_SPLIT_V,
  WIP_SPLIT_WMB,
  WIP_BASEPAIR_V,
  WIP_BASEPAIR_WMB,
  WIP_EXTEND_UNPAIRED,
  VPL_SPLIT_VP,
  VPR_SPLIT_VP_WIP,
  VPR_SPLIT_VP_BASEPAIR,
  VP_WI_CASE1,
  VP_WI_CASE2,
  VP_WI_CASE3,
  VP_STACK,
  VP_INTERNAL_LOOP,
  VP_WIP_VP_LEFT,
  VP_VP_WIP_RIGHT,
  VP_WIP_VPR,
  VP_VPL_WIP,
  WMBW_SPLIT_WMBP_WI,
  WMBP_SPLIT_BE_WMBP_VP,
  WMBP_SPLIT_BE_WMBW_VP,
  WMBP_DIRECT_VP,
  WMBP_SPLIT_BE_WI_VP,
  WMB_SPLIT_BE_WMBP_WI,
  WMB_DIRECT_WMBP,
  WMB_EMPTY,
  BE_BASE_SAMEPAIR,
  BE_STACK,
  BE_INTERNAL_LOOP,
  BE_WIP_WIP,
  BE_WIP_BASEPAIR,
  BE_BASEPAIR_WIP,
};

struct RuleSpec {
  RuleId id;
  NonTerminal parent;
  std::array<RuleChild, 3> children{};
  unsigned char child_count = 0;
};

inline constexpr WmbpStepPrefixKind classify_wmbp_step_prefix_kind(
    RuleId rule) {
  switch (rule) {
    case RuleId::WMBP_DIRECT_VP:
      return WmbpStepPrefixKind::Empty;
    case RuleId::WMBP_SPLIT_BE_WMBP_VP:
      return WmbpStepPrefixKind::PrefixWmbp;
    case RuleId::WMBP_SPLIT_BE_WMBW_VP:
      return WmbpStepPrefixKind::PrefixTail;
    case RuleId::WMBP_SPLIT_BE_WI_VP:
      return WmbpStepPrefixKind::PrefixWi;
    default:
      return WmbpStepPrefixKind::Unknown;
  }
}

inline constexpr std::array<RuleId, 57> all_rule_ids() {
  return {{
      RuleId::V_HAIRPIN,
      RuleId::V_INTERNAL,
      RuleId::V_VM,
      RuleId::WI_BASE_SINGLE,
      RuleId::WI_SPLIT_V,
      RuleId::WI_SPLIT_WMB,
      RuleId::WI_EXTEND_UNPAIRED,
      RuleId::W_EXTEND_UNPAIRED,
      RuleId::W_SPLIT_V,
      RuleId::W_SPLIT_WMB,
      RuleId::W_ROOT_WMB_DECOMPOSITION,
      RuleId::VM_SPLIT_WM_WMv,
      RuleId::VM_SPLIT_WM_WMp,
      RuleId::VM_SPLIT_WMp_BASE,
      RuleId::VM_SCALE2,
      RuleId::WMv_STEM_V,
      RuleId::WMp_STEM_WMB,
      RuleId::WMv_EXTEND_UNPAIRED,
      RuleId::WMp_EXTEND_UNPAIRED,
      RuleId::WM_START_V,
      RuleId::WM_START_WMB,
      RuleId::WM_SPLIT_V,
      RuleId::WM_SPLIT_WMB,
      RuleId::WM_EXTEND_UNPAIRED,
      RuleId::WIP_BASE_V,
      RuleId::WIP_BASE_WMB,
      RuleId::WIP_SPLIT_V,
      RuleId::WIP_SPLIT_WMB,
      RuleId::WIP_BASEPAIR_V,
      RuleId::WIP_BASEPAIR_WMB,
      RuleId::WIP_EXTEND_UNPAIRED,
      RuleId::VPL_SPLIT_VP,
      RuleId::VPR_SPLIT_VP_WIP,
      RuleId::VPR_SPLIT_VP_BASEPAIR,
      RuleId::VP_WI_CASE1,
      RuleId::VP_WI_CASE2,
      RuleId::VP_WI_CASE3,
      RuleId::VP_STACK,
      RuleId::VP_INTERNAL_LOOP,
      RuleId::VP_WIP_VP_LEFT,
      RuleId::VP_VP_WIP_RIGHT,
      RuleId::VP_WIP_VPR,
      RuleId::VP_VPL_WIP,
      RuleId::WMBW_SPLIT_WMBP_WI,
      RuleId::WMBP_SPLIT_BE_WMBP_VP,
      RuleId::WMBP_SPLIT_BE_WMBW_VP,
      RuleId::WMBP_DIRECT_VP,
      RuleId::WMBP_SPLIT_BE_WI_VP,
      RuleId::WMB_SPLIT_BE_WMBP_WI,
      RuleId::WMB_DIRECT_WMBP,
      RuleId::WMB_EMPTY,
      RuleId::BE_BASE_SAMEPAIR,
      RuleId::BE_STACK,
      RuleId::BE_INTERNAL_LOOP,
      RuleId::BE_WIP_WIP,
      RuleId::BE_WIP_BASEPAIR,
      RuleId::BE_BASEPAIR_WIP,
  }};
}

inline constexpr std::string_view rule_name(RuleId rule) {
  switch (rule) {
    case RuleId::V_HAIRPIN: return "V_HAIRPIN";
    case RuleId::V_INTERNAL: return "V_INTERNAL";
    case RuleId::V_VM: return "V_VM";
    case RuleId::WI_BASE_SINGLE: return "WI_BASE_SINGLE";
    case RuleId::WI_SPLIT_V: return "WI_SPLIT_V";
    case RuleId::WI_SPLIT_WMB: return "WI_SPLIT_WMB";
    case RuleId::WI_EXTEND_UNPAIRED: return "WI_EXTEND_UNPAIRED";
    case RuleId::W_EXTEND_UNPAIRED: return "W_EXTEND_UNPAIRED";
    case RuleId::W_SPLIT_V: return "W_SPLIT_V";
    case RuleId::W_SPLIT_WMB: return "W_SPLIT_WMB";
    case RuleId::W_ROOT_WMB_DECOMPOSITION: return "W_ROOT_WMB_DECOMPOSITION";
    case RuleId::VM_SPLIT_WM_WMv: return "VM_SPLIT_WM_WMv";
    case RuleId::VM_SPLIT_WM_WMp: return "VM_SPLIT_WM_WMp";
    case RuleId::VM_SPLIT_WMp_BASE: return "VM_SPLIT_WMp_BASE";
    case RuleId::VM_SCALE2: return "VM_SCALE2";
    case RuleId::WMv_STEM_V: return "WMv_STEM_V";
    case RuleId::WMp_STEM_WMB: return "WMp_STEM_WMB";
    case RuleId::WMv_EXTEND_UNPAIRED: return "WMv_EXTEND_UNPAIRED";
    case RuleId::WMp_EXTEND_UNPAIRED: return "WMp_EXTEND_UNPAIRED";
    case RuleId::WM_START_V: return "WM_START_V";
    case RuleId::WM_START_WMB: return "WM_START_WMB";
    case RuleId::WM_SPLIT_V: return "WM_SPLIT_V";
    case RuleId::WM_SPLIT_WMB: return "WM_SPLIT_WMB";
    case RuleId::WM_EXTEND_UNPAIRED: return "WM_EXTEND_UNPAIRED";
    case RuleId::WIP_BASE_V: return "WIP_BASE_V";
    case RuleId::WIP_BASE_WMB: return "WIP_BASE_WMB";
    case RuleId::WIP_SPLIT_V: return "WIP_SPLIT_V";
    case RuleId::WIP_SPLIT_WMB: return "WIP_SPLIT_WMB";
    case RuleId::WIP_BASEPAIR_V: return "WIP_BASEPAIR_V";
    case RuleId::WIP_BASEPAIR_WMB: return "WIP_BASEPAIR_WMB";
    case RuleId::WIP_EXTEND_UNPAIRED: return "WIP_EXTEND_UNPAIRED";
    case RuleId::VPL_SPLIT_VP: return "VPL_SPLIT_VP";
    case RuleId::VPR_SPLIT_VP_WIP: return "VPR_SPLIT_VP_WIP";
    case RuleId::VPR_SPLIT_VP_BASEPAIR: return "VPR_SPLIT_VP_BASEPAIR";
    case RuleId::VP_WI_CASE1: return "VP_WI_CASE1";
    case RuleId::VP_WI_CASE2: return "VP_WI_CASE2";
    case RuleId::VP_WI_CASE3: return "VP_WI_CASE3";
    case RuleId::VP_STACK: return "VP_STACK";
    case RuleId::VP_INTERNAL_LOOP: return "VP_INTERNAL_LOOP";
    case RuleId::VP_WIP_VP_LEFT: return "VP_WIP_VP_LEFT";
    case RuleId::VP_VP_WIP_RIGHT: return "VP_VP_WIP_RIGHT";
    case RuleId::VP_WIP_VPR: return "VP_WIP_VPR";
    case RuleId::VP_VPL_WIP: return "VP_VPL_WIP";
    case RuleId::WMBW_SPLIT_WMBP_WI: return "WMBW_SPLIT_WMBP_WI";
    case RuleId::WMBP_SPLIT_BE_WMBP_VP: return "WMBP_SPLIT_BE_WMBP_VP";
    case RuleId::WMBP_SPLIT_BE_WMBW_VP: return "WMBP_SPLIT_BE_WMBW_VP";
    case RuleId::WMBP_DIRECT_VP: return "WMBP_DIRECT_VP";
    case RuleId::WMBP_SPLIT_BE_WI_VP: return "WMBP_SPLIT_BE_WI_VP";
    case RuleId::WMB_SPLIT_BE_WMBP_WI: return "WMB_SPLIT_BE_WMBP_WI";
    case RuleId::WMB_DIRECT_WMBP: return "WMB_DIRECT_WMBP";
    case RuleId::WMB_EMPTY: return "WMB_EMPTY";
    case RuleId::BE_BASE_SAMEPAIR: return "BE_BASE_SAMEPAIR";
    case RuleId::BE_STACK: return "BE_STACK";
    case RuleId::BE_INTERNAL_LOOP: return "BE_INTERNAL_LOOP";
    case RuleId::BE_WIP_WIP: return "BE_WIP_WIP";
    case RuleId::BE_WIP_BASEPAIR: return "BE_WIP_BASEPAIR";
    case RuleId::BE_BASEPAIR_WIP: return "BE_BASEPAIR_WIP";
  }
  return "UNKNOWN";
}

inline constexpr RuleSpec rule_spec(RuleId rule) {
  switch (rule) {
    case RuleId::V_HAIRPIN:
      return {rule, NonTerminal::V, {}, 0};
    case RuleId::V_INTERNAL:
      return {rule, NonTerminal::V, {{{NonTerminal::V, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::V_VM:
      return {rule, NonTerminal::V, {{{NonTerminal::VM, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WI_BASE_SINGLE:
      return {rule, NonTerminal::WI, {}, 0};
    case RuleId::WI_SPLIT_V:
      return {rule, NonTerminal::WI, {{{NonTerminal::WI, true, false}, {NonTerminal::V, false, false}}}, 2};
    case RuleId::WI_SPLIT_WMB:
      return {rule, NonTerminal::WI, {{{NonTerminal::WI, true, false}, {NonTerminal::WMB, false, false}}}, 2};
    case RuleId::WI_EXTEND_UNPAIRED:
      return {rule, NonTerminal::WI, {{{NonTerminal::WI, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::W_EXTEND_UNPAIRED:
      return {rule, NonTerminal::W, {{{NonTerminal::W, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::W_SPLIT_V:
      return {rule, NonTerminal::W, {{{NonTerminal::W, true, false}, {NonTerminal::V, false, false}}}, 2};
    case RuleId::W_SPLIT_WMB:
      return {rule, NonTerminal::W, {{{NonTerminal::W, true, false}, {NonTerminal::WMB, false, false}}}, 2};
    case RuleId::W_ROOT_WMB_DECOMPOSITION:
      return {rule,
              NonTerminal::W,
              {{{NonTerminal::W, true, false},
                {NonTerminal::WMB, false, false},
                {NonTerminal::W, false, false}}},
              3};
    case RuleId::VM_SPLIT_WM_WMv:
      return {rule, NonTerminal::VM, {{{NonTerminal::WM, true, false}, {NonTerminal::WMv, false, false}}}, 2};
    case RuleId::VM_SPLIT_WM_WMp:
      return {rule, NonTerminal::VM, {{{NonTerminal::WM, true, false}, {NonTerminal::WMp, false, false}}}, 2};
    case RuleId::VM_SPLIT_WMp_BASE:
      return {rule, NonTerminal::VM, {{{NonTerminal::WMp, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::VM_SCALE2:
      return {rule, NonTerminal::VM, {{{NonTerminal::VM, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WMv_STEM_V:
      return {rule, NonTerminal::WMv, {{{NonTerminal::V, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WMp_STEM_WMB:
      return {rule, NonTerminal::WMp, {{{NonTerminal::WMB, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WMv_EXTEND_UNPAIRED:
      return {rule, NonTerminal::WMv, {{{NonTerminal::WMv, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WMp_EXTEND_UNPAIRED:
      return {rule, NonTerminal::WMp, {{{NonTerminal::WMp, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WM_START_V:
      return {rule, NonTerminal::WM, {{{NonTerminal::V, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WM_START_WMB:
      return {rule, NonTerminal::WM, {{{NonTerminal::WMB, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WM_SPLIT_V:
      return {rule, NonTerminal::WM, {{{NonTerminal::WM, true, false}, {NonTerminal::V, false, false}}}, 2};
    case RuleId::WM_SPLIT_WMB:
      return {rule, NonTerminal::WM, {{{NonTerminal::WM, true, false}, {NonTerminal::WMB, false, false}}}, 2};
    case RuleId::WM_EXTEND_UNPAIRED:
      return {rule, NonTerminal::WM, {{{NonTerminal::WM, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WIP_BASE_V:
      return {rule, NonTerminal::WIP, {{{NonTerminal::V, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WIP_BASE_WMB:
      return {rule, NonTerminal::WIP, {{{NonTerminal::WMB, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WIP_SPLIT_V:
      return {rule, NonTerminal::WIP, {{{NonTerminal::V, true, false}, {NonTerminal::WIP, false, false}}}, 2};
    case RuleId::WIP_SPLIT_WMB:
      return {rule, NonTerminal::WIP, {{{NonTerminal::WMB, true, false}, {NonTerminal::WIP, false, false}}}, 2};
    case RuleId::WIP_BASEPAIR_V:
      return {rule, NonTerminal::WIP, {{{NonTerminal::V, true, true}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WIP_BASEPAIR_WMB:
      return {rule, NonTerminal::WIP, {{{NonTerminal::WMB, true, true}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WIP_EXTEND_UNPAIRED:
      return {rule, NonTerminal::WIP, {{{NonTerminal::WIP, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::VPL_SPLIT_VP:
      return {rule, NonTerminal::VPL, {{{NonTerminal::VP, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::VPR_SPLIT_VP_WIP:
      return {rule, NonTerminal::VPR, {{{NonTerminal::VP, true, false}, {NonTerminal::WIP, false, false}}}, 2};
    case RuleId::VPR_SPLIT_VP_BASEPAIR:
      return {rule, NonTerminal::VPR, {{{NonTerminal::VP, true, true}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::VP_WI_CASE1:
    case RuleId::VP_WI_CASE2:
      return {rule, NonTerminal::VP, {{{NonTerminal::WI, true, false}, {NonTerminal::WI, false, false}}}, 2};
    case RuleId::VP_WI_CASE3:
      return {rule,
              NonTerminal::VP,
              {{{NonTerminal::WI, true, false},
                {NonTerminal::WI, false, false},
                {NonTerminal::WI, false, false}}},
              3};
    case RuleId::VP_STACK:
    case RuleId::VP_INTERNAL_LOOP:
      return {rule, NonTerminal::VP, {{{NonTerminal::VP, true, true}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::VP_WIP_VP_LEFT:
      return {rule, NonTerminal::VP, {{{NonTerminal::WIP, true, false}, {NonTerminal::VP, false, true}}}, 2};
    case RuleId::VP_VP_WIP_RIGHT:
      return {rule, NonTerminal::VP, {{{NonTerminal::VP, true, false}, {NonTerminal::WIP, false, true}}}, 2};
    case RuleId::VP_WIP_VPR:
      return {rule, NonTerminal::VP, {{{NonTerminal::WIP, true, false}, {NonTerminal::VPR, false, true}}}, 2};
    case RuleId::VP_VPL_WIP:
      return {rule, NonTerminal::VP, {{{NonTerminal::VPL, true, false}, {NonTerminal::WIP, false, true}}}, 2};
    case RuleId::WMBW_SPLIT_WMBP_WI:
      return {rule, NonTerminal::WMBW, {{{NonTerminal::WMBP, true, false}, {NonTerminal::WI, false, false}}}, 2};
    case RuleId::WMBP_SPLIT_BE_WMBP_VP:
      return {rule,
              NonTerminal::WMBP,
              {{{NonTerminal::BE, true, true},
                {NonTerminal::WMBP, false, false},
                {NonTerminal::VP, false, false}}},
              3};
    case RuleId::WMBP_SPLIT_BE_WMBW_VP:
      return {rule,
              NonTerminal::WMBP,
              {{{NonTerminal::BE, true, true},
                {NonTerminal::WMBW, false, false},
                {NonTerminal::VP, false, false}}},
              3};
    case RuleId::WMBP_DIRECT_VP:
      return {rule,
              NonTerminal::WMBP,
              {{{NonTerminal::VP, true, false}, {NonTerminal::V, false, false}}},
              1};
    case RuleId::WMBP_SPLIT_BE_WI_VP:
      return {rule,
              NonTerminal::WMBP,
              {{{NonTerminal::BE, true, true},
                {NonTerminal::WI, false, false},
                {NonTerminal::VP, false, false}}},
              3};
    case RuleId::WMB_SPLIT_BE_WMBP_WI:
      return {rule,
              NonTerminal::WMB,
              {{{NonTerminal::BE, true, true},
                {NonTerminal::WMBP, false, false},
                {NonTerminal::WI, false, false}}},
              3};
    case RuleId::WMB_DIRECT_WMBP:
      return {rule, NonTerminal::WMB, {{{NonTerminal::WMBP, true, false}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::WMB_EMPTY:
      return {rule, NonTerminal::WMB, {}, 0};
    case RuleId::BE_BASE_SAMEPAIR:
      return {rule, NonTerminal::BE, {}, 0};
    case RuleId::BE_STACK:
    case RuleId::BE_INTERNAL_LOOP:
      return {rule, NonTerminal::BE, {{{NonTerminal::BE, true, true}, {NonTerminal::V, false, false}}}, 1};
    case RuleId::BE_WIP_WIP:
      return {rule,
              NonTerminal::BE,
              {{{NonTerminal::WIP, true, false},
                {NonTerminal::BE, false, true},
                {NonTerminal::WIP, false, true}}},
              3};
    case RuleId::BE_WIP_BASEPAIR:
      return {rule, NonTerminal::BE, {{{NonTerminal::WIP, true, false}, {NonTerminal::BE, false, true}}}, 2};
    case RuleId::BE_BASEPAIR_WIP:
      return {rule, NonTerminal::BE, {{{NonTerminal::BE, true, true}, {NonTerminal::WIP, false, true}}}, 2};
  }
  return {rule, NonTerminal::W, {}, 0};
}

}  // namespace cparty::scfg2::grammar

#endif  // CPARTY_SCFG2_GRAMMAR_RULES_HH

#ifndef CPARTY_SCFG2_DEDUCTIONS_EXPAND_BASIC_HH
#define CPARTY_SCFG2_DEDUCTIONS_EXPAND_BASIC_HH

#include "scfg2/deductions/deduction.hh"

#include <type_traits>

namespace cparty::scfg2::deductions {

inline constexpr ItemKey make_item(grammar::NonTerminal nt,
                                   cand_pos_t i,
                                   cand_pos_t j,
                                   cand_pos_t ip = -1,
                                   cand_pos_t jp = -1) {
  return ItemKey{nt, i, j, ip, jp};
}

inline constexpr Deduction make_zero_child_deduction(const ItemKey &parent,
                                                     grammar::RuleId rule,
                                                     const grammar::RuleSplit &split) {
  return Deduction{parent, rule, split, {}, 0};
}

namespace detail {

// Runtime and fixed-rule entry points share this one expansion body.
// An integral_constant selector lets the compiler specialize only requested
// rules without forcing the complete runtime switch into every caller.
template <class RuleSelector>
inline constexpr Deduction expand_basic_rule_impl(
    const ItemKey &parent,
    RuleSelector rule_selector,
    const grammar::RuleSplit &split) {
  const grammar::RuleId rule = rule_selector;
  Deduction out{};
  out.parent = parent;
  out.rule = rule;
  out.split = split;


  switch (rule) {
    case grammar::RuleId::V_HAIRPIN:
    case grammar::RuleId::WI_BASE_SINGLE:
    case grammar::RuleId::WMB_EMPTY:
    case grammar::RuleId::BE_BASE_SAMEPAIR:
      return make_zero_child_deduction(parent, rule, split);

    case grammar::RuleId::V_INTERNAL:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::V, split.k, split.l);
      return out;

    case grammar::RuleId::W_EXTEND_UNPAIRED:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::W, parent.i, parent.j - 1);
      return out;

    case grammar::RuleId::W_SPLIT_V:
      if (split.k > parent.i) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::W, parent.i, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      }
      return out;

    case grammar::RuleId::W_SPLIT_WMB:
      if (split.k > parent.i) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::W, parent.i, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      }
      return out;

    case grammar::RuleId::W_ROOT_WMB_DECOMPOSITION:
      out.child_count = 3;
      out.children[0].item = make_item(grammar::NonTerminal::W, parent.i, split.p - 1);
      out.children[1].item = make_item(grammar::NonTerminal::WMB, split.p, split.k);
      out.children[2].item = make_item(grammar::NonTerminal::W, split.k + 1, parent.j);
      return out;

    case grammar::RuleId::WI_SPLIT_V:
      if (split.k > parent.i) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::WI, parent.i, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      }
      return out;

    case grammar::RuleId::WI_SPLIT_WMB:
      if (split.k > parent.i) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::WI, parent.i, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      }
      return out;

    case grammar::RuleId::WI_EXTEND_UNPAIRED:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WI, parent.i, parent.j - 1);
      return out;

    case grammar::RuleId::V_VM:
      out.child_count = 1;
      // Exact V->VM reuses the same (i,j) VM span; the closing pair energy is
      // already carried by the V-local weight.
      out.children[0].item = make_item(grammar::NonTerminal::VM, parent.i, parent.j);
      return out;

    case grammar::RuleId::VM_SPLIT_WM_WMv:
      if (split.k > parent.i + 1) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::WM, parent.i + 1, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::WMv, split.k, parent.j - 1);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::WMv, split.k, parent.j - 1);
      }
      return out;

    case grammar::RuleId::VM_SPLIT_WM_WMp:
      if (split.k > parent.i + 1) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::WM, parent.i + 1, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::WMp, split.k, parent.j - 1);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::WMp, split.k, parent.j - 1);
      }
      return out;

    case grammar::RuleId::VM_SPLIT_WMp_BASE:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WMp, split.k, parent.j - 1);
      return out;

    case grammar::RuleId::WM_START_V:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      return out;

    case grammar::RuleId::WM_SPLIT_V:
      if (split.k > parent.i) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::WM, parent.i, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      }
      return out;

    case grammar::RuleId::WM_START_WMB:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      return out;

    case grammar::RuleId::WM_SPLIT_WMB:
      if (split.k > parent.i) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::WM, parent.i, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      }
      return out;

    case grammar::RuleId::WM_EXTEND_UNPAIRED:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WM, parent.i, parent.j - 1);
      return out;

    case grammar::RuleId::WMv_STEM_V:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::V, parent.i, parent.j);
      return out;

    case grammar::RuleId::WMv_EXTEND_UNPAIRED:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WMv, parent.i, parent.j - 1);
      return out;

    case grammar::RuleId::WMp_STEM_WMB:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WMB, parent.i, parent.j);
      return out;

    case grammar::RuleId::WMp_EXTEND_UNPAIRED:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WMp, parent.i, parent.j - 1);
      return out;

    case grammar::RuleId::WIP_BASE_V:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::V, parent.i, parent.j);
      return out;

    case grammar::RuleId::WIP_BASE_WMB:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WMB, parent.i, parent.j);
      return out;

    case grammar::RuleId::WIP_SPLIT_V:
      if (split.k > parent.i) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::WIP, parent.i, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      }
      return out;

    case grammar::RuleId::WIP_SPLIT_WMB:
      if (split.k > parent.i) {
        out.child_count = 2;
        out.children[0].item = make_item(grammar::NonTerminal::WIP, parent.i, split.k - 1);
        out.children[1].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      } else {
        out.child_count = 1;
        out.children[0].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      }
      return out;

    case grammar::RuleId::WIP_BASEPAIR_V:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::V, split.k, parent.j);
      return out;

    case grammar::RuleId::WIP_BASEPAIR_WMB:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WMB, split.k, parent.j);
      return out;

    case grammar::RuleId::WIP_EXTEND_UNPAIRED:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WIP, parent.i, parent.j - 1);
      return out;

    case grammar::RuleId::VPL_SPLIT_VP:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::VP, split.k, parent.j);
      return out;

    case grammar::RuleId::VPR_SPLIT_VP_WIP:
      out.child_count = 2;
      out.children[0].item = make_item(grammar::NonTerminal::VP, parent.i, split.k);
      out.children[1].item = make_item(grammar::NonTerminal::WIP, split.k + 1, parent.j);
      return out;

    case grammar::RuleId::VPR_SPLIT_VP_BASEPAIR:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::VP, parent.i, split.k);
      return out;

    case grammar::RuleId::VP_WI_CASE1:
      out.child_count = 0;
      if (parent.i + 1 <= split.p - 1) {
        out.children[out.child_count++].item =
            make_item(grammar::NonTerminal::WI, parent.i + 1, split.p - 1);
      }
      if (split.q + 1 <= parent.j - 1) {
        out.children[out.child_count++].item =
            make_item(grammar::NonTerminal::WI, split.q + 1, parent.j - 1);
      }
      return out;

    case grammar::RuleId::VP_WI_CASE2:
      out.child_count = 0;
      if (parent.i + 1 <= split.p - 1) {
        out.children[out.child_count++].item =
            make_item(grammar::NonTerminal::WI, parent.i + 1, split.p - 1);
      }
      if (split.q + 1 <= parent.j - 1) {
        out.children[out.child_count++].item =
            make_item(grammar::NonTerminal::WI, split.q + 1, parent.j - 1);
      }
      return out;

    case grammar::RuleId::VP_WI_CASE3:
      out.child_count = 0;
      if (parent.i + 1 <= split.p - 1) {
        out.children[out.child_count++].item =
            make_item(grammar::NonTerminal::WI, parent.i + 1, split.p - 1);
      }
      if (split.q + 1 <= split.k - 1) {
        out.children[out.child_count++].item =
            make_item(grammar::NonTerminal::WI, split.q + 1, split.k - 1);
      }
      if (split.l + 1 <= parent.j - 1) {
        out.children[out.child_count++].item =
            make_item(grammar::NonTerminal::WI, split.l + 1, parent.j - 1);
      }
      return out;

    case grammar::RuleId::VP_STACK:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::VP, parent.i + 1, parent.j - 1);
      return out;

    case grammar::RuleId::VP_INTERNAL_LOOP:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::VP, split.k, split.l);
      return out;

    case grammar::RuleId::VP_WIP_VP_LEFT:
      out.child_count = 2;
      out.children[0].item = make_item(grammar::NonTerminal::WIP, parent.i + 1, split.k - 1);
      out.children[1].item = make_item(grammar::NonTerminal::VP, split.k, parent.j - 1);
      return out;

    case grammar::RuleId::VP_VP_WIP_RIGHT:
      out.child_count = 2;
      out.children[0].item = make_item(grammar::NonTerminal::VP, parent.i + 1, split.k);
      out.children[1].item = make_item(grammar::NonTerminal::WIP, split.k + 1, parent.j - 1);
      return out;

    case grammar::RuleId::VP_WIP_VPR:
      out.child_count = 2;
      out.children[0].item = make_item(grammar::NonTerminal::WIP, parent.i + 1, split.k - 1);
      out.children[1].item = make_item(grammar::NonTerminal::VPR, split.k, parent.j - 1);
      return out;

    case grammar::RuleId::VP_VPL_WIP:
      out.child_count = 2;
      out.children[0].item = make_item(grammar::NonTerminal::VPL, parent.i + 1, split.k);
      out.children[1].item = make_item(grammar::NonTerminal::WIP, split.k + 1, parent.j - 1);
      return out;

    case grammar::RuleId::BE_STACK: {
      const cand_pos_t next_ip = (parent.i == parent.ip && parent.j == parent.jp) ? (parent.i + 1) : parent.ip;
      const cand_pos_t next_jp = (parent.i == parent.ip && parent.j == parent.jp) ? (parent.j - 1) : parent.jp;
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::BE, parent.i + 1, parent.j - 1, next_ip, next_jp);
      return out;
    }

    case grammar::RuleId::BE_INTERNAL_LOOP:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::BE, split.k, split.l, parent.ip, parent.jp);
      return out;

    case grammar::RuleId::BE_WIP_WIP:
      out.child_count = 3;
      out.children[0].item = make_item(grammar::NonTerminal::WIP, parent.i + 1, split.k - 1);
      out.children[1].item = make_item(grammar::NonTerminal::BE, split.k, split.l, parent.ip, parent.jp);
      out.children[2].item = make_item(grammar::NonTerminal::WIP, split.l + 1, parent.j - 1);
      return out;

    case grammar::RuleId::BE_WIP_BASEPAIR:
      out.child_count = 2;
      out.children[0].item = make_item(grammar::NonTerminal::WIP, parent.i + 1, split.k - 1);
      out.children[1].item = make_item(grammar::NonTerminal::BE, split.k, split.l, parent.ip, parent.jp);
      return out;

    case grammar::RuleId::BE_BASEPAIR_WIP:
      out.child_count = 2;
      out.children[0].item = make_item(grammar::NonTerminal::BE, split.k, split.l, parent.ip, parent.jp);
      out.children[1].item = make_item(grammar::NonTerminal::WIP, split.l + 1, parent.j - 1);
      return out;

    case grammar::RuleId::WMBW_SPLIT_WMBP_WI:
      out.child_count = 2;
      out.children[0].item = make_item(grammar::NonTerminal::WMBP, parent.i, split.k);
      out.children[1].item = make_item(grammar::NonTerminal::WI, split.k + 1, parent.j);
      return out;

    case grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP:
      out.child_count = 3;
      out.children[0].item = make_item(grammar::NonTerminal::BE, split.p, split.q, split.k, split.l);
      out.children[1].item = make_item(grammar::NonTerminal::WMBP, parent.i, split.k - 1);
      out.children[2].item = make_item(grammar::NonTerminal::VP, split.k, parent.j);
      return out;

    case grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP:
      out.child_count = 3;
      out.children[0].item = make_item(grammar::NonTerminal::BE, split.p, split.q, split.k, split.l);
      out.children[1].item = make_item(grammar::NonTerminal::WMBW, parent.i, split.k - 1);
      out.children[2].item = make_item(grammar::NonTerminal::VP, split.k, parent.j);
      return out;

    case grammar::RuleId::WMBP_DIRECT_VP:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::VP, parent.i, parent.j);
      return out;

    case grammar::RuleId::WMBP_SPLIT_BE_WI_VP:
      out.child_count = 0;
      out.children[out.child_count++].item =
          make_item(grammar::NonTerminal::BE, parent.i, split.p, split.k, split.l);
      if (split.k + 1 <= split.p - 1) {
        out.children[out.child_count++].item =
            make_item(grammar::NonTerminal::WI, split.k + 1, split.p - 1);
      }
      out.children[out.child_count++].item = make_item(grammar::NonTerminal::VP, split.p, parent.j);
      return out;

    case grammar::RuleId::WMB_SPLIT_BE_WMBP_WI:
      out.child_count = 0;
      out.children[out.child_count++].item =
          make_item(grammar::NonTerminal::BE, split.p, parent.j, split.q, split.l);
      out.children[out.child_count++].item = make_item(grammar::NonTerminal::WMBP, parent.i, split.k);
      if (split.k + 1 <= split.l - 1) {
        out.children[out.child_count++].item =
            make_item(grammar::NonTerminal::WI, split.k + 1, split.l - 1);
      }
      return out;

    case grammar::RuleId::WMB_DIRECT_WMBP:
      out.child_count = 1;
      out.children[0].item = make_item(grammar::NonTerminal::WMBP, parent.i, parent.j);
      return out;

    default:
      return out;
  }
}

}  // namespace detail

inline constexpr Deduction expand_basic_rule(const ItemKey &parent,
                                             grammar::RuleId rule,
                                             const grammar::RuleSplit &split) {
  return detail::expand_basic_rule_impl(parent, rule, split);
}

template <grammar::RuleId Rule>
inline constexpr Deduction expand_basic_rule(
    const ItemKey &parent,
    const grammar::RuleSplit &split) {
  return detail::expand_basic_rule_impl(
      parent, std::integral_constant<grammar::RuleId, Rule>{}, split);
}

}  // namespace cparty::scfg2::deductions

#endif  // CPARTY_SCFG2_DEDUCTIONS_EXPAND_BASIC_HH

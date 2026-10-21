#ifndef CPARTY_SCFG2_DEDUCTIONS_GENERATE_EXACT_BASIC_HH
#define CPARTY_SCFG2_DEDUCTIONS_GENERATE_EXACT_BASIC_HH

#include "scfg2/deductions/exact_context.hh"
#include "scfg2/deductions/generate_basic.hh"
#include "scfg2/deductions/split_positions.hh"
#include "scfg2/deductions/vp_geometry.hh"
#include "scfg2/deductions/wmb_geometry.hh"

#include <ViennaRNA/params/constants.h>
#include <algorithm>
#include <utility>
#include <vector>

namespace cparty::scfg2::deductions {

// Optional per-rule views for numeric exact evaluation. A null view
// preserves the original split range; an empty selection emits no splits.
struct ExactNumericSplitViews {
  const std::vector<cand_pos_t> *wm_start_v = nullptr, *wm_start_wmb = nullptr;
  const std::vector<cand_pos_t> *wm_split_v = nullptr, *wm_split_wmb = nullptr;
  const std::vector<cand_pos_t> *vm_split_wmv = nullptr, *vm_split_wmp = nullptr;
  const std::vector<cand_pos_t> *vm_base = nullptr;
  const std::vector<cand_pos_t> *wi_split_v = nullptr, *wi_split_wmb = nullptr;
  const std::vector<cand_pos_t> *wip_split_v = nullptr, *wip_split_wmb = nullptr;
  const std::vector<cand_pos_t> *wip_basepair_v = nullptr, *wip_basepair_wmb = nullptr;
  const std::vector<cand_pos_t> *vpr_split_wip = nullptr, *vpr_basepair = nullptr;
  const std::vector<cand_pos_t> *wmbp_vp_starts = nullptr;

  const std::vector<cand_pos_t> *for_rule(
      grammar::RuleId rule, const std::vector<cand_pos_t> *fallback) const {
    using R = grammar::RuleId;
    const std::vector<cand_pos_t> *view = nullptr;
    switch (rule) {
      case R::WM_START_V: view = wm_start_v; break;
      case R::WM_START_WMB: view = wm_start_wmb; break;
      case R::WM_SPLIT_V: view = wm_split_v; break;
      case R::WM_SPLIT_WMB: view = wm_split_wmb; break;
      case R::VM_SPLIT_WM_WMv: view = vm_split_wmv; break;
      case R::VM_SPLIT_WM_WMp: view = vm_split_wmp; break;
      case R::VM_SPLIT_WMp_BASE: view = vm_base; break;
      case R::WI_SPLIT_V: view = wi_split_v; break;
      case R::WI_SPLIT_WMB: view = wi_split_wmb; break;
      case R::WIP_SPLIT_V: view = wip_split_v; break;
      case R::WIP_SPLIT_WMB: view = wip_split_wmb; break;
      case R::WIP_BASEPAIR_V: view = wip_basepair_v; break;
      case R::WIP_BASEPAIR_WMB: view = wip_basepair_wmb; break;
      case R::VPR_SPLIT_VP_WIP: view = vpr_split_wip; break;
      case R::VPR_SPLIT_VP_BASEPAIR: view = vpr_basepair; break;
      case R::WMBP_SPLIT_BE_WMBP_VP: case R::WMBP_SPLIT_BE_WMBW_VP:
      case R::WMBP_SPLIT_BE_WI_VP: view = wmbp_vp_starts; break;
      default: break;
    }
    return view ? view : fallback;
  }
};

// Coordinate restricted by the optional one-dimensional selection below.
// Fixed-coordinate and MAXLOOP-bounded rules are not restricted by this view.
inline std::optional<cand_pos_t> exact_linear_split_coordinate(const Deduction &d) {
  using R = grammar::RuleId;
  switch (d.rule) {
    case R::WMBP_SPLIT_BE_WI_VP: return d.split.p;
    case R::WI_SPLIT_V: case R::WI_SPLIT_WMB:
    case R::W_SPLIT_V: case R::W_SPLIT_WMB: case R::W_ROOT_WMB_DECOMPOSITION:
    case R::WM_START_V: case R::WM_START_WMB:
    case R::WM_SPLIT_V: case R::WM_SPLIT_WMB:
    case R::VM_SPLIT_WM_WMv: case R::VM_SPLIT_WM_WMp: case R::VM_SPLIT_WMp_BASE:
    case R::WIP_SPLIT_V: case R::WIP_SPLIT_WMB:
    case R::WIP_BASEPAIR_V: case R::WIP_BASEPAIR_WMB:
    case R::VPL_SPLIT_VP: case R::VPR_SPLIT_VP_WIP: case R::VPR_SPLIT_VP_BASEPAIR:
    case R::VP_WIP_VP_LEFT: case R::VP_WIP_VPR:
    case R::VP_VP_WIP_RIGHT: case R::VP_VPL_WIP:
    case R::BE_INTERNAL_LOOP: case R::BE_WIP_WIP:
    case R::BE_WIP_BASEPAIR: case R::BE_BASEPAIR_WIP:
    case R::WMBW_SPLIT_WMBP_WI:
    case R::WMBP_SPLIT_BE_WMBP_VP: case R::WMBP_SPLIT_BE_WMBW_VP:
    case R::WMB_SPLIT_BE_WMBP_WI: return d.split.k;
    default: return std::nullopt;
  }
}

namespace {

inline cand_pos_t principal_root_wmb_open(const PairabilityContext &ctx,
                                          cand_pos_t i,
                                          cand_pos_t k) {
  cand_pos_t open = ctx.pair_square(k);
  if (open > 0 && open < k) {
    return open;
  }
  open = ctx.pair_any(k);
  if (open == i && open < k) {
    return open;
  }
  return -1;
}

inline cand_pos_t leftmost_covered_root_wmb_open(const PairabilityContext &ctx,
                                                 cand_pos_t i,
                                                 cand_pos_t k,
                                                 cand_pos_t open) {
  if (open < 0) {
    return open;
  }
  if (ctx.parent_index(open) >= 0) {
    return open;
  }
  for (cand_pos_t pos = open - 1; pos >= i; --pos) {
    const cand_pos_t close = ctx.pair_square(pos);
    if (close > pos && close <= k && ctx.parent_index(pos) >= 0) {
      return pos;
    }
    if (pos == i) {
      break;
    }
  }
  return open;
}

inline cand_pos_t leftmost_crossing_root_wmb_open(const PairabilityContext &ctx,
                                                  cand_pos_t i,
                                                  cand_pos_t k,
                                                  cand_pos_t open) {
  if (open < 0) {
    return -1;
  }
  for (cand_pos_t pos = i; pos <= k; ++pos) {
    const cand_pos_t close = ctx.pair_any(pos);
    if (close <= pos) {
      continue;
    }
    const bool crosses_from_left = pos < open && open < close && close < k;
    const bool crosses_from_right = open < pos && pos < k && k < close;
    if (crosses_from_left || crosses_from_right) {
      return pos;
    }
  }
  return -1;
}

inline cand_pos_t root_wmb_open(const PairabilityContext &ctx,
                                cand_pos_t i,
                                cand_pos_t k) {
  const cand_pos_t open = principal_root_wmb_open(ctx, i, k);
  if (!(open > 0 && open < k)) {
    return -1;
  }
  // The final opening is minimized with candidates in [i,k]. At the lower
  // boundary no earlier opening exists, even for dual carriers.
  if (open == i) return open;
  cand_pos_t p = leftmost_covered_root_wmb_open(ctx, i, k, open);
  const cand_pos_t crossing_open = leftmost_crossing_root_wmb_open(ctx, i, k, open);
  if (crossing_open > 0 && (p < 0 || crossing_open < p)) {
    p = crossing_open;
  }
  return p;
}

template <class QueryContext>
inline bool exact_v_pairing_state_allows_impl(const QueryContext &ctx,
                                              cand_pos_t i,
                                              cand_pos_t j) {
  const cand_pos_t pair_i = ctx.pair_partner(i);
  const cand_pos_t pair_j = ctx.pair_partner(j);
  const bool unpaired = pair_i < -1 && pair_j < -1;
  const bool paired = pair_i == j && pair_j == i;
  return paired || unpaired;
}

inline bool exact_v_pairing_state_allows(const PairabilityContext &ctx,
                                         cand_pos_t i,
                                         cand_pos_t j) {
  return exact_v_pairing_state_allows_impl(ctx, i, j);
}

inline bool exact_v_hairpin_allows(const PairabilityContext &ctx,
                                   cand_pos_t i,
                                   cand_pos_t j) {
  return exact_v_pairing_state_allows(ctx, i, j) &&
         ctx.unpaired_prefix(j - 1) >= (j - i - 1);
}

inline bool exact_vp_open_endpoint_allows(const PairabilityContext &ctx,
                                          cand_pos_t i,
                                          cand_pos_t j) {
  return ctx.is_unpaired(i) && ctx.is_unpaired(j);
}

inline bool exact_vp_parent_pair_admissible(const PairabilityContext &ctx,
                                            cand_pos_t i,
                                            cand_pos_t j) {
  // Match the current constrained-inside VP parent contract: VP spans are not
  // formed directly on scaffold round pairs.
  return ctx.pair_round(i) != j && ctx.can_pair(i, j);
}

inline bool exact_vp_public_endpoint_allows(const PairabilityContext &ctx,
                                            cand_pos_t i,
                                            cand_pos_t j) {
  // Legacy compute_VP distinguishes an ordinary free base (-2) from the
  // pair-forbidden `x` sentinel (-1).  VP creates a new crossing pair, so both
  // endpoints must be ordinary free bases rather than merely non-fixed.
  return ctx.pair_round(i) < -1 &&
         ctx.pair_round(j) < -1 &&
         (j - i) >= (ctx.turn() + 1);
}

inline bool exact_round_open_family_active(const PairabilityContext &ctx,
                                           cand_pos_t i,
                                           cand_pos_t j) {
  return j != i && (j - i) >= (ctx.turn() + 1) && !ctx.weakly_closed(i, j);
}

inline bool exact_round_weakly_closed_family_active(const PairabilityContext &ctx,
                                                    cand_pos_t i,
                                                    cand_pos_t j) {
  return i == j || ctx.weakly_closed(i, j);
}

inline bool exact_v_family_active(const PairabilityContext &ctx,
                                  cand_pos_t i,
                                  cand_pos_t j) {
  // Match the current baseline exact fill gate for round-lane V cells:
  // the closing pair must be physically admissible and the span must live on
  // the weakly-closed round family before V is materialized.
  return exact_round_weakly_closed_family_active(ctx, i, j) &&
         exact_v_pairing_state_allows(ctx, i, j) &&
         ctx.can_pair(i, j);
}

inline bool exact_round_pseudoloop_region_active(const PairabilityContext &ctx,
                                                 cand_pos_t i,
                                                 cand_pos_t j) {
  if ((j - i - 1) <= ctx.turn()) {
    return false;
  }
  const cand_pos_t pair_i = ctx.pair_round(i);
  const cand_pos_t pair_j = ctx.pair_round(j);
  // Ordinary free endpoints (-2) are permitted. Every other endpoint must
  // have its fixed partner inside [i,j]; the forbidden -1 marker fails this
  // test on a nucleotide span. This is only an endpoint test, not closure.
  const auto contained = [i, j](cand_pos_t partner) {
    return partner < -1 || (i <= partner && partner <= j);
  };
  return contained(pair_i) && contained(pair_j);
}

// Completed pseudoloops cannot leave a scaffold pair crossing their boundary.
// This does not apply to WMBP/WMBW prefixes, which are allowed to stay open.
inline bool exact_completed_wmb_parent_active(const PairabilityContext &ctx,
                                              cand_pos_t i, cand_pos_t j) {
  return exact_round_pseudoloop_region_active(ctx, i, j) &&
         (ctx.wmb_closure_policy() == WmbClosurePolicy::LegacyCompatible ||
          ctx.weakly_closed(i, j));
}

struct ExactDirectChildGeometry {
  cand_pos_t count = 0;
  cand_pos_t child_i = -1;
  cand_pos_t child_j = -1;
  cand_pos_t left_gap = -1;
  cand_pos_t right_gap = -1;
};

inline ExactDirectChildGeometry exact_direct_child_geometry(const PairabilityContext &ctx,
                                                            cand_pos_t i,
                                                            cand_pos_t j) {
  ExactDirectChildGeometry geom;
  std::vector<std::pair<cand_pos_t, cand_pos_t>> local_pairs;
  local_pairs.reserve(static_cast<size_t>(std::max<cand_pos_t>(0, j - i + 1)));
  for (cand_pos_t pos = i; pos <= j; ++pos) {
    const cand_pos_t round_partner = ctx.pair_round(pos);
    if (round_partner > pos && round_partner <= j) {
      local_pairs.emplace_back(pos, round_partner);
    }
    const cand_pos_t square_partner = ctx.pair_square(pos);
    if (square_partner > pos && square_partner <= j) {
      local_pairs.emplace_back(pos, square_partner);
    }
  }
  std::sort(local_pairs.begin(), local_pairs.end());
  local_pairs.erase(std::unique(local_pairs.begin(), local_pairs.end()), local_pairs.end());

  for (const auto &[left, right] : local_pairs) {
    if (left <= i || right >= j) {
      continue;
    }
    bool enclosed_by_other = false;
    for (const auto &[outer_left, outer_right] : local_pairs) {
      if (outer_left == i && outer_right == j) {
        continue;
      }
      if (outer_left < left && right < outer_right &&
          outer_left >= i && outer_right <= j) {
        enclosed_by_other = true;
        break;
      }
    }
    if (enclosed_by_other) {
      continue;
    }
    ++geom.count;
    if (geom.count == 1) {
      geom.child_i = left;
      geom.child_j = right;
      geom.left_gap = left - i - 1;
      geom.right_gap = j - right - 1;
    }
  }
  if (geom.count != 1) {
    geom.child_i = -1;
    geom.child_j = -1;
    geom.left_gap = -1;
    geom.right_gap = -1;
  }
  return geom;
}

inline cand_pos_t exact_vp_min_border_left(const PairabilityContext &ctx,
                                           cand_pos_t i,
                                           cand_pos_t j) {
  const cand_pos_t b = ctx.border_b(i, j);
  const cand_pos_t bp_cap = ctx.border_Bp(i, j);
  if (b < 0) {
    return bp_cap;
  }
  if (bp_cap < 0) {
    return b;
  }
  return std::min(b, bp_cap);
}

inline bool exact_vp_recursive_split_allows(
    const ItemKey &parent, grammar::RuleId rule, cand_pos_t k,
    const PairabilityContext &ctx) {
  // Legacy compute_VP cases 6--9 sum the full border range. A single
  // maximal G pair inside the interval does not determine the split point.
  switch (rule) {
    case grammar::RuleId::VP_WIP_VP_LEFT:
    case grammar::RuleId::VP_WIP_VPR:
      return k > parent.i && k < exact_vp_min_border_left(ctx, parent.i, parent.j);
    case grammar::RuleId::VP_VP_WIP_RIGHT:
    case grammar::RuleId::VP_VPL_WIP:
      return k > std::max(ctx.border_B(parent.i, parent.j), ctx.border_bp(parent.i, parent.j)) &&
             k < parent.j;
    default:
      return false;
  }
}

}  // namespace

inline bool exact_wmbp_item_active(const ItemKey &item,
                                   const PairabilityContext &ctx) {
  return item.nonterminal == grammar::NonTerminal::WMBP &&
         exact_round_pseudoloop_region_active(ctx, item.i, item.j);
}

inline bool exact_wmv_wmp_stem_parent_allows(
    const ItemKey &parent,
    const PairabilityContext &ctx) {
  return parent.j - parent.i - 1 >= ctx.turn();
}

inline bool exact_wmv_wmp_extend_parent_allows(
    const ItemKey &parent,
    const PairabilityContext &ctx) {
  return ctx.is_unpaired(parent.j) &&
         parent.j - parent.i >= ctx.turn();
}

inline bool exact_wm_split_parent_allows(const ItemKey &parent) {
  return parent.j - parent.i + 1 >= 4;
}

inline bool exact_wm_start_split_allows(
    const ItemKey &parent,
    cand_pos_t k,
    const PairabilityContext &ctx) {
  return exact_wm_split_parent_allows(parent) &&
         ctx.can_pair_left_span(parent.i, k);
}

inline bool exact_wm_extend_parent_allows(
    const ItemKey &parent,
    const PairabilityContext &ctx) {
  return ctx.is_unpaired(parent.j) &&
         exact_wm_split_parent_allows(parent);
}

inline cand_pos_t exact_vm_binary_begin_k(const ItemKey &parent) {
  return parent.i + 2;
}

inline bool exact_wip_start_split_allows(
    const ItemKey &parent, cand_pos_t k, const PairabilityContext &ctx) {
  // V(k,j) includes the legal minimal hairpin j-k-1 == TURN.
  return k > parent.i && k <= parent.j - ctx.turn() - 1 &&
         ctx.can_pair_left_span(parent.i, k);
}

inline bool exact_vpl_split_allows(
    const ItemKey &parent, cand_pos_t k, const PairabilityContext &ctx) {
  const auto border = std::min(static_cast<cand_pos_tu>(ctx.border_b(parent.i, parent.j)),
                               static_cast<cand_pos_tu>(ctx.border_Bp(parent.i, parent.j)));
  return k > parent.i && k < parent.j && static_cast<cand_pos_tu>(k) < border &&
         ctx.can_pair_left_span(parent.i, k);
}

inline cand_pos_t exact_vm_base_begin_k(const ItemKey &parent) {
  return parent.i + 1;
}

inline cand_pos_t exact_vm_end_k(const ItemKey &parent,
                                 const PairabilityContext &ctx) {
  return parent.j - ctx.turn() - 1;
}

inline bool exact_vm_base_split_allows(
    const ItemKey &parent,
    cand_pos_t k,
    const PairabilityContext &ctx) {
  return ctx.can_pair_left_span(parent.i + 1, k);
}

// Shared VP child identity for the aggregated WMBP families and numeric
// split selection. Every selected start is strictly inside the parent span.
inline ItemKey exact_wmbp_vp_child(const ItemKey &parent, cand_pos_t start) {
  return make_item(grammar::NonTerminal::VP, start, parent.j);
}

struct ExactWmbpBeVpSplitScan {
  std::vector<std::pair<cand_pos_t, WmbpSharedBeVpGeometry>> shared;
  std::vector<std::pair<cand_pos_t, WmbpWiVpAnchors>> wi;
  std::size_t l_visits = 0;
  std::size_t geometry_visits = 0;
  std::size_t overlap_count = 0;
};

inline void observe_exact_wmbp_be_vp_split_position(
    ExactWmbpBeVpSplitScan *out,
    const ItemKey &parent,
    const PairabilityContext &ctx,
    cand_pos_t l,
    bool shared_active,
    bool wi_active) {
  if (!out) return;
  ++out->l_visits;
  if (l <= parent.i || l >= parent.j || l + ctx.turn() > parent.j) return;
  const cand_pos_t scaffold_parent = ctx.parent_index(l);
  if (scaffold_parent < parent.i || scaffold_parent >= parent.j) return;
  ++out->geometry_visits;
  const auto shared = shared_active
      ? wmbp_shared_be_vp_geometry(ctx, parent, l)
      : std::optional<WmbpSharedBeVpGeometry>{};
  const auto wi = wi_active
      ? wmbp_wi_vp_geometry(ctx, parent, l)
      : std::optional<WmbpWiVpAnchors>{};
  if (shared.has_value()) out->shared.emplace_back(l, *shared);
  if (wi.has_value()) out->wi.emplace_back(l, *wi);
  if (shared.has_value() && wi.has_value()) ++out->overlap_count;
}

inline ExactWmbpBeVpSplitScan scan_exact_wmbp_be_vp_split_geometries(
    const ItemKey &parent,
    const PairabilityContext &ctx,
    bool scan_shared,
    bool scan_wi) {
  ExactWmbpBeVpSplitScan out;
  const bool shared_active =
      scan_shared && wmbp_right_boundary_admissible(ctx, parent.j);
  const bool wi_active =
      scan_wi && wmbp_wi_vp_parent_geometry_allows(ctx, parent);
  if (!shared_active && !wi_active) return out;
  for (cand_pos_t l = parent.i + 1; l < parent.j; ++l) {
    observe_exact_wmbp_be_vp_split_position(
        &out, parent, ctx, l, shared_active, wi_active);
  }
  return out;
}

inline ExactWmbpBeVpSplitScan
scan_exact_wmbp_be_vp_split_geometries_for_vp_starts(
    const ItemKey &parent,
    const PairabilityContext &ctx,
    const std::vector<cand_pos_t> &vp_starts,
    bool scan_shared,
    bool scan_wi) {
  ExactWmbpBeVpSplitScan out;
  const bool shared_active =
      scan_shared && wmbp_right_boundary_admissible(ctx, parent.j);
  const bool wi_active =
      scan_wi && wmbp_wi_vp_parent_geometry_allows(ctx, parent);
  if (!shared_active && !wi_active) return out;
  for (const cand_pos_t l : vp_starts) {
    observe_exact_wmbp_be_vp_split_position(
        &out, parent, ctx, l, shared_active, wi_active);
  }
  return out;
}

struct ExactWmbpBeVpFamilyDeductions {
  std::vector<Deduction> wmbp;
  std::vector<Deduction> wmbw;
  std::vector<Deduction> wi;
  std::size_t l_visits = 0;
  std::size_t overlap_count = 0;
};

inline ExactWmbpBeVpFamilyDeductions
generate_exact_wmbp_be_vp_family_deductions_from_scan(
    const ItemKey &parent,
    const PairabilityContext &ctx,
    const ExactWmbpBeVpSplitScan &scan,
    bool emit_wmbp,
    bool emit_wmbw) {
  ExactWmbpBeVpFamilyDeductions out;
  out.l_visits = scan.l_visits;
  out.overlap_count = scan.overlap_count;
  for (const auto &[l, geometry] : scan.shared) {
    grammar::RuleSplit split{};
    split.k = l;
    split.l = geometry.be_ip;
    split.p = geometry.be_i;
    split.q = geometry.be_j;
    const auto append_mode = [&](grammar::RuleId rule,
                                 grammar::NonTerminal left,
                                 std::vector<Deduction> *target) {
      Deduction ded;
      ded.parent = parent;
      ded.rule = rule;
      ded.split = split;
      ded.child_count = 3;
      ded.children[0].item = make_item(grammar::NonTerminal::BE,
                                       split.p,
                                       split.q,
                                       split.l,
                                       ctx.pair_partner(split.l));
      ded.children[1].item = make_item(left, parent.i, split.k - 1);
      ded.children[2].item =
          exact_wmbp_vp_child(parent, split.k);
      if (is_valid(ded)) target->push_back(ded);
    };
    if (emit_wmbp) {
      append_mode(grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP,
                  grammar::NonTerminal::WMBP,
                  &out.wmbp);
    }
    if (emit_wmbw) {
      append_mode(grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP,
                  grammar::NonTerminal::WMBW,
                  &out.wmbw);
    }
  }
  for (const auto &[l, geometry] : scan.wi) {
    grammar::RuleSplit split{};
    split.k = geometry.be_ip;
    split.l = geometry.be_jp;
    split.p = l;
    Deduction ded;
    ded.parent = parent;
    ded.rule = grammar::RuleId::WMBP_SPLIT_BE_WI_VP;
    ded.split = split;
    ded.child_count = 0;
    ded.children[ded.child_count++].item = make_item(
        grammar::NonTerminal::BE,
        parent.i,
        ctx.pair_partner(parent.i),
        split.k,
        split.l);
    if (split.k + 1 <= split.p - 1) {
      ded.children[ded.child_count++].item =
          make_item(grammar::NonTerminal::WI, split.k + 1, split.p - 1);
    }
    ded.children[ded.child_count++].item =
        exact_wmbp_vp_child(parent, split.p);
    if (is_valid(ded)) out.wi.push_back(ded);
  }
  return out;
}

inline ExactWmbpBeVpFamilyDeductions
generate_exact_wmbp_be_vp_family_deductions(
    const ItemKey &parent,
    const PairabilityContext &ctx,
    bool emit_wmbp,
    bool emit_wmbw,
    bool emit_wi) {
  const auto scan = scan_exact_wmbp_be_vp_split_geometries(
      parent, ctx, emit_wmbp || emit_wmbw, emit_wi);
  return generate_exact_wmbp_be_vp_family_deductions_from_scan(
      parent, ctx, scan, emit_wmbp, emit_wmbw);
}

inline ExactWmbpBeVpFamilyDeductions
generate_exact_wmbp_be_vp_family_deductions_for_vp_starts(
    const ItemKey &parent,
    const PairabilityContext &ctx,
    const std::vector<cand_pos_t> &vp_starts,
    bool emit_wmbp,
    bool emit_wmbw,
    bool emit_wi) {
  const auto scan = scan_exact_wmbp_be_vp_split_geometries_for_vp_starts(
      parent, ctx, vp_starts, emit_wmbp || emit_wmbw, emit_wi);
  return generate_exact_wmbp_be_vp_family_deductions_from_scan(
      parent, ctx, scan, emit_wmbp, emit_wmbw);
}

inline std::vector<Deduction>
generate_exact_wmbp_shared_be_vp_deductions(
    const ItemKey &parent,
    const PairabilityContext &ctx,
    bool emit_wmbp,
    bool emit_wmbw) {
  auto family = generate_exact_wmbp_be_vp_family_deductions(
      parent, ctx, emit_wmbp, emit_wmbw, false);
  family.wmbp.insert(family.wmbp.end(),
                     std::make_move_iterator(family.wmbw.begin()),
                     std::make_move_iterator(family.wmbw.end()));
  return family.wmbp;
}

inline void append_exact_wmbw_split_for_wmbp_end(
    std::vector<Deduction> &out,
    const ItemKey &parent,
    const PairabilityContext &ctx,
    cand_pos_t wmbp_end) {
  if (wmbp_end <= parent.i || wmbp_end >= parent.j ||
      !wmbw_parent_geometry_allows(ctx, parent, wmbp_end)) {
    return;
  }
  grammar::RuleSplit split{};
  split.k = wmbp_end;
  const auto deduction = expand_basic_rule(
      parent, grammar::RuleId::WMBW_SPLIT_WMBP_WI, split);
  if (is_valid(deduction)) out.push_back(deduction);
}

inline void append_exact_wmb_split_for_wmbp_end(
    std::vector<Deduction> &out,
    const ItemKey &parent,
    const PairabilityContext &ctx,
    cand_pos_t wmbp_end) {
  if (!exact_completed_wmb_parent_active(ctx, parent.i, parent.j)) return;
  if (!wmb_split_be_wmbp_wi_parent_geometry_allows(ctx, parent)) return;
  const cand_pos_t bp_j = ctx.pair_partner(parent.j);
  if (wmbp_end <= bp_j || wmbp_end >= parent.j) return;
  const auto geometry =
      wmb_split_be_wmbp_wi_geometry(ctx, parent, wmbp_end);
  if (!geometry.has_value()) return;
  grammar::RuleSplit split{};
  split.k = wmbp_end;
  split.l = geometry->be_jp;
  split.p = bp_j;
  split.q = geometry->be_ip;
  const auto deduction = expand_basic_rule(
      parent, grammar::RuleId::WMB_SPLIT_BE_WMBP_WI, split);
  if (is_valid(deduction)) out.push_back(deduction);
}

// Membership tests for single splits, used by beam event builders that
// probe one (parent, child) pair at a time.  Each returns true iff
// generate_exact_basic_deductions(parent, rule, ctx) emits split (k, l)
// before the is_valid filter; scfg2_split_membership_test checks this
// against the generators.

template <class QueryContext>
inline bool exact_v_internal_split_allows_impl(const ItemKey &parent,
                                               cand_pos_t k,
                                               cand_pos_t l,
                                               const QueryContext &ctx) {
  if (parent.j - parent.i - 1 < ctx.turn() ||
      !exact_v_pairing_state_allows_impl(ctx, parent.i, parent.j)) {
    return false;
  }
  const cand_pos_t max_k = std::min(parent.j - ctx.turn() - 2, parent.i + MAXLOOP + 1);
  if (k < parent.i + 1 || k > max_k) return false;
  if (ctx.unpaired_prefix(k - 1) < (k - parent.i - 1)) return false;
  const cand_pos_t min_l =
      std::max(k + ctx.turn() + 1 + MAXLOOP + 2, k + parent.j - parent.i) - MAXLOOP - 2;
  if (l < min_l || l > parent.j - 1) return false;
  if (ctx.unpaired_prefix(parent.j - 1) < (parent.j - l - 1)) return false;
  return ctx.can_pair(k, l);
}

inline bool exact_v_internal_split_allows(const ItemKey &parent,
                                          cand_pos_t k,
                                          cand_pos_t l,
                                          const PairabilityContext &ctx) {
  // Keep the generic context's original short-circuit and query trace.
  return exact_v_internal_split_allows_impl(parent, k, l, ctx);
}

// A query view for one certified, in-bounds child of a known pure context.
// The V agenda's explicit standard-context opt-in establishes that contract;
// ordinary/custom contexts continue through the public wrapper above.
template <class Context>
class ExactVInternalChildQueryCache {
 public:
  ExactVInternalChildQueryCache(const Context &ctx,
                               cand_pos_t k, cand_pos_t l,
                               cand_pos_t left_up)
      : ctx_(ctx), k_(k), l_(l), turn_(ctx.turn()), left_up_(left_up),
        pairable_(ctx.can_pair(k, l)) {}

  cand_pos_t turn() const { return turn_; }
  cand_pos_t pair_partner(cand_pos_t pos) const { return ctx_.pair_partner(pos); }
  cand_pos_t unpaired_prefix(cand_pos_t pos) const {
    return pos == k_ - 1 ? left_up_ : ctx_.unpaired_prefix(pos);
  }
  bool can_pair(cand_pos_t k, cand_pos_t l) const {
    return k == k_ && l == l_ ? pairable_ : ctx_.can_pair(k, l);
  }

 private:
  const Context &ctx_;
  cand_pos_t k_, l_, turn_, left_up_;
  bool pairable_;
};

inline bool exact_vp_internal_loop_split_allows(const ItemKey &parent,
                                                cand_pos_t k,
                                                cand_pos_t l,
                                                const PairabilityContext &ctx) {
  if (parent.j - parent.i < ctx.turn() ||
      !exact_vp_open_endpoint_allows(ctx, parent.i, parent.j)) {
    return false;
  }
  const cand_pos_t min_borders =
      std::min(std::max(ctx.border_Bp(parent.i, parent.j), cand_pos_t(0)),
               std::max(ctx.border_b(parent.i, parent.j), cand_pos_t(0)));
  const cand_pos_t edge_i = std::min(parent.i + MAXLOOP + 1, parent.j - ctx.turn() - 1);
  const cand_pos_t max_k = std::min(min_borders > 0 ? min_borders - 1 : edge_i, edge_i);
  if (k < parent.i + 1 || k > max_k) return false;
  if (!ctx.is_unpaired(k) || !ctx.is_empty_region(parent.i, k)) return false;
  cand_pos_t max_borders =
      std::max(ctx.border_bp(parent.i, parent.j), ctx.border_B(parent.i, parent.j));
  max_borders = std::max(max_borders, k + parent.j - parent.i - MAXLOOP - 2);
  if (l <= max_borders || l > parent.j - 1) return false;
  if (k == parent.i + 1 && l == parent.j - 1) return false;
  return ctx.is_unpaired(l) && ctx.is_empty_region(l, parent.j) && ctx.can_pair(k, l);
}

inline bool exact_be_split_allows(const ItemKey &parent,
                                  grammar::RuleId rule,
                                  cand_pos_t k,
                                  cand_pos_t l,
                                  const PairabilityContext &ctx) {
  if (!(parent.i >= 1 && parent.i < parent.j && parent.i <= parent.ip &&
        parent.ip < parent.jp && parent.jp <= parent.j &&
        ctx.pair_partner(parent.i) == parent.j &&
        ctx.pair_partner(parent.ip) == parent.jp)) {
    return false;
  }
  if (k < parent.i + 1 || k > parent.ip) return false;
  if (ctx.pair_partner(k) != l || !(l >= parent.jp && l < parent.j)) return false;
  switch (rule) {
    case grammar::RuleId::BE_INTERNAL_LOOP:
      return ctx.is_empty_region(parent.i, k) && ctx.is_empty_region(l, parent.j);
    case grammar::RuleId::BE_WIP_WIP:
      return ctx.weakly_closed(parent.i + 1, k - 1) &&
             ctx.weakly_closed(l + 1, parent.j - 1);
    case grammar::RuleId::BE_WIP_BASEPAIR:
      return ctx.weakly_closed(parent.i + 1, k - 1) && ctx.is_empty_region(l, parent.j);
    case grammar::RuleId::BE_BASEPAIR_WIP:
      return ctx.is_empty_region(parent.i, k) && ctx.weakly_closed(l + 1, parent.j - 1);
    default:
      return false;
  }
}

// Append preserves the caller's existing rule stream. Complex families may
// still materialize a temporary buffer; ordinary split rules write directly.
inline void append_exact_deduction_buffer(std::vector<Deduction> &out,
                                           std::vector<Deduction> generated) {
  if (out.empty() && out.capacity() < generated.size()) out = std::move(generated);
  else out.insert(out.end(), std::make_move_iterator(generated.begin()),
                  std::make_move_iterator(generated.end()));
}

template <grammar::RuleId Rule>
inline void append_exact_basic_deduction(std::vector<Deduction> &out,
                                          const ItemKey &parent,
                                          const grammar::RuleSplit &split) {
  const auto deduction = expand_basic_rule<Rule>(parent, split);
  if (is_valid(deduction)) out.push_back(deduction);
}

inline void append_exact_basic_deductions(
    std::vector<Deduction> &out,
    const ItemKey &parent,
    grammar::RuleId rule,
    const PairabilityContext &ctx,
    const std::vector<cand_pos_t> *split_positions = nullptr) {
  const auto push_if_valid = [&](const grammar::RuleSplit &split) {
    const auto ded = expand_basic_rule(parent, rule, split);
    if (is_valid(ded)) {
      out.push_back(ded);
    }
  };
  if (split_positions &&
      (rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP ||
       rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP ||
       rule == grammar::RuleId::WMBP_SPLIT_BE_WI_VP)) {
    auto family = generate_exact_wmbp_be_vp_family_deductions_for_vp_starts(
        parent, ctx, *split_positions,
        rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP,
        rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP,
        rule == grammar::RuleId::WMBP_SPLIT_BE_WI_VP);
    if (rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP)
      append_exact_deduction_buffer(out, std::move(family.wmbp));
    else if (rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP)
      append_exact_deduction_buffer(out, std::move(family.wmbw));
    else
      append_exact_deduction_buffer(out, std::move(family.wi));
    return;
  }
  switch (rule) {
    case grammar::RuleId::W_EXTEND_UNPAIRED:
    case grammar::RuleId::WI_EXTEND_UNPAIRED:
    case grammar::RuleId::WIP_EXTEND_UNPAIRED:
      if (ctx.is_unpaired(parent.j)) {
        append_exact_deduction_buffer(out, generate_basic_deductions_single(parent, rule));
        return;
      }
      return;
    case grammar::RuleId::WM_EXTEND_UNPAIRED:
      if (exact_wm_extend_parent_allows(parent, ctx)) {
        append_exact_deduction_buffer(out, generate_basic_deductions_single(parent, rule));
        return;
      }
      return;
    case grammar::RuleId::WMv_EXTEND_UNPAIRED:
      if (exact_wmv_wmp_extend_parent_allows(parent, ctx)) {
        append_exact_deduction_buffer(out, generate_basic_deductions_single(parent, rule));
        return;
      }
      return;
    case grammar::RuleId::WI_SPLIT_V:
      for (const cand_pos_t k : ExactSplitPositions(parent.i, parent.j - ctx.turn() - 1, split_positions)) {
        if (!ctx.can_pair(k, parent.j)) {
          continue;
        }
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::WI_SPLIT_V>(out, parent, split);
      }
      return;
    case grammar::RuleId::W_SPLIT_V:
      if (!ctx.weakly_closed(1, parent.j)) {
        return;
      }
      for (const cand_pos_t k : ExactSplitPositions(parent.i, parent.j - 1, split_positions)) {
        if (!ctx.weakly_closed(1, k - 1)) {
          continue;
        }
        if (!ctx.can_pair(k, parent.j)) {
          continue;
        }
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::W_SPLIT_V>(out, parent, split);
      }
      return;
    case grammar::RuleId::WM_START_V:
      if (!exact_wm_split_parent_allows(parent)) {
        return;
      }
      for (const cand_pos_t k : ExactSplitPositions(parent.i, parent.j, split_positions)) {
        if (!exact_wm_start_split_allows(parent, k, ctx)) {
          continue;
        }
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::WM_START_V>(out, parent, split);
      }
      return;
    case grammar::RuleId::WM_SPLIT_V:
      if (!exact_wm_split_parent_allows(parent)) {
        return;
      }
      for (const cand_pos_t k : ExactSplitPositions(parent.i + 1, parent.j, split_positions)) {
        grammar::RuleSplit split{};
        split.k = k;
        const auto ded = expand_basic_rule<grammar::RuleId::WM_SPLIT_V>(parent, split);
        if (is_valid(ded)) out.push_back(ded);
      }
      return;
    case grammar::RuleId::WM_START_WMB:
      if (!exact_wm_split_parent_allows(parent)) {
        return;
      }
      for (const cand_pos_t k : ExactSplitPositions(parent.i, parent.j, split_positions)) {
        if (!exact_wm_start_split_allows(parent, k, ctx)) {
          continue;
        }
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::WM_START_WMB>(out, parent, split);
      }
      return;
    case grammar::RuleId::WM_SPLIT_WMB:
      if (!exact_wm_split_parent_allows(parent)) {
        return;
      }
      for (const cand_pos_t k : ExactSplitPositions(parent.i + 1, parent.j, split_positions)) {
        grammar::RuleSplit split{};
        split.k = k;
        const auto ded = expand_basic_rule<grammar::RuleId::WM_SPLIT_WMB>(parent, split);
        if (is_valid(ded)) out.push_back(ded);
      }
      return;
    case grammar::RuleId::W_SPLIT_WMB:
      if (!ctx.weakly_closed(1, parent.j)) {
        return;
      }
      for (const cand_pos_t k : ExactSplitPositions(parent.i, parent.j - 1, split_positions)) {
        if (!ctx.weakly_closed(1, k - 1)) {
          continue;
        }
        if (!(k == parent.i || ctx.weakly_closed(k, parent.j))) {
          continue;
        }
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::W_SPLIT_WMB>(out, parent, split);
      }
      return;
    case grammar::RuleId::W_ROOT_WMB_DECOMPOSITION:
      if (parent.i != 1) {
        return;
      }
      for (const cand_pos_t k : ExactSplitPositions(parent.i, parent.j, split_positions)) {
        // `square` is only the crossing carrier here: the root split is still a
        // single-kernel deduction, but this family needs carrier geometry from
        // the scaffold to recover the legacy root/WMB decomposition shape.
        const cand_pos_t p = root_wmb_open(ctx, parent.i, k);
        if (!(p > 0 && p < k)) {
          continue;
        }
        if (!ctx.weakly_closed(parent.i, p - 1)) {
          continue;
        }
        if (!(k == parent.j || ctx.weakly_closed(k + 1, parent.j))) {
          continue;
        }
        grammar::RuleSplit split{};
        split.k = k;
        split.p = p;
        append_exact_basic_deduction<grammar::RuleId::W_ROOT_WMB_DECOMPOSITION>(out, parent, split);
      }
      return;
    case grammar::RuleId::VM_SPLIT_WM_WMv:
      // Current exact VM case1 geometry always carries an explicit left WM
      // child for this rule family. The k=i+1 edge therefore dies through the
      // empty left span there, so scfg2 should not reopen it as a no-left-WM
      // deduction.
      for (const cand_pos_t k : ExactSplitPositions(exact_vm_binary_begin_k(parent),
                                                   exact_vm_end_k(parent, ctx), split_positions)) {
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::VM_SPLIT_WM_WMv>(out, parent, split);
      }
      return;
    case grammar::RuleId::VM_SPLIT_WM_WMp:
      for (const cand_pos_t k : ExactSplitPositions(exact_vm_binary_begin_k(parent),
                                                   exact_vm_end_k(parent, ctx), split_positions)) {
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::VM_SPLIT_WM_WMp>(out, parent, split);
      }
      return;
    case grammar::RuleId::VM_SPLIT_WMp_BASE:
      for (const cand_pos_t k : ExactSplitPositions(exact_vm_base_begin_k(parent),
                                                   exact_vm_end_k(parent, ctx), split_positions)) {
        if (!exact_vm_base_split_allows(parent, k, ctx)) {
          continue;
        }
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::VM_SPLIT_WMp_BASE>(out, parent, split);
      }
      return;
    case grammar::RuleId::WIP_SPLIT_V:
      for (const cand_pos_t k : ExactSplitPositions(parent.i + 1, parent.j - ctx.turn() - 1, split_positions)) {
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::WIP_SPLIT_V>(out, parent, split);
      }
      return;
    case grammar::RuleId::WIP_SPLIT_WMB:
      for (const cand_pos_t k : ExactSplitPositions(parent.i + 1, parent.j - ctx.turn() - 1, split_positions)) {
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::WIP_SPLIT_WMB>(out, parent, split);
      }
      return;
    case grammar::RuleId::WIP_BASEPAIR_V:
      for (const cand_pos_t k : ExactSplitPositions(parent.i + 1, parent.j - ctx.turn() - 1, split_positions)) {
        if (!exact_wip_start_split_allows(parent, k, ctx)) {
          continue;
        }
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::WIP_BASEPAIR_V>(out, parent, split);
      }
      return;
    case grammar::RuleId::WIP_BASEPAIR_WMB:
      for (const cand_pos_t k : ExactSplitPositions(parent.i + 1, parent.j - ctx.turn() - 1, split_positions)) {
        if (!exact_wip_start_split_allows(parent, k, ctx)) {
          continue;
        }
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::WIP_BASEPAIR_WMB>(out, parent, split);
      }
      return;
    case grammar::RuleId::VPL_SPLIT_VP:
      // Match legacy/current VPL banding: valid left-open VP anchors stop
      // before the first relevant b/Bp boundary. Negative boundaries stay
      // effectively unbounded via the unsigned sentinel convention.
      {
        const cand_pos_tu min_boundary =
            std::min(static_cast<cand_pos_tu>(ctx.border_b(parent.i, parent.j)),
                     static_cast<cand_pos_tu>(ctx.border_Bp(parent.i, parent.j)));
        const cand_pos_t end_k = min_boundary <= static_cast<cand_pos_tu>(parent.j)
            ? static_cast<cand_pos_t>(min_boundary) - 1 : parent.j - 1;
        for (const cand_pos_t k : ExactSplitPositions(parent.i + 1, end_k, split_positions)) {
          if (!exact_vpl_split_allows(parent, k, ctx)) {
            continue;
          }
          grammar::RuleSplit split{};
          split.k = k;
          append_exact_basic_deduction<grammar::RuleId::VPL_SPLIT_VP>(out, parent, split);
        }
      }
      return;
    case grammar::RuleId::VPR_SPLIT_VP_WIP:
      // Match current/legacy BandMaxBpRange: right-open public VP tails only
      // start strictly outside the blocking round-pair boundaries.
      {
        const cand_pos_t max_bp =
            std::max(ctx.border_B(parent.i, parent.j), ctx.border_bp(parent.i, parent.j));
        const cand_pos_t begin_k = std::max(static_cast<cand_pos_t>(parent.i + 1),
                                            static_cast<cand_pos_t>(max_bp + 1));
        for (const cand_pos_t k : ExactSplitPositions(begin_k, parent.j - 1, split_positions)) {
          grammar::RuleSplit split{};
          split.k = k;
          const auto ded = expand_basic_rule<grammar::RuleId::VPR_SPLIT_VP_WIP>(parent, split);
          if (is_valid(ded)) {
            out.push_back(ded);
          }
        }
      }
      return;
    case grammar::RuleId::VPR_SPLIT_VP_BASEPAIR:
      {
        const cand_pos_t max_bp =
            std::max(ctx.border_B(parent.i, parent.j), ctx.border_bp(parent.i, parent.j));
        const cand_pos_t begin_k = std::max(static_cast<cand_pos_t>(parent.i + 1),
                                            static_cast<cand_pos_t>(max_bp + 1));
        for (const cand_pos_t k : ExactSplitPositions(begin_k, parent.j - 1, split_positions)) {
          if (!ctx.can_pair_right_span(k, parent.j)) {
            continue;
          }
          grammar::RuleSplit split{};
          split.k = k;
          append_exact_basic_deduction<grammar::RuleId::VPR_SPLIT_VP_BASEPAIR>(out, parent, split);
        }
      }
      return;
    case grammar::RuleId::WIP_BASE_WMB:
      push_if_valid(grammar::RuleSplit{});
      return;
    case grammar::RuleId::VP_WI_CASE1:
    case grammar::RuleId::VP_WI_CASE2:
    case grammar::RuleId::VP_WI_CASE3: {
      if (!exact_vp_open_endpoint_allows(ctx, parent.i, parent.j)) {
        return;
      }
      const auto geometry = observe_vp_wi_geometry(ctx, parent.i, parent.j);
      if (const auto split = geometry.split_for(rule)) push_if_valid(*split);
      return;
    }
    case grammar::RuleId::VP_STACK:
      if (parent.j - parent.i >= ctx.turn() &&
          parent.i + 1 < parent.j &&
          ctx.is_unpaired(parent.i + 1) &&
          ctx.is_unpaired(parent.j - 1) &&
          ctx.can_pair(parent.i + 1, parent.j - 1)) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::VP_INTERNAL_LOOP:
      if (parent.j - parent.i < ctx.turn() ||
          !exact_vp_open_endpoint_allows(ctx, parent.i, parent.j)) {
        return;
      }
      {
        const auto domain = vp_internal_loop_domain(ctx, parent.i, parent.j, MAXLOOP);
        for (cand_pos_t k = domain.first_left; k <= domain.last_left; ++k) {
          if (!ctx.is_unpaired(k) || !ctx.is_empty_region(parent.i, k)) {
            continue;
          }
          const cand_pos_t first_right = domain.first_right(k);
          for (cand_pos_t l = domain.last_right; l >= first_right; --l) {
            if (k == parent.i + 1 && l == parent.j - 1) {
              continue;
            }
            if (!ctx.is_unpaired(l) || !ctx.is_empty_region(l, parent.j) || !ctx.can_pair(k, l)) {
              continue;
            }
            grammar::RuleSplit split{};
            split.k = k;
            split.l = l;
            append_exact_basic_deduction<grammar::RuleId::VP_INTERNAL_LOOP>(out, parent, split);
          }
        }
      }
      return;
    case grammar::RuleId::VP_WIP_VP_LEFT:
    case grammar::RuleId::VP_WIP_VPR: {
      const cand_pos_t min_bp_j = exact_vp_min_border_left(ctx, parent.i, parent.j);
      for (const cand_pos_t k : ExactSplitPositions(parent.i + 1, min_bp_j - 1, split_positions)) {
        grammar::RuleSplit split{};
        split.k = k;
        push_if_valid(split);
      }
      return;
    }
    case grammar::RuleId::VP_VP_WIP_RIGHT:
    case grammar::RuleId::VP_VPL_WIP: {
      const cand_pos_t max_i_bp = std::max(ctx.border_B(parent.i, parent.j), ctx.border_bp(parent.i, parent.j));
      for (const cand_pos_t k : ExactSplitPositions(max_i_bp + 1, parent.j - 1, split_positions)) {
        grammar::RuleSplit split{};
        split.k = k;
        push_if_valid(split);
      }
      return;
    }
    case grammar::RuleId::BE_BASE_SAMEPAIR:
      if (parent.i == parent.ip && parent.j == parent.jp && parent.i < parent.j) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::BE_STACK:
      if (parent.i >= 1 &&
          parent.i < parent.j &&
          parent.i <= parent.ip &&
          parent.ip < parent.jp &&
          parent.jp <= parent.j &&
          ctx.pair_partner(parent.i) == parent.j &&
          ctx.pair_partner(parent.ip) == parent.jp &&
          ctx.pair_partner(parent.i + 1) == parent.j - 1) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::BE_INTERNAL_LOOP:
    case grammar::RuleId::BE_WIP_WIP:
    case grammar::RuleId::BE_WIP_BASEPAIR:
    case grammar::RuleId::BE_BASEPAIR_WIP:
      if (!(parent.i >= 1 &&
            parent.i < parent.j &&
            parent.i <= parent.ip &&
            parent.ip < parent.jp &&
            parent.jp <= parent.j &&
            ctx.pair_partner(parent.i) == parent.j &&
            ctx.pair_partner(parent.ip) == parent.jp)) {
        return;
      }
      for (const cand_pos_t k : ExactSplitPositions(parent.i + 1, parent.ip, split_positions)) {
        const cand_pos_t l = ctx.pair_partner(k);
        if (!(l >= parent.jp && l < parent.j)) {
          continue;
        }
        const bool empty_left = ctx.is_empty_region(parent.i, k);
        const bool empty_right = ctx.is_empty_region(l, parent.j);
        const bool wc_left = ctx.weakly_closed(parent.i + 1, k - 1);
        const bool wc_right = ctx.weakly_closed(l + 1, parent.j - 1);
        grammar::RuleSplit split{};
        split.k = k;
        split.l = l;
        switch (rule) {
          case grammar::RuleId::BE_INTERNAL_LOOP:
            if (empty_left && empty_right) {
              push_if_valid(split);
            }
            break;
          case grammar::RuleId::BE_WIP_WIP:
            if (wc_left && wc_right) {
              push_if_valid(split);
            }
            break;
          case grammar::RuleId::BE_WIP_BASEPAIR:
            if (wc_left && empty_right) {
              push_if_valid(split);
            }
            break;
          case grammar::RuleId::BE_BASEPAIR_WIP:
            if (empty_left && wc_right) {
              push_if_valid(split);
            }
            break;
          default:
            break;
        }
      }
      return;
    case grammar::RuleId::WMBW_SPLIT_WMBP_WI:
      for (const cand_pos_t l : ExactSplitPositions(parent.i + 1, parent.j - 1, split_positions)) {
        append_exact_wmbw_split_for_wmbp_end(out, parent, ctx, l);
      }
      return;
    case grammar::RuleId::WMBP_DIRECT_VP:
      push_if_valid(grammar::RuleSplit{});
      return;
    case grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP:
      append_exact_deduction_buffer(out, generate_exact_wmbp_shared_be_vp_deductions(
          parent, ctx, true, false));
      return;
    case grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP:
      append_exact_deduction_buffer(out, generate_exact_wmbp_shared_be_vp_deductions(
          parent, ctx, false, true));
      return;
    case grammar::RuleId::WMBP_SPLIT_BE_WI_VP:
      append_exact_deduction_buffer(out, generate_exact_wmbp_be_vp_family_deductions(
                 parent, ctx, false, false, true).wi);
      return;
    case grammar::RuleId::WMB_DIRECT_WMBP:
      if (!exact_completed_wmb_parent_active(ctx, parent.i, parent.j)) return;
      push_if_valid(grammar::RuleSplit{});
      return;
    case grammar::RuleId::WMB_EMPTY:
      if (!exact_completed_wmb_parent_active(ctx, parent.i, parent.j)) return;
      if (parent.i < parent.j) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::WMB_SPLIT_BE_WMBP_WI:
      if (!exact_completed_wmb_parent_active(ctx, parent.i, parent.j)) return;
      if (wmb_split_be_wmbp_wi_parent_geometry_allows(ctx, parent)) {
        const cand_pos_t bp_j = ctx.pair_partner(parent.j);
        for (const cand_pos_t l : ExactSplitPositions(bp_j + 1, parent.j - 1, split_positions)) {
          append_exact_wmb_split_for_wmbp_end(out, parent, ctx, l);
        }
      }
      return;
    case grammar::RuleId::V_INTERNAL:
      if (parent.j - parent.i - 1 < ctx.turn() ||
          !exact_v_pairing_state_allows(ctx, parent.i, parent.j)) {
        return;
      }
      {
        const cand_pos_t max_k = std::min(parent.j - ctx.turn() - 2, parent.i + MAXLOOP + 1);
        for (cand_pos_t k = parent.i + 1; k <= max_k; ++k) {
          if (ctx.unpaired_prefix(k - 1) < (k - parent.i - 1)) {
            continue;
          }
          const cand_pos_t min_l =
              std::max(k + ctx.turn() + 1 + MAXLOOP + 2, k + parent.j - parent.i) - MAXLOOP - 2;
          for (cand_pos_t l = parent.j - 1; l >= min_l; --l) {
            if (ctx.unpaired_prefix(parent.j - 1) < (parent.j - l - 1)) {
              continue;
            }
            if (!ctx.can_pair(k, l)) {
              continue;
            }
            grammar::RuleSplit split{};
            split.k = k;
            split.l = l;
            append_exact_basic_deduction<grammar::RuleId::V_INTERNAL>(out, parent, split);
          }
        }
      }
      return;
    case grammar::RuleId::V_HAIRPIN:
      if (exact_v_hairpin_allows(ctx, parent.i, parent.j)) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::V_VM:
      if (exact_v_pairing_state_allows(ctx, parent.i, parent.j)) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::WMv_STEM_V:
    case grammar::RuleId::WMp_STEM_WMB:
      if (exact_wmv_wmp_stem_parent_allows(parent, ctx)) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::WIP_BASE_V:
      if (parent.i <= parent.j) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::WMp_EXTEND_UNPAIRED:
      if (exact_wmv_wmp_extend_parent_allows(parent, ctx)) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::WI_BASE_SINGLE:
      if (parent.i == parent.j && ctx.is_unpaired(parent.i)) {
        push_if_valid(grammar::RuleSplit{});
      }
      return;
    case grammar::RuleId::WI_SPLIT_WMB:
      for (const cand_pos_t k : ExactSplitPositions(parent.i, parent.j - ctx.turn() - 1, split_positions)) {
        grammar::RuleSplit split{};
        split.k = k;
        append_exact_basic_deduction<grammar::RuleId::WI_SPLIT_WMB>(out, parent, split);
      }
      return;
    default:
      append_exact_deduction_buffer(out, generate_basic_deductions(parent, rule, {}));
      return;
  }
}

inline std::vector<Deduction> generate_exact_basic_deductions(
    const ItemKey &parent,
    grammar::RuleId rule,
    const PairabilityContext &ctx,
    const std::vector<cand_pos_t> *split_positions = nullptr) {
  std::vector<Deduction> out;
  append_exact_basic_deductions(out, parent, rule, ctx, split_positions);
  return out;
}


template <class PairabilityContext>
inline bool exact_vp_public_item_active(const ItemKey &item,
                                        const PairabilityContext &ctx) {
  return item.nonterminal == grammar::NonTerminal::VP &&
         exact_vp_parent_pair_admissible(ctx, item.i, item.j) &&
         exact_round_open_family_active(ctx, item.i, item.j) &&
         exact_vp_public_endpoint_allows(ctx, item.i, item.j);
}

// Shared parent-side omission gate for exact generators and sparse schedules.
// False means the item-level exact dispatcher returns no deductions without
// inspecting child chart values. True does not imply a nonzero result.
template <class PairabilityContext>
inline bool exact_basic_item_parent_active(const ItemKey &item,
                                           const PairabilityContext &ctx) {
  using grammar::NonTerminal;
  switch (item.nonterminal) {
    case NonTerminal::V:
      return exact_v_family_active(ctx, item.i, item.j);
    case NonTerminal::WI:
    case NonTerminal::WIP:
      return exact_round_weakly_closed_family_active(ctx, item.i, item.j);
    case NonTerminal::VPL:
      return exact_round_open_family_active(ctx, item.i, item.j) &&
             ctx.pair_round(item.j) < 0;
    case NonTerminal::VPR: {
      const cand_pos_t pair_j = ctx.pair_round(item.j);
      return exact_round_open_family_active(ctx, item.i, item.j) &&
             pair_j < item.j;
    }
    case NonTerminal::VP:
      return exact_vp_public_item_active(item, ctx);
    case NonTerminal::WMBW:
      return exact_round_pseudoloop_region_active(ctx, item.i, item.j);
    case NonTerminal::WMB:
      return exact_completed_wmb_parent_active(ctx, item.i, item.j);
    case NonTerminal::WMBP:
      return exact_wmbp_item_active(item, ctx);
    default:
      return true;
  }
}

template <class PairabilityContext>
inline std::vector<Deduction> generate_exact_vp_case3_deductions_for_item(
    const ItemKey &item,
    const PairabilityContext &ctx) {
  if (!exact_vp_public_item_active(item, ctx)) {
    return {};
  }
  return generate_exact_basic_deductions(
      item, grammar::RuleId::VP_WI_CASE3, ctx);
}

template <class PairabilityContext>
inline void append_exact_basic_deductions_for_item(
    std::vector<Deduction> &out,
    const ItemKey &item,
    const PairabilityContext &ctx,
    const std::vector<cand_pos_t> *split_positions = nullptr,
    const ExactNumericSplitViews *numeric_splits = nullptr) {
  using grammar::NonTerminal;
  using grammar::RuleId;

  auto append_rule = [&](std::vector<Deduction> &out, RuleId rule) {
    append_exact_basic_deductions(out, item, rule, ctx,
        numeric_splits ? numeric_splits->for_rule(rule, split_positions) : split_positions);
  };

  switch (item.nonterminal) {
    case NonTerminal::V: {
      if (!exact_v_family_active(ctx, item.i, item.j)) {
        return;
      }
      append_rule(out, RuleId::V_HAIRPIN);
      append_rule(out, RuleId::V_INTERNAL);
      append_rule(out, RuleId::V_VM);
      return;
    }
    case NonTerminal::WMv: {
      append_rule(out, RuleId::WMv_STEM_V);
      append_rule(out, RuleId::WMv_EXTEND_UNPAIRED);
      return;
    }
    case NonTerminal::WMp: {
      append_rule(out, RuleId::WMp_STEM_WMB);
      append_rule(out, RuleId::WMp_EXTEND_UNPAIRED);
      return;
    }
    case NonTerminal::WI: {
      if (!exact_round_weakly_closed_family_active(ctx, item.i, item.j)) {
        return;
      }
      append_rule(out, RuleId::WI_BASE_SINGLE);
      append_rule(out, RuleId::WI_SPLIT_V);
      append_rule(out, RuleId::WI_SPLIT_WMB);
      append_rule(out, RuleId::WI_EXTEND_UNPAIRED);
      return;
    }
    case NonTerminal::WM: {
      append_rule(out, RuleId::WM_START_V);
      append_rule(out, RuleId::WM_START_WMB);
      append_rule(out, RuleId::WM_SPLIT_V);
      append_rule(out, RuleId::WM_SPLIT_WMB);
      append_rule(out, RuleId::WM_EXTEND_UNPAIRED);
      return;
    }
    case NonTerminal::VM: {
      append_rule(out, RuleId::VM_SPLIT_WM_WMv);
      append_rule(out, RuleId::VM_SPLIT_WM_WMp);
      append_rule(out, RuleId::VM_SPLIT_WMp_BASE);
      return;
    }
    case NonTerminal::WIP: {
      if (!exact_round_weakly_closed_family_active(ctx, item.i, item.j)) {
        return;
      }
      append_rule(out, RuleId::WIP_BASE_V);
      append_rule(out, RuleId::WIP_BASE_WMB);
      append_rule(out, RuleId::WIP_SPLIT_V);
      append_rule(out, RuleId::WIP_SPLIT_WMB);
      append_rule(out, RuleId::WIP_BASEPAIR_V);
      append_rule(out, RuleId::WIP_BASEPAIR_WMB);
      append_rule(out, RuleId::WIP_EXTEND_UNPAIRED);
      return;
    }
    case NonTerminal::VPL:
      if (!exact_round_open_family_active(ctx, item.i, item.j) ||
          ctx.pair_round(item.j) >= 0) {
        return;
      }
      append_rule(out, RuleId::VPL_SPLIT_VP);
      return;
    case NonTerminal::VPR: {
      const cand_pos_t pair_j = ctx.pair_round(item.j);
      if (!exact_round_open_family_active(ctx, item.i, item.j) ||
          !(pair_j < item.j)) {
        return;
      }
      append_rule(out, RuleId::VPR_SPLIT_VP_WIP);
      append_rule(out, RuleId::VPR_SPLIT_VP_BASEPAIR);
      return;
    }
    case NonTerminal::VP: {
      if (!exact_vp_public_item_active(item, ctx)) {
        return;
      }
      // These fixed-coordinate rules observe the same immutable scaffold
      // geometry. Share that observation while retaining the single-rule
      // expansion, admission and CASE1/CASE2/CASE3 stream order.
      if (exact_vp_open_endpoint_allows(ctx, item.i, item.j)) {
        const auto geometry = observe_vp_wi_geometry(ctx, item.i, item.j);
        for (const auto rule : {RuleId::VP_WI_CASE1, RuleId::VP_WI_CASE2,
                                RuleId::VP_WI_CASE3}) {
          if (const auto split = geometry.split_for(rule)) {
            const auto deduction = expand_basic_rule(item, rule, *split);
            if (is_valid(deduction)) out.push_back(deduction);
          }
        }
      }
      append_rule(out, RuleId::VP_STACK);
      append_rule(out, RuleId::VP_INTERNAL_LOOP);
      append_rule(out, RuleId::VP_WIP_VP_LEFT);
      append_rule(out, RuleId::VP_VP_WIP_RIGHT);
      append_rule(out, RuleId::VP_WIP_VPR);
      append_rule(out, RuleId::VP_VPL_WIP);
      return;
    }
    case NonTerminal::BE: {
      append_rule(out, RuleId::BE_BASE_SAMEPAIR);
      append_rule(out, RuleId::BE_STACK);
      append_rule(out, RuleId::BE_INTERNAL_LOOP);
      append_rule(out, RuleId::BE_WIP_WIP);
      append_rule(out, RuleId::BE_WIP_BASEPAIR);
      append_rule(out, RuleId::BE_BASEPAIR_WIP);
      return;
    }
    case NonTerminal::WMBW:
      if (!exact_round_pseudoloop_region_active(ctx, item.i, item.j)) {
        return;
      }
      append_rule(out, RuleId::WMBW_SPLIT_WMBP_WI);
      return;
    case NonTerminal::WMBP: {
      if (!exact_wmbp_item_active(item, ctx)) {
        return;
      }
      const auto *vp_starts = numeric_splits && numeric_splits->wmbp_vp_starts
          ? numeric_splits->wmbp_vp_starts : split_positions;
      auto family = vp_starts
          ? generate_exact_wmbp_be_vp_family_deductions_for_vp_starts(
                item, ctx, *vp_starts, true, true, true)
          : generate_exact_wmbp_be_vp_family_deductions(
                item, ctx, true, true, true);
      append_exact_deduction_buffer(out, std::move(family.wmbp));
      append_exact_deduction_buffer(out, std::move(family.wmbw));
      append_rule(out, RuleId::WMBP_DIRECT_VP);
      append_exact_deduction_buffer(out, std::move(family.wi));
      return;
    }
    case NonTerminal::WMB: {
      if (!exact_completed_wmb_parent_active(ctx, item.i, item.j)) {
        return;
      }
      append_rule(out, RuleId::WMB_SPLIT_BE_WMBP_WI);
      append_rule(out, RuleId::WMB_DIRECT_WMBP);
      append_rule(out, RuleId::WMB_EMPTY);
      return;
    }
    case NonTerminal::W: {
      if (item.i == item.j) {
        out.push_back(make_zero_child_deduction(
            item, RuleId::W_EXTEND_UNPAIRED, grammar::RuleSplit{}));
        return;
      }
      append_rule(out, RuleId::W_EXTEND_UNPAIRED);
      append_rule(out, RuleId::W_SPLIT_V);
      append_rule(out, RuleId::W_SPLIT_WMB);
      if (item.i == 1) {
        append_rule(out, RuleId::W_ROOT_WMB_DECOMPOSITION);
      }
      return;
    }
    default:
      return;
  }
}

template <class PairabilityContext>
inline std::vector<Deduction> generate_exact_basic_deductions_for_item(
    const ItemKey &item, const PairabilityContext &ctx,
    const std::vector<cand_pos_t> *split_positions = nullptr,
    const ExactNumericSplitViews *numeric_splits = nullptr) {
  std::vector<Deduction> out;
  append_exact_basic_deductions_for_item(out, item, ctx, split_positions, numeric_splits);
  return out;
}

template <class PairabilityContext>
inline std::vector<Deduction>
generate_exact_wmbp_consumer_deductions_for_retained_ends(
    const ItemKey &item,
    const PairabilityContext &ctx,
    const std::vector<cand_pos_t> &retained_wmbp_ends) {
  using grammar::NonTerminal;
  using grammar::RuleId;
  if (!exact_round_pseudoloop_region_active(ctx, item.i, item.j)) return {};

  std::vector<Deduction> out;
  if (item.nonterminal == NonTerminal::WMBW) {
    for (const cand_pos_t end : retained_wmbp_ends) {
      append_exact_wmbw_split_for_wmbp_end(out, item, ctx, end);
    }
    return out;
  }
  if (item.nonterminal != NonTerminal::WMB) return {};

  for (const cand_pos_t end : retained_wmbp_ends) {
    append_exact_wmb_split_for_wmbp_end(out, item, ctx, end);
  }
  const auto append_rule = [&](RuleId rule) {
    const auto generated = generate_exact_basic_deductions(item, rule, ctx);
    out.insert(out.end(), generated.begin(), generated.end());
  };
  if (std::binary_search(retained_wmbp_ends.begin(),
                         retained_wmbp_ends.end(),
                         item.j)) {
    append_rule(RuleId::WMB_DIRECT_WMBP);
  }
  append_rule(RuleId::WMB_EMPTY);
  return out;
}

}  // namespace cparty::scfg2::deductions

#endif  // CPARTY_SCFG2_DEDUCTIONS_GENERATE_EXACT_BASIC_HH

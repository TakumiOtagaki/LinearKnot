#ifndef CPARTY_SCFG2_DEDUCTIONS_WMB_GEOMETRY_HH
#define CPARTY_SCFG2_DEDUCTIONS_WMB_GEOMETRY_HH

#include "scfg2/deductions/exact_context.hh"
#include "scfg2/deductions/item_key.hh"

#include <optional>

namespace cparty::scfg2::deductions {

struct WmbpSharedBeVpGeometry {
  cand_pos_t split_l = -1;
  cand_pos_t be_i = -1;
  cand_pos_t be_j = -1;
  cand_pos_t be_ip = -1;
};

struct WmbpWiVpAnchors {
  cand_pos_t parent_i = -1;
  cand_pos_t be_ip = -1;
  cand_pos_t be_jp = -1;
  cand_pos_t vp_i = -1;
  cand_pos_t vp_j = -1;
};

struct WmbSplitBeWmbpWiAnchors {
  cand_pos_t be_jp = -1;
  cand_pos_t be_ip = -1;
};

inline bool wmbp_right_boundary_admissible(const PairabilityContext &ctx,
                                           cand_pos_t j) {
  const cand_pos_t pair_j = ctx.pair_partner(j);
  return pair_j < 0 || pair_j < j;
}

inline bool wmbw_parent_geometry_allows(const PairabilityContext &ctx,
                                        const ItemKey &parent,
                                        cand_pos_t split_l) {
  const cand_pos_t cover = ctx.parent_index(split_l);
  const bool shared_cover = cover >= 0 && cover == ctx.parent_index(parent.j);
  return shared_cover && ctx.pair_partner(split_l) < 0 &&
         ctx.pair_partner(parent.j) < parent.j;
}

// Compatibility support, separate from the paper's border geometry: a BE
// opening left of its WMBP parent was unavailable in the oracle's fill order.
// Keep this ownership restriction under every executor schedule.
inline bool wmbp_be_source_available(const ItemKey &parent, cand_pos_t be_i) {
  return be_i >= parent.i;
}

inline bool wmbp_shared_be_vp_parent_geometry_allows(const ItemKey &parent,
                                                     cand_pos_t split_l,
                                                     cand_pos_t parent_l,
                                                     cand_pos_t b_ij) {
  if (b_ij > 0) {
    return split_l < b_ij;
  }
  return b_ij < 0 && (parent_l >= 0 && parent_l < parent.j);
}

inline bool wmbp_shared_be_vp_child_geometry_allows(const PairabilityContext &ctx,
                                                    const ItemKey &parent,
                                                    cand_pos_t split_l,
                                                    cand_pos_t be_ip,
                                                    cand_pos_t split_be_open,
                                                    cand_pos_t parent_l) {
  return be_ip >= 0 &&
         split_l > be_ip &&
         split_be_open > 0 &&
         split_l < split_be_open &&
         parent_l >= parent.i &&
         parent_l < parent.j &&
         split_l + ctx.turn() <= parent.j;
}

inline bool wmbp_shared_be_vp_anchor_allows(const PairabilityContext &ctx,
                                            cand_pos_t be_j,
                                            cand_pos_t split_be_open,
                                            WmbpSharedBeVpGeometry *out) {
  const cand_pos_t be_i = ctx.pair_partner(be_j);
  const cand_pos_t be_open = ctx.pair_partner(split_be_open);
  if (!(be_i >= 1 && be_open >= 1)) {
    return false;
  }
  if (out) {
    out->split_l = split_be_open;
    out->be_i = be_i;
    out->be_j = be_j;
    out->be_ip = be_open;
  }
  return true;
}

inline std::optional<WmbpSharedBeVpGeometry> wmbp_shared_be_vp_geometry(
    const PairabilityContext &ctx,
    const ItemKey &parent,
    cand_pos_t split_l) {
  const cand_pos_t parent_l = ctx.parent_index(split_l);
  const cand_pos_t b_ij = ctx.border_b(parent.i, parent.j);
  if (!wmbp_shared_be_vp_parent_geometry_allows(parent, split_l, parent_l, b_ij)) {
    return std::nullopt;
  }
  const cand_pos_t be_ip = ctx.border_bp(parent.i, split_l);
  const cand_pos_t be_j = ctx.border_B(split_l, parent.j);
  const cand_pos_t split_be_open = ctx.border_Bp(split_l, parent.j);
  if (!wmbp_shared_be_vp_child_geometry_allows(
          ctx, parent, split_l, be_ip, split_be_open, parent_l)) {
    return std::nullopt;
  }
  WmbpSharedBeVpGeometry geometry;
  if (!wmbp_shared_be_vp_anchor_allows(ctx, be_j, split_be_open, &geometry)) {
    return std::nullopt;
  }
  if (!wmbp_be_source_available(parent, geometry.be_i)) {
    return std::nullopt;
  }
  return geometry;
}

inline bool wmbp_wi_vp_parent_geometry_allows(const PairabilityContext &ctx,
                                              const ItemKey &parent) {
  return wmbp_right_boundary_admissible(ctx, parent.j) &&
         ctx.pair_partner(parent.i) >= 0;
}

inline bool wmbp_wi_vp_child_geometry_allows(const PairabilityContext &ctx,
                                             const ItemKey &parent,
                                             cand_pos_t split_l,
                                             cand_pos_t be_ip,
                                             cand_pos_t parent_l) {
  return be_ip >= 0 &&
         be_ip < parent.j &&
         split_l + ctx.turn() <= parent.j &&
         parent_l >= parent.i &&
         parent_l < parent.j;
}

inline std::optional<WmbpWiVpAnchors> wmbp_wi_vp_anchor(
    const PairabilityContext &ctx,
    cand_pos_t parent_i,
    cand_pos_t be_ip,
    cand_pos_t vp_i,
    cand_pos_t vp_j) {
  const cand_pos_t be_jp = ctx.pair_partner(be_ip);
  if (be_jp < 0) {
    return std::nullopt;
  }
  return WmbpWiVpAnchors{
      .parent_i = parent_i,
      .be_ip = be_ip,
      .be_jp = be_jp,
      .vp_i = vp_i,
      .vp_j = vp_j,
  };
}

inline std::optional<WmbpWiVpAnchors> wmbp_wi_vp_geometry(
    const PairabilityContext &ctx,
    const ItemKey &parent,
    cand_pos_t split_l) {
  const cand_pos_t be_ip = ctx.border_bp(parent.i, split_l);
  const cand_pos_t parent_l = ctx.parent_index(split_l);
  if (!wmbp_wi_vp_child_geometry_allows(ctx, parent, split_l, be_ip, parent_l)) {
    return std::nullopt;
  }
  return wmbp_wi_vp_anchor(ctx, parent.i, be_ip, split_l, parent.j);
}

inline bool wmb_split_be_wmbp_wi_parent_geometry_allows(const PairabilityContext &ctx,
                                                        const ItemKey &parent) {
  return ctx.pair_partner(parent.j) >= 0 &&
         parent.j > ctx.pair_partner(parent.j) &&
         ctx.pair_partner(parent.j) > parent.i;
}

inline std::optional<WmbSplitBeWmbpWiAnchors> wmb_split_be_wmbp_wi_anchor(
    const PairabilityContext &ctx,
    cand_pos_t split_be_open) {
  const cand_pos_t inner_i = ctx.pair_partner(split_be_open);
  if (inner_i < 1) {
    return std::nullopt;
  }
  return WmbSplitBeWmbpWiAnchors{
      .be_jp = split_be_open,
      .be_ip = inner_i,
  };
}

inline std::optional<WmbSplitBeWmbpWiAnchors> wmb_split_be_wmbp_wi_geometry(
    const PairabilityContext &ctx,
    const ItemKey &parent,
    cand_pos_t split_l) {
  const cand_pos_t split_be_open = ctx.border_Bp(split_l, parent.j);
  // Nucleotides use the inclusive 1..n domain.  In particular, Bp(l,j)==n
  // is a valid right band endpoint, not an exterior sentinel.
  if (split_be_open < 0 || split_be_open > ctx.sequence_length()) {
    return std::nullopt;
  }
  return wmb_split_be_wmbp_wi_anchor(ctx, split_be_open);
}

}  // namespace cparty::scfg2::deductions

#endif  // CPARTY_SCFG2_DEDUCTIONS_WMB_GEOMETRY_HH

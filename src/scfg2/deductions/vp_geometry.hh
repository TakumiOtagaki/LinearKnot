#ifndef CPARTY_SCFG2_DEDUCTIONS_VP_GEOMETRY_HH
#define CPARTY_SCFG2_DEDUCTIONS_VP_GEOMETRY_HH

#include "scfg2/deductions/exact_context.hh"
#include "scfg2/grammar/rules.hh"

#include <algorithm>
#include <optional>

namespace cparty::scfg2::deductions {

// Read-only geometry for SI eq. (v), Cases 1--3. Parent-index and border
// sentinels are part of the current support contract; see README.md.
struct VpWiGeometry {
  cand_pos_t left_cover, right_cover;
  cand_pos_t Bp, B, b, bp;

  std::optional<grammar::RuleSplit> split_for(grammar::RuleId rule) const {
    const bool left_band = Bp >= 0 && B >= 0;
    const bool right_band = b >= 0 && bp >= 0;
    grammar::RuleSplit split{};
    switch (rule) {
      case grammar::RuleId::VP_WI_CASE1:
        if (!(left_cover > 0 && right_cover >= 0 && right_cover < left_cover &&
              left_band && bp < 0)) return std::nullopt;
        split.p = Bp; split.q = B;
        return split;
      case grammar::RuleId::VP_WI_CASE2:
        if (!(left_cover >= 0 && right_cover > 0 && left_cover < right_cover &&
              right_band && Bp < 0)) return std::nullopt;
        split.p = b; split.q = bp;
        return split;
      case grammar::RuleId::VP_WI_CASE3:
        if (!(left_cover > 0 && right_cover > 0 && left_band && right_band))
          return std::nullopt;
        split.p = Bp; split.q = B;
        split.k = b; split.l = bp;
        return split;
      default:
        return std::nullopt;
    }
  }
};

inline VpWiGeometry observe_vp_wi_geometry(const PairabilityContext &ctx,
                                          cand_pos_t i, cand_pos_t j) {
  return {ctx.parent_index(i), ctx.parent_index(j),
          ctx.border_Bp(i, j), ctx.border_B(i, j),
          ctx.border_b(i, j), ctx.border_bp(i, j)};
}

// SI eq. (v), Case 5, plus the existing search/sentinel contract. Derive
// bounds from the two unpaired gaps rather than from a shifted edge formula.
// The VP search uses u_left + u_right < MAXLOOP (an exclusive budget).
struct VpInternalLoopDomain {
  cand_pos_t first_left, last_left;
  cand_pos_t last_right, right_border;
  cand_pos_t unpaired_limit;

  cand_pos_t first_right(cand_pos_t left) const {
    const cand_pos_t left_gap = left - first_left;
    const cand_pos_t largest_right_gap = unpaired_limit - 1 - left_gap;
    return std::max(right_border + 1, last_right - largest_right_gap);
  }
};

inline VpInternalLoopDomain vp_internal_loop_domain(
    const PairabilityContext &ctx, cand_pos_t i, cand_pos_t j,
    cand_pos_t unpaired_limit) {
  const auto Bp = ctx.border_Bp(i, j);
  const auto b = ctx.border_b(i, j);
  cand_pos_t last_left = std::min(i + unpaired_limit, j - ctx.turn() - 1);
  // Preserve the context's missing/zero-border convention. Applying just
  // one positive border here would alter existing deduction support.
  if (Bp > 0 && b > 0) last_left = std::min(last_left, std::min(Bp, b) - 1);
  return {i + 1, last_left, j - 1,
          std::max(ctx.border_bp(i, j), ctx.border_B(i, j)), unpaired_limit};
}

}  // namespace cparty::scfg2::deductions

#endif

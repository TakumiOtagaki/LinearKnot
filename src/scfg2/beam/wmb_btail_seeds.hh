#ifndef CPARTY_SCFG2_BEAM_WMB_BTAIL_SEEDS_HH
#define CPARTY_SCFG2_BEAM_WMB_BTAIL_SEEDS_HH

#include "scfg2/deductions/item_key.hh"
#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

struct WmbBtailSeedCounts { std::size_t direct = 0, joined = 0; };

// Shared deduction surface of the derived BTAIL factor. No weights or
// products are owned here: each seed is WMBP(i,k), with optional WI(k+1,t-1).
template <class Context, class WmbpStartsFn, class EmitFn>
WmbBtailSeedCounts for_each_wmb_btail_seed(
    cand_pos_t t, const Context &ctx, const std::vector<cand_pos_t> &wi_starts,
    WmbpStartsFn &&wmbp_starts, EmitFn &&emit) {
  WmbBtailSeedCounts counts;
  if (t < 2 || ctx.pair_partner(t) < 1 || ctx.pair_partner(t) >= t) return counts;
  const cand_pos_t right = t - 1;
  if (ctx.border_Bp(right, t) == t) {
    for (const cand_pos_t i : wmbp_starts(right)) {
      ++counts.direct;
      emit(i, right, false);
    }
  }
  for (const cand_pos_t wi_i : wi_starts) {
    const cand_pos_t k = wi_i - 1;
    if (k < 1 || k >= right || ctx.border_Bp(k, t) != t) continue;
    for (const cand_pos_t i : wmbp_starts(k)) {
      ++counts.joined;
      emit(i, k, true);
    }
  }
  return counts;
}

}  // namespace cparty::scfg2::beam
#endif

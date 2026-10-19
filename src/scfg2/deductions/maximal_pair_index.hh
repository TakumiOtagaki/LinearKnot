#ifndef CPARTY_SCFG2_DEDUCTIONS_MAXIMAL_PAIR_INDEX_HH
#define CPARTY_SCFG2_DEDUCTIONS_MAXIMAL_PAIR_INDEX_HH

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/geometry/static_scaffold_geometry.hh"

#include <limits>
#include <vector>

namespace cparty::scfg2::deductions {

// O(n) setup/storage, O(1) geometry query for a static laminar scaffold.
// RMQ selects the first completed pair; boundary-gap LCAs lift it to the
// first maximal contained pair. A suffix RMQ detects a second maximal pair.
// See scfg2_constant_time_geometry_20261009.md for the reduction proof.
// Dual-partner or crossing contexts are unusable: callers retain the scan.
class MaximalPairIndex {
 public:
  MaximalPairIndex() = default;

  explicit MaximalPairIndex(const PairabilityContext &ctx) {
    const cand_pos_t n = ctx.sequence_length();
    if (n < 0) return;
    std::vector<cand_pos_t> partners(static_cast<std::size_t>(n) + 1, -2);
    std::vector<cand_pos_t> closing(partners.size(), std::numeric_limits<cand_pos_t>::max());
    for (cand_pos_t pos = 1; pos <= n; ++pos) {
      const cand_pos_t round = ctx.pair_round(pos), square = ctx.pair_square(pos);
      if (round > 0 && square > 0) return;
      const cand_pos_t pair = round > 0 ? round : square;
      partners[pos] = pair;
      if (pair > pos) closing[pos] = pair;
    }
    scaffold_ = geometry::StaticScaffoldGeometry(std::move(partners));
    if (!scaffold_.usable()) return;
    closing_min_ = geometry::StaticRangeMinimum<cand_pos_t>(closing);
    usable_ = true;
  }

  bool usable() const { return usable_; }

  ExactDirectChildGeometry geometry(cand_pos_t i, cand_pos_t j) const {
    ExactDirectChildGeometry geom;
    geom.count = 0;
    assert(usable_ && i >= 1 && i < j && j <= scaffold_.size());
    // Pairs touching one endpoint enclose/exclude their contents, while (i,j)
    // itself is ignored by the scan. Preserve that asymmetry before indexing.
    cand_pos_t lo = i, hi = j;
    const cand_pos_t pj = scaffold_.partner(j);
    if (pj > i && pj < j) hi = pj;
    const cand_pos_t pi = scaffold_.partner(i);
    if (pi > i && pi < hi) lo = pi;
    const cand_pos_t first_completed = completed_pair(lo, hi);
    if (first_completed < 0) return geom;

    cand_pos_t first = first_completed;
    const cand_pos_t parent = scaffold_.parent(first_completed);
    // If even the immediate parent is not fully contained, no higher ancestor
    // can be contained. Wide scaffolds therefore need no LCA/child queries.
    if (parent > lo && scaffold_.partner(parent) < hi) {
      // Ancestors shared with either boundary gap cross that boundary. The
      // deeper LCA is the last excluded ancestor on the pair's ancestor path.
      const cand_pos_t left = scaffold_.lca(scaffold_.gap_after(lo), first_completed);
      const cand_pos_t right = scaffold_.lca(first_completed, scaffold_.gap_before(hi));
      const cand_pos_t excluded = scaffold_.depth(left) >= scaffold_.depth(right) ? left : right;
      first = scaffold_.child_below(excluded, first_completed);
    }
    const cand_pos_t close = scaffold_.partner(first);
    if (completed_pair(close, hi) >= 0) {
      geom.count = 2;  // capped: consumers only distinguish zero/one/multiple
      return geom;
    }
    geom.count = 1;
    geom.child_i = first;
    geom.child_j = close;
    geom.left_gap = first - i - 1;
    geom.right_gap = j - close - 1;
    return geom;
  }

  std::size_t storage_bytes() const {
    return scaffold_.storage_bytes() + closing_min_.storage_bytes();
  }

 private:
  // The opening whose closing endpoint is smallest inside (lo,hi), if any.
  cand_pos_t completed_pair(cand_pos_t lo, cand_pos_t hi) const {
    if (hi - lo <= 1) return -1;
    const cand_pos_t opening = static_cast<cand_pos_t>(closing_min_.argmin(lo + 1, hi - 1));
    const cand_pos_t close = scaffold_.partner(opening);
    return close > opening && close < hi ? opening : -1;
  }

  bool usable_ = false;
  geometry::StaticScaffoldGeometry scaffold_;
  geometry::StaticRangeMinimum<cand_pos_t> closing_min_;
};

}  // namespace cparty::scfg2::deductions
#endif

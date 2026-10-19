#ifndef CPARTY_SCFG2_BEAM_WMBP_BE_VP_SPLIT_CHILD_INDEX_HH
#define CPARTY_SCFG2_BEAM_WMBP_BE_VP_SPLIT_CHILD_INDEX_HH

#include "scfg2/deductions/generate_exact_basic.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

namespace cparty::scfg2::beam {

struct WmbpBeVpIndexedSplitQuery {
  std::vector<cand_pos_t> positions;
  std::size_t parent_bucket_visits = 0;
};

// Scaffold-only index for the common parent-side conditions of the shared-BE
// and WI WMBP split families. Exact border/anchor admissibility remains in the
// deduction geometry functions.
class ExactWmbpBeVpSplitChildIndex {
 public:
  template <class Context>
  static ExactWmbpBeVpSplitChildIndex build(cand_pos_t n,
                                             const Context &ctx) {
    ExactWmbpBeVpSplitChildIndex out;
    for (cand_pos_t l = 1; l <= n; ++l) {
      const cand_pos_t parent = ctx.parent_index(l);
      if (parent >= 1) out.positions_by_parent_[parent].push_back(l);
    }
    return out;
  }

  template <class Context>
  WmbpBeVpIndexedSplitQuery query(const deductions::ItemKey &parent,
                                  const Context &ctx) const {
    WmbpBeVpIndexedSplitQuery out;
    const cand_pos_t lower_l = parent.i + 1;
    const cand_pos_t upper_l = parent.j - ctx.turn();
    if (upper_l < lower_l) return out;
    for (auto bucket = positions_by_parent_.lower_bound(parent.i);
         bucket != positions_by_parent_.end() && bucket->first < parent.j;
         ++bucket) {
      ++out.parent_bucket_visits;
      const auto begin =
          std::lower_bound(bucket->second.begin(), bucket->second.end(), lower_l);
      const auto end =
          std::upper_bound(bucket->second.begin(), bucket->second.end(), upper_l);
      out.positions.insert(out.positions.end(), begin, end);
    }
    std::sort(out.positions.begin(), out.positions.end());
    return out;
  }

 private:
  std::map<cand_pos_t, std::vector<cand_pos_t>> positions_by_parent_;
};

struct IndexedWmbpBeVpSplitScan {
  deductions::ExactWmbpBeVpSplitScan geometry;
  std::size_t parent_bucket_visits = 0;
};

template <class Context>
IndexedWmbpBeVpSplitScan scan_exact_wmbp_be_vp_split_geometries_from_index(
    const deductions::ItemKey &parent,
    const Context &ctx,
    const ExactWmbpBeVpSplitChildIndex &index,
    bool scan_shared,
    bool scan_wi) {
  IndexedWmbpBeVpSplitScan out;
  const bool shared_active =
      scan_shared && deductions::wmbp_right_boundary_admissible(ctx, parent.j);
  const bool wi_active =
      scan_wi && deductions::wmbp_wi_vp_parent_geometry_allows(ctx, parent);
  if (!shared_active && !wi_active) return out;
  auto query = index.query(parent, ctx);
  out.parent_bucket_visits = query.parent_bucket_visits;
  for (const cand_pos_t l : query.positions) {
    ++out.geometry.l_visits;
    ++out.geometry.geometry_visits;
    const auto shared = shared_active
        ? deductions::wmbp_shared_be_vp_geometry(ctx, parent, l)
        : std::optional<deductions::WmbpSharedBeVpGeometry>{};
    const auto wi = wi_active
        ? deductions::wmbp_wi_vp_geometry(ctx, parent, l)
        : std::optional<deductions::WmbpWiVpAnchors>{};
    if (shared.has_value()) out.geometry.shared.emplace_back(l, *shared);
    if (wi.has_value()) out.geometry.wi.emplace_back(l, *wi);
    if (shared.has_value() && wi.has_value()) {
      ++out.geometry.overlap_count;
    }
  }
  return out;
}

inline bool same_wmbp_be_vp_split_geometry_scans(
    const deductions::ExactWmbpBeVpSplitScan &lhs,
    const deductions::ExactWmbpBeVpSplitScan &rhs) {
  if (lhs.shared.size() != rhs.shared.size() ||
      lhs.wi.size() != rhs.wi.size() ||
      lhs.overlap_count != rhs.overlap_count) {
    return false;
  }
  for (std::size_t idx = 0; idx < lhs.shared.size(); ++idx) {
    const auto &[lhs_l, lhs_g] = lhs.shared[idx];
    const auto &[rhs_l, rhs_g] = rhs.shared[idx];
    if (lhs_l != rhs_l || lhs_g.split_l != rhs_g.split_l ||
        lhs_g.be_i != rhs_g.be_i || lhs_g.be_j != rhs_g.be_j ||
        lhs_g.be_ip != rhs_g.be_ip) {
      return false;
    }
  }
  for (std::size_t idx = 0; idx < lhs.wi.size(); ++idx) {
    const auto &[lhs_l, lhs_g] = lhs.wi[idx];
    const auto &[rhs_l, rhs_g] = rhs.wi[idx];
    if (lhs_l != rhs_l || lhs_g.parent_i != rhs_g.parent_i ||
        lhs_g.be_ip != rhs_g.be_ip || lhs_g.be_jp != rhs_g.be_jp ||
        lhs_g.vp_i != rhs_g.vp_i || lhs_g.vp_j != rhs_g.vp_j) {
      return false;
    }
  }
  return true;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMBP_BE_VP_SPLIT_CHILD_INDEX_HH

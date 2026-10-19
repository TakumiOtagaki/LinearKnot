#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_BORDER_PROFILE_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_BORDER_PROFILE_HH

#include "scfg2/beam/vp_case3_tail_column_plan.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase3TailBorderProfile {
  std::size_t source_count = 0;
  std::size_t tail_count = 0;
  std::size_t geometry_count = 0;
  std::size_t contiguous_i_run_count = 0;
  std::size_t constant_p_segment_count = 0;
  std::size_t exact_segment_count = 0;
  std::size_t discontinuous_tail_count = 0;
  std::size_t multi_p_tail_count = 0;
  std::size_t maximum_runs_per_tail = 0;
  std::size_t maximum_geometries_per_tail = 0;
};

inline VpCase3TailBorderProfile profile_vp_case3_tail_borders(
    const ExactVpCase3TailColumnPlan &plan) {
  struct GeometryInfo {
    cand_pos_t p = -1;
    unsigned carrier_mask = 0;
  };
  std::unordered_map<VpCase3TailKey,
                     std::map<cand_pos_t, GeometryInfo>,
                     VpCase3TailKeyHash>
      geometries;
  for (const auto &source : plan.sources()) {
    auto &info = geometries[source.tail][source.deduction.parent.i];
    if (info.p == -1) info.p = source.deduction.split.p;
    const auto nt = source.deduction.parent.nonterminal;
    const unsigned bit = nt == grammar::NonTerminal::VP ? 2u : 0u;
    info.carrier_mask |= bit;
  }

  VpCase3TailBorderProfile profile;
  profile.source_count = plan.source_count();
  profile.tail_count = geometries.size();
  for (const auto &[tail, points] : geometries) {
    (void)tail;
    profile.geometry_count += points.size();
    profile.maximum_geometries_per_tail =
        std::max(profile.maximum_geometries_per_tail, points.size());
    std::size_t runs = 0;
    std::size_t p_segments = 0;
    std::size_t exact_segments = 0;
    std::set<cand_pos_t> distinct_p;
    cand_pos_t previous_i = -2;
    cand_pos_t previous_p = -2;
    unsigned previous_mask = 0;
    for (const auto &[i, info] : points) {
      const bool new_run = i != previous_i + 1;
      if (new_run) ++runs;
      if (new_run || info.p != previous_p) ++p_segments;
      if (new_run || info.p != previous_p ||
          info.carrier_mask != previous_mask) {
        ++exact_segments;
      }
      distinct_p.insert(info.p);
      previous_i = i;
      previous_p = info.p;
      previous_mask = info.carrier_mask;
    }
    profile.contiguous_i_run_count += runs;
    profile.constant_p_segment_count += p_segments;
    profile.exact_segment_count += exact_segments;
    if (runs > 1) ++profile.discontinuous_tail_count;
    if (distinct_p.size() > 1) ++profile.multi_p_tail_count;
    profile.maximum_runs_per_tail =
        std::max(profile.maximum_runs_per_tail, runs);
  }
  return profile;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_BORDER_PROFILE_HH

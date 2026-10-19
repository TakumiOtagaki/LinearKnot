#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_ATTACHMENT_SEGMENTS_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_ATTACHMENT_SEGMENTS_HH

#include "scfg2/beam/vp_case3_tail_column_plan.hh"

#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

inline constexpr unsigned vp_carrier_bit = 2u;

inline constexpr unsigned vp_case3_carrier_bit(grammar::NonTerminal nt) {
  return nt == grammar::NonTerminal::VP ? vp_carrier_bit : 0u;
}

struct VpCase3TailAttachmentSegment {
  VpCase3TailKey tail{};
  cand_pos_t i_begin = -1;
  cand_pos_t i_end = -1;
  cand_pos_t p = -1;
  unsigned carrier_mask = 0;
};

inline std::vector<VpCase3TailAttachmentSegment>
build_vp_case3_tail_attachment_segments(
    const ExactVpCase3TailColumnPlan &plan) {
  struct GeometryInfo {
    cand_pos_t p = -1;
    unsigned carrier_mask = 0;
  };
  std::unordered_map<VpCase3TailKey,
                     std::map<cand_pos_t, GeometryInfo>,
                     VpCase3TailKeyHash>
      by_tail;
  for (const auto &source : plan.sources()) {
    auto &info = by_tail[source.tail][source.deduction.parent.i];
    if (info.p == -1) info.p = source.deduction.split.p;
    info.carrier_mask |=
        vp_case3_carrier_bit(source.deduction.parent.nonterminal);
  }

  std::vector<VpCase3TailAttachmentSegment> out;
  for (const auto &[tail, geometries] : by_tail) {
    VpCase3TailAttachmentSegment current{};
    bool active = false;
    for (const auto &[i, info] : geometries) {
      if (!active || i != current.i_end + 1 || info.p != current.p ||
          info.carrier_mask != current.carrier_mask) {
        if (active) out.push_back(current);
        current = VpCase3TailAttachmentSegment{
            tail, i, i, info.p, info.carrier_mask};
        active = true;
      } else {
        current.i_end = i;
      }
    }
    if (active) out.push_back(current);
  }
  return out;
}

inline bool vp_case3_tail_segments_reconstruct_sources(
    const ExactVpCase3TailColumnPlan &plan,
    const std::vector<VpCase3TailAttachmentSegment> &segments) {
  using Identity =
      std::tuple<unsigned, cand_pos_t, cand_pos_t, cand_pos_t, cand_pos_t,
                 cand_pos_t, cand_pos_t>;
  std::set<Identity> expected;
  for (const auto &source : plan.sources()) {
    expected.emplace(
        vp_case3_carrier_bit(source.deduction.parent.nonterminal),
        source.deduction.parent.i,
        source.tail.j,
        source.deduction.split.p,
        source.tail.q,
        source.tail.k,
        source.tail.l);
  }
  std::set<Identity> reconstructed;
  for (const auto &segment : segments) {
    for (cand_pos_t i = segment.i_begin; i <= segment.i_end; ++i) {
      for (const unsigned bit : {vp_carrier_bit}) {
        if ((segment.carrier_mask & bit) == 0) continue;
        reconstructed.emplace(bit,
                              i,
                              segment.tail.j,
                              segment.p,
                              segment.tail.q,
                              segment.tail.k,
                              segment.tail.l);
      }
    }
  }
  return expected.size() == plan.source_count() &&
         expected == reconstructed;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_ATTACHMENT_SEGMENTS_HH

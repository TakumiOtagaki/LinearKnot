#ifndef CPARTY_SCFG2_BEAM_VIRTUAL_CANDIDATE_HH
#define CPARTY_SCFG2_BEAM_VIRTUAL_CANDIDATE_HH

#include "scfg2/deductions/deduction.hh"

#include <array>
#include <cstddef>
#include <optional>
#include <string_view>

namespace cparty::scfg2::beam {

// Diagnostic-only families used to measure beam-ready factorizations before
// they become grammar nonterminals. They do not participate in admissibility,
// chart identity, local weights, or schedule order.
enum class VirtualBeamFamily : unsigned char {
  VpCase3Boundary,
  VpCase3LeftTail,
  VpCase3Left,
  VpCase3Other,
  BeCore,
  BeRightTail,
  BeLeftRightTail,
};

inline constexpr std::string_view virtual_beam_family_name(
    VirtualBeamFamily family) {
  switch (family) {
    case VirtualBeamFamily::VpCase3Boundary: return "VP_CASE3_BOUNDARY";
    case VirtualBeamFamily::VpCase3LeftTail: return "VP_CASE3_LTAIL";
    case VirtualBeamFamily::VpCase3Left: return "VP_CASE3_LEFT";
    case VirtualBeamFamily::VpCase3Other: return "VP_CASE3_OTHER";
    case VirtualBeamFamily::BeCore: return "BE_CORE";
    case VirtualBeamFamily::BeRightTail: return "BE_RTAIL";
    case VirtualBeamFamily::BeLeftRightTail: return "BE_LRTAIL";
  }
  return "UNKNOWN";
}

enum class VpCase3Partition : unsigned char {
  Boundary,
  LeftNonempty,
  Other,
};

struct VirtualBeamBucketKey {
  VirtualBeamFamily family = VirtualBeamFamily::VpCase3Boundary;
  cand_pos_t right_boundary = -1;
};

inline constexpr bool operator==(const VirtualBeamBucketKey &lhs,
                                 const VirtualBeamBucketKey &rhs) {
  return lhs.family == rhs.family &&
         lhs.right_boundary == rhs.right_boundary;
}

// FullIdentity deliberately keeps the complete source item and split tuple.
// A future beam bucket may be coarse, but observation must not lose the data
// needed to recover the original direct derivation exactly once.
struct VirtualBeamFullIdentity {
  VirtualBeamFamily family = VirtualBeamFamily::VpCase3Boundary;
  deductions::ItemKey source_parent{};
  grammar::RuleSplit source_split{};
};

inline bool operator==(const VirtualBeamFullIdentity &lhs,
                       const VirtualBeamFullIdentity &rhs) {
  return lhs.family == rhs.family &&
         deductions::ItemKeyEq{}(lhs.source_parent, rhs.source_parent) &&
         lhs.source_split.k == rhs.source_split.k &&
         lhs.source_split.l == rhs.source_split.l &&
         lhs.source_split.p == rhs.source_split.p &&
         lhs.source_split.q == rhs.source_split.q;
}

struct VirtualBeamCandidate {
  VirtualBeamBucketKey bucket{};
  VirtualBeamFullIdentity full_identity{};
  bool parent_facing = false;
};

struct VirtualBeamCandidateSet {
  std::array<VirtualBeamCandidate, 3> candidates{};
  std::size_t size = 0;
};

inline constexpr std::optional<VpCase3Partition> vp_case3_partition(
    const deductions::Deduction &deduction) {
  if (deduction.rule != grammar::RuleId::VP_WI_CASE3) {
    return std::nullopt;
  }
  const bool left_empty = deduction.parent.i + 1 > deduction.split.p - 1;
  const bool middle_empty = deduction.split.q + 1 > deduction.split.k - 1;
  if (left_empty && middle_empty) {
    return VpCase3Partition::Boundary;
  }
  if (!left_empty) {
    return VpCase3Partition::LeftNonempty;
  }
  return VpCase3Partition::Other;
}

inline constexpr VirtualBeamFamily primary_vp_case3_family(
    VpCase3Partition partition) {
  switch (partition) {
    case VpCase3Partition::Boundary:
      return VirtualBeamFamily::VpCase3Boundary;
    case VpCase3Partition::LeftNonempty:
      return VirtualBeamFamily::VpCase3Left;
    case VpCase3Partition::Other:
      return VirtualBeamFamily::VpCase3Other;
  }
  return VirtualBeamFamily::VpCase3Other;
}

inline constexpr VirtualBeamCandidate make_virtual_candidate(
    const deductions::Deduction &deduction,
    VirtualBeamFamily family,
    bool parent_facing,
    cand_pos_t right_boundary) {
  return VirtualBeamCandidate{
      VirtualBeamBucketKey{family, right_boundary},
      VirtualBeamFullIdentity{family, deduction.parent, deduction.split},
      parent_facing};
}

inline constexpr VirtualBeamCandidate make_virtual_candidate(
    const deductions::Deduction &deduction,
    VirtualBeamFamily family,
    bool parent_facing) {
  return make_virtual_candidate(
      deduction, family, parent_facing, deduction.parent.j);
}

inline constexpr VirtualBeamCandidateSet virtual_beam_candidates_for_deduction(
    const deductions::Deduction &deduction) {
  VirtualBeamCandidateSet out{};
  if (deduction.rule == grammar::RuleId::BE_WIP_WIP) {
    // The recursive child ends at split.l. The right-tail and complete frame
    // end at the source parent boundary. All three retain the anchored BE
    // parent and full split in their diagnostic full identity.
    out.candidates[out.size++] = make_virtual_candidate(
        deduction, VirtualBeamFamily::BeCore, false, deduction.split.l);
    out.candidates[out.size++] = make_virtual_candidate(
        deduction, VirtualBeamFamily::BeRightTail, false);
    out.candidates[out.size++] = make_virtual_candidate(
        deduction, VirtualBeamFamily::BeLeftRightTail, true);
    return out;
  }

  const auto partition = vp_case3_partition(deduction);
  if (!partition.has_value()) {
    return out;
  }

  if (*partition == VpCase3Partition::LeftNonempty) {
    // LTAIL is the private algebraic factor; LEFT is the parent-facing kept
    // object. Both retain the full original split while they are diagnostic.
    out.candidates[out.size++] = make_virtual_candidate(
        deduction, VirtualBeamFamily::VpCase3LeftTail, false);
  }
  out.candidates[out.size++] = make_virtual_candidate(
      deduction, primary_vp_case3_family(*partition), true);
  return out;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VIRTUAL_CANDIDATE_HH

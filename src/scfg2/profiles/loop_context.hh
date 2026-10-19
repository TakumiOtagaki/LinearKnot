#ifndef CPARTY_SCFG2_PROFILES_LOOP_CONTEXT_HH
#define CPARTY_SCFG2_PROFILES_LOOP_CONTEXT_HH

#include "scfg2/deductions/deduction.hh"

#include <array>
#include <cstddef>
#include <string_view>

namespace cparty::scfg2::profiles {

// Structural context of one base in G u G' (CapR contexts plus PseudoLoop).
// Every unpaired base is emitted by exactly one deduction in a derivation, and
// that deduction's rule fixes the context.  Paired bases are Stem.
enum class LoopContext : unsigned char {
  Stem,
  Hairpin,
  Bulge,
  Internal,
  Multi,
  Exterior,
  PseudoLoop,
};

inline constexpr std::size_t loop_context_count = 7;
inline constexpr std::size_t unpaired_loop_context_count = loop_context_count - 1;

inline constexpr std::string_view loop_context_name(LoopContext context) {
  switch (context) {
    case LoopContext::Stem: return "Stem";
    case LoopContext::Hairpin: return "Hairpin";
    case LoopContext::Bulge: return "Bulge";
    case LoopContext::Internal: return "Internal";
    case LoopContext::Multi: return "Multi";
    case LoopContext::Exterior: return "Exterior";
    case LoopContext::PseudoLoop: return "PseudoLoop";
  }
  return "?";
}

// A run [lo, hi] of unpaired bases emitted by one deduction.
struct UnpairedRun {
  cand_pos_t lo = 0;
  cand_pos_t hi = -1;
  LoopContext context = LoopContext::Exterior;
};

struct UnpairedEmission {
  std::array<UnpairedRun, 2> runs{};
  std::size_t count = 0;
  // True for rules that add one base to the right of a same-state chain
  // X(i,j-1) -> X(i,j), and for the one-base chain starts W(i,i) and
  // WI(i,i).  Their runs are pieces of a longer maximal run.
  bool chain_step = false;
};

namespace detail {

inline void add_run(UnpairedEmission &out,
                    cand_pos_t lo,
                    cand_pos_t hi,
                    LoopContext context) {
  if (lo > hi) return;
  out.runs[out.count++] = UnpairedRun{lo, hi, context};
}

// A two-sided loop closed by an outer and an inner pair.
inline void add_two_sided_loop(UnpairedEmission &out,
                               cand_pos_t i,
                               cand_pos_t j,
                               cand_pos_t k,
                               cand_pos_t l) {
  const bool left = k - i - 1 > 0;
  const bool right = j - l - 1 > 0;
  if (!left && !right) return;
  const LoopContext context =
      left && right ? LoopContext::Internal : LoopContext::Bulge;
  add_run(out, i + 1, k - 1, context);
  add_run(out, l + 1, j - 1, context);
}

}  // namespace detail

// The unpaired bases a deduction emits itself (not through its children).
// `is_unpaired(pos)` tells whether pos is unpaired in the scaffold G.
template <class IsUnpairedFn>
UnpairedEmission unpaired_emission(const deductions::Deduction &deduction,
                                   IsUnpairedFn &&is_unpaired) {
  using grammar::RuleId;
  UnpairedEmission out;
  const cand_pos_t i = deduction.parent.i;
  const cand_pos_t j = deduction.parent.j;
  const auto &s = deduction.split;
  switch (deduction.rule) {
    case RuleId::V_HAIRPIN:
      detail::add_run(out, i + 1, j - 1, LoopContext::Hairpin);
      break;
    case RuleId::V_INTERNAL:
    case RuleId::VP_INTERNAL_LOOP:
    case RuleId::BE_INTERNAL_LOOP:
      detail::add_two_sided_loop(out, i, j, s.k, s.l);
      break;

    case RuleId::W_EXTEND_UNPAIRED:
      // W(i,i) is a zero-child chain start; it emits i only when i is free.
      if (deduction.child_count == 0 && !is_unpaired(j)) break;
      detail::add_run(out, j, j, LoopContext::Exterior);
      out.chain_step = true;
      break;

    case RuleId::VM_SPLIT_WMp_BASE:
      detail::add_run(out, i + 1, s.k - 1, LoopContext::Multi);
      break;
    case RuleId::WM_START_V:
    case RuleId::WM_START_WMB:
      detail::add_run(out, i, s.k - 1, LoopContext::Multi);
      break;
    case RuleId::WM_EXTEND_UNPAIRED:
    case RuleId::WMv_EXTEND_UNPAIRED:
    case RuleId::WMp_EXTEND_UNPAIRED:
      detail::add_run(out, j, j, LoopContext::Multi);
      out.chain_step = true;
      break;

    // Pseudoloop: bases charged by the pseudoknot penalties (PUP inside WI,
    // c' in band-spanning multiloops and band borders).
    case RuleId::WI_BASE_SINGLE:
    case RuleId::WI_EXTEND_UNPAIRED:
    case RuleId::WIP_EXTEND_UNPAIRED:
      detail::add_run(out, j, j, LoopContext::PseudoLoop);
      out.chain_step = true;
      break;
    case RuleId::WIP_BASEPAIR_V:
    case RuleId::WIP_BASEPAIR_WMB:
    case RuleId::VPL_SPLIT_VP:
      detail::add_run(out, i, s.k - 1, LoopContext::PseudoLoop);
      break;
    case RuleId::VPR_SPLIT_VP_BASEPAIR:
      detail::add_run(out, s.k + 1, j, LoopContext::PseudoLoop);
      break;
    case RuleId::BE_WIP_BASEPAIR:
      detail::add_run(out, s.l + 1, j - 1, LoopContext::PseudoLoop);
      break;
    case RuleId::BE_BASEPAIR_WIP:
      detail::add_run(out, i + 1, s.k - 1, LoopContext::PseudoLoop);
      break;

    default:
      break;
  }
  return out;
}

}  // namespace cparty::scfg2::profiles

#endif  // CPARTY_SCFG2_PROFILES_LOOP_CONTEXT_HH

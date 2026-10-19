#ifndef CPARTY_SCFG2_BEAM_BE_FRAME_REUSE_PROFILE_HH
#define CPARTY_SCFG2_BEAM_BE_FRAME_REUSE_PROFILE_HH

#include "scfg2/deductions/deduction.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <tuple>

namespace cparty::scfg2::beam {

// Observation-only identities for the exact factorization
//
//   WIP(left) * BE(core) * WIP(right)
//
// of BE_WIP_WIP. They retain all anchored coordinates needed to recover the
// original source exactly; the coarser right-boundary buckets are measured
// separately below.
using BeCoreIdentity =
    std::tuple<cand_pos_t, cand_pos_t, cand_pos_t, cand_pos_t>;
using BeRightTailIdentity =
    std::tuple<cand_pos_t, cand_pos_t, cand_pos_t, cand_pos_t, cand_pos_t>;
using BeFrameIdentity = std::tuple<cand_pos_t,
                                   cand_pos_t,
                                   cand_pos_t,
                                   cand_pos_t,
                                   cand_pos_t,
                                   cand_pos_t>;

struct BeFrameReuseProfile {
  std::size_t source_count = 0;
  std::size_t core_count = 0;
  std::size_t reused_core_count = 0;
  std::size_t maximum_sources_per_core = 0;
  std::size_t maximum_right_tails_per_core = 0;
  std::size_t right_tail_count = 0;
  std::size_t reused_right_tail_count = 0;
  std::size_t maximum_sources_per_right_tail = 0;
  std::size_t frame_count = 0;
  std::size_t core_column_count = 0;
  std::size_t multi_core_column_count = 0;
  std::size_t maximum_cores_per_column = 0;
  std::size_t right_tail_column_count = 0;
  std::size_t multi_right_tail_column_count = 0;
  std::size_t maximum_right_tails_per_column = 0;
  std::size_t frame_column_count = 0;
  std::size_t multi_frame_column_count = 0;
  std::size_t maximum_frames_per_column = 0;
  std::size_t recursive_source_count = 0;
  std::size_t shared_core_count = 0;
  std::size_t reused_shared_core_count = 0;
  std::size_t maximum_recursive_sources_per_core = 0;
  std::size_t right_shell_source_count = 0;
  std::size_t shared_right_tail_count = 0;
  std::size_t reused_shared_right_tail_count = 0;
  std::size_t maximum_right_shell_sources_per_tail = 0;
};

class BeFrameReuseProfiler {
 public:
  void observe(const deductions::Deduction &deduction) {
    const bool recursive =
        deduction.rule == grammar::RuleId::BE_INTERNAL_LOOP ||
        deduction.rule == grammar::RuleId::BE_WIP_WIP ||
        deduction.rule == grammar::RuleId::BE_WIP_BASEPAIR ||
        deduction.rule == grammar::RuleId::BE_BASEPAIR_WIP;
    if (!recursive) return;
    const auto &parent = deduction.parent;
    const auto &split = deduction.split;
    const BeCoreIdentity core{split.k, split.l, parent.ip, parent.jp};
    ++recursive_sources_per_core_[core];
    const bool has_right_shell =
        deduction.rule == grammar::RuleId::BE_WIP_WIP ||
        deduction.rule == grammar::RuleId::BE_BASEPAIR_WIP;
    if (has_right_shell) {
      const BeRightTailIdentity right_tail{
          parent.j, split.k, split.l, parent.ip, parent.jp};
      ++right_shell_sources_per_tail_[right_tail];
    }
    if (deduction.rule != grammar::RuleId::BE_WIP_WIP) return;
    const BeRightTailIdentity right_tail{
        parent.j, split.k, split.l, parent.ip, parent.jp};
    const BeFrameIdentity frame{parent.i,
                                parent.j,
                                split.k,
                                split.l,
                                parent.ip,
                                parent.jp};
    ++sources_per_core_[core];
    ++sources_per_right_tail_[right_tail];
    right_tails_per_core_[core].insert(right_tail);
    cores_by_column_[split.l].insert(core);
    right_tails_by_column_[parent.j].insert(right_tail);
    frames_by_column_[parent.j].insert(frame);
  }

  BeFrameReuseProfile profile() const {
    BeFrameReuseProfile out;
    for (const auto &[core, count] : sources_per_core_) {
      (void)core;
      out.source_count += count;
      if (count > 1) ++out.reused_core_count;
      out.maximum_sources_per_core =
          std::max(out.maximum_sources_per_core, count);
    }
    out.core_count = sources_per_core_.size();
    for (const auto &[core, tails] : right_tails_per_core_) {
      (void)core;
      out.maximum_right_tails_per_core =
          std::max(out.maximum_right_tails_per_core, tails.size());
    }
    out.right_tail_count = sources_per_right_tail_.size();
    for (const auto &[tail, count] : sources_per_right_tail_) {
      (void)tail;
      if (count > 1) ++out.reused_right_tail_count;
      out.maximum_sources_per_right_tail =
          std::max(out.maximum_sources_per_right_tail, count);
    }
    profile_columns(cores_by_column_,
                    &out.core_column_count,
                    &out.multi_core_column_count,
                    &out.maximum_cores_per_column);
    profile_columns(right_tails_by_column_,
                    &out.right_tail_column_count,
                    &out.multi_right_tail_column_count,
                    &out.maximum_right_tails_per_column);
    // A frame contains its closing column, so these sets are disjoint.
    out.frame_count = profile_columns(frames_by_column_,
                                      &out.frame_column_count,
                                      &out.multi_frame_column_count,
                                      &out.maximum_frames_per_column);
    for (const auto &[core, count] : recursive_sources_per_core_) {
      (void)core;
      out.recursive_source_count += count;
      if (count > 1) ++out.reused_shared_core_count;
      out.maximum_recursive_sources_per_core =
          std::max(out.maximum_recursive_sources_per_core, count);
    }
    out.shared_core_count = recursive_sources_per_core_.size();
    for (const auto &[tail, count] : right_shell_sources_per_tail_) {
      (void)tail;
      out.right_shell_source_count += count;
      if (count > 1) ++out.reused_shared_right_tail_count;
      out.maximum_right_shell_sources_per_tail =
          std::max(out.maximum_right_shell_sources_per_tail, count);
    }
    out.shared_right_tail_count = right_shell_sources_per_tail_.size();
    return out;
  }

 private:
  template <class Identity>
  static std::size_t profile_columns(
      const std::map<cand_pos_t, std::set<Identity>> &columns,
      std::size_t *column_count,
      std::size_t *multi_column_count,
      std::size_t *maximum_per_column) {
    *column_count = columns.size();
    std::size_t total = 0;
    for (const auto &[j, identities] : columns) {
      (void)j;
      total += identities.size();
      if (identities.size() > 1) ++*multi_column_count;
      *maximum_per_column = std::max(*maximum_per_column, identities.size());
    }
    return total;
  }

  std::map<BeCoreIdentity, std::size_t> sources_per_core_;
  std::map<BeRightTailIdentity, std::size_t> sources_per_right_tail_;
  std::map<BeCoreIdentity, std::set<BeRightTailIdentity>>
      right_tails_per_core_;
  std::map<cand_pos_t, std::set<BeCoreIdentity>> cores_by_column_;
  std::map<cand_pos_t, std::set<BeRightTailIdentity>>
      right_tails_by_column_;
  std::map<cand_pos_t, std::set<BeFrameIdentity>> frames_by_column_;
  std::map<BeCoreIdentity, std::size_t> recursive_sources_per_core_;
  std::map<BeRightTailIdentity, std::size_t> right_shell_sources_per_tail_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_FRAME_REUSE_PROFILE_HH

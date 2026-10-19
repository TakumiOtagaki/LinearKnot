#ifndef CPARTY_SCFG2_BEAM_BE_LEFT_RIGHT_FRAME_PROFILE_HH
#define CPARTY_SCFG2_BEAM_BE_LEFT_RIGHT_FRAME_PROFILE_HH

#include "scfg2/beam/be_frame_cache.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <tuple>

namespace cparty::scfg2::beam {

enum class BeRightShellFrameMode : unsigned char {
  RightOnly,
  LeftRight,
};

struct BeRightShellFrameKey {
  BeRightTailKey tail{};
  cand_pos_t i = -1;
  BeRightShellFrameMode mode = BeRightShellFrameMode::RightOnly;
};

struct BeRightShellFrameKeyHash {
  std::size_t operator()(const BeRightShellFrameKey &key) const noexcept {
    std::size_t out = BeRightTailKeyHash{}(key.tail);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.i + 29);
    out = (out * 1315423911u) ^
          static_cast<std::size_t>(key.mode);
    return out;
  }
};

inline BeRightShellFrameKey be_right_shell_frame_key(
    const deductions::Deduction &deduction) {
  return BeRightShellFrameKey{
      BeRightTailKey{deduction.parent.j, be_core_key(deduction)},
      deduction.parent.i,
      deduction.rule == grammar::RuleId::BE_WIP_WIP
          ? BeRightShellFrameMode::LeftRight
          : BeRightShellFrameMode::RightOnly};
}

inline constexpr bool operator==(const BeRightShellFrameKey &lhs,
                                 const BeRightShellFrameKey &rhs) {
  return lhs.tail == rhs.tail && lhs.i == rhs.i && lhs.mode == rhs.mode;
}

inline auto be_right_shell_frame_tie_key(const BeRightShellFrameKey &key) {
  return std::make_tuple(key.i,
                         key.tail.j,
                         key.tail.core.k,
                         key.tail.core.l,
                         key.tail.core.ip,
                         key.tail.core.jp,
                         static_cast<unsigned char>(key.mode));
}

struct BeRightShellFrameProfile {
  std::size_t source_count = 0;
  std::size_t right_only_count = 0;
  std::size_t left_right_count = 0;
  std::size_t frame_count = 0;
  std::size_t tail_count = 0;
  std::size_t multi_frame_tail_count = 0;
  std::size_t maximum_frames_per_tail = 0;
  std::size_t column_count = 0;
  std::size_t multi_frame_column_count = 0;
  std::size_t maximum_frames_per_column = 0;
};

class BeRightShellFrameProfiler {
 public:
  template <class Chart>
  void observe(const deductions::Deduction &deduction,
               const Chart &chart,
               double local_weight,
               double direct_contribution) {
    if (!recursive_be_rule_has_right_wip(deduction.rule)) return;
    const BeRightShellFrameKey frame =
        be_right_shell_frame_key(deduction);
    const auto &tail = frame.tail;
    const auto mode = frame.mode;
    double factor =
        chart.get(deductions::ItemKey{grammar::NonTerminal::BE,
                                      tail.core.k,
                                      tail.core.l,
                                      tail.core.ip,
                                      tail.core.jp}) *
        chart.get(deductions::ItemKey{grammar::NonTerminal::WIP,
                                      tail.core.l + 1,
                                      tail.j - 1,
                                      -1,
                                      -1});
    if (mode == BeRightShellFrameMode::LeftRight) {
      factor *= chart.get(deduction.children[0].item);
      ++left_right_count_;
    } else {
      ++right_only_count_;
    }
    frames_per_tail_[tail].insert(frame);
    frames_by_column_[tail.j].insert(frame);
    const double reconstructed = local_weight * factor;
    const double tolerance =
        std::max(1e-14,
                 std::max(std::fabs(direct_contribution),
                          std::fabs(reconstructed)) *
                     1e-12);
    if (std::fabs(direct_contribution - reconstructed) <= tolerance) {
      ++factor_identity_count_;
    }
    ++source_count_;
  }

  BeRightShellFrameProfile profile() const {
    BeRightShellFrameProfile out;
    out.source_count = source_count_;
    out.right_only_count = right_only_count_;
    out.left_right_count = left_right_count_;
    out.tail_count = frames_per_tail_.size();
    for (const auto &[tail, frames] : frames_per_tail_) {
      (void)tail;
      if (frames.size() > 1) ++out.multi_frame_tail_count;
      out.maximum_frames_per_tail =
          std::max(out.maximum_frames_per_tail, frames.size());
    }
    out.column_count = frames_by_column_.size();
    for (const auto &[j, frames] : frames_by_column_) {
      (void)j;
      // The closing column is part of every frame's identity.
      out.frame_count += frames.size();
      if (frames.size() > 1) ++out.multi_frame_column_count;
      out.maximum_frames_per_column =
          std::max(out.maximum_frames_per_column, frames.size());
    }
    return out;
  }

  bool all_factor_identities_hold() const {
    return factor_identity_count_ == source_count_;
  }

 private:
  struct FrameLess {
    bool operator()(const BeRightShellFrameKey &lhs,
                    const BeRightShellFrameKey &rhs) const {
      return be_right_shell_frame_tie_key(lhs) <
             be_right_shell_frame_tie_key(rhs);
    }
  };
  struct TailLess {
    static auto be_right_tail_tie(const BeRightTailKey &key) {
      return std::make_tuple(key.j,
                             key.core.k,
                             key.core.l,
                             key.core.ip,
                             key.core.jp);
    }

    bool operator()(const BeRightTailKey &lhs,
                    const BeRightTailKey &rhs) const {
      return be_right_tail_tie(lhs) < be_right_tail_tie(rhs);
    }
  };

  std::map<BeRightTailKey,
           std::set<BeRightShellFrameKey, FrameLess>,
           TailLess>
      frames_per_tail_;
  std::map<cand_pos_t, std::set<BeRightShellFrameKey, FrameLess>>
      frames_by_column_;
  std::size_t source_count_ = 0;
  std::size_t right_only_count_ = 0;
  std::size_t left_right_count_ = 0;
  std::size_t factor_identity_count_ = 0;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_LEFT_RIGHT_FRAME_PROFILE_HH

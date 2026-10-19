#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_RIGHT_REUSE_PROFILE_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_RIGHT_REUSE_PROFILE_HH

#include "scfg2/beam/vp_case3_tail_cache.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <tuple>

namespace cparty::scfg2::beam {

struct VpCase3RightReuseProfile {
  std::size_t tail_count = 0;
  std::size_t right_view_count = 0;
  std::size_t nonempty_right_view_count = 0;
  std::size_t empty_right_view_count = 0;
  std::size_t reused_right_view_count = 0;
  std::size_t maximum_tails_per_right_view = 0;
  std::size_t column_count = 0;
  std::size_t multi_right_view_column_count = 0;
  std::size_t maximum_right_views_per_column = 0;
};

inline VpCase3RightReuseProfile profile_vp_case3_right_reuse(
    const ExactVpCase3TailCache &tails) {
  using RightIdentity = std::tuple<cand_pos_t, cand_pos_t>;
  std::map<RightIdentity, std::size_t> tails_per_right;
  std::map<cand_pos_t, std::set<cand_pos_t>> right_l_by_column;
  tails.for_each_tail_key([&](const VpCase3TailKey &key) {
    ++tails_per_right[RightIdentity{key.j, key.l}];
    right_l_by_column[key.j].insert(key.l);
  });

  VpCase3RightReuseProfile out;
  out.tail_count = tails.tail_count();
  out.right_view_count = tails_per_right.size();
  for (const auto &[identity, count] : tails_per_right) {
    const auto [j, l] = identity;
    if (l + 1 <= j - 1) {
      ++out.nonempty_right_view_count;
    } else {
      ++out.empty_right_view_count;
    }
    if (count > 1) ++out.reused_right_view_count;
    out.maximum_tails_per_right_view =
        std::max(out.maximum_tails_per_right_view, count);
  }
  out.column_count = right_l_by_column.size();
  for (const auto &[j, right_l] : right_l_by_column) {
    (void)j;
    if (right_l.size() > 1) ++out.multi_right_view_column_count;
    out.maximum_right_views_per_column =
        std::max(out.maximum_right_views_per_column, right_l.size());
  }
  return out;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_RIGHT_REUSE_PROFILE_HH

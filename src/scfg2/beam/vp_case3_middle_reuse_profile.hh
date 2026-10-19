#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_MIDDLE_REUSE_PROFILE_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_MIDDLE_REUSE_PROFILE_HH

#include "scfg2/beam/vp_case3_tail_cache.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <tuple>

namespace cparty::scfg2::beam {

struct VpCase3MiddleReuseProfile {
  std::size_t tail_count = 0;
  std::size_t middle_view_count = 0;
  std::size_t nonempty_middle_view_count = 0;
  std::size_t empty_middle_view_count = 0;
  std::size_t reused_middle_view_count = 0;
  std::size_t maximum_tails_per_middle_view = 0;
  std::size_t column_count = 0;
  std::size_t multi_middle_view_column_count = 0;
  std::size_t maximum_middle_views_per_column = 0;
};

inline VpCase3MiddleReuseProfile profile_vp_case3_middle_reuse(
    const ExactVpCase3TailCache &tails) {
  using MiddleIdentity = std::tuple<cand_pos_t, cand_pos_t>;
  std::map<MiddleIdentity, std::size_t> tails_per_middle;
  std::map<cand_pos_t, std::set<MiddleIdentity>> middle_by_column;
  tails.for_each_tail_key([&](const VpCase3TailKey &key) {
    const MiddleIdentity identity{key.q, key.k};
    ++tails_per_middle[identity];
    middle_by_column[key.j].insert(identity);
  });

  VpCase3MiddleReuseProfile out;
  out.tail_count = tails.tail_count();
  out.middle_view_count = tails_per_middle.size();
  for (const auto &[identity, count] : tails_per_middle) {
    const auto [q, k] = identity;
    if (q + 1 <= k - 1) {
      ++out.nonempty_middle_view_count;
    } else {
      ++out.empty_middle_view_count;
    }
    if (count > 1) ++out.reused_middle_view_count;
    out.maximum_tails_per_middle_view =
        std::max(out.maximum_tails_per_middle_view, count);
  }
  out.column_count = middle_by_column.size();
  for (const auto &[j, middle_views] : middle_by_column) {
    (void)j;
    if (middle_views.size() > 1) ++out.multi_middle_view_column_count;
    out.maximum_middle_views_per_column =
        std::max(out.maximum_middle_views_per_column, middle_views.size());
  }
  return out;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_MIDDLE_REUSE_PROFILE_HH

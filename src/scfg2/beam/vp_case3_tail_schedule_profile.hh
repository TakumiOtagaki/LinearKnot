#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_SCHEDULE_PROFILE_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_SCHEDULE_PROFILE_HH

#include "scfg2/beam/vp_case3_tail_cache.hh"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase3TailScheduleProfile {
  std::size_t tail_count = 0;
  std::size_t column_count = 0;
  std::size_t topological_tail_count = 0;
  std::size_t committable_column_count = 0;
  std::size_t noncommittable_column_count = 0;
  std::size_t missing_schedule_item_count = 0;
  std::ptrdiff_t maximum_column_overlap = 0;
};

struct VpCase3TailColumnReadiness {
  std::ptrdiff_t latest_tail_ready = -1;
  std::ptrdiff_t earliest_parent =
      std::numeric_limits<std::ptrdiff_t>::max();
};

inline VpCase3TailScheduleProfile profile_vp_case3_tail_schedule(
    const ExactVpCase3TailCache &cache,
    const std::vector<deductions::ItemKey> &schedule) {
  std::unordered_map<deductions::ItemKey,
                     std::ptrdiff_t,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      position;
  for (std::size_t idx = 0; idx < schedule.size(); ++idx) {
    position.emplace(schedule[idx], static_cast<std::ptrdiff_t>(idx));
  }

  VpCase3TailScheduleProfile profile;
  profile.tail_count = cache.tail_count();
  std::unordered_map<cand_pos_t, VpCase3TailColumnReadiness> columns;
  cache.for_each_tail_key([&](const VpCase3TailKey &key) {
    std::ptrdiff_t ready = -1;
    const auto observe_child = [&](cand_pos_t left, cand_pos_t right) {
      if (left > right) return true;
      const deductions::ItemKey child{
          grammar::NonTerminal::WI, left, right, -1, -1};
      const auto it = position.find(child);
      if (it == position.end()) return false;
      ready = std::max(ready, it->second);
      return true;
    };
    if (!observe_child(key.q + 1, key.k - 1)) {
      ++profile.missing_schedule_item_count;
    }
    if (!observe_child(key.l + 1, key.j - 1)) {
      ++profile.missing_schedule_item_count;
    }

    bool tail_is_topological = true;
    auto &column = columns[key.j];
    column.latest_tail_ready = std::max(column.latest_tail_ready, ready);
    for (const auto &attachment : cache.source_attachments(key)) {
      const auto parent_it = position.find(attachment.source_parent);
      if (parent_it == position.end()) {
        ++profile.missing_schedule_item_count;
        tail_is_topological = false;
        continue;
      }
      column.earliest_parent =
          std::min(column.earliest_parent, parent_it->second);
      if (ready >= parent_it->second) tail_is_topological = false;
    }
    if (tail_is_topological) ++profile.topological_tail_count;
  });

  profile.column_count = columns.size();
  for (const auto &[j, column] : columns) {
    (void)j;
    if (column.latest_tail_ready < column.earliest_parent) {
      ++profile.committable_column_count;
    } else {
      ++profile.noncommittable_column_count;
      profile.maximum_column_overlap = std::max(
          profile.maximum_column_overlap,
          column.latest_tail_ready - column.earliest_parent + 1);
    }
  }
  return profile;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_SCHEDULE_PROFILE_HH

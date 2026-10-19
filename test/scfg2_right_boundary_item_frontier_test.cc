#include "scfg2/beam/right_boundary_item_frontier.hh"
#include "scfg2/beam/wmb_item_selection.hh"

#include <cassert>
#include <cstdlib>
#include <vector>

int main() {
  using cparty::scfg2::beam::RightBoundaryTopWidthItemFrontier;
  using cparty::scfg2::deductions::ItemKey;
  using cparty::scfg2::grammar::NonTerminal;

  const ItemKey i12{NonTerminal::WI, 1, 2, -1, -1};
  const ItemKey i22{NonTerminal::WI, 2, 2, -1, -1};
  const ItemKey i32{NonTerminal::WI, 3, 2, -1, -1};
  const ItemKey i13{NonTerminal::WI, 1, 3, -1, -1};

  RightBoundaryTopWidthItemFrontier top_two(2);
  top_two.observe(i22, 5.0);
  top_two.observe(i12, 5.0);
  top_two.observe(i32, 1.0);
  top_two.observe(i13, 9.0);
  std::vector<ItemKey> retained;
  top_two.flush_before(3, [&](const ItemKey &item) {
    retained.push_back(item);
  });
  assert(retained.size() == 2);
  assert(retained[0].i == 1);
  assert(retained[1].i == 2);
  assert(top_two.candidate_count() == 4);
  assert(top_two.selected_count() == 2);
  assert(top_two.pruned_count() == 1);
  assert(top_two.respects_width());

  top_two.flush_before(4, [&](const ItemKey &item) {
    retained.push_back(item);
  });
  assert(retained.size() == 3);
  assert(retained.back().i == 1 && retained.back().j == 3);
  assert(top_two.selected_count() == 3);
  assert(top_two.pruned_count() == 1);

  auto keep_all = RightBoundaryTopWidthItemFrontier::keep_all();
  keep_all.observe(i22, 1.0);
  keep_all.observe(i12, 2.0);
  std::size_t keep_all_count = 0;
  keep_all.flush_before(3, [&](const ItemKey &) { ++keep_all_count; });
  assert(keep_all_count == 2);
  assert(keep_all.selected_count() == keep_all.candidate_count());
  assert(keep_all.pruned_count() == 0);

  RightBoundaryTopWidthItemFrontier partitioned(1);
  partitioned.observe(i12, 3.0);
  partitioned.observe(i22, 2.0);
  partitioned.observe(i32, 1.0);
  std::vector<ItemKey> partition_retained;
  std::vector<ItemKey> partition_pruned;
  partitioned.flush_column_partition(
      2,
      [&](const ItemKey &item) { partition_retained.push_back(item); },
      [&](const ItemKey &item) { partition_pruned.push_back(item); });
  assert(partition_retained.size() == 1);
  assert(partition_retained.front().i == 1);
  assert(partition_pruned.size() == 2);
  assert(partitioned.selected_count() == 1);
  assert(partitioned.pruned_count() == 2);

  // A semantic chart item can be completed more than once as independent
  // deduction phases contribute to it.  The beam is over chart items, not
  // completion notifications: a lower-priority duplicate must neither spend
  // another slot nor prune the retained occurrence of the same item.
  RightBoundaryTopWidthItemFrontier duplicate_completion(1);
  duplicate_completion.observe(i12, 5.0);
  duplicate_completion.observe(i12, 4.0);
  duplicate_completion.observe(i22, 3.0);
  std::vector<ItemKey> duplicate_retained;
  std::vector<ItemKey> duplicate_pruned;
  duplicate_completion.flush_column_partition(
      2,
      [&](const ItemKey &item) { duplicate_retained.push_back(item); },
      [&](const ItemKey &item) { duplicate_pruned.push_back(item); });
  assert(duplicate_retained.size() == 1);
  assert(duplicate_retained.front().i == 1);
  assert(duplicate_pruned.size() == 1);
  assert(duplicate_pruned.front().i == 2);
  assert(duplicate_completion.candidate_count() == 2);
  assert(duplicate_completion.selected_count() == 1);
  assert(duplicate_completion.pruned_count() == 1);

  const ItemKey wmb12{NonTerminal::WMB, 1, 2, -1, -1};
  const ItemKey wmb22{NonTerminal::WMB, 2, 2, -1, -1};
  const ItemKey wmb32{NonTerminal::WMB, 3, 2, -1, -1};
  cparty::scfg2::beam::TopWidthWmbItemSelection wmb_views(1);
  wmb_views.observe_completed(wmb12, 5.0, 1.0);
  wmb_views.observe_completed(wmb22, 4.0, 9.0);
  wmb_views.observe_completed(wmb32, 3.0, 2.0);
  struct ZeroedChart {
    std::vector<ItemKey> items;
    void set(const ItemKey &item, double) { items.push_back(item); }
  } zeroed;
  wmb_views.commit_column(2, zeroed, 0.0);
  assert(wmb_views.internal_starts(2).size() == 1);
  assert(wmb_views.internal_starts(2).front() == 1);
  assert(wmb_views.root_suffix_starts(2).size() == 1);
  assert(wmb_views.root_suffix_starts(2).front() == 2);
  assert(wmb_views.candidate_count() == 3);
  assert(wmb_views.selected_count() == 2);
  assert(zeroed.items.size() == 1);
  assert(zeroed.items.front().i == 3);
  return EXIT_SUCCESS;
}

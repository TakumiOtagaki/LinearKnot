#include "scfg2/beam/right_boundary_item_frontier.hh"
#include "scfg2/beam/wmb_item_selection.hh"
#include "scfg2/engine/chart.hh"
#include <array>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <cassert>

namespace beam = cparty::scfg2::beam;
namespace sc = cparty::scfg2;
using Item = sc::deductions::ItemKey;
using NT = sc::grammar::NonTerminal;
using Actual = beam::RightBoundaryTopWidthItemFrontier;
using Chart = sc::engine::Chart<double>;
using Reference = Actual;

void existing_contracts() {
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
  return;
}

std::size_t checks = 0, changed_comparison_counts = 0;
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
bool same(const Item &a, const Item &b) { return sc::deductions::ItemKeyEq{}(a, b); }
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
void counts(const Actual &actual, const Reference &reference) {
  require(actual.candidate_count() == reference.candidate_count() &&
          actual.selected_count() == reference.selected_count() &&
          actual.pruned_count() == reference.pruned_count() &&
          actual.respects_width() == reference.respects_width(), "frontier counters differ");
  ++checks;
}
void keys(const std::vector<Item> &a, const std::vector<Item> &b) {
  require(a.size() == b.size(), "callback item count differs");
  for (std::size_t p = 0; p < a.size(); ++p)
    require(same(a[p], b[p]), "callback item/order differs");
}
void sort_keys(std::vector<Item> &items) {
  std::sort(items.begin(), items.end(), [](const auto &a, const auto &b) {
    return beam::right_boundary_item_tie_key(a) < beam::right_boundary_item_tie_key(b);
  });
}
void charts(const Chart &a, const Chart &b) {
  const auto &left = a.values();
  const auto &right = b.values();
  require(left.size() == right.size() && left.bucket_count() == right.bucket_count(),
          "canonical chart size/buckets differ");
  auto x = left.begin(), y = right.begin();
  for (; x != left.end(); ++x, ++y)
    require(same(x->first, y->first) && same(x->second, y->second),
            "canonical chart iteration/value bits differ after prune");
  for (std::size_t bucket = 0; bucket < left.bucket_count(); ++bucket) {
    require(left.bucket_size(bucket) == right.bucket_size(bucket), "chart bucket size differs");
    auto p = left.begin(bucket), q = right.begin(bucket);
    for (; p != left.end(bucket); ++p, ++q)
      require(same(p->first, q->first), "chart bucket chain differs");
  }
}
constexpr double zero = -std::numeric_limits<double>::infinity();
bool discard(const double &value) { return value == zero; }

// Callback observations include frontier counters before their end-of-flush
// update and the chart's completed state. The direct contract permits only
// pruned ordering to differ, so intermediate erase order is not compared.
void flush(Actual &actual, Reference &reference, cand_pos_t j, Chart &a, Chart &b,
           bool exact_pruned_order) {
  std::vector<Item> ar, ap, br, bp;
  std::vector<std::array<std::size_t, 3>> ac, bc;
  actual.flush_column_partition(j,
      [&](const Item &item) {
        ar.push_back(item);
        ac.push_back({actual.candidate_count(), actual.selected_count(), actual.pruned_count()});
      }, [&](const Item &item) {
        ap.push_back(item);
        ac.push_back({actual.candidate_count(), actual.selected_count(), actual.pruned_count()});
        a.set(item, zero);
      });
  reference.flush_column_partition(j,
      [&](const Item &item) {
        br.push_back(item);
        bc.push_back({reference.candidate_count(), reference.selected_count(), reference.pruned_count()});
      }, [&](const Item &item) {
        bp.push_back(item);
        bc.push_back({reference.candidate_count(), reference.selected_count(), reference.pruned_count()});
        b.set(item, zero);
      });
  keys(ar, br);
  require(ac == bc, "callback-time frontier observations differ");
  if (!exact_pruned_order) { sort_keys(ap); sort_keys(bp); }
  keys(ap, bp);
  counts(actual, reference);
  charts(a, b);
  changed_comparison_counts += actual.selection_comparisons() != reference.selection_comparisons();
}

void observe(Actual &actual, Reference &reference, Chart &a, Chart &b,
             const Item &item, double priority, double stored) {
  a.set(item, stored); b.set(item, stored);
  actual.observe(item, priority); reference.observe(item, priority);
  counts(actual, reference);
}

void trace(std::size_t width, bool forced, bool direct, unsigned seed) {
  Actual actual(width);
  Reference reference(width);
  constexpr cand_pos_t n = 128;
  Chart a(zero, discard, 0, n), b(zero, discard, 0, n);
  if (direct) require(actual.enable_numeric_canonical_storage(NT::WI, n), "opt-in failed");
  if (forced) {
    auto items = std::make_shared<beam::ItemKeySet>();
    for (cand_pos_t j = 8; j <= n; j += 8)
      for (cand_pos_t i = 1; i <= j; i += 11) items->insert({NT::WI, i, j, -1, -1});
    items->insert({NT::BE, 1, n, 2, n - 1});
    actual.set_forced(items); reference.set_forced(items);
  }
  std::mt19937 random(seed);
  const double special[]{-0.0, +0.0, zero, std::numeric_limits<double>::infinity()};
  for (cand_pos_t j = 8; j <= n; j += 8) {
    for (int step = 0; step < 3 * j; ++step) {
      const cand_pos_t i = 1 + random() % j;
      const double priority = step % 9 == 0 ? special[(step / 9) % 4]
          : double(static_cast<int>(random() % 17) - 8);
      observe(actual, reference, a, b, {NT::WI, i, j, -1, -1}, priority, double(i + j));
    }
    flush(actual, reference, j - 1, a, b, !direct); // Other column is absent.

    flush(actual, reference, j, a, b, !direct);
    if (direct) require(actual.numeric_canonical_storage_enabled(), "ordinary direct trace fell back");
  }
}

void fallback(int kind) {
  Actual actual(7); Reference reference(7);
  Chart a(zero, discard, 0, 64), b(zero, discard, 0, 64);
  require(actual.enable_numeric_canonical_storage(NT::WI, 64), "fallback fixture opt-in failed");
  for (cand_pos_t i = 1; i <= 31; ++i) {
    const Item key{NT::WI, i, 40, -1, -1};
    observe(actual, reference, a, b, key, double(i % 3), double(i));
    observe(actual, reference, a, b, key, double(i % 3 + 1), double(i));
  }
  Item extra{NT::WI, 12, 40, -1, -1};
  double priority = 3.0;
  switch (kind) {
    case 0: extra.nonterminal = NT::BE; extra.ip = 14; extra.jp = 25; break;
    case 1: extra.ip = 17; break;
    case 2: extra.i = -1; break;
    case 3: extra.j = 65; break;
    case 4: extra.j = 41; break;
    case 5: priority = std::numeric_limits<double>::quiet_NaN(); break; // Existing finite key.
    case 6: extra.nonterminal = static_cast<NT>(255); break;
  }
  observe(actual, reference, a, b, extra, priority, 9.0);
  require(!actual.numeric_canonical_storage_enabled(), "unsupported input failed to restore old map");
  flush(actual, reference, 40, a, b, true);
  flush(actual, reference, 41, a, b, true);
  flush(actual, reference, 65, a, b, true);
  require(actual.selection_comparisons() == reference.selection_comparisons(),
          "fallback did not recover old selection comparison protocol");
}


void moved_from_reuse() {
  for (bool assignment : {false, true}) {
    Actual source(2);
    Reference source_reference(2);
    Chart a(zero, discard, 0, 64), b(zero, discard, 0, 64);
    require(source.enable_numeric_canonical_storage(NT::WI, 64), "move fixture opt-in failed");
    for (cand_pos_t i = 1; i <= 9; ++i)
      observe(source, source_reference, a, b, {NT::WI, i, 24, -1, -1}, double(i), double(i));
    Actual moved(1);
    Reference moved_reference(1);
    if (assignment) {
      moved = std::move(source);
      moved_reference = std::move(source_reference);
    } else {
      Actual constructed(std::move(source));
      Reference reference_constructed(std::move(source_reference));
      moved = std::move(constructed);
      moved_reference = std::move(reference_constructed);
    }
    flush(moved, moved_reference, 24, a, b, false);
    // Default moved-from containers are reusable. The direct source detects
    // its moved-away index and safely rejoins the original map path.
    Chart reused(zero, discard, 0, 64), reused_reference(zero, discard, 0, 64);
    observe(source, source_reference, reused, reused_reference,
            {NT::WI, 3, 12, -1, -1}, 4.0, 3.0);
    flush(source, source_reference, 12, reused, reused_reference, true);
  }
}

void before_and_domains() {
  Actual actual(2); Reference reference(2);
  require(!actual.enable_numeric_canonical_storage(NT::BE, 64), "BE enabled primary storage");
  require(!actual.enable_numeric_canonical_storage(NT::WI, -1), "negative extent enabled storage");
  require(actual.enable_numeric_canonical_storage(NT::WI, 64), "flush-before opt-in failed");
  for (cand_pos_t i = 1; i <= 9; ++i) {
    const Item item{NT::WI, i, 12, -1, -1};
    actual.observe(item, double(i % 2)); reference.observe(item, double(i % 2));
  }
  std::vector<Item> a, b;
  actual.flush_before(12, [&](const Item &item) { a.push_back(item); });
  reference.flush_before(12, [&](const Item &item) { b.push_back(item); });
  keys(a, b); counts(actual, reference);
  actual.flush_before(13, [&](const Item &item) { a.push_back(item); });
  reference.flush_before(13, [&](const Item &item) { b.push_back(item); });
  keys(a, b); counts(actual, reference);

  Actual pending(2);
  pending.observe({NT::WI, 1, 2, -1, -1}, 0.0);
  require(!pending.enable_numeric_canonical_storage(NT::WI, 64), "pending map allowed destructive opt-in");
}

void copy_different_domains() {
  Actual source(7), target(2);
  Reference reference(7);
  Chart a(zero, discard, 0, 128), b(zero, discard, 0, 128);
  require(source.enable_numeric_canonical_storage(NT::V, 128), "copy source opt-in failed");
  require(target.enable_numeric_canonical_storage(NT::WI, 64), "copy target opt-in failed");
  target.observe({NT::WI, 1, 40, -1, -1}, 1.0);
  for (cand_pos_t i = 1; i <= 23; ++i)
    observe(source, reference, a, b, {NT::V, i, 80, -1, -1}, double(i % 3), double(i));
  target = source;
  target = target;
  auto copy = target;
  Chart copied(a), copied_reference(b);
  flush(copy, reference, 80, copied, copied_reference, false);
  require(source.numeric_canonical_storage_enabled() && target.numeric_canonical_storage_enabled(),
          "copy changed source storage mode");
  Reference target_reference(7);
  for (cand_pos_t i = 1; i <= 23; ++i) target_reference.observe({NT::V, i, 80, -1, -1}, double(i % 3));
  flush(target, target_reference, 80, a, b, false);
}

void wmb_union_storage() {
  for (const std::size_t width : {std::size_t(0), std::size_t(1), std::size_t(2)}) {
    beam::TopWidthWmbItemSelection reference(width), actual(width);
    actual.enable_numeric_canonical_storage(16);
    Chart a(zero, discard, 0, 16), b(zero, discard, 0, 16);
    for (const cand_pos_t j : {cand_pos_t(12), cand_pos_t(16)}) {
      for (cand_pos_t i = 1; i <= 6; ++i) {
        const Item item{NT::WMB, i, j, -1, -1};
        a.set(item, double(i)); b.set(item, double(i));
        reference.observe_ranked(item, true, double(7 - i), double(i));
        actual.observe_ranked(item, true, double(7 - i), double(i));
      }
      actual.commit_column(j, a, zero);
      reference.commit_column(j, b, zero);
      require(actual.internal_starts(j) == reference.internal_starts(j) &&
              actual.root_suffix_starts(j) == reference.root_suffix_starts(j),
              "WMB two-view retained callback order differs");
      require(actual.candidate_count() == reference.candidate_count() &&
              actual.selected_count() == reference.selected_count(), "WMB two-view counts differ");
      for (cand_pos_t i = 1; i <= 6; ++i)
        require(actual.root_suffix_contains(i, j) == reference.root_suffix_contains(i, j),
                "WMB root membership differs");
      charts(a, b);
    }
  }
}

int main() {
  try {
    existing_contracts();
    for (const auto width : {std::size_t(0), std::size_t(1), std::size_t(7),
                             std::size_t(100), std::numeric_limits<std::size_t>::max()})
      for (const bool forced : {false, true}) trace(width, forced, true, 20261021);
    for (int kind = 0; kind < 7; ++kind) fallback(kind);
    before_and_domains();
    copy_different_domains();
    moved_from_reuse();
    wmb_union_storage();
    std::cout << "scfg2_right_boundary_item_frontier=ok checks=" << checks
              << " comparison_count_changes=" << changed_comparison_counts << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}

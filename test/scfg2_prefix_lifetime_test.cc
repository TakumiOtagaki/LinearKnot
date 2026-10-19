#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case12_prefix_interval.hh"
#include "scfg2/beam/vp_case3_nonleft_prefix_interval.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace cparty::scfg2;
using C12 = beam::VpCase12PrefixRecord;
using C3 = beam::VpCase3NonleftPrefixRecord;
std::size_t checks = 0, restored_case2 = 0, restored_case3 = 0;
std::size_t largest_pending = 0, largest_finite_pending = 0, max_frames = 0;

void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
template <class Record> bool before(const Record &a, const Record &b) {
  return a.priority != b.priority ? a.priority > b.priority : a.geometry < b.geometry;
}
template <class Record>
std::vector<Record> top(std::vector<Record> records, std::size_t width) {
  std::sort(records.begin(), records.end(), before<Record>);
  if (records.size() > width) records.resize(width);
  return records;
}
template <class Record>
void same(const std::vector<Record> &a, const std::vector<Record> &b) {
  require(a.size() == b.size(), "prefix support differs from interval-stabbing oracle");
  for (std::size_t i = 0; i < a.size(); ++i)
    require(!(a[i].geometry < b[i].geometry) && !(b[i].geometry < a[i].geometry) &&
                a[i].priority == b[i].priority,
            "prefix identity, tie order, or first-registration value changed");
  ++checks;
}
template <class Record>
void count_restored(const std::vector<Record> &active,
                    std::set<decltype(Record{}.geometry)> &ever,
                    std::set<decltype(Record{}.geometry)> &previous,
                    std::size_t &restored) {
  std::set<decltype(Record{}.geometry)> current;
  for (const auto &record : active) {
    if (!previous.count(record.geometry) && ever.count(record.geometry)) ++restored;
    ever.insert(record.geometry);
    current.insert(record.geometry);
  }
  previous = std::move(current);
}

// This oracle keeps every first registration, selects by interval membership,
// and has no pending agenda, frame stack, or expiry/restore simulation.
template <class S>
void check(std::string g, std::size_t width, bool ties, bool bounded_notifications = false) {
  std::string seq(g.size(), 'A');
  for (std::size_t i = 0; i < g.size(); ++i)
    seq[i] = g[i] == '(' ? 'G' : g[i] == ')' ? 'C' : (i % 2 ? 'U' : 'A');
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(g.size());
  const auto left = beam::VpCase12DirectNonemptyLeftQuery::build(n, ctx);
  const auto empty = beam::VpCase12DirectEmptyLeftQuery::build(n, ctx, width);
  const auto middle = beam::VpCase3NonleftChildEventIndex::build(n, ctx, 30);
  beam::VpCase12PrefixIntervalAgenda a12(width);
  beam::VpCase3NonleftPrefixIntervalAgenda a3(width);
  std::map<beam::VpCase12PrefixGeometry, double> all12;
  std::map<beam::VpCase3NonleftPrefixGeometry, double> all3;
  std::set<beam::VpCase12PrefixGeometry> ever2, previous2;
  std::set<beam::VpCase3NonleftPrefixGeometry> ever3, previous3;
  std::map<beam::VpCase12PrefixGeometry, deductions::ItemKey> trigger12;
  std::map<beam::VpCase3NonleftPrefixGeometry, deductions::ItemKey> trigger3;
  std::vector<deductions::ItemKey> old_wi;
  std::size_t duplicate12 = 0, duplicate3 = 0;

  const auto storage = [&] {
    for (const auto s : {a12.storage_snapshot(), a3.storage_snapshot()}) {
      largest_pending = std::max(largest_pending, s.largest_pending_bucket);
      if (width != std::numeric_limits<std::size_t>::max())
        largest_finite_pending = std::max(largest_finite_pending, s.largest_pending_bucket);
      max_frames = std::max(max_frames, s.saved_frames);
      require(width == std::numeric_limits<std::size_t>::max() ||
                  s.largest_pending_bucket <= width,
              "pending activation group exceeds its beam width");
      require(s.dedup_keys <= old_wi.size() + g.size() + 1,
              "dedup storage exceeds WI triggers plus activation columns");
      if (bounded_notifications) {
        const auto n = g.size() + 2;
        require(s.pending_records <= 2 * n * width &&
                    s.dedup_keys <= n * (width + 1) &&
                    s.active_records + s.saved_frame_records <= 2 * n * width,
                "agenda live storage exceeds the retained-WI O(nb) budget");
      }
    }
  };
  for (cand_pos_t j = 1; j <= n + 1; ++j) {
    if (j <= n) {
      for (const auto &p : empty.query_activation(j, ctx))
        all12.emplace(p, S::priority(S::one()));
      for (const auto &p : middle.query_empty_middle_prefixes(j, ctx, width))
        all3.emplace(p, S::priority(S::one()));
    }
    std::vector<cand_pos_t> right;
    for (cand_pos_t i = 1; i < j; ++i) right.push_back(i);
    a12.template prepare_column<S>(j, right, empty, ctx);
    a3.template prepare_column<S>(j, middle, ctx);

    std::map<cand_pos_t, std::vector<C12>> groups;
    std::vector<C12> expected1, expected2;
    std::vector<C3> expected3;
    for (const auto &[p, priority] : all12) {
      if (p.target_begin > j || j > p.target_end) continue;
      if (p.rule == grammar::RuleId::VP_WI_CASE1) groups[p.q].push_back({p, priority});
      else expected2.push_back({p, priority});
    }
    for (auto &[q, records] : groups) {
      (void)q;
      auto kept = top(std::move(records), width);
      expected1.insert(expected1.end(), kept.begin(), kept.end());
    }
    for (const auto &[p, priority] : all3)
      if (p.k <= j && j <= p.target_close) expected3.push_back({p, priority});
    expected2 = top(std::move(expected2), width);
    expected3 = top(std::move(expected3), width);
    same(a12.active_case1_prefixes(), expected1);
    same(a12.active_case2_prefixes(), expected2);
    same(a3.active_prefixes(), expected3);
    count_restored(expected2, ever2, previous2, restored_case2);
    count_restored(expected3, ever3, previous3, restored_case3);
    storage();
    if (j > n) break;

    // A second preparation must not register the empty-child event again.
    if (j % 3 == 0) {
      duplicate12 += empty.query_activation(j, ctx).size();
      duplicate3 += middle.query_empty_middle_prefixes(j, ctx, width).size();
      a12.template prepare_column<S>(j, right, empty, ctx);
      a3.template prepare_column<S>(j, middle, ctx);
      same(a12.active_case1_prefixes(), expected1);
      same(a12.active_case2_prefixes(), expected2);
      same(a3.active_prefixes(), expected3);
    }

    for (cand_pos_t i = 1; i <= j; ++i) {
      if (bounded_notifications && static_cast<std::size_t>(i) > width) break;
      const deductions::ItemKey wi{grammar::NonTerminal::WI, i, j, -1, -1};
      const double raw = ties ? 1. : 2. + j + i % 3;
      const auto value = S::from_weight(raw);
      a12.template observe_left_wi<S>(wi, S::zero(), left, ctx);
      a3.template observe_middle_wi<S>(wi, S::zero(), middle, ctx);
      const auto p12 = left.query_completed_left(wi, ctx, width);
      const auto p3 = middle.query_middle_prefixes(wi, ctx, width);
      for (const auto &p : p12) {
        require(p.target_begin > j, "CASE1/2 readiness contract failed");
        if (p.rule == grammar::RuleId::VP_WI_CASE2)
          require(p.target_begin == j + 1, "CASE2 does not activate next column");
        require(trigger12.emplace(p, wi).second, "CASE1/2 geometry has multiple WI triggers");
        require(all12.emplace(p, S::priority(value)).second, "empty/nonempty CASE1/2 overlap");
      }
      for (const auto &p : p3) {
        require(p.k == j + 1, "CASE3 does not activate next column");
        require(trigger3.emplace(p, wi).second, "CASE3 geometry has multiple WI triggers");
        require(all3.emplace(p, S::priority(value)).second, "empty/nonempty CASE3 overlap");
      }
      a12.template observe_left_wi<S>(wi, value, left, ctx);
      a3.template observe_middle_wi<S>(wi, value, middle, ctx);
      // Re-notification with a larger value must preserve the first value.
      a12.template observe_left_wi<S>(wi, S::from_weight(raw * 100.), left, ctx);
      a3.template observe_middle_wi<S>(wi, S::from_weight(raw * 100.), middle, ctx);
      duplicate12 += p12.size(); duplicate3 += p3.size();
      if (!p12.empty() || !p3.empty()) old_wi.push_back(wi);
    }
    // Include old notifications whose prefixes have activated or expired.
    if (!old_wi.empty() && j % 5 == 0) {
      const auto wi = old_wi.front();
      a12.template observe_left_wi<S>(wi, S::from_weight(10000.), left, ctx);
      a3.template observe_middle_wi<S>(wi, S::from_weight(10000.), middle, ctx);
      duplicate12 += left.query_completed_left(wi, ctx, width).size();
      duplicate3 += middle.query_middle_prefixes(wi, ctx, width).size();
    }
    storage();
  }
  require(a12.stats().duplicate_prefixes == duplicate12 &&
              a3.stats().duplicate_prefixes == duplicate3,
          "duplicate registration accounting changed");
  if (width != std::numeric_limits<std::size_t>::max()) {
    const auto budget = [&](std::size_t candidates) {
      return 32 * (candidates + 1) * static_cast<std::size_t>(std::ceil(std::log2(width + 2)));
    };
    require(a12.pending_comparisons() <= budget(a12.stats().prefix_candidates) &&
                a3.pending_comparisons() <= budget(a3.stats().prefix_candidates),
            "pending heap work exceeds logarithmic comparison budget");
  }
  require(a12.storage_snapshot().pending_records == 0 &&
              a3.storage_snapshot().pending_records == 0 &&
              a12.storage_snapshot().saved_frames == 0 &&
              a3.storage_snapshot().saved_frames == 0,
          "agenda has pending prefixes or frames after final expiry");
}
} // namespace

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    const std::vector<std::string> inputs{
      "......", ".(....)..(....)..",
      ".(....)..(..(....)..(....)..(....)..)..",
      ".(..(....)..(..(....)..(....)..)..(....)..)..",
      "....(...)((....)..)((.(.((...))))(.(....)((...))))"};
    for (const auto &g : inputs)
      for (auto width : {std::size_t(1), std::size_t(2), std::size_t(5),
                         std::numeric_limits<std::size_t>::max()})
        for (bool ties : {false, true}) {
          check<semirings::SumProductSemiring>(g, width, ties);
          check<semirings::LogMaxSemiring>(g, width, ties);
          if (width != std::numeric_limits<std::size_t>::max()) {
            check<semirings::SumProductSemiring>(g, width, ties, true);
            check<semirings::LogMaxSemiring>(g, width, ties, true);
          }
        }
    require(restored_case2 > 0 && restored_case3 > 0 && max_frames >= 2,
            "nested restoration witness is vacuous");
    std::cout << "scfg2_prefix_lifetime=ok checks=" << checks
              << " restored_case2=" << restored_case2 << " restored_case3=" << restored_case3
              << " max_frames=" << max_frames << " largest_pending=" << largest_pending
              << " largest_finite_pending=" << largest_finite_pending << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

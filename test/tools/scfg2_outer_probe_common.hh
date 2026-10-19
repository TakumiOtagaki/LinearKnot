// Observation and independent witnesses for the shared outer implementation.
#ifndef SCFG2_OUTER_PROBE_COMMON_HH
#define SCFG2_OUTER_PROBE_COMMON_HH

#include "scfg2/beam/outer_frontier.hh"

namespace sc = cparty::scfg2;
using Item = sc::deductions::ItemKey;
using NT = sc::grammar::NonTerminal;
using Sum = sc::semirings::LogSumExpSemiring;
using Max = sc::semirings::LogMaxSemiring;
using Chart = sc::engine::Chart<double>;
using Clock = std::chrono::steady_clock;
constexpr double zero = -std::numeric_limits<double>::infinity();
using ItemMap = std::unordered_map<Item, double, sc::deductions::ItemKeyHash,
                                 sc::deductions::ItemKeyEq>;
inline constexpr std::size_t outer_probe_families =
    static_cast<std::size_t>(NT::BE) + 1;
using ItemSet = std::unordered_set<Item, sc::deductions::ItemKeyHash,
                                  sc::deductions::ItemKeyEq>;

using sc::beam::outer_detail::Outer;
using sc::beam::outer_detail::MaxPool;
using sc::beam::outer_detail::less_key;
using sc::beam::outer_detail::seconds;
using sc::beam::outer_detail::better_link;
using sc::beam::outer_detail::update_outer;

// Independent batch reference checks every prefix of a deterministic stream,
// including eviction/re-entry, duplicates, ties and mandatory items.
inline void check_max_pool() {
  for (bool prefilter : {false, true}) for (std::size_t width : {0, 1, 3, 7}) {
    MaxPool pool;
    MaxPool::Entries all;
    ItemMap previous_expected;
    std::size_t expected_evicted = 0;
    std::uint32_t rng = 17423;
    for (int step = 0; step < 200; ++step) {
      rng = rng * 1664525u + 1013904223u;
      const int id = 1 + (rng % 23);
      const Item key{NT::V, id, 50, -1, -1};
      rng = rng * 1664525u + 1013904223u;
      const double value = double(rng % 17) - 8;
      const bool mandatory = id % 5 == 0;
      const Outer::Link link{{step % 2 ? NT::W : NT::VM,
                              1 + int((rng >> 8) % 13), 60, -1, -1}, value};
      std::size_t previous_novel = 0;
      for (const auto &[item, score] : previous_expected) {
        (void)score;
        if (item.i % 5 != 0) ++previous_novel;
      }
      if (!mandatory && !previous_expected.count(key) && previous_novel == width)
        ++expected_evicted;
      const auto old = all.find(key);
      if (old == all.end()) all.emplace(key, MaxPool::Entry{value, link, mandatory});
      else if (value > old->second.value ||
               (value == old->second.value && better_link(link, old->second.link)))
        old->second = {value, link, mandatory};
      pool.offer(key, value, link, mandatory, width, prefilter);
      std::vector<std::pair<Item, double>> sorted;
      ItemMap expected;
      for (const auto &[item, entry] : all) {
        if (entry.mandatory) expected.emplace(item, entry.value);
        else sorted.emplace_back(item, entry.value);
      }
      std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
        if (a.second != b.second) return a.second > b.second;
        return less_key(a.first, b.first);
      });
      for (std::size_t k = 0; k < std::min(width, sorted.size()); ++k)
        expected.emplace(sorted[k]);
      if (pool.entries.size() != expected.size()) throw std::runtime_error("max pool support mismatch");
      for (const auto &[item, score] : expected) {
        const auto &entry = pool.entries.at(item);
        const auto &want = all.at(item);
        if (entry.value != score || entry.mandatory != want.mandatory ||
            !sc::deductions::ItemKeyEq{}(entry.link.parent, want.link.parent) ||
            entry.link.edge != want.link.edge)
          throw std::runtime_error("max pool value/link/ownership mismatch");
      }
      if (pool.evicted != expected_evicted || pool.novel.size() != std::min(width, sorted.size()))
        throw std::runtime_error("max pool eviction/rank count mismatch");
      std::size_t rank = 0;
      for (const auto &entry : pool.novel) {
        if (!sc::deductions::ItemKeyEq{}(entry.key, sorted[rank].first) ||
            entry.value != sorted[rank].second)
          throw std::runtime_error("max pool ordered rank mismatch");
        ++rank;
      }
      previous_expected = std::move(expected);
      // Node reuse must not introduce stored iterators with copy/move aliases.
      if (step % 31 == 0) {
        auto copied = pool;
        MaxPool assigned;
        assigned = copied;
        auto moved = std::move(assigned);
        pool = std::move(moved);
      }
    }
  }
}

#endif

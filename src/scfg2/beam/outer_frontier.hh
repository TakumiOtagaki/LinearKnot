// Shared data structures for bounded completion-context ranking.
#ifndef CPARTY_SCFG2_BEAM_OUTER_FRONTIER_HH
#define CPARTY_SCFG2_BEAM_OUTER_FRONTIER_HH

#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam::outer_detail {

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

inline double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
inline bool less_key(const Item &a, const Item &b) {
  return std::tie(a.i, a.j, a.ip, a.jp) < std::tie(b.i, b.j, b.ip, b.jp);
}
struct Outer {
  struct Link { Item parent; double edge; };
  using Links = std::unordered_map<Item, Link, sc::deductions::ItemKeyHash,
                                   sc::deductions::ItemKeyEq>;
  Chart values{zero};
  Links links;
  std::size_t visited = 0, generated = 0, eligible = 0, selected = 0;
  std::size_t max_parent_candidates = 0;
  std::size_t peak_pending = 0, peak_guide_states = 0, evicted = 0;
  std::array<std::size_t, sc::grammar::all_rule_ids().size()> generated_rules{};
  std::array<std::size_t, sc::grammar::all_rule_ids().size()> eligible_rules{};
  double elapsed = 0;
};

inline bool better_link(const Outer::Link &a, const Outer::Link &b) {
  if (a.parent.nonterminal != b.parent.nonterminal)
    return a.parent.nonterminal < b.parent.nonterminal;
  return less_key(a.parent, b.parent);
}

inline void update_outer(Outer &out, const Item &child, double value,
                  const Outer::Link &link) {
  const double old = out.values.get(child);
  if (value < old) return;
  const auto found = out.links.find(child);
  if (value == old && found != out.links.end() && !better_link(link, found->second)) return;
  out.values.set(child, value);
  out.links.insert_or_assign(child, link);
}

// A max pool can discard subthreshold partial values: later larger offers
// suffice for re-entry. A sum pool could not forget those earlier summands.
struct MaxPool {
  struct Entry { double value; Outer::Link link; bool mandatory; };
  struct Rank { Item key; double value; };
  struct Better {
    bool operator()(const Rank &a, const Rank &b) const {
      if (a.value != b.value) return a.value > b.value;
      return less_key(a.key, b.key);
    }
  };
  using Entries = std::unordered_map<Item, Entry, sc::deductions::ItemKeyHash,
                                    sc::deductions::ItemKeyEq>;
  Entries entries;
  std::set<Rank, Better> novel;
  std::size_t evicted = 0;

  void offer(const Item &key, double value, const Outer::Link &link,
              bool mandatory, std::size_t width, bool prefilter = false) {
    const auto old = entries.find(key);
    // A new max candidate below the full frontier would be inserted and
    // immediately evicted. Preserve the tie rule and eviction count while
    // avoiding two node allocations and their hash/tree mutations.
    if (prefilter && old == entries.end() && !mandatory && novel.size() == width &&
        (width == 0 || !Better{}({key, value}, *std::prev(novel.end())))) {
      ++evicted;
      return;
    }
    if (old != entries.end()) {
      if (old->second.mandatory != mandatory) throw std::runtime_error("pool ownership changed");
      if (value < old->second.value ||
          (value == old->second.value && !better_link(link, old->second.link))) return;
      // An equal-score predecessor update leaves the rank unchanged. Avoid
      // erasing/reallocating its tree node just to replace this link.
      if (value == old->second.value) {
        old->second.link = link;
        return;
      }
      if (!mandatory) novel.erase({key, old->second.value});
      old->second = Entry{value, link, mandatory};
    } else {
      entries.insert_or_assign(key, Entry{value, link, mandatory});
    }
    if (!mandatory) {
      novel.insert({key, value});
      if (novel.size() > width) {
        const auto worst = std::prev(novel.end());
        entries.erase(worst->key);
        novel.erase(worst);
        ++evicted;
      }
    }
  }
};

} // namespace cparty::scfg2::beam::outer_detail

#endif

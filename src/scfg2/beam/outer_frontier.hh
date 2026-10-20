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
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
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
  static constexpr std::size_t no_heap_index = std::numeric_limits<std::size_t>::max();
  struct Entry {
    double value;
    Outer::Link link;
    bool mandatory;
    std::size_t heap_index = no_heap_index;
  };
  struct Rank { Item key; double value; };
  struct Better {
    bool operator()(const Rank &a, const Rank &b) const {
      if (a.value != b.value) return a.value > b.value;
      return less_key(a.key, b.key);
    }
  };
  using Entries = std::unordered_map<Item, Entry, sc::deductions::ItemKeyHash,
                                    sc::deductions::ItemKeyEq>;

  // A worst-first indexed heap replaces one allocation per rank-tree update.
  // Entry pointers survive unordered_map rehashes. Only pool copy/move needs
  // to rebind them, and swaps update their indices without another map probe.
  class NovelRanks {
   public:
    using const_iterator = std::vector<Rank>::const_iterator;
    std::size_t size() const { return heap_.size(); }
    bool empty() const { return heap_.empty(); }
    std::size_t capacity() const { return heap_.capacity(); }
    std::size_t view_capacity() const { return view_.capacity(); }

    // Preserve the diagnostic sorted-range API. Production asks only for
    // size, cutoff and clear, so it never builds this view. Any rank mutation
    // invalidates view iterators; repeated observations reuse the same view.
    const_iterator begin() const { return ranked().begin(); }
    const_iterator end() const { return ranked().end(); }
    void clear() {
      // Production moves entries to its completed-parent map first. Do not
      // dereference the old pointers here; release both histories entirely.
      std::vector<Cell>{}.swap(heap_);
      std::vector<Rank>{}.swap(view_);
      view_valid_ = true;
    }

   private:
    friend struct MaxPool;
    struct Cell { Rank rank; Entry *entry; };
    const Rank &worst() const { return heap_.front().rank; }
    void swap_cells(std::size_t a, std::size_t b) {
      std::swap(heap_[a], heap_[b]);
      heap_[a].entry->heap_index = a;
      heap_[b].entry->heap_index = b;
    }
    void sift_up(std::size_t position) {
      while (position) {
        const auto parent = (position - 1) / 2;
        if (!Better{}(heap_[parent].rank, heap_[position].rank)) break;
        swap_cells(parent, position);
        position = parent;
      }
    }
    void sift_down(std::size_t position) {
      while (position < heap_.size() / 2) {
        auto child = position * 2 + 1;
        if (child + 1 < heap_.size() &&
            Better{}(heap_[child].rank, heap_[child + 1].rank)) ++child;
        if (!Better{}(heap_[position].rank, heap_[child].rank)) break;
        swap_cells(position, child);
        position = child;
      }
    }
    void insert(const Rank &rank, Entry &entry) {
      const auto position = heap_.size();
      heap_.push_back({rank, &entry});
      entry.heap_index = position;
      sift_up(position);
      view_valid_ = false;
    }
    void improve(Entry &entry, double value) {
      const auto position = entry.heap_index;
      heap_[position].rank.value = value;
      // Scores only increase. A better rank moves away from the worst root.
      sift_down(position);
      view_valid_ = false;
    }
    void pop_erased_worst() {
      // The canonical entry has already been erased, matching the original
      // map mutation order. Never touch its now-dangling Entry pointer.
      if (heap_.size() == 1) heap_.pop_back();
      else {
        heap_.front() = heap_.back();
        heap_.pop_back();
        heap_.front().entry->heap_index = 0;
        sift_down(0);
      }
      view_valid_ = false;
    }
    void rebind(Entries &entries) {
      for (std::size_t p = 0; p < heap_.size(); ++p) {
        auto &entry = entries.at(heap_[p].rank.key);
        heap_[p].entry = &entry;
        entry.heap_index = p;
      }
    }
    const std::vector<Rank> &ranked() const {
      if (!view_valid_) {
        view_.clear();
        view_.reserve(heap_.size());
        for (const auto &cell : heap_) view_.push_back(cell.rank);
        std::sort(view_.begin(), view_.end(), Better{});
        view_valid_ = true;
      }
      return view_;
    }
    std::vector<Cell> heap_;
    mutable std::vector<Rank> view_;
    mutable bool view_valid_ = true;
  };

  Entries entries;
  NovelRanks novel;
  std::size_t evicted = 0;

  MaxPool() = default;
  MaxPool(const MaxPool &other)
      : entries(other.entries), novel(other.novel), evicted(other.evicted) {
    novel.rebind(entries);
  }
  MaxPool &operator=(const MaxPool &other) {
    if (this != &other) {
      entries = other.entries;
      novel = other.novel;
      evicted = other.evicted;
      novel.rebind(entries);
    }
    return *this;
  }
  MaxPool(MaxPool &&other)
      : entries(std::move(other.entries)), novel(std::move(other.novel)),
        evicted(other.evicted) {
    novel.rebind(entries);
    other.novel.clear();
  }
  MaxPool &operator=(MaxPool &&other) {
    if (this != &other) {
      entries = std::move(other.entries);
      novel = std::move(other.novel);
      evicted = other.evicted;
      novel.rebind(entries);
      other.novel.clear();
    }
    return *this;
  }

  void offer(const Item &key, double value, const Outer::Link &link,
              bool mandatory, std::size_t width, bool prefilter = false) {
    const auto old = entries.find(key);
    // A new max candidate below the full frontier would be inserted and
    // immediately evicted. Preserve the tie rule and eviction count while
    // avoiding a map-node allocation and the matching heap update.
    if (prefilter && old == entries.end() && !mandatory && novel.size() == width &&
        (width == 0 || !Better{}({key, value}, novel.worst()))) {
      ++evicted;
      return;
    }
    if (old != entries.end()) {
      if (old->second.mandatory != mandatory) throw std::runtime_error("pool ownership changed");
      if (value < old->second.value ||
          (value == old->second.value && !better_link(link, old->second.link))) return;
      // An equal-score predecessor update leaves the rank unchanged. Avoid
      // reordering its heap cell just to replace this link.
      if (value == old->second.value) {
        old->second.link = link;
        return;
      }
      const auto heap_index = old->second.heap_index;
      old->second = Entry{value, link, mandatory, heap_index};
      if (!mandatory) novel.improve(old->second, value);
    } else {
      const auto inserted = entries.insert_or_assign(key, Entry{value, link, mandatory});
      if (!mandatory) novel.insert({key, value}, inserted.first->second);
    }
    if (!mandatory) {
      if (novel.size() > width) {
        const auto worst = novel.worst();
        entries.erase(worst.key);
        novel.pop_erased_worst();
        ++evicted;
      }
    }
  }
};

} // namespace cparty::scfg2::beam::outer_detail

#endif

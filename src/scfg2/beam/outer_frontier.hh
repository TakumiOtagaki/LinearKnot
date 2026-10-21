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
#include <new>
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

  // Mirror only for closed internal callers. The canonical Entries remains
  // the sole authority for presence, ownership, values and iteration order.
  class EntryLookup {
   public:
    Entry *find(const Item &key) const {
      if (slots_.empty()) return nullptr;
      const auto p = position(slots_, key);
      return slots_[p].entry;
    }
    void set(const Item &key, Entry *entry) {
      if (slots_.empty()) rehash(8);
      auto p = position(slots_, key);
      if (slots_[p].entry) { slots_[p].entry = entry; return; }
      if (size_ >= slots_.size() / 2) {
        if (slots_.size() > slots_.max_size() / 2)
          throw std::length_error("pool lookup capacity overflow");
        rehash(slots_.size() * 2);
        p = position(slots_, key);
      }
      slots_[p] = {key, entry};
      ++size_;
    }
    void erase(const Item &key) {
      if (slots_.empty()) return;
      const auto p = position(slots_, key);
      if (!slots_[p].entry) return;
      const auto mask = slots_.size() - 1;
      auto hole = p;
      for (auto cursor = (hole + 1) & mask; slots_[cursor].entry;
           cursor = (cursor + 1) & mask) {
        const auto home = hash(slots_[cursor].key) & mask;
        if (((hole - home) & mask) < ((cursor - home) & mask)) {
          slots_[hole] = slots_[cursor];
          hole = cursor;
        }
      }
      slots_[hole] = {};
      if (!--size_) { clear(); return; }
      auto capacity = slots_.size();
      while (capacity > 8 && size_ < capacity / 8) capacity /= 2;
      if (capacity != slots_.size()) rehash(capacity);
    }
    void clear() noexcept { std::vector<Slot>{}.swap(slots_); size_ = 0; }

   private:
    struct Slot { Item key{}; Entry *entry = nullptr; };
    static std::size_t hash(const Item &key) {
      auto value = sc::deductions::ItemKeyHash{}(key);
      // Keep all five fields, including BE anchors and noncanonical primary
      // coordinates; mix upper bits before using a power-of-two slot count.
      value ^= value >> (std::numeric_limits<std::size_t>::digits / 2);
      return value;
    }
    static std::size_t position(const std::vector<Slot> &slots, const Item &key) {
      const auto mask = slots.size() - 1;
      auto p = hash(key) & mask;
      while (slots[p].entry && !sc::deductions::ItemKeyEq{}(slots[p].key, key))
        p = (p + 1) & mask;
      return p;
    }
    void rehash(std::size_t capacity) {
      std::vector<Slot> replacement(capacity);
      for (const auto &slot : slots_) {
        if (!slot.entry) continue;
        replacement[position(replacement, slot.key)] = slot;
      }
      slots_.swap(replacement);
    }
    std::vector<Slot> slots_;
    std::size_t size_ = 0;
  };

  // Only a validated canonical primary domain may use this representation.
  // Within a fixed (family,end) or (family,span) pool, start determines the
  // complete key. Every map Entry pointer remains authoritative and stable.
  class PrimaryEntryLookup {
   public:
    Entry *find(cand_pos_t start) const {
      if (slots_.empty()) return nullptr;
      return slots_[position(slots_, start)].entry;
    }
    void set(cand_pos_t start, Entry *entry) {
      if (slots_.empty()) rehash(8);
      auto p = position(slots_, start);
      if (slots_[p].entry) { slots_[p].entry = entry; return; }
      if (size_ >= slots_.size() / 2) {
        if (slots_.size() > slots_.max_size() / 2)
          throw std::length_error("primary pool lookup capacity overflow");
        rehash(slots_.size() * 2);
        p = position(slots_, start);
      }
      slots_[p] = {start, entry};
      ++size_;
    }
    void erase(cand_pos_t start) {
      if (slots_.empty()) return;
      const auto p = position(slots_, start);
      if (!slots_[p].entry) return;
      const auto mask = slots_.size() - 1;
      auto hole = p;
      for (auto cursor = (hole + 1) & mask; slots_[cursor].entry;
           cursor = (cursor + 1) & mask) {
        const auto home = hash(slots_[cursor].start) & mask;
        if (((hole - home) & mask) < ((cursor - home) & mask)) {
          slots_[hole] = slots_[cursor];
          hole = cursor;
        }
      }
      slots_[hole] = {};
      if (!--size_) { clear(); return; }
      auto capacity = slots_.size();
      while (capacity > 8 && size_ < capacity / 8) capacity /= 2;
      if (capacity != slots_.size()) rehash(capacity);
    }
    void clear() noexcept { std::vector<Slot>{}.swap(slots_); size_ = 0; }

   private:
    struct Slot { cand_pos_t start = 0; Entry *entry = nullptr; };
    static std::size_t hash(cand_pos_t start) {
      auto value = static_cast<std::uint64_t>(start) * UINT64_C(0x9e3779b97f4a7c15);
      value ^= value >> 32;
      return static_cast<std::size_t>(value);
    }
    static std::size_t position(const std::vector<Slot> &slots, cand_pos_t start) {
      const auto mask = slots.size() - 1;
      auto p = hash(start) & mask;
      while (slots[p].entry && slots[p].start != start) p = (p + 1) & mask;
      return p;
    }
    void rehash(std::size_t capacity) {
      std::vector<Slot> replacement(capacity);
      for (const auto &slot : slots_) {
        if (!slot.entry) continue;
        replacement[position(replacement, slot.start)] = slot;
      }
      slots_.swap(replacement);
    }
    std::vector<Slot> slots_;
    std::size_t size_ = 0;
  };

  Entries entries;
  NovelRanks novel;
  std::size_t evicted = 0;

  // Trusted internal lifecycle only: after opt-in all insertions/erasures
  // must go through offer(). Disable before ANY direct entries mutation,
  // assignment, swap or move-out. Public/default pools keep ordinary find.
  // A MaxPool copy starts with no mirror; a move disables the source mirror.
  void enable_trusted_index() noexcept {
    disable_trusted_index();
    try {
      for (auto &[key, entry] : entries) lookup_.set(key, &entry);
      trusted_lookup_ = true;
    } catch (const std::bad_alloc &) { disable_trusted_index(); }
      catch (const std::length_error &) { disable_trusted_index(); }
  }
  void disable_trusted_index() noexcept {
    trusted_lookup_ = false;
    lookup_.clear();
    primary_lookup_.clear();
    primary_domain_ = PrimaryDomain::None;
  }

  void enable_trusted_primary_end_index(NT family, cand_pos_t end) noexcept {
    enable_trusted_primary_index(family, end, PrimaryDomain::End);
  }
  void enable_trusted_primary_span_index(NT family, cand_pos_t span) noexcept {
    enable_trusted_primary_index(family, span, PrimaryDomain::Span);
  }

  MaxPool() = default;
  MaxPool(const MaxPool &other)
      : entries(other.entries), novel(other.novel), evicted(other.evicted) {
    novel.rebind(entries);
  }
  MaxPool &operator=(const MaxPool &other) {
    if (this != &other) {
      disable_trusted_index();
      entries = other.entries;
      novel = other.novel;
      evicted = other.evicted;
      novel.rebind(entries);
    }
    return *this;
  }
  MaxPool(MaxPool &&other)
      : entries((other.disable_trusted_index(), std::move(other.entries))),
        novel(std::move(other.novel)),
        evicted(other.evicted) {
    novel.rebind(entries);
    other.disable_trusted_index();
    other.novel.clear();
  }
  MaxPool &operator=(MaxPool &&other) {
    if (this != &other) {
      disable_trusted_index();
      other.disable_trusted_index();
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
    Entry *old = find_entry(key);
    // A new max candidate below the full frontier would be inserted and
    // immediately evicted. Preserve the tie rule and eviction count while
    // avoiding a map-node allocation and the matching heap update.
    if (prefilter && old == nullptr && !mandatory && novel.size() == width &&
        (width == 0 || !Better{}({key, value}, novel.worst()))) {
      ++evicted;
      return;
    }
    if (old != nullptr) {
      if (old->mandatory != mandatory) throw std::runtime_error("pool ownership changed");
      if (value < old->value ||
          (value == old->value && !better_link(link, old->link))) return;
      // An equal-score predecessor update leaves the rank unchanged. Avoid
      // reordering its heap cell just to replace this link.
      if (value == old->value) {
        old->link = link;
        return;
      }
      const auto heap_index = old->heap_index;
      *old = Entry{value, link, mandatory, heap_index};
      if (!mandatory) novel.improve(*old, value);
    } else {
      const auto inserted = entries.insert_or_assign(key, Entry{value, link, mandatory});
      if (trusted_lookup_) {
        try {
          if (primary_domain_ == PrimaryDomain::None)
            lookup_.set(key, &inserted.first->second);
          else primary_lookup_.set(key.i, &inserted.first->second);
        }
        catch (const std::bad_alloc &) { disable_trusted_index(); }
        catch (const std::length_error &) { disable_trusted_index(); }
      }
      if (!mandatory) novel.insert({key, value}, inserted.first->second);
    }
    if (!mandatory) {
      if (novel.size() > width) {
        const auto worst = novel.worst();
        if (trusted_lookup_) {
          try {
            if (primary_domain_ == PrimaryDomain::None) lookup_.erase(worst.key);
            else primary_lookup_.erase(worst.key.i);
          }
          catch (const std::bad_alloc &) { disable_trusted_index(); }
          catch (const std::length_error &) { disable_trusted_index(); }
        }
        entries.erase(worst.key);
        novel.pop_erased_worst();
        ++evicted;
      }
    }
  }

 private:
  enum class PrimaryDomain { None, End, Span };
  bool primary_domain_matches(const Item &key) const {
    return key.nonterminal == primary_family_ && key.ip == -1 && key.jp == -1 &&
        key.i >= 1 && key.j >= key.i &&
        (primary_domain_ == PrimaryDomain::End ? key.j == primary_extent_
                                             : key.j - key.i == primary_extent_);
  }
  Entry *find_entry(const Item &key) {
    if (trusted_lookup_) {
      if (primary_domain_ == PrimaryDomain::None) return lookup_.find(key);
      if (primary_domain_matches(key)) return primary_lookup_.find(key.i);
      // A malformed or differently grouped offer/query must never alias a
      // live start. The map Entry pointer remains valid when metadata clears.
      disable_trusted_index();
    }
    const auto found = entries.find(key);
    return found == entries.end() ? nullptr : &found->second;
  }
  void enable_trusted_primary_index(NT family, cand_pos_t extent,
                                     PrimaryDomain domain) noexcept {
    disable_trusted_index();
    if (static_cast<std::size_t>(family) >= static_cast<std::size_t>(NT::BE) ||
        (domain == PrimaryDomain::End ? extent < 1 : extent < 0)) {
      enable_trusted_index();
      return;
    }
    primary_domain_ = domain;
    primary_family_ = family;
    primary_extent_ = extent;
    try {
      for (auto &[key, entry] : entries) {
        if (!primary_domain_matches(key)) { disable_trusted_index(); return; }
        primary_lookup_.set(key.i, &entry);
      }
      trusted_lookup_ = true;
    } catch (const std::bad_alloc &) { disable_trusted_index(); }
      catch (const std::length_error &) { disable_trusted_index(); }
  }
  EntryLookup lookup_;
  PrimaryEntryLookup primary_lookup_;
  bool trusted_lookup_ = false;
  PrimaryDomain primary_domain_ = PrimaryDomain::None;
  NT primary_family_ = NT::W;
  cand_pos_t primary_extent_ = 0;
};

} // namespace cparty::scfg2::beam::outer_detail

#endif

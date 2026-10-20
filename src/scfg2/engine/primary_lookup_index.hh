#ifndef CPARTY_SCFG2_ENGINE_PRIMARY_LOOKUP_INDEX_HH
#define CPARTY_SCFG2_ENGINE_PRIMARY_LOOKUP_INDEX_HH

#include "scfg2/deductions/item_key.hh"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cparty::scfg2::engine {

// Numeric point-lookup cache only. The owning chart retains the authoritative
// keys, values, presence flags and iteration order. Each primary (family,j)
// row hashes its positive start coordinate; no n-by-n index is allocated.
class PrimaryLookupIndex {
 public:
  using key_type = deductions::ItemKey;

  explicit PrimaryLookupIndex(cand_pos_t n) : n_(n) {
    if (n < 0) throw std::invalid_argument("negative primary lookup length");
    const auto length = static_cast<std::size_t>(n);
    if (length > rows_.max_size() / primary_families)
      throw std::length_error("primary lookup dimensions overflow");
    rows_.resize(primary_families * length);
  }

  PrimaryLookupIndex(const PrimaryLookupIndex &other)
      : n_(other.n_), rows_(other.rows_), size_(other.size_),
        active_rows_(other.active_rows_) {
    // vector copies need not preserve spare capacity. Count the actual
    // allocations of this independent copy, rather than copying byte totals.
    for (const auto &row : rows_) allocated_capacity_ += row.slots.capacity();
  }
  PrimaryLookupIndex &operator=(const PrimaryLookupIndex &other) {
    if (this != &other) {
      PrimaryLookupIndex copy(other);
      swap(copy);
    }
    return *this;
  }
  PrimaryLookupIndex(PrimaryLookupIndex &&other) noexcept
      : n_(std::exchange(other.n_, 0)), rows_(std::move(other.rows_)),
        size_(std::exchange(other.size_, 0)),
        active_rows_(std::exchange(other.active_rows_, 0)),
        allocated_capacity_(std::exchange(other.allocated_capacity_, 0)) {
    std::vector<Row>{}.swap(other.rows_);
  }
  PrimaryLookupIndex &operator=(PrimaryLookupIndex &&other) noexcept {
    if (this == &other) return *this;
    n_ = std::exchange(other.n_, 0);
    rows_ = std::move(other.rows_);
    size_ = std::exchange(other.size_, 0);
    active_rows_ = std::exchange(other.active_rows_, 0);
    allocated_capacity_ = std::exchange(other.allocated_capacity_, 0);
    std::vector<Row>{}.swap(other.rows_);
    return *this;
  }

  cand_pos_t length() const { return n_; }
  std::size_t size() const { return size_; }
  std::size_t row_count() const { return rows_.size(); }
  std::size_t active_rows() const { return active_rows_; }
  std::size_t allocated_capacity() const { return allocated_capacity_; }
  std::size_t row_capacity(const key_type &key) const {
    const auto index = row_index(key);
    return index ? rows_[*index].slots.capacity() : 0;
  }
  // O(1) dynamic storage observation, including unused allocated slots and
  // row headers. The PrimaryLookupIndex object itself is not included.
  std::size_t memory_bytes() const {
    return rows_.capacity() * sizeof(Row) + allocated_capacity_ * sizeof(Slot);
  }

  // Unsupported full keys use the owning chart's general lookup. Supported
  // but absent keys return a value, preserving any default-value bit pattern.
  std::optional<double> get(const key_type &key, double default_value) const {
    const auto index = row_index(key);
    if (!index) return std::nullopt;
    const auto &row = rows_[*index];
    if (row.slots.empty()) return default_value;
    const auto position = find(row, key.i);
    return row.slots[position].start == key.i
        ? row.slots[position].value : default_value;
  }

  void set(const key_type &key, double value) {
    const auto index = row_index(key);
    if (!index) return;
    auto &row = rows_[*index];
    if (row.slots.empty()) rehash(row, minimum_capacity);
    auto position = find(row, key.i);
    if (row.slots[position].start == key.i) {
      row.slots[position].value = value;
      return;
    }
    // At most half the logical slots are occupied, so every probe terminates
    // and expected lookup/erase cost stays bounded independently of n.
    if (row.count >= row.slots.size() / 2) {
      if (row.slots.size() > row.slots.max_size() / 2)
        throw std::length_error("primary lookup row capacity overflow");
      rehash(row, row.slots.size() * 2);
      position = find(row, key.i);
    }
    row.slots[position] = Slot{key.i, value};
    if (!row.count) ++active_rows_;
    ++row.count;
    ++size_;
  }

  void erase(const key_type &key) {
    const auto index = row_index(key);
    if (!index) return;
    auto &row = rows_[*index];
    if (row.slots.empty()) return;
    const auto position = find(row, key.i);
    if (row.slots[position].start != key.i) return;

    // Backward-shift deletion closes the hole within its probe cluster.
    // No tombstones survive repeated prune/reinsert cycles, including a
    // cluster which wraps from the last slot back to zero.
    const auto mask = row.slots.size() - 1;
    auto hole = position;
    for (auto cursor = (hole + 1) & mask; row.slots[cursor].start;
         cursor = (cursor + 1) & mask) {
      const auto home = bucket(row.slots[cursor].start, mask);
      if (((hole - home) & mask) < ((cursor - home) & mask)) {
        row.slots[hole] = row.slots[cursor];
        hole = cursor;
      }
    }
    row.slots[hole] = Slot{};
    --row.count;
    --size_;
    if (!row.count) {
      allocated_capacity_ -= row.slots.capacity();
      std::vector<Slot>{}.swap(row.slots);
      --active_rows_;
      return;
    }

    // Old columns must not retain their candidate-phase peak capacity.
    // Geometric shrinking maintains capacity <= 8*live slots per nonempty
    // row, with an eight-slot minimum, rather than accumulating O(n^2).
    auto capacity = row.slots.size();
    while (capacity > minimum_capacity && row.count < capacity / 8)
      capacity /= 2;
    if (capacity != row.slots.size()) rehash(row, capacity);
  }

  void swap(PrimaryLookupIndex &other) noexcept {
    using std::swap;
    swap(n_, other.n_);
    rows_.swap(other.rows_);
    swap(size_, other.size_);
    swap(active_rows_, other.active_rows_);
    swap(allocated_capacity_, other.allocated_capacity_);
  }

 private:
  static constexpr std::size_t primary_families =
      static_cast<std::size_t>(grammar::NonTerminal::BE);
  static constexpr std::size_t minimum_capacity = 8;
  struct Slot {
    cand_pos_t start = 0;  // Zero is not a supported nucleotide coordinate.
    double value = 0.0;
  };
  struct Row {
    std::vector<Slot> slots;
    std::size_t count = 0;
  };

  std::optional<std::size_t> row_index(const key_type &key) const {
    const auto family = static_cast<std::size_t>(key.nonterminal);
    if (family >= primary_families || key.ip != -1 || key.jp != -1 ||
        key.i < 1 || key.j < key.i || key.j > n_) return std::nullopt;
    return family * static_cast<std::size_t>(n_) + key.j - 1;
  }
  static std::size_t bucket(cand_pos_t start, std::size_t mask) {
    // Integer-only hashing; the row already identifies family and end.
    auto hash = static_cast<std::uint64_t>(start) * UINT64_C(0x9e3779b97f4a7c15);
    hash ^= hash >> 32;
    return static_cast<std::size_t>(hash) & mask;
  }
  static std::size_t find(const Row &row, cand_pos_t start) {
    const auto mask = row.slots.size() - 1;
    auto position = bucket(start, mask);
    while (row.slots[position].start && row.slots[position].start != start)
      position = (position + 1) & mask;
    return position;
  }
  void rehash(Row &row, std::size_t capacity) {
    std::vector<Slot> replacement(capacity);
    const auto mask = capacity - 1;
    for (const auto &slot : row.slots) {
      if (!slot.start) continue;
      auto position = bucket(slot.start, mask);
      while (replacement[position].start) position = (position + 1) & mask;
      replacement[position] = slot;
    }
    const auto old_capacity = row.slots.capacity();
    row.slots.swap(replacement);
    allocated_capacity_ -= old_capacity;
    allocated_capacity_ += row.slots.capacity();
  }

  cand_pos_t n_;
  std::vector<Row> rows_;
  std::size_t size_ = 0;
  std::size_t active_rows_ = 0;
  std::size_t allocated_capacity_ = 0;
};

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_PRIMARY_LOOKUP_INDEX_HH

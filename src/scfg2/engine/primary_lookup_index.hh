#ifndef CPARTY_SCFG2_ENGINE_PRIMARY_LOOKUP_INDEX_HH
#define CPARTY_SCFG2_ENGINE_PRIMARY_LOOKUP_INDEX_HH

#include "scfg2/deductions/item_key.hh"

#include <cstddef>
#include <cstdint>
#include <new>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cparty::scfg2::engine {

// Numeric point-lookup cache only. The owning chart retains the authoritative
// keys, values, presence flags and iteration order. Each primary (family,j)
// row hashes its positive start coordinate, or uses a bounded direct interval
// after selection. No n-by-n index is allocated.
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
    for (const auto &row : rows_) {
      allocated_capacity_ += row.starts.capacity();
      allocated_value_capacity_ += row.values.capacity();
    }
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
        allocated_capacity_(std::exchange(other.allocated_capacity_, 0)),
        allocated_value_capacity_(std::exchange(other.allocated_value_capacity_, 0)) {
    std::vector<Row>{}.swap(other.rows_);
  }
  PrimaryLookupIndex &operator=(PrimaryLookupIndex &&other) noexcept {
    if (this == &other) return *this;
    n_ = std::exchange(other.n_, 0);
    rows_ = std::move(other.rows_);
    size_ = std::exchange(other.size_, 0);
    active_rows_ = std::exchange(other.active_rows_, 0);
    allocated_capacity_ = std::exchange(other.allocated_capacity_, 0);
    allocated_value_capacity_ = std::exchange(other.allocated_value_capacity_, 0);
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
    return index ? rows_[*index].starts.capacity() : 0;
  }
  // O(1) dynamic storage observation, including unused allocated slots and
  // row headers. The PrimaryLookupIndex object itself is not included.
  std::size_t memory_bytes() const {
    return rows_.capacity() * sizeof(Row) + allocated_capacity_ * sizeof(cand_pos_t) +
        allocated_value_capacity_ * sizeof(double);
  }

  // Unsupported full keys use the owning chart's general lookup. Supported
  // but absent keys return a value, preserving any default-value bit pattern.
  std::optional<double> get(const key_type &key, double default_value) const {
    const auto index = row_index(key);
    if (!index) return std::nullopt;
    const auto &row = rows_[*index];
    if (row.starts.empty()) return default_value;
    if (row.first_start > 0) {
      if (key.i < row.first_start) return default_value;
      const auto position = static_cast<std::size_t>(key.i - row.first_start);
      return position < row.starts.size() && row.starts[position]
          ? row.values[position] : default_value;
    }
    const auto position = find(row, key.i);
    return row.starts[position] == key.i
        ? row.values[position] : default_value;
  }

  void set(const key_type &key, double value) {
    const auto index = row_index(key);
    if (!index) return;
    auto &row = rows_[*index];
    make_mutable(row);
    if (row.starts.empty()) rehash(row, minimum_capacity);
    auto position = find(row, key.i);
    if (row.starts[position] == key.i) {
      row.values[position] = value;
      return;
    }
    // At most half the logical slots are occupied, so every probe terminates
    // and expected lookup/erase cost stays bounded independently of n.
    if (row.count >= row.starts.size() / 2) {
      if (row.starts.size() > row.starts.max_size() / 2 ||
          row.values.size() > row.values.max_size() / 2)
        throw std::length_error("primary lookup row capacity overflow");
      rehash(row, row.starts.size() * 2);
      position = find(row, key.i);
    }
    row.starts[position] = key.i;
    row.values[position] = value;
    if (!row.count) ++active_rows_;
    ++row.count;
    ++size_;
  }

  void erase(const key_type &key) {
    const auto index = row_index(key);
    if (!index) return;
    auto &row = rows_[*index];
    if (row.starts.empty()) return;
    make_mutable(row);
    const auto position = find(row, key.i);
    if (row.starts[position] != key.i) return;

    // Backward-shift deletion closes the hole within its probe cluster.
    // No tombstones survive repeated prune/reinsert cycles, including a
    // cluster which wraps from the last slot back to zero.
    const auto mask = row.starts.size() - 1;
    auto hole = position;
    for (auto cursor = (hole + 1) & mask; row.starts[cursor];
         cursor = (cursor + 1) & mask) {
      const auto home = bucket(row.starts[cursor], mask);
      if (((hole - home) & mask) < ((cursor - home) & mask)) {
        row.starts[hole] = row.starts[cursor];
        row.values[hole] = row.values[cursor];
        hole = cursor;
      }
    }
    row.starts[hole] = 0;
    row.values[hole] = 0.0;
    --row.count;
    --size_;
    if (!row.count) {
      allocated_capacity_ -= row.starts.capacity();
      allocated_value_capacity_ -= row.values.capacity();
      std::vector<cand_pos_t>{}.swap(row.starts);
      std::vector<double>{}.swap(row.values);
      --active_rows_;
      return;
    }

    // Old columns must not retain their candidate-phase peak capacity.
    // Geometric shrinking maintains capacity <= 8*live slots per nonempty
    // row, with an eight-slot minimum, rather than accumulating O(n^2).
    auto capacity = row.starts.size();
    while (capacity > minimum_capacity && row.count < capacity / 8)
      capacity /= 2;
    if (capacity != row.starts.size()) rehash(row, capacity);
  }

  // Explicit post-selection observation boundary, attempted at most once
  // per row. A direct interval is never larger than this row's current hash
  // table, whose post-prune capacity is bounded by eight times its live keys.
  // No chart record, query order or arithmetic changes representation here.
  void seal_row(grammar::NonTerminal family, cand_pos_t end) {
    const auto index = row_index(key_type{family, 1, end, -1, -1});
    if (!index) return;
    auto &row = rows_[*index];
    if (row.seal_attempted) return;
    row.seal_attempted = true;
    if (!row.count) return;
    cand_pos_t first = end, last = 0;
    for (const auto start : row.starts) {
      if (!start) continue;
      if (start < first) first = start;
      if (start > last) last = start;
    }
    const auto extent = static_cast<std::size_t>(last) -
        static_cast<std::size_t>(first) + 1;
    if (extent > row.starts.size()) return;
    try {
      std::vector<cand_pos_t> starts(extent, 0);
      std::vector<double> values(extent, 0.0);
      for (std::size_t p = 0; p < row.starts.size(); ++p) {
        const auto start = row.starts[p];
        if (!start) continue;
        const auto position = static_cast<std::size_t>(start - first);
        starts[position] = start; // Presence is independent of every value bit.
        values[position] = row.values[p];
      }
      const auto old_capacity = row.starts.capacity();
      const auto old_value_capacity = row.values.capacity();
      row.starts.swap(starts);
      row.values.swap(values);
      allocated_capacity_ -= old_capacity;
      allocated_capacity_ += row.starts.capacity();
      allocated_value_capacity_ -= old_value_capacity;
      allocated_value_capacity_ += row.values.capacity();
      row.first_start = first;
    } catch (const std::bad_alloc &) {
    } catch (const std::length_error &) {
    }
  }

  void swap(PrimaryLookupIndex &other) noexcept {
    using std::swap;
    swap(n_, other.n_);
    rows_.swap(other.rows_);
    swap(size_, other.size_);
    swap(active_rows_, other.active_rows_);
    swap(allocated_capacity_, other.allocated_capacity_);
    swap(allocated_value_capacity_, other.allocated_value_capacity_);
  }

 private:
  static constexpr std::size_t primary_families =
      static_cast<std::size_t>(grammar::NonTerminal::BE);
  static constexpr std::size_t minimum_capacity = 8;
  struct Row {
    // Probe keys without bringing the value/padding of every collided slot
    // into cache. Both arrays have the same logical size; hashed rows use powers of two.
    // Zero is not a supported nucleotide coordinate.
    std::vector<cand_pos_t> starts;
    std::vector<double> values;
    std::size_t count = 0;
    cand_pos_t first_start = 0; // Positive only for the direct interval mode.
    bool seal_attempted = false;
  };

  void make_mutable(Row &row) {
    if (!row.first_start) return;
    auto capacity = minimum_capacity;
    while (row.count > capacity / 2) {
      if (capacity > row.starts.max_size() / 2 ||
          capacity > row.values.max_size() / 2)
        throw std::length_error("primary lookup row capacity overflow");
      capacity *= 2;
    }
    // rehash accepts direct or hashed input and commits only after both
    // replacement allocations succeed. Keep seal_attempted: later updates
    // stay hashed, so each row pays at most one promotion and one demotion.
    rehash(row, capacity);
  }

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
    const auto mask = row.starts.size() - 1;
    auto position = bucket(start, mask);
    while (row.starts[position] && row.starts[position] != start)
      position = (position + 1) & mask;
    return position;
  }
  void rehash(Row &row, std::size_t capacity) {
    // Allocate and populate both replacements before either owner changes.
    // Failure of the second allocation leaves the current index untouched.
    std::vector<cand_pos_t> replacement_starts(capacity, 0);
    std::vector<double> replacement_values(capacity, 0.0);
    const auto mask = capacity - 1;
    for (std::size_t p = 0; p < row.starts.size(); ++p) {
      const auto start = row.starts[p];
      if (!start) continue;
      auto position = bucket(start, mask);
      while (replacement_starts[position]) position = (position + 1) & mask;
      replacement_starts[position] = start;
      replacement_values[position] = row.values[p];
    }
    const auto old_capacity = row.starts.capacity();
    const auto old_value_capacity = row.values.capacity();
    row.starts.swap(replacement_starts);
    row.values.swap(replacement_values);
    allocated_capacity_ -= old_capacity;
    allocated_capacity_ += row.starts.capacity();
    allocated_value_capacity_ -= old_value_capacity;
    allocated_value_capacity_ += row.values.capacity();
    row.first_start = 0;
  }

  cand_pos_t n_;
  std::vector<Row> rows_;
  std::size_t size_ = 0;
  std::size_t active_rows_ = 0;
  // The public slot-capacity observation counts key slots. Account for the
  // two actual allocation capacities separately in memory_bytes().
  std::size_t allocated_capacity_ = 0;
  std::size_t allocated_value_capacity_ = 0;
};

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_PRIMARY_LOOKUP_INDEX_HH

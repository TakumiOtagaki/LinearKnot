#ifndef CPARTY_SCFG2_BEAM_COLUMN_STORE_HH
#define CPARTY_SCFG2_BEAM_COLUMN_STORE_HH

#include "scfg2/deductions/item_key.hh"

#include <cstddef>
#include <memory>
#include <new>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

// Authoritative storage for one private agenda column. Public chart identity,
// pruning, local weights and semiring checkpoints remain with their owners.
// Canonical keys use start-indexed slots; every other complete key retains the
// old map behavior through the sparse fallback. No iterator feeds execution.
template <class Value, grammar::NonTerminal Family, std::size_t InitialCapacityHint = 0>
class ColumnStore {
 public:
  using key_type = deductions::ItemKey;
  using mapped_type = Value;
  using map_type = std::unordered_map<key_type, Value,
      deductions::ItemKeyHash, deductions::ItemKeyEq>;
  using value_type = typename map_type::value_type;
  using Entry = value_type;

  ColumnStore() = default;
  ColumnStore(const ColumnStore &) = default;
  ColumnStore &operator=(const ColumnStore &other) {
    if (this != &other) {
      ColumnStore copy(other);
      swap(copy);
    }
    return *this;
  }
  ColumnStore(ColumnStore &&other) noexcept
      : column_(std::exchange(other.column_, -1)),
        direct_ready_(std::exchange(other.direct_ready_, false)),
        slots_(std::move(other.slots_)), touched_(std::move(other.touched_)),
        fallback_(std::move(other.fallback_)) {
    std::vector<std::optional<Entry>>{}.swap(other.slots_);
    std::vector<cand_pos_t>{}.swap(other.touched_);
    other.fallback_.clear();
  }
  ColumnStore &operator=(ColumnStore &&other) noexcept {
    if (this != &other) {
      ColumnStore moved(std::move(other));
      swap(moved);
    }
    return *this;
  }

  // Grow only at a column boundary. Returned Entry/Value pointers are stable
  // throughout insertion and lookup within that column, until clear/reset.
  // Empty optionals construct no Values, including graph/tape identities.
  void reset_column(cand_pos_t column) {
    clear();
    direct_ready_ = false;
    if (column >= 0) {
      const auto extent = static_cast<std::size_t>(column) + 1;
      try {
        if (extent > slots_.capacity()) {
          auto capacity = slots_.capacity() ? slots_.capacity() : 1;
          while (capacity < extent) {
            if (capacity > slots_.max_size() / 2) { capacity = extent; break; }
            capacity *= 2;
          }
          slots_.reserve(capacity);
        }
        if (extent > slots_.size()) slots_.resize(extent);
        // Match the geometric slot capacity, rather than reserving every j.
        // No bookkeeping allocation is then needed during this column.
        touched_.reserve(slots_.capacity());
        direct_ready_ = true;
      } catch (const std::bad_alloc &) {
      } catch (const std::length_error &) {
      }
    }
    column_ = column;
  }
  void begin_column(cand_pos_t column) { reset_column(column); }
  void clear() {
    for (const auto start : touched_) slots_[start].reset();
    touched_.clear();
    fallback_.clear();
  }
  cand_pos_t column() const { return column_; }
  std::size_t size() const { return touched_.size() + fallback_.size(); }
  bool empty() const { return size() == 0; }
  static constexpr std::nullptr_t end() { return nullptr; }

  Entry *find(const key_type &key) {
    if (canonical(key)) {
      auto &slot = slots_[key.i];
      return slot ? std::addressof(*slot) : nullptr;
    }
    const auto found = fallback_.find(key);
    return found == fallback_.end() ? nullptr : std::addressof(*found);
  }
  const Entry *find(const key_type &key) const {
    if (canonical(key)) {
      const auto &slot = slots_[key.i];
      return slot ? std::addressof(*slot) : nullptr;
    }
    const auto found = fallback_.find(key);
    return found == fallback_.end() ? nullptr : std::addressof(*found);
  }

  template <class... Args>
  std::pair<Entry *, bool> try_emplace(const key_type &key, Args &&...args) {
    if (!canonical(key)) {
      auto [found, inserted] = fallback_.try_emplace(key, std::forward<Args>(args)...);
      if (inserted) reserve_initial_capacity(found->second);
      return {std::addressof(*found), inserted};
    }
    auto &slot = slots_[key.i];
    if (slot) return {std::addressof(*slot), false};
    // Record capacity first. A throwing Value constructor leaves this slot
    // absent and rolls back the touched record, without a default Value.
    touched_.push_back(key.i);
    try {
      slot.emplace(std::piecewise_construct, std::forward_as_tuple(key),
          std::forward_as_tuple(std::forward<Args>(args)...));
    } catch (...) {
      touched_.pop_back();
      throw;
    }
    reserve_initial_capacity(slot->second);
    return {std::addressof(*slot), true};
  }
  Value &operator[](const key_type &key) { return try_emplace(key).first->second; }

  void swap(ColumnStore &other) noexcept {
    using std::swap;
    swap(column_, other.column_);
    swap(direct_ready_, other.direct_ready_);
    slots_.swap(other.slots_);
    touched_.swap(other.touched_);
    fallback_.swap(other.fallback_);
  }

 private:
  // Optional allocation advice for explicitly opted-in vector-like Values.
  // Generic semiring Values never enter this path or gain a reserve contract.
  static void reserve_initial_capacity(Value &value) {
    if constexpr (InitialCapacityHint != 0) {
      if (value.empty() && value.capacity() == 0) {
        try { value.reserve(InitialCapacityHint); }
        catch (const std::bad_alloc &) {}
        catch (const std::length_error &) {}
      }
    }
  }
  bool canonical(const key_type &key) const {
    return direct_ready_ && key.nonterminal == Family && key.j == column_ && key.ip == -1 &&
        key.jp == -1 && key.i >= 1 && key.i <= column_;
  }
  cand_pos_t column_ = -1;
  bool direct_ready_ = false;
  std::vector<std::optional<Entry>> slots_;
  std::vector<cand_pos_t> touched_;
  map_type fallback_;
};

}  // namespace cparty::scfg2::beam
#endif

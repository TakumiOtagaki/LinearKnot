#ifndef CPARTY_SCFG2_ENGINE_CHART_HH
#define CPARTY_SCFG2_ENGINE_CHART_HH

#include "scfg2/deductions/item_key.hh"
#include "scfg2/engine/primary_lookup_index.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <unordered_map>
#include <cstdint>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace cparty::scfg2::engine {

template <class Value>
class Chart {
 public:
  using value_type = Value;
  using key_type = deductions::ItemKey;
  using map_type = std::unordered_map<key_type, Value, deductions::ItemKeyHash, deductions::ItemKeyEq>;
  using discard_fn = bool (*)(const Value &);

  // An opt-in predicate may discard semiring zero entries. The default
  // preserves exact callers' historical contains()/zero-entry contract.
  explicit Chart(Value default_value, discard_fn discard = nullptr,
                 cand_pos_t dense_length = 0, cand_pos_t primary_lookup_length = 0)
      : default_value_(default_value), discard_(discard) {
    if (dense_length < 0) throw std::invalid_argument("negative dense chart length");
    if (primary_lookup_length < 0) throw std::invalid_argument("negative primary lookup length");
    if (dense_length && primary_lookup_length)
      throw std::invalid_argument("dense storage and primary lookup are mutually exclusive");
    if (dense_length) dense_.emplace(dense_length, default_value);
    if constexpr (std::is_same_v<Value, double>) {
      if (primary_lookup_length) {
        try { primary_lookup_.emplace(primary_lookup_length); }
        catch (const std::bad_alloc &) {}
        catch (const std::length_error &) {}
      }
    } else if (primary_lookup_length) {
      throw std::invalid_argument("primary lookup requires double values");
    }
  }

  Chart(const Chart &other)
      : default_value_(other.default_value_), discard_(other.discard_),
        values_(other.values_), dense_(other.dense_), materialized_(other.materialized_) {
    copy_primary_lookup(other);
  }
  Chart &operator=(const Chart &other) {
    if (this == &other) return *this;
    primary_lookup_.reset();
    default_value_ = other.default_value_;
    discard_ = other.discard_;
    values_ = other.values_;
    dense_ = other.dense_;
    materialized_ = other.materialized_;
    copy_primary_lookup(other);
    return *this;
  }
  Chart(Chart &&other)
      : default_value_(other.default_value_), discard_(other.discard_),
        values_(std::move(other.values_)), dense_(std::move(other.dense_)),
        materialized_(std::move(other.materialized_)),
        primary_lookup_(std::move(other.primary_lookup_)) {
    other.dense_.reset();
    other.materialized_.reset();
    other.primary_lookup_.reset();
  }
  Chart &operator=(Chart &&other) {
    if (this == &other) return *this;
    primary_lookup_.reset();
    default_value_ = other.default_value_;
    discard_ = other.discard_;
    values_ = std::move(other.values_);
    dense_ = std::move(other.dense_);
    materialized_ = std::move(other.materialized_);
    primary_lookup_ = std::move(other.primary_lookup_);
    other.dense_.reset();
    other.materialized_.reset();
    other.primary_lookup_.reset();
    return *this;
  }

  // Exact numeric runs opt into primary-span triangles. BE and noncanonical
  // keys keep their full sparse identity; default/beam charts remain sparse.
  cand_pos_t dense_length() const { return dense_ ? dense_->n : 0; }
  // An exact-only observation index: a bit is set precisely when the dense
  // primary value differs from logarithmic zero. Sparse/beam charts and
  // non-log defaults allocate no rows and always request the caller's fallback.
  cand_pos_t live_primary_index_length() const {
    return dense_ && !dense_->live_starts.empty() ? dense_->n : 0;
  }
  std::size_t live_primary_index_bytes() const {
    return dense_ ? (dense_->live_starts.capacity() + dense_->live_ends.capacity()) *
                        sizeof(std::uint64_t) : 0;
  }

  template <class Fn>
  bool for_each_live_primary_start_union(
      grammar::NonTerminal first_family, cand_pos_t first_end,
      grammar::NonTerminal second_family, cand_pos_t second_end,
      cand_pos_t first, cand_pos_t last, Fn &&fn) const {
    const auto *left = live_primary_row(first_family, first_end, false);
    const auto *right = live_primary_row(second_family, second_end, false);
    if (!left || !right || first < 1 || last > dense_->n) return false;
    if (first > last) return true;
    for_each_live_bit(left, right, first, last, std::forward<Fn>(fn));
    return true;
  }

  template <class Fn>
  bool for_each_live_primary_end(grammar::NonTerminal family, cand_pos_t start,
                                 cand_pos_t first, cand_pos_t last, Fn &&fn) const {
    const auto *row = live_primary_row(family, start, true);
    if (!row || first < 1 || last > dense_->n) return false;
    if (first > last) return true;
    for_each_live_bit(row, row, first, last, std::forward<Fn>(fn));
    return true;
  }
  cand_pos_t primary_lookup_length() const { return primary_lookup_ ? primary_lookup_->length() : 0; }
  std::size_t primary_lookup_bytes() const { return primary_lookup_ ? primary_lookup_->memory_bytes() : 0; }
  // Optional numeric point-mirror compaction after a complete selection.
  // Generic values and charts without a mirror retain their existing path.
  void seal_primary_lookup_row(grammar::NonTerminal family, cand_pos_t end) {
    if constexpr (std::is_same_v<Value, double>) {
      if (primary_lookup_) primary_lookup_->seal_row(family, end);
    }
  }
  std::size_t size() const { return values_.size() + (dense_ ? dense_->count : 0); }

  Value get(const key_type &key) const {
    if (const auto index = dense_index(key)) {
      // Missing cells are initialized/reset to default_value_. Presence is
      // needed for contains()/iteration, but reading it adds no information.
      if constexpr (std::is_same_v<Value, double>) return dense_->values[*index];
      else return dense_->present[*index] ? dense_->values[*index] : default_value_;
    }
    if constexpr (std::is_same_v<Value, double>) {
      if (primary_lookup_)
        if (const auto value = primary_lookup_->get(key, default_value_)) return *value;
    }
    const auto it = values_.find(key);
    if (it == values_.end()) {
      return default_value_;
    }
    return it->second;
  }

  void set(const key_type &key, Value value) {
    const auto index = dense_index(key);
    if (semirings::ValueBoundary<Value>::discard(value) ||
        (discard_ && discard_(value))) {
      if (index) {
        if (dense_->present[*index]) --dense_->count;
        dense_->present[*index] = false;
        dense_->values[*index] = default_value_;
        dense_->set_live(key, default_value_);
      } else {
        values_.erase(key);
      }
      update_primary_lookup(key, default_value_, false);
      if (materialized_) materialized_->erase(key);
    } else {
      // A representation boundary must run once, including when the legacy
      // map view has been requested (factorized values can acquire identity).
      value = semirings::checkpoint_value(value);
      if (index) {
        if (!dense_->present[*index]) ++dense_->count;
        dense_->present[*index] = true;
        dense_->values[*index] = value;
        dense_->set_live(key, value);
      } else {
        values_[key] = value;
      }
      update_primary_lookup(key, value, true);
      if (materialized_) (*materialized_)[key] = value;
    }
  }

  bool contains(const key_type &key) const {
    if (const auto index = dense_index(key)) return dense_->present[*index];
    return values_.find(key) != values_.end();
  }

  template <class Fn>
  void for_each_entry(Fn &&fn) const {
    if (dense_) {
      std::size_t index = 0;
      for (std::size_t family = 0; family < primary_families; ++family) {
        for (cand_pos_t j = 1; j <= dense_->n; ++j) {
          for (cand_pos_t i = 1; i <= j; ++i, ++index) {
            if (dense_->present[index]) {
              fn(key_type{static_cast<grammar::NonTerminal>(family), i, j, -1, -1},
                 dense_->values[index]);
            }
          }
        }
      }
    }
    for (const auto &[key, value] : values_) fn(key, value);
  }

  // Compatibility/diagnostic view. Numeric consumers use for_each_entry to
  // avoid duplicating dense storage. A dense view follows set() mutations;
  // references to that optional view expire when the chart is reassigned.
  const map_type &values() const {
    if (!dense_) return values_;
    if (!materialized_) {
      materialized_.emplace();
      materialized_->reserve(size());
      for_each_entry([&](const key_type &key, const Value &value) {
        materialized_->emplace(key, value);
      });
    }
    return *materialized_;
  }

 private:
  // The canonical map retains its mutation and iteration order. This optional
  // numeric index only mirrors checkpointed values. If its allocation fails,
  // disable it before any subsequent lookup can observe stale data.
  void update_primary_lookup(const key_type &key, const Value &value, bool present) {
    if constexpr (std::is_same_v<Value, double>) {
      if (!primary_lookup_) return;
      try {
        if (present) primary_lookup_->set(key, value);
        else primary_lookup_->erase(key);
      } catch (const std::bad_alloc &) { primary_lookup_.reset(); }
        catch (const std::length_error &) { primary_lookup_.reset(); }
    }
  }

  void copy_primary_lookup(const Chart &other) {
    if (!other.primary_lookup_) return;
    try { primary_lookup_.emplace(*other.primary_lookup_); }
    catch (const std::bad_alloc &) {}
    catch (const std::length_error &) {}
  }

  static constexpr std::size_t primary_families =
      static_cast<std::size_t>(grammar::NonTerminal::BE);
  struct DensePrimary {
    cand_pos_t n;
    std::size_t stride;
    std::vector<Value> values;
    std::vector<unsigned char> present;
    std::size_t live_words = 0;
    std::vector<std::uint64_t> live_starts;
    std::vector<std::uint64_t> live_ends;
    std::size_t count = 0;
    DensePrimary(cand_pos_t length, const Value &zero) : n(length), stride(0) {
      const auto extent = static_cast<std::size_t>(length);
      if (extent > std::numeric_limits<std::size_t>::max() / (extent + 1))
        throw std::length_error("dense chart dimensions overflow");
      stride = extent * (extent + 1) / 2;
      if (stride > values.max_size() / primary_families ||
          stride > present.max_size() / primary_families)
        throw std::length_error("dense chart dimensions overflow");
      values.assign(primary_families * stride, zero);
      present.assign(primary_families * stride, false);
      if constexpr (std::is_same_v<Value, double>) {
        if (zero == -std::numeric_limits<double>::infinity()) {
          live_words = (extent + 63) / 64;
          if (extent > live_starts.max_size() / primary_families / live_words) {
            live_words = 0;
            return;
          }
          const auto cells = primary_families * extent * live_words;
          try {
            live_starts.assign(cells, 0);
            live_ends.assign(cells, 0);
          } catch (const std::bad_alloc &) {
            disable_live_index();
          } catch (const std::length_error &) {
            disable_live_index();
          }
        }
      }
    }

    DensePrimary(const DensePrimary &other)
        : n(other.n), stride(other.stride), values(other.values),
          present(other.present), count(other.count) {
      // The dense cells/presence remain authoritative. An optional support
      // copy must not make an otherwise successful chart copy fail.
      if (!other.live_starts.empty()) {
        try {
          live_starts = other.live_starts;
          live_ends = other.live_ends;
          live_words = other.live_words;
        } catch (const std::bad_alloc &) {
          disable_live_index();
        } catch (const std::length_error &) {
          disable_live_index();
        }
      }
    }
    DensePrimary &operator=(const DensePrimary &other) {
      if (this == &other) return *this;
      DensePrimary copy(other);
      *this = std::move(copy);
      return *this;
    }
    DensePrimary(DensePrimary &&) = default;
    DensePrimary &operator=(DensePrimary &&) = default;

    void disable_live_index() {
      std::vector<std::uint64_t>().swap(live_starts);
      std::vector<std::uint64_t>().swap(live_ends);
      live_words = 0;
    }

    void set_live(const key_type &key, const Value &value) {
      if constexpr (std::is_same_v<Value, double>) {
        if (live_starts.empty()) return;
        const auto family = static_cast<std::size_t>(key.nonterminal);
        const auto begin = static_cast<std::size_t>(key.i - 1);
        const auto end = static_cast<std::size_t>(key.j - 1);
        auto &start_word = live_starts[(family * n + end) * live_words + begin / 64];
        auto &end_word = live_ends[(family * n + begin) * live_words + end / 64];
        const std::uint64_t start_bit = std::uint64_t{1} << (begin % 64);
        const std::uint64_t end_bit = std::uint64_t{1} << (end % 64);
        // NaN, positive infinity and every finite number remain live, exactly
        // as LogSumExp/LogMax::is_zero; this is never a finite cutoff.
        if (value != -std::numeric_limits<double>::infinity()) {
          start_word |= start_bit;
          end_word |= end_bit;
        } else {
          start_word &= ~start_bit;
          end_word &= ~end_bit;
        }
      }
    }
  };

  const std::uint64_t *live_primary_row(grammar::NonTerminal nt,
                                       cand_pos_t fixed, bool end_row) const {
    const auto family = static_cast<std::size_t>(nt);
    if (!dense_ || dense_->live_starts.empty() || family >= primary_families ||
        fixed < 1 || fixed > dense_->n) return nullptr;
    const auto offset = (family * static_cast<std::size_t>(dense_->n) + fixed - 1) *
                        dense_->live_words;
    return (end_row ? dense_->live_ends : dense_->live_starts).data() + offset;
  }

  template <class Fn>
  static void for_each_live_bit(const std::uint64_t *left,
                                const std::uint64_t *right,
                                cand_pos_t first, cand_pos_t last, Fn &&fn) {
    const auto first_bit = static_cast<std::size_t>(first - 1);
    const auto last_bit = static_cast<std::size_t>(last - 1);
    const auto first_word = first_bit / 64;
    const auto last_word = last_bit / 64;
    for (auto word = first_word; word <= last_word; ++word) {
      auto bits = left[word] | right[word];
      if (word == first_word) bits &= ~std::uint64_t{0} << (first_bit % 64);
      if (word == last_word) bits &= ~std::uint64_t{0} >> (63 - last_bit % 64);
      while (bits) {
        const auto bit = static_cast<unsigned>(__builtin_ctzll(bits));
        fn(static_cast<cand_pos_t>(word * 64 + bit + 1));
        bits &= bits - 1;
      }
    }
  }

  std::optional<std::size_t> dense_index(const key_type &key) const {
    const auto family = static_cast<std::size_t>(key.nonterminal);
    if (!dense_ || family >= primary_families || key.ip != -1 || key.jp != -1 ||
        key.i < 1 || key.j < key.i || key.j > dense_->n) return std::nullopt;
    const auto end = static_cast<std::size_t>(key.j);
    return family * dense_->stride + end * (end - 1) / 2 + key.i - 1;
  }

  Value default_value_;
  discard_fn discard_;
  map_type values_;
  std::optional<DensePrimary> dense_;
  mutable std::optional<map_type> materialized_;
  std::optional<PrimaryLookupIndex> primary_lookup_;
};

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_CHART_HH

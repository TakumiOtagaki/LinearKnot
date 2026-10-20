#ifndef CPARTY_SCFG2_ENGINE_CHART_HH
#define CPARTY_SCFG2_ENGINE_CHART_HH

#include "scfg2/deductions/item_key.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <unordered_map>
#include <limits>
#include <optional>
#include <stdexcept>
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
                 cand_pos_t dense_length = 0)
      : default_value_(default_value), discard_(discard) {
    if (dense_length < 0) throw std::invalid_argument("negative dense chart length");
    if (dense_length) dense_.emplace(dense_length, default_value);
  }

  Chart(const Chart &) = default;
  Chart &operator=(const Chart &) = default;
  Chart(Chart &&other)
      : default_value_(other.default_value_), discard_(other.discard_),
        values_(std::move(other.values_)), dense_(std::move(other.dense_)),
        materialized_(std::move(other.materialized_)) {
    other.dense_.reset();
    other.materialized_.reset();
  }
  Chart &operator=(Chart &&other) {
    if (this == &other) return *this;
    default_value_ = other.default_value_;
    discard_ = other.discard_;
    values_ = std::move(other.values_);
    dense_ = std::move(other.dense_);
    materialized_ = std::move(other.materialized_);
    other.dense_.reset();
    other.materialized_.reset();
    return *this;
  }

  // Exact numeric runs opt into primary-span triangles. BE and noncanonical
  // keys keep their full sparse identity; default/beam charts remain sparse.
  cand_pos_t dense_length() const { return dense_ ? dense_->n : 0; }
  std::size_t size() const { return values_.size() + (dense_ ? dense_->count : 0); }

  Value get(const key_type &key) const {
    if (const auto index = dense_index(key)) {
      return dense_->present[*index] ? dense_->values[*index] : default_value_;
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
      } else {
        values_.erase(key);
      }
      if (materialized_) materialized_->erase(key);
    } else {
      // A representation boundary must run once, including when the legacy
      // map view has been requested (factorized values can acquire identity).
      value = semirings::checkpoint_value(value);
      if (index) {
        if (!dense_->present[*index]) ++dense_->count;
        dense_->present[*index] = true;
        dense_->values[*index] = value;
      } else {
        values_[key] = value;
      }
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
  static constexpr std::size_t primary_families =
      static_cast<std::size_t>(grammar::NonTerminal::BE);
  struct DensePrimary {
    cand_pos_t n;
    std::size_t stride;
    std::vector<Value> values;
    std::vector<unsigned char> present;
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
    }
  };

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
};

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_CHART_HH

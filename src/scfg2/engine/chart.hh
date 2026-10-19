#ifndef CPARTY_SCFG2_ENGINE_CHART_HH
#define CPARTY_SCFG2_ENGINE_CHART_HH

#include "scfg2/deductions/item_key.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <unordered_map>

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
  explicit Chart(Value default_value, discard_fn discard = nullptr)
      : default_value_(default_value), discard_(discard) {}

  Value get(const key_type &key) const {
    const auto it = values_.find(key);
    if (it == values_.end()) {
      return default_value_;
    }
    return it->second;
  }

  void set(const key_type &key, Value value) {
    if (semirings::ValueBoundary<Value>::discard(value) ||
        (discard_ && discard_(value))) {
      values_.erase(key);
    } else {
      values_[key] = semirings::checkpoint_value(value);
    }
  }

  bool contains(const key_type &key) const {
    return values_.find(key) != values_.end();
  }

  const map_type &values() const { return values_; }

 private:
  Value default_value_;
  discard_fn discard_;
  map_type values_;
};

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_CHART_HH

#ifndef CPARTY_SCFG2_SEMIRINGS_VALUE_BOUNDARY_HH
#define CPARTY_SCFG2_SEMIRINGS_VALUE_BOUNDARY_HH

#include <type_traits>
#include <utility>

namespace cparty::scfg2::semirings {

// A representation boundary, with no effect on numeric semirings or the
// arithmetic-tape oracle. Compact outside values seal an item/private sum
// here so that later consumers share its identity rather than its summands.
template <class Value, class = void>
struct ValueBoundary {
  static Value checkpoint(Value value) { return value; }
  static bool discard(const Value &) { return false; }
};

template <class Value>
struct ValueBoundary<Value, std::void_t<decltype(
    std::declval<const Value &>().checkpointed())>> {
  static Value checkpoint(const Value &value) { return value.checkpointed(); }
  static bool discard(const Value &value) { return value.discard_chart_entry(); }
};

template <class Value>
Value checkpoint_value(const Value &value) {
  return ValueBoundary<Value>::checkpoint(value);
}

}  // namespace cparty::scfg2::semirings

#endif

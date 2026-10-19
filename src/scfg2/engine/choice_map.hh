#ifndef CPARTY_SCFG2_ENGINE_CHOICE_MAP_HH
#define CPARTY_SCFG2_ENGINE_CHOICE_MAP_HH

#include "scfg2/deductions/deduction.hh"

#include <unordered_map>

namespace cparty::scfg2::engine {

template <class Value>
struct ChoiceRecord {
  deductions::Deduction deduction{};
  Value contribution{};
  bool has_value = false;
};

template <class Value>
class ChoiceMap {
 public:
  using key_type = deductions::ItemKey;
  using map_type = std::unordered_map<key_type,
                                      ChoiceRecord<Value>,
                                      deductions::ItemKeyHash,
                                      deductions::ItemKeyEq>;

  template <class BetterFn>
  void observe(const key_type &parent,
               const deductions::Deduction &deduction,
               Value contribution,
               BetterFn &&better) {
    auto &slot = values_[parent];
    if (!slot.has_value || better(contribution, slot.contribution)) {
      slot.deduction = deduction;
      slot.contribution = contribution;
      slot.has_value = true;
    }
  }

  template <class BetterFn, class TieBeforeFn>
  void observe_with_tie_break(const key_type &parent,
                              const deductions::Deduction &deduction,
                              Value contribution,
                              BetterFn &&better,
                              TieBeforeFn &&tie_before) {
    auto &slot = values_[parent];
    const bool tied = slot.has_value &&
        !better(contribution, slot.contribution) &&
        !better(slot.contribution, contribution);
    if (!slot.has_value || better(contribution, slot.contribution) ||
        (tied && tie_before(deduction, slot.deduction))) {
      slot.deduction = deduction;
      slot.contribution = contribution;
      slot.has_value = true;
    }
  }

  const ChoiceRecord<Value> *find(const key_type &key) const {
    const auto it = values_.find(key);
    return it == values_.end() ? nullptr : &it->second;
  }

  const map_type &values() const { return values_; }
  void erase(const key_type &key) { values_.erase(key); }

 private:
  map_type values_;
};

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_CHOICE_MAP_HH

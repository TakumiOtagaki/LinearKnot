#ifndef CPARTY_SCFG2_ENGINE_COMPACT_CHOICE_MAP_HH
#define CPARTY_SCFG2_ENGINE_COMPACT_CHOICE_MAP_HH

#include "scfg2/deductions/expand_basic.hh"
#include "scfg2/engine/choice_map.hh"

#include <optional>
#include <unordered_map>

namespace cparty::scfg2::engine {

// Store a selected deduction's descriptor, not a second copy of its parent
// and all three children. Decode through the shared grammar expansion.
// Runtime/private deductions which cannot round-trip exactly retain their
// full representation in the exception table; even inactive child slots are
// preserved because the historical tie key includes them.
template <class Value>
class CompactChoiceMap {
 public:
  using key_type = deductions::ItemKey;
  struct Record {
    Value contribution{};
    grammar::RuleSplit split{};
    grammar::RuleId rule{};
    bool fallback = false;
  };
  using map_type = std::unordered_map<key_type, Record,
      deductions::ItemKeyHash, deductions::ItemKeyEq>;
  using fallback_type = std::unordered_map<key_type, deductions::Deduction,
      deductions::ItemKeyHash, deductions::ItemKeyEq>;

  template <class BetterFn, class TieBeforeFn>
  void observe_with_tie_break(const key_type &parent,
                              const deductions::Deduction &deduction,
                              Value contribution, BetterFn &&better,
                              TieBeforeFn &&tie_before) {
    auto [found, inserted] = values_.try_emplace(parent);
    auto &slot = found->second;
    if (inserted || better(contribution, slot.contribution) ||
        (!better(slot.contribution, contribution) &&
         tie_before(deduction, decode(parent, slot)))) {
      const auto expanded = canonical(parent, deduction.rule, deduction.split);
      const bool fallback = !same_deduction(expanded, deduction);
      if (fallback) fallbacks_[parent] = deduction;
      else if (slot.fallback) fallbacks_.erase(parent);
      slot = Record{contribution, deduction.split, deduction.rule, fallback};
    }
  }

  // Return an owning value: decoding never invalidates an earlier read and
  // does not maintain a growing decoded-record cache.
  std::optional<ChoiceRecord<Value>> find(const key_type &key) const {
    const auto found = values_.find(key);
    if (found == values_.end()) return std::nullopt;
    return ChoiceRecord<Value>{decode(key, found->second),
                               found->second.contribution, true};
  }

  void erase(const key_type &key) {
    const auto found = values_.find(key);
    if (found == values_.end()) return;
    if (found->second.fallback) fallbacks_.erase(key);
    values_.erase(found);
  }

  const map_type &values() const { return values_; }
  const fallback_type &fallback_values() const { return fallbacks_; }

 private:
  static deductions::Deduction canonical(const key_type &parent,
      grammar::RuleId rule, const grammar::RuleSplit &split) {
    return deductions::expand_basic_rule(parent, rule, split);
  }

  static bool same_deduction(const deductions::Deduction &a,
                            const deductions::Deduction &b) {
    const deductions::ItemKeyEq equal;
    if (!equal(a.parent, b.parent) || a.child_count != b.child_count) return false;
    for (std::size_t c = 0; c < a.children.size(); ++c)
      if (!equal(a.children[c].item, b.children[c].item)) return false;
    return true;
  }

  deductions::Deduction decode(const key_type &key, const Record &record) const {
    if (record.fallback) return fallbacks_.at(key);
    return canonical(key, record.rule, record.split);
  }

  map_type values_;
  fallback_type fallbacks_;
};

}  // namespace cparty::scfg2::engine
#endif

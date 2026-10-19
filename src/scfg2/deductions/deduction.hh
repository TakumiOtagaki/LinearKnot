#ifndef CPARTY_SCFG2_DEDUCTIONS_DEDUCTION_HH
#define CPARTY_SCFG2_DEDUCTIONS_DEDUCTION_HH

#include "scfg2/deductions/item_key.hh"
#include "scfg2/grammar/rules.hh"

#include <array>
#include <cstddef>

namespace cparty::scfg2::deductions {

struct ChildRef {
  ItemKey item;
};

struct Deduction {
  ItemKey parent;
  grammar::RuleId rule = grammar::RuleId::W_EXTEND_UNPAIRED;
  grammar::RuleSplit split{};
  std::array<ChildRef, 3> children{};
  std::size_t child_count = 0;
};

inline constexpr bool is_valid(const ChildRef &child) {
  return is_valid_item_key(child.item);
}

inline constexpr bool is_valid(const Deduction &deduction) {
  if (!is_valid_item_key(deduction.parent)) {
    return false;
  }
  if (deduction.child_count > deduction.children.size()) {
    return false;
  }
  for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
    if (!is_valid(deduction.children[idx])) {
      return false;
    }
  }
  return true;
}

}  // namespace cparty::scfg2::deductions

#endif  // CPARTY_SCFG2_DEDUCTIONS_DEDUCTION_HH

#ifndef CPARTY_SCFG2_DEDUCTIONS_GENERATE_BASIC_HH
#define CPARTY_SCFG2_DEDUCTIONS_GENERATE_BASIC_HH

#include "scfg2/deductions/expand_basic.hh"

#include <vector>

namespace cparty::scfg2::deductions {

inline std::vector<Deduction> generate_basic_deductions(
    const ItemKey &parent,
    grammar::RuleId rule,
    const std::vector<grammar::RuleSplit> &splits) {
  std::vector<Deduction> out;
  out.reserve(splits.size());
  for (const auto &split : splits) {
    Deduction deduction = expand_basic_rule(parent, rule, split);
    if (is_valid(deduction)) {
      out.push_back(deduction);
    }
  }
  return out;
}

inline std::vector<Deduction> generate_basic_deductions_single(
    const ItemKey &parent,
    grammar::RuleId rule) {
  return generate_basic_deductions(parent, rule, std::vector<grammar::RuleSplit>{{}});
}

}  // namespace cparty::scfg2::deductions

#endif  // CPARTY_SCFG2_DEDUCTIONS_GENERATE_BASIC_HH

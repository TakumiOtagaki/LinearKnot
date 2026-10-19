#ifndef CPARTY_SCFG2_DEDUCTIONS_BASIC_SPLIT_ENUMERATOR_HH
#define CPARTY_SCFG2_DEDUCTIONS_BASIC_SPLIT_ENUMERATOR_HH

#include "scfg2/deductions/item_key.hh"
#include "scfg2/grammar/rules.hh"

#include <vector>

namespace cparty::scfg2::deductions {

inline std::vector<grammar::RuleSplit> enumerate_basic_splits(
    const ItemKey &parent,
    grammar::RuleId rule) {
  std::vector<grammar::RuleSplit> out;
  switch (rule) {
    case grammar::RuleId::W_SPLIT_V:
      for (cand_pos_t k = parent.i; k <= parent.j; ++k) {
        grammar::RuleSplit split{};
        split.k = k;
        out.push_back(split);
      }
      return out;
    case grammar::RuleId::W_EXTEND_UNPAIRED:
    case grammar::RuleId::V_HAIRPIN:
      return {grammar::RuleSplit{}};
    default:
      return out;
  }
}

}  // namespace cparty::scfg2::deductions

#endif  // CPARTY_SCFG2_DEDUCTIONS_BASIC_SPLIT_ENUMERATOR_HH

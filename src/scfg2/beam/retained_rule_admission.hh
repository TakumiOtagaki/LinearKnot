#ifndef CPARTY_SCFG2_BEAM_RETAINED_RULE_ADMISSION_HH
#define CPARTY_SCFG2_BEAM_RETAINED_RULE_ADMISSION_HH

#include "scfg2/deductions/deduction.hh"

#include <cstdint>
#include <unordered_map>

namespace cparty::scfg2::beam {

// Records which rule families actually contributed, rather than their split
// lists or products. Candidate observation precedes public pruning. On the
// next family/column the previous barrier has committed, so only surviving
// parents are copied from scratch into the immutable admission table.
// This does not distinguish splits within one rule: such families need a
// retained-child generator, temporal cutoff, or explicit private-cell view.
template <class Semiring>
class RetainedRuleAdmission {
 public:
  template <class Chart, class LocalWeight>
  void operator()(const deductions::ItemKey &parent,
                   const deductions::Deduction &source,
                   const typename Semiring::value_type &candidate,
                   const Chart &chart, LocalWeight) {
    if (!has_column_ || parent.nonterminal != family_ || parent.j != column_) {
      finish(chart);
      family_ = parent.nonterminal;
      column_ = parent.j;
      has_column_ = true;
    }
    if (!Semiring::is_zero(candidate)) pending_[parent] |= bit(source.rule);
  }

  template <class Chart>
  void finish(const Chart &chart) {
    for (const auto &[parent, rules] : pending_)
      if (!Semiring::is_zero(chart.get(parent))) admitted_[parent] = rules;
    pending_.clear();
    has_column_ = false;
  }

  std::uint64_t rules_for(const deductions::ItemKey &parent) const {
    const auto found = admitted_.find(parent);
    return found == admitted_.end() ? 0 : found->second;
  }
  static bool has(std::uint64_t rules, grammar::RuleId rule) {
    return rules & bit(rule);
  }
  bool has(const deductions::ItemKey &parent, grammar::RuleId rule) const {
    return has(rules_for(parent), rule);
  }
  std::size_t size() const { return admitted_.size(); }
  std::size_t pending_size() const { return pending_.size(); }

 private:
  static constexpr std::uint64_t bit(grammar::RuleId rule) {
    static_assert(grammar::all_rule_ids().size() <= 64);
    return std::uint64_t{1} << static_cast<unsigned>(rule);
  }
  using Map = std::unordered_map<deductions::ItemKey, std::uint64_t,
      deductions::ItemKeyHash, deductions::ItemKeyEq>;
  Map pending_, admitted_;
  grammar::NonTerminal family_{};
  cand_pos_t column_ = -1;
  bool has_column_ = false;
};

}  // namespace cparty::scfg2::beam
#endif

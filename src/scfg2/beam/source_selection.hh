#ifndef CPARTY_SCFG2_BEAM_SOURCE_SELECTION_HH
#define CPARTY_SCFG2_BEAM_SOURCE_SELECTION_HH

#include "scfg2/deductions/deduction.hh"
#include "scfg2/engine/traceback.hh"

#include <vector>

namespace cparty::scfg2::beam {

inline bool same_source_deduction(const deductions::Deduction &lhs,
                                  const deductions::Deduction &rhs) {
  return deductions::ItemKeyEq{}(lhs.parent, rhs.parent) &&
         lhs.rule == rhs.rule && lhs.split.k == rhs.split.k &&
         lhs.split.l == rhs.split.l && lhs.split.p == rhs.split.p &&
         lhs.split.q == rhs.split.q;
}

template <class Value>
void collect_selected_source_deductions(
    const engine::DeductionTraceNode<Value> &node,
    std::vector<deductions::Deduction> *out) {
  if (!out) return;
  std::vector<const engine::DeductionTraceNode<Value> *> pending{&node};
  while (!pending.empty()) {
    const auto *current = pending.back();
    pending.pop_back();
    if (current->chosen_deduction.has_value())
      out->push_back(*current->chosen_deduction);
    for (auto child = current->children.rbegin();
         child != current->children.rend(); ++child)
      pending.push_back(&*child);
  }
}

inline bool source_deduction_is_selected(
    const deductions::Deduction &deduction,
    const std::vector<deductions::Deduction> &selected) {
  for (const auto &candidate : selected) {
    if (same_source_deduction(deduction, candidate)) return true;
  }
  return false;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_SOURCE_SELECTION_HH

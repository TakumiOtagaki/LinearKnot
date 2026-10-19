#ifndef CPARTY_SCFG2_DEDUCTIONS_RETAINED_START_SPLITS_HH
#define CPARTY_SCFG2_DEDUCTIONS_RETAINED_START_SPLITS_HH

#include "scfg2/deductions/generate_exact_basic.hh"

namespace cparty::scfg2::deductions {

// Start aggregates only select the parent. Once that parent is admitted,
// reify its original deductions over the retained seed positions, with no
// interval-length scan and no saved private prefix branch list.
template <class Emit, class StartsFn>
void for_each_retained_start_split(
    const ItemKey &parent, const PairabilityContext &ctx,
    StartsFn &&starts, Emit &&emit) {
  using NT = grammar::NonTerminal;
  using R = grammar::RuleId;
  const auto append = [&](R rule, cand_pos_t k) {
    grammar::RuleSplit split{};
    split.k = k;
    const auto d = expand_basic_rule(parent, rule, split);
    if (is_valid(d)) emit(d);
  };
  if (parent.nonterminal == NT::VM) {
    for (auto k : starts(NT::WMp, parent.j - 1))
      if (k >= exact_vm_base_begin_k(parent) && k <= exact_vm_end_k(parent, ctx) &&
          exact_vm_base_split_allows(parent, k, ctx)) append(R::VM_SPLIT_WMp_BASE, k);
    return;
  }
  if (parent.nonterminal == NT::VPL) {
    for (auto k : starts(NT::VP, parent.j))
      if (exact_vpl_split_allows(parent, k, ctx)) append(R::VPL_SPLIT_VP, k);
    return;
  }
  // The forward prefix builds the union of seed positions first, then reads
  // both families at every seed. A root-only WMB value can therefore enter
  // here when V supplied the same seed position. Preserve that boundary.
  std::vector<cand_pos_t> seeds = starts(NT::V, parent.j);
  const auto &wmb = starts(NT::WMB, parent.j);
  seeds.insert(seeds.end(), wmb.begin(), wmb.end());
  std::sort(seeds.begin(), seeds.end());
  seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());
  for (auto k : seeds) {
    if (parent.nonterminal == NT::WM) {
      if (k >= parent.i && k <= parent.j && exact_wm_start_split_allows(parent, k, ctx)) {
        append(R::WM_START_V, k);
        append(R::WM_START_WMB, k);
      }
    } else if (parent.nonterminal == NT::WIP) {
      if (k == parent.i) {
        append(R::WIP_BASE_V, k);
        append(R::WIP_BASE_WMB, k);
      } else if (exact_wip_start_split_allows(parent, k, ctx)) {
        append(R::WIP_BASEPAIR_V, k);
        append(R::WIP_BASEPAIR_WMB, k);
      }
    }
  }
}

}  // namespace cparty::scfg2::deductions
#endif

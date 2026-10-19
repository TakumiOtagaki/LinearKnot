#ifndef CPARTY_SCFG2_ENGINE_RUNTIME_CHILD_ITEM_HH
#define CPARTY_SCFG2_ENGINE_RUNTIME_CHILD_ITEM_HH

#include "scfg2/deductions/deduction.hh"
#include <cstddef>

namespace cparty::scfg2::engine {

// Shared child lookup; the grammar directly names every consumed state.
struct RuntimeDeductionChildItem {
  deductions::ItemKey operator()(const deductions::Deduction &source,
                                  std::size_t c) const {
    auto child = source.children[c].item;
    return child;
  }
};

}  // namespace cparty::scfg2::engine
#endif

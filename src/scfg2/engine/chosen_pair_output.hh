#ifndef CPARTY_SCFG2_ENGINE_CHOSEN_PAIR_OUTPUT_HH
#define CPARTY_SCFG2_ENGINE_CHOSEN_PAIR_OUTPUT_HH

#include "scfg2/engine/choice_map.hh"
#include "scfg2/engine/traceback_pairs.hh"

#include <unordered_set>
#include <vector>

namespace cparty::scfg2::engine {

struct ChosenPairOutput {
  std::vector<TracebackPair> pairs;
  std::vector<TracebackPair> generated_pairs;
  std::size_t visited_items = 0;
  std::size_t repeated_items = 0;
  std::size_t child_occurrences = 0;
};

// Observation-only pair consumer of immutable selected choices. A second
// occurrence of an item has exactly the same descendants and emissions, so
// it cannot add anything to the historical stable-unique pair output. DFS
// first visits preserve that output's order without constructing a tree.
// This does not change the full occurrence-tree contract in traceback.hh.
template <class Choices, class Context>
ChosenPairOutput collect_chosen_pair_output(
    const deductions::ItemKey &root, const Choices &choices,
    const Context &ctx) {
  ChosenPairOutput out;
  std::unordered_set<deductions::ItemKey, deductions::ItemKeyHash,
                     deductions::ItemKeyEq> visited;
  std::vector<deductions::ItemKey> pending{root};
  while (!pending.empty()) {
    const auto item = pending.back();
    pending.pop_back();
    if (!visited.insert(item).second) {
      ++out.repeated_items;
      continue;
    }
    ++out.visited_items;
    const auto choice = choices.find(item);
    if (!choice || !choice->has_value) continue;
    // Append helpers deduplicate their argument with a scan: keep that list
    // local and fixed-size, rather than passing the growing global output.
    std::vector<TracebackPair> local;
    append_deduction_output_pairs(choice->deduction, ctx, &local);
    out.pairs.insert(out.pairs.end(), local.begin(), local.end());
    local.clear();
    append_deduction_generated_pairs(choice->deduction, ctx, &local);
    out.generated_pairs.insert(out.generated_pairs.end(), local.begin(), local.end());
    out.child_occurrences += choice->deduction.child_count;
    for (std::size_t c = choice->deduction.child_count; c > 0; --c)
      pending.push_back(choice->deduction.children[c - 1].item);
  }
  stable_unique_traceback_pairs(&out.pairs);
  stable_unique_traceback_pairs(&out.generated_pairs);
  return out;
}

}  // namespace cparty::scfg2::engine

#endif

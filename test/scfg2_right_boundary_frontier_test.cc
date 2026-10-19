#include "scfg2/beam/right_boundary_item_frontier.hh"

#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
  using namespace cparty::scfg2;
  using namespace beam;
  std::size_t comparisons = 0;
  try {
    for (std::size_t width : {0, 1, 2, 20, 65}) for (bool forced : {false, true}) {
      RightBoundaryTopWidthItemFrontier frontier(width);
      constexpr cand_pos_t count = 10000, j = count + 1;
      const deductions::ItemKey force{grammar::NonTerminal::WI, count, j, -1, -1};
      auto keys = std::make_shared<ItemKeySet>();
      if (forced) keys->insert(force);
      frontier.set_forced(keys);
      std::vector<RankedRightBoundaryItem> expected, mandatory;
      for (cand_pos_t i = 1; i <= count; ++i) {
        const deductions::ItemKey item{grammar::NonTerminal::WI, i, j, -1, -1};
        double value = (i * 37) % 31;
        frontier.observe(item, value);
        if (i % 7 == 0) { value += 100.; frontier.observe(item, value); }
        (keys->count(item) ? mandatory : expected).push_back({item, value});
      }
      std::sort(expected.begin(), expected.end(), right_boundary_item_rank_before);
      expected.resize(std::min(width, expected.size()));
      expected.insert(expected.end(), mandatory.begin(), mandatory.end());
      std::sort(expected.begin(), expected.end(), right_boundary_item_rank_before);
      std::vector<deductions::ItemKey> actual, pruned;
      frontier.flush_column_partition(j, [&](const auto &key) { actual.push_back(key); },
                                      [&](const auto &key) { pruned.push_back(key); });
      if (actual.size() != expected.size() || actual.size() + pruned.size() != count ||
          frontier.candidate_count() != count || frontier.selected_count() != actual.size() ||
          frontier.pruned_count() != pruned.size() || !frontier.respects_width())
        throw std::runtime_error("frontier partition or forced-width contract changed");
      for (std::size_t i = 0; i < actual.size(); ++i)
        if (!deductions::ItemKeyEq{}(actual[i], expected[i].item))
          throw std::runtime_error("heap selection differs from independent full-sort oracle");
      const auto budget = 32 * count * static_cast<std::size_t>(std::ceil(std::log2(width + 2)));
      if (frontier.selection_comparisons() > budget)
        throw std::runtime_error("frontier comparison budget exceeded");
      comparisons += frontier.selection_comparisons();
    }
    std::cout << "scfg2_right_boundary_frontier=ok cases=10 candidates=100000 comparisons="
              << comparisons << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  return 0;
}

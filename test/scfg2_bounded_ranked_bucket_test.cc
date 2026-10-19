#include "scfg2/beam/bounded_ranked_bucket.hh"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

namespace {
using cparty::scfg2::beam::BoundedRankedBucket;
struct Cell { int score, coordinate, identity; };
struct Rank {
  bool operator()(const Cell &a, const Cell &b) const {
    return a.score != b.score ? a.score > b.score : a.coordinate < b.coordinate;
  }
};
void require(bool ok, const char *message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
void compare(const std::vector<Cell> &actual, const std::vector<Cell> &expected) {
  require(actual.size() == expected.size(), "retained count differs from sort oracle");
  for (std::size_t x = 0; x < actual.size(); ++x)
    require(actual[x].score == expected[x].score &&
            actual[x].coordinate == expected[x].coordinate &&
            actual[x].identity == expected[x].identity, "rank/order differs from stable sort oracle");
}
}  // namespace

int main() {
  std::mt19937 rng(20261009);
  std::size_t checked = 0;
  // Every insertion is compared against independent full sorting/truncation.
  // Repeated primary keys exercise FIFO ties at and across the beam cutoff.
  for (std::size_t width : {0u, 1u, 2u, 7u, 16u, 64u}) {
    BoundedRankedBucket<Cell, Rank> bucket;
    std::vector<Cell> reference;
    for (int x = 0; x < 3000; ++x) {
      const Cell cell{static_cast<int>(rng() % 13), static_cast<int>(rng() % 7), x};
      bucket.insert(cell, width);
      reference.push_back(cell);
      std::stable_sort(reference.begin(), reference.end(), Rank{});
      reference.resize(std::min(reference.size(), width));
      compare(bucket.ranked(), reference);
      const auto before = bucket.stats();
      compare(bucket.ranked(), reference);
      require(bucket.stats().view_builds == before.view_builds &&
              bucket.stats().view_comparisons == before.view_comparisons &&
              bucket.stats().view_cells == before.view_cells, "cached read did work again");
      ++checked;
    }
  }
  // Width transitions, including KeepAll, preserve current retained history.
  BoundedRankedBucket<Cell, Rank> changing;
  std::vector<Cell> reference;
  for (int x = 0; x < 300; ++x) {
    const std::size_t width = x < 70 ? std::numeric_limits<std::size_t>::max() :
                              x < 150 ? 3 : x < 220 ? 0 : 12;
    const Cell cell{x % 9, x % 4, x};
    changing.insert(cell, width);
    reference.push_back(cell);
    std::stable_sort(reference.begin(), reference.end(), Rank{});
    if (reference.size() > width) reference.resize(width);
    compare(changing.ranked(), reference);
    ++checked;
  }
  // b^2 improving candidates per column force replacements. There is just
  // one output phase per column; increasing b cannot hide a full sort per
  // insertion inside the counted logarithmic update budget.
  for (std::size_t width : {1u, 4u, 16u, 64u, 256u}) {
    BoundedRankedBucket<Cell, Rank> bucket;
    const std::size_t columns = 16, per_column = width * width;
    for (std::size_t col = 0; col < columns; ++col) {
      for (std::size_t x = 0; x < per_column; ++x) {
        const int score = static_cast<int>(col * per_column + x);
        bucket.insert(Cell{score, 0, score}, width);
      }
      const auto &view = bucket.ranked();
      require(view.size() == width, "scaling beam width exceeded or lost");
      for (std::size_t x = 0; x < width; ++x)
        require(view[x].score == static_cast<int>((col + 1) * per_column - 1 - x),
                "scaling retained maximum incorrect");
    }
    std::size_t lg = 0;
    for (std::size_t x = width + 2; x > 1; x >>= 1) ++lg;
    const auto &stats = bucket.stats();
    require(stats.update_comparisons <= columns * per_column * (6 * lg + 6),
            "update comparisons exceeded logarithmic budget");
    require(stats.view_builds == columns && stats.view_cells == columns * width,
            "view materialized outside output phases");
    require(stats.view_comparisons <= columns * width * (8 * lg + 8),
            "output sorting exceeded per-column budget");
    std::cout << "b=" << width << " insertions=" << columns * per_column
              << " update_comparisons=" << stats.update_comparisons
              << " view_comparisons=" << stats.view_comparisons << '\n';
  }
  std::cout << "differential_insertions=" << checked << "\nscfg2_bounded_ranked_bucket=ok\n";
}

#include "scfg2/geometry/static_range_minimum.hh"

#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

namespace {
using cparty::scfg2::geometry::PlusMinusRangeMinimum;
using cparty::scfg2::geometry::StaticRangeMinimum;
std::size_t queries = 0, max_reads = 0;

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void check_walk(const std::vector<int> &walk) {
  PlusMinusRangeMinimum index(walk);
  for (std::size_t lo = 0; lo < walk.size(); ++lo) {
    std::size_t first = lo, last = lo;
    for (std::size_t hi = lo; hi < walk.size(); ++hi) {
      if (walk[hi] < walk[first]) first = last = hi;
      else if (walk[hi] == walk[first]) last = hi;
      std::size_t first_reads = 0, last_reads = 0;
      require(index.argmin(lo, hi, PlusMinusRangeMinimum::Tie::First, &first_reads) == first,
              "leftmost +/-1 RMQ differs from scan");
      require(index.argmin(lo, hi, PlusMinusRangeMinimum::Tie::Last, &last_reads) == last,
              "rightmost +/-1 RMQ differs from scan");
      require(first_reads <= 4 && last_reads <= 4, "RMQ exceeded four table reads");
      max_reads = std::max({max_reads, first_reads, last_reads});
      ++queries;
    }
  }
}

void check_values(const std::vector<int> &values) {
  StaticRangeMinimum<int> index(values);
  for (std::size_t lo = 0; lo < values.size(); ++lo) {
    std::size_t first = lo;
    for (std::size_t hi = lo; hi < values.size(); ++hi) {
      if (values[hi] < values[first]) first = hi;
      std::size_t reads = 0;
      require(index.argmin(lo, hi, &reads) == first, "Cartesian RMQ differs from scan");
      require(reads <= 4, "Cartesian RMQ exceeded four table reads");
      max_reads = std::max(max_reads, reads);
      ++queries;
    }
  }
}
}  // namespace

int main() {
  for (std::size_t n = 1; n <= 12; ++n) {
    for (std::size_t mask = 0; mask < (std::size_t{1} << (n - 1)); ++mask) {
      std::vector<int> walk(n);
      for (std::size_t x = 1; x < n; ++x)
        walk[x] = walk[x - 1] + ((mask >> (x - 1)) & 1 ? 1 : -1);
      check_walk(walk);
    }
  }
  std::mt19937 rng(20261009);
  for (int trial = 0; trial < 60; ++trial) {
    const std::size_t n = 1 + rng() % 300;
    std::vector<int> values(n), walk(n);
    for (std::size_t x = 0; x < n; ++x) {
      values[x] = static_cast<int>(rng() % 19) - 9;
      if (x) walk[x] = walk[x - 1] + (rng() % 2 ? 1 : -1);
    }
    check_values(values);
    check_walk(walk);
  }
  // Monotone/constant Cartesian trees have depth n. They exercise iterative
  // traversal, long inputs, tie stability, and storage rather than wall time.
  for (std::size_t n : {4096u, 16384u, 65536u, 262144u}) {
    for (int shape = 0; shape < 3; ++shape) {
      std::vector<int> values(n);
      for (std::size_t x = 0; x < n; ++x)
        values[x] = shape == 0 ? static_cast<int>(x) :
                    shape == 1 ? -static_cast<int>(x) : 7;
      StaticRangeMinimum<int> index(values);
      require(index.storage_bytes() <= 256 * n, "RMQ storage exceeded linear witness bound");
      for (int q = 0; q < 2000; ++q) {
        std::size_t lo = rng() % n, hi = rng() % n;
        if (lo > hi) std::swap(lo, hi);
        std::size_t reads = 0;
        require(index.argmin(lo, hi, &reads) == (shape == 1 ? hi : lo),
                "long monotone/constant RMQ mismatch");
        require(reads <= 4, "long RMQ exceeded four table reads");
        max_reads = std::max(max_reads, reads);
        ++queries;
      }
      std::cout << "n=" << n << " shape=" << shape
                << " storage_bytes=" << index.storage_bytes() << '\n';
    }
  }
  std::cout << "queries=" << queries << " max_table_reads=" << max_reads
            << "\nscfg2_static_range_minimum=ok\n";
}

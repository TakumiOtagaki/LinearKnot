#ifndef CPARTY_SCFG2_GEOMETRY_STATIC_RANGE_MINIMUM_HH
#define CPARTY_SCFG2_GEOMETRY_STATIC_RANGE_MINIMUM_HH

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cparty::scfg2::geometry {

// Bender--Farach-Colton blocking for a +/-1 depth walk. With block width
// k = floor(log2(m))/2, the universal micro table has O(2^k k^2) = O(m)
// entries and the block-minimum sparse table has O((m/k) log(m/k)) = O(m).
// Both construction and storage are linear; an inclusive query reads at most
// four table entries. Keep both tie directions for Euler-walk child navigation.
class PlusMinusRangeMinimum {
 public:
  enum class Tie { First, Last };

  PlusMinusRangeMinimum() = default;
  explicit PlusMinusRangeMinimum(std::vector<int> depth) { build(std::move(depth)); }

  void build(std::vector<int> depth) {
    depth_ = std::move(depth);
    types_.clear();
    micro_.clear();
    macro_.clear();
    logs_.clear();
    if (depth_.empty()) return;
    std::size_t lg = 0;
    for (std::size_t m = depth_.size(); m > 1; m >>= 1) ++lg;
    width_ = std::max<std::size_t>(1, lg / 2);
    const std::size_t type_count = std::size_t{1} << (width_ - 1);
    micro_.resize(type_count * width_ * width_);
    std::vector<int> local(width_);
    for (std::size_t type = 0; type < type_count; ++type) {
      local[0] = 0;
      for (std::size_t x = 1; x < width_; ++x)
        local[x] = local[x - 1] + ((type >> (x - 1)) & 1 ? 1 : -1);
      for (std::size_t lo = 0; lo < width_; ++lo) {
        std::size_t first = lo, last = lo;
        for (std::size_t hi = lo; hi < width_; ++hi) {
          if (local[hi] < local[first]) first = last = hi;
          else if (local[hi] == local[first]) last = hi;
          micro_[micro_offset(type, lo, hi)] = {
              static_cast<std::uint8_t>(first), static_cast<std::uint8_t>(last)};
        }
      }
    }
    const std::size_t blocks = (depth_.size() + width_ - 1) / width_;
    types_.assign(blocks, 0);
    std::vector<Extrema> minima(blocks);
    for (std::size_t block = 0; block < blocks; ++block) {
      const std::size_t start = block * width_;
      const std::size_t end = std::min(depth_.size(), start + width_);
      for (std::size_t x = start + 1; x < end; ++x) {
        assert(depth_[x] == depth_[x - 1] + 1 || depth_[x] == depth_[x - 1] - 1);
        if (depth_[x] > depth_[x - 1]) types_[block] |= std::size_t{1} << (x - start - 1);
      }
      minima[block] = within_block(block, 0, end - start - 1);
    }
    logs_.assign(blocks + 1, 0);
    for (std::size_t x = 2; x <= blocks; ++x) logs_[x] = logs_[x / 2] + 1;
    macro_.push_back(std::move(minima));
    for (std::size_t level = 1; level <= logs_[blocks]; ++level) {
      const std::size_t half = std::size_t{1} << (level - 1);
      std::vector<Extrema> row(blocks - 2 * half + 1);
      for (std::size_t x = 0; x < row.size(); ++x)
        row[x] = combine(macro_[level - 1][x], macro_[level - 1][x + half]);
      macro_.push_back(std::move(row));
    }
  }

  std::size_t argmin(std::size_t lo, std::size_t hi, Tie tie = Tie::First,
                     std::size_t *table_reads = nullptr) const {
    assert(lo <= hi && hi < depth_.size());
    const std::size_t left = lo / width_, right = hi / width_;
    if (left == right) {
      if (table_reads) ++*table_reads;
      const Extrema result = within_block(left, lo % width_, hi % width_);
      return tie == Tie::First ? result.first : result.last;
    }
    if (table_reads) *table_reads += 2;
    Extrema result = combine(within_block(left, lo % width_, width_ - 1),
                             within_block(right, 0, hi % width_));
    if (right > left + 1) {
      const std::size_t count = right - left - 1, level = logs_[count];
      const std::size_t span = std::size_t{1} << level;
      if (table_reads) *table_reads += 2;
      result = combine(result, combine(macro_[level][left + 1],
                                       macro_[level][right - span]));
    }
    return tie == Tie::First ? result.first : result.last;
  }

  // Logical storage in bytes (excluding vector headers/capacity). Used by
  // scaling witnesses, independent of elapsed-time noise.
  std::size_t storage_bytes() const {
    std::size_t bytes = depth_.size() * sizeof(int) +
        (types_.size() + logs_.size()) * sizeof(std::size_t) +
        micro_.size() * sizeof(Offsets);
    for (const auto &row : macro_) bytes += row.size() * sizeof(Extrema);
    return bytes;
  }

 private:
  struct Offsets { std::uint8_t first = 0, last = 0; };
  struct Extrema { std::size_t first = 0, last = 0; };

  std::size_t micro_offset(std::size_t type, std::size_t lo, std::size_t hi) const {
    return (type * width_ + lo) * width_ + hi;
  }
  Extrema within_block(std::size_t block, std::size_t lo, std::size_t hi) const {
    const Offsets offset = micro_[micro_offset(types_[block], lo, hi)];
    return {block * width_ + offset.first, block * width_ + offset.last};
  }
  Extrema combine(Extrema a, Extrema b) const {
    if (depth_[a.first] < depth_[b.first]) return a;
    if (depth_[b.first] < depth_[a.first]) return b;
    return {std::min(a.first, b.first), std::max(a.last, b.last)};
  }

  std::size_t width_ = 1;
  std::vector<int> depth_;
  std::vector<std::size_t> types_, logs_;
  std::vector<Offsets> micro_;
  std::vector<std::vector<Extrema>> macro_;
};

// General static RMQ, reduced to +/-1 RMQ by a Cartesian tree. Strictly
// greater values are popped, so ties select the leftmost input minimum.
// Iterative construction/traversal also handles monotone arrays of arbitrary
// depth without recursion or a fixed number of sparse-table levels.
template <class Value>
class StaticRangeMinimum {
 public:
  StaticRangeMinimum() = default;
  explicit StaticRangeMinimum(const std::vector<Value> &values) {
    const std::size_t n = values.size();
    if (!n) return;
    std::vector<std::size_t> left(n, n), right(n, n), stack;
    stack.reserve(n);
    for (std::size_t x = 0; x < n; ++x) {
      std::size_t popped = n;
      while (!stack.empty() && values[stack.back()] > values[x]) {
        popped = stack.back();
        stack.pop_back();
      }
      if (!stack.empty()) right[stack.back()] = x;
      left[x] = popped;
      stack.push_back(x);
    }
    const std::size_t root = stack.front();
    first_.resize(n);
    euler_.reserve(2 * n - 1);
    std::vector<int> depths;
    depths.reserve(2 * n - 1);
    struct Frame { std::size_t node; unsigned next; int depth; };
    std::vector<Frame> walk{{root, 0, 0}};
    auto append = [&](std::size_t node, int depth) {
      euler_.push_back(node);
      depths.push_back(depth);
    };
    first_[root] = 0;
    append(root, 0);
    while (!walk.empty()) {
      Frame &frame = walk.back();
      if (frame.next < 2) {
        const std::size_t child = frame.next++ == 0 ? left[frame.node] : right[frame.node];
        if (child == n) continue;
        const int depth = frame.depth + 1;
        first_[child] = euler_.size();
        append(child, depth);
        walk.push_back({child, 0, depth});
      } else {
        walk.pop_back();
        if (!walk.empty()) append(walk.back().node, walk.back().depth);
      }
    }
    rmq_.build(std::move(depths));
  }

  std::size_t argmin(std::size_t lo, std::size_t hi,
                     std::size_t *table_reads = nullptr) const {
    assert(lo <= hi && hi < first_.size());
    std::size_t a = first_[lo], b = first_[hi];
    if (a > b) std::swap(a, b);
    return euler_[rmq_.argmin(a, b, PlusMinusRangeMinimum::Tie::First, table_reads)];
  }
  std::size_t storage_bytes() const {
    return (first_.size() + euler_.size()) * sizeof(std::size_t) + rmq_.storage_bytes();
  }

 private:
  std::vector<std::size_t> first_, euler_;
  PlusMinusRangeMinimum rmq_;
};

}  // namespace cparty::scfg2::geometry
#endif

#ifndef CPARTY_SCFG2_GEOMETRY_STATIC_SCAFFOLD_GEOMETRY_HH
#define CPARTY_SCFG2_GEOMETRY_STATIC_SCAFFOLD_GEOMETRY_HH

#include "scfg2/types.hh"
#include "scfg2/geometry/static_range_minimum.hh"
#include "scfg2/geometry/scaffold_input.hh"

#include <limits>
#include <stdexcept>
#include <string_view>

namespace cparty::scfg2::geometry {

// Observation-only index of a static laminar scaffold. Closing endpoints have
// the outer parent, exactly as in legacy CParty. O(n) setup/storage; all public
// queries take O(1). No grammar, weight, or ownership decisions live here.
class StaticScaffoldGeometry {
 public:
  StaticScaffoldGeometry() = default;
  explicit StaticScaffoldGeometry(std::vector<cand_pos_t> partners)
      : partner_(std::move(partners)) { build(); }

  StaticScaffoldGeometry(std::string_view structure, cand_pos_t n) {
    if (n < 0 || static_cast<std::size_t>(n) != structure.size())
      throw std::invalid_argument("scaffold length mismatch");
    std::string error;
    if (!validate_scaffold_input(structure, &error))
      throw std::invalid_argument(error);
    partner_.assign(static_cast<std::size_t>(n) + 1, -2);
    std::vector<cand_pos_t> stack;
    for (cand_pos_t pos = 1; pos <= n; ++pos) {
      const char ch = structure[static_cast<std::size_t>(pos - 1)];
      if (ch == 'x' || ch == 'X') partner_[pos] = -1;
      if (ch == '(') stack.push_back(pos);
      if (ch == ')') {
        if (stack.empty()) throw std::invalid_argument("unbalanced scaffold");
        const cand_pos_t opening = stack.back();
        stack.pop_back();
        partner_[pos] = opening;
        partner_[opening] = pos;
      }
    }
    if (!stack.empty()) throw std::invalid_argument("unbalanced scaffold");
    build();
  }

  bool usable() const { return usable_; }
  cand_pos_t size() const { return n_; }
  cand_pos_t partner(cand_pos_t pos) const {
    return pos >= 1 && pos <= n_ ? partner_[pos] : -1;
  }
  cand_pos_t parent(cand_pos_t pos) const {
    return valid_node(pos) ? parent_[pos] : -1;
  }
  cand_pos_t unpaired_prefix(cand_pos_t pos) const {
    return pos >= 1 && pos <= n_ ? up_[pos] : 0;
  }
  const std::vector<cand_pos_t> &unpaired_counts() const { return up_; }
  int depth(cand_pos_t node) const { return level_[node]; }

  cand_pos_t lca(cand_pos_t a, cand_pos_t b) const {
    std::size_t lo = first_[a], hi = first_[b];
    if (lo > hi) std::swap(lo, hi);
    return euler_[rmq_.argmin(lo, hi)];
  }

  // The last visit to ancestor before node's first visit is immediately
  // followed by the child on that path. This removes level-ancestor search.
  cand_pos_t child_below(cand_pos_t ancestor, cand_pos_t node) const {
    assert(ancestor != node && lca(ancestor, node) == ancestor);
    const std::size_t visit = rmq_.argmin(first_[ancestor], first_[node],
                                         PlusMinusRangeMinimum::Tie::Last);
    return euler_[visit + 1];
  }

  // Stack top just after/before the specified coordinate. Opening endpoints
  // belong to the following gap, closing endpoints to the preceding gap.
  cand_pos_t gap_after(cand_pos_t pos) const {
    if (!pos) return 0;
    return partner_[pos] > pos ? pos : parent_[pos];
  }
  cand_pos_t gap_before(cand_pos_t pos) const {
    if (pos > n_) return 0;
    return partner_[pos] > 0 && partner_[pos] < pos ? partner_[pos] : parent_[pos];
  }

  cand_pos_t bp(cand_pos_t i, cand_pos_t l) const {
    if (!band_base(l)) return -2;
    return parent_[l] < i ? -1 : parent_[l];
  }
  cand_pos_t Bp(cand_pos_t l, cand_pos_t j) const {
    if (!band_base(l)) return -2;
    return partner_[parent_[l]] > j ? -1 : partner_[parent_[l]];
  }
  cand_pos_t b(cand_pos_t i, cand_pos_t l) const {
    if (!band_base(l)) return -2;
    if (parent_[l] < i) return -1;
    const cand_pos_t ancestor = lca(i, l);
    if (i == ancestor) return i;
    const std::size_t visit = rmq_.argmin(first_[i], first_[l],
                                         PlusMinusRangeMinimum::Tie::Last);
    const cand_pos_t child = euler_[visit + 1];
    return child < l && partner_[child] > l ? child : -100;
  }
  cand_pos_t B(cand_pos_t l, cand_pos_t j) const {
    if (!band_base(l)) return -2;
    if (partner_[parent_[l]] > j) return -1;
    const std::size_t visit = rmq_.argmin(first_[l], first_[j]);
    const cand_pos_t ancestor = euler_[visit];
    if (j == ancestor) return j;
    const cand_pos_t child = euler_[visit - 1];
    return child < j && partner_[child] > l ? partner_[child] : -100;
  }
  bool weakly_closed(cand_pos_t i, cand_pos_t j) const {
    if (i == 1 && j == 0) return true;
    if (!valid_node(i) || !valid_node(j)) return false;
    if ((partner_[i] > 0 && partner_[i] < i) || partner_[j] > j) return false;
    return parent_[i] == parent_[j] && level_[i] == level_[j];
  }

  std::size_t storage_bytes() const {
    return (partner_.size() + parent_.size() + up_.size() + euler_.size()) *
        sizeof(cand_pos_t) + level_.size() * sizeof(int) +
        first_.size() * sizeof(std::size_t) + rmq_.storage_bytes();
  }

 private:
  bool valid_node(cand_pos_t pos) const { return usable_ && pos >= 0 && pos <= n_; }
  bool band_base(cand_pos_t pos) const {
    return pos >= 1 && pos <= n_ && partner_[pos] < 0 && parent_[pos] > 0;
  }

  void build() {
    if (partner_.empty() || partner_.size() - 1 >
        static_cast<std::size_t>(std::numeric_limits<cand_pos_t>::max())) return;
    n_ = static_cast<cand_pos_t>(partner_.size() - 1);
    const std::size_t size = partner_.size();
    parent_.assign(size, -1);
    up_.assign(size, 0);
    std::vector<cand_pos_t> stack{0}, first_child(size, -1), last_child(size, -1),
        next_sibling(size, -1);
    cand_pos_t count = 0;
    for (cand_pos_t pos = 1; pos <= n_; ++pos) {
      const cand_pos_t pair = partner_[pos];
      if (pair > 0 && (pair > n_ || pair == pos || partner_[pair] != pos)) return;
      if (pair > 0 && pair < pos) {
        if (stack.size() == 1 || stack.back() != pair) return;  // crossing
        stack.pop_back();
        count = 0;
      }
      const cand_pos_t outer = stack.back();
      parent_[pos] = outer;
      if (last_child[outer] < 0) first_child[outer] = pos;
      else next_sibling[last_child[outer]] = pos;
      last_child[outer] = pos;
      if (pair > pos) {
        stack.push_back(pos);
        count = 0;
      }
      up_[pos] = count++;
    }
    if (stack.size() != 1) return;
    first_.resize(size);
    level_.resize(size);
    euler_.reserve(2 * size - 1);
    std::vector<int> depths;
    depths.reserve(2 * size - 1);
    struct Frame { cand_pos_t node, next; };
    std::vector<Frame> walk{{0, first_child[0]}};
    first_[0] = 0;
    level_[0] = 0;
    euler_.push_back(0);
    depths.push_back(0);
    while (!walk.empty()) {
      Frame &frame = walk.back();
      if (frame.next >= 0) {
        const cand_pos_t child = frame.next;
        frame.next = next_sibling[child];
        first_[child] = euler_.size();
        level_[child] = level_[frame.node] + 1;
        euler_.push_back(child);
        depths.push_back(level_[child]);
        walk.push_back({child, first_child[child]});
      } else {
        walk.pop_back();
        if (!walk.empty()) {
          euler_.push_back(walk.back().node);
          depths.push_back(level_[walk.back().node]);
        }
      }
    }
    rmq_.build(std::move(depths));
    usable_ = true;
  }

  cand_pos_t n_ = 0;
  bool usable_ = false;
  std::vector<cand_pos_t> partner_, parent_, up_, euler_;
  std::vector<int> level_;
  std::vector<std::size_t> first_;
  PlusMinusRangeMinimum rmq_;
};

}  // namespace cparty::scfg2::geometry
#endif

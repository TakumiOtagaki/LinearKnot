#ifndef CPARTY_SCFG2_BEAM_V_RETAINED_AGENDA_HH
#define CPARTY_SCFG2_BEAM_V_RETAINED_AGENDA_HH

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/beam/memory_observation.hh"

#include <algorithm>
#include <bitset>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct VRetainedAgendaStats {
  std::size_t hairpin_parent_probes = 0;
  std::size_t hairpin_sources = 0;
  std::size_t internal_parent_probes = 0;
  std::size_t internal_sources = 0;
  std::size_t vm_parent_probes = 0;
  std::size_t vm_sources = 0;
};

// Source-owned V agenda. INTERNAL parents are emitted forward from retained V
// children, while V_VM parents are emitted from the retained same-column VM
// view. Hairpins have no semantic child and therefore remain an explicit
// all-start seed scan; isolating that scan is intentional because replacing it
// requires a separate left-to-right seed-state design.
template <class Context>
class VRetainedAgenda {
 public:
  void build_column(cand_pos_t j,
                    const std::vector<cand_pos_t> &retained_vm_starts,
                    std::size_t max_hairpin_unpaired,
                    Context &ctx) {
    column_ = j;
    parent_starts_.clear();

    cand_pos_t min_hairpin_i = 1;
    if (max_hairpin_unpaired != std::numeric_limits<std::size_t>::max()) {
      const auto capped = static_cast<cand_pos_t>(max_hairpin_unpaired);
      min_hairpin_i = std::max<cand_pos_t>(1, j - capped - 1);
    }
    for (cand_pos_t i = j; i >= min_hairpin_i; --i) {
      ++stats_.hairpin_parent_probes;
      const deductions::ItemKey parent{
          grammar::NonTerminal::V, i, j, -1, -1};
      append_generated(parent, grammar::RuleId::V_HAIRPIN, nullptr, ctx);
    }
    // The hairpin cap is a pruning heuristic for G' pairs.  A scaffold pair
    // closing at j must stay formable whatever its loop length.
    const cand_pos_t scaffold_i = ctx.pair_partner(j);
    if (scaffold_i > 0 && scaffold_i < min_hairpin_i) {
      ++stats_.hairpin_parent_probes;
      const deductions::ItemKey parent{
          grammar::NonTerminal::V, scaffold_i, j, -1, -1};
      append_generated(parent, grammar::RuleId::V_HAIRPIN, nullptr, ctx);
    }
    for (const cand_pos_t i : retained_vm_starts) {
      ++stats_.vm_parent_probes;
      const deductions::ItemKey parent{
          grammar::NonTerminal::V, i, j, -1, -1};
      append_generated(parent, grammar::RuleId::V_VM, nullptr, ctx);
    }

    const auto found = parent_starts_by_column_.find(j);
    if (found != parent_starts_by_column_.end()) {
      parent_starts_.insert(
          parent_starts_.end(), found->second.begin(), found->second.end());
      parent_starts_by_column_.erase(found);
    }
    std::sort(parent_starts_.begin(),
              parent_starts_.end(),
              std::greater<cand_pos_t>());
    parent_starts_.erase(
        std::unique(parent_starts_.begin(), parent_starts_.end()),
        parent_starts_.end());
  }

  void observe_committed_v_column(
      cand_pos_t child_j,
      const std::vector<cand_pos_t> &retained_child_starts,
      cand_pos_t n,
      Context &ctx) {
    for (const cand_pos_t child_i : retained_child_starts) {
      const deductions::ItemKey child{
          grammar::NonTerminal::V, child_i, child_j, -1, -1};
      if (child_i <= 1) continue;
      // The shared split guard requires left <= unpaired_prefix(child_i-1)
      // and left+right <= MAXLOOP. Exclude only the necessarily rejected part
      // of the rectangular window; the final shared predicate still runs.
      const cand_pos_t max_left = std::min<cand_pos_t>(
          MAXLOOP, ctx.unpaired_prefix(child_i - 1));
      const cand_pos_t min_parent_i =
          std::max<cand_pos_t>(1, child_i - max_left - 1);
      for (cand_pos_t parent_i = min_parent_i;
           parent_i < child_i;
           ++parent_i) {
        const cand_pos_t left = child_i - parent_i - 1;
        const cand_pos_t max_parent_j =
            std::min<cand_pos_t>(n, child_j + MAXLOOP + 1 - left);
        for (cand_pos_t parent_j = child_j + 1;
             parent_j <= max_parent_j;
             ++parent_j) {
          ++stats_.internal_parent_probes;
          const deductions::ItemKey parent{
              grammar::NonTerminal::V, parent_i, parent_j, -1, -1};
          append_generated(
              parent, grammar::RuleId::V_INTERNAL, &child, ctx);
        }
      }
    }
  }

  const std::vector<cand_pos_t> &parent_starts() const {
    return parent_starts_;
  }

  std::vector<deductions::Deduction> take_for(
      const deductions::ItemKey &parent) {
    const auto found = deductions_by_parent_.find(parent);
    if (found == deductions_by_parent_.end()) return {};
    auto out = std::move(found->second.sources);
    deductions_by_parent_.erase(found);
    return out;
  }

  cand_pos_t column() const { return column_; }
  const VRetainedAgendaStats &stats() const { return stats_; }

  // Parent-side admission may skip take_for(). Such candidates are still
  // expired after the complete V phase; no later column can consume them.
  void release_column() {
    for (const cand_pos_t i : parent_starts_)
      deductions_by_parent_.erase({grammar::NonTerminal::V, i, column_, -1, -1});
    parent_starts_.clear();
  }
  std::size_t pending_parent_count() const { return deductions_by_parent_.size(); }
  std::size_t pending_column_count() const { return parent_starts_by_column_.size(); }

  ContainerMemory memory_snapshot(const char *name) const {
    auto out = hash_memory(name, deductions_by_parent_);
    for (const auto &[parent, cell] : deductions_by_parent_) {
      (void)parent;
      out.vector_size_bytes += cell.sources.size() * sizeof(deductions::Deduction);
      out.vector_capacity_bytes += cell.sources.capacity() * sizeof(deductions::Deduction);
    }
    out.add(hash_vector_memory("", parent_starts_by_column_));
    out.add(vector_memory("", parent_starts_));
    return out;
  }

 private:
  void append_generated(const deductions::ItemKey &parent,
                        grammar::RuleId rule,
                        const deductions::ItemKey *required_child,
                        Context &ctx) {
    if (rule == grammar::RuleId::V_INTERNAL && required_child != nullptr) {
      // One (parent, child) probe: test the split directly instead of
      // regenerating every internal loop of the parent (O(L^2) per probe).
      if (deductions::exact_v_internal_split_allows(
              parent, required_child->i, required_child->j, ctx)) {
        grammar::RuleSplit split{};
        split.k = required_child->i;
        split.l = required_child->j;
        auto source = deductions::expand_basic_rule<grammar::RuleId::V_INTERNAL>(
            parent, split);
        if (deductions::is_valid(source))
          append_source(parent, std::move(source), required_child);
      }
      return;
    }
    for (auto source : deductions::generate_exact_basic_deductions(parent, rule, ctx))
      append_source(parent, std::move(source), required_child);
  }

  void append_source(const deductions::ItemKey &parent,
                     deductions::Deduction source,
                     const deductions::ItemKey *required_child) {
    if (required_child != nullptr &&
        (source.child_count != 1 ||
         !deductions::ItemKeyEq{}(source.children[0].item, *required_child))) return;
    const auto rule = source.rule;
    auto [found, inserted] = deductions_by_parent_.try_emplace(parent);
    auto &cell = found->second;
    if (rule == grammar::RuleId::V_INTERNAL) {
      // The shared split guard bounds each wing by MAXLOOP. Keep the vector's
      // first-arrival order, but test repeated (k,l) in constant space/time
      // per queued parent instead of rescanning its accumulated sources.
      const auto left = source.split.k - parent.i - 1;
      const auto right = parent.j - source.split.l - 1;
      if (left < 0 || left > MAXLOOP || right < 0 || right > MAXLOOP)
        throw std::logic_error("V internal source exceeds the shared loop window");
      const auto slot = static_cast<std::size_t>(left) * loop_width + right;
      if (cell.internal_splits.test(slot)) return;
      cell.internal_splits.set(slot);
    } else {
      const bool duplicate = std::any_of(
          cell.sources.begin(), cell.sources.end(), [&](const auto &old) {
            return old.rule == source.rule && old.split.k == source.split.k &&
                   old.split.l == source.split.l;
          });
      if (duplicate) return;
    }
    cell.sources.push_back(std::move(source));
    if (inserted) {
      if (parent.j == column_) {
        parent_starts_.push_back(parent.i);
      } else {
        parent_starts_by_column_[parent.j].push_back(parent.i);
      }
    }
    if (rule == grammar::RuleId::V_HAIRPIN) {
      ++stats_.hairpin_sources;
    } else if (rule == grammar::RuleId::V_INTERNAL) {
      ++stats_.internal_sources;
    } else {
      ++stats_.vm_sources;
    }
  }

  static constexpr std::size_t loop_width = MAXLOOP + 1;
  struct ParentSources {
    std::vector<deductions::Deduction> sources;
    // Released together with the existing parent source vector in take_for.
    std::bitset<loop_width * loop_width> internal_splits;
  };
  cand_pos_t column_ = -1;
  VRetainedAgendaStats stats_{};
  std::unordered_map<deductions::ItemKey,
                     ParentSources,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      deductions_by_parent_;
  std::unordered_map<cand_pos_t, std::vector<cand_pos_t>>
      parent_starts_by_column_;
  std::vector<cand_pos_t> parent_starts_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_V_RETAINED_AGENDA_HH

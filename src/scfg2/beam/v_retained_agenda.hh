#ifndef CPARTY_SCFG2_BEAM_V_RETAINED_AGENDA_HH
#define CPARTY_SCFG2_BEAM_V_RETAINED_AGENDA_HH

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/beam/memory_observation.hh"

#include <algorithm>
#include <array>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::replay {
template <bool LogDomain> class WFinalInsideContext;
}

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
  VRetainedAgenda() = default;
  VRetainedAgenda(const VRetainedAgenda &) = default;
  VRetainedAgenda(VRetainedAgenda &&) noexcept = default;
  VRetainedAgenda &operator=(VRetainedAgenda &&) noexcept = default;
  VRetainedAgenda &operator=(const VRetainedAgenda &other) {
    if (this != &other) {
      // Keep parent lists, map owners and numeric pending cells together if
      // any allocation fails while copying either representation.
      VRetainedAgenda copy(other);
      *this = std::move(copy);
    }
    return *this;
  }
  // Explicit internal opt-in before the first column/event. A failed optional
  // allocation leaves the old packed map untouched. No empty ParentSources
  // (or their vectors) are constructed in the preallocated cell metadata.
  bool enable_numeric_pending(cand_pos_t n) {
    if (n <= 0 || column_ != -1 || pending_parent_count() ||
        !parent_starts_by_column_.empty() || !parent_starts_.empty()) return false;
    const auto stride = static_cast<std::size_t>(n) + 1;
    try {
      if (stride > std::vector<std::optional<ParentSources>>{}.max_size() / ring_columns ||
          stride > std::vector<cand_pos_t>{}.max_size() / ring_columns ||
          stride > std::vector<unsigned char>{}.max_size() / ring_columns) return false;
      NumericPending replacement;
      replacement.n = n;
      replacement.stride = stride;
      replacement.cells.resize(stride * ring_columns);
      replacement.touched.resize(stride * ring_columns);
      replacement.marked.resize(stride * ring_columns, false);
      numeric_pending_.emplace(std::move(replacement));
      return true;
    } catch (const std::bad_alloc &) { return false; }
      catch (const std::length_error &) { return false; }
  }

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
      append_generated<grammar::RuleId::V_HAIRPIN>(parent, nullptr, ctx);
    }
    // The hairpin cap is a pruning heuristic for G' pairs.  A scaffold pair
    // closing at j must stay formable whatever its loop length.
    const cand_pos_t scaffold_i = ctx.pair_partner(j);
    if (scaffold_i > 0 && scaffold_i < min_hairpin_i) {
      ++stats_.hairpin_parent_probes;
      const deductions::ItemKey parent{
          grammar::NonTerminal::V, scaffold_i, j, -1, -1};
      append_generated<grammar::RuleId::V_HAIRPIN>(parent, nullptr, ctx);
    }
    for (const cand_pos_t i : retained_vm_starts) {
      ++stats_.vm_parent_probes;
      const deductions::ItemKey parent{
          grammar::NonTerminal::V, i, j, -1, -1};
      append_generated<grammar::RuleId::V_VM>(parent, nullptr, ctx);
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
      const cand_pos_t left_up = ctx.unpaired_prefix(child_i - 1);
      const cand_pos_t max_left = std::min<cand_pos_t>(MAXLOOP, left_up);
      std::optional<deductions::ExactVInternalChildQueryCache<Context>> child_queries;
      if constexpr (std::is_same_v<std::decay_t<Context>,
                                  replay::WFinalInsideContext<true>>) {
        // Ring enablement is itself gated by numeric lookup + no memory
        // observation at the standard runner boundary. Recheck this event's
        // complete domain before the unchecked base-array pairing query.
        if (numeric_pending_ && numeric_pending_->n > 0 &&
            n == numeric_pending_->n && child_i >= 1 && child_i < child_j &&
            child_j <= n && n == ctx.sequence_length())
          child_queries.emplace(ctx, child_i, child_j, left_up);
      }
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
          append_generated<grammar::RuleId::V_INTERNAL>(
              parent, &child, ctx, child_queries ? &*child_queries : nullptr);
        }
      }
    }
  }

  const std::vector<cand_pos_t> &parent_starts() const {
    return parent_starts_;
  }

  std::vector<deductions::Deduction> take_for(
      const deductions::ItemKey &parent) {
    // The packed index represents only canonical V parents. Other complete
    // ItemKeys must keep the old full-key map's absent-key behavior.
    if (parent.nonterminal != grammar::NonTerminal::V ||
        parent.ip != -1 || parent.jp != -1 ||
        !deductions::is_valid_primary_span(parent)) return {};
    if (const auto index = ring_cell(parent)) {
      std::vector<deductions::Deduction> out;
      expand_sources(parent, *numeric_pending_->cells[*index], out);
      erase_ring_cell(*index);
      return out;
    }
    if (numeric_pending_ && numeric_pending_->n && !numeric_pending_->mixed) return {};
    const auto found = deductions_by_parent_.find(parent_key(parent.i, parent.j));
    if (found == deductions_by_parent_.end()) return {};
    std::vector<deductions::Deduction> out;
    expand_sources(parent, found->second, out);
    deductions_by_parent_.erase(found);
    return out;
  }

  cand_pos_t column() const { return column_; }
  const VRetainedAgendaStats &stats() const { return stats_; }

  // Parent-side admission may skip take_for(). Such candidates are still
  // expired after the complete V phase; no later column can consume them.
  void release_column() {
    for (const cand_pos_t i : parent_starts_) {
      const deductions::ItemKey parent{grammar::NonTerminal::V, i, column_, -1, -1};
      if (const auto index = ring_cell(parent)) erase_ring_cell(*index);
      else if (!numeric_pending_ || !numeric_pending_->n || numeric_pending_->mixed)
        deductions_by_parent_.erase(parent_key(i, column_));
    }
    if (numeric_pending_ && numeric_pending_->n && column_ >= 1 &&
        column_ <= numeric_pending_->n) {
      const auto row = static_cast<std::size_t>(column_) % ring_columns;
      if (numeric_pending_->columns[row].stamp == column_ &&
          !numeric_pending_->columns[row].live) clear_ring_row(row);
    }
    parent_starts_.clear();
  }
  std::size_t pending_parent_count() const {
    return deductions_by_parent_.size() + (numeric_pending_ ? numeric_pending_->live : 0);
  }
  std::size_t pending_column_count() const { return parent_starts_by_column_.size(); }

  ContainerMemory memory_snapshot(const char *name) const {
    auto out = hash_memory(name, deductions_by_parent_);
    for (const auto &[parent, cell] : deductions_by_parent_) {
      (void)parent;
      out.vector_size_bytes += cell.sources.size() * sizeof(CompactSource);
      out.vector_capacity_bytes += cell.sources.capacity() * sizeof(CompactSource);
    }
    if (numeric_pending_) {
      out.add(vector_memory("", numeric_pending_->cells));
      out.add(vector_memory("", numeric_pending_->touched));
      out.add(vector_memory("", numeric_pending_->marked));
      for (const auto &cell : numeric_pending_->cells) {
        if (!cell) continue;
        out.vector_size_bytes += cell->sources.size() * sizeof(CompactSource);
        out.vector_capacity_bytes += cell->sources.capacity() * sizeof(CompactSource);
      }
    }
    out.add(hash_vector_memory("", parent_starts_by_column_));
    out.add(vector_memory("", parent_starts_));
    return out;
  }

 private:
  template <grammar::RuleId Rule>
  void append_generated(const deductions::ItemKey &parent,
                        const deductions::ItemKey *required_child,
                        Context &ctx,
                        const deductions::ExactVInternalChildQueryCache<Context> *child_queries = nullptr) {
    static_assert(Rule == grammar::RuleId::V_HAIRPIN ||
                  Rule == grammar::RuleId::V_INTERNAL || Rule == grammar::RuleId::V_VM);
    grammar::RuleSplit split{};
    if constexpr (Rule == grammar::RuleId::V_INTERNAL) {
      // One (parent, child) probe: test the split directly instead of
      // regenerating every internal loop of the parent (O(L^2) per probe).
      if (required_child == nullptr) return;
      const bool allowed = child_queries
          ? deductions::exact_v_internal_split_allows_impl(
                parent, required_child->i, required_child->j, *child_queries)
          : deductions::exact_v_internal_split_allows(
                parent, required_child->i, required_child->j, ctx);
      if (!allowed) return;
      split.k = required_child->i;
      split.l = required_child->j;
    } else if constexpr (Rule == grammar::RuleId::V_HAIRPIN) {
      if (!deductions::exact_v_hairpin_allows(ctx, parent.i, parent.j)) return;
    } else {
      if (!deductions::exact_v_pairing_state_allows(ctx, parent.i, parent.j)) return;
    }
    // Validate the same shared expansion before queuing it, including the
    // required child's complete identity. No full Deduction crosses the
    // append boundary or survives in the future-column queue.
    const auto source = deductions::expand_basic_rule<Rule>(parent, split);
    if (!deductions::is_valid(source)) return;
    if (required_child != nullptr &&
        (source.child_count != 1 ||
         !deductions::ItemKeyEq{}(source.children[0].item, *required_child))) return;
    append_source(parent, Rule, split.k, split.l);
  }

  // Only the three fixed V schemas enter this private queue. INTERNAL wing
  // lengths are bounded by the unchanged shared guard; absolute coordinates
  // are restored only by take_for(), through the unchanged shared expansion.
  struct CompactSource {
    std::uint8_t rule;
    std::uint8_t left = 0;
    std::uint8_t right = 0;
  };
  static_assert(MAXLOOP >= 0 &&
                    MAXLOOP <= std::numeric_limits<std::uint8_t>::max(),
                "compact V wings require an unsigned byte");
  static_assert(static_cast<unsigned>(grammar::RuleId::V_HAIRPIN) <= 255 &&
                    static_cast<unsigned>(grammar::RuleId::V_INTERNAL) <= 255 &&
                    static_cast<unsigned>(grammar::RuleId::V_VM) <= 255,
                "compact V rule tags require an unsigned byte");
  static_assert(sizeof(CompactSource) == 3,
                "compact V source must use three bytes");

  static std::uint64_t parent_key(cand_pos_t i, cand_pos_t j) {
    static_assert(std::numeric_limits<cand_pos_t>::digits <= 32,
                  "packed V parent coordinates require at most 32 bits");
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) << 32) |
        static_cast<std::uint32_t>(j);
  }

  void append_source(const deductions::ItemKey &parent,
                     grammar::RuleId rule, int k, int l) {
    auto [storage, inserted] = source_storage(parent);
    auto &cell = *storage;
    CompactSource source{static_cast<std::uint8_t>(rule)};
    if (rule == grammar::RuleId::V_INTERNAL) {
      // The shared split guard bounds each wing by MAXLOOP. Keep the vector's
      // first-arrival order, but test repeated (k,l) in constant space/time
      // per queued parent instead of rescanning its accumulated sources.
      const auto left = k - parent.i - 1;
      const auto right = parent.j - l - 1;
      if (left < 0 || left > MAXLOOP || right < 0 || right > MAXLOOP || left + right > MAXLOOP)
        throw std::logic_error("V internal source exceeds the shared loop window");
      // Row left contains MAXLOOP + 1 - left legal right wings.
      // Its prefix starts after left * (2 * loop_width + 1 - left) / 2 bits.
      const auto row = static_cast<std::size_t>(left);
      const auto slot = row * (2 * loop_width + 1 - row) / 2 + right;
      if (cell.internal_splits.test(slot)) return;
      cell.internal_splits.set(slot);
      source.left = static_cast<std::uint8_t>(left);
      source.right = static_cast<std::uint8_t>(right);
    } else {
      const bool duplicate = std::any_of(
          cell.sources.begin(), cell.sources.end(), [&](const auto &old) {
            // Both non-INTERNAL schemas have k/l == -1; their zero wings
            // preserve the old (rule,k,l) identity without storing sentinels.
            return old.rule == source.rule && old.left == source.left &&
                old.right == source.right;
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
    std::vector<CompactSource> sources;
    // Released together with the existing parent source vector in take_for.
    std::bitset<loop_width * (loop_width + 1) / 2> internal_splits;
  };
  static constexpr std::size_t ring_columns = MAXLOOP + 2;
  struct RingColumn {
    cand_pos_t stamp = -1;
    std::size_t live = 0;
    std::size_t touched = 0;
  };
  struct NumericPending {
    cand_pos_t n = 0;
    std::size_t stride = 0;
    std::array<RingColumn, ring_columns> columns{};
    std::vector<std::optional<ParentSources>> cells;
    std::vector<cand_pos_t> touched;
    std::vector<unsigned char> marked;
    std::size_t live = 0;
    bool mixed = false;

    NumericPending() = default;
    NumericPending(const NumericPending &) = default;
    NumericPending &operator=(const NumericPending &other) {
      if (this != &other) {
        NumericPending copy(other);
        *this = std::move(copy);
      }
      return *this;
    }
    NumericPending(NumericPending &&other) noexcept
        : n(std::exchange(other.n, 0)), stride(std::exchange(other.stride, 0)),
          columns(other.columns), cells(std::move(other.cells)),
          touched(std::move(other.touched)), marked(std::move(other.marked)),
          live(std::exchange(other.live, 0)), mixed(std::exchange(other.mixed, false)) {
      other.columns = {};
    }
    NumericPending &operator=(NumericPending &&other) noexcept {
      if (this == &other) {
        // The enclosing agenda's ordinary maps/vectors may empty themselves
        // on self-move. Keep no pending owner behind that moved outer state.
        NumericPending empty;
        return *this = std::move(empty);
      }
      n = std::exchange(other.n, 0);
      stride = std::exchange(other.stride, 0);
      columns = other.columns;
      cells = std::move(other.cells);
      touched = std::move(other.touched);
      marked = std::move(other.marked);
      live = std::exchange(other.live, 0);
      mixed = std::exchange(other.mixed, false);
      other.columns = {};
      return *this;
    }
  };

  static void expand_sources(const deductions::ItemKey &parent,
                            const ParentSources &cell,
                            std::vector<deductions::Deduction> &out) {
    out.reserve(cell.sources.size());
    // Expansion is pure; preserve every shared-schema field/arrival position.
    for (const auto &source : cell.sources) {
      grammar::RuleSplit split{};
      const auto rule = static_cast<grammar::RuleId>(source.rule);
      if (rule == grammar::RuleId::V_INTERNAL) {
        split.k = parent.i + 1 + source.left;
        split.l = parent.j - 1 - source.right;
      } // Hairpin/VM retain RuleSplit's exact k/l/p/q == -1 defaults.
      switch (rule) {
        case grammar::RuleId::V_HAIRPIN:
          out.push_back(deductions::expand_basic_rule<grammar::RuleId::V_HAIRPIN>(parent, split));
          break;
        case grammar::RuleId::V_INTERNAL:
          out.push_back(deductions::expand_basic_rule<grammar::RuleId::V_INTERNAL>(parent, split));
          break;
        case grammar::RuleId::V_VM:
          out.push_back(deductions::expand_basic_rule<grammar::RuleId::V_VM>(parent, split));
          break;
        default:
          throw std::logic_error("unexpected compact V source rule");
      }
    }
  }

  bool ring_key_supported(const deductions::ItemKey &parent) const {
    return numeric_pending_ && numeric_pending_->n > 0 &&
        parent.nonterminal == grammar::NonTerminal::V && parent.ip == -1 && parent.jp == -1 &&
        parent.i >= 1 && parent.j >= parent.i && parent.j <= numeric_pending_->n;
  }
  std::optional<std::size_t> ring_cell(const deductions::ItemKey &parent) const {
    if (!ring_key_supported(parent)) return std::nullopt;
    const auto &pending = *numeric_pending_;
    const auto row = static_cast<std::size_t>(parent.j) % ring_columns;
    if (pending.columns[row].stamp != parent.j) return std::nullopt;
    const auto index = row * pending.stride + static_cast<std::size_t>(parent.i);
    return pending.cells[index] ? std::optional<std::size_t>{index} : std::nullopt;
  }
  void erase_ring_cell(std::size_t index) {
    auto &pending = *numeric_pending_;
    pending.cells[index].reset(); // Release this source vector's full capacity.
    --pending.columns[index / pending.stride].live;
    --pending.live;
  }
  void clear_ring_row(std::size_t row) {
    auto &pending = *numeric_pending_;
    const auto offset = row * pending.stride;
    auto &column = pending.columns[row];
    // Call only once the row has no live owners. Marks deduplicate touched
    // starts even if a consumed parent is subsequently reintroduced.
    for (std::size_t t = 0; t < column.touched; ++t) {
      const auto index = offset + static_cast<std::size_t>(pending.touched[offset + t]);
      pending.cells[index].reset();
      pending.marked[index] = false;
    }
    column = {};
  }
  std::pair<ParentSources *, bool> source_storage(const deductions::ItemKey &parent) {
    if (numeric_pending_ && numeric_pending_->n > 0) {
      auto &pending = *numeric_pending_;
      if (const auto index = ring_cell(parent)) return {&*pending.cells[*index], false};
      // Out-of-domain/colliding events can coexist safely in the old map.
      // Once such an event occurred, keep an existing complete parent in
      // that owner even if its geometry later fits the current ring window.
      if (pending.mixed) {
        const auto found = deductions_by_parent_.find(parent_key(parent.i, parent.j));
        if (found != deductions_by_parent_.end()) return {&found->second, false};
      }
      if (ring_key_supported(parent) && column_ >= 1 && parent.j >= column_ &&
          parent.j - column_ <= MAXLOOP + 1) {
        const auto row = static_cast<std::size_t>(parent.j) % ring_columns;
        auto &column = pending.columns[row];
        if (column.stamp == parent.j || !column.live) {
          if (column.stamp != parent.j) { clear_ring_row(row); column.stamp = parent.j; }
          const auto offset = row * pending.stride;
          const auto index = offset + static_cast<std::size_t>(parent.i);
          pending.cells[index].emplace();
          ++column.live;
          ++pending.live;
          if (!pending.marked[index]) {
            pending.marked[index] = true;
            pending.touched[offset + column.touched++] = parent.i;
          }
          return {&*pending.cells[index], true};
        }
      }
      pending.mixed = true;
    }
    auto [found, inserted] = deductions_by_parent_.try_emplace(parent_key(parent.i, parent.j));
    return {&found->second, inserted};
  }
  cand_pos_t column_ = -1;
  VRetainedAgendaStats stats_{};
  std::unordered_map<std::uint64_t, ParentSources> deductions_by_parent_;
  std::unordered_map<cand_pos_t, std::vector<cand_pos_t>>
      parent_starts_by_column_;
  std::vector<cand_pos_t> parent_starts_;
  std::optional<NumericPending> numeric_pending_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_V_RETAINED_AGENDA_HH

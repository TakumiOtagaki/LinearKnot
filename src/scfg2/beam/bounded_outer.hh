// Bounded max completion contexts from pilot inside values.
#ifndef CPARTY_SCFG2_BEAM_BOUNDED_OUTER_HH
#define CPARTY_SCFG2_BEAM_BOUNDED_OUTER_HH

#include "scfg2/beam/outer_frontier.hh"
#include "scfg2/geometry/static_scaffold_geometry.hh"
#include <cstdint>
#include <memory>
#include <new>
#include <stdexcept>

namespace cparty::scfg2::beam {
namespace outer_detail {
struct Work {
  std::size_t input_records = 0, guide_records = 0, index_records = 0;
  std::size_t comparisons = 0, build_moves = 0, undo_moves = 0, cursor_steps = 0;
  std::size_t query_outputs = 0, proposals = 0, max_proposals = 0, checks = 0;
  double validation_seconds = 0;
  std::size_t weight_evaluations = 0, child_alpha_lookups = 0;
  std::size_t leaf_deductions = 0, zero_sibling_deductions = 0;
  double preparation_seconds = 0, proposal_seconds = 0, generation_seconds = 0;
  double propagation_seconds = 0, finalization_seconds = 0;
  std::array<double, outer_probe_families> family_generation{}, family_propagation{};
};
struct Record { Item key; double alpha; };

// Reusable position deduplication. Only words marked by this parent are sorted,
// enumerated and reset; clearing a parent never scans the O(n)-sized metadata.
class ProposalPositionBits {
 public:
  explicit ProposalPositionBits(int n)
      : upper_(std::int64_t(n) + 1), word_extent_(word_count_for(n)) {}

  bool enabled() const { return enabled_; }
  std::size_t word_count() const { return words_.size(); }
  std::size_t active_word_count() const { return touched_.size(); }
  std::size_t touched_capacity() const { return touched_.capacity(); }

  void sort_unique(std::vector<int> &positions) {
    if (positions.size() < 2) return;
    if (!enabled_) { reference_sort_unique(positions); return; }
    // Validate before changing either positions or the bitset. Custom contexts
    // may produce other sentinels; retain their original sort/unique behavior.
    for (const int position : positions) {
      if (position < -2 || std::int64_t(position) > upper_) {
        reference_sort_unique(positions);
        return;
      }
    }
    // Allocate once, and reserve before marking, so allocation failure cannot
    // leave dirty words. Distinct touched words <= min(word_extent, raw_count).
    // This index is optional: allocation failure retains the original sorter.
    try {
      if (words_.empty()) words_.resize(word_extent_, 0);
      touched_.reserve(std::min(words_.size(), positions.size()));
    } catch (const std::bad_alloc &) {
      disable(); reference_sort_unique(positions); return;
    } catch (const std::length_error &) {
      disable(); reference_sort_unique(positions); return;
    }
    for (const int position : positions) {
      const auto shifted = static_cast<std::size_t>(std::int64_t(position) + 2);
      const auto id = shifted / 64;
      auto &word = words_[id];
      if (!word) touched_.push_back(id);
      word |= std::uint64_t{1} << (shifted % 64);
    }
    std::sort(touched_.begin(), touched_.end());
    positions.clear();
    // The original vector's capacity is sufficient: unique_count <= raw_count.
    for (const auto id : touched_) {
      auto bits = words_[id];
      words_[id] = 0;
      while (bits) {
        const auto offset = static_cast<unsigned>(__builtin_ctzll(bits));
        positions.push_back(static_cast<int>(std::int64_t(id * 64 + offset) - 2));
        bits &= bits - 1;
      }
    }
    touched_.clear();
  }

 private:
  static void reference_sort_unique(std::vector<int> &positions) {
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
  }
  void disable() {
    enabled_ = false;
    std::vector<std::uint64_t>{}.swap(words_);
    std::vector<std::size_t>{}.swap(touched_);
  }
  static std::size_t word_count_for(int n) {
    if (n < 0) throw std::runtime_error("negative proposal position extent");
    return (static_cast<std::size_t>(n) + 4 + 63) / 64;
  }
  std::int64_t upper_;
  std::size_t word_extent_;
  bool enabled_ = true;
  std::vector<std::uint64_t> words_;
  std::vector<std::size_t> touched_;
};

// The same constant-size undo design as probe_scfg2_vm_outer_budget.py.
class PrefixTop {
 public:
  const std::vector<Record> &records;
  struct Undo { int position; std::size_t evicted; };
  std::vector<std::size_t> top;
  std::vector<Undo> undo;
  std::size_t cursor, width;
  int last_upper = std::numeric_limits<int>::max();
  Work &work;
  bool better(std::size_t a, std::size_t b) const {
    if (records[a].alpha != records[b].alpha) return records[a].alpha > records[b].alpha;
    return less_key(records[a].key, records[b].key);
  }
  PrefixTop(const std::vector<Record> &records, std::size_t width, Work &work)
      : records(records), cursor(records.size()), width(width), work(work) {
    work.index_records += records.size();
    for (std::size_t id = 0; id < records.size(); ++id) {
      if (id && records[id - 1].key.j >= records[id].key.j)
        throw std::runtime_error("prefix records not in increasing end order");
      std::size_t p = 0;
      while (p < top.size()) {
        ++work.comparisons;
        if (better(id, top[p])) break;
        ++p;
      }
      if (p >= width) { undo.push_back({-1, records.size()}); continue; }
      work.build_moves += top.size() - p + 1;
      top.insert(top.begin() + p, id);
      std::size_t evicted = records.size();
      if (top.size() > width) { evicted = top.back(); top.pop_back(); }
      undo.push_back({static_cast<int>(p), evicted});
    }
  }
  const std::vector<std::size_t> &query(int upper) {
    if (upper > last_upper) throw std::runtime_error("nonmonotone prefix query");
    last_upper = upper;
    while (cursor && records[cursor - 1].key.j > upper) {
      const auto &u = undo[--cursor];
      ++work.cursor_steps;
      if (u.position < 0) continue;
      const auto p = static_cast<std::size_t>(u.position);
      if (top.at(p) != cursor) throw std::runtime_error("prefix undo mismatch");
      work.undo_moves += top.size() - p;
      top.erase(top.begin() + p);
      if (u.evicted != records.size()) { top.push_back(u.evicted); ++work.undo_moves; }
    }
    work.query_outputs += top.size();
    return top;
  }
};

class Proposals {
 public:
  static constexpr std::size_t families = outer_probe_families;
  int n;
  std::size_t b;
  Work &work;
  Chart alpha{zero};
  std::vector<std::array<std::vector<Record>, families>> right, left;
  std::vector<std::array<int, families>> best_left, best_right;
  std::unordered_map<std::size_t, std::unique_ptr<PrefixTop>> indices;
  sc::geometry::StaticScaffoldGeometry scaffold;
  ProposalPositionBits position_bits;

  template <class Context>
  Proposals(const Chart &pilot, int n, std::size_t b, Context &ctx, Work &work)
      : n(n), b(b), work(work), right(n + 1), left(n + 1),
        best_left(n + 1), best_right(n + 1), position_bits(n) {
    if (!b) throw std::runtime_error("zero proposal width");
    std::vector<int> partners(n + 1, -2);
    for (int p = 1; p <= n; ++p) {
      if (ctx.pair_square(p) > 0) throw std::runtime_error("bounded guide requires laminar G");
      partners[p] = ctx.pair_partner(p);
    }
    scaffold = sc::geometry::StaticScaffoldGeometry(std::move(partners));
    if (!scaffold.usable()) throw std::runtime_error("invalid guide scaffold");
    // Bound effective guide width, including forced exceptions. Keep root W
    // separately; all other families retain at most b pilot values per column.
    std::vector<std::array<MaxPool, families>> pools(n + 1);
    for (const auto &[key, value] : pilot.values()) {
      ++work.input_records;
      if (value == zero) continue;
      pools.at(key.j).at(static_cast<std::size_t>(key.nonterminal)).offer(
          key, value, {key, 0}, key.nonterminal == NT::W && key.i == 1, b);
    }
    for (int j = 1; j <= n; ++j) {
      best_left[j].fill(-1); best_right[j].fill(-1);
      for (std::size_t f = 0; f < families; ++f) {
        auto &column = right[j][f];
        for (const auto &[key, entry] : pools[j][f].entries) {
          alpha.set(key, entry.value);
          column.push_back({key, entry.value});
        }
        std::sort(column.begin(), column.end(), [](const auto &a, const auto &b) {
          return less_key(a.key, b.key);
        });
        for (std::size_t p = 0; p < column.size(); ++p) {
          const auto &r = column[p];
          auto &best = best_right[j][f];
          if (best < 0 || r.alpha > column[best].alpha) best = p;
          if (f == static_cast<std::size_t>(NT::BE)) continue;
          auto &row = left[r.key.i][f];
          auto &row_best = best_left[r.key.i][f];
          if (row_best < 0 || r.alpha > row[row_best].alpha) row_best = row.size();
          row.push_back(r); // right columns are read in increasing end order
        }
      }
    }
    work.guide_records = alpha.values().size();
    if (work.guide_records > families * std::size_t(n) * (b + 1))
      throw std::runtime_error("guide alpha width bound violated");
  }
  const std::vector<Record> &column(NT nt, int end) const {
    static const std::vector<Record> empty;
    return end >= 1 && end <= n ? right[end][static_cast<std::size_t>(nt)] : empty;
  }
  // The guide alpha records are immutable after construction. Reuse the
  // sorted proposal columns for point reads instead of probing the duplicate
  // global chart. Full keys distinguish BE and noncanonical primary records.
  double get_alpha(const Item &key) const {
    const auto family = static_cast<std::size_t>(key.nonterminal);
    if (family >= families || key.j < 1 || key.j > n) return zero;
    const auto &records = right[key.j][family];
    const auto found = std::lower_bound(records.begin(), records.end(), key,
        [](const Record &record, const Item &target) {
          return less_key(record.key, target);
        });
    return found != records.end() && sc::deductions::ItemKeyEq{}(found->key, key)
        ? found->alpha : zero;
  }
  void from_right(std::vector<int> &out, NT nt, int end, int offset = 0) {
    for (const auto &r : column(nt, end)) out.push_back(r.key.i + offset);
  }
  void from_left(std::vector<int> &out, const Item &parent, NT nt,
                 int start, int upper, int offset) {
    if (start < 1 || start > n) return;
    const auto family = static_cast<std::size_t>(nt);
    const auto channel = static_cast<std::size_t>(parent.nonterminal) * families + family;
    const auto id = channel * (n + 1) + start;
    auto &index = indices[id];
    if (!index) index = std::make_unique<PrefixTop>(left[start][family], 2 * b + 1, work);
    for (auto p : index->query(upper)) out.push_back(index->records[p].key.j + offset);
  }
  // No range search or hidden scan: one global raw-alpha hint, then the shared
  // geometry/weight checks decide whether it is usable for this parent.
  void best_left_hint(std::vector<int> &out, NT nt, int start, int offset) {
    if (start < 1 || start > n) return;
    const auto f = static_cast<std::size_t>(nt);
    const int best = best_left[start][f];
    if (best >= 0) out.push_back(left[start][f][best].key.j + offset);
  }
  void best_right_hint(std::vector<int> &out, NT nt, int end, int offset) {
    if (end < 1 || end > n) return;
    const auto f = static_cast<std::size_t>(nt);
    const int best = best_right[end][f];
    if (best >= 0) out.push_back(right[end][f][best].key.i + offset);
  }
  void coverage(std::vector<int> &out, int lo, int hi) const {
    if (lo > hi) return;
    const auto length = static_cast<std::size_t>(hi - lo + 1);
    const auto count = std::min(b, length);
    for (std::size_t p = 0; p < count; ++p) {
      out.push_back(lo + p); // nearby completions
      out.push_back(lo + (count == 1 ? 0 : p * (length - 1) / (count - 1)));
    }
  }
  template <class Context>
  void fill(const Item &a, Context &ctx, std::vector<int> &out) {
    out.clear();
    const int i = a.i, j = a.j;
    const auto prefix = [&](NT nt, int start, int offset, int upper) {
      from_left(out, a, nt, start, upper, offset);
    };
    switch (a.nonterminal) {
      case NT::W: case NT::WI: case NT::WM: case NT::WIP:
        from_right(out, NT::V, j); from_right(out, NT::WMB, j);
        prefix(a.nonterminal, i, 1, j - 1);
        coverage(out, i, j - 1);
        break;
      case NT::VM:
        from_right(out, NT::WMv, j - 1); from_right(out, NT::WMp, j - 1);
        prefix(NT::WM, i + 1, 1, j - ctx.turn() - 2);
        coverage(out, i + 1, j - ctx.turn() - 1);
        break;
      case NT::VP: from_right(out, NT::VP, j - 1); from_right(out, NT::VPR, j - 1);
        from_right(out, NT::WIP, j - 1, -1);
        prefix(NT::WIP, i + 1, 1, j - 2);
        prefix(NT::VP, i + 1, 0, j - 2);
        prefix(NT::VPL, i + 1, 0, j - 2);
        coverage(out, i + 1, j - 1);
        break;
      case NT::VPL: {
        from_right(out, NT::VP, j);
        const auto boundary = std::min(static_cast<cand_pos_tu>(ctx.border_b(i, j)),
                                       static_cast<cand_pos_tu>(ctx.border_Bp(i, j)));
        coverage(out, i + 1, boundary <= static_cast<cand_pos_tu>(j)
            ? static_cast<int>(boundary) - 1 : j - 1);
        break;
      }
      case NT::VPR:
        from_right(out, NT::WIP, j, -1);
        prefix(NT::VP, i, 0, j - 1);
        coverage(out, std::max(i + 1, std::max(ctx.border_B(i, j), ctx.border_bp(i, j)) + 1), j - 1);
        break;
      case NT::WMBW:
        from_right(out, NT::WI, j, -1);
        prefix(NT::WMBP, i, 0, j - 1);
        coverage(out, i + 1, j - 1);
        break;
      case NT::WMBP:
        from_right(out, NT::VP, j);
        prefix(NT::WMBP, i, 1, j - 2);
        prefix(NT::WMBW, i, 1, j - 2);
        for (const auto &be : column(NT::BE, ctx.pair_partner(i))) if (be.key.i == i) {
          out.push_back(be.key.ip + 1); // empty-WI completion
          best_left_hint(out, NT::WI, be.key.ip + 1, 1);
        }
        coverage(out, i + 1, j - 1);
        break;
      case NT::WMB:
        prefix(NT::WMBP, i, 0, j - 1);
        for (const auto &be : column(NT::BE, j)) {
          out.push_back(be.key.jp - 1); // empty-WI completion
          best_right_hint(out, NT::WI, be.key.jp - 1, -1);
        }
        coverage(out, std::max(i + 1, ctx.pair_partner(j) + 1), j - 1);
        break;
      case NT::BE:
        if (a.ip > i && ctx.pair_partner(i) == j && ctx.pair_partner(a.ip) == a.jp)
          out.push_back(scaffold.child_below(i, a.ip));
        break;
      default: break; // fixed-coordinate and MAXLOOP-bounded families
    }
    work.proposals += out.size();
    work.max_proposals = std::max(work.max_proposals, out.size());
    if (out.size() > 64 * b + 32) throw std::runtime_error("per-parent proposal budget exceeded");
    position_bits.sort_unique(out);
  }
  void check_budget() const {
    const auto q = 2 * b + 1;
    if (work.index_records > families * work.guide_records ||
        work.comparisons > work.index_records * q ||
        work.build_moves > work.index_records * (q + 1) ||
        work.undo_moves > work.index_records * (q + 1) ||
        work.cursor_steps > work.index_records)
      throw std::runtime_error("prefix structural work budget violated");
  }
};

using Signature = std::tuple<sc::grammar::RuleId, int, int, int, int>;
inline Signature signature(const sc::deductions::Deduction &d) {
  return {d.rule, d.split.k, d.split.l, d.split.p, d.split.q};
}
struct Result { Outer outer; Work work; };
} // namespace outer_detail

// Optimize=false preserves the original product-of-siblings evaluation as
// an independent numerical reference for the constant-factor changes below.
// TrackLinks=false validates a predecessor when its child is finalized, then
// releases the link with that pending pool. Higher groups are already final
// because every offer requires target < g. TrackLinks=true retains the original
// linked guide and post-pass validation for observations/differential tests.
template <class Context, bool Optimize = true, bool TrackLinks = true>
outer_detail::Result bounded_outer(const outer_detail::Chart &pilot, const outer_detail::Item &root,
    Context &ctx, std::size_t proposal_width, std::size_t frontier_width,
    bool validate = false, bool profile = false, bool numeric_lookup_index = false) {
  using namespace outer_detail;
  const auto start = Clock::now();
  Result result;
  auto &out = result.outer;
  if (numeric_lookup_index) out.values = Chart(zero, nullptr, 0, root.j);
  Proposals proposals(pilot, root.j, proposal_width, ctx, result.work);
  constexpr auto families = outer_probe_families;
  std::array<int, families> phase{};
  int ordinal = 0;
  for (auto nt : sc::engine::basic_exact_span_family_order())
    phase.at(static_cast<std::size_t>(nt)) = ordinal++;
  const auto group = [&](const Item &key) {
    if (key.i < 1 || key.j > root.j || key.j < key.i)
      throw std::runtime_error("invalid bounded target span");
    return (key.j - key.i) * families + phase.at(static_cast<std::size_t>(key.nonterminal));
  };
  std::vector<MaxPool> pools(families * root.j);
  pools.at(group(root)).offer(root, 0, {root, 0}, true, frontier_width);
  std::size_t pending = 1;
  out.peak_pending = out.peak_guide_states = 1;
  sc::engine::AllLocalWeightModel<Context> weights(&ctx);
  sc::engine::RuntimeDeductionChildItem child_item;
  // The shared generator only reads these positions during this parent.
  // Reuse O(p) scratch capacity instead of allocating a new vector each time.
  std::vector<int> positions;
  std::vector<sc::deductions::Deduction> ds;
  std::size_t validated_links = 0;
  if (profile) result.work.preparation_seconds = seconds(start);
  for (int g = static_cast<int>(pools.size()) - 1; g >= 0; --g) {
    auto parents = std::move(pools[g].entries);
    pools[g].novel.clear();
    pending -= parents.size();
    for (const auto &[parent, entry] : parents) {
      out.values.set(parent, entry.value);
      if (!sc::deductions::ItemKeyEq{}(parent, root)) {
        if constexpr (TrackLinks) out.links.emplace(parent, entry.link);
        else {
          const auto &link = entry.link;
          const double parent_value = out.values.get(link.parent);
          if (group(link.parent) <= group(parent) || parent_value == zero ||
              std::fabs(parent_value + link.edge - entry.value) > 1e-9)
            throw std::runtime_error("bounded guide disconnected predecessor");
          ++validated_links;
        }
      }
    }
    for (const auto &record : parents) {
      const auto &parent = record.first;
      const auto &entry = record.second;
      ++out.visited;
      auto stamp = profile ? Clock::now() : start;
      proposals.fill(parent, ctx, positions);
      if (profile) {
        const auto now = Clock::now();
        result.work.proposal_seconds += std::chrono::duration<double>(now - stamp).count();
        stamp = now;
      }
      ds.clear();
      sc::deductions::append_exact_basic_deductions_for_item(ds, parent, ctx, &positions);
      if (profile) {
        const auto now = Clock::now();
        const double elapsed = std::chrono::duration<double>(now - stamp).count();
        result.work.generation_seconds += elapsed;
        result.work.family_generation.at(static_cast<std::size_t>(parent.nonterminal)) += elapsed;
        stamp = now;
      }
      if (ds.size() > 128 * (proposal_width + 1) + 2 * (MAXLOOP + 1) * (MAXLOOP + 1) + 16)
        throw std::runtime_error("bounded per-parent deduction budget violated");
      if (validate) {
        const auto check_start = Clock::now();
        std::vector<Signature> got, want;
        for (const auto &d : ds) got.push_back(signature(d));
        for (const auto &d : sc::deductions::generate_exact_basic_deductions_for_item(parent, ctx)) {
          const auto coordinate = sc::deductions::exact_linear_split_coordinate(d);
          if (!coordinate || std::binary_search(positions.begin(), positions.end(), *coordinate))
            want.push_back(signature(d));
        }
        if (got != want) throw std::runtime_error("bounded/full selected deduction mismatch");
        ++result.work.checks;
        result.work.validation_seconds += seconds(check_start);
        if (profile) stamp = Clock::now();
      }
      out.generated += ds.size();
      std::size_t count = 0;
      const auto offer = [&](const sc::deductions::Deduction &d,
                             const Item &child, double edge, bool mandatory) {
        if (edge == zero) return;
        const auto target = group(child);
        if (!std::isfinite(edge) || target >= static_cast<std::size_t>(g))
          throw std::runtime_error("bounded guide weight/order violation");
        ++count; ++out.eligible; ++out.selected;
        ++out.eligible_rules.at(static_cast<std::size_t>(d.rule));
        auto &pool = pools[target];
        const auto old = pool.entries.size();
        pool.offer(child, entry.value + edge, {parent, edge}, mandatory, frontier_width, Optimize);
        pending -= old; pending += pool.entries.size();
        out.peak_pending = std::max(out.peak_pending, pending);
        out.peak_guide_states = std::max(out.peak_guide_states, out.values.values().size() + pending);
      };
      for (const auto &d : ds) {
        ++out.generated_rules.at(static_cast<std::size_t>(d.rule));
        sc::engine::require_legacy_span_schedule_support(d);
        if constexpr (Optimize) {
          if (!d.child_count) { ++result.work.leaf_deductions; continue; }
          constexpr auto arity = std::tuple_size<decltype(d.children)>::value;
          std::array<Item, arity> children;
          std::array<double, arity> alphas;
          std::size_t zeros = 0;
          for (std::size_t c = 0; c < d.child_count; ++c) {
            children[c] = child_item(d, c);
            alphas[c] = proposals.get_alpha(children[c]);
            ++result.work.child_alpha_lookups;
            if (alphas[c] == zero) ++zeros;
          }
          // With two missing children, every one-hole context has a zero
          // sibling. One missing child is allowed: it is the hole itself.
          if (zeros > 1) { ++result.work.zero_sibling_deductions; continue; }
          ++result.work.weight_evaluations;
          const double local = weights.log_weight_for(d);
          if (local == zero) continue;
          for (std::size_t c = 0; c < d.child_count; ++c) {
            if (zeros == 1 && alphas[c] != zero) continue;
            double edge = local;
            // Keep the original addition order, including rounding.
            for (std::size_t k = 0; k < d.child_count; ++k)
              if (c != k) edge += alphas[k];
            offer(d, children[c], edge, alphas[c] != zero);
          }
        } else {
          ++result.work.weight_evaluations;
          const double local = weights.log_weight_for(d);
          if (local == zero) continue;
          for (std::size_t c = 0; c < d.child_count; ++c) {
            double edge = local;
            for (std::size_t k = 0; k < d.child_count; ++k)
              if (c != k) {
                ++result.work.child_alpha_lookups;
                edge += proposals.get_alpha(child_item(d, k));
              }
            if (edge == zero) continue;
            const auto child = child_item(d, c);
            ++result.work.child_alpha_lookups;
            offer(d, child, edge, proposals.get_alpha(child) != zero);
          }
        }
      }
      out.max_parent_candidates = std::max(out.max_parent_candidates, count);
      if (profile) {
        const double elapsed = seconds(stamp);
        result.work.propagation_seconds += elapsed;
        result.work.family_propagation.at(static_cast<std::size_t>(parent.nonterminal)) += elapsed;
      }
    }
  }
  const auto final_start = profile ? Clock::now() : start;
  for (const auto &pool : pools) out.evicted += pool.evicted;
  proposals.check_budget();
  const auto bound = result.work.guide_records + families * root.j * frontier_width;
  const auto predecessor_count = TrackLinks ? out.links.size() : validated_links;
  if (pending || out.peak_guide_states > bound || predecessor_count + 1 != out.values.values().size())
    throw std::runtime_error("bounded guide state/connectivity budget violated");
  if constexpr (TrackLinks) for (const auto &[child, link] : out.links) {
    const double parent_value = out.values.get(link.parent);
    if (group(link.parent) <= group(child) || parent_value == zero ||
        std::fabs(parent_value + link.edge - out.values.get(child)) > 1e-9)
      throw std::runtime_error("bounded guide disconnected predecessor");
  }
  out.elapsed = seconds(start) - result.work.validation_seconds;
  if (profile) result.work.finalization_seconds = seconds(final_start);
  return result;
}
} // namespace cparty::scfg2::beam
#endif

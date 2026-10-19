#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_NONLEFT_PREFIX_INTERVAL_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_NONLEFT_PREFIX_INTERVAL_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/vp_case3_nonleft_last_child_index.hh"
#include "scfg2/beam/prefix_agenda_storage.hh"
#include "scfg2/beam/bounded_prefix_pending.hh"
#include "scfg2/semirings/semiring.hh"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <set>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase3NonleftPrefixRecord {
  VpCase3NonleftPrefixGeometry geometry{};
  double priority = 0.0;
};

struct VpCase3NonleftPrefixRank {
  bool operator()(const VpCase3NonleftPrefixRecord &lhs,
                  const VpCase3NonleftPrefixRecord &rhs) const {
    return lhs.priority != rhs.priority ? lhs.priority > rhs.priority : lhs.geometry < rhs.geometry;
  }
};

struct VpCase3NonleftPrefixIntervalStats {
  std::size_t middle_queries = 0;
  std::size_t boundary_probes = 0;
  std::size_t prefix_candidates = 0;
  std::size_t active_prefix_probes = 0;
  std::size_t emitted_sources = 0;
  std::size_t duplicate_prefixes = 0;
  std::size_t duplicate_sources = 0;
  std::size_t border_mismatches = 0;
  std::size_t max_active_prefixes = 0;
};

// Exact interval agenda for the non-left CASE3 normalization
// boundary * middle-WI -> prefix, prefix * right-WI -> VP.  Construction of
// prefix geometry is child-driven; final parents are recovered only from
// prefixes active in the current right-sibling scaffold interval.
class VpCase3NonleftPrefixIntervalAgenda {
 public:
  explicit VpCase3NonleftPrefixIntervalAgenda(
      std::size_t width = std::numeric_limits<std::size_t>::max(),
      grammar::NonTerminal carrier = grammar::NonTerminal::VP)
      : width_(width), carrier_(carrier) {}

  template <class Semiring, class Context>
  void observe_middle_wi(
      const deductions::ItemKey &middle,
      typename Semiring::value_type value,
      const VpCase3NonleftChildEventIndex &index,
      const Context &ctx) {
    if (Semiring::is_zero(value)) return;
    ++stats_.middle_queries;
    const auto geometries = index.query_middle_prefixes(middle, ctx, width_);
    const bool duplicate = !geometries.empty() && !registered_wi_.insert(middle).second;
    register_prefixes(geometries, Semiring::priority(value), duplicate);
  }

  template <class Context>
  void observe_middle_wi(
      const deductions::ItemKey &middle,
      double value,
      const VpCase3NonleftChildEventIndex &index,
      const Context &ctx) {
    observe_middle_wi<semirings::SumProductSemiring>(
        middle, value, index, ctx);
  }

  template <class Semiring, class Context>
  void prepare_column(cand_pos_t j,
                      const VpCase3NonleftChildEventIndex &index,
                      Context &ctx) {
    while (!frames_.empty() && frames_.back().expiry == j) {
      selected_active_ = std::move(frames_.back().previous);
      frames_.pop_back();
    }

    if (static_cast<std::size_t>(j) >= registered_empty_.size()) registered_empty_.resize(j + 1);
    const bool duplicate_empty = registered_empty_[j] != 0;
    registered_empty_[j] = 1;
    register_prefixes(
        index.query_empty_middle_prefixes(j, ctx, width_),
        Semiring::priority(Semiring::one()), duplicate_empty);
    auto incoming = pending_.take(j);
    if (!incoming.records.empty()) {
      auto candidates = std::move(incoming.records);
      const cand_pos_t expiry = candidates.front().geometry.target_close + 1;
      frames_.push_back(ActiveFrame{expiry, selected_active_});
      candidates.insert(candidates.end(),
                        selected_active_.begin(),
                        selected_active_.end());
      sort_records(candidates);
      const std::size_t selected = std::min(width_, candidates.size());
      candidates.resize(selected);
      selected_active_ = std::move(candidates);
    }
    stats_.max_active_prefixes =
        std::max(stats_.max_active_prefixes, selected_active_.size());

    sources_.clear();
    std::set<std::pair<cand_pos_t, cand_pos_t>> seen;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    for (const auto &record : selected_active_) {
      const auto &geometry = record.geometry;
      ++stats_.active_prefix_probes;
      const deductions::ItemKey parent{
          carrier_, geometry.i, j, -1, -1};
      for (auto source :
           deductions::generate_exact_vp_case3_deductions_for_item(
               parent, ctx)) {
        engine::require_legacy_span_schedule_support(source);
        const auto partition = vp_case3_partition(source);
        if (!partition.has_value() ||
            *partition == VpCase3Partition::LeftNonempty ||
            weights::is_zero(weight.local_weight_for(source))) {
          continue;
        }
        if (source.split.p != geometry.p ||
            source.split.q != geometry.q ||
            source.split.k != geometry.k) {
          ++stats_.border_mismatches;
          continue;
        }
        if (!seen.emplace(source.parent.i, source.parent.j).second) {
          ++stats_.duplicate_sources;
          continue;
        }
        sources_[source.parent].push_back(std::move(source));
        ++stats_.emitted_sources;
      }
    }
  }

  template <class Context>
  void prepare_column(cand_pos_t j,
                      const VpCase3NonleftChildEventIndex &index,
                      Context &ctx) {
    prepare_column<semirings::SumProductSemiring>(j, index, ctx);
  }

  std::vector<deductions::Deduction> take_for(
      const deductions::ItemKey &parent) {
    const auto found = sources_.find(parent);
    if (found == sources_.end()) return {};
    auto out = std::move(found->second);
    sources_.erase(found);
    return out;
  }

  std::vector<cand_pos_t> parent_starts() const {
    std::vector<cand_pos_t> out;
    out.reserve(sources_.size());
    for (const auto &[parent, sources] : sources_) {
      if (!sources.empty()) out.push_back(parent.i);
    }
    return out;
  }

  const VpCase3NonleftPrefixIntervalStats &stats() const { return stats_; }

  const std::vector<VpCase3NonleftPrefixRecord> &active_prefixes() const {
    return selected_active_;
  }
  PrefixAgendaStorage storage_snapshot() const {
    PrefixAgendaStorage out;
    pending_.append_storage(out);
    out.dedup_keys = registered_wi_.size() +
        std::count(registered_empty_.begin(), registered_empty_.end(), 1);
    out.active_records = selected_active_.size();
    out.saved_frames = frames_.size();
    for (const auto &frame : frames_) out.saved_frame_records += frame.previous.size();
    return out;
  }
  std::size_t pending_comparisons() const { return pending_.comparisons(); }

  ContainerMemory memory_snapshot(const char *name) const {
    auto out = pending_.memory_snapshot(name);
    out.add(vector_memory("", frames_));
    for (const auto &frame : frames_) out.add(vector_memory("", frame.previous));
    out.add(vector_memory("", selected_active_));
    out.add(hash_memory("", registered_wi_));
    out.add(vector_memory("", registered_empty_));
    out.add(hash_vector_memory("", sources_));
    return out;
  }

 private:
  void register_prefixes(
      const std::vector<VpCase3NonleftPrefixGeometry> &geometries,
      double priority, bool duplicate) {
    stats_.boundary_probes += geometries.size();
    std::vector<VpCase3NonleftPrefixGeometry> selected = geometries;
    std::sort(selected.begin(), selected.end());
    const std::size_t selected_count = std::min(width_, selected.size());
    selected.resize(selected_count);
    for (const auto &geometry : selected) {
      ++stats_.prefix_candidates;
      if (duplicate) {
        ++stats_.duplicate_prefixes;
        continue;
      }
      pending_.insert(geometry.k, VpCase3NonleftPrefixRecord{geometry, priority}, width_);
    }
  }

  static void sort_records(std::vector<VpCase3NonleftPrefixRecord> &records) {
    std::sort(records.begin(), records.end(), VpCase3NonleftPrefixRank{});
  }

  struct ActiveFrame {
    cand_pos_t expiry = -1;
    std::vector<VpCase3NonleftPrefixRecord> previous;
  };

  BoundedPrefixPending<VpCase3NonleftPrefixRecord, VpCase3NonleftPrefixRank> pending_;
  std::vector<ActiveFrame> frames_;
  std::vector<VpCase3NonleftPrefixRecord> selected_active_;
  // q,k identify the one nonempty middle WI. Empty middle has q=k-1.
  std::unordered_set<deductions::ItemKey, deductions::ItemKeyHash, deductions::ItemKeyEq> registered_wi_;
  std::vector<unsigned char> registered_empty_;
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      sources_;
  VpCase3NonleftPrefixIntervalStats stats_{};
  std::size_t width_ = std::numeric_limits<std::size_t>::max();
  grammar::NonTerminal carrier_ = grammar::NonTerminal::VP;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_NONLEFT_PREFIX_INTERVAL_HH

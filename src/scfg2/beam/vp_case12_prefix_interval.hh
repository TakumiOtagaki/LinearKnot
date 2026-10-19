#ifndef CPARTY_SCFG2_BEAM_VP_CASE12_PREFIX_INTERVAL_HH
#define CPARTY_SCFG2_BEAM_VP_CASE12_PREFIX_INTERVAL_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/vp_case12_last_child_index.hh"
#include "scfg2/beam/prefix_agenda_storage.hh"
#include "scfg2/beam/bounded_prefix_pending.hh"
#include "scfg2/semirings/semiring.hh"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase12PrefixRecord {
  VpCase12PrefixGeometry geometry{};
  double priority = 0.0;
};

struct VpCase12PrefixRank {
  bool operator()(const VpCase12PrefixRecord &lhs, const VpCase12PrefixRecord &rhs) const {
    return lhs.priority != rhs.priority ? lhs.priority > rhs.priority : lhs.geometry < rhs.geometry;
  }
};

struct VpCase12PrefixIntervalStats {
  std::size_t prefix_candidates = 0;
  std::size_t emitted_sources = 0;
  std::size_t duplicate_prefixes = 0;
  std::size_t duplicate_sources = 0;
  std::size_t border_mismatches = 0;
  std::size_t max_case1_group = 0;
  std::size_t max_case2_active = 0;
};

// CASE1 keeps a bounded left-prefix bucket per q.  At final column j it joins
// only q values exposed by retained WI(q+1,j-1), plus the empty-right q=j-1.
// CASE2 has a common right factor for fixed j and can therefore keep one
// top-width active set across its nested right-child band intervals.
class VpCase12PrefixIntervalAgenda {
 public:
  explicit VpCase12PrefixIntervalAgenda(
      std::size_t width = std::numeric_limits<std::size_t>::max(),
      grammar::NonTerminal carrier = grammar::NonTerminal::VP)
      : width_(width), carrier_(carrier) {}

  template <class Semiring, class Context>
  void observe_left_wi(
      const deductions::ItemKey &left,
      typename Semiring::value_type value,
      const VpCase12DirectNonemptyLeftQuery &query,
      Context &ctx) {
    if (Semiring::is_zero(value)) return;
    const auto geometries = query.query_completed_left(left, ctx, width_);
    const bool duplicate = !geometries.empty() && !registered_wi_.insert(left).second;
    register_prefixes(geometries, Semiring::priority(value), duplicate);
  }

  template <class Context>
  void observe_left_wi(
      const deductions::ItemKey &left,
      double value,
      const VpCase12DirectNonemptyLeftQuery &query,
      Context &ctx) {
    observe_left_wi<semirings::SumProductSemiring>(left, value, query, ctx);
  }

  template <class Semiring, class Context>
  void prepare_column(
      cand_pos_t j,
      const std::vector<cand_pos_t> &retained_right_starts,
      const VpCase12DirectEmptyLeftQuery &empty_query,
      Context &ctx) {
    expire_case1(j);
    while (!case2_frames_.empty() && case2_frames_.back().expiry == j) {
      case2_active_ = std::move(case2_frames_.back().previous);
      case2_frames_.pop_back();
    }

    if (static_cast<std::size_t>(j) >= registered_empty_.size()) registered_empty_.resize(j + 1);
    const bool duplicate_empty = registered_empty_[j] != 0;
    registered_empty_[j] = 1;
    register_prefixes(empty_query.query_activation(j, ctx),
                      Semiring::priority(Semiring::one()), duplicate_empty);
    activate_case1(j);
    activate_case2(j);

    stats_.max_case2_active =
        std::max(stats_.max_case2_active, case2_active_.size());
    sources_.clear();
    std::set<std::pair<grammar::RuleId, cand_pos_t>> seen_parent;
    engine::AllLocalWeightModel<Context> weight(&ctx);

    std::set<cand_pos_t> final_qs;
    if (j > 1) final_qs.insert(j - 1);
    for (const cand_pos_t start : retained_right_starts) {
      if (start >= 1 && start <= j - 1) final_qs.insert(start - 1);
    }
    for (const cand_pos_t q : final_qs) {
      if (q < 0 || static_cast<std::size_t>(q) >= case1_active_.size()) continue;
      for (const auto &record : case1_active_[q]) {
        emit(record, j, ctx, weight, seen_parent);
      }
    }

    const cand_pos_t case2_q = ctx.parent_index(j);
    const bool case2_right_available = final_qs.count(case2_q) != 0;
    if (case2_right_available) {
      for (const auto &record : case2_active_) {
        emit(record, j, ctx, weight, seen_parent);
      }
    }
  }

  template <class Context>
  void prepare_column(
      cand_pos_t j,
      const std::vector<cand_pos_t> &retained_right_starts,
      const VpCase12DirectEmptyLeftQuery &empty_query,
      Context &ctx) {
    prepare_column<semirings::SumProductSemiring>(
        j, retained_right_starts, empty_query, ctx);
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

  const VpCase12PrefixIntervalStats &stats() const { return stats_; }

  std::vector<VpCase12PrefixRecord> active_case1_prefixes() const {
    std::vector<VpCase12PrefixRecord> out;
    for (const auto &records : case1_active_) {
      out.insert(out.end(), records.begin(), records.end());
    }
    return out;
  }
  const std::vector<VpCase12PrefixRecord> &active_case2_prefixes() const {
    return case2_active_;
  }
  PrefixAgendaStorage storage_snapshot() const {
    PrefixAgendaStorage out;
    case1_pending_.append_storage(out);
    case2_pending_.append_storage(out);
    out.dedup_keys = registered_wi_.size() +
        std::count(registered_empty_.begin(), registered_empty_.end(), 1);
    out.active_records = case2_active_.size();
    for (const auto &records : case1_active_) {
      out.active_records += records.size();
    }
    out.saved_frames = case2_frames_.size();
    for (const auto &frame : case2_frames_) out.saved_frame_records += frame.previous.size();
    return out;
  }
  std::size_t pending_comparisons() const {
    return case1_pending_.comparisons() + case2_pending_.comparisons();
  }

  ContainerMemory memory_snapshot(const char *name) const {
    auto out = case1_pending_.memory_snapshot(name);
    out.add(case2_pending_.memory_snapshot(""));
    out.add(nested_vector_memory("", case1_active_));
    out.add(nested_vector_memory("", case1_expiry_));
    out.add(vector_memory("", case2_frames_));
    for (const auto &frame : case2_frames_) out.add(vector_memory("", frame.previous));
    out.add(vector_memory("", case2_active_));
    out.add(hash_memory("", registered_wi_));
    out.add(vector_memory("", registered_empty_));
    out.add(hash_vector_memory("", sources_));
    return out;
  }

 private:
  void register_prefixes(
      const std::vector<VpCase12PrefixGeometry> &geometries,
      double priority, bool duplicate) {
    for (const auto &geometry : geometries) {
      ++stats_.prefix_candidates;
      if (duplicate) {
        ++stats_.duplicate_prefixes;
        continue;
      }
      auto &pending = geometry.rule == grammar::RuleId::VP_WI_CASE1
          ? case1_pending_ : case2_pending_;
      pending.insert(geometry.target_begin, VpCase12PrefixRecord{geometry, priority}, width_);
    }
  }

  static void sort_records(std::vector<VpCase12PrefixRecord> &records) {
    std::sort(records.begin(), records.end(), VpCase12PrefixRank{});
  }

  void activate_case1(cand_pos_t j) {
    auto incoming = case1_pending_.take(j);
    if (incoming.records.empty()) return;
    // Every CASE1 geometry has target_begin=q+1: exactly one group per column.
    const cand_pos_t q = j - 1;
    const cand_pos_t expiry = incoming.records.front().geometry.target_end + 1;
    if (static_cast<std::size_t>(q) >= case1_active_.size()) case1_active_.resize(q + 1);
    if (static_cast<std::size_t>(expiry) >= case1_expiry_.size()) case1_expiry_.resize(expiry + 1);
    case1_active_[q] = std::move(incoming.records);
    case1_expiry_[expiry].push_back(q);
    stats_.max_case1_group = std::max(stats_.max_case1_group, case1_active_[q].size());
  }

  void activate_case2(cand_pos_t j) {
    auto incoming = case2_pending_.take(j);
    auto candidates = std::move(incoming.records);
    if (!candidates.empty()) {
      const cand_pos_t expiry = candidates.front().geometry.target_end + 1;
      case2_frames_.push_back(Case2Frame{expiry, case2_active_});
      candidates.insert(candidates.end(),
                        case2_active_.begin(),
                        case2_active_.end());
      sort_records(candidates);
      const std::size_t selected = std::min(width_, candidates.size());
      candidates.resize(selected);
      case2_active_ = std::move(candidates);
    }
  }

  void expire_case1(cand_pos_t j) {
    if (j < 0 || static_cast<std::size_t>(j) >= case1_expiry_.size()) return;
    for (const cand_pos_t q : case1_expiry_[j]) case1_active_[q].clear();
    case1_expiry_[j].clear();
  }

  template <class Context>
  void emit(const VpCase12PrefixRecord &record,
            cand_pos_t j,
            Context &ctx,
            engine::AllLocalWeightModel<Context> &weight,
            std::set<std::pair<grammar::RuleId, cand_pos_t>> &seen_parent) {
    const auto &geometry = record.geometry;
    const deductions::ItemKey parent{
        carrier_, geometry.i, j, -1, -1};
    for (auto source : deductions::generate_exact_basic_deductions(
             parent, geometry.rule, ctx)) {
      engine::require_legacy_span_schedule_support(source);
      if (weights::is_zero(weight.local_weight_for(source))) {
        continue;
      }
      const bool matches = source.split.p == geometry.p &&
          (geometry.rule == grammar::RuleId::VP_WI_CASE2 ||
           source.split.q == geometry.q);
      if (!matches) {
        ++stats_.border_mismatches;
        continue;
      }
      if (!seen_parent.emplace(source.rule, source.parent.i).second) {
        ++stats_.duplicate_sources;
        continue;
      }
      sources_[source.parent].push_back(std::move(source));
      ++stats_.emitted_sources;
    }
  }

  struct Case2Frame {
    cand_pos_t expiry = -1;
    std::vector<VpCase12PrefixRecord> previous;
  };

  BoundedPrefixPending<VpCase12PrefixRecord, VpCase12PrefixRank> case1_pending_, case2_pending_;
  std::vector<std::vector<VpCase12PrefixRecord>> case1_active_;
  std::vector<std::vector<cand_pos_t>> case1_expiry_;
  std::vector<Case2Frame> case2_frames_;
  std::vector<VpCase12PrefixRecord> case2_active_;
  // A geometry uniquely identifies its nonempty WI trigger. Empty-trigger
  // geometries are disjoint and registered once per activation column.
  std::unordered_set<deductions::ItemKey, deductions::ItemKeyHash, deductions::ItemKeyEq> registered_wi_;
  std::vector<unsigned char> registered_empty_;
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      sources_;
  VpCase12PrefixIntervalStats stats_{};
  std::size_t width_ = std::numeric_limits<std::size_t>::max();
  grammar::NonTerminal carrier_ = grammar::NonTerminal::VP;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE12_PREFIX_INTERVAL_HH

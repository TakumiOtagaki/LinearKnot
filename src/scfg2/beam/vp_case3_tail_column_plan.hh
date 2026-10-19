#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_COLUMN_PLAN_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_COLUMN_PLAN_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/vp_case3_tail_cache.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <array>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase3PlannedSource {
  VpCase3TailKey tail{};
  deductions::Deduction deduction{};
};

class ExactVpCase3TailColumnPlan {
 public:
  static ExactVpCase3TailColumnPlan from_sources(
      std::vector<VpCase3PlannedSource> sources) {
    ExactVpCase3TailColumnPlan plan;
    std::unordered_set<VpCase3TailKey, VpCase3TailKeyHash> seen;
    for (auto &source : sources) {
      if (seen.insert(source.tail).second) {
        plan.columns_[source.tail.j].push_back(source.tail);
      }
      plan.sources_.push_back(std::move(source));
    }
    return plan;
  }

  template <class Context>
  static ExactVpCase3TailColumnPlan discover(cand_pos_t n,
                                             const Context &ctx) {
    ExactVpCase3TailColumnPlan plan;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    std::unordered_set<VpCase3TailKey,
                       VpCase3TailKeyHash>
        seen;
    for (cand_pos_t j = 1; j <= n; ++j) {
      for (cand_pos_t i = j; i >= 1; --i) {
        for (const auto nt : {grammar::NonTerminal::VP}) {
          const deductions::ItemKey parent{nt, i, j, -1, -1};
          const auto deductions =
              deductions::generate_exact_vp_case3_deductions_for_item(
                  parent, ctx);
          for (const auto &deduction : deductions) {
            engine::require_legacy_span_schedule_support(deduction);
            if (vp_case3_partition(deduction) !=
                    VpCase3Partition::LeftNonempty ||
                weights::is_zero(weight.local_weight_for(deduction))) {
              continue;
            }
            const VpCase3TailKey key{j,
                                     deduction.split.q,
                                     deduction.split.k,
                                     deduction.split.l};
            plan.sources_.push_back(VpCase3PlannedSource{key, deduction});
            if (seen.insert(key).second) plan.columns_[j].push_back(key);
          }
        }
      }
    }
    return plan;
  }

  template <class Semiring, class Chart>
  std::size_t materialize_column(cand_pos_t j,
                                 const Chart &chart,
                                 ExactVpCase3TailCache &cache) const {
    const auto it = columns_.find(j);
    if (it == columns_.end()) return 0;
    std::size_t count = 0;
    for (const auto &key : it->second) {
      count += cache.template materialize_tail<Semiring>(key, chart) ? 1 : 0;
    }
    return count;
  }

  std::size_t tail_count() const {
    std::size_t count = 0;
    for (const auto &[j, keys] : columns_) {
      (void)j;
      count += keys.size();
    }
    return count;
  }

  std::size_t source_count() const { return sources_.size(); }
  std::size_t column_count() const { return columns_.size(); }
  const auto &sources() const { return sources_; }

 private:
  std::unordered_map<cand_pos_t, std::vector<VpCase3TailKey>> columns_;
  std::vector<VpCase3PlannedSource> sources_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_TAIL_COLUMN_PLAN_HH

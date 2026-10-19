#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_WI_FRONTIER_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_WI_FRONTIER_HH

#include "scfg2/beam/right_boundary_item_frontier.hh"
#include "scfg2/beam/vp_case3_left_child_index.hh"
#include "scfg2/beam/vp_case3_right_endpoint_producer_registry.hh"

#include <cstddef>
#include <limits>

namespace cparty::scfg2::beam {

enum class VpCase3WiFrontierScope {
  EligibleCase3Producers,
  AllWiItems,
};

// Producer frontier. Candidates are bucketed by the completed WI right
// boundary and are not queried until the next target column begins.
// KeepAll is represented by an unbounded width; finite width suppresses the
// child query entirely and therefore has no direct-CASE3 fallback.
class VpCase3WiProducerFrontier {
 public:
  static VpCase3WiProducerFrontier keep_all(
      VpCase3WiFrontierScope scope =
          VpCase3WiFrontierScope::EligibleCase3Producers) {
    return VpCase3WiProducerFrontier(
        std::numeric_limits<std::size_t>::max(), scope);
  }

  explicit VpCase3WiProducerFrontier(
      std::size_t width,
      VpCase3WiFrontierScope scope =
          VpCase3WiFrontierScope::EligibleCase3Producers)
      : frontier_(width), scope_(scope) {}

  template <class Semiring, class Chart, class Context>
  void observe_completed_child(
      const deductions::ItemKey &item,
      const Chart &chart,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    if (item.nonterminal != grammar::NonTerminal::WI) return;
    if (scope_ == VpCase3WiFrontierScope::EligibleCase3Producers &&
        !index.recognizes_completed_wi(item, ctx)) {
      return;
    }
    frontier_.observe(item, Semiring::priority(chart.get(item)));
  }

  template <class Chart, class Context>
  void observe_completed_child(
      const deductions::ItemKey &item,
      const Chart &chart,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    observe_completed_child<semirings::SumProductSemiring>(
        item, chart, index, ctx);
  }

  template <class Chart, class Context>
  void flush_before(
      cand_pos_t target_j,
      const Chart &chart,
      PendingVpCase3LeftCandidates &pending,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    frontier_.flush_before(target_j, [&](const deductions::ItemKey &item) {
      pending.observe_completed_child(
          item, chart, index, ctx);
    });
  }

  template <class Semiring, class Chart, class Context>
  void flush_before_registering(
      cand_pos_t target_j,
      const Chart &chart,
      TopWidthVpCase3PrecedingSiblingProducerFrontier &registry,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    frontier_.flush_before(target_j, [&](const deductions::ItemKey &item) {
      registry.observe_completed_wi(
          item, Semiring::priority(chart.get(item)), index, ctx);
    });
  }

  template <class Chart, class Context>
  void flush_before_registering(
      cand_pos_t target_j,
      const Chart &chart,
      TopWidthVpCase3PrecedingSiblingProducerFrontier &registry,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    flush_before_registering<semirings::SumProductSemiring>(
        target_j, chart, registry, index, ctx);
  }

  std::size_t candidate_count() const { return frontier_.candidate_count(); }
  std::size_t selected_count() const { return frontier_.selected_count(); }
  std::size_t pruned_count() const { return frontier_.pruned_count(); }
  bool respects_width() const { return frontier_.respects_width(); }

 private:
  RightBoundaryTopWidthItemFrontier frontier_;
  VpCase3WiFrontierScope scope_ =
      VpCase3WiFrontierScope::EligibleCase3Producers;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_WI_FRONTIER_HH

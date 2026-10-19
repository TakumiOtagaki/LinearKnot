#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_right_endpoint_producer_registry.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace {

using cparty::scfg2::beam::VpCase3RegisteredLeftProducer;
using cparty::scfg2::deductions::ItemKey;
using Identity = std::tuple<cand_pos_t,
                            cand_pos_t,
                            cand_pos_t,
                            cand_pos_t,
                            cand_pos_t,
                            cand_pos_t>;

double priority_for(const ItemKey &item) {
  return 1.0 + static_cast<double>(item.i) * 0.01 +
      static_cast<double>(item.j) * 0.0001;
}

Identity identity(
    const cparty::scfg2::beam::VpCase3RightEndpointProducerMatch &match) {
  return {match.right_path.ancestor,
          match.right_path.child,
          match.producer.item.i,
          match.producer.item.j,
          match.producer.ancestor,
          match.producer.left_child};
}

void check_width(
    std::size_t width,
    cand_pos_t n,
    const std::vector<ItemKey> &schedule,
    const cparty::scfg2::beam::ExactVpCase3BorderEventChildIndex &index,
    const cparty::scfg2::replay::WFinalExactInsideContext &ctx) {
  cparty::scfg2::beam::TopWidthVpCase3PrecedingSiblingProducerFrontier
      frontier(width);
  using ByChild =
      std::map<cand_pos_t, std::vector<VpCase3RegisteredLeftProducer>>;
  std::map<cand_pos_t, ByChild> offline;
  std::map<cand_pos_t, std::vector<ItemKey>> wi_by_right;
  for (const auto &item : schedule) {
    if (item.nonterminal == cparty::scfg2::grammar::NonTerminal::WI) {
      wi_by_right[item.j].push_back(item);
    }
  }

  for (cand_pos_t j = 1; j <= n; ++j) {
    const auto selected = frontier.query_target_endpoint(j, index, ctx);
    std::set<Identity> actual;
    for (const auto &match : selected) actual.insert(identity(match));

    std::set<Identity> expected;
    for (const auto &right_path : index.paths_for_endpoint(j)) {
      std::vector<VpCase3RegisteredLeftProducer> compatible;
      const auto ancestor = offline.find(right_path.ancestor);
      if (ancestor != offline.end()) {
        for (auto child = ancestor->second.begin();
             child != ancestor->second.lower_bound(right_path.child);
             ++child) {
          compatible.insert(
              compatible.end(), child->second.begin(), child->second.end());
        }
      }
      std::sort(
          compatible.begin(),
          compatible.end(),
          cparty::scfg2::beam::vp_case3_registered_producer_rank_before);
      if (compatible.size() > width) compatible.resize(width);
      for (const auto &producer : compatible) {
        expected.emplace(right_path.ancestor,
                         right_path.child,
                         producer.item.i,
                         producer.item.j,
                         producer.ancestor,
                         producer.left_child);
      }
    }
    assert(actual == expected);

    const auto completed = wi_by_right.find(j);
    if (completed == wi_by_right.end()) continue;
    for (const auto &item : completed->second) {
      const auto factor = priority_for(item);
      frontier.observe_completed_wi(item, factor, index, ctx);
      if (!index.recognizes_completed_wi(item, ctx)) continue;
      const cand_pos_t endpoint = item.i - 1;
      for (const auto &path : index.paths_for_endpoint(endpoint)) {
        offline[path.ancestor][path.child].push_back(
            VpCase3RegisteredLeftProducer{item,
                                          factor,
                                          endpoint,
                                          item.j + 1,
                                          path.ancestor,
                                          path.child});
      }
    }
  }
  assert(frontier.late_registration_count() == 0);
  assert(frontier.nonmonotone_query_count() == 0);
  assert(frontier.maximum_prefix_size() <= width);
  assert(frontier.maximum_returned_per_path() <= width);
  // Each static scaffold child is passed by one ancestor's monotone cursor.
  assert(frontier.visited_child_count() <= static_cast<std::size_t>(n));
}

}  // namespace

int main(int argc, char **argv) {
  assert(argc == 1 || argc == 3 || argc == 4);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 4 ? argv[3] : "params/rna_Turner04.par",
         1);
  std::string sequence = argc >= 3
      ? argv[1]
      : "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  std::string scaffold = argc >= 3
      ? argv[2]
      : "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  assert(sequence.size() == scaffold.size());
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto schedule =
      cparty::scfg2::engine::build_pkfree_basic_exact_right_boundary_schedule(
          n, ctx);
  const auto index =
      cparty::scfg2::beam::ExactVpCase3BorderEventChildIndex::build(n, ctx);
  for (const auto width : {std::size_t(1), std::size_t(2), std::size_t(5),
                           std::size_t(10)}) {
    check_width(width, n, schedule, index, ctx);
  }
  std::cout << "scfg2_vp_case3_endpoint_prefix_frontier=ok\n";
  return EXIT_SUCCESS;
}

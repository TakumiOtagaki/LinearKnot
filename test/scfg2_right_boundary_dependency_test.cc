#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/right_boundary_dependency_profile.hh"
#include "scfg2/beam/vp_wmbp_online_schedule.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  std::string sequence =
      "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  std::string scaffold =
      "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule =
      cparty::scfg2::engine::build_pkfree_basic_exact_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);

  cparty::scfg2::beam::RightBoundaryDependencyProfile profile;
  for (const auto &item : schedule) {
    const auto deductions =
        cparty::scfg2::deductions::generate_exact_basic_deductions_for_item(
            item, ctx);
    for (const auto &deduction : deductions) profile.observe(deduction);
  }

  assert(profile.edge_count > 0);
  assert(profile.prior_column_count > 0);
  assert(profile.same_column_inward_count > 0);
  assert(profile.same_span_phase_count > 0);
  assert(profile.future_column_count == 0);
  assert(profile.same_column_outward_count == 0);
  assert(profile.phase_graph_is_acyclic());
  assert(profile.same_column_family_graph_is_acyclic());
  for (std::size_t from = 0;
       from < cparty::scfg2::beam::nonterminal_count;
       ++from) {
    for (std::size_t to = 0;
         to < cparty::scfg2::beam::nonterminal_count;
         ++to) {
      if (from == to || !profile.same_column_family_edges[from][to]) {
        continue;
      }
      const auto child =
          static_cast<cparty::scfg2::grammar::NonTerminal>(from);
      const auto parent =
          static_cast<cparty::scfg2::grammar::NonTerminal>(to);
      assert(cparty::scfg2::beam::vp_wmbp_online_column_phase_rank(child) <
             cparty::scfg2::beam::vp_wmbp_online_column_phase_rank(parent));
    }
  }
  std::cout << "scfg2_right_boundary_dependency=ok edges="
            << profile.edge_count
            << " prior_column=" << profile.prior_column_count
            << " same_column_inward=" << profile.same_column_inward_count
            << " same_span_phase=" << profile.same_span_phase_count << '\n';
  return EXIT_SUCCESS;
}

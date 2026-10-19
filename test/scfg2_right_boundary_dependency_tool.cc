#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/right_boundary_dependency_profile.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    std::cerr << "usage: scfg2_right_boundary_dependency_tool <seq> <G> "
                 "[param_file]\n";
    return 2;
  }
  try {
    std::string sequence = argv[1];
    std::string scaffold = argv[2];
    if (sequence.size() != scaffold.size()) {
      throw std::runtime_error("sequence and scaffold lengths differ");
    }
    if (argc == 4) {
      setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[3], 1);
    } else if (!std::getenv("CPARTY_FIXED_ENERGY_PARAM_FILE")) {
      setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    }
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
    std::cout << "n\tedges\tprior_column\tfuture_column"
                 "\tsame_column_inward\tsame_column_outward"
                 "\tsame_span_phase\tphase_acyclic"
                 "\tcolumn_family_acyclic\n"
              << sequence.size() << '\t' << profile.edge_count << '\t'
              << profile.prior_column_count << '\t'
              << profile.future_column_count << '\t'
              << profile.same_column_inward_count << '\t'
              << profile.same_column_outward_count << '\t'
              << profile.same_span_phase_count << '\t'
              << (profile.phase_graph_is_acyclic() ? 1 : 0) << '\t'
              << (profile.same_column_family_graph_is_acyclic() ? 1 : 0)
              << '\n';
    for (std::size_t from = 0;
         from < cparty::scfg2::beam::nonterminal_count;
         ++from) {
      for (std::size_t to = 0;
           to < cparty::scfg2::beam::nonterminal_count;
           ++to) {
        if (!profile.phase_edges[from][to]) continue;
        std::cerr << "same_span_phase="
                  << cparty::scfg2::grammar::nonterminal_name(
                         static_cast<cparty::scfg2::grammar::NonTerminal>(from))
                  << "->"
                  << cparty::scfg2::grammar::nonterminal_name(
                         static_cast<cparty::scfg2::grammar::NonTerminal>(to))
                  << '\n';
      }
    }
    for (std::size_t from = 0;
         from < cparty::scfg2::beam::nonterminal_count;
         ++from) {
      for (std::size_t to = 0;
           to < cparty::scfg2::beam::nonterminal_count;
           ++to) {
        if (from == to || !profile.same_column_family_edges[from][to]) {
          continue;
        }
        std::cerr << "same_column_family="
                  << cparty::scfg2::grammar::nonterminal_name(
                         static_cast<cparty::scfg2::grammar::NonTerminal>(from))
                  << "->"
                  << cparty::scfg2::grammar::nonterminal_name(
                         static_cast<cparty::scfg2::grammar::NonTerminal>(to))
                  << '\n';
      }
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "scfg2_right_boundary_dependency_tool failed: "
              << error.what() << '\n';
    return 1;
  }
}

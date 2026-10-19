#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/beam/virtual_candidate.hh"

#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>

namespace {

using cparty::scfg2::grammar::RuleId;
using cparty::scfg2::deductions::ItemKey;

constexpr std::array<RuleId, 9> vp_rules{
    RuleId::VP_WI_CASE1,
    RuleId::VP_WI_CASE2,
    RuleId::VP_WI_CASE3,
    RuleId::VP_STACK,
    RuleId::VP_INTERNAL_LOOP,
    RuleId::VP_WIP_VP_LEFT,
    RuleId::VP_VP_WIP_RIGHT,
    RuleId::VP_WIP_VPR,
    RuleId::VP_VPL_WIP,
};

std::size_t rule_index(RuleId rule) {
  for (std::size_t idx = 0; idx < vp_rules.size(); ++idx) {
    if (vp_rules[idx] == rule) return idx;
  }
  return vp_rules.size();
}

template <class Semiring>
int run(const std::string &sequence,
        const std::string &scaffold,
        std::size_t width) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto widths = width == std::numeric_limits<std::size_t>::max()
      ? cparty::scfg2::replay::WFinalVpWmbpOnlineWidths{}
      : cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(width);
  std::array<std::size_t, vp_rules.size()> sources{};
  std::array<std::size_t, 3> case3_partitions{};
  std::array<std::unordered_set<
                 ItemKey,
                 cparty::scfg2::deductions::ItemKeyHash,
                 cparty::scfg2::deductions::ItemKeyEq>,
             vp_rules.size()>
      parents;
  auto result = cparty::scfg2::replay::run_w_final_vp_wmbp_online_observed<
      Semiring>(
      ctx,
      static_cast<cand_pos_t>(sequence.size()),
      widths,
      [&](const ItemKey &item,
          const cparty::scfg2::deductions::Deduction &deduction,
          double candidate,
          const auto &,
          double) {
        if (item.nonterminal !=
                cparty::scfg2::grammar::NonTerminal::VP ||
            candidate == Semiring::zero()) {
          return;
        }
        const auto idx = rule_index(deduction.rule);
        assert(idx < vp_rules.size());
        ++sources[idx];
        parents[idx].insert(item);
        if (deduction.rule == RuleId::VP_WI_CASE3) {
          const auto partition =
              cparty::scfg2::beam::vp_case3_partition(deduction);
          assert(partition.has_value());
          ++case3_partitions[static_cast<std::size_t>(*partition)];
        }
      });
  std::size_t total = 0;
  std::cout << "rule\tpositive_sources\tpositive_parents\n";
  for (std::size_t idx = 0; idx < vp_rules.size(); ++idx) {
    total += sources[idx];
    std::cout << cparty::scfg2::grammar::rule_name(vp_rules[idx]) << '\t'
              << sources[idx] << '\t' << parents[idx].size() << '\n';
  }
  std::cout << "VP_WI_CASE3_BOUNDARY\t" << case3_partitions[0] << "\t0\n"
            << "VP_WI_CASE3_LEFT_NONEMPTY\t" << case3_partitions[1]
            << "\t0\n"
            << "VP_WI_CASE3_OTHER\t" << case3_partitions[2] << "\t0\n";
  std::cout << "TOTAL\t" << total << '\t'
            << result.positive_items[static_cast<std::size_t>(
                   cparty::scfg2::grammar::NonTerminal::VP)]
            << '\n';
  return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3 || argc > 6) {
    std::cerr << "usage: " << argv[0]
              << " SEQUENCE SCAFFOLD [WIDTH|all] [inside|max] [PARAM_FILE]\n";
    return EXIT_FAILURE;
  }
  const std::string sequence = argv[1];
  const std::string scaffold = argv[2];
  const std::string width_arg = argc >= 4 ? argv[3] : "all";
  const std::string mode = argc >= 5 ? argv[4] : "inside";
  const std::size_t width = width_arg == "all"
      ? std::numeric_limits<std::size_t>::max()
      : std::stoull(width_arg);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc >= 6 ? argv[5] : "params/rna_Turner04.par",
         1);
  if (mode == "inside") {
    return run<cparty::scfg2::semirings::SumProductSemiring>(
        sequence, scaffold, width);
  }
  if (mode == "max") {
    return run<cparty::scfg2::semirings::MaxProductSemiring>(
        sequence, scaffold, width);
  }
  return EXIT_FAILURE;
}

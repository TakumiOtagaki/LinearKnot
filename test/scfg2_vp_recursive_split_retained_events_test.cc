#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_recursive_split_retained_events.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

bool close(double lhs, double rhs) {
  const double scale = std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  return std::fabs(lhs - rhs) <= 1e-10 * scale;
}

template <class Semiring>
std::array<std::size_t, 4> check_events(cparty::scfg2::weights::RuntimeSnapshot &owner) {
  using namespace cparty::scfg2;
  replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const auto schedule =
      engine::build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  const auto chart = replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      [](const Deduction &) { return true; },
      [](const ItemKey &, const Deduction &,
         typename Semiring::value_type, const auto &, double) {});
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);

  using Starts = std::vector<std::vector<cand_pos_t>>;
  const auto collect_starts = [&](grammar::NonTerminal nt) {
    Starts out(static_cast<std::size_t>(n + 1));
    for (cand_pos_t j = 1; j <= n; ++j) {
      for (cand_pos_t i = 1; i <= j; ++i) {
        if (chart.get(ItemKey{nt, i, j, -1, -1}) != Semiring::zero()) {
          out[static_cast<std::size_t>(j)].push_back(i);
        }
      }
    }
    return out;
  };
  const auto vp = collect_starts(grammar::NonTerminal::VP);
  const auto wip = collect_starts(grammar::NonTerminal::WIP);
  const auto vpl = collect_starts(grammar::NonTerminal::VPL);
  const auto vpr = collect_starts(grammar::NonTerminal::VPR);
  const auto at = [](const Starts &starts,
                     cand_pos_t right) -> const std::vector<cand_pos_t> & {
    if (right < 0 || static_cast<std::size_t>(right) >= starts.size()) {
      static const std::vector<cand_pos_t> empty;
      return empty;
    }
    return starts[static_cast<std::size_t>(right)];
  };

  using Totals = std::unordered_map<ItemKey,
                                    typename Semiring::value_type,
                                    deductions::ItemKeyHash,
                                    deductions::ItemKeyEq>;
  using SourceId =
      std::tuple<grammar::RuleId, cand_pos_t, cand_pos_t, cand_pos_t>;
  Totals driven;
  std::set<SourceId> driven_ids;
  beam::VpRecursiveSplitRetainedEvents events;
  for (cand_pos_t j = 1; j <= n; ++j) {
    events.prepare_column(
        j,
        [&](cand_pos_t right) -> const auto & { return at(vp, right); },
        [&](cand_pos_t right) -> const auto & { return at(wip, right); },
        [&](cand_pos_t right) -> const auto & { return at(vpl, right); },
        [&](cand_pos_t right) -> const auto & { return at(vpr, right); },
        ctx);
    for (cand_pos_t i = 1; i <= j; ++i) {
      const ItemKey parent{grammar::NonTerminal::VP, i, j, -1, -1};
      for (const auto &source : events.take_for(parent)) {
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, source, weight.weight_for(source));
        assert(term != Semiring::zero());
        assert(driven_ids.emplace(
            source.rule, source.parent.i, source.parent.j, source.split.k)
                   .second);
        auto [found, inserted] = driven.emplace(parent, Semiring::zero());
        found->second = Semiring::plus(found->second, term);
      }
    }
  }
  if (events.stats().duplicate_sources != 0) {
    std::cerr << "recursive VP duplicate sources="
              << events.stats().duplicate_sources << " probes="
              << events.stats().join_probes[0] << ','
              << events.stats().join_probes[1] << ','
              << events.stats().join_probes[2] << ','
              << events.stats().join_probes[3] << '\n';
  }
  assert(events.stats().duplicate_sources == 0);

  const std::array<grammar::RuleId, 4> rules{
      grammar::RuleId::VP_WIP_VP_LEFT,
      grammar::RuleId::VP_VP_WIP_RIGHT,
      grammar::RuleId::VP_WIP_VPR,
      grammar::RuleId::VP_VPL_WIP};
  std::set<SourceId> exact_ids;
  std::set<SourceId> generator_ids;
  std::array<std::size_t, 4> counts{};
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = 1; i <= j; ++i) {
      const ItemKey parent{grammar::NonTerminal::VP, i, j, -1, -1};
      auto aggregate = Semiring::zero();
      for (const auto rule : rules) {
        for (const auto &source : deductions::generate_exact_basic_deductions(parent, rule, ctx)) {
          if (!engine::legacy_span_schedule_supports(source)) continue;
          const auto term = replay::runtime_deduction_contribution<Semiring>(
              chart, source, weight.weight_for(source));
          if (term != Semiring::zero()) generator_ids.emplace(rule, i, j, source.split.k);
        }
        // Independent reference: patched legacy compute_VP cases 6--9 sum
        // the complete border range. Do not ask the exact generator for this
        // range: a shared shortcut there and in events can hide lost terms.
        const bool left = rule == grammar::RuleId::VP_WIP_VP_LEFT ||
                          rule == grammar::RuleId::VP_WIP_VPR;
        const cand_pos_t begin = left ? i + 1
            : std::max(ctx.border_B(i, j), ctx.border_bp(i, j)) + 1;
        const cand_pos_t end = left
            ? static_cast<cand_pos_t>(std::min(
                  static_cast<cand_pos_tu>(ctx.border_b(i, j)),
                  static_cast<cand_pos_tu>(ctx.border_Bp(i, j))))
            : j;
        for (cand_pos_t k = begin; k < end; ++k) {
          grammar::RuleSplit split{};
          split.k = k;
          const auto source = deductions::expand_basic_rule(parent, rule, split);
          if (!deductions::is_valid(source)) continue;
          if (!engine::legacy_span_schedule_supports(source)) continue;
          const double local = weight.weight_for(source);
          if (local == 0.0) continue;
          const auto term = replay::runtime_deduction_contribution<Semiring>(
              chart, source, local);
          if (term == Semiring::zero()) continue;
          aggregate = Semiring::plus(aggregate, term);
          exact_ids.emplace(rule, i, j, source.split.k);
          ++counts[beam::vp_recursive_split_rule_index(rule)];
        }
      }
      const auto found = driven.find(parent);
      const auto actual = found == driven.end()
          ? Semiring::zero()
          : found->second;
      if (!close(aggregate, actual)) {
        std::cerr << "recursive VP event mismatch parent=" << i << ',' << j
                  << " exact=" << aggregate << " driven=" << actual << '\n';
        assert(false);
      }
    }
  }
  assert(driven_ids == exact_ids);
  assert(generator_ids == exact_ids);
  assert(events.stats().emitted_sources == counts);
  return counts;
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
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  auto inside =
      check_events<cparty::scfg2::semirings::SumProductSemiring>(owner);
  const auto viterbi =
      check_events<cparty::scfg2::semirings::MaxProductSemiring>(owner);
  assert(inside == viterbi);
  if (argc == 1) {
    // An active parent with exactly one contained maximal G pair; k=10
    // contributes to WIP*VPR although that pair opens at 11. Reflect the
    // input as well to exercise the right-oriented VPL*WIP loss.
    for (const bool mirror : {false, true}) {
      std::string seq = "GGGAAAAACGGAAAAACAAAACCAC";
      std::string g = "(.........(.....)....)...";
      if (mirror) {
        std::reverse(seq.begin(), seq.end());
        std::reverse(g.begin(), g.end());
        for (auto &c : g) {
          if (c == '(') c = ')';
          else if (c == ')') c = '(';
        }
      }
      cparty::scfg2::weights::RuntimeSnapshot witness(
      seq, g);
      const auto sum = check_events<cparty::scfg2::semirings::SumProductSemiring>(witness);
      const auto best = check_events<cparty::scfg2::semirings::MaxProductSemiring>(witness);
      assert(sum == best);
      for (std::size_t r = 0; r < inside.size(); ++r) inside[r] += sum[r];
    }
    assert(std::accumulate(inside.begin(), inside.end(), std::size_t{0}) > 0);
  }
  std::cout << "vp_wip_vp_left\tvp_vp_wip_right\t"
               "vp_wip_vpr\tvp_vpl_wip\n"
            << inside[0] << '\t' << inside[1] << '\t'
            << inside[2] << '\t' << inside[3] << '\n';
  return EXIT_SUCCESS;
}

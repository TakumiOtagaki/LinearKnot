#include "scfg2/weights/runtime_snapshot.hh"
#include "tools/scfg2_pilot_guided_inside.hh"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {
void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
bool close(double a, double b) {
  return a == b || std::fabs(a - b) < 1e-9;
}
void same_chart(const Chart &a, const Chart &b) {
  for (const auto &[key, value] : a.values())
    require(close(value, b.get(key)), "rank-only run changed a chart value");
  for (const auto &[key, value] : b.values())
    require(close(value, a.get(key)), "rank-only run changed chart support");
}
void below_exact(const Chart &beam, const Chart &exact) {
  for (const auto &[key, value] : beam.values())
    require(value <= exact.get(key) + 1e-8, "finite inside exceeds exact inside");
}

template <class Context>
void same_outer_reference(const PilotGuidedInsideRun &run, Context &ctx, const Item &root) {
  const auto reference = bounded_outer<Context, false>(
      run.pilot, root, ctx, run.pilot_width, 4 * run.pilot_width);
  const auto &a = run.guide.outer;
  const auto &b = reference.outer;
  require(a.values.values().size() == b.values.values().size(), "outer optimization changed support");
  for (const auto &[key, value] : a.values.values())
    require(value == b.values.get(key), "outer optimization changed a score bitwise");
  require(a.links.size() == b.links.size(), "outer optimization changed predecessor support");
  for (const auto &[key, link] : a.links) {
    const auto &other = b.links.at(key);
    require(sc::deductions::ItemKeyEq{}(link.parent, other.parent) && link.edge == other.edge,
            "outer optimization changed the chosen predecessor");
  }
  require(a.generated_rules == b.generated_rules && a.eligible_rules == b.eligible_rules,
          "outer optimization changed generated/eligible rules");
  // Production needs scores only. Finalizing a child must validate its link
  // against an already-final parent without retaining a second global map.
  const auto scores = sc::beam::bounded_outer<Context, true, false>(
      run.pilot, root, ctx, run.pilot_width, 4 * run.pilot_width);
  const auto &c = scores.outer;
  require(c.links.empty(), "score-only outer retained predecessor links");
  require(c.values.values().size() == a.values.values().size(), "score-only outer changed support");
  for (const auto &[key, value] : a.values.values())
    require(value == c.values.get(key), "score-only outer changed a score bitwise");
  require(a.visited == c.visited && a.generated == c.generated &&
          a.eligible == c.eligible && a.selected == c.selected &&
          a.max_parent_candidates == c.max_parent_candidates &&
          a.peak_pending == c.peak_pending && a.peak_guide_states == c.peak_guide_states &&
          a.evicted == c.evicted && a.generated_rules == c.generated_rules &&
          a.eligible_rules == c.eligible_rules,
          "score-only outer changed traversal/selection work");
}

void rank_contract() {
  for (const auto &[b, p] : std::vector<std::pair<std::size_t, std::size_t>>{
           {1, 1}, {2, 2}, {4, 2}, {5, 3}, {20, 5}, {100, 10}, {200, 15}})
    require(sqrt_pilot_width(b) == p, "sqrt pilot width is incorrect");
  for (const auto &[b, p] : std::vector<std::pair<std::size_t, std::size_t>>{
           {20, 8}, {27, 9}, {64, 16}, {100, 22}, {125, 25}})
    require(pilot_width_for(b, PilotExponent::TwoThirds) == p,
            "two-thirds pilot rounding is incorrect");
  // Independent small-integer inequalities, without floating point or the
  // fixed-word multiplication used by pilot_width_for.
  for (std::size_t b = 1; b <= 512; ++b) {
    const auto int_power = [](std::uint64_t x, unsigned e) {
      std::uint64_t product = 1;
      while (e--) product *= x;
      return product;
    };
    const auto p = pilot_width_for(b, PilotExponent::TwoThirds);
    require(int_power(p - 1, 3) < int_power(b, 2) &&
            int_power(p, 3) >= int_power(b, 2), "rational pilot width is not the ceiling");
  }
  Chart guide(zero), pilot(zero);
  PilotRankingCounts counts;
  const PilotOuterRanking rank{guide, pilot, counts};
  const Item item{NT::V, 7, 11, -1, -1};
  require(rank(item, -2, zero) == zero, "missing guide did not preserve prefix fallback");
  require(rank(item, -2, -9) == -9, "missing guide changed an existing finite rank");
  guide.set(item, -3);
  require(rank(item, -2, zero) == -5, "pilot-absent target cannot use a context");
  pilot.set(item, 100);
  require(rank(item, -2, zero) == -5, "target pilot alpha contaminated the rank");
  require(counts.total(counts.fallback) == 2 && counts.total(counts.novel) == 1,
          "fallback/discovery observation failed");
}

void witness(std::string sequence, std::string scaffold, sc::WmbClosurePolicy policy) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  sc::replay::WFinalLogWeightContext ctx(owner, policy);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const Item root{NT::W, 1, n, -1, -1};
  const auto no_observation = [](const auto &, const auto &, auto, const auto &, auto) {};
  const auto exact = sc::replay::run_w_final_exact_dp_schedule<Sum>(
      sc::engine::build_pkfree_basic_exact_schedule(n, ctx), ctx,
      [](const auto &) { return true; }, no_observation);
  require(std::isfinite(exact.get(root)), "exact witness root is unreachable");

  auto keep_all = sc::replay::WFinalVpWmbpOnlineWidths{};
  std::size_t calls = 0;
  const auto arbitrary = [&](const Item &key, double alpha, double) {
    ++calls;
    return alpha + 1000 * key.i - 300 * key.j;
  };
  const auto saturated = sc::replay::run_w_final_vp_wmbp_online_ranked_inside(
      ctx, n, keep_all, arbitrary);
  same_chart(exact, saturated.chart);
  require(calls > 0, "ranking experiment never reached a public selector");

  auto widths = sc::replay::WFinalVpWmbpOnlineWidths::uniform(20);
  widths.forced_items = sc::replay::scaffold_only_derivation_items(ctx, n, widths);
  for (std::size_t b : {2, 20}) {
    auto finite = sc::replay::WFinalVpWmbpOnlineWidths::uniform(b);
    finite.forced_items = widths.forced_items;
    const auto baseline = sc::replay::run_w_final_vp_wmbp_online<Sum>(ctx, n, finite);
    Chart empty(zero);
    PilotRankingCounts fallback_counts;
    const auto fallback = sc::replay::run_w_final_vp_wmbp_online_ranked_inside(
        ctx, n, finite, PilotOuterRanking{empty, empty, fallback_counts});
    same_chart(baseline.chart, fallback.chart);
    require(fallback_counts.total(fallback_counts.fallback) > 0,
            "empty guide fallback was not exercised");
    const auto guided = run_pilot_guided_inside(ctx, n, finite, b);
    const auto repeated = run_pilot_guided_inside(ctx, n, finite, b, 4, PilotExponent::TwoThirds);
    require(guided.pilot_width == pilot_width_for(b, PilotExponent::TwoThirds),
            "default pilot exponent is not two thirds");
    same_chart(guided.inside.chart, repeated.inside.chart);
    same_outer_reference(guided, ctx, root);
    require(std::isfinite(guided.inside.chart.get(root)), "guided root is unreachable");
    require(guided.ranks.total(guided.ranks.covered) > 0, "bounded guide is never used");
    if (scaffold == "((...)(...))") {
      require(guided.pilot.get({NT::W, 1, 6, -1, -1}) == zero,
              "forced multiloop no longer exercises the zero exterior prefix");
      require(std::isfinite(guided.guide.outer.values.get({NT::WMv, 7, 11, -1, -1})) &&
              std::isfinite(guided.guide.outer.values.get({NT::V, 7, 11, -1, -1})),
              "zero exterior prefix hides a valid multiloop context");
    }
    below_exact(baseline.chart, exact);
    below_exact(guided.pilot, exact);
    below_exact(guided.inside.chart, exact);
    if (b == 20) {
      const auto alternative = run_pilot_guided_inside(ctx, n, finite, b, 4, PilotExponent::Half);
      same_outer_reference(alternative, ctx, root);
      require(std::isfinite(alternative.inside.chart.get(root)), "q sweep disconnected the root");
      below_exact(alternative.pilot, exact);
      below_exact(alternative.inside.chart, exact);
    }
    for (cand_pos_t i = 1; i <= n; ++i) {
      const auto j = ctx.pair_partner(i);
      if (j > i)
        require(std::isfinite(guided.inside.chart.get({NT::V, i, j, -1, -1})),
                "guided ranking lost a forced scaffold pair");
    }
    std::cout << "n=" << n << " b=" << b << " p=" << guided.pilot_width
              << " covered=" << guided.ranks.total(guided.ranks.covered)
              << " fallback=" << guided.ranks.total(guided.ranks.fallback)
              << " novel=" << guided.ranks.total(guided.ranks.novel) << '\n';
  }
}

void changed_selection() {
  std::string sequence = "GCGCGCGCGCGCGCGCGCGCGCGC";
  std::string scaffold(sequence.size(), '.');
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  sc::replay::WFinalLogWeightContext ctx(owner);
  auto widths = sc::replay::WFinalVpWmbpOnlineWidths::uniform(2);
  widths.forced_items = sc::replay::scaffold_only_derivation_items(ctx, sequence.size(), widths);
  const auto baseline = sc::replay::run_w_final_vp_wmbp_online<Sum>(ctx, sequence.size(), widths);
  const auto changed = sc::replay::run_w_final_vp_wmbp_online_ranked_inside(
      ctx, sequence.size(), widths,
      [](const Item &key, double alpha, double existing) {
        return key.nonterminal == NT::V ? alpha + 1000 * key.i : existing;
      });
  std::size_t different = 0;
  for (const auto &[key, value] : baseline.chart.values())
    if (!close(value, changed.chart.get(key))) ++different;
  for (const auto &[key, value] : changed.chart.values())
    if (!close(value, baseline.chart.get(key))) ++different;
  require(different > 0, "public rank hook does not affect finite selection");
}
} // namespace

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    rank_contract();
    check_max_pool();
    bounded_outer_detail::check_prefix_top();
    for (auto policy : {sc::WmbClosurePolicy::RequireWeakClosure, sc::WmbClosurePolicy::LegacyCompatible}) {
      witness("GGAAACGAAACC", "((...)(...))", policy);
      witness("GCGCGCGCGCGCGCGCGCGCGCGC", "........................", policy);
      witness("GCGCGCGCGCGCGC", "(.(......)...)", policy);
    }
    changed_selection();
    std::cout << "scfg2_pilot_guided_inside=ok\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}

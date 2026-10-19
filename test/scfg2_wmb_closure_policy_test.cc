#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/centroid.hh"
#include "scfg2/api/mea.hh"
#include "scfg2/api/profiles.hh"
#include "scfg2/api/sampling.hh"
#include "scfg2/profiles/loop_context.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace {
using namespace cparty::scfg2;
using NT = grammar::NonTerminal;
using S = semirings::LogSumExpSemiring;
using Key = deductions::ItemKey;
void require(bool ok, const char *message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
bool close(double a, double b, double tol = 1e-10) {
  return a == b || std::fabs(a - b) <= tol;
}

void partial_prefix() {
  std::string seq = "GCGCGCGCGCGCGC", g = "(.(......)...)";
  cparty::scfg2::weights::RuntimeSnapshot pf(
      seq, g);
  replay::WFinalLogWeightContext ctx(pf);
  const auto chart = replay::run_w_final_exact_dp_schedule<S>(
      engine::build_pkfree_basic_exact_schedule(seq.size(), ctx), ctx,
      [](const auto &) { return true; },
      [](const auto &, const auto &, auto, const auto &, auto) {});
  const Key open{NT::WMB, 2, 7, -1, -1};
  require(!ctx.weakly_closed(2, 7) && !ctx.weakly_closed(8, 13) &&
          ctx.weakly_closed(2, 13), "closure witness geometry changed");
  require(std::isfinite(chart.get({NT::WMBP, 2, 7, -1, -1})),
          "open WMBP prefix was discarded");
  require(S::is_zero(chart.get(open)) &&
          S::is_zero(chart.get({NT::WMB, 8, 13, -1, -1})) &&
          std::isfinite(chart.get({NT::WMB, 2, 13, -1, -1})),
          "completed WMB domain is incorrect");
  require(!deductions::exact_basic_item_parent_active(open, ctx) &&
          deductions::generate_exact_basic_deductions_for_item(open, ctx).empty(),
          "public WMB parent gate bypassed");
  for (auto rule : {grammar::RuleId::WMB_DIRECT_WMBP,
                    grammar::RuleId::WMB_EMPTY,
                    grammar::RuleId::WMB_SPLIT_BE_WMBP_WI})
    require(deductions::generate_exact_basic_deductions(open, rule, ctx).empty(),
            "family-level WMB gate bypassed");
  require(deductions::exact_basic_item_parent_active(
              {NT::WMBW, 2, 7, -1, -1}, ctx),
          "partial WMBW incorrectly requires closure");
  for (std::size_t width : {20, 1000}) {
    const auto run = replay::run_w_final_vp_wmbp_online<S>(
        ctx, seq.size(), replay::WFinalVpWmbpOnlineWidths::uniform(width));
    require(std::isfinite(run.chart.get({NT::WMBP, 2, 7, -1, -1})),
            "beam lost required open prefix");
    for (const auto &[key, value] : run.chart.values())
      if (key.nonterminal == NT::WMB && !S::is_zero(value))
        require(ctx.weakly_closed(key.i, key.j), "beam materialized open WMB");
    require(close(run.chart.get({NT::W, 1, 14, -1, -1}),
                  chart.get({NT::W, 1, 14, -1, -1})), "saturated beam Z differs");
  }
  const auto evaluated = api::evaluate_energy(seq, g, ".[....][....].");
  // Independent legacy target-chain term with closure gate and patch0009;
  // correcting the ordinary closing mismatch changes 26.24 to 25.14.
  require(evaluated.ok && evaluated.reachable && close(evaluated.energy, 25.14),
          "valid canonical chain energy changed");
}

// Independent event query on the generic executor: forbid every deduction
// owning base 7, paired or unpaired. A nonzero root proves omitted ownership.
double omitted_mass(std::string seq, std::string g, WmbClosurePolicy policy) {
  cparty::scfg2::weights::RuntimeSnapshot pf(
      seq, g);
  replay::WFinalLogWeightContext ctx(pf, policy);
  const auto schedule = engine::build_pkfree_basic_exact_schedule(seq.size(), ctx);
  const auto no_observation = [](const auto &, const auto &, auto, const auto &, auto) {};
  const auto chart = replay::run_w_final_exact_dp_schedule<S>(
      schedule, ctx, [](const auto &) { return true; }, no_observation);
  const auto missing = replay::run_w_final_exact_dp_schedule<S>(
      schedule, ctx, [&](const auto &d) {
        const auto &key = d.parent;
        if ((key.nonterminal == NT::V || key.nonterminal == NT::VP) &&
            (key.i == 7 || key.j == 7)) return false;
        const auto emission = profiles::unpaired_emission(
            d, [&](cand_pos_t i) { return ctx.is_unpaired(i); });
        for (std::size_t r = 0; r < emission.count; ++r)
          if (emission.runs[r].lo <= 7 && 7 <= emission.runs[r].hi) return false;
        return true;
      }, no_observation);
  const Key root{NT::W, 1, static_cast<cand_pos_t>(seq.size()), -1, -1};
  if (policy == WmbClosurePolicy::RequireWeakClosure)
    require(S::is_zero(missing.get(root)), "omitted event remains reachable in log domain");
  return std::exp(missing.get(root) - chart.get(root));
}

void consumer_support() {
  const std::string seq = "GCAUAUUAAUAUUCUGCUAUAAAUAUAGGUGCACUUAGGGCAUGC";
  const std::string g = "((.(((....))(.(.....))).((....((.(....)))))))";
  for (auto policy : {WmbClosurePolicy::RequireWeakClosure,
                      WmbClosurePolicy::LegacyCompatible}) {
    api::ExactSession session(seq, g, policy);
    const auto z = session.inside();
    const auto outside = session.outside();
    const auto best = session.viterbi();
    const auto evaluated = session.evaluate_energy(best.generated_pairs);
    const auto probability = session.probability(best.generated_pairs);
    const auto profile = api::exact_structural_profiles(seq, g, {1, 2, 5}, policy);
    const auto sampled = api::exact_sample(seq, g, 3, 123, policy);
    const auto mea = api::exact_mea(seq, g, 1., policy);
    const auto centroid = api::exact_centroid(seq, g, 1., policy);
    require(z.ok && outside.ok && best.ok && evaluated.ok && probability.ok &&
            profile.ok && sampled.ok && mea.ok && centroid.ok, "exact consumer failed");
    require(z.wmb_closure_policy == policy && outside.wmb_closure_policy == policy &&
            best.wmb_closure_policy == policy && evaluated.wmb_closure_policy == policy &&
            probability.wmb_closure_policy == policy && profile.wmb_closure_policy == policy &&
            sampled.wmb_closure_policy == policy && mea.wmb_closure_policy == policy &&
            centroid.wmb_closure_policy == policy, "exact result lost policy provenance");
    for (double log_z : {outside.log_partition_weight, profile.log_partition_weight,
                         sampled.log_partition_weight, mea.log_partition_weight,
                         centroid.log_partition_weight})
      require(close(log_z, z.log_partition_weight), "exact consumers disagree on Z");
    require(best.log_weight <= z.log_partition_weight + 1e-10 &&
            close(best.energy, evaluated.energy) && probability.probability <= 1. + 1e-12,
            "exact cross-semiring consistency failed");
    const double omitted = omitted_mass(seq, g, policy);
    const double gap = outside.unpaired_probabilities[7] - profile.accessibility[0][7];
    if (policy == WmbClosurePolicy::RequireWeakClosure) {
      require(omitted == 0., "base 7 can still be omitted");
      for (std::size_t i = 1; i <= seq.size(); ++i) {
        require(close(profile.accessibility[0][i], outside.unpaired_probabilities[i]),
                "corrected accessibility disagrees with paired complement");
        require(close(std::accumulate(profile.contexts[i].begin(), profile.contexts[i].end(), 0.), 1.),
                "corrected contexts do not normalize");
      }
    } else {
      require(close(omitted, 2.32677792562853e-7, 1e-15) && close(gap, omitted, 1e-13),
              "legacy omitted-base regression changed");
    }
    std::cout << std::setprecision(17) << "policy=" << wmb_closure_policy_name(policy)
              << " logZ=" << z.log_partition_weight << " omitted=" << omitted << '\n';
    for (std::size_t width : {20, 1000}) {
      api::BeamConfig config;
      config.beam_size = width;
      config.wmb_closure_policy = policy;
      api::BeamSession beam_session(seq, g, config);
      const auto bz = beam_session.inside();
      const auto bo = beam_session.outside();
      const auto bv = beam_session.viterbi();
      const auto bp = beam_session.probability(bv.generated_pairs);
      const auto bc = api::beam_structural_profiles(seq, g, {1, 2, 5}, config);
      const auto bs = api::beam_sample(seq, g, 3, 123, config);
      const auto bm = api::beam_mea(seq, g, 1., config);
      const auto bg = api::beam_centroid(seq, g, 1., config);
      require(bz.ok && bo.ok && bv.ok && bp.ok && bc.ok && bs.ok && bm.ok && bg.ok,
              "beam consumer failed");
      require(bz.wmb_closure_policy == policy && bo.wmb_closure_policy == policy &&
              bv.wmb_closure_policy == policy && bp.wmb_closure_policy == policy &&
              bc.wmb_closure_policy == policy && bs.wmb_closure_policy == policy &&
              bm.wmb_closure_policy == policy && bg.wmb_closure_policy == policy,
              "beam result lost policy provenance");
      for (double log_z : {bo.log_partition_weight, bc.log_partition_weight,
                           bs.log_partition_weight, bm.log_partition_weight, bg.log_partition_weight})
        require(close(log_z, bz.log_partition_weight), "beam consumers disagree on retained Z");
      if (width == 1000) require(close(bz.log_partition_weight, z.log_partition_weight),
                                 "saturated beam differs from exact");
      if (policy == WmbClosurePolicy::RequireWeakClosure)
        for (std::size_t i = 1; i <= seq.size(); ++i)
          require(close(bc.accessibility[0][i], bo.unpaired_probabilities[i]),
                  "beam accessibility disagrees with paired complement");
    }
  }
}
}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  partial_prefix();
  consumer_support();
  std::cout << "scfg2_wmb_closure_policy=ok\n";
}

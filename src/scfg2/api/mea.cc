#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/mea.hh"
#include "scfg2/api/centroid.hh"

#include "scfg2/api/detail.hh"
#include "scfg2/api/beam_detail.hh"
#include "scfg2/beam/retained_outside.hh"
#include "scfg2/decoding/inside_pair_max.hh"
#include "scfg2/decoding/pair_posterior.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <chrono>
#include <exception>

namespace cparty::scfg2::api {
namespace {

template <class Run>
PosteriorDecodingResult decode_with(const std::string &sequence, const std::string &scaffold,
                   double gamma, decoding::PairPosteriorObjective objective,
                   std::size_t beam_size, WmbClosurePolicy policy,
                   Run &&run_inside) {
  const auto started = std::chrono::steady_clock::now();
  PosteriorDecodingResult out;
  out.wmb_closure_policy = policy;
  out.gamma = gamma;
  out.beam_size = beam_size;
  out.scaffold_structure = scaffold;
  if (!std::isfinite(gamma) || gamma < 0.0) {
    out.error = "decoding gamma must be finite and nonnegative";
  } else if (detail::validate_session_inputs(sequence, scaffold, &out.error)) {
    try {
      std::string seq = sequence, g = scaffold;
      weights::RuntimeSnapshot owner(seq, g);
      replay::WFinalLogWeightContext ctx(owner, policy);
      const auto n = static_cast<cand_pos_t>(sequence.size());
      using S = semirings::FactorizedLogSumExpSemiring;
      semirings::FactorizedInsideGraph graph;
      const auto chart = [&] {
        S::Scope scope(graph);
        return run_inside(ctx, n);
      }();
      const auto root = chart.get({grammar::NonTerminal::W, 1, n, -1, -1});
      out.log_partition_weight = root.log_value;
      if (S::is_zero(root)) {
        out.ok = true;  // Valid call, unreachable retained ensemble.
      } else {
        const auto marginal = engine::collect_outside_marginals(
            chart, graph.posteriors(root.node), n, ctx);
        if (!marginal.ok) throw std::runtime_error(marginal.error);
        const auto decoding_started = std::chrono::steady_clock::now();
        const auto decoded = decoding::decode_pair_posterior(
            graph, root.node, chart, marginal, n, ctx, gamma, objective);
        out.expected_accuracy = decoded.expected_accuracy;
        out.expected_true_positive_pairs = decoded.expected_true_positive_pairs;
        out.expected_true_negative_pairs = decoded.expected_true_negative_pairs;
        out.generated_structure = decoded.generated_structure;
        out.union_structure = detail::compose_union_structure(scaffold, decoded.generated_structure);
        // The decoder already ordered these pairs in linear time.
        for (const auto &pair : decoded.pairs) out.generated_pairs.push_back({pair.i, pair.j});
        out.decoding_states = graph.states().size();
        out.decoding_terms = decoded.terms_evaluated;
        out.decoding_factors = decoded.factors_evaluated;
        out.decoding_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - decoding_started).count();
        out.ok = out.reachable = true;
      }
    } catch (const std::exception &e) {
      out.ok = out.reachable = false;
      out.error = e.what();
    }
  }
  out.elapsed_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  return out;
}

PosteriorDecodingResult beam_decode(const std::string &sequence, const std::string &scaffold,
                   double gamma, BeamConfig config,
                   decoding::PairPosteriorObjective objective) {
  std::string config_error;
  if (!detail::validate_beam_config(config, &config_error)) {
    PosteriorDecodingResult out;
    out.wmb_closure_policy = config.wmb_closure_policy;
    out.gamma = gamma;
    out.error = config_error;
    return out;
  }
  const auto started = std::chrono::steady_clock::now();
  PosteriorDecodingResult out;
  out.wmb_closure_policy = config.wmb_closure_policy;
  out.gamma = gamma;
  out.beam_size = config.beam_size;
  out.scaffold_structure = scaffold;
  if (!std::isfinite(gamma) || gamma < 0.0) {
    out.error = "decoding gamma must be finite and nonnegative";
  } else if (detail::validate_session_inputs(sequence, scaffold, &out.error)) {
    try {
      std::string seq = sequence, g = scaffold;
      weights::RuntimeSnapshot owner(seq, g);
      replay::WFinalLogWeightContext ctx(owner, config.wmb_closure_policy);
      const auto n = static_cast<cand_pos_t>(sequence.size());
      const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
      const auto prepared = detail::prepare_beam_ranking(ctx, n, config);
      // Numeric inside/outside has no saved edge forest. Release those work
      // tables before replaying the identical alpha-ranked schedule below.
      const auto marginal = [&] {
        using L = semirings::LogSumExpSemiring;
        beam::RetainedRuleAdmission<L> rules;
        const auto numeric = replay::run_w_final_vp_wmbp_online_observed<L>(
            ctx, n, prepared.widths, rules, prepared.ranking);
        out.log_partition_weight = numeric.chart.get(root);
        if (L::is_zero(out.log_partition_weight)) return engine::OutsideMarginalSummary{};
        rules.finish(numeric.chart);
        beam::RetainedOutsideProgram<L, replay::WFinalLogWeightContext> outside(
            n, numeric.chart, ctx, rules, prepared.widths, prepared.ranking);
        const auto reverse = outside.outside(root);
        return engine::collect_outside_marginals_by_item(numeric.chart,
            [&](const auto &item, double alpha) {
              return std::exp(alpha + reverse.beta.get(item) - out.log_partition_weight);
            }, n, ctx);
      }();
      if (semirings::LogSumExpSemiring::is_zero(out.log_partition_weight)) {
        out.ok = true;
      } else {
        if (!marginal.ok) throw std::runtime_error(marginal.error);
        const auto decoding_started = std::chrono::steady_clock::now();
        const decoding::PairPosteriorScoring scoring(marginal, n, gamma, objective);
        using S = semirings::InsidePairMaxSemiring;
        auto storage = std::make_shared<semirings::PairMaxStorage>();
        const auto selected = [&] {
          S::Scope scope(storage);
          return replay::run_w_final_vp_wmbp_online_observed_transformed<S>(
              ctx, n, prepared.widths,
              [](const auto &, const auto &, auto, const auto &, auto) {},
              [&](const deductions::Deduction &d, auto value) {
                const auto feature = scoring.feature(d.parent, ctx, true);
                return feature ? S::with_pair_reward(value, d.parent, scoring.reward(*feature)) : value;
              }, prepared.ranking);
        }();
        if (selected.chart.get(root).log_value != out.log_partition_weight)
          throw std::runtime_error("posterior replay changed inside alpha");
        const auto decoded = decoding::decode_inside_pair_max(selected.chart.get(root),
            marginal, n, ctx, gamma, objective);
        out.expected_accuracy = decoded.expected_accuracy;
        out.expected_true_positive_pairs = decoded.expected_true_positive_pairs;
        out.expected_true_negative_pairs = decoded.expected_true_negative_pairs;
        out.generated_structure = decoded.generated_structure;
        out.union_structure = detail::compose_union_structure(scaffold, decoded.generated_structure);
        for (const auto &pair : decoded.pairs) out.generated_pairs.push_back({pair.i, pair.j});
        out.decoding_states = storage->peak_nodes;
        out.decoding_terms = decoded.terms_evaluated;
        out.decoding_factors = decoded.factors_evaluated;
        out.decoding_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - decoding_started).count();
        out.ok = out.reachable = true;
      }
    } catch (const std::exception &e) {
      out.ok = out.reachable = false;
      out.error = e.what();
    }
  }
  out.elapsed_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  return out;
}

PosteriorDecodingResult exact_decode(const std::string &sequence, const std::string &scaffold,
                                    double gamma, decoding::PairPosteriorObjective objective, WmbClosurePolicy policy) {
  return decode_with(sequence, scaffold, gamma, objective, 0, policy,
                 [](auto &ctx, cand_pos_t n) {
    return replay::run_w_final_exact_dp_schedule<semirings::FactorizedLogSumExpSemiring>(
        engine::build_pkfree_basic_exact_schedule(n, ctx), ctx,
        [](const deductions::Deduction &) { return true; },
        [](const auto &, const auto &, const auto &, const auto &, auto) {});
  });
}

}  // namespace

MeaResult beam_mea(const std::string &sequence, const std::string &scaffold,
                   double gamma, BeamConfig config) {
  return beam_decode(sequence, scaffold, gamma, config, decoding::PairPosteriorObjective::Mea);
}

MeaResult exact_mea(const std::string &sequence, const std::string &scaffold, double gamma, WmbClosurePolicy policy) {
  return exact_decode(sequence, scaffold, gamma, decoding::PairPosteriorObjective::Mea, policy);
}

CentroidResult beam_centroid(const std::string &sequence, const std::string &scaffold,
                             double gamma, BeamConfig config) {
  return beam_decode(sequence, scaffold, gamma, config, decoding::PairPosteriorObjective::GammaCentroid);
}

CentroidResult exact_centroid(const std::string &sequence, const std::string &scaffold, double gamma, WmbClosurePolicy policy) {
  return exact_decode(sequence, scaffold, gamma, decoding::PairPosteriorObjective::GammaCentroid, policy);
}

}  // namespace cparty::scfg2::api

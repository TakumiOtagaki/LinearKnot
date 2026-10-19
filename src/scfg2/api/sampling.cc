#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/sampling.hh"

#include "scfg2/api/detail.hh"
#include "scfg2/api/beam_detail.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/sampling/tape_sampler.hh"
#include "scfg2/sampling/factorized_sampler.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <chrono>
#include <cmath>
#include <exception>
#include <random>

namespace cparty::scfg2::api {
namespace {

template <class Semiring, class Graph, class Sampler, class RunFn>
SampleStreamResult sample_with(const std::string &sequence,
                         const std::string &scaffold,
                         std::size_t count,
                         std::uint64_t seed,
                         WmbClosurePolicy policy,
                         const SampleConsumer &consumer, RunFn &&run_inside) {
  const auto started = std::chrono::steady_clock::now();
  SampleStreamResult out;
  if (!consumer) out.error = "sample consumer must be callable";
  else if (detail::validate_session_inputs(sequence, scaffold, &out.error)) {
    try {
      std::string seq = sequence;
      std::string g = scaffold;
      weights::RuntimeSnapshot pf(seq, g);
      replay::WFinalLogWeightContext ctx(pf, policy);
      const auto n = static_cast<cand_pos_t>(sequence.size());
      Graph graph;
      const auto chart = [&] {
        typename Semiring::Scope scope(graph);
        return run_inside(ctx, n);
      }();
      const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
      const auto root_value = chart.get(root);
      out.log_partition_weight = root_value.log_value;
      if (!std::isfinite(root_value.log_value)) {
        out.error = "no structure is reachable";
      } else {
        out.ensemble_free_energy =
            pf.partition_free_energy_from_log_weight(root_value.log_value);
        Sampler sampler(graph, root_value.node, chart, ctx);
        std::mt19937_64 rng(seed);
        for (std::size_t s = 0; s < count; ++s) {
          const auto pairs = sampler.sample(rng);
          SampledStructure structure;
          structure.generated_pairs = detail::to_base_pairs(pairs);
          structure.union_structure = detail::compose_union_structure(scaffold,
              detail::render_generated_structure(static_cast<int>(n), pairs, ctx));
          ++out.samples_emitted;
          if (!consumer(std::move(structure))) {
            out.stopped_by_consumer = true;
            break;
          }
        }
        out.ok = true;
      }
    } catch (const std::exception &e) {
      out.error = e.what();
    } catch (...) {
      out.error = "sampling consumer or kernel threw a non-standard exception";
    }
  }
  out.wmb_closure_policy = policy;
  out.elapsed_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  return out;
}

template <class RunEach>
SampleResult collect_samples(RunEach &&run_each) {
  SampleResult out;
  const auto summary = run_each([&](SampledStructure &&structure) {
    out.generated_pairs.push_back(std::move(structure.generated_pairs));
    out.union_structures.push_back(std::move(structure.union_structure));
    return true;
  });
  if (!summary.ok) out = SampleResult{};  // Historical all-or-error collection.
  out.ok = summary.ok;
  out.error = summary.error;
  out.wmb_closure_policy = summary.wmb_closure_policy;
  out.beam_size = summary.beam_size;
  out.elapsed_seconds = summary.elapsed_seconds;
  if (summary.ok) {
    out.log_partition_weight = summary.log_partition_weight;
    out.ensemble_free_energy = summary.ensemble_free_energy;
  }
  return out;
}

}  // namespace

SampleStreamResult exact_sample_each(const std::string &sequence,
                          const std::string &scaffold_structure,
                          std::size_t count,
                          std::uint64_t seed, const SampleConsumer &consumer, WmbClosurePolicy policy) {
  return sample_with<semirings::TapeLogSumExpSemiring,
                     semirings::LogSumExpTape, sampling::TapeStructureSampler>(
                     sequence, scaffold_structure, count, seed,
                     policy, consumer,
                     [](auto &ctx, cand_pos_t n) {
                       const auto schedule =
                           engine::build_pkfree_basic_exact_schedule(n, ctx);
                       return replay::run_w_final_exact_dp_schedule_with_contribution<
                           semirings::TapeLogSumExpSemiring>(
                           schedule, ctx,
                           [](const deductions::Deduction &) { return true; },
                           [](const deductions::ItemKey &,
                              const deductions::Deduction &, auto,
                              const auto &, auto) {},
                           [](const deductions::Deduction &deduction,
                              const auto &chart, auto local) {
                             using S = semirings::TapeLogSumExpSemiring;
                             // A zero child makes the entire candidate dead.
                             // Check before appending its local/product tape
                             // nodes; keep every live operation in its original
                             // order, including node-based pair ownership.
                             for (std::size_t c = 0; c < deduction.child_count; ++c) {
                               if (S::is_zero(replay::runtime_child_total_for_deduction(
                                       chart, deduction, c))) return S::zero();
                             }
                             return replay::runtime_deduction_contribution<S>(
                                 chart, deduction, local);
                           });
                     });
}

SampleStreamResult beam_sample_each(const std::string &sequence,
                         const std::string &scaffold_structure,
                         std::size_t count,
                         std::uint64_t seed,
                         const SampleConsumer &consumer, BeamConfig config) {
  std::string config_error;
  if (!detail::validate_beam_config(config, &config_error)) {
    SampleStreamResult out;
    out.wmb_closure_policy = config.wmb_closure_policy;
    out.error = config_error;
    return out;
  }
  auto out = sample_with<semirings::FactorizedLogSumExpSemiring,
                         semirings::FactorizedInsideGraph,
                         sampling::FactorizedStructureSampler>(
      sequence, scaffold_structure, count, seed,
      config.wmb_closure_policy, consumer,
      [&](auto &ctx, cand_pos_t n) {
        return detail::run_beam_inside<semirings::FactorizedLogSumExpSemiring>(
            ctx, n, config).chart;
      });
  out.beam_size = config.beam_size;
  return out;
}

SampleResult exact_sample(const std::string &sequence, const std::string &scaffold,
                          std::size_t count, std::uint64_t seed, WmbClosurePolicy policy) {
  return collect_samples([&](const SampleConsumer &consumer) {
    return exact_sample_each(sequence, scaffold, count, seed, consumer, policy);
  });
}

SampleResult beam_sample(const std::string &sequence, const std::string &scaffold,
                         std::size_t count, std::uint64_t seed, BeamConfig config) {
  return collect_samples([&](const SampleConsumer &consumer) {
    return beam_sample_each(sequence, scaffold, count, seed, consumer, config);
  });
}

}  // namespace cparty::scfg2::api

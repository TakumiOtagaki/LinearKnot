#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/subopt.hh"

#include "scfg2/api/detail.hh"
#include "scfg2/api/beam_detail.hh"
#include "scfg2/engine/factorized_pair_owners.hh"
#include "scfg2/replay/lazy_kbest.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <chrono>
#include <cmath>
#include <unordered_set>

namespace cparty::scfg2::api {
namespace {
template <class Run>
SuboptStreamResult enumerate_with(const std::string &sequence, const std::string &scaffold,
                            SuboptConfig config, std::size_t beam_size,
                            WmbClosurePolicy policy, const SuboptConsumer &consumer, Run &&run_inside) {
  const auto started = std::chrono::steady_clock::now();
  SuboptStreamResult out;
  out.wmb_closure_policy = policy;
  out.beam_size = beam_size;
  if (!consumer) out.error = "subopt consumer must be callable";
  else if (config.top_k == 0) out.error = "top_k must be positive";
  else if (config.delta_energy &&
           (!std::isfinite(*config.delta_energy) || *config.delta_energy < 0.0))
    out.error = "delta energy must be finite and nonnegative";
  else if (detail::validate_session_inputs(sequence, scaffold, &out.error)) {
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
      out.graph_states = graph.states().size();
      out.graph_terms = graph.terms().size();
      out.graph_factors = graph.factors().size();
      if (!S::is_zero(root)) {
        replay::LazyKBest lazy(graph);
        engine::FactorizedPairOwners pair_owners(graph, chart, ctx);
        const auto best_log = lazy.best_log_weight(root.node);
        out.mfe_energy = owner.partition_free_energy_from_log_weight(best_log);
        const double rt = owner.partition_free_energy_from_log_weight(-1.0);
        std::unordered_set<std::string> seen_pair_maps;
        for (std::size_t rank = 0; rank < config.top_k; ++rank) {
          if (!lazy.ensure(root.node, rank)) break;
          const auto score = lazy.derivation(root.node, rank).log_weight;
          if (config.delta_energy && rt * (best_log - score) > *config.delta_energy + 1e-9) {
            out.stop_reason = SuboptStopReason::DeltaEnergy;
            break;
          }
          SuboptStructure structure;
          structure.log_weight = score;
          structure.energy = owner.partition_free_energy_from_log_weight(score);
          structure.generated_structure.assign(sequence.size(), '.');
          std::vector<cand_pos_t> partner(sequence.size() + 1, 0);
          lazy.visit_derivation(root.node, rank, [&](std::uint32_t node) {
            const auto *pair = pair_owners.find(node);
            if (!pair) return;
            if (partner[pair->i] == pair->j) return;  // public wrapper aliases
            if (partner[pair->i] || partner[pair->j])
              throw std::logic_error("subopt derivation reuses a paired endpoint");
            partner[pair->i] = pair->j;
            partner[pair->j] = pair->i;
            const bool round = pair->lane == engine::TracebackPairLane::Round;
            structure.generated_structure[pair->i - 1] = round ? '(' : '[';
            structure.generated_structure[pair->j - 1] = round ? ')' : ']';
          });
          // Linear pair ordering/rendering and a lane-independent identity.
          // Fixed-width binary endpoints avoid decimal-index log(n) output.
          std::string key;
          for (cand_pos_t i = 1; i <= n; ++i) {
            if (partner[i] <= i) continue;
            structure.generated_pairs.push_back({static_cast<int>(i), static_cast<int>(partner[i])});
            for (const auto endpoint : {i, partner[i]})
              for (unsigned shift = 0; shift < 32; shift += 8)
                key.push_back(static_cast<char>((static_cast<std::uint32_t>(endpoint) >> shift) & 255));
          }
          if (!seen_pair_maps.insert(std::move(key)).second)
            throw std::logic_error("subopt found multiple derivations for one physical pair map (grammar ambiguity)");
          structure.union_structure = detail::compose_union_structure(scaffold, structure.generated_structure);
          ++out.structures_emitted;
          if (!consumer(std::move(structure))) {
            out.stop_reason = SuboptStopReason::Consumer;
            break;
          }
          if (out.structures_emitted == config.top_k) out.stop_reason = SuboptStopReason::TopK;
        }
        out.candidates_created = lazy.candidates_created();
        out.derivations_cached = lazy.derivations_cached();
        out.reachable = true;
      }
      out.ok = true;
    } catch (const std::exception &e) {
      out.ok = out.reachable = false;
      out.error = e.what();
    } catch (...) {
      out.ok = out.reachable = false;
      out.error = "subopt consumer or kernel threw a non-standard exception";
    }
  }
  out.elapsed_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  return out;
}
template <class RunEach>
SuboptResult collect_subopt(RunEach &&run_each) {
  SuboptResult out;
  const auto summary = run_each([&](SuboptStructure &&structure) {
    out.structures.push_back(std::move(structure));
    return true;
  });
  if (!summary.ok) out.structures.clear();
  out.ok = summary.ok;
  out.reachable = summary.reachable;
  out.error = summary.error;
  out.wmb_closure_policy = summary.wmb_closure_policy;
  out.beam_size = summary.beam_size;
  out.stop_reason = summary.stop_reason;
  out.log_partition_weight = summary.log_partition_weight;
  out.mfe_energy = summary.mfe_energy;
  out.graph_states = summary.graph_states;
  out.graph_terms = summary.graph_terms;
  out.graph_factors = summary.graph_factors;
  out.candidates_created = summary.candidates_created;
  out.derivations_cached = summary.derivations_cached;
  out.elapsed_seconds = summary.elapsed_seconds;
  return out;
}
}  // namespace

SuboptStreamResult exact_subopt_each(const std::string &sequence, const std::string &scaffold,
                          const SuboptConsumer &consumer, SuboptConfig config, WmbClosurePolicy policy) {
  return enumerate_with(sequence, scaffold, config, 0, policy, consumer, [](auto &ctx, cand_pos_t n) {
    return replay::run_w_final_exact_dp_schedule<semirings::FactorizedLogSumExpSemiring>(
        engine::build_pkfree_basic_exact_schedule(n, ctx), ctx,
        [](const deductions::Deduction &) { return true; },
        [](const auto &, const auto &, const auto &, const auto &, auto) {});
  });
}

SuboptStreamResult beam_subopt_each(const std::string &sequence, const std::string &scaffold,
                         const SuboptConsumer &consumer, SuboptConfig config, BeamConfig beam) {
  std::string config_error;
  if (!detail::validate_beam_config(beam, &config_error)) {
    SuboptStreamResult out;
    out.wmb_closure_policy = beam.wmb_closure_policy;
    out.error = config_error;
    return out;
  }
  return enumerate_with(sequence, scaffold, config, beam.beam_size, beam.wmb_closure_policy, consumer,
                         [&](auto &ctx, cand_pos_t n) {
    return detail::run_beam_inside<semirings::FactorizedLogSumExpSemiring>(
        ctx, n, beam).chart;
  });
}
SuboptResult exact_subopt(const std::string &sequence, const std::string &scaffold,
                          SuboptConfig config, WmbClosurePolicy policy) {
  return collect_subopt([&](const SuboptConsumer &consumer) {
    return exact_subopt_each(sequence, scaffold, consumer, config, policy);
  });
}

SuboptResult beam_subopt(const std::string &sequence, const std::string &scaffold,
                         SuboptConfig config, BeamConfig beam) {
  return collect_subopt([&](const SuboptConsumer &consumer) {
    return beam_subopt_each(sequence, scaffold, consumer, config, beam);
  });
}
}  // namespace cparty::scfg2::api

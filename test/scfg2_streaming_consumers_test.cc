#include "scfg2/api/sampling.hh"
#include "scfg2/api/subopt.hh"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

using namespace cparty::scfg2;
namespace {
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
bool same_pairs(const std::vector<api::BasePair> &a, const std::vector<api::BasePair> &b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (a[i].i != b[i].i || a[i].j != b[i].j) return false;
  return true;
}

api::SuboptConfig top(std::size_t count) {
  api::SuboptConfig config;
  config.top_k = count;
  return config;
}

void samples(bool exact, api::BeamConfig config) {
  const std::string seq = "GCGCGCGCGCGCGC", g = "(.(......)...)";
  const auto batch = exact ? api::exact_sample(seq, g, 32, 73, config.wmb_closure_policy)
                           : api::beam_sample(seq, g, 32, 73, config);
  require(batch.ok, "sample batch failed");
  const auto run = [&](std::size_t count, const api::SampleConsumer &consumer) {
    return exact ? api::exact_sample_each(seq, g, count, 73, consumer, config.wmb_closure_policy)
                 : api::beam_sample_each(seq, g, count, 73, consumer, config);
  };
  std::size_t received = 0;
  const auto full = run(32, [&](api::SampledStructure &&s) {
    require(s.union_structure == batch.union_structures.at(received) &&
            same_pairs(s.generated_pairs, batch.generated_pairs.at(received)), "sample draw order changed");
    ++received;
    return true;
  });
  require(full.ok && full.samples_emitted == 32 && received == 32 && !full.stopped_by_consumer &&
          full.log_partition_weight == batch.log_partition_weight &&
          full.ensemble_free_energy == batch.ensemble_free_energy &&
          full.beam_size == batch.beam_size && full.wmb_closure_policy == config.wmb_closure_policy,
          "sample summary differs");
  received = 0;
  // No count-sized reserve or eager request of all outputs is possible here.
  const auto stopped = run(std::numeric_limits<std::size_t>::max(), [&](api::SampledStructure &&s) {
    require(s.union_structure == batch.union_structures.at(received), "sample prefix changed with count");
    return ++received < 3;
  });
  require(stopped.ok && stopped.samples_emitted == 3 && received == 3 && stopped.stopped_by_consumer,
          "sample consumer stop generated extra output");
  const auto zero = run(0, [&](auto &&) { throw std::logic_error("unexpected zero-count callback"); return true; });
  require(zero.ok && zero.samples_emitted == 0 && !zero.stopped_by_consumer, "zero samples failed");
  const auto empty = run(1, {});
  require(!empty.ok && empty.samples_emitted == 0, "empty sample consumer was accepted");
  received = 0;
  const auto error = run(32, [&](auto &&) {
    if (++received == 2) throw std::runtime_error("sample sink failure");
    return true;
  });
  require(!error.ok && error.error == "sample sink failure" && error.samples_emitted == 2 && received == 2,
          "sample error lost delivered prefix metadata");
  const auto foreign = run(1, [](auto &&) -> bool { throw 7; });
  require(!foreign.ok && foreign.samples_emitted == 1 && !foreign.error.empty(), "nonstandard sample error escaped");
}

void subopt(bool exact, api::BeamConfig beam) {
  const std::string seq = "GCGCGCGCGCGCGC", g = "(.(......)...)";
  auto config = top(10);
  const auto batch = exact ? api::exact_subopt(seq, g, config, beam.wmb_closure_policy)
                           : api::beam_subopt(seq, g, config, beam);
  require(batch.ok && batch.structures.size() >= 3, "subopt batch lacks witness ranks");
  const auto run = [&](api::SuboptConfig c, const api::SuboptConsumer &consumer) {
    return exact ? api::exact_subopt_each(seq, g, consumer, c, beam.wmb_closure_policy)
                 : api::beam_subopt_each(seq, g, consumer, c, beam);
  };
  std::size_t received = 0;
  const auto full = run(config, [&](api::SuboptStructure &&s) {
    const auto &want = batch.structures.at(received++);
    require(s.log_weight == want.log_weight && s.energy == want.energy &&
            s.union_structure == want.union_structure && s.generated_structure == want.generated_structure &&
            same_pairs(s.generated_pairs, want.generated_pairs), "subopt rank/tie order changed");
    return true;
  });
  require(full.ok && full.reachable && full.structures_emitted == batch.structures.size() &&
          full.stop_reason == batch.stop_reason && full.log_partition_weight == batch.log_partition_weight &&
          full.mfe_energy == batch.mfe_energy && full.graph_states == batch.graph_states &&
          full.graph_terms == batch.graph_terms && full.graph_factors == batch.graph_factors &&
          full.candidates_created == batch.candidates_created && full.derivations_cached == batch.derivations_cached,
          "subopt summary changed");
  received = 0;
  const auto stopped = run(top(std::numeric_limits<std::size_t>::max()), [&](api::SuboptStructure &&s) {
    require(s.union_structure == batch.structures.at(received).union_structure, "subopt prefix changed with K");
    return ++received < 3;
  });
  const auto top3 = exact ? api::exact_subopt(seq, g, top(3), beam.wmb_closure_policy)
                          : api::beam_subopt(seq, g, top(3), beam);
  require(stopped.ok && received == 3 && stopped.structures_emitted == 3 &&
          stopped.stop_reason == api::SuboptStopReason::Consumer &&
          stopped.candidates_created == top3.candidates_created &&
          stopped.derivations_cached == top3.derivations_cached, "consumer stop requested rank K+1");
  const auto empty = run(config, {});
  require(!empty.ok && empty.structures_emitted == 0, "empty subopt consumer was accepted");
  received = 0;
  const auto error = run(config, [&](auto &&) {
    if (++received == 2) throw std::runtime_error("subopt sink failure");
    return true;
  });
  require(!error.ok && error.error == "subopt sink failure" && error.structures_emitted == 2 && received == 2,
          "subopt error lost delivered prefix metadata");
  const auto foreign = run(config, [](auto &&) -> bool { throw 7; });
  require(!foreign.ok && foreign.structures_emitted == 1 && !foreign.error.empty(), "nonstandard subopt error escaped");
  config.delta_energy = 0.;
  const auto band = exact ? api::exact_subopt(seq, g, config, beam.wmb_closure_policy)
                          : api::beam_subopt(seq, g, config, beam);
  received = 0;
  const auto band_stream = run(config, [&](auto &&s) {
    require(s.log_weight == band.structures.at(received++).log_weight, "stream energy band changed"); return true;
  });
  require(band_stream.ok && received == band.structures.size() && band_stream.stop_reason == band.stop_reason,
          "stream band termination changed");
}
} // namespace

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", SCFG2_TEST_PARAM_FILE, 1);
    api::BeamConfig beam{3};
    for (auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible}) {
      beam.wmb_closure_policy = policy;
      samples(true, beam); subopt(true, beam);
      for (auto ranking : {api::BeamRanking::Baseline, api::BeamRanking::Outer}) {
        beam.ranking = ranking;
        samples(false, beam); subopt(false, beam);
      }
    }
    std::size_t calls = 0;
    const auto exhausted = api::beam_subopt_each("AAAAAA", "......", [&](auto &&s) {
      require(s.union_structure == "......", "exhaustion witness structure changed"); ++calls; return true;
    }, top(100), api::BeamConfig{1});
    require(exhausted.ok && exhausted.reachable && calls == 1 && exhausted.structures_emitted == 1 &&
            exhausted.stop_reason == api::SuboptStopReason::Exhausted, "stream exhaustion failed");
    const auto invalid = api::beam_sample_each("AAAAAA", "......", 100, 73,
        [&](auto &&) { ++calls; return true; }, api::BeamConfig{0});
    require(!invalid.ok && invalid.samples_emitted == 0 && calls == 1, "invalid beam delivered output");
    const auto invalid_subopt = api::exact_subopt_each("AAAAAA", "......",
        [&](auto &&) { ++calls; return true; }, top(0));
    require(!invalid_subopt.ok && invalid_subopt.structures_emitted == 0 && calls == 1, "invalid subopt delivered output");
    // The legacy domain has a known duplicate derivation. Stream a valid
    // prefix, then fail; the collecting API must still clear its outputs.
    const std::string seq = "GCGCGCGCGCGCGC", g = "(.(......)...)";
    const auto legacy = api::exact_subopt(seq, g, top(100000), WmbClosurePolicy::LegacyCompatible);
    calls = 0;
    const auto streamed_legacy = api::exact_subopt_each(seq, g, [&](auto &&) { ++calls; return true; },
        top(100000), WmbClosurePolicy::LegacyCompatible);
    require(!legacy.ok && legacy.structures.empty() && !streamed_legacy.ok &&
            streamed_legacy.error == legacy.error && streamed_legacy.structures_emitted == calls && calls > 0,
            "streaming changed ambiguity detection or batch rollback");
    std::cout << "scfg2_streaming_consumers=ok modes=6 sampling_and_subopt stop_prefix_errors=ok\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}

// LinearKnot sampling: exact arithmetic-tape and beam factorized-forest sampling.
// Fixed-seed exact draws also match an unfiltered arithmetic-tape reference.
//   1. Exact: each frequent sampled structure's frequency matches its exact
//      probability; sampled pair frequencies match exact outside bpp.
//   2. KeepAll beam: sampled pair frequencies match exact outside bpp.
//   3. Finite beam: sampled pair frequencies match beam_outside bpp, and
//      sampled window-unpaired frequencies match beam accessibility (an
//      independent check of finite-beam windows with w > 1).
// Frequency comparisons allow 5 binomial standard errors plus 2/N.

#include "scfg2/api/exact.hh"
#include "scfg2/api/profiles.hh"
#include "scfg2/api/sampling.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/sampling/tape_sampler.hh"
#include "scfg2/weights/runtime_snapshot.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace {

using namespace cparty::scfg2::api;
using PairKey = std::pair<int, int>;

bool fail(const std::string &label, const std::string &message) {
  std::cerr << "[" << label << "] " << message << "\n";
  return false;
}

// Preserve the old local-then-children recording path as a reference. The
// public exact sampler skips candidates with a zero child before recording;
// this helper intentionally retains their otherwise unreachable tape nodes.
template <class Context, class Observer>
auto unfiltered_exact_tape(Context &ctx, cand_pos_t n,
    cparty::scfg2::semirings::LogSumExpTape &tape, Observer observe) {
  using namespace cparty::scfg2;
  using T = semirings::TapeLogSumExpSemiring;
  T::Scope scope(tape);
  return replay::run_w_final_exact_dp_schedule<T>(
      engine::build_pkfree_basic_exact_schedule(n, ctx), ctx,
      [](const deductions::Deduction &) { return true; }, observe);
}

bool check_exact_tape_reference() {
  using namespace cparty::scfg2;
  struct Fixture {
    std::string label, sequence, scaffold;
    bool long_be = false;
    bool crossing = false;
  };
  std::string long_sequence(42, 'A');
  for (std::size_t i = 0; i < long_sequence.size(); ++i)
    long_sequence[i] = "GCUA"[i % 4];
  for (const int i : {1, 18}) long_sequence[i - 1] = 'G';
  for (const int i : {26, 42}) long_sequence[i - 1] = 'C';
  const std::vector<Fixture> fixtures{
      {"empty-G", "GAGAAACAC", "........."},
      {"VP-split", "GGGAAAAACGGAAAAACAAAACCAC", "(.........(.....)....)...", false, true},
      {"nested-G", "GCGCGCGCGCGCGC", "(.(......)...)"},
      {"long-BE", long_sequence,
       "(" + std::string(16, '.') + "(" + std::string(7, '.') + ")" +
           std::string(15, '.') + ")", true}};
  constexpr std::size_t draws = 128;
  for (const auto policy : {WmbClosurePolicy::RequireWeakClosure,
                            WmbClosurePolicy::LegacyCompatible}) {
    for (const auto &fixture : fixtures) {
      weights::RuntimeSnapshot owner(fixture.sequence, fixture.scaffold);
      replay::WFinalLogWeightContext ctx(owner, policy);
      const auto n = static_cast<cand_pos_t>(fixture.sequence.size());
      semirings::LogSumExpTape tape;
      bool saw_long_be = false;
      const auto chart = unfiltered_exact_tape(ctx, n, tape,
          [&](const auto &, const deductions::Deduction &d, auto value,
              const auto &, auto) {
            if (d.rule == grammar::RuleId::BE_INTERNAL_LOOP &&
                d.split.k - d.parent.i - 1 + d.parent.j - d.split.l - 1 > 30 &&
                std::isfinite(value.log_value)) saw_long_be = true;
          });
      const auto root = chart.get({grammar::NonTerminal::W, 1, n, -1, -1});
      if (!std::isfinite(root.log_value) || (fixture.long_be && !saw_long_be))
        return fail(fixture.label, "reference lost its reachable ensemble or long-BE contribution");
      sampling::TapeStructureSampler sampler(tape, root.node, chart, ctx);
      bool saw_crossing = false;
      for (const std::uint64_t seed : {7, 73, 20261018}) {
        const auto actual = exact_sample(fixture.sequence, fixture.scaffold, draws, seed, policy);
        if (!actual.ok || actual.log_partition_weight != root.log_value ||
            actual.generated_pairs.size() != draws || actual.union_structures.size() != draws)
          return fail(fixture.label, "exact sampling summary differs from unfiltered tape");
        std::mt19937_64 rng(seed);
        for (std::size_t s = 0; s < draws; ++s) {
          const auto expected = sampler.sample(rng);
          const auto &pairs = actual.generated_pairs[s];
          if (pairs.size() != expected.size())
            return fail(fixture.label, "fixed-seed sample pair count changed");
          std::string structure = fixture.scaffold;
          for (std::size_t p = 0; p < pairs.size(); ++p) {
            const auto &pair = expected[p];
            if (pairs[p].i != pair.i || pairs[p].j != pair.j)
              return fail(fixture.label, "fixed-seed ordered sample pairs changed");
            const bool square = pair.lane == engine::TracebackPairLane::Square;
            saw_crossing = saw_crossing || square;
            structure[pair.i - 1] = square ? '[' : '(';
            structure[pair.j - 1] = square ? ']' : ')';
          }
          if (actual.union_structures[s] != structure)
            return fail(fixture.label, "fixed-seed sample rendering changed");
        }
      }
      if (fixture.crossing && !saw_crossing)
        return fail(fixture.label, "fixed-seed control did not exercise a crossing pair");
    }
  }
  std::cout << "exact_unfiltered_tape_controls=ok fixtures=4 closures=2 seeds=3 draws=128\n";
  return true;
}

bool check_exact_empty_and_unreachable() {
  using namespace cparty::scfg2;
  const std::string sequence(32, 'A'), scaffold(32, '.');
  for (const auto policy : {WmbClosurePolicy::RequireWeakClosure,
                            WmbClosurePolicy::LegacyCompatible}) {
    // The unfiltered reference establishes that these valid inputs have,
    // respectively, one unit-weight structure and no reachable structure.
    for (const bool reachable : {true, false}) {
      const std::string seq = reachable ? sequence : "AAAAAA";
      const std::string g = reachable ? scaffold : "(....)";
      weights::RuntimeSnapshot owner(seq, g);
      replay::WFinalLogWeightContext ctx(owner, policy);
      const auto n = static_cast<cand_pos_t>(seq.size());
      semirings::LogSumExpTape tape;
      const auto chart = unfiltered_exact_tape(ctx, n, tape,
          [](const auto &, const auto &, auto, const auto &, auto) {});
      const auto reference = chart.get({grammar::NonTerminal::W, 1, n, -1, -1});
      if (reachable ? reference.log_value != 0.0
                    : !semirings::TapeLogSumExpSemiring::is_zero(reference))
        return fail("exact-boundary", "unexpected unfiltered reference reachability");
      for (const std::size_t count : {0, 1}) {
        std::size_t received = 0;
        bool empty_structure = true;
        const auto streamed = exact_sample_each(seq, g, count, 73,
            [&](SampledStructure &&sample) {
              ++received;
              empty_structure = empty_structure && sample.generated_pairs.empty() &&
                  sample.union_structure == g;
              return true;
            }, policy);
        if (reachable) {
          if (!streamed.ok || streamed.log_partition_weight != 0.0 ||
              streamed.ensemble_free_energy != 0.0 || streamed.samples_emitted != count ||
              received != count || !empty_structure || streamed.stopped_by_consumer)
            return fail("exact-boundary", "all-unpaired sampling or zero-count summary changed");
        } else if (streamed.ok || streamed.error != "no structure is reachable" ||
                   streamed.samples_emitted != 0 || received != 0 ||
                   streamed.log_partition_weight != reference.log_value) {
          return fail("exact-boundary", "unreachable sampling delivered a structure or changed its error");
        }
      }
    }
  }
  return true;
}

bool within(double frequency, double p, std::size_t n) {
  const double sigma = std::sqrt(std::max(p * (1.0 - p), 0.0) / n);
  return std::fabs(frequency - p) <= 5.0 * sigma + 2.0 / n;
}

std::map<PairKey, double> sample_pair_frequency(const SampleResult &samples) {
  std::map<PairKey, double> out;
  for (const auto &pairs : samples.generated_pairs) {
    for (const auto &pair : pairs) out[{pair.i, pair.j}] += 1.0;
  }
  for (auto &[_, value] : out) value /= samples.generated_pairs.size();
  return out;
}

bool check_pairs(const std::string &label, const SampleResult &samples,
                 const std::vector<PairProbability> &want,
                 const std::string &scaffold) {
  auto got = sample_pair_frequency(samples);
  const std::size_t n = samples.generated_pairs.size();
  for (const auto &pair : want) {
    const bool scaffold_pair =
        pair.i >= 1 && scaffold[pair.i - 1] == '(' && scaffold[pair.j - 1] == ')';
    if (scaffold_pair) continue;  // G pairs are not generated
    const double f = got[{pair.i, pair.j}];
    if (!within(f, pair.union_probability, n)) {
      return fail(label, "pair " + std::to_string(pair.i) + "-" +
                             std::to_string(pair.j) + " sampled " +
                             std::to_string(f) + " outside " +
                             std::to_string(pair.union_probability));
    }
    got.erase({pair.i, pair.j});
  }
  for (const auto &[key, f] : got) {
    if (f > 5.0 / n) {
      return fail(label, "sampled pair " + std::to_string(key.first) + "-" +
                             std::to_string(key.second) + " has no outside mass");
    }
  }
  return true;
}

bool check_windows(const std::string &label, const SampleResult &samples,
                   const StructuralProfileResult &profile,
                   const std::string &scaffold) {
  const int n = static_cast<int>(scaffold.size());
  const std::size_t count = samples.generated_pairs.size();
  for (std::size_t w_idx = 0; w_idx < profile.access_lengths.size(); ++w_idx) {
    const int w = profile.access_lengths[w_idx];
    std::vector<double> freq(n + 2, 0.0);
    for (const auto &structure : samples.union_structures) {
      for (int s = 1; s + w - 1 <= n; ++s) {
        bool free = true;
        for (int t = s; t < s + w && free; ++t) free = structure[t - 1] == '.';
        if (free) freq[s] += 1.0 / count;
      }
    }
    for (int s = 1; s + w - 1 <= n; ++s) {
      const double p = profile.accessibility[w_idx][s];
      if (!within(freq[s], p, count)) {
        return fail(label, "window w=" + std::to_string(w) + " s=" +
                               std::to_string(s) + " sampled " +
                               std::to_string(freq[s]) + " profile " +
                               std::to_string(p));
      }
    }
  }
  return true;
}

char complement(char c) {
  switch (c) {
    case 'A': return 'U';
    case 'U': return 'A';
    case 'G': return 'C';
    default: return 'G';
  }
}

std::string random_scaffold(std::mt19937 &rng, std::string &sequence) {
  const int n = static_cast<int>(sequence.size());
  std::string g(n, '.');
  auto rand_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, std::max(lo, hi))(rng);
  };
  const int hairpins = rand_int(1, 2);
  const int width = n / hairpins;
  for (int h = 0; h < hairpins; ++h) {
    const int a = h * width + rand_int(0, 2);
    const int b = (h + 1) * width - 1 - rand_int(0, 2);
    for (int t = 0; t < rand_int(1, 3); ++t) {
      if (b - a - 2 * t < 6) break;
      g[a + t] = '(';
      g[b - t] = ')';
      sequence[b - t] = complement(sequence[a + t]);
    }
  }
  return g;
}

}  // namespace

int main(int argc, char **argv) {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  if (!check_exact_tape_reference() || !check_exact_empty_and_unreachable())
    return EXIT_FAILURE;
  const auto invalid_width = beam_sample("AAAAAA", "......", 2, 7, BeamConfig{0});
  if (invalid_width.ok || invalid_width.error != "beam size must be positive") {
    fail("public-api", "zero beam width was not rejected");
    return EXIT_FAILURE;
  }
  const auto unit_width = beam_sample("AAAAAA", "......", 2, 7, BeamConfig{1});
  if (!unit_width.ok || unit_width.union_structures != std::vector<std::string>(2, "......")) {
    fail("public-api", "positive beam width lost the all-unpaired ensemble");
    return EXIT_FAILURE;
  }
  const int cases = argc > 1 ? std::atoi(argv[1]) : 6;
  const std::size_t draws = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 4000;
  std::mt19937 rng(20261008);
  std::uniform_int_distribution<int> base(0, 3);
  std::uniform_int_distribution<int> length(16, 40);
  bool ok = true;
  int checked = 0;
  std::size_t pk_samples = 0;
  for (int c = 0; c < cases; ++c) {
    std::string sequence(length(rng), 'A');
    for (auto &ch : sequence) ch = "ACGU"[base(rng)];
    const std::string scaffold = random_scaffold(rng, sequence);
    const std::string label = "case" + std::to_string(c) + " " + sequence + " " + scaffold;

    ExactSession session(sequence, scaffold);
    const auto outside = session.outside();
    if (!outside.ok || !outside.reachable) continue;

    // 1. Exact.
    const auto exact = exact_sample(sequence, scaffold, draws, 7 + c);
    if (!exact.ok) {
      ok = fail(label, "exact_sample: " + exact.error) && ok;
      continue;
    }
    if (std::fabs(exact.log_partition_weight - outside.log_partition_weight) > 1e-9) {
      ok = fail(label, "exact sampler log Z differs") && ok;
    }
    ++checked;
    for (const auto &structure : exact.union_structures) {
      if (structure.find('[') != std::string::npos) ++pk_samples;
    }
    ok = check_pairs(label + " exact", exact, outside.pair_probabilities, scaffold) && ok;
    std::map<std::string, double> structure_freq;
    std::map<std::string, std::vector<BasePair>> structure_pairs;
    for (std::size_t s = 0; s < exact.union_structures.size(); ++s) {
      structure_freq[exact.union_structures[s]] += 1.0 / draws;
      structure_pairs[exact.union_structures[s]] = exact.generated_pairs[s];
    }
    for (const auto &[structure, f] : structure_freq) {
      if (f < 0.02) continue;
      const auto prob = session.probability(structure_pairs[structure]);
      if (!prob.ok || !prob.reachable) {
        ok = fail(label, "sampled structure not reachable: " + structure) && ok;
        continue;
      }
      if (!within(f, prob.probability, draws)) {
        ok = fail(label, "structure " + structure + " sampled " +
                             std::to_string(f) + " exact " +
                             std::to_string(prob.probability)) && ok;
      }
    }

    // 2. KeepAll beam.
    BeamConfig keep_all;
    keep_all.beam_size = std::numeric_limits<std::size_t>::max();
    keep_all.max_hairpin_unpaired = std::numeric_limits<std::size_t>::max();
    keep_all.max_case3_ancestor_depth = std::numeric_limits<std::size_t>::max();
    const auto full = beam_sample(sequence, scaffold, draws, 11 + c, keep_all);
    if (!full.ok) {
      ok = fail(label, "keep-all beam_sample: " + full.error) && ok;
    } else {
      ok = check_pairs(label + " keep-all", full, outside.pair_probabilities, scaffold) && ok;
    }

    // 3. Finite beams.
    for (const std::size_t b : {3, 10}) {
      BeamConfig config;
      config.beam_size = b;
      const auto samples = beam_sample(sequence, scaffold, draws, 13 + c, config);
      const auto beam_out = beam_outside(sequence, scaffold, config);
      if (!samples.ok || !beam_out.ok || !beam_out.reachable) continue;
      const std::string tag = label + " b=" + std::to_string(b);
      ok = check_pairs(tag, samples, beam_out.pair_probabilities, scaffold) && ok;
      const auto profile =
          beam_structural_profiles(sequence, scaffold, {1, 2, 4, 7}, config);
      if (!profile.ok) {
        ok = fail(tag, "profiles: " + profile.error) && ok;
        continue;
      }
      ok = check_windows(tag, samples, profile, scaffold) && ok;
    }
  }
  std::cout << "checked_cases=" << checked
            << " exact_samples_with_crossing_pairs=" << pk_samples << "\n";
  if (!ok || checked == 0) return EXIT_FAILURE;
  std::cout << "scfg2_pk_sampling=ok\n";
  return EXIT_SUCCESS;
}

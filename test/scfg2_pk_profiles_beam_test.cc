// Beam PK-CapR / PK-Raccess.
//   1. KeepAll widths: beam_profiles() equals exact_profiles() (same DAG).
//   2. Finite widths: every value is a probability, each base's contexts
//      sum to one, and accessibility does not grow with the window length.

#include "scfg2/api/beam.hh"
#include "scfg2/api/profiles.hh"
#include "scfg2/profiles/beam_profiles.hh"
#include "scfg2/profiles/exact_profiles.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

using cparty::scfg2::profiles::ProfileResult;
using cparty::scfg2::profiles::loop_context_count;

bool fail(const std::string &label, const std::string &message) {
  std::cerr << "[" << label << "] " << message << "\n";
  return false;
}

bool same(const std::string &label, const ProfileResult &got,
          const ProfileResult &want, double tol) {
  if (got.accessibility_by_loop.size() != want.accessibility_by_loop.size())
    return fail(label, "loop breakdown availability differs");
  if (std::fabs(got.log_partition_weight - want.log_partition_weight) > 1e-9) {
    return fail(label, "log Z differs");
  }
  for (std::size_t pos = 1; pos < want.context.size(); ++pos) {
    for (std::size_t c = 0; c < loop_context_count; ++c) {
      if (std::fabs(got.context[pos][c] - want.context[pos][c]) > tol) {
        std::string row;
        for (std::size_t d = 0; d < loop_context_count; ++d) {
          row += " " + std::to_string(d) + ":" + std::to_string(got.context[pos][d]) +
                 "/" + std::to_string(want.context[pos][d]);
        }
        std::cerr << "  row" << row << "\n";
        return fail(label, "context " + std::to_string(c) + " at " +
                               std::to_string(pos) + ": beam " +
                               std::to_string(got.context[pos][c]) + " exact " +
                               std::to_string(want.context[pos][c]));
      }
    }
  }
  for (std::size_t w = 0; w < want.accessibility.size(); ++w) {
    for (std::size_t s = 1; s < want.accessibility[w].size(); ++s) {
      if (std::fabs(got.accessibility[w][s] - want.accessibility[w][s]) > tol) {
        return fail(label, "accessibility w=" +
                               std::to_string(want.window_lengths[w]) + " s=" +
                               std::to_string(s) + ": beam " +
                               std::to_string(got.accessibility[w][s]) +
                               " exact " +
                               std::to_string(want.accessibility[w][s]));
      }
      if (!want.accessibility_by_loop.empty())
        for (std::size_t c = 0; c < cparty::scfg2::profiles::unpaired_loop_context_count; ++c)
          if (std::fabs(got.accessibility_by_loop[w][s][c] -
                        want.accessibility_by_loop[w][s][c]) > tol)
            return fail(label, "loop accessibility differs");
    }
  }
  return true;
}

bool sane(const std::string &label, const ProfileResult &got) {
  constexpr double kTol = 1e-9;
  for (std::size_t pos = 1; pos < got.context.size(); ++pos) {
    double sum = 0.0;
    for (std::size_t c = 0; c < loop_context_count; ++c) {
      const double p = got.context[pos][c];
      if (!(p >= -kTol && p <= 1.0 + kTol)) {
        return fail(label, "context out of [0,1] at " + std::to_string(pos));
      }
      sum += p;
    }
    if (std::fabs(sum - 1.0) > kTol) return fail(label, "contexts do not sum to 1");
  }
  for (std::size_t w = 0; w < got.accessibility.size(); ++w) {
    for (std::size_t s = 1; s < got.accessibility[w].size(); ++s) {
      const double p = got.accessibility[w][s];
      if (!(p >= -kTol && p <= 1.0 + kTol)) {
        return fail(label, "accessibility out of [0,1]");
      }
      if (!got.accessibility_by_loop.empty()) {
        double sum = 0.0;
        for (std::size_t c = 0; c < cparty::scfg2::profiles::unpaired_loop_context_count; ++c) {
          const double q = got.accessibility_by_loop[w][s][c];
          if (!(q >= -kTol && q <= p + kTol))
            return fail(label, "loop accessibility outside [0,total]");
          if (got.window_lengths[w] == 1 &&
              std::fabs(q - got.context[s][c + 1]) > kTol)
            return fail(label, "length-one loop accessibility differs from CapR");
          sum += q;
        }
        if (std::fabs(sum - p) > kTol)
          return fail(label, "loop probabilities do not sum to total");
      }
      if (w > 0 && got.window_lengths[w] > got.window_lengths[w - 1] &&
          s < got.accessibility[w - 1].size() &&
          p > got.accessibility[w - 1][s] + kTol) {
        return fail(label, "accessibility grows with window length at " +
                               std::to_string(s));
      }
    }
    // A single base's accessibility is its unpaired probability.
    if (got.window_lengths[w] == 1) {
      for (std::size_t s = 1; s < got.accessibility[w].size(); ++s) {
        const double stem = got.context[s][0];
        if (std::fabs(got.accessibility[w][s] - (1.0 - stem)) > kTol) {
          return fail(label, "P_1 != 1 - P_stem at " + std::to_string(s));
        }
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

// Two or three hairpin stems, sometimes inside an outer stem.
std::string random_scaffold(std::mt19937 &rng, std::string &sequence) {
  const int n = static_cast<int>(sequence.size());
  std::string g(n, '.');
  auto rand_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, std::max(lo, hi))(rng);
  };
  auto place = [&](int a, int b, int len) {
    for (int t = 0; t < len; ++t) {
      if (b - a - 2 * t < 4) return;
      g[a + t] = '(';
      g[b - t] = ')';
      sequence[b - t] = complement(sequence[a + t]);
    }
  };
  int lo = 0, hi = n - 1;
  if (rand_int(0, 2) == 0) {
    place(0, n - 1, 2);
    lo = 3;
    hi = n - 4;
  }
  const int hairpins = rand_int(2, 3);
  const int width = (hi - lo + 1) / hairpins;
  for (int h = 0; h < hairpins; ++h) {
    const int a = lo + h * width + rand_int(0, 2);
    const int b = lo + (h + 1) * width - 1 - rand_int(0, 2);
    if (b - a >= 6) place(a, b, rand_int(2, 4));
  }
  return g;
}

bool check_streamed_whole_runs() {
  using namespace cparty::scfg2::profiles;
  constexpr cand_pos_t n = 257;
  const std::vector<cand_pos_t> windows{0, 1, 2, 19, 128, n, n + 10, 2};
  ProfileAccumulator stored(n, true), streamed(n, windows, true);
  std::mt19937 rng(20261009);
  for (int k = 0; k < 20000; ++k) {
    UnpairedEmission emission;
    emission.count = 1 + k % 2;
    for (std::size_t r = 0; r < emission.count; ++r) {
      const auto lo = static_cast<cand_pos_t>(1 + rng() % n);
      const auto hi = static_cast<cand_pos_t>(lo + rng() % (n - lo + 1));
      emission.runs[r] = UnpairedRun{lo, hi, static_cast<LoopContext>(1 + rng() % 6)};
    }
    const double probability = 1e-7 * (1 + rng() % 5);
    stored.add_whole_runs(emission, probability);
    streamed.add_whole_runs(emission, probability);
  }
  if (stored.stored_whole_runs() != 20000 || streamed.stored_whole_runs() != 0)
    return fail("streamed-runs", "regenerated emissions were saved as branch records");
  if (!same("streamed-runs", streamed.finish(0.0, windows), stored.finish(0.0, windows), 1e-12)) return false;
  if (streamed.finish(0.0, {1}).ok)
    return fail("streamed-runs", "changed windows were silently accepted");
  return true;
}

bool check_chain_window_prefixes() {
  using namespace cparty::scfg2::profiles;
  for (const bool mirrored : {false, true}) {
    ProfileAccumulator acc(40, true);
    std::vector<std::pair<cand_pos_t, ProfileAccumulator::ChainStep>> steps;
    for (cand_pos_t q = 0; q < 40; ++q) {
      if (q == 21) continue;
      const cand_pos_t pos = mirrored ? 40 - q : q + 1;
      auto &step = acc.chain_step(7, 1, pos, LoopContext::Exterior, mirrored);
      step.log_entry = q == 22 ? -std::numeric_limits<double>::infinity() : -2.0 - 0.2*q;
      step.log_step = q == 31 ? -std::numeric_limits<double>::infinity() : -0.06*(1 + q%3);
      step.has_prev = q != 0 && q != 13;
      step.log_state = -0.1*q;
      step.state_adjoint = 0.1 + 0.01*q;
      steps.push_back({mirrored ? -pos : pos, step});
    }
    const std::vector<cand_pos_t> windows{1, 2, 9, 17, 40};
    const auto result = acc.finish(0.0, windows);
    for (std::size_t wi = 0; wi < windows.size(); ++wi) {
      const auto w = windows[wi];
      std::vector<double> want(41);
      for (std::size_t s = 0; s < steps.size(); ++s) {
        if (!std::isfinite(steps[s].second.log_entry)) continue;
        double log_mass = steps[s].second.log_entry;
        bool linked = s + w <= steps.size();
        for (cand_pos_t t = 1; t < w && linked; ++t) {
          const auto &next = steps[s+t];
          linked = next.first == steps[s].first + t && next.second.has_prev;
          log_mass += next.second.log_step;
        }
        if (!linked) continue;
        const auto &last = steps[s+w-1];
        const auto pos = mirrored ? -last.first : steps[s].first;
        want[pos] += std::exp(log_mass - last.second.log_state) * last.second.state_adjoint;
      }
      for (cand_pos_t pos = 1; pos <= 40; ++pos) {
        if (std::fabs(want[pos] - result.accessibility[wi][pos]) > 1e-12)
          return fail("chain-prefix", "window posterior differs from direct product");
        for (std::size_t c = 0; c < unpaired_loop_context_count; ++c)
          if (std::fabs(result.accessibility_by_loop[wi][pos][c] -
                        (c == 4 ? want[pos] : 0.0)) > 1e-12)
            return fail("chain-prefix", "loop posterior differs from direct product");
      }
    }
    if (result.chain_window_steps != steps.size() * windows.size())
      return fail("chain-prefix", "work depends on window length");
  }
  const cand_pos_t n = 4096;
  ProfileAccumulator all_free(n, true);
  for (cand_pos_t t = 1; t <= n; ++t) {
    auto &step = all_free.chain_step(7, 1, t, LoopContext::Exterior);
    step.log_entry = step.log_state = -0.2*t;
    step.log_step = -0.2;
    step.has_prev = t > 1;
    step.state_adjoint = 1.0;
    step.emit_probability = 1.0;
  }
  const auto result = all_free.finish(0.0, {1, n/2, n});
  if (result.chain_window_steps != 3 * std::size_t(n))
    return fail("long-chain-prefix", "long windows exceed linear step count");
  for (std::size_t wi = 0; wi < result.window_lengths.size(); ++wi)
    for (cand_pos_t s = 1; s + result.window_lengths[wi] - 1 <= n; ++s)
      if (std::fabs(result.accessibility[wi][s] - 1.0) > 1e-12)
        return fail("long-chain-prefix", "all-free window probability differs from one");
      else if (std::fabs(result.accessibility_by_loop[wi][s][4] - 1.0) > 1e-12)
        return fail("long-chain-prefix", "long window lost its exterior label");
  return true;
}

}  // namespace

int main(int argc, char **argv) {
  using namespace cparty::scfg2;
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  const auto invalid_width = api::beam_structural_profiles(
      "AAAAAA", "......", {1, 2}, api::BeamConfig{0});
  if (invalid_width.ok || invalid_width.error != "beam size must be positive") {
    fail("public-api", "zero beam width was not rejected");
    return EXIT_FAILURE;
  }
  const auto unit_width = api::beam_structural_profiles(
      "AAAAAA", "......", {1, 2}, api::BeamConfig{1});
  if (!unit_width.ok || std::fabs(unit_width.accessibility[1][1] - 1.0) > 1e-12) {
    fail("public-api", "positive beam width lost the all-unpaired ensemble");
    return EXIT_FAILURE;
  }
  if (!unit_width.accessibility_by_loop.empty()) {
    fail("public-api", "opt-out call allocated loop windows");
    return EXIT_FAILURE;
  }
  const auto unit_by_loop = api::beam_structural_profiles(
      "AAAAAA", "......", {1, 2}, api::BeamConfig{1}, true);
  if (!unit_by_loop.ok || unit_by_loop.accessibility_by_loop.size() != 2 ||
      std::fabs(unit_by_loop.accessibility_by_loop[1][1][4] - 1.0) > 1e-12) {
    fail("public-api", "loop accessibility missing from API");
    return EXIT_FAILURE;
  }
  const int count = argc > 1 ? std::atoi(argv[1]) : 10;
  const int max_length = argc > 2 ? std::atoi(argv[2]) : 50;
  const std::vector<cand_pos_t> windows = {1, 2, 4, 8};
  std::mt19937 rng(20261008);
  std::uniform_int_distribution<int> base(0, 3);
  std::uniform_int_distribution<int> length(20, max_length);
  bool ok = check_chain_window_prefixes() && check_streamed_whole_runs();
  // Band-wing witnesses exercise WIP/VPL starts, VPR rightward prefixes,
  // mirrored leftward prefixes and nested WI/WMB emissions. Compare independent
  // reverse engines on each retained support; KeepAll also matches exact.
  struct Fixture { const char *name; const char *sequence; const char *scaffold; };
  const Fixture fixtures[] = {
    {"gallery", "GGGGGAACCAACCC", "((.....))....."},
    {"be-left-branch-right-run", "GAGAAACAAGAAGAAAAAAAAACAAACAAAAC", "(.(...)..(............)...)....."},
    {"be-right-branch-left-run", "CAAAACAAACAAAAAAAAAGAAGAACAAAGAG", ".....(...(............)..(...).)"},
    {"vp-left-branch-right-run", "GAGAAACAAGAAGAAAAAAAAACAAACAAAAC", "..(...).....(..................)"},
    {"vp-right-branch-left-run", "CAAAACAAACAAAAAAAAAGAAGAACAAAGAG", "(..................).....(...).."},
    {"be-wing-with-nested-wmb", "GAGGGGGAACCAACCCAAGAAGAAAAAAAAACAAACAAAAC", "(.((.....)).......(............)...)....."},
    {"wi-gap-with-nested-wing", "GGGGGGAGAAACAAGAAGAAAAAAAAACAAACAAAACACCAACCC", "((...(.(...)..(............)...)......))....."}
  };
  for (const auto &fixture : fixtures) {
    for (const auto policy : {WmbClosurePolicy::RequireWeakClosure,
                              WmbClosurePolicy::LegacyCompatible}) {
      const auto exact = profiles::exact_profiles(
          fixture.sequence, fixture.scaffold, windows, policy, true);
      if (!exact.ok) { ok = fail(fixture.name, exact.error) && ok; continue; }
      for (const auto ranking : {beam::RankingMode::Baseline, beam::RankingMode::Outer}) {
        for (const std::size_t b : {0, 2, 5, 20}) {
          // Default low-level widths have no beam or geometry limits.
          const auto widths = b ? replay::WFinalVpWmbpOnlineWidths::uniform(b)
                                : replay::WFinalVpWmbpOnlineWidths{};
          const auto numeric = profiles::beam_profiles(
              fixture.sequence, fixture.scaffold, widths, windows, policy,
              ranking, beam::PilotExponent::TwoThirds, true);
          const auto taped = profiles::beam_profiles_taped(
              fixture.sequence, fixture.scaffold, widths, windows, policy,
              ranking, beam::PilotExponent::TwoThirds, true);
          const auto label = std::string(fixture.name) + " b=" + std::to_string(b);
          if (!numeric.ok || !taped.ok)
            ok = fail(label, numeric.error + " / " + taped.error) && ok;
          else {
            ok = same(label, numeric, taped, 1e-9) && sane(label, numeric) && ok;
            if (b == 0) ok = same(label + " exact", numeric, exact, 1e-9) && ok;
          }
        }
      }
    }
  }
  // A missing aggregate rule must not silently acquire a context.
  try {
    (void)profiles::detail::aggregated_run_context(grammar::RuleId::V_HAIRPIN);
    ok = fail("aggregate-context", "accepted a non-aggregate rule") && ok;
  } catch (const std::logic_error &) {}
  try {
    (void)profiles::fixed_rule_context(grammar::RuleId::V_INTERNAL);
    ok = fail("fixed-context", "accepted a geometry-dependent rule") && ok;
  } catch (const std::logic_error &) {}
  for (int c = 0; c < count; ++c) {
    std::string sequence(length(rng), 'A');
    for (auto &ch : sequence) ch = "ACGU"[base(rng)];
    const std::string scaffold = random_scaffold(rng, sequence);
    const std::string label = "case" + std::to_string(c) + " " + sequence + " " + scaffold;
    const auto exact = profiles::exact_profiles(sequence, scaffold, windows,
        WmbClosurePolicy::RequireWeakClosure, true);
    if (!exact.ok) {
      std::cout << "[" << label << "] skip: " << exact.error << "\n";
      continue;
    }
    const auto keep_all = profiles::beam_profiles(
        sequence, scaffold, replay::WFinalVpWmbpOnlineWidths{}, windows,
        WmbClosurePolicy::RequireWeakClosure, beam::RankingMode::Baseline,
        beam::PilotExponent::TwoThirds, true);
    if (!keep_all.ok) {
      ok = fail(label, "keep-all beam: " + keep_all.error) && ok;
      continue;
    }
    if (keep_all.reverse_records == 0 || keep_all.reverse_terms != 0 || keep_all.reverse_factors != 0)
      ok = fail(label, "regenerated profiles saved a branch forest") && ok;
    ok = same(label + " keep-all", keep_all, exact, 1e-9) && ok;
    ok = sane(label + " exact", exact) && ok;
    for (const std::size_t b : {2, 5, 20}) {
      const auto finite = profiles::beam_profiles(
          sequence, scaffold, replay::WFinalVpWmbpOnlineWidths::uniform(b),
          windows, WmbClosurePolicy::RequireWeakClosure,
          beam::RankingMode::Baseline, beam::PilotExponent::TwoThirds, true);
      if (!finite.ok) continue;
      ok = sane(label + " b=" + std::to_string(b), finite) && ok;
      // Unpaired probability from emissions equals the pair-based value of
      // the existing beam outside on the same widths.
      api::BeamConfig config;
      config.beam_size = b;
      config.ranking = api::BeamRanking::Baseline; // compare the baseline low-level profile above
      const auto outside = api::beam_outside(sequence, scaffold, config);
      if (!outside.ok) {
        ok = fail(label, "beam_outside: " + outside.error) && ok;
        continue;
      }
      for (std::size_t pos = 1; pos < finite.context.size(); ++pos) {
        const double unpaired = 1.0 - finite.context[pos][0];
        if (std::fabs(unpaired - outside.unpaired_probabilities[pos]) > 1e-9) {
          ok = fail(label + " b=" + std::to_string(b),
                    "unpaired at " + std::to_string(pos) + ": emissions " +
                        std::to_string(unpaired) + " pairs " +
                        std::to_string(outside.unpaired_probabilities[pos])) &&
               ok;
          break;
        }
      }
    }
  }
  if (!ok) return EXIT_FAILURE;
  std::cout << "scfg2_pk_profiles_beam=ok\n";
  return EXIT_SUCCESS;
}

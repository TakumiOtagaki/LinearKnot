#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/beam.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {
using namespace cparty::scfg2;

struct Witness {
  const char *file;
  int j;
  double legacy_energy;
  double legacy_target_count;
};

bool check(const Witness &w) {
  std::ifstream input(std::string("test/fixtures/scfg2_wmbp_leftward_be_20261008/") + w.file);
  std::string sequence, scaffold;
  std::getline(input, sequence);
  std::getline(input, scaffold);
  if (sequence.empty() || sequence.size() != scaffold.size()) return false;
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const auto fail = [&](const std::string &what) {
    std::cerr << w.file << ": " << what << '\n';
    return false;
  };

  // Goldens are direct patched-legacy inside and an external V-marker
  // measurement at MFE seed 0 / pf_scale=1, independent of SCFG2 scheduling.
  // Remeasured with multiloop closing patch0009; see the 2026-10-17 note.
  api::ExactSession exact(sequence, scaffold);
  const auto inside = exact.inside();
  if (!inside.ok || !std::isfinite(inside.ensemble_free_energy) ||
      std::fabs(inside.ensemble_free_energy - w.legacy_energy) > 2e-7)
    return fail("inside disagrees with patched legacy: " +
                std::to_string(inside.ensemble_free_energy));

  using S = semirings::FactorizedLogSumExpSemiring;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  semirings::FactorizedInsideGraph graph;
  const auto chart = [&] {
    S::Scope scope(graph);
    return replay::run_w_final_exact_inside_schedule<S>(
        owner, engine::build_pkfree_basic_exact_schedule(n, ctx));
  }();
  const auto root = chart.get({grammar::NonTerminal::W, 1, n, -1, -1});
  const auto target = chart.get({grammar::NonTerminal::V, 27, w.j, -1, -1});
  if (S::is_zero(root) || S::is_zero(target)) return fail("missing audit state");
  const auto posterior = graph.posteriors(root.node);
  if (std::fabs(posterior[target.node] - w.legacy_target_count) > 1e-9)
    return fail("target derivative disagrees with external legacy marker");

  // A count below one alone can conceal a rare duplicated derivation.
  // Check the maximum count over every positive term in the exact forest.
  std::vector<unsigned char> maximum(graph.states().size());
  for (std::size_t node = 0; node < maximum.size(); ++node) {
    if (node == target.node) { maximum[node] = 1; continue; }
    const auto &state = graph.states()[node];
    for (std::size_t t = state.first_term;
         t < std::size_t(state.first_term) + state.term_count; ++t) {
      const auto &term = graph.terms()[t];
      unsigned count = 0;
      for (std::size_t f = term.first_factor;
           f < std::size_t(term.first_factor) + term.factor_count; ++f)
        count = std::min(2u, count + maximum[graph.factors()[f]]);
      maximum[node] = std::max<unsigned>(maximum[node], count);
    }
  }
  if (maximum[root.node] != 1) return fail("positive derivation repeats target V");

  const auto outside = exact.outside();
  if (!outside.ok) return fail("exact outside: " + outside.error);
  // A saturated finite beam exercises retained joins and shared prefix
  // caches; it must match exact with the same seed/geometry universe.
  const api::BeamConfig config{std::size_t(n) + 1, std::size_t(n), std::size_t(n)};
  api::BeamSession beam(sequence, scaffold, config);
  const auto beam_inside = beam.inside();
  const auto beam_outside = beam.outside();
  if (!beam_inside.ok || !beam_outside.ok ||
      std::fabs(beam_inside.log_partition_weight - S::priority(root)) > 1e-9 ||
      std::fabs(beam_outside.log_partition_weight - S::priority(root)) > 1e-9)
    return fail("saturated beam inside/outside: " + beam_outside.error);
  for (const auto *pairs : {&outside.pair_probabilities, &beam_outside.pair_probabilities}) {
    bool found = false;
    for (const auto &pair : *pairs) {
      if (pair.i != 27 || pair.j != w.j) continue;
      found = true;
      if (std::fabs(pair.union_probability - w.legacy_target_count) > 1e-9 ||
          pair.square_probability != 0.0)
        return fail("public target marginal differs from legacy marker");
    }
    if (!found) return fail("public target marginal absent");
  }

  const auto best = exact.viterbi();
  const auto beam_best = beam.viterbi();
  if (!best.ok || !beam_best.ok || !std::isfinite(best.energy) ||
      std::fabs(best.energy - beam_best.energy) > 1e-9 ||
      best.weight > inside.partition_weight * (1.0 + 1e-9))
    return fail("shared Viterbi surface disagrees with saturated beam/inside");
  return true;
}
}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  const Witness witnesses[] = {
      {"stem20.txt", 69, -37.78103825737757, 0.999499243324801},
      {"stem35.txt", 99, -87.2855382007357, 0.9995019586334782},
      {"stem50.txt", 129, -136.79001320869293, 0.9995046337224037},
  };
  for (const auto &w : witnesses) if (!check(w)) return EXIT_FAILURE;
  std::cout << "scfg2_wmbp_be_availability=ok\n";
}

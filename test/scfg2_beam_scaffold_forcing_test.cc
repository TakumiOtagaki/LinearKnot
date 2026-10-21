// A finite beam must never lose a scaffold pair.  Every structure contains
// all pairs of G, so the beam keeps one scaffold-only derivation (G' empty)
// outside the width.  For random scaffolds, including hairpin loops longer
// than the beam hairpin cap and wide multiloops, and for very small widths,
// check that the root is reachable and that every scaffold pair has a
// nonzero V item.

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/run_term_observer.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

char complement(char c) {
  switch (c) {
    case 'A': return 'U';
    case 'U': return 'A';
    case 'G': return 'C';
    default: return 'G';
  }
}

void stem(std::string &g, std::string &seq, int a, int b, int len) {
  for (int t = 0; t < len && b - a - 2 * t >= 4; ++t) {
    g[a + t] = '(';
    g[b - t] = ')';
    seq[b - t] = complement(seq[a + t]);
  }
}

// An outer stem closing a multiloop of `branches` hairpins; one hairpin has
// a loop of `long_loop` bases (beyond the beam hairpin cap of 30).
std::string scaffold(std::mt19937 &rng, std::string &seq, int branches, int long_loop) {
  const int n = static_cast<int>(seq.size());
  std::string g(n, '.');
  stem(g, seq, 0, n - 1, 3);
  int pos = 5;
  for (int k = 0; k < branches; ++k) {
    const int gap = std::uniform_int_distribution<int>(2, 9)(rng);
    const int loop = k == 0 ? long_loop : std::uniform_int_distribution<int>(3, 12)(rng);
    const int len = std::uniform_int_distribution<int>(2, 5)(rng);
    const int a = pos + gap, b = a + 2 * len + loop - 1;
    if (b >= n - 6) break;
    stem(g, seq, a, b, len);
    pos = b + 1;
  }
  return g;
}

void check_forced_order(const cparty::scfg2::beam::ItemKeySet &actual,
                        const cparty::scfg2::beam::ItemKeySet &reference) {
  const cparty::scfg2::deductions::ItemKeyEq equal;
  if (actual.size() != reference.size() || actual.bucket_count() != reference.bucket_count()) {
    std::cerr << "forcing changed forced item count/buckets\n";
    std::exit(EXIT_FAILURE);
  }
  auto expected = reference.begin();
  for (const auto &item : actual) {
    if (expected == reference.end() || !equal(item, *expected++)) {
      std::cerr << "forcing changed forced full keys/iteration order\n";
      std::exit(EXIT_FAILURE);
    }
  }
  for (std::size_t bucket = 0; bucket < actual.bucket_count(); ++bucket) {
    if (actual.bucket_size(bucket) != reference.bucket_size(bucket)) {
      std::cerr << "forcing changed forced bucket occupancy\n";
      std::exit(EXIT_FAILURE);
    }
    auto expected_item = reference.begin(bucket);
    for (auto item = actual.begin(bucket); item != actual.end(bucket); ++item)
      if (!equal(*item, *expected_item++)) {
        std::cerr << "forcing changed forced local bucket order\n";
        std::exit(EXIT_FAILURE);
      }
  }
}

// Memory-observed forcing retains the original complete choice store. A
// no-op observer is therefore an independent public reference for the plain
// private forcing optimization; every other context/width setting is identical.
template <class Context>
auto forcing_items_with_observer_check(
    Context &ctx, cand_pos_t n,
    const cparty::scfg2::replay::WFinalVpWmbpOnlineWidths &widths,
    std::size_t &checked) {
  namespace sc = cparty::scfg2;
  const auto plain = sc::replay::scaffold_only_derivation_items(ctx, n, widths);
  auto observed_widths = widths;
  observed_widths.memory_observer = [](const sc::beam::MemoryObservation &) {};
  const auto observed = sc::replay::scaffold_only_derivation_items(ctx, n, observed_widths);
  check_forced_order(*plain, *observed);
  ++checked;
  return plain;
}

// One ordered stream covers all observer callback kinds, including inactive
// alias fields. Copy values during the callback; never retain borrowed aliases.
struct RunTerms : cparty::scfg2::beam::RunTermObserver<double> {
  using Rule = cparty::scfg2::grammar::RuleId;
  using Event = std::tuple<int, Rule, cand_pos_t, cand_pos_t,
                           std::uint64_t, std::uint64_t, bool, bool>;
  std::vector<Event> events;
  static std::uint64_t bits(double value) {
    static_assert(sizeof(value) == sizeof(std::uint64_t));
    std::uint64_t out;
    std::memcpy(&out, &value, sizeof(out));
    return out;
  }
  void on_run_term(Rule rule, cand_pos_t lo, cand_pos_t hi, const double &term) override {
    events.emplace_back(0, rule, lo, hi, bits(term), 0, false, false);
  }
  void on_vpr_run_cell(cand_pos_t i, cand_pos_t t, const double &total, bool carried) override {
    events.emplace_back(1, Rule{}, i, t, bits(total), 0, false, carried);
  }
  void on_left_run_unit(Rule rule, cand_pos_t j, cand_pos_t base,
                        const double &unit, const double *alias, bool linked) override {
    events.emplace_back(2, rule, j, base, bits(unit), alias ? bits(*alias) : 0,
                        alias != nullptr, linked);
  }
};

void numeric_forcing_checks(cparty::scfg2::replay::WFinalLogWeightContext &ctx,
                            cand_pos_t n,
                            cparty::scfg2::replay::WFinalVpWmbpOnlineWidths widths,
                            std::size_t &checked) {
  namespace sc = cparty::scfg2;
  widths.numeric_lookup_index = false;
  const auto reference = sc::replay::scaffold_only_derivation_items(ctx, n, widths);
  auto numeric_widths = widths;
  numeric_widths.numeric_lookup_index = true;
  const auto actual = sc::replay::scaffold_only_derivation_items(ctx, n, numeric_widths);
  check_forced_order(*actual, *reference);
  auto observed_widths = numeric_widths;
  bool observed_full_family = false;
  observed_widths.memory_observer = [&](const sc::beam::MemoryObservation &event) {
    if (std::strcmp(event.event, "forward_end") != 0) return;
    for (const auto &family : event.families)
      if (family.nonterminal == sc::grammar::NonTerminal::WI && family.live_items)
        observed_full_family = true;
  };
  const auto observed = sc::replay::scaffold_only_derivation_items(ctx, n, observed_widths);
  check_forced_order(*observed, *reference);
  if (reference->empty() || !observed_full_family) {
    std::cerr << "numeric forcing/memory fallback witness is vacuous\n";
    std::exit(EXIT_FAILURE);
  }

  // Install the same collector for both lookup settings. An active observer
  // must preserve the complete forcing path and its callback values/order.
  RunTerms recorder;
  std::shared_ptr<const sc::beam::ItemKeySet> old_terms, new_terms;
  {
    sc::beam::ScopedRunTermObserver<double> scope(&recorder);
    old_terms = sc::replay::scaffold_only_derivation_items(ctx, n, widths);
  }
  auto expected_events = std::move(recorder.events);
  recorder.events.clear();
  {
    sc::beam::ScopedRunTermObserver<double> scope(&recorder);
    new_terms = sc::replay::scaffold_only_derivation_items(ctx, n, numeric_widths);
  }
  check_forced_order(*old_terms, *reference);
  check_forced_order(*new_terms, *reference);
  if (expected_events.empty() || recorder.events != expected_events) {
    std::cerr << "numeric forcing changed active run observer event order/value bits\n";
    std::exit(EXIT_FAILURE);
  }
  ++checked;
}

}  // namespace

int main(int argc, char **argv) {
  using namespace cparty::scfg2;
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  const int cases = argc > 1 ? std::atoi(argv[1]) : 12;
  // --no-force: run without the scaffold-only derivation, to show the test
  // detects the failure it guards against.
  const bool force = !(argc > 2 && std::string(argv[2]) == "--no-force");
  if (force) {
    // V stays KeepAll while WM is finite. The scaffold forcing probe uses
    // a huge finite width; reserving that many prefix slots used to throw.
    std::string seq = "GGGAAACCC", g = "(((...)))";
    cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
    replay::WFinalExactInsideContext ctx(owner);
    auto mixed = replay::WFinalVpWmbpOnlineWidths{};
    mixed.wm = 20;
    const auto reference = replay::run_w_final_vp_wmbp_online<semirings::LogSumExpSemiring>(
        ctx, seq.size(), replay::WFinalVpWmbpOnlineWidths{});
    const auto actual = replay::run_w_final_vp_wmbp_online<semirings::LogSumExpSemiring>(
        ctx, seq.size(), mixed);
    const deductions::ItemKey root{grammar::NonTerminal::W, 1,
        static_cast<cand_pos_t>(seq.size()), -1, -1};
    if (!std::isfinite(actual.chart.get(root)) ||
        std::fabs(actual.chart.get(root) - reference.chart.get(root)) > 1e-12) {
      std::cerr << "mixed KeepAll/finite scaffold forcing changed the saturated root\n";
      return EXIT_FAILURE;
    }
  }
  std::mt19937 rng(20261009);
  int checked = 0;
  double total_seconds = 0.0;
  int failures = 0;
  std::size_t maximum_forced_per_family_column = 0;
  std::size_t forcing_observer_checks = 0;
  for (int c = 0; c < cases; ++c) {
    const int n = std::uniform_int_distribution<int>(120, 320)(rng);
    std::string seq(n, 'A');
    for (auto &ch : seq) ch = "ACGU"[std::uniform_int_distribution<int>(0, 3)(rng)];
    std::string g = scaffold(rng, seq, std::uniform_int_distribution<int>(2, 6)(rng),
                             std::uniform_int_distribution<int>(31, 45)(rng));
    cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
    replay::WFinalExactInsideContext ctx(owner);
    if (force) {
      // The standard forced set is a single G'=empty traceback. Its pure
      // pk-free branches partition the sequence, so no family occurs twice
      // at the same right boundary. This is stronger than an arbitrary
      // caller-supplied forced set and bounds each forced family's width by b+1.
      const auto forcing_widths = replay::WFinalVpWmbpOnlineWidths::uniform(20);
      const auto forced_items = forcing_items_with_observer_check(
          ctx, n, forcing_widths, forcing_observer_checks);
      if (c == 0) {
        // Production log contexts use both maintained closure policies. The
        // same nonempty branching/long-loop scaffold also exercises retained
        // fixed VM choices and the unreachable WI/WIP/BE choice boundary.
        for (const auto policy : {WmbClosurePolicy::LegacyCompatible,
                                  WmbClosurePolicy::RequireWeakClosure}) {
          replay::WFinalLogWeightContext log_ctx(owner, policy);
          numeric_forcing_checks(log_ctx, n, forcing_widths, forcing_observer_checks);
        }
      }
      std::map<std::pair<grammar::NonTerminal, cand_pos_t>, std::size_t> counts;
      for (const auto &item : *forced_items) {
        using NT = grammar::NonTerminal;
        if (item.nonterminal != NT::W && item.nonterminal != NT::V &&
            item.nonterminal != NT::VM && item.nonterminal != NT::WM &&
            item.nonterminal != NT::WMv) {
          std::cerr << "scaffold-only traceback contains a PK family\n";
          return EXIT_FAILURE;
        }
        const auto count = ++counts[{item.nonterminal, item.j}];
        maximum_forced_per_family_column = std::max(maximum_forced_per_family_column, count);
        if (count > 1) {
          std::cerr << "multiple scaffold-only forced items in one family/column G=" << g << '\n';
          return EXIT_FAILURE;
        }
      }
    }
    for (const std::size_t b : {1, 2, 5, 20}) {
      auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(b);
      widths.force_scaffold_derivation = force;
      const auto started = std::chrono::steady_clock::now();
      const auto run = replay::run_w_final_vp_wmbp_online<semirings::LogSumExpSemiring>(
          ctx, n, widths);
      total_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - started).count();
      const auto value = [&](grammar::NonTerminal nt, cand_pos_t i, cand_pos_t j) {
        return run.chart.get(deductions::ItemKey{nt, i, j, -1, -1});
      };
      if (!std::isfinite(value(grammar::NonTerminal::W, 1, n))) {
        ++failures;
        if (!force) continue;
        std::cerr << "root unreachable b=" << b << " G=" << g << "\n";
        return EXIT_FAILURE;
      }
      for (cand_pos_t i = 1; i <= n; ++i) {
        const cand_pos_t j = ctx.pair_partner(i);
        if (j > i && !std::isfinite(value(grammar::NonTerminal::V, i, j))) {
          ++failures;
          if (!force) break;
          std::cerr << "scaffold pair (" << i << "," << j << ") lost at b=" << b
                    << " G=" << g << "\n";
          return EXIT_FAILURE;
        }
      }
      ++checked;
    }
  }
  std::cout << "checked_runs=" << checked << " failures=" << failures
            << " beam_seconds=" << total_seconds
            << " max_forced_per_family_column=" << maximum_forced_per_family_column
            << " forcing_observer_checks=" << forcing_observer_checks
            << "\nscfg2_beam_scaffold_forcing=ok\n";
  return EXIT_SUCCESS;
}

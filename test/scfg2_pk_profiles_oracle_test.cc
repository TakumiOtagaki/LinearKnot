// Oracle for PK-CapR contexts and PK-Raccess accessibility.
//
// Small inputs are solved by enumerating every derivation of the exact
// grammar.  For each derivation we record its generated pairs and the context
// each deduction assigns to the bases it emits.  The test checks
//   1. labels: every base unpaired in G u G' is emitted exactly once, and no
//      paired base is emitted;
//   2. contexts: exact_profiles() equals the derivation-weighted context
//      frequencies;
//   3. accessibility: exact_profiles() equals the probability that a window
//      is unpaired, computed from the pair sets alone.
//   4. loop accessibility: the entire window has one emitted loop label;
//      its weighted frequency equals the context-resolved result.

#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/profiles/exact_profiles.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;
using cparty::scfg2::profiles::LoopContext;

constexpr char kNone = '-';

std::map<std::string, int> &rule_coverage() {
  static std::map<std::string, int> coverage;
  return coverage;
}

std::array<bool, cparty::scfg2::profiles::unpaired_loop_context_count> &window_loop_coverage() {
  static std::array<bool, cparty::scfg2::profiles::unpaired_loop_context_count> coverage{};
  return coverage;
}

char context_code(LoopContext context) {
  return "SHBIMEP"[static_cast<int>(context)];
}

struct Summary {
  std::vector<cparty::scfg2::engine::TracebackPair> pairs;
  std::string emitted;  // one-based; kNone when not emitted
  double weight = 0.0;
};

using SummaryMap = std::map<std::string, Summary>;

std::string summary_key(const Summary &summary) {
  std::string key = summary.emitted + "|";
  auto pairs = summary.pairs;
  std::sort(pairs.begin(), pairs.end(), [](const auto &a, const auto &b) {
    return std::tie(a.i, a.j) < std::tie(b.i, b.j);
  });
  for (const auto &pair : pairs) {
    key += std::to_string(pair.i) + "," + std::to_string(pair.j) + ";";
  }
  return key;
}

void merge(SummaryMap *map, Summary summary) {
  const auto key = summary_key(summary);
  auto [it, inserted] = map->try_emplace(key, summary);
  if (!inserted) it->second.weight += summary.weight;
}

ItemKey runtime_child_item(const Deduction &deduction, std::size_t idx) {
  ItemKey child = deduction.children[idx].item;
  if (deduction.rule == cparty::scfg2::grammar::RuleId::WMBP_DIRECT_VP &&
      child.nonterminal == cparty::scfg2::grammar::NonTerminal::VP) {
    child.nonterminal = cparty::scfg2::grammar::NonTerminal::VP;
  }
  return child;
}

bool fail(const std::string &label, const std::string &message) {
  std::cerr << "[" << label << "] " << message << "\n";
  return false;
}

bool check_case(const std::string &label,
                std::string sequence,
                std::string scaffold) {
  using namespace cparty::scfg2;
  const cand_pos_t n = static_cast<cand_pos_t>(sequence.size());
  const std::vector<cand_pos_t> windows = {1, 2, 3, 5, n, n + 3, 2};

  const auto profile = profiles::exact_profiles(sequence, scaffold, windows,
      WmbClosurePolicy::RequireWeakClosure, true);
  if (!profile.ok) return fail(label, "exact_profiles: " + profile.error);
  const auto total_only = profiles::exact_profiles(sequence, scaffold, windows);
  if (!total_only.ok || !total_only.accessibility_by_loop.empty() ||
      total_only.accessibility != profile.accessibility ||
      total_only.context != profile.context)
    return fail(label, "loop breakdown changed default profiles");

  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  storage::exact_store::reset();
  const auto is_unpaired = [&](cand_pos_t pos) { return ctx.is_unpaired(pos); };
  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);
  std::unordered_map<ItemKey, SummaryMap, deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      chart;
  for (const auto &item : schedule) {
    SummaryMap groups;
    for (const auto &deduction :
         deductions::generate_exact_basic_deductions_for_item(item, ctx)) {
      if (!engine::legacy_span_schedule_supports(deduction)) continue;
      const double local = weight.weight_for(deduction);
      if (!(local > 0.0)) continue;
      Summary initial;
      initial.weight = local;
      initial.emitted.assign(static_cast<std::size_t>(n) + 1, kNone);
      engine::append_deduction_generated_pairs(deduction, ctx, &initial.pairs);
      const auto emission = profiles::unpaired_emission(deduction, is_unpaired);
      for (std::size_t r = 0; r < emission.count; ++r) {
        const auto &run = emission.runs[r];
        for (cand_pos_t pos = run.lo; pos <= run.hi; ++pos) {
          initial.emitted[static_cast<std::size_t>(pos)] = context_code(run.context);
        }
      }
      SummaryMap partial;
      merge(&partial, std::move(initial));
      for (std::size_t c = 0; c < deduction.child_count && !partial.empty(); ++c) {
        const auto found = chart.find(runtime_child_item(deduction, c));
        if (found == chart.end() || found->second.empty()) {
          partial.clear();
          break;
        }
        SummaryMap joined;
        for (const auto &[_, prefix] : partial) {
          for (const auto &[__, child] : found->second) {
            Summary combined = prefix;
            for (std::size_t pos = 1; pos <= static_cast<std::size_t>(n); ++pos) {
              if (child.emitted[pos] == kNone) continue;
              if (combined.emitted[pos] != kNone) {
                fail(label, "base " + std::to_string(pos) + " emitted twice");
                return false;
              }
              combined.emitted[pos] = child.emitted[pos];
            }
            for (const auto &pair : child.pairs) {
              engine::append_unique_traceback_pair(&combined.pairs, pair);
            }
            combined.weight *= child.weight;
            merge(&joined, std::move(combined));
          }
        }
        partial = std::move(joined);
      }
      if (!partial.empty() && emission.count > 0) {
        ++rule_coverage()[std::string(grammar::rule_name(deduction.rule))];
      }
      for (auto &[_, summary] : partial) merge(&groups, std::move(summary));
    }
    chart[item] = std::move(groups);
  }
  const ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
  const auto &derivations = chart[root];
  double z = 0.0;
  for (const auto &[_, summary] : derivations) z += summary.weight;
  if (!(z > 0.0)) return fail(label, "enumeration found no derivation");
  if (std::fabs(std::log(z) - profile.log_partition_weight) > 1e-9) {
    return fail(label, "enumerated log Z " + std::to_string(std::log(z)) +
                           " != exact log Z " +
                           std::to_string(profile.log_partition_weight));
  }

  std::vector<std::array<double, profiles::loop_context_count>> context(
      static_cast<std::size_t>(n) + 1);
  std::vector<std::vector<double>> access(windows.size(),
                                          std::vector<double>(n + 1, 0.0));
  using LoopRow = std::array<double, profiles::unpaired_loop_context_count>;
  std::vector<std::vector<LoopRow>> loop_access(windows.size(),
                                               std::vector<LoopRow>(n + 1));
  for (const auto &[_, summary] : derivations) {
    const double p = summary.weight / z;
    std::vector<bool> paired(static_cast<std::size_t>(n) + 1, false);
    for (cand_pos_t pos = 1; pos <= n; ++pos) {
      if (!ctx.is_unpaired(pos)) paired[pos] = true;  // paired in G
    }
    for (const auto &pair : summary.pairs) paired[pair.i] = paired[pair.j] = true;
    for (cand_pos_t pos = 1; pos <= n; ++pos) {
      const char code = summary.emitted[pos];
      if (paired[pos] != (code == kNone)) {
        std::string pairs;
        for (const auto &pair : summary.pairs) {
          pairs += std::to_string(pair.i) + "-" + std::to_string(pair.j) + " ";
        }
        return fail(label, "base " + std::to_string(pos) +
                               (paired[pos] ? " paired but emitted"
                                            : " unpaired but not emitted") +
                               " emitted=" + summary.emitted.substr(1) +
                               " pairs=" + pairs);
      }
      const std::string codes = "SHBIMEP";
      const std::size_t c = paired[pos] ? 0 : codes.find(code);
      context[pos][c] += p;
    }
    for (std::size_t w_idx = 0; w_idx < windows.size(); ++w_idx) {
      const cand_pos_t w = windows[w_idx];
      for (cand_pos_t s = 1; s + w - 1 <= n; ++s) {
        bool free = true;
        for (cand_pos_t t = s; t < s + w; ++t) free = free && !paired[t];
        if (free) {
          access[w_idx][s] += p;
          const char loop = summary.emitted[s];
          for (cand_pos_t t = s; t < s + w; ++t)
            if (summary.emitted[t] != loop)
              return fail(label, "unpaired window crosses loop labels");
          const auto c = std::string("HBIMEP").find(loop);
          if (c >= profiles::unpaired_loop_context_count)
            return fail(label, "unpaired window has no loop label");
          if (w > 1 && p > 0.0) window_loop_coverage()[c] = true;
          loop_access[w_idx][s][c] += p;
        }
      }
    }
  }

  constexpr double kTol = 1e-9;
  for (cand_pos_t pos = 1; pos <= n; ++pos) {
    for (std::size_t c = 0; c < profiles::loop_context_count; ++c) {
      const double got = profile.context[pos][c];
      const double want = context[pos][c];
      if (std::fabs(got - want) > kTol) {
        return fail(label, "context " +
                               std::string(profiles::loop_context_name(
                                   static_cast<LoopContext>(c))) +
                               " at " + std::to_string(pos) + ": got " +
                               std::to_string(got) + " want " +
                               std::to_string(want));
      }
    }
  }
  for (std::size_t w_idx = 0; w_idx < windows.size(); ++w_idx) {
    const cand_pos_t w = windows[w_idx];
    for (cand_pos_t s = 1; s + w - 1 <= n; ++s) {
      const double got = profile.accessibility[w_idx][s];
      const double want = access[w_idx][s];
      if (std::fabs(got - want) > kTol) {
        return fail(label, "accessibility w=" + std::to_string(w) + " s=" +
                               std::to_string(s) + ": got " +
                               std::to_string(got) + " want " +
                               std::to_string(want));
      }
      double sum = 0.0;
      for (std::size_t c = 0; c < profiles::unpaired_loop_context_count; ++c) {
        const double value = profile.accessibility_by_loop[w_idx][s][c];
        if (std::fabs(value - loop_access[w_idx][s][c]) > kTol)
          return fail(label, "loop accessibility differs from enumerated frequency");
        if (w == 1 && std::fabs(value - profile.context[s][c + 1]) > kTol)
          return fail(label, "length-one loop accessibility differs from CapR");
        sum += value;
      }
      if (std::fabs(sum - got) > kTol)
        return fail(label, "loop accessibility does not sum to total");
    }
  }
  std::size_t with_pseudoloop = 0;
  for (const auto &[_, summary] : derivations) {
    if (summary.emitted.find('P') != std::string::npos) ++with_pseudoloop;
  }
  std::cout << "[" << label << "] ok derivations=" << derivations.size()
            << " with_pseudoloop=" << with_pseudoloop << "\n";
  return true;
}

std::string random_scaffold(std::mt19937 &rng, std::string &sequence) {
  const int n = static_cast<int>(sequence.size());
  std::string g(n, '.');
  auto complement = [](char c) {
    switch (c) {
      case 'A': return 'U';
      case 'U': return 'A';
      case 'G': return 'C';
      default: return 'G';
    }
  };
  auto rand_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, std::max(lo, hi))(rng);
  };
  auto place_stem = [&](int a, int b, int len) {
    for (int t = 0; t < len; ++t) {
      if (b - a - 2 * t < 4 || g[a + t] != '.' || g[b - t] != '.') return t;
      bool inner_free = true;
      for (int x = a + t + 1; x < b - t; ++x) inner_free = inner_free && g[x] == '.';
      if (!inner_free) return t;
      g[a + t] = '(';
      g[b - t] = ')';
      sequence[b - t] = complement(sequence[a + t]);
    }
    return len;
  };
  int lo = 0, hi = n - 1;
  if (rand_int(0, 2) == 0) {
    const int len = place_stem(0, n - 1, rand_int(1, 2));
    lo = len + 1;
    hi = n - 2 - len;
  }
  const int hairpins = rand_int(1, 2);
  const int width = (hi - lo + 1) / hairpins;
  for (int h = 0; h < hairpins; ++h) {
    const int a0 = lo + h * width;
    const int b0 = a0 + width - 1;
    const int a = rand_int(a0, a0 + 2);
    const int b = rand_int(std::max(a + 5, b0 - 2), b0);
    if (b <= hi) place_stem(a, b, rand_int(1, 3));
  }
  return g;
}

}  // namespace

int main(int argc, char **argv) {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  bool ok = true;
  ok &= check_case("wmb-witness", "GCGCGCGCGCGCGC", "(.(......)...)");
  ok &= check_case("two-hairpins", "GGGAAACCCAGGGAAACCC", "(((...))).........."
  );
  ok &= check_case("empty-g", "GGGAAAUCCCAGG", ".............");
  ok &= check_case("gallery", "GGGGGAACCAACCC", "((.....)).....");
  ok &= check_case("fixed-bulge", "GAAGAAACC", "(..(...))");
  ok &= check_case("fixed-internal", "GAAGAAACAAC", "(..(...)..)");
  ok &= check_case("fixed-multi-exterior", "AAGAAGAAACAAGAAACAACAA",
                   "..(..(...)..(...)..)..");
  // Minimal three-base hairpin at the end of a BE wing's WIP start.
  // Enumerate the completed structures, pair emissions and unpaired runs.
  ok &= check_case("wip-minimal-hairpin", "GAGAAACAAGAAGAAAAAAAAACCAAAAC",
                   "(.(...)..(............)).....");
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_DirksPierce09.par", 1);
  ok &= check_case("wip-minimal-hairpin-dp09", "GAGAAACAAGAAGAAAAAAAAACCAAAAC",
                   "(.(...)..(............)).....");
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  const int count = argc > 1 ? std::atoi(argv[1]) : 12;
  std::mt19937 rng(20261007);
  std::uniform_int_distribution<int> base(0, 3);
  const int max_length = argc > 2 ? std::atoi(argv[2]) : 16;
  std::uniform_int_distribution<int> length(12, max_length);
  for (int c = 0; c < count; ++c) {
    std::string sequence(length(rng), 'A');
    for (auto &ch : sequence) ch = "ACGU"[base(rng)];
    const std::string scaffold = random_scaffold(rng, sequence);
    ok &= check_case("random" + std::to_string(c) + " " + sequence + " " + scaffold,
                     sequence, scaffold);
  }
  for (const auto &[rule, hits] : rule_coverage()) {
    std::cout << "coverage " << rule << " " << hits << "\n";
  }
  for (const bool covered : window_loop_coverage())
    if (!covered) ok = fail("loop-coverage", "missing a loop type for windows longer than one");
  if (!ok) return EXIT_FAILURE;
  std::cout << "scfg2_pk_profiles_oracle=ok\n";
  return EXIT_SUCCESS;
}

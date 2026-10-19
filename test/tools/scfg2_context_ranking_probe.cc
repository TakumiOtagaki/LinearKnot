#include "scfg2/weights/runtime_snapshot.hh"
// Offline ranking experiment on an immutable exact chart. This is not a
// production beam or a completion estimator available during forward pruning.
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/engine/outside.hh"
#include "scfg2/beam/wmb_btail_seeds.hh"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>

namespace sc = cparty::scfg2;
using Item = sc::deductions::ItemKey;
using Deduction = sc::deductions::Deduction;
using NT = sc::grammar::NonTerminal;
using LogSemiring = sc::semirings::LogSumExpSemiring;
using Clock = std::chrono::steady_clock;
constexpr std::size_t top_count = 16;
constexpr std::size_t mode_count = 10;
constexpr std::size_t oracle_mode = 8;
const double zero = LogSemiring::zero();
const std::array<const char *, mode_count> modes{
    "raw", "w_prefix", "zero_fallback", "max_context", "top4_contexts",
    "top16_contexts", "shortest_max", "shortest_sum", "oracle",
    "oracle_support_raw"};

double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

struct Guide {
  double sum = zero;
  double shortest_sum = zero;
  double shortest_max = zero;
  double log_paths = zero;
  int shortest = 1000000000;
  int longest = -1;
  std::array<double, top_count> top{};
  std::size_t size = 0;
  std::size_t incoming = 0;
  Guide() { top.fill(zero); }
};

using Guides = std::unordered_map<Item, Guide, sc::deductions::ItemKeyHash,
                                sc::deductions::ItemKeyEq>;

// Merge two sorted streams, preserving equal-weight distinct spine paths.
void add_top(const Guide &parent, double edge, Guide &child) {
  const auto old = child.top;
  const auto old_size = child.size;
  std::size_t a = 0, b = 0, size = 0;
  while (size < top_count && (a < old_size || b < parent.size)) {
    const double x = a < old_size ? old[a] : zero;
    const double y = b < parent.size ? parent.top[b] + edge : zero;
    if (a < old_size && (b == parent.size || x >= y)) {
      child.top[size++] = x;
      ++a;
    } else {
      child.top[size++] = y;
      ++b;
    }
  }
  child.size = size;
}

double top_sum(const Guide &guide, std::size_t k) {
  double value = zero;
  for (std::size_t i = 0; i < std::min(k, guide.size); ++i)
    value = LogSemiring::plus(value, guide.top[i]);
  return value;
}

bool key_less(const Item &a, const Item &b) {
  return std::tie(a.i, a.j, a.ip, a.jp) <
         std::tie(b.i, b.j, b.ip, b.jp);
}

struct Row {
  Item item;
  std::string family;
  double alpha = zero;
  double posterior = 0;
  bool prefix_zero = false;
  int shortest = -1, longest = -1;
  double log_paths = zero;
  std::array<double, mode_count> score{};
};

std::array<double, mode_count> scores(double alpha, double prefix,
                                     const Guide &guide) {
  return {alpha, alpha + prefix, alpha + (prefix == zero ? 0 : prefix),
          alpha + top_sum(guide, 1), alpha + top_sum(guide, 4),
          alpha + top_sum(guide, 16), alpha + guide.shortest_max,
          alpha + guide.shortest_sum, alpha + guide.sum,
          guide.sum == zero ? zero : alpha};
}

// A derived factor has no public item identity. These sums observe complete
// WMB terms, grouped by the actual BTAIL key (i,t); they do not add a new rule.
struct Tail {
  double alpha = zero;
  double full_mass = zero;
  std::map<std::tuple<int, int, int>, double> contexts;
};

std::ofstream output(const std::string &path) {
  std::ofstream stream(path);
  if (!stream) throw std::runtime_error("cannot write " + path);
  stream << std::setprecision(17);
  return stream;
}

void report_rows(const std::vector<Row> &rows, const std::vector<int> &widths,
                 const std::string &prefix) {
  std::map<std::pair<std::string, int>, std::vector<std::size_t>> buckets;
  for (std::size_t r = 0; r < rows.size(); ++r)
    buckets[{rows[r].family, rows[r].item.j}].push_back(r);
  auto table = output(prefix + ".ranking.tsv");
  table << "family\twidth\tmode\tcompetitive_buckets\tcandidates"
           "\texpected_occurrences\tretained_occurrences\toracle_retained"
           "\tretained_fraction\tzero_prefix_occurrences\n";
  auto items = output(prefix + ".items.tsv");
  items << "family\ti\tj\tip\tjp\talpha\tposterior\tprefix_zero"
           "\tshortest_hops\tlongest_hops\tlog10_spine_paths";
  for (const auto *mode : modes) items << '\t' << mode << "_rank";
  items << '\n';
  std::vector<std::array<int, mode_count>> ranks(rows.size());
  for (const auto &[bucket, indices] : buckets) {
    for (std::size_t mode = 0; mode < mode_count; ++mode) {
      auto sorted = indices;
      std::sort(sorted.begin(), sorted.end(), [&](auto a, auto b) {
        if (rows[a].score[mode] != rows[b].score[mode])
          return rows[a].score[mode] > rows[b].score[mode];
        return key_less(rows[a].item, rows[b].item);
      });
      for (std::size_t r = 0; r < sorted.size(); ++r)
        ranks[sorted[r]][mode] = static_cast<int>(r + 1);
    }
  }
  for (const int width : widths) {
    struct Stats {
      std::size_t buckets = 0, candidates = 0;
      double total = 0, oracle = 0, zero_prefix = 0;
      std::array<double, mode_count> retained{};
    };
    std::map<std::string, Stats> stats;
    for (const auto &[bucket, indices] : buckets) {
      if (indices.size() <= static_cast<std::size_t>(width)) continue;
      double total = 0;
      for (auto r : indices) total += rows[r].posterior;
      if (total == 0) continue;
      auto &st = stats[bucket.first];
      ++st.buckets;
      st.candidates += indices.size();
      st.total += total;
      for (auto r : indices) {
        const auto q = rows[r].posterior;
        if (rows[r].prefix_zero) st.zero_prefix += q;
        for (std::size_t mode = 0; mode < mode_count; ++mode)
          if (ranks[r][mode] <= width) st.retained[mode] += q;
      }
      st.oracle = st.retained[oracle_mode];
    }
    for (const auto &[family, st] : stats)
      for (std::size_t mode = 0; mode < mode_count; ++mode) {
        if (family == "BTAIL" && (mode == 6 || mode == 7)) continue;
        table << family << '\t' << width << '\t' << modes[mode] << '\t'
              << st.buckets << '\t' << st.candidates << '\t' << st.total
              << '\t' << st.retained[mode] << '\t' << st.oracle << '\t'
              << st.retained[mode] / st.total << '\t' << st.zero_prefix << '\n';
      }
  }
  for (std::size_t r = 0; r < rows.size(); ++r) {
    const auto &row = rows[r];
    if (row.posterior < 1e-8) continue;
    items << row.family << '\t' << row.item.i << '\t' << row.item.j << '\t'
          << row.item.ip << '\t' << row.item.jp << '\t' << row.alpha << '\t'
          << row.posterior << '\t' << row.prefix_zero << '\t' << row.shortest
          << '\t' << row.longest << '\t' << row.log_paths / std::log(10.0);
    for (std::size_t mode = 0; mode < mode_count; ++mode)
      items << '\t' << (row.family == "BTAIL" && (mode == 6 || mode == 7)
                            ? -1 : ranks[r][mode]);
    items << '\n';
  }
}

int main(int argc, char **argv) {
  try {
    if (argc < 5 || argc > 6 ||
        (argc == 6 && std::string(argv[5]) != "time-scalar"))
      throw std::runtime_error(
          "usage: probe INPUT given|empty OUTPUT_PREFIX WIDTHS_CSV [time-scalar]");
    std::ifstream input(argv[1]);
    std::string line, sequence, scaffold;
    while (std::getline(input, line)) {
      if (line.empty() || line.front() == '>') continue;
      if (sequence.empty()) sequence = line;
      else { scaffold = line; break; }
    }
    if (std::string(argv[2]) == "empty") scaffold.assign(sequence.size(), '.');
    else if (std::string(argv[2]) != "given")
      throw std::runtime_error("mode must be given or empty");
    if (sequence.empty() || scaffold.size() != sequence.size() ||
        sequence.size() > 512)
      throw std::runtime_error("expected a sequence/scaffold of 1..512 nt");
    std::vector<int> widths;
    std::stringstream width_stream(argv[4]);
    while (std::getline(width_stream, line, ',')) {
      const int width = std::stoi(line);
      if (width < 1) throw std::runtime_error("width must be positive");
      widths.push_back(width);
    }
    if (widths.empty()) throw std::runtime_error("no widths");
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    const auto total_start = Clock::now();
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    sc::replay::WFinalLogWeightContext ctx(owner);
    const auto n = static_cast<cand_pos_t>(sequence.size());
    const auto inside_start = Clock::now();
    const auto schedule = sc::engine::build_pkfree_basic_exact_schedule(n, ctx);
    const auto inside = sc::replay::run_w_final_exact_dp_schedule<LogSemiring>(
        schedule, ctx, [](const auto &) { return true; },
        [](const auto &, const auto &, const auto &, const auto &, auto) {});
    const double inside_seconds = seconds(inside_start);
    const Item root{NT::W, 1, n, -1, -1};
    const double log_z = inside.get(root);
    if (!std::isfinite(log_z)) throw std::runtime_error("root is unreachable");
    std::cerr << "inside n=" << n << " log_z=" << log_z
              << " seconds=" << inside_seconds << '\n';
    sc::engine::AllLocalWeightModel<decltype(ctx)> weights(&ctx);
    sc::engine::RuntimeDeductionChildItem child_item;
    Guides guides;
    auto &initial = guides[root];
    initial.sum = initial.shortest_sum = initial.shortest_max = 0;
    initial.log_paths = 0;
    initial.shortest = initial.longest = 0;
    initial.top[0] = 0;
    initial.size = 1;
    std::unordered_set<Item, sc::deductions::ItemKeyHash,
                       sc::deductions::ItemKeyEq> processed;
    std::map<std::pair<int, int>, Tail> tails;
    std::size_t generated = 0, live_terms = 0, edges = 0;
    std::size_t bound_checks = 0;
    double generator_seconds = 0;
    const auto outside_start = Clock::now();
    for (auto cursor = schedule.rbegin(); cursor != schedule.rend(); ++cursor) {
      const auto found = guides.find(*cursor);
      if (found == guides.end()) continue;
      const Guide parent = found->second;
      if (!processed.insert(*cursor).second)
        throw std::runtime_error("duplicate reverse item");
      const auto generation_start = Clock::now();
      const auto deductions =
          sc::deductions::generate_exact_basic_deductions_for_item(*cursor, ctx);
      generator_seconds += seconds(generation_start);
      generated += deductions.size();
      for (const auto &ded : deductions) {
        if (!sc::engine::legacy_span_schedule_supports(ded)) continue;
        const double local = weights.log_weight_for(ded);
        if (local == zero) continue;
        std::array<Item, 4> children{};
        std::array<double, 4> alpha{};
        double full = local;
        bool live = true;
        if (ded.child_count > children.size())
          throw std::runtime_error("unexpected rule arity");
        for (std::size_t c = 0; c < ded.child_count; ++c) {
          children[c] = child_item(ded, c);
          alpha[c] = inside.get(children[c]);
          if (alpha[c] == zero) { live = false; break; }
          full += alpha[c];
        }
        if (!live) continue;
        ++live_terms;
        if (ded.rule == sc::grammar::RuleId::WMB_SPLIT_BE_WMBP_WI) {
          Item be;
          for (std::size_t c = 0; c < ded.child_count; ++c)
            if (children[c].nonterminal == NT::BE) be = children[c];
          if (be.nonterminal != NT::BE) throw std::runtime_error("missing BE");
          auto &tail = tails[{ded.parent.i, be.jp}];
          tail.full_mass = LogSemiring::plus(tail.full_mass, parent.sum + full);
          const auto context = std::make_tuple(be.i, be.j, be.ip);
          const double value = parent.sum + local + inside.get(be);
          const auto previous = tail.contexts.emplace(context, value);
          if (!previous.second && std::fabs(previous.first->second - value) > 1e-9)
            throw std::runtime_error("BTAIL local factor depends on split");
        }
        for (std::size_t c = 0; c < ded.child_count; ++c) {
          if (processed.count(children[c]))
            throw std::runtime_error("reverse schedule is not topological");
          double edge = local;
          for (std::size_t d = 0; d < ded.child_count; ++d)
            if (d != c) edge += alpha[d];
          auto &child = guides[children[c]];
          ++child.incoming;
          ++edges;
          child.sum = LogSemiring::plus(child.sum, parent.sum + edge);
          child.log_paths = LogSemiring::plus(child.log_paths, parent.log_paths);
          child.longest = std::max(child.longest, parent.longest + 1);
          const int hops = parent.shortest + 1;
          if (hops < child.shortest) {
            child.shortest = hops;
            child.shortest_sum = parent.shortest_sum + edge;
            child.shortest_max = parent.shortest_max + edge;
          } else if (hops == child.shortest) {
            child.shortest_sum = LogSemiring::plus(child.shortest_sum,
                                        parent.shortest_sum + edge);
            child.shortest_max = std::max(child.shortest_max,
                                         parent.shortest_max + edge);
          }
          add_top(parent, edge, child);
        }
      }
    }
    const double outside_seconds = seconds(outside_start);
    std::size_t max_incoming = 0;
    int max_shortest = 0, max_longest = 0;
    double max_log_paths = 0, max_bound_error = 0;
    for (const auto &[item, guide] : guides) {
      max_incoming = std::max(max_incoming, guide.incoming);
      max_shortest = std::max(max_shortest, guide.shortest);
      max_longest = std::max(max_longest, guide.longest);
      max_log_paths = std::max(max_log_paths, guide.log_paths / std::log(10.0));
      const std::array<double, 4> lower{
          top_sum(guide, 1), top_sum(guide, 4), top_sum(guide, 16),
          guide.shortest_sum};
      for (double value : lower) {
        const double error = value - guide.sum;
        max_bound_error = std::max(max_bound_error, error);
        if (error > 1e-9) throw std::runtime_error("guide exceeds full outside");
        ++bound_checks;
      }
      if (lower[0] > lower[1] + 1e-9 || lower[1] > lower[2] + 1e-9 ||
          guide.shortest_max > guide.shortest_sum + 1e-9)
        throw std::runtime_error("truncated guide is not monotone");
      if (guide.log_paths <= std::log(double(top_count)) + 1e-10 &&
          std::fabs(lower[2] - guide.sum) > 1e-9)
        throw std::runtime_error("top-16 lost one of at most sixteen paths");
    }
    double reference_error = 0;
    double scalar_sum_seconds = -1, scalar_max_seconds = -1;
    if (argc == 6) {
      // Same generated positive-support graph for both scalar operations.
      // Exact inside is supplied unchanged, even to the max-spine pass.
      const auto generate_live = [&](const Item &item) {
        auto ds = sc::deductions::generate_exact_basic_deductions_for_item(item, ctx);
        ds.erase(std::remove_if(ds.begin(), ds.end(), [&](const auto &d) {
          if (!sc::engine::legacy_span_schedule_supports(d)) return true;
          for (std::size_t c = 0; c < d.child_count; ++c)
            if (inside.get(child_item(d, c)) == zero) return true;
          return false;
        }), ds.end());
        return ds;
      };
      const auto scalar_start = Clock::now();
      {
        const auto beta = sc::engine::run_regenerated_outside<LogSemiring>(
            schedule, inside, root, generate_live,
            [&](const auto &d) { return weights.local_weight_for(d); }, child_item);
        scalar_sum_seconds = seconds(scalar_start);
        for (const auto &[item, guide] : guides)
          if (std::fabs(beta.get(item) - guide.sum) > 1e-9)
            throw std::runtime_error("scalar sum differs from experiment");
      }
      const auto max_start = Clock::now();
      {
        const auto beta = sc::engine::run_regenerated_outside<
            sc::semirings::LogMaxSemiring>(
            schedule, inside, root, generate_live,
            [&](const auto &d) { return weights.local_weight_for(d); }, child_item);
        scalar_max_seconds = seconds(max_start);
        for (const auto &[item, guide] : guides)
          if (std::fabs(beta.get(item) - top_sum(guide, 1)) > 1e-9)
            throw std::runtime_error("scalar max differs from experiment");
      }
    }
    if (n <= 32) {
      const auto reference = sc::engine::run_regenerated_outside<LogSemiring>(
          schedule, inside, root,
          [&](const Item &item) {
            auto ds = sc::deductions::generate_exact_basic_deductions_for_item(item, ctx);
            ds.erase(std::remove_if(ds.begin(), ds.end(), [](const auto &d) {
              return !sc::engine::legacy_span_schedule_supports(d);
            }), ds.end());
            return ds;
          }, [&](const auto &d) { return weights.local_weight_for(d); }, child_item);
      for (const auto &[item, value] : inside.values()) {
        if (value == zero) continue;
        const auto observed = guides.find(item);
        const double beta = observed == guides.end() ? zero : observed->second.sum;
        const double expected = reference.get(item);
        if (beta == zero && expected == zero) continue;
        if (!std::isfinite(beta) || !std::isfinite(expected))
          throw std::runtime_error("outside support differs from shared engine");
        reference_error = std::max(reference_error, std::fabs(beta - expected));
      }
      if (reference_error > 1e-9)
        throw std::runtime_error("outside differs from shared engine");
    }
    std::vector<Row> rows;
    for (const auto &[item, value] : inside.values()) {
      if (value == zero) continue;
      const auto found = guides.find(item);
      const Guide guide = found == guides.end() ? Guide{} : found->second;
      const double prefix = item.i <= 1 ? 0 : inside.get({NT::W, 1, item.i - 1, -1, -1});
      Row row;
      row.item = item;
      row.family = sc::grammar::nonterminal_name(item.nonterminal);
      row.alpha = value;
      row.posterior = std::exp(value + guide.sum - log_z);
      row.prefix_zero = prefix == zero;
      if (found != guides.end()) {
        row.shortest = guide.shortest;
        row.longest = guide.longest;
      }
      row.log_paths = guide.log_paths;
      row.score = scores(value, prefix, guide);
      rows.push_back(row);
    }
    const auto public_positive = rows.size();
    // Observe the production factorization using the shared BTAIL seed helper.
    std::vector<std::vector<int>> wi_starts(n + 1), wmbp_starts(n + 1);
    for (const auto &[item, value] : inside.values()) {
      if (value == zero) continue;
      if (item.nonterminal == NT::WI) wi_starts[item.j].push_back(item.i);
      if (item.nonterminal == NT::WMBP) wmbp_starts[item.j].push_back(item.i);
    }
    for (int t = 2; t <= n; ++t)
      sc::beam::for_each_wmb_btail_seed(t, ctx, wi_starts[t - 1],
          [&](int k) -> const auto & { return wmbp_starts[k]; },
          [&](int i, int k, bool has_wi) {
            const double value = inside.get({NT::WMBP, i, k, -1, -1}) +
                (has_wi ? inside.get({NT::WI, k + 1, t - 1, -1, -1}) : 0);
            auto &tail = tails[{i, t}];
            tail.alpha = LogSemiring::plus(tail.alpha, value);
          });
    double tail_factor_error = 0;
    for (const auto &[key, tail] : tails) {
      if (tail.alpha == zero) {
        if (tail.full_mass != zero) throw std::runtime_error("missing BTAIL seeds");
        continue;
      }
      double sum = zero, best = zero;
      std::vector<double> contexts;
      for (const auto &[context, value] : tail.contexts) {
        sum = LogSemiring::plus(sum, value);
        best = std::max(best, value);
        contexts.push_back(value);
      }
      if (sum != zero) {
        tail_factor_error = std::max(tail_factor_error,
            std::fabs(tail.alpha + sum - tail.full_mass));
        if (tail_factor_error > 1e-8)
          throw std::runtime_error("BTAIL factorization differs from complete terms");
      }
      std::sort(contexts.begin(), contexts.end(), std::greater<double>());
      Guide completion;
      completion.sum = sum;
      completion.shortest_sum = sum;  // Not a hop experiment for BTAIL.
      completion.shortest_max = best;
      completion.size = std::min(top_count, contexts.size());
      for (std::size_t c = 0; c < completion.size; ++c) completion.top[c] = contexts[c];
      const double prefix = key.first <= 1 ? 0 :
          inside.get({NT::W, 1, key.first - 1, -1, -1});
      Row row;
      row.item = {NT::WMBP, key.first, key.second, -1, -1};
      row.family = "BTAIL";
      row.alpha = tail.alpha;
      row.posterior = std::exp(tail.full_mass - log_z);
      row.prefix_zero = prefix == zero;
      row.score = scores(tail.alpha, prefix, completion);
      // BTAIL top-K means BE completion alternatives with exact WMB outside,
      // not top-K root spines. Keep unsupported shortest columns impossible.
      row.score[6] = row.score[7] = zero;
      rows.push_back(row);
    }
    std::sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) {
      if (a.family != b.family) return a.family < b.family;
      return key_less(a.item, b.item);
    });
    report_rows(rows, widths, argv[3]);
    auto summary = output(std::string(argv[3]) + ".summary.tsv");
    summary << "n\tlog_z\tscheduled_items\tpositive_items\treachable_items"
               "\tgenerated_deductions\tlive_terms\toutside_edges"
               "\tmax_incoming\tmax_shortest_hops\tmax_longest_hops"
               "\tmax_log10_spine_paths\tinside_seconds\tguide_seconds"
               "\tgeneration_seconds\ttotal_seconds\tbound_checks"
               "\tmax_bound_error\treference_error\ttail_factor_error"
               "\tscalar_sum_seconds\tscalar_max_seconds\n";
    summary << n << '\t' << log_z << '\t' << schedule.size() << '\t'
            << public_positive << '\t' << guides.size() << '\t'
            << generated << '\t' << live_terms << '\t' << edges << '\t'
            << max_incoming << '\t' << max_shortest << '\t' << max_longest
            << '\t' << max_log_paths << '\t' << inside_seconds << '\t'
            << outside_seconds << '\t' << generator_seconds << '\t'
            << seconds(total_start) << '\t' << bound_checks << '\t'
            << max_bound_error << '\t' << reference_error << '\t'
            << tail_factor_error << '\t' << scalar_sum_seconds << '\t'
            << scalar_max_seconds << '\n';
    std::cerr << "guide vertices=" << guides.size() << " edges=" << edges
              << " seconds=" << outside_seconds << " checks=" << bound_checks << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}

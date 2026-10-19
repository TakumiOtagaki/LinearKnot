#include "scfg2/weights/runtime_snapshot.hh"
// Diagnostic only: exact alpha/beta evaluate ranking, while each experimental
// outer uses ONLY a finite-beam pilot's sibling alphas. Full exact candidate
// enumeration remains here deliberately; this is not a linear-time runtime.
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/engine/outside.hh"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include "scfg2_outer_probe_common.hh"

std::vector<int> integers(const std::string &csv, bool allow_zero = false) {
  std::vector<int> out;
  std::stringstream stream(csv);
  std::string token;
  while (std::getline(stream, token, ',')) {
    const int value = std::stoi(token);
    if (value < (allow_zero ? 0 : 1)) throw std::runtime_error("invalid width/cap");
    out.push_back(value);
  }
  if (out.empty()) throw std::runtime_error("empty width/cap list");
  return out;
}
std::ofstream output(const std::string &path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << std::setprecision(17);
  return out;
}
void check_le(double a, double b, const char *message) {
  if (a == zero) return;
  if (!std::isfinite(a) || b == zero || a > b + 1e-8)
    throw std::runtime_error(message);
}

// Observation only: one raw-alpha reserve per (family,j), before public
// pruning. The final pilot chart is always included as a context backbone.
struct CandidateReserve {
  struct Candidate { Item key; double alpha; };
  struct Better {
    bool operator()(const Candidate &a, const Candidate &b) const {
      if (a.alpha != b.alpha) return a.alpha > b.alpha;
      return less_key(a.key, b.key);
    }
  };
  static constexpr std::size_t families = static_cast<std::size_t>(NT::BE) + 1;
  std::size_t limit;
  std::vector<std::array<std::set<Candidate, Better>, families>> buckets;
  CandidateReserve(int n, std::size_t limit) : limit(limit), buckets(n + 1) {}
  template <class... Args> void operator()(const Args &...) {}
  void after_item(const Item &key, const Chart &chart) {
    const double alpha = chart.get(key);
    if (alpha == zero) return;
    auto &bucket = buckets.at(key.j).at(static_cast<std::size_t>(key.nonterminal));
    bucket.insert({key, alpha});
    if (bucket.size() > limit) bucket.erase(std::prev(bucket.end()));
  }
  ItemSet allowed(const Chart &pilot, std::size_t width) const {
    ItemSet out;
    for (const auto &[key, alpha] : pilot.values())
      if (alpha != zero) out.insert(key);
    for (const auto &column : buckets)
      for (const auto &bucket : column) {
        std::size_t count = 0;
        for (const auto &candidate : bucket) {
          if (count++ == width) break;
          out.insert(candidate.key);
        }
      }
    return out;
  }
};

// Target alpha is not read to score its context. For A <- B C, the A -> B
// edge is local(A <- B C) * pilot_alpha(C), even if B was pruned by the pilot.
// A cap keeps the best distinct targets PER PARENT AND TARGET FAMILY by this
// edge weight, before any child outer is known. Equal targets collapse by max.
// cap=0 is unlimited. No exact alpha/beta enters this function.
template <class Context>
Outer pilot_outer(const std::vector<Item> &schedule, const Chart &pilot,
                  const Item &root, Context &ctx, int cap, bool strict,
                  bool preserve_pilot = false, const ItemSet *allowed = nullptr) {
  const auto start = Clock::now();
  Outer result;
  result.values.set(root, 0);
  sc::engine::AllLocalWeightModel<Context> weights(&ctx);
  sc::engine::RuntimeDeductionChildItem child_item;
  std::unordered_set<Item, sc::deductions::ItemKeyHash,
                     sc::deductions::ItemKeyEq> processed;
  for (auto cursor = schedule.rbegin(); cursor != schedule.rend(); ++cursor) {
    processed.insert(*cursor);
    const double parent_outer = result.values.get(*cursor);
    if (parent_outer == zero) continue;
    ++result.visited;
    ItemMap edges;
    const auto deductions =
        sc::deductions::generate_exact_basic_deductions_for_item(*cursor, ctx);
    result.generated += deductions.size();
    for (const auto &d : deductions) {
      ++result.generated_rules.at(static_cast<std::size_t>(d.rule));
      if (!sc::engine::legacy_span_schedule_supports(d)) continue;
      const double local = weights.log_weight_for(d);
      if (local == zero) continue;
      for (std::size_t c = 0; c < d.child_count; ++c) {
        const Item target = child_item(d, c);
        if (allowed && !allowed->count(target)) continue;
        if (strict && pilot.get(target) == zero) continue;
        double edge = local;
        for (std::size_t k = 0; k < d.child_count; ++k)
          if (c != k) edge += pilot.get(child_item(d, k));
        if (edge == zero) continue;
        if (!std::isfinite(edge)) throw std::runtime_error("invalid edge");
        if (processed.count(target)) throw std::runtime_error("non-topological edge");
        ++result.eligible;
        ++result.eligible_rules.at(static_cast<std::size_t>(d.rule));
        const auto found = edges.find(target);
        if (found == edges.end()) edges.emplace(target, edge);
        else found->second = std::max(found->second, edge);
      }
    }
    result.max_parent_candidates = std::max(result.max_parent_candidates, edges.size());
    std::map<NT, std::vector<std::pair<Item, double>>> groups;
    for (const auto &edge : edges) groups[edge.first.nonterminal].push_back(edge);
    for (auto &[family, candidates] : groups) {
      (void)family;
      std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b) {
        if (a.second != b.second) return a.second > b.second;
        return less_key(a.first, b.first);
      });
      std::size_t novel = 0;
      for (const auto &[child, edge] : candidates)
        if (!cap || (preserve_pilot && pilot.get(child) != zero) || novel++ < static_cast<std::size_t>(cap)) {
          ++result.selected;
          update_outer(result, child, parent_outer + edge, {*cursor, edge});
        }
    }
  }
  if (result.visited != result.values.values().size())
    throw std::runtime_error("outer target missing from schedule");
  result.elapsed = seconds(start);
  return result;
}

// Reverse family/span buckets are topological, unlike a naive right-end
// family quota (W can read same-family, same-right-end suffix W). Derive only
// bucket order from the shared exact schedule; no exact values enter here.
template <class Context>
Outer span_outer(const std::vector<Item> &schedule, const Chart &pilot,
                 const Item &root, Context &ctx, std::size_t width) {
  const auto start = Clock::now();
  constexpr auto families = CandidateReserve::families;
  std::vector<std::array<int, families>> group_at(root.j + 1);
  for (auto &row : group_at) row.fill(-1);
  std::vector<MaxPool> pools;
  for (const auto &item : schedule) {
    const auto span = item.j - item.i + 1;
    auto &group = group_at.at(span).at(static_cast<std::size_t>(item.nonterminal));
    if (group < 0) { group = pools.size(); pools.emplace_back(); }
  }
  const auto group_for = [&](const Item &key) {
    const auto span = key.j - key.i + 1;
    if (span < 1 || span > root.j) throw std::runtime_error("invalid span target");
    const int group = group_at.at(span).at(static_cast<std::size_t>(key.nonterminal));
    if (group < 0) throw std::runtime_error("unscheduled span family");
    return group;
  };
  Outer result;
  pools.at(group_for(root)).offer(root, 0, {root, 0}, true, width);
  std::size_t pending = 1;
  result.peak_pending = result.peak_guide_states = 1;
  sc::engine::AllLocalWeightModel<Context> weights(&ctx);
  sc::engine::RuntimeDeductionChildItem child_item;
  for (int g = static_cast<int>(pools.size()) - 1; g >= 0; --g) {
    auto parents = std::move(pools[g].entries);
    pools[g].novel.clear();
    pending -= parents.size();
    for (const auto &[parent, entry] : parents) {
      result.values.set(parent, entry.value);
      if (!sc::deductions::ItemKeyEq{}(parent, root)) result.links.emplace(parent, entry.link);
    }
    for (const auto &[parent, entry] : parents) {
      ++result.visited;
      const auto deductions = sc::deductions::generate_exact_basic_deductions_for_item(parent, ctx);
      result.generated += deductions.size();
      std::size_t parent_candidates = 0;
      for (const auto &d : deductions) {
        ++result.generated_rules.at(static_cast<std::size_t>(d.rule));
        if (!sc::engine::legacy_span_schedule_supports(d)) continue;
        const double local = weights.log_weight_for(d);
        if (local == zero) continue;
        for (std::size_t c = 0; c < d.child_count; ++c) {
          double edge = local;
          for (std::size_t k = 0; k < d.child_count; ++k)
            if (c != k) edge += pilot.get(child_item(d, k));
          if (edge == zero) continue;
          if (!std::isfinite(edge)) throw std::runtime_error("invalid span edge");
          const Item child = child_item(d, c);
          const int target = group_for(child);
          if (target >= g) throw std::runtime_error("span bucket dependency violation");
          ++parent_candidates; ++result.eligible;
          ++result.eligible_rules.at(static_cast<std::size_t>(d.rule));
          ++result.selected;
          auto &pool = pools[target];
          const auto previous_size = pool.entries.size();
          pool.offer(child, entry.value + edge, {parent, edge}, pilot.get(child) != zero, width);
          pending -= previous_size;
          pending += pool.entries.size();
          result.peak_pending = std::max(result.peak_pending, pending);
          result.peak_guide_states = std::max(result.peak_guide_states, result.values.values().size() + pending);
        }
      }
      result.max_parent_candidates = std::max(result.max_parent_candidates, parent_candidates);
    }
  }
  for (const auto &pool : pools) result.evicted += pool.evicted;
  if (pending) throw std::runtime_error("unprocessed span candidates");
  const auto bound = pilot.values().size() + families * std::size_t(root.j) * width;
  if (result.peak_guide_states > bound) throw std::runtime_error("span state budget violated");
  for (const auto &[child, link] : result.links) {
    if (group_for(link.parent) <= group_for(child) || result.values.get(link.parent) == zero)
      throw std::runtime_error("disconnected span path");
    if (std::fabs(result.values.get(link.parent) + link.edge - result.values.get(child)) > 1e-9)
      throw std::runtime_error("invalid span path weight");
  }
  if (result.links.size() + 1 != result.values.values().size()) throw std::runtime_error("missing span predecessor");
  result.elapsed = seconds(start);
  return result;
}

#include "scfg2_bounded_outer.hh"

void bounded_header(std::ostream &out) {
  out << "n\tpilot_width\tbounded_factor\tinput_records\tguide_records\tindex_records"
         "\tbuild_comparisons\tbuild_moves\tundo_moves\tcursor_steps\tquery_outputs"
         "\tproposals\tmax_proposals\tvalidation_checks\tvalidation_seconds"
         "\touter_items\tgenerated_deductions\teligible_edges\tpeak_guide_states\touter_seconds\n";
}
void bounded_row(std::ostream &out, int n, int p, int factor,
                 const bounded_outer_detail::Result &result) {
  const auto &w = result.work;
  const auto &g = result.outer;
  out << n << '\t' << p << '\t' << factor << '\t' << w.input_records << '\t'
      << w.guide_records << '\t' << w.index_records << '\t' << w.comparisons << '\t'
      << w.build_moves << '\t' << w.undo_moves << '\t' << w.cursor_steps << '\t'
      << w.query_outputs << '\t' << w.proposals << '\t' << w.max_proposals << '\t'
      << w.checks << '\t' << w.validation_seconds << '\t' << g.values.values().size()
      << '\t' << g.generated << '\t' << g.eligible << '\t' << g.peak_guide_states
      << '\t' << g.elapsed << '\n';
}

struct Row {
  Item key;
  double alpha, q, prefix;
};

// Causal control only: uses the already generated UNBOUNDED pilot-outer tree.
// Closing a reserve under these predecessors retains its best paths exactly,
// but does not provide a cheap or proven linear construction of that reserve.
ItemSet close_best_paths(ItemSet allowed, const Outer &all, const Item &root) {
  ItemSet complete;
  complete.insert(root);
  std::vector<Item> seeds(allowed.begin(), allowed.end());
  for (auto cursor : seeds) {
    if (all.values.get(cursor) == zero) continue;
    while (complete.insert(cursor).second) {
      allowed.insert(cursor);
      const auto found = all.links.find(cursor);
      if (found == all.links.end()) throw std::runtime_error("missing best-path parent");
      cursor = found->second.parent;
    }
  }
  allowed.insert(root);
  return allowed;
}

void report_bridge_paths(std::ofstream &out, const std::vector<Row> &rows,
                          const ItemSet &allowed, const Outer &restricted,
                          const Outer &all, int pilot_width, int factor,
                          const Chart &pilot) {
  for (const bool discovery : {false, true}) {
    std::vector<const Row *> lost;
    for (const auto &row : rows) {
      if (row.q == 0 || all.values.get(row.key) == zero ||
          bool(allowed.count(row.key)) == discovery) continue;
      if (all.values.get(row.key) > restricted.values.get(row.key) + 1e-8)
        lost.push_back(&row);
    }
    std::sort(lost.begin(), lost.end(), [](const auto *a, const auto *b) {
      if (a->q != b->q) return a->q > b->q;
      if (a->key.nonterminal != b->key.nonterminal) return a->key.nonterminal < b->key.nonterminal;
      return less_key(a->key, b->key);
    });
    if (lost.size() > 8) lost.resize(8);
    for (const auto *row : lost) {
      std::vector<Item> path;
      Item cursor = row->key;
      for (;;) {
        path.push_back(cursor);
        const auto link = all.links.find(cursor);
        if (link == all.links.end()) break;
        cursor = link->second.parent;
        if (path.size() > all.values.values().size()) throw std::runtime_error("cyclic best path");
      }
      std::reverse(path.begin(), path.end());
      for (std::size_t step = 0; step < path.size(); ++step) {
        const auto &item = path[step];
        out << pilot_width << '\t' << factor << '\t' << (discovery ? "undiscovered_target" : "missing_bridge")
            << '\t' << sc::grammar::nonterminal_name(row->key.nonterminal) << '\t' << row->key.i << '\t' << row->key.j
            << '\t' << row->q << '\t' << step << '\t' << sc::grammar::nonterminal_name(item.nonterminal)
            << '\t' << item.i << '\t' << item.j << '\t' << item.ip << '\t' << item.jp
            << '\t' << (pilot.get(item) != zero) << '\t' << bool(allowed.count(item))
            << '\t' << restricted.values.get(item) << '\t' << all.values.get(item) << '\n';
      }
    }
  }
}

int main(int argc, char **argv) {
  try {
    if (argc < 7 || argc > 11) throw std::runtime_error(
        "usage: probe INPUT given|empty OUTPUT_PREFIX RANK_WIDTHS PILOT_WIDTHS EDGE_CAPS(0=all) [RESERVE_FACTORS [SPAN_FACTORS [BOUNDED_FACTORS [bounded-only]]]]");
    check_max_pool();
    bounded_outer_detail::check_prefix_top();
    std::ifstream input(argv[1]);
    std::string line, sequence, scaffold;
    while (std::getline(input, line)) {
      if (line.empty() || line[0] == '>') continue;
      if (sequence.empty()) sequence = line;
      else { scaffold = line; break; }
    }
    if (std::string(argv[2]) == "empty") scaffold.assign(sequence.size(), '.');
    else if (std::string(argv[2]) != "given") throw std::runtime_error("invalid scaffold mode");
    if (sequence.empty() || sequence.size() > 512 || scaffold.size() != sequence.size())
      throw std::runtime_error("expected sequence/scaffold of 1..512 nt");
    const auto rank_widths = integers(argv[4]);
    const auto pilot_widths = integers(argv[5]);
    const auto reserve_factors = argc >= 8 && std::string(argv[7]).size() ? integers(argv[7]) : std::vector<int>{};
    const auto span_factors = argc >= 9 && std::string(argv[8]).size() ? integers(argv[8]) : std::vector<int>{};
    const auto bounded_factors = argc >= 10 && std::string(argv[9]).size() ? integers(argv[9]) : std::vector<int>{};
    const bool bounded_only = argc == 11;
    if (bounded_only && (std::string(argv[10]) != "bounded-only" || bounded_factors.empty()))
      throw std::runtime_error("bounded-only requires bounded factors");
    auto caps = integers(argv[6], true);
    if (std::find(caps.begin(), caps.end(), 0) == caps.end()) caps.push_back(0);
    std::sort(caps.begin(), caps.end());
    caps.erase(std::unique(caps.begin(), caps.end()), caps.end());
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    const auto total_start = Clock::now();
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    sc::replay::WFinalLogWeightContext ctx(owner);
    const auto n = static_cast<cand_pos_t>(sequence.size());
    const Item root{NT::W, 1, n, -1, -1};
    if (bounded_only) {
      // No exact schedule/chart/unrestricted outer is constructed on this path.
      auto summary = output(std::string(argv[3]) + ".summary.tsv");
      auto work = output(std::string(argv[3]) + ".bounded.tsv");
      bounded_header(summary); bounded_header(work);
      auto empty_ranking = output(std::string(argv[3]) + ".ranking.tsv");
      empty_ranking << "mode\n";
      for (int p : pilot_widths) {
        auto widths = sc::replay::WFinalVpWmbpOnlineWidths::uniform(p);
        const auto run = sc::replay::run_w_final_vp_wmbp_online<Sum>(ctx, n, widths);
        for (int factor : bounded_factors) {
          const auto result = bounded_outer(run.chart, root, ctx, p, std::size_t(p) * factor);
          bounded_row(summary, n, p, factor, result);
          bounded_row(work, n, p, factor, result);
        }
      }
      return 0;
    }
    const auto schedule = sc::engine::build_pkfree_basic_exact_schedule(n, ctx);
    const auto exact = sc::replay::run_w_final_exact_dp_schedule<Sum>(
        schedule, ctx, [](const auto &) { return true; },
        [](const auto &, const auto &, const auto &, const auto &, auto) {});
    const double log_z = exact.get(root);
    if (!std::isfinite(log_z)) throw std::runtime_error("unreachable exact root");
    // Evaluation truth only: one regeneration supplies both exact beta and
    // the hindsight max-spine baseline. Never supplied to a pilot estimator.
    Chart beta(zero), max_exact(zero);
    beta.set(root, 0); max_exact.set(root, 0);
    sc::engine::AllLocalWeightModel<decltype(ctx)> weights(&ctx);
    sc::engine::RuntimeDeductionChildItem child_item;
    for (auto cursor = schedule.rbegin(); cursor != schedule.rend(); ++cursor) {
      if (beta.get(*cursor) == zero) continue;
      for (const auto &d : sc::deductions::generate_exact_basic_deductions_for_item(*cursor, ctx)) {
        if (!sc::engine::legacy_span_schedule_supports(d)) continue;
        bool live = true;
        for (std::size_t c = 0; c < d.child_count; ++c)
          if (exact.get(child_item(d, c)) == zero) live = false;
        if (!live) continue;
        const auto local = [&](const auto &source) { return weights.local_weight_for(source); };
        sc::engine::propagate_outside_deduction<Sum>(d, exact, beta.get(*cursor), &beta, local, child_item);
        sc::engine::propagate_outside_deduction<Max>(d, exact, max_exact.get(*cursor), &max_exact, local, child_item);
      }
    }
    std::vector<Row> rows;
    std::map<std::pair<NT, int>, std::vector<std::size_t>> buckets;
    for (const auto &[item, alpha] : exact.values()) {
      if (alpha == zero) continue;
      buckets[{item.nonterminal, item.j}].push_back(rows.size());
      rows.push_back({item, alpha, std::exp(alpha + beta.get(item) - log_z),
          item.i <= 1 ? 0 : exact.get({NT::W, 1, item.i - 1, -1, -1})});
    }
    auto ranking = output(std::string(argv[3]) + ".ranking.tsv");
    ranking << "pilot_width\tedge_cap\tpreserve_pilot\treserve_factor\tspan_factor\tpath_closure\tfamily\twidth\tmode\tcompetitive_buckets"
               "\tcandidates\texpected_occurrences\tretained_occurrences\toracle_retained"
               "\tcovered_occurrences\tcovered_pilot_omitted_occurrences\tpilot_omitted_occurrences\tallowed_occurrences\tbounded_factor\n";
    auto summary = output(std::string(argv[3]) + ".summary.tsv");
    summary << "n\tlog_z\tpilot_width\tpilot_log_z\tpilot_items\tpilot_seconds\tedge_cap"
               "\tpreserve_pilot\treserve_factor\tspan_factor\tpath_closure\tallowed_items\touter_items\touter_pilot_omitted_items\tgenerated_deductions"
               "\teligible_edges\tpropagation_offers\tmax_parent_candidates\touter_seconds"
               "\tstrict_outer_items\tstrict_outer_seconds\treference_error\ttotal_seconds"
               "\tpeak_pending\tpeak_guide_states\tfrontier_evictions\tbounded_factor\n";
    auto items = output(std::string(argv[3]) + ".items.tsv");
    items << "pilot_width\tedge_cap\tpreserve_pilot\treserve_factor\tspan_factor\tpath_closure\tfamily\ti\tj\tip\tjp\tposterior\tpilot_present\tcandidate_allowed\touter\touter_rank\tbounded_factor\n";
    auto bridges = output(std::string(argv[3]) + ".bridge_paths.tsv");
    bridges << "pilot_width\treserve_factor\tcategory\tseed_family\tseed_i\tseed_j\tseed_q"
               "\tstep\tfamily\ti\tj\tip\tjp\tpilot_present\treserve_present\trestricted_outer\tfull_outer\n";
    auto rule_counts = output(std::string(argv[3]) + ".rules.tsv");
    rule_counts << "pilot_width\tedge_cap\tpreserve_pilot\treserve_factor\tspan_factor\tpath_closure"
                   "\trule\tgenerated_deductions\teligible_occurrences\tbounded_factor\n";
    auto bounded_work = output(std::string(argv[3]) + ".bounded.tsv");
    bounded_header(bounded_work);
    for (const int pilot_width : pilot_widths) {
      const auto pilot_start = Clock::now();
      auto widths = sc::replay::WFinalVpWmbpOnlineWidths::uniform(pilot_width);
      const int max_factor = reserve_factors.empty() ? 1 :
          *std::max_element(reserve_factors.begin(), reserve_factors.end());
      CandidateReserve reserve(n, std::size_t(pilot_width) * max_factor);
      const auto pilot_run = sc::replay::run_w_final_vp_wmbp_online_observed<Sum>(ctx, n, widths, reserve);
      const auto &pilot = pilot_run.chart;
      const double pilot_seconds = seconds(pilot_start);
      if (!std::isfinite(pilot.get(root))) throw std::runtime_error("unreachable pilot root");
      for (const auto &[item, alpha] : pilot.values())
        check_le(alpha, exact.get(item), "pilot alpha exceeds exact");
      auto strict = pilot_outer(schedule, pilot, root, ctx, 0, true);
      auto all = pilot_outer(schedule, pilot, root, ctx, 0, false);
      double reference_error = 0;
      if (n <= 32) {
        const auto reference = sc::engine::run_regenerated_outside<Max>(
            schedule, pilot, root, [&](const auto &item) {
              auto ds = sc::deductions::generate_exact_basic_deductions_for_item(item, ctx);
              ds.erase(std::remove_if(ds.begin(), ds.end(), [](const auto &d) {
                return !sc::engine::legacy_span_schedule_supports(d);
              }), ds.end());
              return ds;
            }, [&](const auto &d) { return weights.local_weight_for(d); }, child_item);
        if (reference.values().size() != all.values.values().size())
          throw std::runtime_error("shared engine support differs");
        for (const auto &[key, value] : reference.values()) {
          check_le(value, all.values.get(key), "shared engine value differs");
          check_le(all.values.get(key), value, "shared engine value differs");
          reference_error = std::max(reference_error, std::fabs(value - all.values.get(key)));
        }
        bool same_inside = true;
        for (const auto &row : rows)
          if (pilot.get(row.key) == zero || std::fabs(pilot.get(row.key) - row.alpha) > 1e-9)
            same_inside = false;
        if (same_inside)
          for (const auto &row : rows) {
            check_le(max_exact.get(row.key), all.values.get(row.key), "exact pilot lost a max context");
            check_le(all.values.get(row.key), max_exact.get(row.key), "exact pilot gained a max context");
          }
      }
      std::vector<std::tuple<int, bool, int, int, bool, int>> variants;
      for (const int cap : caps) {
        variants.emplace_back(cap, false, 0, 0, false, 0);
        if (cap) variants.emplace_back(cap, true, 0, 0, false, 0);
      }
      for (const int factor : reserve_factors) {
        variants.emplace_back(0, true, factor, 0, false, 0);
        variants.emplace_back(0, true, factor, 0, true, 0);
      }
      for (const int factor : span_factors) variants.emplace_back(0, true, 0, factor, false, 0);
      for (const int factor : bounded_factors) variants.emplace_back(0, false, 0, 0, false, factor);
      std::map<int, Chart> plain_by_cap;
      for (const auto &[cap, preserve_pilot, reserve_factor, span_factor, path_closure, bounded_factor] : variants) {
        ItemSet allowed;
        if (reserve_factor) {
          allowed = reserve.allowed(pilot, std::size_t(pilot_width) * reserve_factor);
          const auto bound = pilot.values().size() + CandidateReserve::families * n * std::size_t(pilot_width) * reserve_factor;
          if (allowed.size() > bound) throw std::runtime_error("reserve state bound violated");
        }
        Outer capped;
        const Outer *outer = &all;
        if (path_closure) {
          const auto closure_start = Clock::now();
          allowed = close_best_paths(std::move(allowed), all, root);
          for (const auto &key : allowed) {
            const double value = all.values.get(key);
            if (value == zero) continue;
            capped.values.set(key, value);
            const auto link = all.links.find(key);
            if (link != all.links.end()) capped.links.emplace(key, link->second);
          }
          capped.elapsed = seconds(closure_start);
          capped.visited = capped.values.values().size();
          if (n <= 32) {
            const auto regenerated = pilot_outer(schedule, pilot, root, ctx, 0, false, false, &allowed);
            if (regenerated.values.values().size() != capped.values.values().size())
              throw std::runtime_error("path closure support mismatch");
            for (const auto &[key, value] : capped.values.values()) {
              check_le(value, regenerated.values.get(key), "path closure lost max");
              check_le(regenerated.values.get(key), value, "path closure gained max");
            }
          }
          outer = &capped;
        } else if (cap || reserve_factor) {
          capped = pilot_outer(schedule, pilot, root, ctx, cap, false, preserve_pilot,
              reserve_factor ? &allowed : nullptr);
          outer = &capped;
        }
        if (reserve_factor && !path_closure)
          report_bridge_paths(bridges, rows, allowed, *outer, all, pilot_width, reserve_factor, pilot);
        if (span_factor) {
          capped = span_outer(schedule, pilot, root, ctx, std::size_t(pilot_width) * span_factor);
          outer = &capped;
        }
        if (bounded_factor) {
          auto result = bounded_outer(pilot, root, ctx, pilot_width,
                                       std::size_t(pilot_width) * bounded_factor, true);
          bounded_row(bounded_work, n, pilot_width, bounded_factor, result);
          capped = std::move(result.outer);
          outer = &capped;
        }
        if (cap && !preserve_pilot) plain_by_cap.emplace(cap, outer->values);
        std::size_t omitted_items = 0;
        for (const auto &[key, value] : outer->values.values()) {
          if (pilot.get(key) == zero) ++omitted_items;
          check_le(value, all.values.get(key), "capped outer exceeds unlimited outer");
          if (exact.get(key) != zero)
            check_le(value, max_exact.get(key), "pilot outer exceeds exact max-spine");
        }
        for (const auto &[key, value] : strict.values.values()) {
          check_le(value, all.values.get(key), "strict outer exceeds one-hole outer");
          if (preserve_pilot) check_le(value, outer->values.get(key), "pilot backbone lost a strict path");
        }
        if (cap && preserve_pilot)
          for (const auto &[key, value] : plain_by_cap.at(cap).values())
            check_le(value, outer->values.get(key), "pilot preservation lost a capped path");
        if (reserve_factor && outer->values.values().size() > allowed.size())
          throw std::runtime_error("outer escaped candidate reserve");
        const std::vector<std::string> modes{
            "raw", "w_prefix", "zero_fallback", "exact_max", "oracle",
            "pilot_strict", "outer", "outer_fallback"};
        constexpr std::size_t mode_count = 8;
        std::vector<std::array<int, mode_count>> ranks(rows.size());
        for (const auto &[bucket, indices] : buckets) {
          (void)bucket;
          const auto score = [&](std::size_t r, std::size_t mode) {
            const auto &row = rows[r];
            const double guide = outer->values.get(row.key);
            switch (mode) {
              case 0: return row.alpha;
              case 1: return row.alpha + row.prefix;
              case 2: return row.alpha + (row.prefix == zero ? 0 : row.prefix);
              case 3: return row.alpha + max_exact.get(row.key);
              case 4: return row.alpha + beta.get(row.key);
              case 5: return row.alpha + strict.values.get(row.key);
              case 6: return row.alpha + guide;
              default: return row.alpha + (guide == zero ? (row.prefix == zero ? 0 : row.prefix) : guide);
            }
          };
          for (std::size_t mode = 0; mode < mode_count; ++mode) {
            auto sorted = indices;
            std::sort(sorted.begin(), sorted.end(), [&](auto a, auto b) {
              const double x = score(a, mode), y = score(b, mode);
              if (x != y) return x > y;
              return less_key(rows[a].key, rows[b].key);
            });
            for (std::size_t r = 0; r < sorted.size(); ++r) ranks[sorted[r]][mode] = r + 1;
          }
        }
        for (const int width : rank_widths) {
          struct Stats {
            std::size_t buckets = 0, candidates = 0;
            double total = 0, covered = 0, recovered = 0, omitted = 0, allowed = 0;
            std::array<double, mode_count> retained{};
          };
          std::map<NT, Stats> stats;
          for (const auto &[bucket, indices] : buckets) {
            if (indices.size() <= static_cast<std::size_t>(width)) continue;
            double total = 0;
            for (auto r : indices) total += rows[r].q;
            if (total == 0) continue;
            auto &st = stats[bucket.first];
            ++st.buckets; st.candidates += indices.size(); st.total += total;
            for (auto r : indices) {
              const auto &row = rows[r];
              const bool covered = outer->values.get(row.key) != zero;
              const bool omitted = pilot.get(row.key) == zero;
              if (!reserve_factor || allowed.count(row.key)) st.allowed += row.q;
              if (covered) st.covered += row.q;
              if (omitted) st.omitted += row.q;
              if (covered && omitted) st.recovered += row.q;
              for (std::size_t mode = 0; mode < mode_count; ++mode)
                if (ranks[r][mode] <= width) st.retained[mode] += row.q;
            }
          }
          for (const auto &[family, st] : stats)
            for (std::size_t mode = 0; mode < mode_count; ++mode)
              ranking << pilot_width << '\t' << cap << '\t' << preserve_pilot << '\t' << reserve_factor << '\t' << span_factor << '\t' << path_closure << '\t' << sc::grammar::nonterminal_name(family)
                      << '\t' << width << '\t' << modes[mode] << '\t' << st.buckets << '\t' << st.candidates
                      << '\t' << st.total << '\t' << st.retained[mode] << '\t' << st.retained[4]
                      << '\t' << st.covered << '\t' << st.recovered << '\t' << st.omitted
                      << '\t' << st.allowed << '\t' << bounded_factor << '\n';
        }
        for (std::size_t r = 0; r < rows.size(); ++r) {
          const auto &row = rows[r];
          if (row.q < 1e-6) continue;
          items << pilot_width << '\t' << cap << '\t' << preserve_pilot << '\t' << reserve_factor << '\t' << span_factor << '\t' << path_closure << '\t' << sc::grammar::nonterminal_name(row.key.nonterminal)
                << '\t' << row.key.i << '\t' << row.key.j << '\t' << row.key.ip << '\t' << row.key.jp
                << '\t' << row.q << '\t' << (pilot.get(row.key) != zero)
                << '\t' << (!reserve_factor || allowed.count(row.key)) << '\t'
                << outer->values.get(row.key) << '\t' << ranks[r][6] << '\t' << bounded_factor << '\n';
        }
        summary << n << '\t' << log_z << '\t' << pilot_width << '\t' << pilot.get(root)
                << '\t' << pilot.values().size() << '\t' << pilot_seconds << '\t' << cap
                << '\t' << preserve_pilot << '\t' << reserve_factor << '\t' << span_factor << '\t' << path_closure << '\t' << allowed.size()
                << '\t' << outer->values.values().size() << '\t' << omitted_items
                << '\t' << outer->generated << '\t' << outer->eligible << '\t' << outer->selected
                << '\t' << outer->max_parent_candidates << '\t' << outer->elapsed
                << '\t' << strict.values.values().size() << '\t' << strict.elapsed
                << '\t' << reference_error << '\t' << seconds(total_start)
                << '\t' << outer->peak_pending << '\t' << outer->peak_guide_states << '\t' << outer->evicted << '\t' << bounded_factor << '\n';
        for (const auto rule : sc::grammar::all_rule_ids()) {
          const auto idx = static_cast<std::size_t>(rule);
          if (!outer->generated_rules[idx]) continue;
          rule_counts << pilot_width << '\t' << cap << '\t' << preserve_pilot << '\t' << reserve_factor
                      << '\t' << span_factor << '\t' << path_closure << '\t' << sc::grammar::rule_name(rule)
                      << '\t' << outer->generated_rules[idx] << '\t' << outer->eligible_rules[idx] << '\t' << bounded_factor << '\n';
        }
        std::cerr << "pilot=" << pilot_width << " cap=" << cap << " preserve_pilot=" << preserve_pilot
                  << " reserve_factor=" << reserve_factor << " span_factor=" << span_factor
                  << " path_closure=" << path_closure << " bounded_factor=" << bounded_factor << " items="
                  << outer->values.values().size() << " edges=" << outer->selected
                  << " seconds=" << outer->elapsed << '\n';
      }
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

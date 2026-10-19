#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_loop_retained_events.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using namespace cparty::scfg2;
using Key = deductions::ItemKey;
using Rule = grammar::RuleId;
using Source = deductions::Deduction;
using Pair = std::pair<cand_pos_t, cand_pos_t>;

void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}

// Independent predicate controls isolate boundary decisions without imposing
// additional pairing or closure policy on the event producer. In particular,
// can_pair does not enforce TURN: the shared generator owns its explicit gates.
struct SyntheticVpContext final : deductions::PairabilityContext {
  struct Borders { cand_pos_t b = -1, bp = -1, B = -1, Bp = -1; };
  cand_pos_t n = 48;
  std::vector<cand_pos_t> partners = std::vector<cand_pos_t>(n + 1, -2);
  std::map<Pair, Borders> borders;
  std::set<Pair> forbidden_pairs, closed_spans;

  cand_pos_t turn() const override { return 3; }
  cand_pos_t sequence_length() const override { return n; }
  bool can_pair(cand_pos_t i, cand_pos_t j) const override {
    return i >= 1 && j <= n && i < j && !forbidden_pairs.count({i, j});
  }
  cand_pos_t pair_partner(cand_pos_t p) const override {
    return p >= 1 && p <= n ? partners[p] : -2;
  }
  bool is_unpaired(cand_pos_t p) const override {
    return p >= 1 && p <= n && partners[p] < 0;
  }
  cand_pos_t unpaired_prefix(cand_pos_t p) const override {
    cand_pos_t count = 0;
    for (; p >= 1 && is_unpaired(p); --p) ++count;
    return count;
  }
  bool is_empty_region(cand_pos_t i, cand_pos_t j) const override {
    if (j <= i) return false;
    for (cand_pos_t p = i + 1; p < j; ++p)
      if (!is_unpaired(p)) return false;
    return true;
  }
  bool weakly_closed(cand_pos_t i, cand_pos_t j) const override {
    return closed_spans.count({i, j});
  }
  Borders border_values(cand_pos_t i, cand_pos_t j) const {
    const auto found = borders.find({i, j});
    return found == borders.end() ? Borders{} : found->second;
  }
  cand_pos_t border_b(cand_pos_t i, cand_pos_t j) const override {
    return border_values(i, j).b;
  }
  cand_pos_t border_bp(cand_pos_t i, cand_pos_t j) const override {
    return border_values(i, j).bp;
  }
  cand_pos_t border_B(cand_pos_t i, cand_pos_t j) const override {
    return border_values(i, j).B;
  }
  cand_pos_t border_Bp(cand_pos_t i, cand_pos_t j) const override {
    return border_values(i, j).Bp;
  }
};

void same_source(const Source &actual, const Source &expected) {
  require(deductions::ItemKeyEq{}(actual.parent, expected.parent),
          "VP event changed parent identity");
  require(std::tie(actual.rule, actual.split.k, actual.split.l,
                   actual.split.p, actual.split.q) ==
              std::tie(expected.rule, expected.split.k, expected.split.l,
                       expected.split.p, expected.split.q),
          "VP event changed source identity or first-arrival order");
  require(actual.child_count == expected.child_count,
          "VP event changed child arity");
  for (std::size_t c = 0; c < actual.child_count; ++c)
    require(deductions::ItemKeyEq{}(actual.children[c].item,
                                  expected.children[c].item),
            "VP event changed child identity");
}

// Expected admission of the sole physical VP state.
struct EventExpectation {
  Pair parent, child;
  Rule rule;
  unsigned carriers;
};

void independent_event_case(SyntheticVpContext &ctx,
                            const std::vector<Pair> &children,
                            const std::vector<EventExpectation> &checks) {
  constexpr grammar::NonTerminal carriers[]{
      grammar::NonTerminal::VP};
  std::map<cand_pos_t, std::vector<cand_pos_t>> columns;
  for (const auto &[i, j] : children) columns[j].push_back(i);
  std::vector<std::pair<cand_pos_t, std::vector<cand_pos_t>>> deliveries(
      columns.begin(), columns.end());

  for (std::size_t carrier_index = 0; carrier_index < 1; ++carrier_index) {
    const auto carrier = carriers[carrier_index];
    for (bool reversed : {false, true}) {
      auto ordered = deliveries;
      if (reversed) {
        std::reverse(ordered.begin(), ordered.end());
        for (auto &[j, starts] : ordered)
          std::reverse(starts.begin(), starts.end());
      }
      beam::VpLoopRetainedEvents<SyntheticVpContext> events(carrier);
      std::map<Pair, std::size_t> arrival;
      const auto deliver = [&](cand_pos_t j, const std::vector<cand_pos_t> &starts) {
        for (const auto i : starts) arrival.emplace(Pair{i, j}, arrival.size());
        events.observe_committed_column(j, starts, ctx.n, ctx);
      };
      for (const auto &[j, starts] : ordered) {
        auto repeated = starts;
        repeated.push_back(starts.front());
        deliver(j, repeated);
        std::reverse(repeated.begin(), repeated.end());
        deliver(j, repeated);
      }
      // Revisit earlier columns after all other first arrivals as well.
      for (auto column = ordered.rbegin(); column != ordered.rend(); ++column)
        deliver(column->first, column->second);

      // Full-parent generation is independent of the producer's reverse
      // window, direct split predicate, and duplicate representation. Only
      // delivered children determine which complete sources should survive.
      std::map<Pair, std::vector<Source>> expected;
      std::map<cand_pos_t, std::vector<std::pair<std::size_t, cand_pos_t>>> parent_order;
      std::size_t stacks = 0, internals = 0;
      for (cand_pos_t j = 1; j <= ctx.n; ++j) {
        for (cand_pos_t i = 1; i <= j; ++i) {
          const Key parent{carrier, i, j, -1, -1};
          std::vector<Source> sources;
          for (const auto &source :
               deductions::generate_exact_basic_deductions_for_item(parent, ctx)) {
            if (source.rule != Rule::VP_STACK && source.rule != Rule::VP_INTERNAL_LOOP)
              continue;
            require(source.child_count == 1, "VP loop reference is not unary");
            const auto &child = source.children[0].item;
            require(child.nonterminal == carrier, "VP loop reference changed carrier");
            if (arrival.count({child.i, child.j})) sources.push_back(source);
          }
          const auto first_arrival = [&](const Source &source) {
            const auto &child = source.children[0].item;
            return arrival.at({child.i, child.j});
          };
          std::sort(sources.begin(), sources.end(), [&](const Source &a, const Source &b) {
            const auto a_order = first_arrival(a), b_order = first_arrival(b);
            if (a_order != b_order) return a_order < b_order;
            return a.rule == Rule::VP_STACK && b.rule != Rule::VP_STACK;
          });
          if (!sources.empty()) {
            parent_order[j].emplace_back(first_arrival(sources.front()), i);
            for (const auto &source : sources) {
              stacks += source.rule == Rule::VP_STACK;
              internals += source.rule == Rule::VP_INTERNAL_LOOP;
            }
            expected.emplace(Pair{i, j}, std::move(sources));
          }
          const auto found = expected.find({i, j});
          const std::vector<Source> empty;
          const auto &want = found == expected.end() ? empty : found->second;
          const auto &actual = events.deductions_for(parent);
          require(actual.size() == want.size(),
                  "VP event source support changed under repeated/reversed deliveries");
          for (std::size_t s = 0; s < actual.size(); ++s) same_source(actual[s], want[s]);
        }
        auto &order = parent_order[j];
        std::sort(order.begin(), order.end());
        std::vector<cand_pos_t> starts;
        for (const auto &[first, i] : order) starts.push_back(i);
        require(events.parent_starts(j) == starts, "VP event parent order changed");
      }
      require(events.stats().stack_sources == stacks &&
                  events.stats().internal_sources == internals,
              "VP event source counts changed under repeated/reversed deliveries");
      // parent_probes intentionally remains free to shrink with the producer.
      for (const auto &check : checks) {
        const Key parent{carrier, check.parent.first, check.parent.second, -1, -1};
        std::size_t count = 0;
        for (const auto &source : events.deductions_for(parent)) {
          const auto &child = source.children[0].item;
          count += source.rule == check.rule &&
                   Pair{child.i, child.j} == check.child;
        }
        const bool present = (check.carriers & (1u << carrier_index)) != 0;
        require(count == (present ? 1u : 0u), "VP boundary witness was not exercised");
      }
    }
  }
}

void independent_event_witnesses() {
  constexpr unsigned all = 1, closed_only = 0;
  SyntheticVpContext free;
  independent_event_case(free,
      {{18, 30}, {19, 30}, {17, 31}, {2, 47}, {2, 46}, {3, 47},
       {31, 36}, {32, 36}, {3, 8}, {2, 3}, {2, 4}, {3, 4}, {2, 2},
       {1, 12}, {20, 48}},
      {{{17, 31}, {18, 30}, Rule::VP_STACK, all},
       {{17, 31}, {18, 30}, Rule::VP_INTERNAL_LOOP, 0},
       {{2, 45}, {18, 30}, Rule::VP_INTERNAL_LOOP, all}, // combined 29
       {{1, 45}, {18, 30}, Rule::VP_INTERNAL_LOOP, 0},   // combined 30
       {{1, 37}, {31, 36}, Rule::VP_INTERNAL_LOOP, all}, // left 29
       {{1, 37}, {32, 36}, Rule::VP_INTERNAL_LOOP, 0},   // left 30
       {{2, 38}, {3, 8}, Rule::VP_INTERNAL_LOOP, all},   // right 29
       {{2, 39}, {3, 8}, Rule::VP_INTERNAL_LOOP, 0},     // right 30
       {{1, 48}, {2, 47}, Rule::VP_STACK, all},
       {{1, 48}, {2, 46}, Rule::VP_INTERNAL_LOOP, all},
       {{1, 48}, {3, 47}, Rule::VP_INTERNAL_LOOP, all},
       {{1, 4}, {2, 3}, Rule::VP_STACK, closed_only},    // parent span == TURN
       {{1, 5}, {2, 4}, Rule::VP_STACK, all},            // parent span == TURN+1
       {{1, 3}, {2, 2}, Rule::VP_STACK, 0},             // below TURN
       {{2, 6}, {3, 4}, Rule::VP_INTERNAL_LOOP, 0},     // k past TURN edge
       {{2, 7}, {3, 4}, Rule::VP_INTERNAL_LOOP, all}}); // k on TURN edge

  const Pair parent{3, 44}, child{18, 30};
  // Strict left/right borders, including each absent-left sentinel separately.
  for (const auto &bounds : std::vector<std::pair<SyntheticVpContext::Borders, bool>>{
           {{18, -1, -1, 18}, false}, {{19, 29, 29, 19}, true},
           {{19, 30, 29, 19}, false}, {{19, 29, 30, 19}, false},
           {{18, -1, -1, -1}, true}, {{-1, -1, -1, 18}, true}}) {
    SyntheticVpContext ctx;
    ctx.borders[parent] = bounds.first;
    independent_event_case(ctx, {child},
        {{parent, child, Rule::VP_INTERNAL_LOOP, bounds.second ? all : 0u}});
  }
  for (bool right : {false, true}) {
    SyntheticVpContext ctx;
    const auto p = right ? 35 : 10;
    ctx.partners[p] = p + 2; ctx.partners[p + 2] = p;
    independent_event_case(ctx, {child, {4, 43}},
        {{parent, child, Rule::VP_INTERNAL_LOOP, 0},
         {parent, {4, 43}, Rule::VP_STACK, all}});
  }
  {
    SyntheticVpContext ctx;
    ctx.partners[18] = 20; ctx.partners[20] = 18;
    independent_event_case(ctx, {child},
        {{parent, child, Rule::VP_INTERNAL_LOOP, 0},
         {{17, 31}, child, Rule::VP_STACK, 0}});
  }
  {
    SyntheticVpContext ctx;
    ctx.partners[parent.first] = -1; // forbidden public endpoint, still unpaired
    independent_event_case(ctx, {child, {4, 43}},
        {{parent, child, Rule::VP_INTERNAL_LOOP, closed_only},
         {parent, {4, 43}, Rule::VP_STACK, closed_only}});
  }
  {
    SyntheticVpContext ctx;
    ctx.partners[parent.first] = parent.second;
    ctx.partners[parent.second] = parent.first;
    independent_event_case(ctx, {child, {4, 43}},
        {{parent, child, Rule::VP_INTERNAL_LOOP, 0},
         {parent, {4, 43}, Rule::VP_STACK, 0}});
  }
  for (const auto &blocked : {parent, child}) {
    SyntheticVpContext ctx;
    ctx.forbidden_pairs.insert(blocked);
    independent_event_case(ctx, {child},
        {{parent, child, Rule::VP_INTERNAL_LOOP, 0}});
  }
  {
    SyntheticVpContext ctx;
    ctx.closed_spans.insert(parent);
    independent_event_case(ctx, {child},
        {{parent, child, Rule::VP_INTERNAL_LOOP, closed_only}});
  }
}

void streaming_lifetime_witness() {
  SyntheticVpContext ctx;
  ctx.n = 120;
  ctx.partners.assign(ctx.n + 1, -2);
  for (const auto carrier : {grammar::NonTerminal::VP}) {
    beam::VpLoopRetainedEvents<SyntheticVpContext> history(carrier), stream(carrier);
    std::size_t released_parents = 0;
    for (cand_pos_t j = 1; j <= ctx.n; ++j) {
      // Compare all sources and their arrival order before consumption.
      require(history.parent_starts(j) == stream.parent_starts(j),
              "column release changed future parent starts");
      for (const auto i : history.parent_starts(j)) {
        const Key parent{carrier, i, j, -1, -1};
        const auto &a = history.deductions_for(parent);
        const auto &b = stream.deductions_for(parent);
        require(a.size() == b.size(), "column release changed future source support");
        for (std::size_t k = 0; k < a.size(); ++k) same_source(a[k], b[k]);
        ++released_parents;
      }
      const std::vector<cand_pos_t> starts = j > 8
          ? std::vector<cand_pos_t>{j - 7, j - 5, j - 3} : std::vector<cand_pos_t>{};
      history.observe_committed_column(j, starts, ctx.n, ctx);
      stream.observe_committed_column(j, starts, ctx.n, ctx);
      stream.release_column(j);
      stream.release_column(j); // Empty/repeated release is harmless.
      require(stream.parent_starts(j).empty(), "consumed start inventory survived");
      for (const auto i : history.parent_starts(j))
        require(stream.deductions_for({carrier, i, j, -1, -1}).empty(),
                "consumed loop source survived");
      require(stream.pending_column_count() <= MAXLOOP,
              "streaming loop events exceeded the fixed lookahead window");
      require(stream.stats().parent_probes == history.stats().parent_probes &&
                  stream.stats().stack_sources == history.stats().stack_sources &&
                  stream.stats().internal_sources == history.stats().internal_sources,
              "column release changed cumulative work");
    }
    require(released_parents > 0 && stream.pending_parent_count() == 0 &&
                stream.pending_column_count() == 0 && history.pending_parent_count() > 0,
            "streaming lifetime witness did not retire historical storage");
  }
}


bool close(double lhs, double rhs) {
  const double scale = std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  return std::fabs(lhs - rhs) <= 1e-10 * scale;
}

template <class Semiring>
std::size_t check_events(cparty::scfg2::weights::RuntimeSnapshot &owner) {
  using namespace cparty::scfg2;
  replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const auto schedule =
      engine::build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  const auto chart = replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      [](const deductions::Deduction &) { return true; },
      [](const deductions::ItemKey &, const deductions::Deduction &,
         typename Semiring::value_type, const auto &, double) {});
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);
  beam::VpLoopRetainedEvents<replay::WFinalExactInsideContext> events;
  for (cand_pos_t j = 1; j <= n; ++j) {
    std::vector<cand_pos_t> starts;
    for (cand_pos_t i = 1; i <= j; ++i) {
      if (chart.get(deductions::ItemKey{
              grammar::NonTerminal::VP, i, j, -1, -1}) != Semiring::zero()) {
        starts.push_back(i);
      }
    }
    events.observe_committed_column(j, starts, n, ctx);
  }

  std::size_t positive_exact = 0;
  std::size_t positive_event = 0;
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = 1; i <= j; ++i) {
      const deductions::ItemKey parent{
          grammar::NonTerminal::VP, i, j, -1, -1};
      auto exact = Semiring::zero();
      auto driven = Semiring::zero();
      for (const auto &source :
           deductions::generate_exact_basic_deductions_for_item(parent, ctx)) {
        if (source.rule != grammar::RuleId::VP_STACK &&
            source.rule != grammar::RuleId::VP_INTERNAL_LOOP) {
          continue;
        }
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, source, weight.weight_for(source));
        exact = Semiring::plus(exact, term);
        positive_exact += term != Semiring::zero() ? 1 : 0;
      }
      for (const auto &source : events.deductions_for(parent)) {
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, source, weight.weight_for(source));
        driven = Semiring::plus(driven, term);
        positive_event += term != Semiring::zero() ? 1 : 0;
      }
      if (!close(exact, driven)) {
        std::cerr << "loop mismatch parent=" << i << ',' << j
                  << " exact=" << exact << " driven=" << driven << '\n';
        assert(false);
      }
    }
  }
  assert(positive_event == positive_exact);
  return positive_exact;
}

}  // namespace

int main(int argc, char **argv) {
  assert(argc == 1 || argc == 3 || argc == 4);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 4 ? argv[3] : "params/rna_Turner04.par",
         1);
  std::string sequence = argc >= 3
      ? argv[1]
      : "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  std::string scaffold = argc >= 3
      ? argv[2]
      : "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  independent_event_witnesses();
  streaming_lifetime_witness();
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  const auto inside =
      check_events<cparty::scfg2::semirings::SumProductSemiring>(owner);
  const auto viterbi =
      check_events<cparty::scfg2::semirings::MaxProductSemiring>(owner);
  if (argc == 1) {
    assert(inside > 0);
    assert(viterbi > 0);
  }
  std::cout << "scfg2_vp_loop_retained_events=ok\n";
  return EXIT_SUCCESS;
}

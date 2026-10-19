#include "scfg2/beam/source_selection.hh"
#include "scfg2/engine/traceback_pairs.hh"
#include "scfg2/engine/chosen_pair_output.hh"

#include <cstdlib>
#include <iostream>
#include <new>
#include <optional>
#include <type_traits>
#include <utility>
#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace {
bool forbid_allocation = false;
struct NoAllocation {
  NoAllocation() { forbid_allocation = true; }
  ~NoAllocation() { forbid_allocation = false; }
};
void require(bool ok, const char *message) {
  if (!ok) {
    std::cerr << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}
}  // namespace

void *operator new(std::size_t size) {
  if (forbid_allocation) throw std::bad_alloc();
  if (void *memory = std::malloc(size ? size : 1)) return memory;
  throw std::bad_alloc();
}
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete(void *memory, std::size_t) noexcept { std::free(memory); }
void *operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory, std::size_t) noexcept { std::free(memory); }

namespace {
using namespace cparty::scfg2;
using Node = engine::DeductionTraceNode<double>;
using NT = grammar::NonTerminal;
using R = grammar::RuleId;
using Key = deductions::ItemKey;

struct Context : deductions::PairabilityContext {
  cand_pos_t turn() const override { return 3; }
  cand_pos_t sequence_length() const override { return 100000; }
  bool can_pair(cand_pos_t, cand_pos_t) const override { return false; }
  cand_pos_t unpaired_prefix(cand_pos_t i) const override { return i; }
  cand_pos_t pair_partner(cand_pos_t) const override { return -2; }
};

void check_chain(const Node &tree, int n) {
  const Node *node = &tree;
  for (int j = n; j >= 1; --j) {
    require(node->item.nonterminal == NT::W && node->item.i == 1 &&
            node->item.j == j && node->chart_value == 1.0 &&
            node->chosen_contribution == 1.0 && node->chosen_deduction &&
            node->chosen_deduction->rule == R::W_EXTEND_UNPAIRED,
            "deep traceback content changed");
    require(node->children.size() == (j > 1 ? 1u : 0u),
            "deep traceback child count changed");
    if (j > 1) node = &node->children[0];
  }
}

void check_deep_chain(int n) {
  engine::Chart<double> chart(0.0);
  engine::ChoiceMap<double> choices;
  for (int j = 1; j <= n; ++j) {
    const Key item{NT::W, 1, j, -1, -1};
    deductions::Deduction d;
    d.parent = item;
    d.rule = R::W_EXTEND_UNPAIRED;
    if (j > 1) {
      d.child_count = 1;
      d.children[0].item = {NT::W, 1, j - 1, -1, -1};
    }
    chart.set(item, 1.0);
    choices.observe(item, d, 1.0, [](double a, double b) { return a > b; });
  }
  std::optional<Node> tree(engine::build_traceback_tree(
      Key{NT::W, 1, n, -1, -1}, chart, choices));
  check_chain(*tree, n);
  Context ctx;
  const auto chosen = engine::collect_chosen_pair_output(
      Key{NT::W, 1, n, -1, -1}, choices, ctx);
  require(chosen.pairs.empty() && chosen.generated_pairs.empty() &&
          chosen.visited_items == static_cast<std::size_t>(n) &&
          chosen.child_occurrences == static_cast<std::size_t>(n - 1),
          "choice-DAG deep chain changed output or repeated work");
  require(engine::collect_traceback_pairs(*tree, ctx).empty() &&
          engine::collect_generated_traceback_pairs(*tree, ctx).empty(),
          "all-unpaired traceback emitted a pair");
  std::vector<deductions::Deduction> sources;
  beam::collect_selected_source_deductions(*tree, &sources);
  require(sources.size() == static_cast<std::size_t>(n),
          "deep source collector lost a deduction");
  for (int k = 0; k < n; ++k)
    require(sources[k].parent.j == n - k, "deep source preorder changed");

  auto copy = *tree;
  check_chain(copy, n);
  copy.children[0].chart_value = 2.0;
  require(tree->children[0].chart_value == 1.0, "traceback copy aliases source");
  copy.children[0].chart_value = 1.0;
  auto assigned = *tree;
  assigned = copy;
  assigned = assigned;
  check_chain(assigned, n);
  auto moved = std::move(copy);
  check_chain(moved, n);
  require(copy.children.empty(), "move construction copied deep children");
  {
    NoAllocation guard;
    assigned = std::move(moved);
    assigned = std::move(assigned);
  }
  check_chain(assigned, n);
  require(moved.children.empty(), "move assignment copied deep children");
  {
    NoAllocation guard;
    assigned.children = std::move(assigned.children.back().children);
  }
  require(assigned.children[0].item.j == n - 2,
          "move assignment from a descendant lost its buffer");
  {
    NoAllocation guard;
    assigned.children.clear();
    tree.reset();
  }
  require(assigned.children.empty(), "deep clear retained children");
}

void check_branching_order() {
  engine::Chart<double> chart(0.0);
  engine::ChoiceMap<double> choices;
  const Key root{NT::W, 1, 30, -1, -1};
  const Key first{NT::V, 5, 20, -1, -1}, second{NT::V, 2, 25, -1, -1};
  const Key missing{NT::V, 7, 18, -1, -1};
  // A synthetic branching feature tree tests traversal and first-occurrence
  // dedup independently of grammar admission. The long W chain above is a
  // valid all-unpaired derivation.
  deductions::Deduction d;
  d.parent = root;
  d.rule = R::W_ROOT_WMB_DECOMPOSITION;
  d.child_count = 3;
  d.children[0].item = first;
  d.children[1].item = second;
  d.children[2].item = first;
  choices.observe(root, d, 3.0, [](double a, double b) { return a > b; });
  for (const auto item : {first, second}) {
    d = {};
    d.parent = item;
    d.rule = R::V_HAIRPIN;
    if (item.i == second.i) {
      d.child_count = 1;
      d.children[0].item = missing;
    }
    choices.observe(item, d, 2.0, [](double a, double b) { return a > b; });
    chart.set(item, 2.0);
  }
  chart.set(root, 3.0);
  chart.set(missing, 4.0);
  auto tree = engine::build_traceback_tree(root, chart, choices);
  require(!tree.children[1].children[0].chosen_deduction &&
          tree.children[1].children[0].chart_value == 4.0,
          "missing-choice leaf changed");
  Context ctx;
  std::vector<engine::TracebackPair> occurrences;
  engine::collect_traceback_pair_occurrences_into(tree, ctx, &occurrences, false);
  require(occurrences.size() == 3 && occurrences[0].i == 5 &&
          occurrences[1].i == 2 && occurrences[2].i == 5,
          "pair occurrence preorder changed");
  const auto pairs = engine::collect_traceback_pairs(tree, ctx);
  const auto generated = engine::collect_generated_traceback_pairs(tree, ctx);
  const auto chosen = engine::collect_chosen_pair_output(root, choices, ctx);
  require(chosen.visited_items == 4 && chosen.repeated_items == 1 &&
          chosen.pairs.size() == pairs.size() &&
          chosen.generated_pairs.size() == generated.size() &&
          std::equal(pairs.begin(), pairs.end(), chosen.pairs.begin(), engine::traceback_pair_equal) &&
          std::equal(generated.begin(), generated.end(), chosen.generated_pairs.begin(), engine::traceback_pair_equal),
          "choice-DAG changed stable first-occurrence output");
  require(pairs.size() == 2 && pairs[0].i == 5 && pairs[1].i == 2 &&
          generated.size() == 2 && generated[0].i == 5 && generated[1].i == 2,
          "first-occurrence pair dedup changed");
  std::vector<deductions::Deduction> sources;
  beam::collect_selected_source_deductions(tree, &sources);
  require(sources.size() == 4 && sources[0].parent.i == 1 &&
          sources[1].parent.i == 5 && sources[2].parent.i == 2 &&
          sources[3].parent.i == 5, "branching source preorder changed");
}

void check_shared_choice_dag(int depth, bool tree_reference) {
  engine::Chart<double> chart(0.0);
  engine::ChoiceMap<double> choices;
  for (int i = 1; i <= depth; ++i) {
    const Key item{NT::V, i, 2 * depth + 10 - i, -1, -1};
    deductions::Deduction d;
    d.parent = item;
    d.rule = R::V_HAIRPIN;
    if (i < depth) {
      d.child_count = 2;
      d.children[0].item = d.children[1].item =
          Key{NT::V, i + 1, 2 * depth + 9 - i, -1, -1};
    }
    choices.observe(item, d, 1.0, [](double a, double b) { return a > b; });
    chart.set(item, 1.0);
  }
  Context ctx;
  const Key root{NT::V, 1, 2 * depth + 9, -1, -1};
  const auto out = engine::collect_chosen_pair_output(root, choices, ctx);
  require(out.visited_items == static_cast<std::size_t>(depth) &&
          out.repeated_items == static_cast<std::size_t>(depth - 1) &&
          out.child_occurrences == 2 * static_cast<std::size_t>(depth - 1) &&
          out.pairs.size() == static_cast<std::size_t>(depth) &&
          out.generated_pairs.size() == out.pairs.size(),
          "shared choice DAG expanded occurrences or dropped pair identity");
  for (int i = 1; i <= depth; ++i)
    require(out.pairs[i - 1].i == i &&
            out.pairs[i - 1].lane == engine::TracebackPairLane::Round &&
            out.generated_pairs[i - 1].lane == engine::TracebackPairLane::Square,
            "choice-DAG pair lane or preorder changed");
  if (tree_reference) {
    const auto tree = engine::build_traceback_tree(root, chart, choices);
    const auto pairs = engine::collect_traceback_pairs(tree, ctx);
    require(pairs.size() == out.pairs.size() &&
            std::equal(pairs.begin(), pairs.end(), out.pairs.begin(), engine::traceback_pair_equal),
            "shared choice DAG differs from independent occurrence-tree output");
  }
}

void check_stable_pair_identity() {
  using Lane = engine::TracebackPairLane;
  std::vector<engine::TracebackPair> pairs{
      {1, 9, Lane::Square}, {1, 9, Lane::Round}, {1, 8, Lane::Round},
      {1, 9, Lane::Square}, {2, 8, Lane::Round}, {1, 9, Lane::Round}};
  const std::vector<engine::TracebackPair> expected{
      {1, 9, Lane::Square}, {1, 9, Lane::Round}, {1, 8, Lane::Round}, {2, 8, Lane::Round}};
  engine::stable_unique_traceback_pairs(&pairs);
  require(pairs.size() == expected.size() &&
          std::equal(pairs.begin(), pairs.end(), expected.begin(), engine::traceback_pair_equal),
          "stable pair hash changed first occurrence or full lane identity");
}

struct CopyFailure {};
struct CountedValue {
  inline static int live = 0, assignments_until_throw = -1;
  CountedValue() { ++live; }
  CountedValue(const CountedValue &) { ++live; }
  ~CountedValue() { --live; }
  CountedValue &operator=(const CountedValue &) {
    if (assignments_until_throw == 0) throw CopyFailure{};
    if (assignments_until_throw > 0) --assignments_until_throw;
    return *this;
  }
};

void check_copy_unwind() {
  {
    engine::DeductionTraceNode<CountedValue> root;
    auto *node = &root;
    for (int i = 0; i < 10000; ++i) {
      node->children.emplace_back();
      node = &node->children.back();
    }
    const int before = CountedValue::live;
    CountedValue::assignments_until_throw = 500;
    bool threw = false;
    try {
      auto copy = root;
      (void)copy;
    } catch (const CopyFailure &) {
      threw = true;
    }
    CountedValue::assignments_until_throw = -1;
    require(threw && CountedValue::live == before,
            "partial iterative copy leaked live values");
  }
  require(CountedValue::live == 0, "copy-unwind source leaked live values");
}
}  // namespace

int main() {
  static_assert(std::is_aggregate_v<Node>);
  static_assert(std::is_nothrow_move_constructible_v<Node>);
  static_assert(std::is_nothrow_move_assignable_v<Node>);
#if defined(__unix__) || defined(__APPLE__)
  rlimit limit{};
  require(getrlimit(RLIMIT_STACK, &limit) == 0, "cannot read stack limit");
  if (limit.rlim_cur == RLIM_INFINITY || limit.rlim_cur > 1024 * 1024) {
    limit.rlim_cur = 1024 * 1024;
    require(setrlimit(RLIMIT_STACK, &limit) == 0, "cannot constrain stack");
  }
#endif
  check_branching_order();
  check_stable_pair_identity();
  check_shared_choice_dag(8, true);
  check_shared_choice_dag(70, false);  // Occurrence tree would have 2^70-1 nodes.
  check_deep_chain(100000);
  check_copy_unwind();
  std::cout << "scfg2_traceback_lifecycle=ok depth=100000 stack_limit=1048576\n";
}

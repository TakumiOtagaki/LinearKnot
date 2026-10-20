#include "scfg2/deductions/expand_basic.hh"
#include "scfg2/deductions/generate_exact_basic.hh"

#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using namespace cparty::scfg2;
using deductions::Deduction;
using deductions::ItemKey;
using grammar::NonTerminal;
using grammar::RuleId;
using grammar::RuleSplit;

constexpr auto rules = grammar::all_rule_ids();

constexpr bool same_item(const ItemKey &lhs, const ItemKey &rhs) {
  return lhs.nonterminal == rhs.nonterminal && lhs.i == rhs.i &&
         lhs.j == rhs.j && lhs.ip == rhs.ip && lhs.jp == rhs.jp;
}

constexpr bool same_deduction(const Deduction &lhs, const Deduction &rhs) {
  if (!same_item(lhs.parent, rhs.parent) || lhs.rule != rhs.rule ||
      lhs.split.k != rhs.split.k || lhs.split.l != rhs.split.l ||
      lhs.split.p != rhs.split.p || lhs.split.q != rhs.split.q ||
      lhs.child_count != rhs.child_count) {
    return false;
  }
  // Unused children retain their initialized fields too.
  for (std::size_t child = 0; child < lhs.children.size(); ++child) {
    if (!same_item(lhs.children[child].item, rhs.children[child].item))
      return false;
  }
  return true;
}

template <RuleId Rule>
constexpr bool constexpr_expansions_match() {
  constexpr std::array<ItemKey, 4> parents{{
      {grammar::rule_spec(Rule).parent, 3, 18, -1, -1},
      {NonTerminal::VP, 3, 18, -1, -1},
      {NonTerminal::BE, 3, 18, 3, 18},
      {NonTerminal::BE, 3, 18, 6, 14},
  }};
  constexpr std::array<RuleSplit, 6> splits{{
      {}, {3, 18, 4, 17}, {4, 17, 5, 16},
      {5, 16, 4, 17}, {11, 16, 5, 9}, {10, 17, 4, 9},
  }};
  for (const auto &parent : parents) {
    for (const auto &split : splits) {
      if (!same_deduction(deductions::expand_basic_rule(parent, Rule, split),
                          deductions::expand_basic_rule<Rule>(parent, split)))
        return false;
    }
  }
  return true;
}

constexpr ItemKey vp_parent{NonTerminal::VP, 3, 18, -1, -1};
constexpr auto vp_stack = deductions::expand_basic_rule<RuleId::VP_STACK>(
    vp_parent, RuleSplit{});
static_assert(vp_stack.child_count == 1 &&
              same_item(vp_stack.children[0].item,
                        ItemKey{NonTerminal::VP, 4, 17, -1, -1}),
              "fixed expansion must preserve constexpr VP recursion");
constexpr auto same_pair_stack = deductions::expand_basic_rule<RuleId::BE_STACK>(
    ItemKey{NonTerminal::BE, 3, 18, 3, 18}, RuleSplit{});
static_assert(same_item(same_pair_stack.children[0].item,
                        ItemKey{NonTerminal::BE, 4, 17, 4, 17}),
              "fixed expansion must preserve constexpr BE anchor movement");
constexpr auto runtime_constexpr = deductions::expand_basic_rule(
    ItemKey{NonTerminal::V, 3, 18, -1, -1},
    RuleId::V_INTERNAL, RuleSplit{5, 16, -1, -1});
static_assert(runtime_constexpr.child_count == 1 &&
              same_item(runtime_constexpr.children[0].item,
                        ItemKey{NonTerminal::V, 5, 16, -1, -1}),
              "runtime entry point must remain constexpr");

constexpr std::array<NonTerminal, 15> carriers{{
    NonTerminal::W, NonTerminal::WI, NonTerminal::V, NonTerminal::VM,
    NonTerminal::WMv, NonTerminal::WMp, NonTerminal::WM, NonTerminal::WIP,
    NonTerminal::VPL, NonTerminal::VPR, NonTerminal::VP, NonTerminal::WMBW,
    NonTerminal::WMBP, NonTerminal::WMB, NonTerminal::BE,
}};
constexpr std::array<std::array<cand_pos_t, 4>, 8> spans{{
    {{1, 1, -1, -1}}, {{3, 4, -1, -1}}, {{3, 18, -1, -1}},
    {{3, 18, 3, 18}}, {{3, 18, 6, 14}}, {{0, 18, -1, -1}},
    {{9, 3, 6, 8}}, {{-1, -1, -1, -1}},
}};
constexpr std::array<RuleSplit, 18> splits{{
    {}, {2, 3, 1, 2}, {3, 18, 4, 17}, {4, 17, 4, 17},
    {5, 16, 5, 16}, {6, 14, 6, 14}, {4, 17, 5, 17},
    {4, 17, 4, 16}, {4, 17, 5, 16}, {4, 5, 5, 16},
    {4, 6, 6, 16}, {10, 17, 4, 9}, {10, 16, 4, 9},
    {11, 17, 4, 9}, {11, 16, 4, 9}, {10, 17, 5, 9},
    {10, 16, 5, 9}, {11, 17, 5, 9},
}};

template <RuleId Rule>
std::size_t check_rule() {
  static_assert(constexpr_expansions_match<Rule>(),
                "fixed and runtime entry points must agree in constexpr use");
  std::size_t checks = 0;
  for (const auto carrier : carriers) {
    for (const auto &span : spans) {
      const ItemKey parent{carrier, span[0], span[1], span[2], span[3]};
      const auto check_split = [&](const RuleSplit &split) {
        const auto runtime = deductions::expand_basic_rule(parent, Rule, split);
        const auto fixed = deductions::expand_basic_rule<Rule>(parent, split);
        if (!same_deduction(runtime, fixed)) {
          std::cerr << "fixed expansion mismatch rule=" << grammar::rule_name(Rule)
                    << " carrier=" << grammar::nonterminal_name(carrier)
                    << " parent=" << parent.i << ',' << parent.j << ','
                    << parent.ip << ',' << parent.jp
                    << " split=" << split.k << ',' << split.l << ','
                    << split.p << ',' << split.q << '\n';
          std::exit(EXIT_FAILURE);
        }
        ++checks;
      };
      for (const auto &split : splits) check_split(split);
      // Complete the eight empty/nonempty wing combinations for VP case 3.
      check_split(RuleSplit{11, 16, 5, 9});
    }
  }
  return checks;
}

template <std::size_t... Index>
std::size_t check_all_rules(std::index_sequence<Index...>) {
  return (check_rule<rules[Index]>() + ...);
}

// Compare the actual generator's complete ordered stream to the runtime
// expansion oracle, including bounded split views and inactive child slots.
std::size_t check_wm_split_streams() {
  struct Context : deductions::PairabilityContext {
    cand_pos_t turn() const override { return 3; }
    cand_pos_t sequence_length() const override { return 32; }
    bool can_pair(cand_pos_t, cand_pos_t) const override { return true; }
    cand_pos_t unpaired_prefix(cand_pos_t pos) const override { return pos; }
  } ctx;
  const Deduction prefix = deductions::expand_basic_rule<RuleId::V_HAIRPIN>(
      ItemKey{NonTerminal::V, 1, 8, -1, -1}, RuleSplit{});
  for (const auto rule : rules) {
    const ItemKey parent{grammar::rule_spec(rule).parent, 1, 24, -1, -1};
    const auto expected = deductions::generate_exact_basic_deductions(parent, rule, ctx);
    std::vector<Deduction> buffer{prefix};
    deductions::append_exact_basic_deductions(buffer, parent, rule, ctx);
    if (buffer.size() != expected.size() + 1 || !same_deduction(buffer.front(), prefix)) {
      std::cerr << "append lost existing candidate prefix\n";
      std::exit(EXIT_FAILURE);
    }
    for (std::size_t c = 0; c < expected.size(); ++c) {
      if (!same_deduction(buffer[c + 1], expected[c])) {
        std::cerr << "append changed single-rule stream\n";
        std::exit(EXIT_FAILURE);
      }
    }
  }
  const std::vector<cand_pos_t> selected{1, 3, 4, 7, 8, 12, 18, 32};
  const std::vector<cand_pos_t> empty;
  std::size_t checks = 0;
  for (const auto rule : {RuleId::WM_SPLIT_V, RuleId::WM_SPLIT_WMB}) {
    for (int i = 0; i <= 18; ++i) {
      for (int j = i; j <= 32; ++j) {
        const ItemKey parent{NonTerminal::WM, i, j, -1, -1};
        for (const auto *view : {static_cast<const std::vector<cand_pos_t> *>(nullptr),
                                &selected, &empty}) {
          std::vector<Deduction> expected;
          if (deductions::exact_wm_split_parent_allows(parent)) {
            for (int k = i + 1; k <= j; ++k) {
              if (view && !std::binary_search(view->begin(), view->end(), k)) continue;
              const auto ded = deductions::expand_basic_rule(parent, rule, RuleSplit{k});
              if (deductions::is_valid(ded)) expected.push_back(ded);
            }
          }
          const auto generated = deductions::generate_exact_basic_deductions(parent, rule, ctx, view);
          if (expected.size() != generated.size()) {
            std::cerr << "WM split candidate count changed\n";
            std::exit(EXIT_FAILURE);
          }
          for (std::size_t c = 0; c < generated.size(); ++c) {
            if (!same_deduction(expected[c], generated[c])) {
              std::cerr << "WM split candidate fields/order changed\n";
              std::exit(EXIT_FAILURE);
            }
            ++checks;
          }
        }
      }
    }
  }
  return checks;
}

// The item dispatcher may share read-only geometry across fixed VP cases;
// the unchanged single-rule entries provide the independent ordered oracle.
std::size_t check_vp_geometry_streams() {
  struct Context : deductions::PairabilityContext {
    deductions::VpWiGeometry geometry{};
    bool pairs = true, closed = false, left_free = true, right_free = true;
    cand_pos_t left_round = -2, right_round = -2;
    mutable std::size_t parent_queries = 0;
    cand_pos_t turn() const override { return 3; }
    cand_pos_t sequence_length() const override { return 32; }
    bool can_pair(cand_pos_t i, cand_pos_t j) const override {
      return pairs && i == 3 && j == 18;
    }
    cand_pos_t unpaired_prefix(cand_pos_t pos) const override { return pos; }
    cand_pos_t pair_round(cand_pos_t pos) const override {
      return pos == 3 ? left_round : pos == 18 ? right_round : -2;
    }
    bool is_unpaired(cand_pos_t pos) const override {
      return pos == 3 ? left_free : pos == 18 ? right_free : true;
    }
    bool weakly_closed(cand_pos_t, cand_pos_t) const override { return closed; }
    bool is_empty_region(cand_pos_t i, cand_pos_t j) const override { return i < j; }
    cand_pos_t parent_index(cand_pos_t pos) const override {
      ++parent_queries;
      return pos == 3 ? geometry.left_cover : geometry.right_cover;
    }
    cand_pos_t border_b(cand_pos_t, cand_pos_t) const override { return geometry.b; }
    cand_pos_t border_bp(cand_pos_t, cand_pos_t) const override { return geometry.bp; }
    cand_pos_t border_B(cand_pos_t, cand_pos_t) const override { return geometry.B; }
    cand_pos_t border_Bp(cand_pos_t, cand_pos_t) const override { return geometry.Bp; }
  };
  constexpr std::array<RuleId, 9> vp_rules{{
      RuleId::VP_WI_CASE1, RuleId::VP_WI_CASE2, RuleId::VP_WI_CASE3,
      RuleId::VP_STACK, RuleId::VP_INTERNAL_LOOP, RuleId::VP_WIP_VP_LEFT,
      RuleId::VP_VP_WIP_RIGHT, RuleId::VP_WIP_VPR, RuleId::VP_VPL_WIP}};
  std::vector<deductions::VpWiGeometry> geometries{
      // Cases 1 and 2, missing/zero borders and invalid fixed splits.
      {2, 1, 8, 12, -1, -1}, {1, 2, -1, -1, 8, 12},
      {-1, -1, -1, -1, -1, -1}, {0, 0, 0, 0, 0, 0},
      {2, 1, 19, 20, -1, -1}, {1, 2, -1, -1, 2, 1},
      {2, 1, 0, 0, -1, -1}, {1, 2, -1, -1, 0, 0},
      {2, 1, -1, 12, 8, -1}, {1, 2, -1, 12, 8, -1}};
  // All eight empty/nonempty wing combinations of Case 3, including its
  // boundary and nonleft variants. Geometry fields are Bp, B, b, bp.
  for (bool left_empty : {false, true})
    for (bool middle_empty : {false, true})
      for (bool right_empty : {false, true})
        geometries.push_back({2, 2, left_empty ? 4 : 5, 9,
                              middle_empty ? 10 : 11, right_empty ? 17 : 16});
  const std::vector<cand_pos_t> selected{1, 4, 5, 9, 10, 11, 16, 18, 32};
  const std::vector<cand_pos_t> empty;
  std::array<std::size_t, 3> case_counts{};
  std::size_t checks = 0, shared_observations = 0;
  const auto fail = [](const char *message) {
    std::cerr << message << '\n';
    std::exit(EXIT_FAILURE);
  };
  for (const auto &geometry : geometries) {
    for (int gate = 0; gate < 8; ++gate) {
      Context reference;
      reference.geometry = geometry;
      reference.pairs = gate != 1;
      reference.closed = gate == 2;
      reference.left_free = gate != 3;
      reference.right_free = gate != 4;
      reference.left_round = gate == 5 ? -1 : gate == 6 ? 18 : -2;
      reference.right_round = gate == 7 ? -1 : -2;
      for (const auto &parent : {ItemKey{NonTerminal::VP, 3, 18, -1, -1},
                                ItemKey{NonTerminal::VP, 3, 18, 4, 17}}) {
        for (const auto *view : {static_cast<const std::vector<cand_pos_t> *>(nullptr),
                                &selected, &empty}) {
          auto shared = reference;
          reference.parent_queries = shared.parent_queries = 0;
          const auto prefix = deductions::expand_basic_rule(
              ItemKey{NonTerminal::V, 1, 8, -1, -1}, RuleId::V_HAIRPIN, RuleSplit{});
          std::vector<Deduction> expected{prefix}, generated{prefix};
          if (deductions::exact_vp_public_item_active(parent, reference)) {
            for (const auto rule : vp_rules)
              deductions::append_exact_basic_deductions(expected, parent, rule, reference, view);
          }
          deductions::append_exact_basic_deductions_for_item(generated, parent, shared, view);
          if (expected.size() != generated.size()) fail("shared VP geometry changed stream size");
          for (std::size_t c = 0; c < generated.size(); ++c) {
            if (!same_deduction(expected[c], generated[c]))
              fail("shared VP geometry changed complete candidate fields/order");
            for (std::size_t rule = 0; rule < case_counts.size(); ++rule)
              if (generated[c].rule == vp_rules[rule]) ++case_counts[rule];
            ++checks;
          }
          if (reference.parent_queries == 6) {
            if (shared.parent_queries != 2) fail("VP geometry was not observed once per parent");
            ++shared_observations;
          } else if (reference.parent_queries != 0 || shared.parent_queries != 0) {
            fail("inactive VP parent observed fixed-case geometry");
          }
        }
      }
    }
  }
  if (!shared_observations) fail("shared VP geometry observation witness is vacuous");
  for (const auto count : case_counts)
    if (!count) fail("shared VP geometry omitted a fixed case witness");
  return checks;
}

}  // namespace

int main() {
  const auto checks = check_all_rules(std::make_index_sequence<rules.size()>{});
  const auto streams = check_wm_split_streams();
  const auto vp_streams = check_vp_geometry_streams();
  std::cout << "scfg2_fixed_rule_expansion=ok rules=" << rules.size()
            << " checks=" << checks << " generated_candidates=" << streams
            << " vp_geometry_candidates=" << vp_streams << '\n';
  return EXIT_SUCCESS;
}

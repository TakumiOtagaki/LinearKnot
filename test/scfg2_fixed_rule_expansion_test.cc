#include "scfg2/deductions/expand_basic.hh"

#include <array>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>

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

}  // namespace

int main() {
  const auto checks = check_all_rules(std::make_index_sequence<rules.size()>{});
  std::cout << "scfg2_fixed_rule_expansion=ok rules=" << rules.size()
            << " checks=" << checks << '\n';
  return EXIT_SUCCESS;
}

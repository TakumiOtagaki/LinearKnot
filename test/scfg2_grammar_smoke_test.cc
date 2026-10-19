#include "scfg2/grammar/nonterminal.hh"
#include "scfg2/grammar/rules.hh"

#include <cassert>
#include <string_view>

namespace {

using cparty::scfg2::grammar::NonTerminal;
using cparty::scfg2::grammar::RuleId;
using cparty::scfg2::grammar::RuleSpec;

void test_nonterminal_names_roundtrip() {
  NonTerminal parsed = NonTerminal::W;
  assert(cparty::scfg2::grammar::parse_nonterminal_name("VP", &parsed));
  assert(parsed == NonTerminal::VP);
  assert(cparty::scfg2::grammar::nonterminal_name(NonTerminal::WMBP) == std::string_view("WMBP"));
  assert(!cparty::scfg2::grammar::parse_nonterminal_name("VP_CASE3_R", &parsed));
  assert(!cparty::scfg2::grammar::parse_nonterminal_name("VP_DIRECT", &parsed));
  assert(!cparty::scfg2::grammar::parse_nonterminal_name("VP_CLOSED", &parsed));
  static_assert(static_cast<unsigned>(NonTerminal::BE) + 1 == 15);
}

void test_rule_names_exist() {
  assert(cparty::scfg2::grammar::rule_name(RuleId::W_SPLIT_V) == std::string_view("W_SPLIT_V"));
  assert(cparty::scfg2::grammar::rule_name(RuleId::BE_WIP_WIP) == std::string_view("BE_WIP_WIP"));
}

void test_all_rule_ids_surface() {
  const auto rules = cparty::scfg2::grammar::all_rule_ids();
  assert(rules.front() == RuleId::V_HAIRPIN);
  assert(rules.back() == RuleId::BE_BASEPAIR_WIP);
  for (const auto rule : rules) {
    assert(cparty::scfg2::grammar::rule_name(rule) != std::string_view("UNKNOWN"));
  }
}

void test_rule_spec_parents_and_arity() {
  const RuleSpec w_split_v = cparty::scfg2::grammar::rule_spec(RuleId::W_SPLIT_V);
  assert(w_split_v.parent == NonTerminal::W);
  assert(w_split_v.child_count == 2);
  assert(w_split_v.children[0].nonterminal == NonTerminal::W);
  assert(w_split_v.children[1].nonterminal == NonTerminal::V);

  const RuleSpec vp_stack = cparty::scfg2::grammar::rule_spec(RuleId::VP_STACK);
  assert(vp_stack.parent == NonTerminal::VP);
  assert(vp_stack.child_count == 1);
  assert(vp_stack.children[0].nonterminal == NonTerminal::VP);
  assert(vp_stack.children[0].uses_outer_span);
  assert(vp_stack.children[0].uses_secondary_span);

  const RuleSpec wmb_empty = cparty::scfg2::grammar::rule_spec(RuleId::WMB_EMPTY);
  assert(wmb_empty.parent == NonTerminal::WMB);
  assert(wmb_empty.child_count == 0);

  const auto direct = cparty::scfg2::grammar::rule_spec(RuleId::WMBP_DIRECT_VP);
  assert(direct.child_count == 1);
  assert(direct.children[0].nonterminal == NonTerminal::VP);
}

}  // namespace

int main() {
  test_nonterminal_names_roundtrip();
  test_rule_names_exist();
  test_all_rule_ids_surface();
  test_rule_spec_parents_and_arity();
  return 0;
}

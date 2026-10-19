#include "scfg2/engine/compact_choice_map.hh"

#include <cstdlib>
#include <iostream>
#include <tuple>

namespace sc = cparty::scfg2;
using D = sc::deductions::Deduction;
void require(bool ok, const char *message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
bool same(const D &a, const D &b) {
  const sc::deductions::ItemKeyEq equal;
  if (!equal(a.parent, b.parent) || a.rule != b.rule || a.child_count != b.child_count ||
      std::tie(a.split.k,a.split.l,a.split.p,a.split.q) !=
      std::tie(b.split.k,b.split.l,b.split.p,b.split.q)) return false;
  for (std::size_t c=0; c<3; ++c)
    if (!equal(a.children[c].item,b.children[c].item)) return false;
  return true;
}
int main() {
  sc::engine::ChoiceMap<double> full;
  sc::engine::CompactChoiceMap<double> compact;
  const auto better = [](double a,double b) { return a > b + 1e-12; };
  const auto tie = [](const D &a,const D &b) { return a.split.k < b.split.k; };
  std::size_t roundtrips=0;
  for (const auto rule : sc::grammar::all_rule_ids()) {
    const auto nt=sc::grammar::rule_spec(rule).parent;
    const sc::deductions::ItemKey parent{nt,1,1000+int(roundtrips),
        nt==sc::grammar::NonTerminal::BE ? 100 : -1,
        nt==sc::grammar::NonTerminal::BE ? 500 : -1};
    auto source=sc::deductions::expand_basic_rule(parent,rule,{20,70,30,80});
    if (rule==sc::grammar::RuleId::WMBP_DIRECT_VP)
      source.children[0].item.nonterminal=sc::grammar::NonTerminal::VP;
    const auto observe = [&](const D &d,double value) {
      full.observe_with_tie_break(parent,d,value,better,tie);
      compact.observe_with_tie_break(parent,d,value,better,tie);
      const auto *expected=full.find(parent);
      const auto actual=compact.find(parent);
      require(actual && same(expected->deduction,actual->deduction) &&
              expected->contribution==actual->contribution,"compact choice differs from full record");
    };
    observe(source,1);
    require(compact.fallback_values().empty(),"canonical rule unexpectedly requires fallback");
    // Noncanonical active or inactive child slots must survive exactly.
    auto exceptional=source;
    exceptional.children[2].item.i+=3;
    observe(exceptional,2);
    require(compact.fallback_values().size()==1,"noncanonical deduction lost its fallback");
    const auto saved=compact.find(parent);
    observe(source,3);
    require(compact.fallback_values().empty(),"replaced fallback leaked");
    require(saved && same(saved->deduction,exceptional),"later reads invalidated an owning result");
    auto near=source;
    near.split.k-=1;
    observe(near,3-1e-13);
    observe(source,3+1e-13);
    compact.erase(parent);
    full.erase(parent);
    require(!compact.find(parent) && compact.fallback_values().empty(),"erase retained a choice");
    ++roundtrips;
  }
  require(sizeof(sc::engine::CompactChoiceMap<double>::Record)*3 <=
          sizeof(sc::engine::ChoiceRecord<double>),"compact record lost its storage bound");
  std::cout << "scfg2_compact_choices=ok rules=" << roundtrips
      << " full_record_bytes=" << sizeof(sc::engine::ChoiceRecord<double>)
      << " compact_record_bytes=" << sizeof(sc::engine::CompactChoiceMap<double>::Record) << '\n';
}

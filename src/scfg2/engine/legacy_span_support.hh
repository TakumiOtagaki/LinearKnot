#ifndef CPARTY_SCFG2_ENGINE_LEGACY_SPAN_SUPPORT_HH
#define CPARTY_SCFG2_ENGINE_LEGACY_SPAN_SUPPORT_HH

#include "scfg2/deductions/deduction.hh"

#include <stdexcept>
#include <string>
#include <tuple>

namespace cparty::scfg2::engine {

inline constexpr int legacy_span_phase(grammar::NonTerminal nt) {
  switch (nt) {
    case grammar::NonTerminal::VM: return 0;
    case grammar::NonTerminal::V: return 1;
    case grammar::NonTerminal::WMv: return 2;
    case grammar::NonTerminal::VP: return 3;
    case grammar::NonTerminal::BE: return 4;
    case grammar::NonTerminal::VPL: return 5;
    case grammar::NonTerminal::VPR: return 6;
    case grammar::NonTerminal::WMBP: return 7;
    case grammar::NonTerminal::WMBW: return 8;
    case grammar::NonTerminal::WMB: return 9;
    case grammar::NonTerminal::WMp: return 10;
    case grammar::NonTerminal::WM: return 11;
    case grammar::NonTerminal::WIP: return 12;
    case grammar::NonTerminal::WI: return 13;
    case grammar::NonTerminal::W: return 14;
  }
  return 15;
}

inline auto legacy_span_schedule_rank(const deductions::ItemKey &item) {
  return std::make_tuple(item.j - item.i + 1,
                         legacy_span_phase(item.nonterminal),
                         item.i,
                         item.ip,
                         item.jp);
}

inline deductions::ItemKey runtime_child_item_for_support(
    const deductions::Deduction &deduction,
    std::size_t child_idx) {
  deductions::ItemKey child = deduction.children[child_idx].item;
  return child;
}

inline bool legacy_span_schedule_supports(
    const deductions::Deduction &deduction) {
  if (deduction.child_count > deduction.children.size()) return false;
  const auto &parent = deduction.parent;
  const auto parent_span = parent.j - parent.i + 1;
  int parent_phase = -1;
  for (std::size_t child_idx = 0;
       child_idx < deduction.child_count;
       ++child_idx) {
    const auto child = runtime_child_item_for_support(deduction, child_idx);
    const auto child_span = child.j - child.i + 1;
    // The common shorter-span edge needs no phase/key comparison. The
    // remaining branches are exactly the rank tuple's lexicographic order;
    // equal ranks are still rejected, including on unknown nonterminals.
    if (child_span < parent_span) continue;
    if (child_span > parent_span) return false;
    if (parent_phase < 0) parent_phase = legacy_span_phase(parent.nonterminal);
    const int child_phase = legacy_span_phase(child.nonterminal);
    if (child_phase < parent_phase) continue;
    if (child_phase > parent_phase) return false;
    if (child.i < parent.i) continue;
    if (child.i > parent.i) return false;
    if (child.ip < parent.ip) continue;
    if (child.ip > parent.ip) return false;
    if (child.jp >= parent.jp) return false;
  }
  return true;
}

// The canonical span rank is a static dependency invariant, not a grammar
// admission policy or an observation of the actual executor's commit order.
// Keep the boolean predicate above for audits. Runtime consumers must report
// a generator/order bug instead of silently deleting its contribution.
// Keep allocation-heavy diagnostic formatting off the valid-edge hot path.
[[gnu::cold, gnu::noinline]] inline void report_legacy_span_schedule_violation(
    const deductions::Deduction &deduction) {
  if (deduction.child_count > deduction.children.size()) {
    throw std::logic_error("internal scfg2 schedule invariant: invalid child count");
  }
  const auto parent_rank = legacy_span_schedule_rank(deduction.parent);
  for (std::size_t c = 0; c < deduction.child_count; ++c) {
    const auto child = runtime_child_item_for_support(deduction, c);
    if (legacy_span_schedule_rank(child) < parent_rank) continue;
    const auto describe = [](const deductions::ItemKey &key) {
      return std::string(grammar::nonterminal_name(key.nonterminal)) + "(" +
          std::to_string(key.i) + "," + std::to_string(key.j) + ";" +
          std::to_string(key.ip) + "," + std::to_string(key.jp) + ")" +
          " span=" + std::to_string(key.j - key.i + 1) +
          " phase=" + std::to_string(legacy_span_phase(key.nonterminal));
    };
    throw std::logic_error(
        "internal scfg2 schedule invariant violated: rule=" +
        std::string(grammar::rule_name(deduction.rule)) +
        " parent=" + describe(deduction.parent) +
        " runtime_child[" + std::to_string(c) + "]=" + describe(child) +
        " declared_child=" + describe(deduction.children[c].item) +
        " split=(k=" + std::to_string(deduction.split.k) +
        ",l=" + std::to_string(deduction.split.l) +
        ",p=" + std::to_string(deduction.split.p) +
        ",q=" + std::to_string(deduction.split.q) +
        "); expected runtime child rank < parent rank");
  }
}

inline void require_legacy_span_schedule_support(
    const deductions::Deduction &deduction) {
  if (!legacy_span_schedule_supports(deduction))
    report_legacy_span_schedule_violation(deduction);
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_LEGACY_SPAN_SUPPORT_HH

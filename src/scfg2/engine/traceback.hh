#ifndef CPARTY_SCFG2_ENGINE_TRACEBACK_HH
#define CPARTY_SCFG2_ENGINE_TRACEBACK_HH

#include "scfg2/engine/chart.hh"
#include "scfg2/engine/choice_map.hh"

#include <optional>
#include <utility>
#include <vector>

namespace cparty::scfg2::engine {

template <class Value>
struct DeductionTraceNode;

// Keep the node's aggregate/vector interface while making ownership of a
// deep derivation independent of the process stack. The teardown link is
// used only during clear(); ordinary copies and moves start with no link.
template <class Value>
class TraceChildren : public std::vector<DeductionTraceNode<Value>> {
 public:
  using Node = DeductionTraceNode<Value>;
  using Base = std::vector<Node>;
  using Base::Base;

  TraceChildren() = default;
  TraceChildren(const TraceChildren &other);
  TraceChildren(TraceChildren &&other) noexcept;
  TraceChildren &operator=(const TraceChildren &other);
  TraceChildren &operator=(TraceChildren &&other) noexcept;
  ~TraceChildren() noexcept;
  void clear() noexcept;

 private:
  void copy_from(const Base &other);
  TraceChildren *teardown_parent_ = nullptr;
};

template <class Value>
struct DeductionTraceNode {
  deductions::ItemKey item;
  Value chart_value{};
  std::optional<deductions::Deduction> chosen_deduction;
  Value chosen_contribution{};
  TraceChildren<Value> children;
};

template <class Value>
TraceChildren<Value>::TraceChildren(const TraceChildren &other) : Base() {
  copy_from(other);
}

template <class Value>
TraceChildren<Value>::TraceChildren(TraceChildren &&other) noexcept
    : Base(std::move(other)) {}

template <class Value>
TraceChildren<Value> &TraceChildren<Value>::operator=(const TraceChildren &other) {
  if (this != &other) {
    TraceChildren copy(other);
    Base::swap(copy);
  }
  return *this;
}

template <class Value>
TraceChildren<Value> &TraceChildren<Value>::operator=(TraceChildren &&other) noexcept {
  if (this != &other) {
    // Detach first so assignment from a descendant stays valid while the
    // previous tree is destroyed, without relying on vector's release order.
    Base detached(std::move(static_cast<Base &>(other)));
    clear();
    Base::operator=(std::move(detached));
  }
  return *this;
}

template <class Value>
TraceChildren<Value>::~TraceChildren() noexcept {
  clear();
}

template <class Value>
void TraceChildren<Value>::clear() noexcept {
  // Descend through the last child, then destroy it only after its own
  // children are empty. No worklist allocation is needed in a destructor.
  if (Base::empty()) return;
  auto *current = this;
  while (true) {
    if (current->empty()) {
      auto *parent = current->teardown_parent_;
      current->teardown_parent_ = nullptr;
      if (!parent) break;
      current = parent;
    } else if (!current->back().children.empty()) {
      auto *child = &current->back().children;
      child->teardown_parent_ = current;
      current = child;
    } else {
      // The child's destructor sees an empty children vector and returns
      // immediately, rather than recursing down the derivation.
      static_cast<Base &>(*current).pop_back();
    }
  }
}

template <class Value>
void TraceChildren<Value>::copy_from(const Base &other) {
  Base::resize(other.size());
  std::vector<std::pair<const Node *, Node *>> pending;
  for (std::size_t i = 0; i < other.size(); ++i)
    pending.emplace_back(&other[i], &(*this)[i]);
  while (!pending.empty()) {
    const auto [source, target] = pending.back();
    pending.pop_back();
    target->item = source->item;
    target->chart_value = source->chart_value;
    target->chosen_deduction = source->chosen_deduction;
    target->chosen_contribution = source->chosen_contribution;
    static_cast<Base &>(target->children).resize(source->children.size());
    for (std::size_t i = 0; i < source->children.size(); ++i)
      pending.emplace_back(&source->children[i], &target->children[i]);
  }
}

template <class Value>
DeductionTraceNode<Value> build_traceback_tree(
    const deductions::ItemKey &root,
    const Chart<Value> &chart,
    const ChoiceMap<Value> &choices) {
  using Node = DeductionTraceNode<Value>;
  Node tree;
  std::vector<std::pair<Node *, deductions::ItemKey>> pending{{&tree, root}};
  while (!pending.empty()) {
    const auto [node, item] = pending.back();
    pending.pop_back();
    node->item = item;
    node->chart_value = chart.get(item);
    const auto *choice = choices.find(item);
    if (!choice || !choice->has_value) continue;
    node->chosen_deduction = choice->deduction;
    node->chosen_contribution = choice->contribution;
    node->children.resize(choice->deduction.child_count);
    for (std::size_t idx = choice->deduction.child_count; idx > 0; --idx) {
      pending.emplace_back(&node->children[idx - 1],
                           choice->deduction.children[idx - 1].item);
    }
  }
  return tree;
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_TRACEBACK_HH

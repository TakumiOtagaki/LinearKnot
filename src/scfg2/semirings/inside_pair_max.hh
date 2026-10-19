#ifndef CPARTY_SCFG2_SEMIRINGS_INSIDE_PAIR_MAX_HH
#define CPARTY_SCFG2_SEMIRINGS_INSIDE_PAIR_MAX_HH

#include "scfg2/deductions/item_key.hh"
#include "scfg2/semirings/semiring.hh"

#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

namespace cparty::scfg2::semirings {

struct PairMaxStorage {
  std::size_t nodes = 0, factors = 0, products = 0, product_factors = 0;
  std::size_t peak_nodes = 0, peak_factors = 0, peak_products = 0;
  std::size_t peak_product_factors = 0, max_arity = 0;
  std::size_t sums = 0, multiplications = 0;
};

struct PairMaxNode {
  double score;
  std::optional<deductions::ItemKey> feature;
  std::vector<std::shared_ptr<const PairMaxNode>> factors;
  std::shared_ptr<PairMaxStorage> storage;
  mutable const PairMaxNode *pending_delete = nullptr;

  // Last-reference release is iterative even for a very deep selected DAG.
  // The intrusive queue needs no allocation during destruction.
  static void destroy(const PairMaxNode *node) noexcept {
    static thread_local const PairMaxNode *pending = nullptr;
    static thread_local bool draining = false;
    node->pending_delete = pending;
    pending = node;
    if (draining) return;
    draining = true;
    while (pending) {
      const auto *next = pending;
      pending = next->pending_delete;
      --next->storage->nodes;
      next->storage->factors -= next->factors.size();
      delete next;
    }
    draining = false;
  }
};

struct PairMaxProduct {
  double reward = 0.0;
  std::optional<deductions::ItemKey> feature;
  std::vector<std::shared_ptr<const PairMaxNode>> factors;
  std::shared_ptr<PairMaxStorage> storage;
  mutable std::weak_ptr<const PairMaxNode> sealed;
  bool sum = false;
  std::size_t counted_factors = 0;
  ~PairMaxProduct() {
    if (storage) {
      --storage->products;
      storage->product_factors -= counted_factors;
    }
  }
};

struct InsidePairMaxValue {
  double log_value = -std::numeric_limits<double>::infinity();
  double score = -std::numeric_limits<double>::infinity();
  std::shared_ptr<const PairMaxNode> node;
  std::shared_ptr<PairMaxProduct> product;

  bool discard_chart_entry() const { return LogSumExpSemiring::is_zero(log_value); }
  InsidePairMaxValue checkpointed() const {
    if (node || !product) return *this;
    if (!product->feature && product->factors.empty()) return {log_value, score, {}, {}};
    if (!product->feature && product->factors.size() == 1)
      return {log_value, score, product->factors.front(), {}};
    auto selected = product->sealed.lock();
    if (!selected) {
      auto &s = *product->storage;
      auto *raw = new PairMaxNode{score, product->feature, product->factors,
                                  product->storage};
      ++s.nodes;
      s.factors += raw->factors.size();
      s.peak_nodes = std::max(s.peak_nodes, s.nodes);
      s.peak_factors = std::max(s.peak_factors, s.factors);
      selected = std::shared_ptr<const PairMaxNode>(raw, PairMaxNode::destroy);
      product->sealed = selected;
    }
    return {log_value, score, std::move(selected), {}};
  }
};

inline bool operator<(const InsidePairMaxValue &a, const InsidePairMaxValue &b) {
  return a.log_value < b.log_value;
}

// Alpha is the unchanged inside semiring. Only the selected posterior trace
// is retained; no list of losing alternatives or permanent allocation pool.
// All scheduler decisions, including representative choices, use alpha.
struct InsidePairMaxSemiring {
  using value_type = InsidePairMaxValue;
  class Scope {
   public:
    explicit Scope(std::shared_ptr<PairMaxStorage> storage) : previous_(active_) {
      if (previous_) throw std::logic_error("nested pair-max scopes are unsupported");
      active_ = std::move(storage);
    }
    ~Scope() { active_ = std::move(previous_); }
    Scope(const Scope &) = delete;
    Scope &operator=(const Scope &) = delete;
   private:
    std::shared_ptr<PairMaxStorage> previous_;
  };

  static value_type zero() { return {}; }
  static value_type one() { return from_log_weight(0.0); }
  static value_type from_weight(double w) { return from_log_weight(std::log(w)); }
  static value_type from_weight(weights::LogWeight w) { return from_log_weight(w.value); }
  static value_type from_log_weight(double w) {
    return LogSumExpSemiring::is_zero(w) ? zero() : value_type{w, 0.0, {}, {}};
  }
  static bool is_zero(const value_type &v) { return v.discard_chart_entry(); }
  static double priority(const value_type &v) { return v.log_value; }
  static double traceback_priority(const value_type &v) { return v.log_value; }
  static bool better(const value_type &a, const value_type &b) {
    return a.log_value > b.log_value;
  }

  static value_type times(const value_type &a, const value_type &b) {
    if (is_zero(a) || is_zero(b)) return zero();
    ++stats().multiplications;
    auto p = product();
    append(*p, a);
    append(*p, b);
    return finish(a.log_value + b.log_value, std::move(p));
  }

  static value_type plus(const value_type &a, const value_type &b) {
    if (is_zero(a)) return b;
    if (is_zero(b)) return a;
    ++stats().sums;
    // Strict improvement preserves the first forward alternative on ties.
    const auto &best = b.score > a.score ? b : a;
    auto p = product();
    if (best.node || (best.product && best.product->sum &&
        (&best == &b || best.product.use_count() != 1))) {
      append(*p, best.checkpointed());
    } else if (best.product) {
      p->reward = best.product->reward;
      p->feature = best.product->feature;
      p->factors = best.product->factors;
    }
    p->sum = true;
    return finish(LogSumExpSemiring::plus(a.log_value, b.log_value), std::move(p));
  }

  // Called by the executor's deduction-value transform, before parent sum.
  // Reward is independent of the local thermodynamic weight and of alpha.
  static value_type with_pair_reward(const value_type &v,
                                     const deductions::ItemKey &item, double reward) {
    if (is_zero(v)) return v;
    auto p = product();
    append(*p, v);
    p->feature = item;
    p->reward = reward;
    return finish(v.log_value, std::move(p));
  }

 private:
  static PairMaxStorage &stats() {
    if (!active_) throw std::logic_error("pair-max operation requires a scope");
    return *active_;
  }
  static std::shared_ptr<PairMaxProduct> product() {
    auto &s = stats();
    auto p = std::make_shared<PairMaxProduct>();
    p->storage = active_;
    ++s.products;
    s.peak_products = std::max(s.peak_products, s.products);
    return p;
  }
  static void append(PairMaxProduct &p, const value_type &v) {
    if (v.node || (v.product && (v.product->sum || v.product->feature))) {
      const auto sealed = v.checkpointed();
      if (sealed.node) p.factors.push_back(sealed.node);
    } else if (v.product) {
      p.factors.insert(p.factors.end(), v.product->factors.begin(), v.product->factors.end());
    }
  }
  static value_type finish(double alpha, std::shared_ptr<PairMaxProduct> p) {
    auto &s = *p->storage;
    s.product_factors += p->factors.size();
    p->counted_factors = p->factors.size();
    s.peak_product_factors = std::max(s.peak_product_factors, s.product_factors);
    s.max_arity = std::max(s.max_arity, p->factors.size());
    double score = p->reward;
    for (const auto &f : p->factors) score += f->score;
    return {alpha, score, {}, std::move(p)};
  }
  inline static thread_local std::shared_ptr<PairMaxStorage> active_;
};

}  // namespace cparty::scfg2::semirings
#endif

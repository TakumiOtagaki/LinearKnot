#ifndef CPARTY_SCFG2_REPLAY_LAZY_KBEST_HH
#define CPARTY_SCFG2_REPLAY_LAZY_KBEST_HH

#include "scfg2/semirings/factorized_logsumexp.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::replay {

// Huang/Chiang lazy k-best on a frozen acyclic AND/OR forest. Scores are
// thermodynamic log weights (larger is better), not inside log masses.
// Factor occurrences are ordered: repeated child IDs remain separate tails.
class LazyKBest {
 public:
  struct Derivation {
    double log_weight;
    std::uint32_t term;
    std::vector<std::size_t> child_ranks;  // zero-based
  };

  explicit LazyKBest(const semirings::FactorizedInsideGraph &graph)
      : graph_(graph), best_(graph.states().size(), negative_infinity()),
        best_term_(graph.states().size(), semirings::invalid_outside_state),
        memo_(graph.states().size()) {
    for (std::size_t node = 0; node < graph.states().size(); ++node) {
      const auto &state = graph.states()[node];
      for (std::size_t t = state.first_term;
           t < std::size_t(state.first_term) + state.term_count; ++t) {
        const auto &term = graph.terms()[t];
        double score = term.local_log_weight;
        for (std::size_t f = 0; f < term.factor_count; ++f) {
          const auto child = graph.factors()[term.first_factor + f];
          if (child >= node) throw std::logic_error("k-best forest is not topological");
          score += best_[child];
        }
        if (std::isnan(score) || score == std::numeric_limits<double>::infinity())
          throw std::logic_error("k-best forest has a nonfinite local score");
        if (score > best_[node]) {
          best_[node] = score;
          best_term_[node] = static_cast<std::uint32_t>(t);
        }
      }
    }
  }

  double best_log_weight(std::uint32_t node) const { return best_.at(node); }
  std::size_t candidates_created() const { return candidates_created_; }
  std::size_t derivations_cached() const { return derivations_cached_; }

  // Explicit dependency stack avoids C++ recursion even for long unary
  // chains. A node expands its previous winner exactly once before the next
  // heap pop; child next ranks are themselves requested only as needed.
  bool ensure(std::uint32_t node, std::size_t rank) {
    if (node >= memo_.size()) throw std::out_of_range("invalid k-best root");
    if (memo_[node].ranked.size() > rank) return true;
    if (memo_[node].exhausted) return false;
    struct Request { std::uint32_t node; std::size_t rank; };
    std::vector<Request> pending{{node, rank}};
    while (!pending.empty()) {
      const auto request = pending.back();
      auto &memo = memo_[request.node];
      if (memo.ranked.empty() && !memo.exhausted) {
        if (best_term_[request.node] == semirings::invalid_outside_state) {
          memo.exhausted = true;
        } else {
          const auto t = best_term_[request.node];
          memo.ranked.push_back({best_[request.node], t,
                                 std::vector<std::size_t>(graph_.terms()[t].factor_count, 0)});
          ++derivations_cached_;
        }
      }
      if (memo.ranked.size() > request.rank || memo.exhausted) {
        pending.pop_back();
        continue;
      }
      if (!memo.initialized) initialize_candidates(request.node);
      const auto &last = memo.ranked.back();
      const auto &term = graph_.terms()[last.term];
      if (memo.next_factor < term.factor_count) {
        const auto f = memo.next_factor;
        const auto child = graph_.factors()[term.first_factor + f];
        const auto next_rank = last.child_ranks[f] + 1;
        const auto &child_memo = memo_[child];
        if (child_memo.ranked.size() <= next_rank && !child_memo.exhausted) {
          pending.push_back({child, next_rank});
          continue;
        }
        if (child_memo.ranked.size() > next_rank) {
          auto ranks = last.child_ranks;
          ranks[f] = next_rank;
          enqueue(request.node, last.term, std::move(ranks));
        }
        ++memo.next_factor;
        continue;
      }
      if (memo.heap.empty()) {
        memo.exhausted = true;
      } else {
        std::pop_heap(memo.heap.begin(), memo.heap.end(), Worse{});
        memo.ranked.push_back(std::move(memo.heap.back()));
        memo.heap.pop_back();
        memo.next_factor = 0;
        ++derivations_cached_;
      }
    }
    return memo_[node].ranked.size() > rank;
  }

  // Reference remains valid until another ensure() extends this node's list.
  const Derivation &derivation(std::uint32_t node, std::size_t rank) {
    if (!ensure(node, rank)) throw std::out_of_range("k-best rank is exhausted");
    return memo_[node].ranked[rank];
  }

  template <class Visit>
  void visit_derivation(std::uint32_t root, std::size_t rank, Visit &&visit) {
    struct Occurrence { std::uint32_t node; std::size_t rank; };
    std::vector<Occurrence> pending{{root, rank}};
    while (!pending.empty()) {
      const auto occurrence = pending.back();
      pending.pop_back();
      const auto &d = derivation(occurrence.node, occurrence.rank);
      visit(occurrence.node);
      const auto &term = graph_.terms()[d.term];
      for (std::size_t f = term.factor_count; f > 0; --f)
        pending.push_back({graph_.factors()[term.first_factor + f - 1], d.child_ranks[f - 1]});
    }
  }

 private:
  struct Key {
    std::uint32_t term;
    std::vector<std::size_t> ranks;
    bool operator==(const Key &other) const { return term == other.term && ranks == other.ranks; }
  };
  struct KeyHash {
    std::size_t operator()(const Key &key) const {
      std::size_t hash = key.term;
      for (const auto rank : key.ranks)
        hash ^= rank + 0x9e3779b9U + (hash << 6) + (hash >> 2);
      return hash;
    }
  };
  struct Worse {
    bool operator()(const Derivation &a, const Derivation &b) const {
      if (a.log_weight != b.log_weight) return a.log_weight < b.log_weight;
      if (a.term != b.term) return a.term > b.term;
      return a.child_ranks > b.child_ranks;
    }
  };
  struct Memo {
    std::vector<Derivation> ranked;
    std::vector<Derivation> heap;
    std::unordered_set<Key, KeyHash> seen;
    std::size_t next_factor = 0;
    bool initialized = false;
    bool exhausted = false;
  };
  static double negative_infinity() { return -std::numeric_limits<double>::infinity(); }

  void enqueue(std::uint32_t node, std::uint32_t t, std::vector<std::size_t> ranks,
                bool update_heap = true) {
    auto &memo = memo_[node];
    if (!memo.seen.insert({t, ranks}).second) return;
    const auto &term = graph_.terms()[t];
    double score = term.local_log_weight;
    for (std::size_t f = 0; f < term.factor_count; ++f) {
      const auto child = graph_.factors()[term.first_factor + f];
      score += ranks[f] == 0 ? best_[child] : memo_[child].ranked[ranks[f]].log_weight;
    }
    if (score == negative_infinity()) return;
    memo.heap.push_back({score, t, std::move(ranks)});
    if (update_heap) std::push_heap(memo.heap.begin(), memo.heap.end(), Worse{});
    ++candidates_created_;
  }

  void initialize_candidates(std::uint32_t node) {
    auto &memo = memo_[node];
    memo.initialized = true;
    const auto &state = graph_.states()[node];
    for (std::size_t t = state.first_term;
         t < std::size_t(state.first_term) + state.term_count; ++t) {
      std::vector<std::size_t> ranks(graph_.terms()[t].factor_count, 0);
      if (t == best_term_[node]) {
        memo.seen.insert({static_cast<std::uint32_t>(t), std::move(ranks)});
      } else {
        enqueue(node, static_cast<std::uint32_t>(t), std::move(ranks), false);
      }
    }
    std::make_heap(memo.heap.begin(), memo.heap.end(), Worse{});
  }

  const semirings::FactorizedInsideGraph &graph_;
  std::vector<double> best_;
  std::vector<std::uint32_t> best_term_;
  std::vector<Memo> memo_;
  std::size_t candidates_created_ = 0;
  std::size_t derivations_cached_ = 0;
};

}  // namespace cparty::scfg2::replay
#endif

#ifndef CPARTY_SCFG2_BEAM_WMBP_WI_PREFIX_ADMISSION_HH
#define CPARTY_SCFG2_BEAM_WMBP_WI_PREFIX_ADMISSION_HH

#include "scfg2/beam/wmbp_wi_prefix_store.hh"
#include "scfg2/beam/public_item_ranking.hh"
#include "scfg2/beam/wmbp_be_vp_prefix_profile.hh"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

// A temporal top-b set is reconstructed by comparing a candidate with its
// worst admitted rank at the original query time. Save one boundary per
// retained BE and per queried VP, never the b selected prefixes per query.
// Child alphas are immutable after commit. BE commit order resolves the
// forward bucket's FIFO ties (one prefix per BE in each VP-start bucket).
class WmbpWiPrefixAdmission {
 public:
  // A single chronological pass over immutable retained alphas reproduces
  // the forward queries. No DP value is recomputed and no inside history is
  // required. BE commit order is the public frontier's priority/key order.
  template <class Semiring, class Chart>
  static WmbpWiPrefixAdmission from_retained_chart(
      cand_pos_t n, const Chart &chart,
      std::size_t history_width, std::size_t prefix_width,
      const PublicItemRanking &ranking = {}) {
    std::vector<std::vector<deductions::ItemKey>> be(n + 1);
    std::vector<std::vector<cand_pos_t>> wi(n + 1), vp(n + 1);
    for (const auto &[item, value] : chart.values()) {
      if (Semiring::is_zero(value)) continue;
      if (item.nonterminal == grammar::NonTerminal::BE) be[item.j].push_back(item);
      if (item.nonterminal == grammar::NonTerminal::WI) wi[item.j].push_back(item.i);
      if (item.nonterminal == grammar::NonTerminal::VP) vp[item.j].push_back(item.i);
    }
    WmbpWiPrefixAdmission out;
    WmbpWiPrefixStore<Semiring> store;
    for (cand_pos_t j = 1; j <= n; ++j) {
      std::sort(be[j].begin(), be[j].end(), [&](const auto &a, const auto &b) {
        const auto priority = [&](const auto &key) {
          const double raw = Semiring::priority(chart.get(key));
          return ranking ? ranking(key, raw, raw) : raw;
        };
        const double x = priority(a);
        const double y = priority(b);
        return x != y ? x > y : item_tie_key(a) < item_tie_key(b);
      });
      store.commit_be_column(be[j], chart, prefix_width);
      out.template observe_be_column<Semiring>(be[j], chart, store);
      out.template observe_vp_column<Semiring>(j, vp[j], store);
      // WI(j) commits after BE/VP/WMBP(j) in the shared online schedule.
      store.commit_wi_column(j, wi[j], chart, history_width);
    }
    return out;
  }

  template <class Semiring, class Chart>
  void observe_be_column(const std::vector<deductions::ItemKey> &items,
                         const Chart &chart,
                         const WmbpWiPrefixStore<Semiring> &store) {
    for (const auto &be : items) {
      if (Semiring::is_zero(chart.get(be))) continue;
      const auto boundary = store.wi_history_boundary(be.ip + 1);
      const auto [it, inserted] = be_.emplace(
          be, BeAdmission{next_be_order_++, boundary});
      (void)it;
      if (!inserted)
        throw std::logic_error("WMBP admission observed a BE twice");
    }
  }

  template <class Semiring>
  void observe_vp_column(cand_pos_t j,
                         const std::vector<cand_pos_t> &starts,
                         const WmbpWiPrefixStore<Semiring> &store) {
    for (const auto i : starts) {
      std::optional<PrefixBoundary> boundary;
      if (const auto *last = store.prefix_boundary(i)) {
        const auto be = be_.find(last->be);
        if (be == be_.end())
          throw std::logic_error("WMBP admission is missing a committed BE");
        boundary = PrefixBoundary{
            Semiring::priority(last->factor), last->be.i, last->be.j,
            be->second.order};
      }
      const auto [it, inserted] = vp_.emplace(
          deductions::ItemKey{grammar::NonTerminal::VP, i, j, -1, -1},
          boundary);
      (void)it;
      if (!inserted)
        throw std::logic_error("WMBP admission observed a VP query twice");
    }
  }

  template <class Semiring, class Chart>
  bool allows(const deductions::Deduction &source, const Chart &chart) const {
    if (source.rule != grammar::RuleId::WMBP_SPLIT_BE_WI_VP) return true;
    const auto &core = source.children[0].item;
    const auto &right = wmbp_be_vp_right_child(source);
    const auto be = be_.find(core);
    const auto vp = vp_.find(right);
    if (be == be_.end() || vp == vp_.end() || !vp->second ||
        core.j > source.parent.j) return false;
    auto factor = chart.get(core);
    if (Semiring::is_zero(factor)) return false;
    if (source.child_count == 3) {
      const auto &wi = source.children[1].item;
      const auto value = chart.get(wi);
      if (Semiring::is_zero(value) || wi.j + 1 >= core.j ||
          !be->second.wi) return false;
      const auto &boundary = *be->second.wi;
      const double priority = Semiring::priority(value);
      if (priority < boundary.priority ||
          (priority == boundary.priority && wi.j > boundary.right))
        return false;
      factor = Semiring::times(factor, value);
    }
    const auto &boundary = *vp->second;
    const double priority = Semiring::priority(factor);
    if (priority != boundary.priority) return priority > boundary.priority;
    if (core.i != boundary.i) return core.i < boundary.i;
    if (core.j != boundary.j) return core.j < boundary.j;
    return be->second.order <= boundary.order;
  }

  std::size_t be_boundary_count() const { return be_.size(); }
  std::size_t vp_boundary_count() const { return vp_.size(); }

 private:
  struct BeAdmission {
    std::size_t order;
    std::optional<WmbpWiHistoryBoundary> wi;
  };
  struct PrefixBoundary {
    double priority;
    cand_pos_t i, j;
    std::size_t order;
  };
  std::unordered_map<deductions::ItemKey, BeAdmission,
                     deductions::ItemKeyHash, deductions::ItemKeyEq> be_;
  std::unordered_map<deductions::ItemKey, std::optional<PrefixBoundary>,
                     deductions::ItemKeyHash, deductions::ItemKeyEq> vp_;
  std::size_t next_be_order_ = 0;
};

}  // namespace cparty::scfg2::beam
#endif

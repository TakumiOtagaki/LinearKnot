#ifndef CPARTY_SCFG2_BEAM_WMB_REGENERATED_BTAIL_CELLS_HH
#define CPARTY_SCFG2_BEAM_WMB_REGENERATED_BTAIL_CELLS_HH

#include "scfg2/beam/wmb_retained_agenda.hh"
#include "scfg2/engine/chart.hh"

namespace cparty::scfg2::beam {

// Rebuild selected BTAIL alphas with the shared seed generator, storing only
// b private cells per closing column. Public WMB split terms and private seed
// products are regenerated and differentiated through that same surface.
template <class Semiring>
class WmbRegeneratedBtailCells {
 public:
  using Value = typename Semiring::value_type;
  using Outside = std::vector<std::unordered_map<cand_pos_t, Value>>;

  template <class Chart, class Context, class StartsFn>
  static WmbRegeneratedBtailCells from_retained_chart(
      cand_pos_t n, const Chart &inside, const Context &ctx,
      StartsFn &&starts, std::size_t width) {
    WmbRegeneratedBtailCells out;
    out.n_ = n;
    for (cand_pos_t t = 2; t <= n; ++t)
      out.store_.build(t, inside, ctx, starts(grammar::NonTerminal::WI, t - 1),
          [&](cand_pos_t k) -> const auto & { return starts(grammar::NonTerminal::WMBP, k); }, width);
    return out;
  }

  const WmbBtailStore<Semiring> &store() const { return store_; }
  Outside empty_outside() const { return Outside(n_ + 1); }
  std::size_t size() const { return store_.stats().selected; }

  // Reverse the shared WMB aggregate: beta_BE += beta_WMB * w * BTAIL,
  // beta_BTAIL += beta_WMB * w * BE. The direct WMBP rule is a separate
  // original public deduction and belongs to the generic outside consumer.
  template <class Chart, class Context>
  void propagate_wmb_splits(cand_pos_t j, const Chart &inside, const Context &ctx,
                            const std::vector<deductions::ItemKey> &be_items,
                            engine::Chart<Value> *public_outside, Outside *outside) const {
    for_each_wmb_be_btail_term(j, ctx, be_items, store_,
        [&](const auto &source, const auto &be, const auto &tail, auto local) {
          if (Semiring::is_zero(inside.get(source.parent))) return;
          const auto beta = public_outside->get(source.parent);
          const auto weighted = Semiring::times(beta, Semiring::from_weight(local));
          const auto to_be = Semiring::times(weighted, tail.total);
          if (!Semiring::is_zero(to_be)) public_outside->set(be,
              Semiring::plus(public_outside->get(be), to_be));
          add(&(*outside)[be.jp], tail.i, Semiring::times(weighted, inside.get(be)));
        });
  }

  // After all WMB consumers at t (and later columns), reverse BTAIL(i,t)
  // into retained WMBP / WI children. No selected branch lists are stored.
  template <class Chart, class Context, class StartsFn>
  void propagate_column(cand_pos_t t, const Chart &inside, const Context &ctx,
                         StartsFn &&starts, Outside *outside,
                         engine::Chart<Value> *public_outside) const {
    if (t < 2 || (*outside)[t].empty()) return;
    for_each_wmb_btail_seed(t, ctx, starts(grammar::NonTerminal::WI, t - 1),
        [&](cand_pos_t k) -> const auto & { return starts(grammar::NonTerminal::WMBP, k); },
        [&](cand_pos_t i, cand_pos_t k, bool has_wi) {
          const auto found = (*outside)[t].find(i);
          if (found == (*outside)[t].end()) return;
          const deductions::ItemKey left{grammar::NonTerminal::WMBP, i, k, -1, -1};
          const deductions::ItemKey right{grammar::NonTerminal::WI, k + 1, t - 1, -1, -1};
          const auto to_left = Semiring::times(found->second,
              has_wi ? inside.get(right) : Semiring::one());
          if (!Semiring::is_zero(to_left)) public_outside->set(left,
              Semiring::plus(public_outside->get(left), to_left));
          if (has_wi) {
            const auto to_right = Semiring::times(found->second, inside.get(left));
            if (!Semiring::is_zero(to_right)) public_outside->set(right,
                Semiring::plus(public_outside->get(right), to_right));
          }
        });
  }

 private:
  static void add(std::unordered_map<cand_pos_t, Value> *column,
                    cand_pos_t i, Value value) {
    if (Semiring::is_zero(value)) return;
    auto [found, inserted] = column->try_emplace(i, Semiring::zero());
    (void)inserted;
    found->second = Semiring::plus(found->second, value);
  }
  cand_pos_t n_ = 0;
  WmbBtailStore<Semiring> store_;
};

}  // namespace cparty::scfg2::beam
#endif

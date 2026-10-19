#ifndef CPARTY_SCFG2_BEAM_VPR_REGENERATED_BASE_CELLS_HH
#define CPARTY_SCFG2_BEAM_VPR_REGENERATED_BASE_CELLS_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/vpr_retained_agenda.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/chart.hh"
#include "scfg2/engine/legacy_span_support.hh"

namespace cparty::scfg2::beam {

// Genuine private recurrence C(i,j), distinct from the public VPR(i,j):
//   C(i,j) = cp * (kept C(i,j-1) + retained VP(i,j-1)).
// The private top-b carry can outlive a pruned public VPR. Rebuild it once
// from immutable retained VP alphas using the shared forward agenda, then
// reverse this two-term recurrence without keeping any branch graph.
template <class Semiring>
class VprRegeneratedBaseCells {
 public:
  using Value = typename Semiring::value_type;
  struct Cell { Value alpha; bool carried; bool public_base; };
  using Column = std::unordered_map<cand_pos_t, Cell>;
  using Outside = std::vector<std::unordered_map<cand_pos_t, Value>>;

  template <class Chart, class Context, class StartsFn>
  static VprRegeneratedBaseCells from_retained_chart(
      cand_pos_t n, const Chart &chart, Context &ctx,
      StartsFn &&starts, std::size_t width) {
    VprRegeneratedBaseCells out;
    out.cells_.resize(n + 1);
    struct Recorder : RunTermObserver<Value> {
      std::vector<Column> *columns;
      explicit Recorder(std::vector<Column> *p) : columns(p) {}
      void on_run_term(grammar::RuleId, cand_pos_t, cand_pos_t, const Value &) override {}
      void on_left_run_unit(grammar::RuleId, cand_pos_t, cand_pos_t,
                            const Value &, const Value *, bool) override {}
      void on_vpr_run_cell(cand_pos_t i, cand_pos_t j,
                            const Value &alpha, bool carried) override {
        (*columns)[j].emplace(i, Cell{alpha, carried, false});
      }
    } recorder(&out.cells_);
    ScopedRunTermObserver<Value> scope(&recorder);
    VprRetainedAgenda<Semiring> agenda;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    for (cand_pos_t j = 1; j <= n; ++j) {
      agenda.build_column(j, chart, ctx, starts(grammar::NonTerminal::WIP, j),
          [&](cand_pos_t k) -> const auto & { return starts(grammar::NonTerminal::VP, k); },
          width);
      for (auto &[i, cell] : out.cells_[j]) {
        const deductions::ItemKey parent{grammar::NonTerminal::VPR, i, j, -1, -1};
        if (Semiring::is_zero(chart.get(parent))) continue;
        for (const auto &source : agenda.deductions_for(parent)) {
          engine::require_legacy_span_schedule_support(source);
          if (agenda.reconstructs_base(source) &&
              !weights::is_zero(weight.local_weight_for(source))) cell.public_base = true;
        }
      }
    }
    return out;
  }

  Outside empty_outside() const { return Outside(cells_.size()); }
  const Column &cells(cand_pos_t j) const { return cells_[j]; }
  std::size_t size() const {
    std::size_t total = 0;
    for (const auto &column : cells_) total += column.size();
    return total;
  }

  void add_public_outside(const deductions::ItemKey &parent, Value beta,
                           Outside *outside) const {
    const auto found = cells_[parent.j].find(parent.i);
    if (found != cells_[parent.j].end() && found->second.public_base)
      add(&(*outside)[parent.j], parent.i, beta);
  }

  // Call in descending column order after all current public VPR beta has
  // arrived. Carry contributions from j+1 have already arrived in C(i,j).
  template <class Chart, class Context>
  void propagate_column(cand_pos_t j, const Chart &inside, const Context &ctx,
                         Outside *outside, engine::Chart<Value> *public_outside) const {
    const auto cp = Semiring::from_weight(weights::local_weight(ctx, ctx.vpr_split_vp_basepair_weight(1)));
    for (const auto &[i, cell] : cells_[j]) {
      const auto found = (*outside)[j].find(i);
      if (found == (*outside)[j].end()) continue;
      const auto beta = Semiring::times(found->second, cp);
      if (cell.carried) add(&(*outside)[j - 1], i, beta);
      const deductions::ItemKey seed{grammar::NonTerminal::VP, i, j - 1, -1, -1};
      if (i < j - 1 && !Semiring::is_zero(inside.get(seed)))
        public_outside->set(seed, Semiring::plus(public_outside->get(seed), beta));
    }
  }

 private:
  static void add(std::unordered_map<cand_pos_t, Value> *column,
                    cand_pos_t i, Value value) {
    if (Semiring::is_zero(value)) return;
    auto [found, inserted] = column->try_emplace(i, Semiring::zero());
    (void)inserted;
    found->second = Semiring::plus(found->second, value);
  }
  std::vector<Column> cells_;
};

}  // namespace cparty::scfg2::beam
#endif

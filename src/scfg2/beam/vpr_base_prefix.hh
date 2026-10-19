#ifndef CPARTY_SCFG2_BEAM_VPR_BASE_PREFIX_HH
#define CPARTY_SCFG2_BEAM_VPR_BASE_PREFIX_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/item_key.hh"

#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

template <class Value>
struct VprBasePrefixCell {
  Value total{};
  Value best_term{};
  cand_pos_t best_k = -1;
};

template <class Semiring>
class VprBasePrefixColumns {
 public:
  using Value = typename Semiring::value_type;

  template <class Chart, class Context>
  void build_column(cand_pos_t j, const Chart &chart, const Context &ctx) {
    std::vector<VprBasePrefixCell<Value>> current(
        static_cast<std::size_t>(j + 1),
        VprBasePrefixCell<Value>{Semiring::zero(), Semiring::zero(), -1});
    if (j > 1 && ctx.can_pair_right_span(j - 1, j)) {
      // Every term VP(i,k) * cp^(j-k) gains one unpaired base per column.
      const Value one_unpaired =
          Semiring::from_weight(weights::local_weight(ctx, ctx.vpr_split_vp_basepair_weight(1)));
      for (cand_pos_t i = 1; i < j; ++i) {
        auto &cell = current[static_cast<std::size_t>(i)];
        if (static_cast<std::size_t>(i) < previous_.size()) {
          cell = previous_[static_cast<std::size_t>(i)];
          cell.total = Semiring::times(cell.total, one_unpaired);
          cell.best_term = Semiring::times(cell.best_term, one_unpaired);
        }
        const cand_pos_t k = j - 1;
        if (k <= i) continue;
        const Value term = Semiring::times(
            chart.get(deductions::ItemKey{
                grammar::NonTerminal::VP, i, k, -1, -1}),
            one_unpaired);
        cell.total = Semiring::plus(cell.total, term);
        if (Semiring::better(term, cell.best_term)) {
          cell.best_term = term;
          cell.best_k = k;
        }
      }
    }
    previous_ = std::move(current);
  }

  const VprBasePrefixCell<Value> &cell(cand_pos_t i) const {
    return previous_[static_cast<std::size_t>(i)];
  }

 private:
  std::vector<VprBasePrefixCell<Value>> previous_;
};

}  // namespace cparty::scfg2::beam

#endif

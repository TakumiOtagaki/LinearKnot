#ifndef CPARTY_SCFG2_BEAM_RETAINED_CHART_VIEWS_HH
#define CPARTY_SCFG2_BEAM_RETAINED_CHART_VIEWS_HH

#include "scfg2/beam/right_boundary_item_frontier.hh"
#include "scfg2/beam/public_item_ranking.hh"
#include "scfg2/beam/vp_wmbp_online_schedule.hh"

#include <array>

namespace cparty::scfg2::beam {

// Immutable, sorted retained-child indices. Reconstruct both WMB consumer
// views from their union: every original winner is in that union, and the
// priority/key order is identical to the forward frontier.
class RetainedChartViews {
 public:
  using Key = deductions::ItemKey;
  using NT = grammar::NonTerminal;
  static constexpr std::size_t family_count = static_cast<std::size_t>(NT::BE) + 1;

  template <class S, class Chart>
  static RetainedChartViews build(cand_pos_t n, const Chart &chart, std::size_t wmb_width,
                                  const PublicItemRanking &ranking = {}) {
    RetainedChartViews out;
    for (auto &family : out.items_) family.resize(n + 1);
    for (auto &family : out.starts_) family.resize(n + 1);
    out.internal_.resize(n + 1);
    out.root_.resize(n + 1);
    for (const auto &[key, value] : chart.values()) {
      if (S::is_zero(value)) continue;
      out.items_[static_cast<std::size_t>(key.nonterminal)][key.j].push_back(key);
    }
    for (std::size_t nt = 0; nt < family_count; ++nt) {
      for (cand_pos_t j = 1; j <= n; ++j) {
        auto &items = out.items_[nt][j];
        std::sort(items.begin(), items.end(), [](const auto &a, const auto &b) {
          return right_boundary_item_tie_key(a) < right_boundary_item_tie_key(b);
        });
        auto &starts = out.starts_[nt][j];
        for (const auto &key : items) starts.push_back(key.i);
        starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
      }
    }
    for (cand_pos_t j = 1; j <= n; ++j) {
      RightBoundaryTopWidthItemFrontier internal(wmb_width), root(wmb_width);
      for (const auto &key : out.items(NT::WMB, j)) {
        const auto value = chart.get(key);
        const auto prefix = key.i == 1 ? S::one()
            : chart.get(Key{NT::W, 1, key.i - 1, -1, -1});
        internal.observe(key, S::priority(value));
        const double raw = S::priority(value);
        const double existing = S::priority(S::times(prefix, value));
        root.observe(key, ranking ? ranking(key, raw, existing) : existing);
      }
      internal.flush_column(j, [&](const auto &key) { out.internal_[j].push_back(key.i); });
      root.flush_column(j, [&](const auto &key) { out.root_[j].push_back(key.i); });
      std::sort(out.internal_[j].begin(), out.internal_[j].end());
      std::sort(out.root_[j].begin(), out.root_[j].end());
      auto &be = out.items_[static_cast<std::size_t>(NT::BE)][j];
      std::sort(be.begin(), be.end(), [&](const auto &a, const auto &b) {
        const auto priority = [&](const auto &key) {
          const double raw = S::priority(chart.get(key));
          return ranking ? ranking(key, raw, raw) : raw;
        };
        return right_boundary_item_rank_before(
            {a, priority(a)}, {b, priority(b)});
      });
    }
    return out;
  }

  const std::vector<cand_pos_t> &starts(NT nt, cand_pos_t j) const {
    if (j < 0 || static_cast<std::size_t>(j) >= starts_[0].size()) return empty_starts();
    return starts_[static_cast<std::size_t>(nt)][j];
  }
  const std::vector<cand_pos_t> &internal_starts(NT nt, cand_pos_t j) const {
    if (nt != NT::WMB) return starts(nt, j);
    return at(internal_, j);
  }
  const std::vector<cand_pos_t> &root_starts(NT nt, cand_pos_t j) const {
    if (nt != NT::WMB) return starts(nt, j);
    return at(root_, j);
  }
  const std::vector<Key> &items(NT nt, cand_pos_t j) const {
    return items_[static_cast<std::size_t>(nt)][j];
  }

 private:
  static const std::vector<cand_pos_t> &empty_starts() {
    static const std::vector<cand_pos_t> empty;
    return empty;
  }
  static const std::vector<cand_pos_t> &at(
      const std::vector<std::vector<cand_pos_t>> &columns, cand_pos_t j) {
    return j < 0 || static_cast<std::size_t>(j) >= columns.size()
        ? empty_starts() : columns[j];
  }
  std::array<std::vector<std::vector<Key>>, family_count> items_;
  std::array<std::vector<std::vector<cand_pos_t>>, family_count> starts_;
  std::vector<std::vector<cand_pos_t>> internal_, root_;
};

}  // namespace cparty::scfg2::beam
#endif

// Pilot inside and bounded completion-context ranking on the shared kernel.
#ifndef CPARTY_SCFG2_BEAM_PILOT_OUTER_HH
#define CPARTY_SCFG2_BEAM_PILOT_OUTER_HH

#include "scfg2/beam/bounded_outer.hh"
#include "scfg2/beam/ranking_config.hh"

namespace cparty::scfg2::beam::outer_detail {

inline std::size_t sqrt_pilot_width(std::size_t width) {
  if (!width || width == std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument("pilot requires a positive finite main width");
  auto p = static_cast<std::size_t>(std::sqrt(static_cast<long double>(width)));
  while (p > width / p) --p;
  while (p + 1 <= width / (p + 1)) ++p;
  return p * p == width ? p : p + 1;
}

using beam::PilotExponent;

inline const char *pilot_exponent_name(PilotExponent exponent) {
  switch (exponent) {
    case PilotExponent::Half: return "0.5";
    case PilotExponent::TwoThirds: return "2/3";
  }
  throw std::invalid_argument("invalid pilot exponent");
}

namespace pilot_width_detail {
// Fixed three-word integers suffice for p^3 and b^2 at every size_t width.
// Compare rational powers exactly, including perfect powers such as b=27.
using Power = std::array<std::uint64_t, 3>;
inline Power power(std::size_t base, unsigned exponent) {
  static_assert(std::numeric_limits<std::size_t>::digits <= 64);
  Power result{};
  result[0] = 1;
  for (unsigned k = 0; k < exponent; ++k) {
    __uint128_t carry = 0;
    for (auto &word : result) {
      const __uint128_t product = __uint128_t(word) * base + carry;
      word = static_cast<std::uint64_t>(product);
      carry = product >> 64;
    }
  }
  return result;
}
inline bool at_least(const Power &a, const Power &b) {
  for (std::size_t k = a.size(); k-- > 0;)
    if (a[k] != b[k]) return a[k] > b[k];
  return true;
}
} // namespace pilot_width_detail

inline std::size_t pilot_width_for(std::size_t width, PilotExponent exponent) {
  if (exponent == PilotExponent::Half) return sqrt_pilot_width(width);
  if (!width || width == std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument("pilot requires a positive finite main width");
  if (exponent != PilotExponent::TwoThirds)
    throw std::invalid_argument("invalid pilot exponent");
  const auto target = pilot_width_detail::power(width, 2);
  std::size_t low = 1, high = width;
  while (low < high) {
    const auto mid = low + (high - low) / 2;
    if (pilot_width_detail::at_least(pilot_width_detail::power(mid, 3), target))
      high = mid;
    else low = mid + 1;
  }
  return low;
}

struct PilotRankingCounts {
  std::array<std::size_t, outer_probe_families> covered{}, fallback{}, novel{};
  std::size_t total(const std::array<std::size_t, outer_probe_families> &counts) const {
    std::size_t result = 0;
    for (auto count : counts) result += count;
    return result;
  }
};

inline double outer_priority(double context, double alpha, double existing) {
  if (context == zero) return existing;
  if (!std::isfinite(context)) throw std::logic_error("invalid outer context");
  return alpha + context;
}

// Target pilot alpha is used only to count new discoveries, never to rank
// or admit a target. Missing guide values keep the exact existing priority.
struct PilotOuterRanking {
  const Chart &guide;
  const Chart &pilot;
  PilotRankingCounts &counts;
  double operator()(const Item &item, double alpha, double existing) const {
    const auto family = static_cast<std::size_t>(item.nonterminal);
    const double context = guide.get(item);
    if (context == zero) {
      ++counts.fallback.at(family);
      return existing;
    }
    ++counts.covered.at(family);
    if (pilot.get(item) == zero) ++counts.novel.at(family);
    return outer_priority(context, alpha, existing);
  }
};

struct PilotGuidedInsideRun {
  sc::replay::WFinalVpWmbpOnlineRun<double> inside;
  Chart pilot;
  Result guide;
  PilotRankingCounts ranks;
  std::size_t pilot_width;
  double pilot_seconds, main_seconds;
};

template <class Context>
PilotGuidedInsideRun run_pilot_guided_inside(
    Context &ctx, cand_pos_t n,
    const sc::replay::WFinalVpWmbpOnlineWidths &main_widths,
    std::size_t main_width, std::size_t outer_factor = 4,
    PilotExponent exponent = PilotExponent::TwoThirds) {
  const auto p = pilot_width_for(main_width, exponent);
  if (!outer_factor || p > std::numeric_limits<std::size_t>::max() / outer_factor)
    throw std::invalid_argument("invalid outer frontier factor");
  auto widths = sc::replay::WFinalVpWmbpOnlineWidths::uniform(p);
  widths.max_hairpin_unpaired = main_widths.max_hairpin_unpaired;
  widths.max_case3_ancestor_depth = main_widths.max_case3_ancestor_depth;
  widths.forced_items = main_widths.forced_items;
  widths.force_scaffold_derivation = main_widths.force_scaffold_derivation;
  const auto pilot_start = Clock::now();
  auto pilot = sc::replay::run_w_final_vp_wmbp_online<Sum>(ctx, n, widths);
  const double pilot_seconds = seconds(pilot_start);
  const Item root{NT::W, 1, n, -1, -1};
  const auto outer_start = Clock::now();
  auto guide = bounded_outer(pilot.chart, root, ctx, p, p * outer_factor);
  // Include destruction of the guide builder's temporary indices / pools.
  // The internal diagnostic clock stops before those scratch objects die.
  guide.outer.elapsed = seconds(outer_start);
  PilotRankingCounts counts;
  const PilotOuterRanking rank{guide.outer.values, pilot.chart, counts};
  const auto main_start = Clock::now();
  auto main = sc::replay::run_w_final_vp_wmbp_online_ranked_inside(
      ctx, n, main_widths, rank);
  const double main_seconds = seconds(main_start);
  return {std::move(main), std::move(pilot.chart), std::move(guide), counts,
          p, pilot_seconds, main_seconds};
}

} // namespace cparty::scfg2::beam::outer_detail
#endif

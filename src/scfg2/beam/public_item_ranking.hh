#ifndef CPARTY_SCFG2_BEAM_PUBLIC_ITEM_RANKING_HH
#define CPARTY_SCFG2_BEAM_PUBLIC_ITEM_RANKING_HH

#include "scfg2/deductions/item_key.hh"
#include <functional>

namespace cparty::scfg2::beam {

// Item, raw semiring priority, existing family-specific priority. This only
// chooses public states; weights and private-factor selection stay unchanged.
using PublicItemRanking = std::function<double(
    const deductions::ItemKey &, double, double)>;

} // namespace cparty::scfg2::beam
#endif

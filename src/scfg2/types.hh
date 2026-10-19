#ifndef CPARTY_SCFG2_TYPES_HH
#define CPARTY_SCFG2_TYPES_HH

#include <cstdint>

// Shared scalar vocabulary. Coordinates are signed because the grammar uses
// negative border sentinels. Public pair coordinates remain 1-based.
using cand_pos_t = std::int_least32_t;
using cand_pos_tu = std::uint_least32_t;
using pair_type = std::int_least16_t;
using pf_t = double;

#endif

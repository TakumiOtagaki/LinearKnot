#ifndef CPARTY_SCFG2_GRAMMAR_SURFACE_KIND_HH
#define CPARTY_SCFG2_GRAMMAR_SURFACE_KIND_HH

namespace cparty::scfg2::grammar {

// WMBP step-prefix labels used by exact replay observation. These are
// grammar/storage vocabulary; they do not define separate scoring semantics.
enum class WmbpStepPrefixKind : unsigned char {
  Unknown = 0,
  Empty = 1,
  PrefixWmbp = 2,
  PrefixTail = 3,
  PrefixWi = 4,
};

}  // namespace cparty::scfg2::grammar

#endif  // CPARTY_SCFG2_GRAMMAR_SURFACE_KIND_HH

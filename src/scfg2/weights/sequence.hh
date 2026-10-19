#ifndef CPARTY_SCFG2_WEIGHTS_SEQUENCE_HH
#define CPARTY_SCFG2_WEIGHTS_SEQUENCE_HH

#include <cctype>
#include <string>

namespace cparty::scfg2::weights {

// Encoded bases and sequence-based special-loop lookups must see the same
// RNA spelling. Keep the caller's original input outside the weight owner.
inline std::string normalize_rna_sequence(std::string sequence) {
  for (char &base : sequence) {
    base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    if (base == 'T') base = 'U';
  }
  return sequence;
}

}  // namespace cparty::scfg2::weights

#endif

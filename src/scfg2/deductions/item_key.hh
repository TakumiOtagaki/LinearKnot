#ifndef CPARTY_SCFG2_DEDUCTIONS_ITEM_KEY_HH
#define CPARTY_SCFG2_DEDUCTIONS_ITEM_KEY_HH

#include "scfg2/types.hh"
#include "scfg2/grammar/nonterminal.hh"

#include <cstddef>

namespace cparty::scfg2::deductions {

struct ItemKey {
  grammar::NonTerminal nonterminal = grammar::NonTerminal::W;
  cand_pos_t i = -1;
  cand_pos_t j = -1;
  cand_pos_t ip = -1;
  cand_pos_t jp = -1;
};

inline constexpr bool uses_secondary_span(grammar::NonTerminal nt) {
  return nt == grammar::NonTerminal::BE;
}

inline constexpr bool is_valid_primary_span(const ItemKey &key) {
  return key.i >= 1 && key.j >= key.i;
}

inline constexpr bool is_valid_item_key(const ItemKey &key) {
  if (!is_valid_primary_span(key)) {
    return false;
  }
  if (!uses_secondary_span(key.nonterminal)) {
    return key.ip == -1 && key.jp == -1;
  }
  return key.ip >= key.i && key.jp > key.ip && key.jp <= key.j;
}

struct ItemKeyHash {
  std::size_t operator()(const ItemKey &key) const noexcept {
    std::size_t h = static_cast<std::size_t>(key.nonterminal);
    h = (h * 1315423911u) ^ static_cast<std::size_t>(key.i + 3);
    h = (h * 1315423911u) ^ static_cast<std::size_t>(key.j + 5);
    h = (h * 1315423911u) ^ static_cast<std::size_t>(key.ip + 7);
    h = (h * 1315423911u) ^ static_cast<std::size_t>(key.jp + 11);
    return h;
  }
};

struct ItemKeyEq {
  bool operator()(const ItemKey &lhs, const ItemKey &rhs) const noexcept {
    return lhs.nonterminal == rhs.nonterminal &&
           lhs.i == rhs.i &&
           lhs.j == rhs.j &&
           lhs.ip == rhs.ip &&
           lhs.jp == rhs.jp;
  }
};

}  // namespace cparty::scfg2::deductions

#endif  // CPARTY_SCFG2_DEDUCTIONS_ITEM_KEY_HH

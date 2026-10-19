#ifndef CPARTY_SCFG2_STORAGE_EXACT_STORE_HH
#define CPARTY_SCFG2_STORAGE_EXACT_STORE_HH

#include "scfg2/types.hh"

#include <cstddef>
#include <functional>

namespace cparty::scfg2::storage {

// Historical observation-carrier identifiers; not the executor grammar.
// Explicit values preserve the recorded diagnostic identities.
enum class NonTerminal : unsigned char {
  V = 0,
  WI = 1,
  W = 2,
  VM = 3,
  WMv = 4,
  WMp = 5,
  WM = 6,
  WIP = 7,
  WIP_VP_L = 8,
  WIP_VP_R = 9,
  WIP_BE_L = 10,
  WIP_BE_R = 11,
  VPL = 12,
  VPR = 13,
  VP = 14,
  VP_CASE3_BOUNDARY = 15,
  VP_CASE3_R = 16,
  VP_CASE3_MR = 17,
  VP_CASE3_LEFT = 18,
  VP_C12_R2 = 19,
  VP_C12_2 = 20,
  VP_INTERNAL_LOCAL = 21,
  WMBW = 22,
  WMBP = 23,
  WMB_BTAIL = 24,
  WMB_BSTEP = 25,
  WMB_NONEMPTY = 26,
  WMB_EMPTY_LOCAL = 27,
  WMB = 28,
  BE = 29,
};

struct ExactItemKey {
  NonTerminal nt = NonTerminal::W;
  cand_pos_t i = -1;
  cand_pos_t j = -1;
  cand_pos_t ip = -1;
  cand_pos_t jp = -1;
  cand_pos_t aux0 = -1;
  cand_pos_t aux1 = -1;
  cand_pos_t aux2 = -1;
  cand_pos_t aux3 = -1;

  bool operator==(const ExactItemKey &other) const;
};

struct ExactItemKeyHash {
  std::size_t operator()(const ExactItemKey &key) const;
};

namespace exact_store {

void reset();
void set_total(const ExactItemKey &item, pf_t total);
void erase(const ExactItemKey &item);
pf_t total(const ExactItemKey &item);
using VisitFn = std::function<void(const ExactItemKey &, pf_t)>;
void for_each_item_span(NonTerminal nt,
                        cand_pos_t i,
                        cand_pos_t j,
                        const VisitFn &visit);

}  // namespace exact_store
}  // namespace cparty::scfg2::storage

#endif  // CPARTY_SCFG2_STORAGE_EXACT_STORE_HH

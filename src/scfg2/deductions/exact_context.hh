#ifndef CPARTY_SCFG2_DEDUCTIONS_EXACT_CONTEXT_HH
#define CPARTY_SCFG2_DEDUCTIONS_EXACT_CONTEXT_HH

#include "scfg2/types.hh"
#include "scfg2/grammar/model_policy.hh"

namespace cparty::scfg2::deductions {

class PairabilityContext {
 public:
  virtual ~PairabilityContext() = default;

  virtual WmbClosurePolicy wmb_closure_policy() const {
    return WmbClosurePolicy::RequireWeakClosure;
  }

  virtual cand_pos_t turn() const = 0;
  virtual cand_pos_t sequence_length() const = 0;
  virtual bool can_pair(cand_pos_t i, cand_pos_t j) const = 0;
  virtual cand_pos_t unpaired_prefix(cand_pos_t pos) const = 0;
  virtual cand_pos_t pair_partner(cand_pos_t pos) const {
    (void)pos;
    return -1;
  }
  virtual cand_pos_t pair_round(cand_pos_t pos) const {
    return pair_partner(pos);
  }
  virtual cand_pos_t pair_square(cand_pos_t pos) const {
    (void)pos;
    return -1;
  }
  virtual cand_pos_t pair_any(cand_pos_t pos) const {
    const cand_pos_t round = pair_round(pos);
    return round >= 0 ? round : pair_square(pos);
  }
  virtual cand_pos_t parent_index(cand_pos_t pos) const {
    (void)pos;
    return -1;
  }
  virtual cand_pos_t border_b(cand_pos_t i, cand_pos_t j) const {
    (void)i;
    (void)j;
    return -1;
  }
  virtual cand_pos_t border_bp(cand_pos_t i, cand_pos_t j) const {
    (void)i;
    (void)j;
    return -1;
  }
  virtual cand_pos_t border_B(cand_pos_t i, cand_pos_t j) const {
    (void)i;
    (void)j;
    return -1;
  }
  virtual cand_pos_t border_Bp(cand_pos_t i, cand_pos_t j) const {
    (void)i;
    (void)j;
    return -1;
  }
  virtual bool is_unpaired(cand_pos_t pos) const {
    (void)pos;
    return false;
  }
  virtual bool is_empty_region(cand_pos_t i, cand_pos_t j) const {
    if (j <= i) return false;
    for (cand_pos_t pos = i + 1; pos < j; ++pos) {
      if (!is_unpaired(pos)) {
        return false;
      }
    }
    return true;
  }

  virtual bool can_pair_left_span(cand_pos_t left, cand_pos_t split) const {
    if (split < left) return false;
    if (split == left) return true;
    return unpaired_prefix(split - 1) >= (split - left);
  }

  // [split+1, right] is free in G.  Legacy CParty checked only
  // [split, right-1] (`up[j-1] >= j-k`), which let the VPR right run end on a
  // G-paired base; the full-run guard is also in oracle patch 0006.
  virtual bool can_pair_right_span(cand_pos_t split, cand_pos_t right) const {
    if (right < split) return false;
    if (right == split) return true;
    return unpaired_prefix(right) >= (right - split);
  }

  virtual bool weakly_closed(cand_pos_t i, cand_pos_t j) const {
    if (j <= i) return false;
    for (cand_pos_t pos = i; pos <= j; ++pos) {
      if (!is_unpaired(pos)) {
        return false;
      }
    }
    return true;
  }
};

}  // namespace cparty::scfg2::deductions

#endif  // CPARTY_SCFG2_DEDUCTIONS_EXACT_CONTEXT_HH

#ifndef CPARTY_SCFG2_DEDUCTIONS_SPLIT_POSITIONS_HH
#define CPARTY_SCFG2_DEDUCTIONS_SPLIT_POSITIONS_HH

#include "scfg2/types.hh"
#include <cstddef>
#include <vector>

namespace cparty::scfg2::deductions {

// A view, not a materialized interval. A null selection visits the original
// integer range. An explicit sorted, unique selection visits only its entries
// in range, including O(selection.size()) work for rejected coordinates.
class ExactSplitPositions {
 public:
  struct Iterator {
    cand_pos_t value, lower, upper;
    const std::vector<cand_pos_t> *selected;
    std::size_t index = 0;
    bool done = false;
    void seek() {
      while (index < selected->size() &&
             ((*selected)[index] < lower || (*selected)[index] > upper)) ++index;
      done = index == selected->size();
      if (!done) value = (*selected)[index];
    }
    cand_pos_t operator*() const { return value; }
    Iterator &operator++() {
      if (selected) { ++index; seek(); }
      else if (value == upper) done = true;
      else ++value;
      return *this;
    }
    bool operator!=(const Iterator &other) const { return done != other.done; }
  };
  ExactSplitPositions(cand_pos_t lower, cand_pos_t upper,
                      const std::vector<cand_pos_t> *selected)
      : lower_(lower), upper_(upper), selected_(selected) {}
  Iterator begin() const {
    Iterator it{lower_, lower_, upper_, selected_};
    if (lower_ > upper_) it.done = true;
    else if (selected_) it.seek();
    return it;
  }
  Iterator end() const { return Iterator{0, 0, 0, nullptr, 0, true}; }
 private:
  cand_pos_t lower_, upper_;
  const std::vector<cand_pos_t> *selected_;
};

}  // namespace cparty::scfg2::deductions
#endif

#ifndef CPARTY_SCFG2_BEAM_RETAINED_CHOICES_HH
#define CPARTY_SCFG2_BEAM_RETAINED_CHOICES_HH

#include "scfg2/engine/choice_map.hh"
#include <optional>
#include <vector>

namespace cparty::scfg2::beam {

// Choice observation precedes pruning. Retire only the just-completed
// family/column's rejected keys when the next family is observed, keeping
// O(nb) surviving choices plus current-family scratch. Used also by the
// scaffold-forcing pre-pass, whose choices must not keep all old candidates.
template <class S, class Choices = engine::ChoiceMap<typename S::value_type>>
class RetainedChoices {
 public:
  explicit RetainedChoices(Choices *choices)
      : choices_(choices) {}
  template <class Chart>
  void before_parent(const deductions::ItemKey &parent, const Chart &chart) {
    if (previous_ && (previous_->nonterminal != parent.nonterminal || previous_->j != parent.j))
      finish(chart);
    if (!previous_ || !deductions::ItemKeyEq{}(*previous_, parent)) pending_.push_back(parent);
    previous_ = parent;
  }
  template <class Chart>
  void finish(const Chart &chart) {
    for (const auto &parent : pending_)
      if (S::is_zero(chart.get(parent))) choices_->erase(parent);
    pending_.clear();
    previous_.reset();
  }
 private:
  Choices *choices_;
  std::optional<deductions::ItemKey> previous_;
  std::vector<deductions::ItemKey> pending_;
};

}  // namespace cparty::scfg2::beam
#endif

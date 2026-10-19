#ifndef CPARTY_SCFG2_BEAM_BE_RIGHT_SHELL_SOURCE_AGENDA_HH
#define CPARTY_SCFG2_BEAM_BE_RIGHT_SHELL_SOURCE_AGENDA_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/be_left_right_frame_selection.hh"
#include "scfg2/deductions/expand_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct BeRightShellSourceAgendaStats {
  std::size_t frame_candidates = 0;
  std::size_t emitted_sources = 0;
  std::size_t rejected_sources = 0;
  std::size_t duplicate_sources = 0;
};

// Reifies selected mode-tagged BE_LRTAIL frames as the original public BE
// deductions. Local-weight ownership remains on that original deduction.
class BeRightShellSourceAgenda {
 public:
  template <class Value, class Selection, class Context>
  void prepare_column(
                      const std::vector<BasicScoredBeRightShellFrame<Value>> &frames,
                      const Selection &selection,
                      Context &ctx) {
    current_.clear();
    engine::AllLocalWeightModel<Context> weight(&ctx);
    for (const auto &frame : frames) {
      ++stats_.frame_candidates;
      if (!selection.allows(frame.key)) continue;
      const auto &key = frame.key;
      const auto rule = key.mode == BeRightShellFrameMode::LeftRight
          ? grammar::RuleId::BE_WIP_WIP
          : grammar::RuleId::BE_BASEPAIR_WIP;
      const deductions::ItemKey parent{
          grammar::NonTerminal::BE,
          key.i,
          key.tail.j,
          key.tail.core.ip,
          key.tail.core.jp};
      grammar::RuleSplit split{};
      split.k = key.tail.core.k;
      split.l = key.tail.core.l;
      auto source = deductions::expand_basic_rule(parent, rule, split);
      engine::require_legacy_span_schedule_support(source);
      if (weights::is_zero(weight.local_weight_for(source))) {
        ++stats_.rejected_sources;
        continue;
      }
      auto &sources = current_[parent];
      const bool duplicate = std::any_of(
          sources.begin(), sources.end(), [&](const auto &old) {
            return old.rule == source.rule && old.split.k == source.split.k &&
                old.split.l == source.split.l;
          });
      if (duplicate) {
        ++stats_.duplicate_sources;
        continue;
      }
      sources.push_back(std::move(source));
      ++stats_.emitted_sources;
    }
  }

  const std::vector<deductions::Deduction> &deductions_for(
      const deductions::ItemKey &parent) const {
    const auto found = current_.find(parent);
    if (found == current_.end()) return empty();
    return found->second;
  }

  std::vector<deductions::ItemKey> parents() const {
    std::vector<deductions::ItemKey> out;
    out.reserve(current_.size());
    for (const auto &[parent, sources] : current_) {
      if (!sources.empty()) out.push_back(parent);
    }
    std::sort(out.begin(), out.end(), [](const auto &lhs, const auto &rhs) {
      return std::tie(lhs.i, lhs.ip, lhs.jp) >
          std::tie(rhs.i, rhs.ip, rhs.jp);
    });
    return out;
  }

  const BeRightShellSourceAgendaStats &stats() const { return stats_; }

 private:
  using Sources = std::unordered_map<deductions::ItemKey,
                                     std::vector<deductions::Deduction>,
                                     deductions::ItemKeyHash,
                                     deductions::ItemKeyEq>;

  static const std::vector<deductions::Deduction> &empty() {
    static const std::vector<deductions::Deduction> value;
    return value;
  }

  Sources current_;
  BeRightShellSourceAgendaStats stats_{};
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_RIGHT_SHELL_SOURCE_AGENDA_HH

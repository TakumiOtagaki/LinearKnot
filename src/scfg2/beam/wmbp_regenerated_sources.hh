#ifndef CPARTY_SCFG2_BEAM_WMBP_REGENERATED_SOURCES_HH
#define CPARTY_SCFG2_BEAM_WMBP_REGENERATED_SOURCES_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/wmbp_be_vp_prefix_events.hh"
#include "scfg2/beam/wmbp_wi_prefix_admission.hh"
#include "scfg2/engine/runtime_child_item.hh"

namespace cparty::scfg2::beam {

// The shared exact generator fixes all anchors from (parent, right VP).
// Each of b retained parents therefore probes b retained right starts, not
// every historical BE x WI prefix. The temporal WI/prefix boundaries above
// reproduce admission before the parent-local prefix selection is repeated.
// This provider owns only O(b) per-parent scratch; the boundary history is
// O(nb). Over all columns its work is O(nb^2 log(b+2)) for finite widths.
// vp_starts is sorted and unique, as in the public retained VP view.
template <class Semiring, class Chart, class Context>
std::vector<deductions::Deduction> regenerate_wmbp_sources(
    const deductions::ItemKey &parent,
    const std::vector<cand_pos_t> &vp_starts,
    const Chart &chart, Context &ctx,
    const WmbpWiPrefixAdmission &admission, std::size_t prefix_width) {
  if (Semiring::is_zero(chart.get(parent))) return {};
  auto sources = generate_wmbp_parent_deductions_for_retained_vp_starts(
      parent, vp_starts, ctx);
  engine::AllLocalWeightModel<Context> weight(&ctx);
  sources.erase(std::remove_if(sources.begin(), sources.end(),
      [&](const auto &source) {
        engine::require_legacy_span_schedule_support(source);
        if (weights::is_zero(weight.local_weight_for(source)) ||
            !admission.template allows<Semiring>(source, chart)) return true;
        for (std::size_t c = 0; c < source.child_count; ++c)
          if (Semiring::is_zero(chart.get(engine::RuntimeDeductionChildItem{}(source, c)))) return true;
        return false;
      }), sources.end());
  TopWidthWmbpBeVpPrefixSelection selection(prefix_width);
  selection.template commit_parent<Semiring>(
      parent, score_wmbp_be_vp_prefixes_from_deductions<Semiring>(
          sources, chart, ctx));
  sources.erase(std::remove_if(sources.begin(), sources.end(),
      [&](const auto &source) { return !selection.allows(source); }), sources.end());
  return sources;
}

}  // namespace cparty::scfg2::beam
#endif

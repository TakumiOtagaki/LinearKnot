#ifndef CPARTY_SCFG2_BEAM_RETAINED_OUTSIDE_HH
#define CPARTY_SCFG2_BEAM_RETAINED_OUTSIDE_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/be_regenerated_sources.hh"
#include "scfg2/beam/retained_chart_views.hh"
#include "scfg2/beam/retained_rule_admission.hh"
#include "scfg2/beam/vm_retained_join_agenda.hh"
#include "scfg2/beam/vp_recursive_split_retained_events.hh"
#include "scfg2/beam/vpr_regenerated_base_cells.hh"
#include "scfg2/beam/wi_retained_join_agenda.hh"
#include "scfg2/beam/wip_retained_join_agenda.hh"
#include "scfg2/beam/wm_retained_join_agenda.hh"
#include "scfg2/beam/wmb_regenerated_btail_cells.hh"
#include "scfg2/beam/wmbp_regenerated_sources.hh"
#include "scfg2/beam/wmbw_retained_join_agenda.hh"
#include "scfg2/beam/w_retained_join_agenda.hh"
#include "scfg2/deductions/retained_start_splits.hh"
#include "scfg2/engine/outside.hh"

namespace cparty::scfg2::beam {

// Shared deduction consumers over immutable retained inside alphas. The
// only persistent numerical private cells are VPR carries and WMB BTAIL.
// All ordinary branches, start-prefix sums and BE/WMBP products are replayed
// for the current column, then discarded.
template <class S, class Context>
class RetainedOutsideProgram {
 public:
  using Value = typename S::value_type;
  using Key = deductions::ItemKey;
  using D = deductions::Deduction;
  using NT = grammar::NonTerminal;
  using R = grammar::RuleId;
  using Chart = engine::Chart<Value>;
  struct Result {
    Chart beta{S::zero()};
    typename VprRegeneratedBaseCells<S>::Outside vpr_beta;
    typename WmbRegeneratedBtailCells<S>::Outside btail_beta;
  };

  template <class Widths>
  RetainedOutsideProgram(cand_pos_t n, const Chart &inside, Context &ctx,
                          const RetainedRuleAdmission<S> &rules, const Widths &widths,
                          const PublicItemRanking &ranking = {})
      : n_(n), inside_(inside), ctx_(ctx), rules_(rules),
        views_(RetainedChartViews::build<S>(n, inside, widths.wmb, ranking)),
        be_(BeRegeneratedSources<S>::from_retained_chart(n, inside, ctx)),
        prefix_(WmbpWiPrefixAdmission::from_retained_chart<S>(
            n, inside, widths.wmbp_wi_history, widths.wmbp_wi_prefix, ranking)),
        tail_width_(widths.be_right_tail), shell_width_(widths.be_right_shell),
        prefix_width_(widths.wmbp_prefix) {
    const auto starts = [&](NT nt, cand_pos_t j) -> const auto & { return views_.starts(nt, j); };
    vpr_ = VprRegeneratedBaseCells<S>::from_retained_chart(n, inside, ctx, starts, widths.vpr);
    btail_ = WmbRegeneratedBtailCells<S>::from_retained_chart(n, inside, ctx, starts, widths.wmb_btail);
  }

  const RetainedChartViews &views() const { return views_; }
  std::size_t private_cell_count() const { return vpr_.size() + btail_.size(); }
  const auto &vpr_cells(cand_pos_t j) const { return vpr_.cells(j); }

  void prepare_column(cand_pos_t j) {
    prepare_column(j, true);
  }

 private:
  // The alpha projection is needed only by explicit source-reconstruction
  // observations. Reverse propagation uses the shared term generator below.
  void prepare_column(cand_pos_t j, bool project_wmb_alpha) {
    const auto starts = [&](NT nt, cand_pos_t k) -> const auto & {
      return views_.internal_starts(nt, k);
    };
    vm_.build_column(j, inside_, ctx_, starts, 0);
    wm_.build_column(j, inside_, ctx_, starts, 0);
    wip_.build_column(j, inside_, ctx_, starts, 0);
    wi_.build_column(j, ctx_, starts);
    wmbw_.build_column(j, ctx_,
        [&](cand_pos_t k) -> const auto & { return views_.starts(NT::WI, k); },
        [&](cand_pos_t k) -> const auto & { return views_.starts(NT::WMBP, k); });
    vpr_binary_.build_column(j, inside_, ctx_, views_.starts(NT::WIP, j),
        [&](cand_pos_t k) -> const auto & { return views_.starts(NT::VP, k); }, 0);
    w_.build_column(j, ctx_, [&](NT nt, cand_pos_t k) -> const auto & {
      return views_.root_starts(nt, k);
    }, [&](cand_pos_t i, cand_pos_t k) {
      const auto &starts = views_.root_starts(NT::WMB, k);
      return std::binary_search(starts.begin(), starts.end(), i);
    });
    const auto recursive = [&](auto &events, NT vp) {
      events.prepare_column(j,
          [&](cand_pos_t k) -> const auto & { return views_.starts(vp, k); },
          [&](cand_pos_t k) -> const auto & { return views_.starts(NT::WIP, k); },
          [&](cand_pos_t k) -> const auto & { return views_.starts(NT::VPL, k); },
          [&](cand_pos_t k) -> const auto & { return views_.starts(NT::VPR, k); }, ctx_);
    };
    recursive(vp_recursive_, NT::VP);
    be_.prepare_column(j, inside_, ctx_, tail_width_, shell_width_);
    wmb_split_alpha_.clear();
    if (!project_wmb_alpha) return;
    for_each_wmb_be_btail_term(j, ctx_, views_.items(NT::BE, j), btail_.store(),
        [&](const auto &source, const auto &be, const auto &tail, auto local) {
          if (S::is_zero(inside_.get(source.parent))) return;
          const auto term = S::times(S::from_weight(local), S::times(inside_.get(be), tail.total));
          auto [found, inserted] = wmb_split_alpha_.try_emplace(source.parent, S::zero());
          (void)inserted;
          found->second = S::plus(found->second, term);
        });
  }

 public:
  // One read per parent after prepare_column (recursive events are moved).
  std::vector<D> deductions_for(const Key &parent) {
    const auto admitted_rules = rules_.rules_for(parent);
    std::vector<D> out;
    const auto append_rule = [&](R rule) {
      if (!rules_.has(admitted_rules, rule)) return;
      auto sources = deductions::generate_exact_basic_deductions(parent, rule, ctx_);
      out.insert(out.end(), sources.begin(), sources.end());
    };
    switch (parent.nonterminal) {
      case NT::BE: out = be_.deductions_for(parent, inside_, ctx_); break;
      case NT::WMBP:
        out = regenerate_wmbp_sources<S>(parent, views_.starts(NT::VP, parent.j),
            inside_, ctx_, prefix_, prefix_width_);
        break;
      case NT::VM: out = vm_.deductions_for(parent); break;
      case NT::WM: out = wm_.deductions_for(parent); break;
      case NT::WIP: out = wip_.deductions_for(parent); break;
      case NT::WI: out = wi_.deductions_for(parent); break;
      case NT::WMBW: out = wmbw_.deductions_for(parent); break;
      case NT::VPR: out = vpr_binary_.deductions_for(parent); break;
      case NT::W: out = w_.deductions_for(parent); break;
      case NT::VPL: break;
      case NT::WMB: append_rule(R::WMB_DIRECT_WMBP); append_rule(R::WMB_EMPTY); break;
      case NT::V:
        append_rule(R::V_HAIRPIN); append_rule(R::V_INTERNAL); append_rule(R::V_VM); break;
      case NT::WMv: append_rule(R::WMv_STEM_V); append_rule(R::WMv_EXTEND_UNPAIRED); break;
      case NT::WMp: append_rule(R::WMp_STEM_WMB); append_rule(R::WMp_EXTEND_UNPAIRED); break;
      case NT::VP:
        append_rule(R::VP_WI_CASE1); append_rule(R::VP_WI_CASE2); append_rule(R::VP_WI_CASE3);
        append_rule(R::VP_STACK); append_rule(R::VP_INTERNAL_LOOP);
        auto sources = vp_recursive_.take_for(parent);
        out.insert(out.end(), sources.begin(), sources.end());
        break;
    }
    bool start = false;
    if (parent.nonterminal == NT::VM) start = rules_.has(admitted_rules, R::VM_SPLIT_WMp_BASE);
    if (parent.nonterminal == NT::WM)
      start = rules_.has(admitted_rules, R::WM_START_V) || rules_.has(admitted_rules, R::WM_START_WMB);
    if (parent.nonterminal == NT::WIP)
      start = rules_.has(admitted_rules, R::WIP_BASE_V) || rules_.has(admitted_rules, R::WIP_BASE_WMB) ||
              rules_.has(admitted_rules, R::WIP_BASEPAIR_V) || rules_.has(admitted_rules, R::WIP_BASEPAIR_WMB);
    if (parent.nonterminal == NT::VPL) start = rules_.has(admitted_rules, R::VPL_SPLIT_VP);
    if (start) deductions::for_each_retained_start_split(parent, ctx_,
        [&](NT nt, cand_pos_t j) -> const auto & { return views_.internal_starts(nt, j); },
        [&](const auto &source) { out.push_back(source); });
    engine::AllLocalWeightModel<Context> weight(&ctx_);
    out.erase(std::remove_if(out.begin(), out.end(), [&](const auto &source) {
      engine::require_legacy_span_schedule_support(source);
      if ((!start_rule(source.rule) && !rules_.has(admitted_rules, source.rule)) ||
          weights::is_zero(weight.local_weight_for(source)))
        return true;
      for (std::size_t c = 0; c < source.child_count; ++c)
        if (S::is_zero(inside_.get(engine::RuntimeDeductionChildItem{}(source, c)))) return true;
      return false;
    }), out.end());
    return out;
  }

  Value private_contribution(const Key &parent) const {
    if (parent.nonterminal == NT::WMB) {
      const auto found = wmb_split_alpha_.find(parent);
      return found == wmb_split_alpha_.end() ? S::zero() : found->second;
    }
    if (parent.nonterminal == NT::VPR) {
      const auto found = vpr_.cells(parent.j).find(parent.i);
      if (found != vpr_.cells(parent.j).end() && found->second.public_base) return found->second.alpha;
    }
    return S::zero();
  }

  struct IgnoreDeduction {
    void operator()(const D &, Value) const {}
  };

  template <class Observer = IgnoreDeduction>
  Result outside(const Key &root, Observer &&observe = {}) {
    Result result;
    result.vpr_beta = vpr_.empty_outside();
    result.btail_beta = btail_.empty_outside();
    result.beta.set(root, S::one());
    engine::AllLocalWeightModel<Context> weight(&ctx_);
    for (cand_pos_t j = n_; j >= 1; --j) {
      prepare_column(j, false);
      for (auto phase = vp_wmbp_online_column_phases.rbegin();
           phase != vp_wmbp_online_column_phases.rend(); ++phase) {
        // Starts are ascending: root W(1,j) precedes same-column suffix W.
        for (const auto &parent : views_.items(*phase, j)) {
          const auto beta = result.beta.get(parent);
          if (S::is_zero(beta)) continue;
          auto sources = deductions_for(parent);
          for (const auto &source : sources) {
            observe(source, beta);
            engine::propagate_outside_deduction<S>(source, inside_, beta, &result.beta,
                [&](const auto &d) { return weight.local_weight_for(d); }, engine::RuntimeDeductionChildItem{});
          }
          if (*phase == NT::VPR) vpr_.add_public_outside(parent, beta, &result.vpr_beta);
        }
        if (*phase == NT::WMB)
          btail_.propagate_wmb_splits(j, inside_, ctx_, views_.items(NT::BE, j),
              &result.beta, &result.btail_beta);
      }
      vpr_.propagate_column(j, inside_, ctx_, &result.vpr_beta, &result.beta);
      btail_.propagate_column(j, inside_, ctx_,
          [&](NT nt, cand_pos_t k) -> const auto & { return views_.starts(nt, k); },
          &result.btail_beta, &result.beta);
    }
    return result;
  }

 private:
  static bool start_rule(R rule) {
    return rule == R::VM_SPLIT_WMp_BASE || rule == R::WM_START_V || rule == R::WM_START_WMB ||
           rule == R::WIP_BASE_V || rule == R::WIP_BASE_WMB || rule == R::WIP_BASEPAIR_V ||
           rule == R::WIP_BASEPAIR_WMB || rule == R::VPL_SPLIT_VP;
  }
  cand_pos_t n_;
  const Chart &inside_;
  Context &ctx_;
  const RetainedRuleAdmission<S> &rules_;
  RetainedChartViews views_;
  BeRegeneratedSources<S> be_;
  WmbpWiPrefixAdmission prefix_;
  std::size_t tail_width_, shell_width_, prefix_width_;
  VprRegeneratedBaseCells<S> vpr_;
  WmbRegeneratedBtailCells<S> btail_;
  VmRetainedJoinAgenda<S> vm_;
  WmRetainedJoinAgenda<S> wm_;
  WipRetainedJoinAgenda<S> wip_;
  WiRetainedJoinAgenda wi_;
  WmbwRetainedJoinAgenda wmbw_;
  VprRetainedAgenda<S> vpr_binary_;
  WRetainedJoinAgenda w_;
  VpRecursiveSplitRetainedEvents vp_recursive_;
  std::unordered_map<Key, Value, deductions::ItemKeyHash, deductions::ItemKeyEq> wmb_split_alpha_;
};

}  // namespace cparty::scfg2::beam
#endif

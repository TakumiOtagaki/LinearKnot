#ifndef CPARTY_SCFG2_REPLAY_W_FINAL_EXACT_INSIDE_HH
#define CPARTY_SCFG2_REPLAY_W_FINAL_EXACT_INSIDE_HH

#include "scfg2/engine/runtime_child_item.hh"

#include "rna_pair_table.hh"
#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/storage/exact_store.hh"
#include "scfg2/deductions/exact_context.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/geometry/static_scaffold_geometry.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/engine/basic_exact_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"
#include "scfg2/engine/viterbi_result.hh"
#include "scfg2/semirings/semiring.hh"
#include "scfg2/weights/context_be.hh"
#include "scfg2/weights/log_loops.hh"
#include "scfg2/weights/context_vp.hh"
#include "scfg2/weights/context_vpl.hh"
#include "scfg2/weights/context_vpr.hh"
#include "scfg2/weights/context_w.hh"
#include "scfg2/weights/context_wi.hh"
#include "scfg2/weights/context_wip.hh"
#include "scfg2/weights/context_wmb.hh"
#include "scfg2/weights/context_wmbp.hh"
#include "scfg2/weights/context_wmbw.hh"

#include <initializer_list>
#include <array>
#include <optional>
#include <algorithm>
#include <type_traits>
#include <utility>
#include <vector>

namespace cparty::scfg2::replay {

template <bool LogDomain>
class WFinalInsideContext final : public engine::BasicExactWeightContext,
                                       public weights::WTransitionContext,
                                       public weights::WITransitionContext,
                                       public weights::WIPTransitionContext,
                                       public weights::VPLTransitionContext,
                                       public weights::VPRTransitionContext,
                                       public weights::VPTransitionContext,
                                       public weights::WMBWTransitionContext,
                                       public weights::WMBPTransitionContext,
                                       public weights::WMBTransitionContext,
                                       public weights::BETransitionContext,
                                       public deductions::PairabilityContext {
 public:
  static constexpr bool log_local_weights = LogDomain;

  explicit WFinalInsideContext(weights::RuntimeSnapshot &owner,
      WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure)
      : owner_(owner), wmb_closure_policy_(policy),
        tree_(owner_.scaffold(), owner_.length()) {}

  double exterior_unpaired_weight() const override {
    return unit();
  }

  WmbClosurePolicy wmb_closure_policy() const override {
    return wmb_closure_policy_;
  }

  double hairpin_weight(cand_pos_t i, cand_pos_t j) const override {
    return hairpin_factor(i, j);
  }

  double internal_weight(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const override {
    return internal_factor(i, j, k, l);
  }

  double vm_closing_weight(cand_pos_t i, cand_pos_t j) const override {
    (void)i;
    (void)j;
    // Exact V->VM carries the multiloop total through the VM child itself.
    // The current exact deduction uses local weight 1 here.
    return unit();
  }

  double split_wm_wmv_weight(cand_pos_t i, cand_pos_t j) const override {
    return product({mlclosing_stem_factor(i, j), mlclosing_factor()});
  }

  double split_wm_wmp_weight(cand_pos_t i, cand_pos_t j) const override {
    return split_wm_wmv_weight(i, j);
  }

  double split_wmp_base_weight(cand_pos_t i,
                               cand_pos_t j,
                               cand_pos_t left_unpaired) const override {
    return product({mlbase_factor(left_unpaired), split_wm_wmp_weight(i, j)});
  }

  double vm_base_unpaired_step_weight() const override {
    return mlbase_factor(1);
  }

  double vm_base_gap_weight(cand_pos_t left_unpaired) const override {
    return mlbase_factor(left_unpaired);
  }

  double vm_wrapper_weight() const override {
    return unit();
  }

  double wmv_stem_weight(cand_pos_t i, cand_pos_t j) const override {
    return mlstem_factor(i, j);
  }

  double wmp_stem_weight(cand_pos_t, cand_pos_t) const override {
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::psm), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::b)});
  }

  double wmv_extend_unpaired_weight() const override {
    return mlbase_factor(1);
  }

  double wmp_extend_unpaired_weight() const override {
    return mlbase_factor(1);
  }

  double start_v_weight(cand_pos_t i, cand_pos_t j, cand_pos_t k) const override {
    return product({mlbase_factor(k - i), mlstem_factor(k, j)});
  }

  double start_wmb_weight(cand_pos_t i, cand_pos_t, cand_pos_t k) const override {
    return product({mlbase_factor(k - i), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::psm), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::b)});
  }

  double split_v_weight(cand_pos_t k, cand_pos_t j) const override {
    return mlstem_factor(k, j);
  }

  double split_wmb_weight(cand_pos_t, cand_pos_t) const override {
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::psm), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::b)});
  }

  double extend_unpaired_weight() const override {
    return mlbase_factor(1);
  }

  double exp_extloop(cand_pos_t i, cand_pos_t j) const override {
    return external_factor(i, j);
  }

  double w_extend_unpaired_weight() const override {
    return exterior_unpaired_weight();
  }

  double w_split_v_weight(cand_pos_t i, cand_pos_t j) const override {
    return exp_extloop(i, j);
  }

  double w_split_wmb_weight() const override {
    return pk_factor(&weights::RuntimeSnapshot::LogPkParameters::ps);
  }

  double w_root_wmb_decomposition_weight() const override {
    return pk_factor(&weights::RuntimeSnapshot::LogPkParameters::ps);
  }

  double wi_base_single_weight() const override {
    return pup_factor(1);
  }

  double wi_split_v_weight() const override {
    return pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pps);
  }

  double wi_split_wmb_weight() const override {
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::psp), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pps)});
  }

  double wi_extend_unpaired_weight() const override {
    return pup_factor(1);
  }

  double wip_base_v_weight() const override {
    return pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp);
  }

  double wip_base_wmb_weight() const override {
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::psm)});
  }

  double wip_split_v_weight() const override {
    return pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp);
  }

  double wip_split_wmb_weight() const override {
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::psm)});
  }

  double wip_basepair_v_weight(cand_pos_t left_unpaired) const override {
    return product({cp_factor(left_unpaired), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp)});
  }

  double wip_basepair_wmb_weight(cand_pos_t left_unpaired) const override {
    return product({cp_factor(left_unpaired), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::psm)});
  }

  double wip_extend_unpaired_weight() const override {
    return cp_factor(1);
  }

  double vpl_split_vp_weight(cand_pos_t left_unpaired) const override {
    return cp_factor(left_unpaired);
  }

  double vpr_split_vp_wip_weight() const override {
    return unit();
  }

  double vpr_split_vp_basepair_weight(cand_pos_t right_unpaired) const override {
    return cp_factor(right_unpaired);
  }

  double vp_wi_case_weight() const override {
    return unit();
  }

  double vp_stack_weight(cand_pos_t i, cand_pos_t j) const override {
    return stack_pk_factor(i, j);
  }

  double vp_internal_loop_weight(cand_pos_t i,
                                 cand_pos_t j,
                                 cand_pos_t k,
                                 cand_pos_t l) const override {
    return internal_pk_factor(i, k, l, j);
  }

  double vp_split_weight() const override {
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::ap), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp)});
  }

  cand_pos_t turn() const override { return TURN; }

  cand_pos_t sequence_length() const override { return owner_.length(); }

  bool can_pair(cand_pos_t i, cand_pos_t j) const override {
    if (i < 1 || j < 1 || i > j) {
      return false;
    }
    return owner_.bases_can_pair(i, j);
  }

  cand_pos_t unpaired_prefix(cand_pos_t pos) const override {
    const auto &up = tree_.unpaired_counts();
    if (pos <= 0 || static_cast<size_t>(pos) >= up.size()) {
      return 0;
    }
    return static_cast<cand_pos_t>(up[static_cast<size_t>(pos)]);
  }

  cand_pos_t pair_partner(cand_pos_t pos) const override {
    return tree_.partner(pos);
  }

  cand_pos_t pair_round(cand_pos_t pos) const override {
    return pair_partner(pos);
  }

  cand_pos_t parent_index(cand_pos_t pos) const override {
    return tree_.parent(pos);
  }

  cand_pos_t border_b(cand_pos_t i, cand_pos_t j) const override { return tree_.b(i, j); }
  cand_pos_t border_bp(cand_pos_t i, cand_pos_t j) const override { return tree_.bp(i, j); }
  cand_pos_t border_B(cand_pos_t i, cand_pos_t j) const override { return tree_.B(i, j); }
  cand_pos_t border_Bp(cand_pos_t i, cand_pos_t j) const override { return tree_.Bp(i, j); }

  bool is_unpaired(cand_pos_t pos) const override {
    return pos >= 1 && pos <= tree_.size() && tree_.partner(pos) < 0;
  }

  bool is_empty_region(cand_pos_t i, cand_pos_t j) const override {
    if (j <= i) return false;
    if (j == i + 1) return true;
    const cand_pos_t up_index = j - 1;
    const auto &up = tree_.unpaired_counts();
    if (up.empty() || up_index <= 0 ||
        static_cast<std::size_t>(up_index) >= up.size()) {
      return false;
    }
    return static_cast<cand_pos_t>(up[static_cast<std::size_t>(up_index)]) >=
           (j - i - 1);
  }

  bool weakly_closed(cand_pos_t i, cand_pos_t j) const override {
    return tree_.weakly_closed(i, j);
  }

  double wmbp_split_be_wmbp_vp_weight(
      const weights::WmbpSharedBeVpAnchors &) const override {
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pb), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pb)});
  }

  double wmbp_split_be_wmbw_vp_weight(
      const weights::WmbpSharedBeVpAnchors &) const override {
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pb), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pb)});
  }

  double wmbp_direct_vp_weight(cand_pos_t i, cand_pos_t j) const override {
    (void)i;
    (void)j;
    // BE/VP children own the scaffold support; this wrapper adds only PB.
    return pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pb);
  }

  double wmbp_split_be_wi_vp_weight(
      const deductions::WmbpWiVpAnchors &anchors) const override {
    if (anchors.be_ip < 1 || anchors.be_jp < anchors.be_ip) {
      return impossible();
    }
    // Anchored-BE support lives on the explicit BE child. This rule adds PB^2.
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pb), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pb)});
  }

  double wmbw_split_wmbp_wi_weight() const override {
    return unit();
  }

  double wmb_split_be_wmbp_wi_weight(
      const weights::WmbSplitBeWmbpWiWeightAnchors &anchors) const override {
    if (anchors.be_j < 1 || anchors.be_jp < 1 || anchors.be_ip < 1) {
      return impossible();
    }
    return pk_factor(&weights::RuntimeSnapshot::LogPkParameters::pb);
  }

  double wmb_direct_wmbp_weight() const override {
    // Keep the direct wrapper local term minimal until WMBP/WMB exact-block
    // ownership is rewritten explicitly. Prior probes showed that a naive
    // symmetric "re-sponsor the right block here" tweak was not sufficient.
    return unit();
  }

  double wmb_empty_weight() const override {
    return impossible();
  }

  double be_base_samepair_weight() const override {
    return unit();
  }

  double be_stack_weight(cand_pos_t i, cand_pos_t j) const override {
    return stack_pk_factor(i, j);
  }

  double be_internal_loop_weight(cand_pos_t i,
                                 cand_pos_t j,
                                 cand_pos_t k,
                                 cand_pos_t l) const override {
    return internal_pk_factor(i, k, l, j);
  }

  double be_split_wip_wip_weight() const override {
    return product({pk_factor(&weights::RuntimeSnapshot::LogPkParameters::ap), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp)});
  }

  double be_split_wip_basepair_weight(cand_pos_t right_unpaired) const override {
    return product({cp_factor(right_unpaired), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::ap), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp)});
  }

  double be_split_basepair_wip_weight(cand_pos_t left_unpaired) const override {
    return product({cp_factor(left_unpaired), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::ap), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp), pk_factor(&weights::RuntimeSnapshot::LogPkParameters::bp)});
  }

 private:

  static constexpr double unit() { return LogDomain ? 0. : 1.; }
  static constexpr double impossible() { return LogDomain ? weights::log_zero : 0.; }
  static double product(std::initializer_list<double> factors) {
    double value = unit();
    for (double factor : factors) {
      if constexpr (LogDomain) value += factor;
      else value *= factor;
    }
    return value;
  }
  // Raw observation is only a projection of the same unscaled log factors.
  // No second parameter domain or numerical scaling lives in the owner.
  static double project_log_factor(double log_factor) {
    if constexpr (LogDomain) return log_factor;
    else return std::exp(log_factor);
  }
  double mlbase_factor(cand_pos_t length) const {
    return project_log_factor(length == 0 ? 0. : length * owner_.log_pf_parameters().expMLbase);
  }
  double cp_factor(cand_pos_t length) const {
    return project_log_factor(length == 0 ? 0. : length * owner_.log_pk_parameters().cp);
  }
  double pup_factor(cand_pos_t length) const {
    return project_log_factor(length == 0 ? 0. : length * owner_.log_pk_parameters().pup);
  }
  double pk_factor(double weights::RuntimeSnapshot::LogPkParameters::*member) const {
    return project_log_factor(owner_.log_pk_parameters().*member);
  }
  double mlclosing_factor() const {
    return project_log_factor(owner_.log_pf_parameters().expMLclosing);
  }
  double hairpin_factor(cand_pos_t i, cand_pos_t j) const {
    const int type = cparty::rna::pair_type_for(owner_.encoded_base(i), owner_.encoded_base(j));
    if (!type) return impossible();
    return project_log_factor(weights::log_hairpin(j-i-1, type, owner_.adjacent_base(i+1),
        owner_.adjacent_base(j-1), owner_.sequence().c_str()+i-1,
        owner_.log_pf_parameters()));
  }
  double internal_log_factor(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
    return weights::log_internal_loop(k-i-1,j-l-1,
        cparty::rna::pair_type_for(owner_.encoded_base(i),owner_.encoded_base(j)),
        cparty::rna::reverse_pair_type(cparty::rna::pair_type_for(owner_.encoded_base(k),owner_.encoded_base(l))),
        owner_.adjacent_base(i+1),owner_.adjacent_base(j-1),
        owner_.adjacent_base(k-1),owner_.adjacent_base(l+1),
        owner_.log_pf_parameters(), owner_.loop_extrapolation());
  }
  double internal_factor(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
    return project_log_factor(internal_log_factor(i,j,k,l));
  }
  double stack_pk_factor(cand_pos_t i, cand_pos_t j) const {
    return i+1 == j-1 ? impossible() : project_log_factor(
        owner_.pseudoknot_model().stack_multiplier * internal_log_factor(i,j,i+1,j-1));
  }
  double internal_pk_factor(cand_pos_t i,cand_pos_t k,cand_pos_t l,cand_pos_t j) const {
    return k == i+1 && l == j-1 ? impossible() : project_log_factor(
        owner_.pseudoknot_model().internal_multiplier * internal_log_factor(i,j,k,l));
  }
  double stem_factor(cand_pos_t i,cand_pos_t j,bool multiloop,bool closing) const {
    const auto &p = owner_.log_pf_parameters();
    const int type = closing ? cparty::rna::pair_type_for(owner_.encoded_base(j),owner_.encoded_base(i))
                             : cparty::rna::pair_type_for(owner_.encoded_base(i),owner_.encoded_base(j));
    int si=-1,sj=-1;
    if (p.model_details.dangles == 1 || p.model_details.dangles == 2) {
      if (closing) {
        // Both neighbors are inside a multiloop's closing pair, even at
        // sequence endpoints. Reverse the pair and neighbor orientation.
        si = owner_.encoded_base(j-1);
        sj = owner_.encoded_base(i+1);
      } else {
        si = i > 1 ? owner_.encoded_base(i-1) : -1;
        sj = j < owner_.length() ? owner_.encoded_base(j+1) : -1;
      }
    }
    return project_log_factor(weights::log_stem(type,si,sj,multiloop,p));
  }
  double mlclosing_stem_factor(cand_pos_t i,cand_pos_t j) const {
    return stem_factor(i,j,true,true);
  }
  double mlstem_factor(cand_pos_t i,cand_pos_t j) const {
    return stem_factor(i,j,true,false);
  }
  double external_factor(cand_pos_t i,cand_pos_t j) const {
    return stem_factor(i,j,false,false);
  }

  weights::RuntimeSnapshot &owner_;
  WmbClosurePolicy wmb_closure_policy_;
  geometry::StaticScaffoldGeometry tree_;
};

// Bounded raw projections are retained for historical observation witnesses.
// Both views read the same log-only parameter snapshot.
using WFinalExactInsideContext = WFinalInsideContext<false>;
using WFinalLogWeightContext = WFinalInsideContext<true>;

// Explicit numeric-exact policy. Reference callers retain the complete chart
// and candidate observations. The fast policy omits zero-child observations
// and VM values whose sole canonical consumer, V(i,j), cannot be generated.
// It must not be used for finite-beam retention or custom child products.
struct ExactExecutionOptions {
  bool dense_primary = false;
  bool skip_zero_children = false;
  bool omit_unused_vm = false;
  // Also permits skipping admission-predicate calls on proved-zero branches.
  // Enable only for a pure predicate and canonical numeric child products.
  bool filter_zero_candidates = false;

  bool reuse_candidate_buffer = false;

  static constexpr ExactExecutionOptions numeric() { return {true, true, true, true, true}; }
};

// Run-local reusable selections. Child identities come from the shared fixed
// expansion, never from a second rule schema. The generator still owns
// all geometric guards and final rule/split order.
struct ExactNumericSplitScratch {
  std::array<std::vector<cand_pos_t>, 16> selected;

  deductions::ExactNumericSplitViews views() const {
    return {&selected[0], &selected[1], &selected[2], &selected[3],
            &selected[4], &selected[5], &selected[6], &selected[7], &selected[8],
            &selected[9], &selected[10], &selected[11], &selected[12], &selected[13], &selected[14], &selected[15]};
  }

  static bool supports_item(grammar::NonTerminal nt) {
    using N = grammar::NonTerminal;
    return nt == N::WM || nt == N::VM || nt == N::WI || nt == N::WIP || nt == N::VPR || nt == N::WMBP;
  }

  template <class Semiring, grammar::RuleId Rule, class Chart>
  static bool children_nonzero(const deductions::ItemKey &parent,
                               cand_pos_t k, const Chart &chart) {
    const auto d = deductions::expand_basic_rule<Rule>(parent, grammar::RuleSplit{k});
    // Keep schedule diagnostics before reading any child, including on branches
    // that will be omitted. The hot predicate can inline for the fixed rule;
    // full diagnostic formatting remains on the failure path.
    if (!engine::legacy_span_schedule_supports(d))
      engine::require_legacy_span_schedule_support(d);
    for (std::size_t c = 0; c < d.child_count; ++c)
      if (Semiring::is_zero(chart.get(d.children[c].item))) return false;
    return true;
  }

  template <class Semiring, class Chart, class Context>
  static bool standard_live_index(const deductions::ItemKey &parent,
                                   const Chart &chart, const Context &ctx) {
    if constexpr (std::is_same_v<std::decay_t<Chart>, engine::Chart<double>> &&
                  std::is_same_v<std::decay_t<Context>, WFinalLogWeightContext> &&
                  (std::is_same_v<Semiring, semirings::LogSumExpSemiring> ||
                   std::is_same_v<Semiring, semirings::LogMaxSemiring>)) {
      const auto n = chart.live_primary_index_length();
      return n && parent.ip == -1 && parent.jp == -1 && parent.i >= 1 &&
             parent.i <= parent.j && parent.j <= n && ctx.turn() >= 0;
    } else {
      return false;
    }
  }

  // Infer the required right child's row and split-coordinate offset from
  // shared expansion. Every indexed split still executes the complete shared
  // expansion and schedule check in children_nonzero. Unsupported contexts,
  // charts, parent keys and child shapes retain the old full-range scan.
  // For these fixed rules on canonical spans, every omitted dependency also
  // precedes its parent: WM's equal-span V/WMB children have earlier phases;
  // WI/WIP exclude their last endpoint; VM's children are inside (i,j).
  // The left child, when present, is strictly shorter as well. The endpoint
  // checks below retain fallback if the shared schema changes that premise.
  template <class Semiring, grammar::RuleId First, grammar::RuleId Second,
            class Chart, class Context, class Fn>
  static void for_each_right_live_start(const deductions::ItemKey &parent,
                                        const Chart &chart, const Context &ctx,
                                        cand_pos_t begin, cand_pos_t end, Fn &&fn) {
    if (begin > end) return;
    if constexpr (std::is_same_v<std::decay_t<Chart>, engine::Chart<double>>) {
      if (standard_live_index<Semiring>(parent, chart, ctx) &&
          begin >= parent.i && end <= parent.j) {
        const auto left_begin = deductions::expand_basic_rule<First>(parent, grammar::RuleSplit{begin});
        const auto left_end = deductions::expand_basic_rule<First>(parent, grammar::RuleSplit{end});
        const auto right_begin = deductions::expand_basic_rule<Second>(parent, grammar::RuleSplit{begin});
        const auto right_end = deductions::expand_basic_rule<Second>(parent, grammar::RuleSplit{end});
        if (left_begin.child_count && left_end.child_count &&
            right_begin.child_count && right_end.child_count &&
            engine::legacy_span_schedule_supports(left_begin) &&
            engine::legacy_span_schedule_supports(left_end) &&
            engine::legacy_span_schedule_supports(right_begin) &&
            engine::legacy_span_schedule_supports(right_end)) {
          const auto &a = left_begin.children[left_begin.child_count - 1].item;
          const auto &b = left_end.children[left_end.child_count - 1].item;
          const auto &c = right_begin.children[right_begin.child_count - 1].item;
          const auto &d = right_end.children[right_end.child_count - 1].item;
          const auto offset = a.i - begin;
          if (a.ip == -1 && a.jp == -1 && b.ip == -1 && b.jp == -1 &&
              c.ip == -1 && c.jp == -1 && d.ip == -1 && d.jp == -1 &&
              a.nonterminal == b.nonterminal && a.j == b.j &&
              c.nonterminal == d.nonterminal && c.j == d.j &&
              b.i - end == offset && c.i - begin == offset && d.i - end == offset &&
              chart.for_each_live_primary_start_union(a.nonterminal, a.j,
                  c.nonterminal, c.j, begin + offset, end + offset,
                  [&](cand_pos_t start) { fn(start - offset); })) return;
        }
      }
    }
    for (cand_pos_t k = begin; k <= end; ++k) fn(k);
  }

  template <class Semiring, grammar::RuleId Rule, class Chart, class Context, class Fn>
  static void for_each_left_live_end(const deductions::ItemKey &parent,
                                     const Chart &chart, const Context &ctx,
                                     cand_pos_t begin, cand_pos_t end, Fn &&fn) {
    if (begin > end) return;
    if constexpr (std::is_same_v<std::decay_t<Chart>, engine::Chart<double>>) {
      if (standard_live_index<Semiring>(parent, chart, ctx) &&
          begin > parent.i && end < parent.j) {
        const auto first = deductions::expand_basic_rule<Rule>(parent, grammar::RuleSplit{begin});
        const auto last = deductions::expand_basic_rule<Rule>(parent, grammar::RuleSplit{end});
        if (first.child_count && last.child_count &&
            engine::legacy_span_schedule_supports(first) &&
            engine::legacy_span_schedule_supports(last)) {
          const auto &a = first.children[0].item;
          const auto &b = last.children[0].item;
          const auto offset = a.j - begin;
          if (a.ip == -1 && a.jp == -1 && b.ip == -1 && b.jp == -1 &&
              a.nonterminal == b.nonterminal && a.i == b.i && b.j - end == offset &&
              chart.for_each_live_primary_end(a.nonterminal, a.i,
                  begin + offset, end + offset,
                  [&](cand_pos_t finish) { fn(finish - offset); })) return;
        }
      }
    }
    for (cand_pos_t k = begin; k <= end; ++k) fn(k);
  }

  template <class Semiring, class Chart, class Context, class Fn>
  static void for_each_wmbp_live_start(const deductions::ItemKey &parent,
                                       const Chart &chart, const Context &ctx,
                                       cand_pos_t begin, cand_pos_t end, Fn &&fn) {
    if (begin > end) return;
    if constexpr (std::is_same_v<std::decay_t<Chart>, engine::Chart<double>>) {
      if (standard_live_index<Semiring>(parent, chart, ctx) &&
          begin > parent.i && end < parent.j) {
        const auto a = deductions::exact_wmbp_vp_child(parent, begin);
        const auto b = deductions::exact_wmbp_vp_child(parent, end);
        const auto offset = a.i - begin;
        if (a.ip == -1 && a.jp == -1 && b.ip == -1 && b.jp == -1 &&
            a.nonterminal == b.nonterminal && a.j == b.j && b.i - end == offset &&
            chart.for_each_live_primary_start_union(a.nonterminal, a.j,
                a.nonterminal, a.j, begin + offset, end + offset,
                [&](cand_pos_t start) { fn(start - offset); })) return;
      }
    }
    for (cand_pos_t k = begin; k <= end; ++k) fn(k);
  }

  template <class Semiring, class Chart, class Context>
  void prepare(const deductions::ItemKey &parent, const Chart &chart, const Context &ctx) {
    for (auto &s : selected) s.clear();
    using R = grammar::RuleId;
    if (parent.nonterminal == grammar::NonTerminal::WM) {
      if (!deductions::exact_wm_split_parent_allows(parent)) return;
      for_each_right_live_start<Semiring, R::WM_START_V, R::WM_START_WMB>(
          parent, chart, ctx, parent.i, parent.j, [&](cand_pos_t k) {
        if (deductions::exact_wm_start_split_allows(parent, k, ctx)) {
          if (children_nonzero<Semiring, R::WM_START_V>(parent, k, chart)) selected[0].push_back(k);
          if (children_nonzero<Semiring, R::WM_START_WMB>(parent, k, chart)) selected[1].push_back(k);
        }
        if (k > parent.i) {
          if (children_nonzero<Semiring, R::WM_SPLIT_V>(parent, k, chart)) selected[2].push_back(k);
          if (children_nonzero<Semiring, R::WM_SPLIT_WMB>(parent, k, chart)) selected[3].push_back(k);
        }
      });
    } else if (parent.nonterminal == grammar::NonTerminal::VM) {
      for_each_right_live_start<Semiring, R::VM_SPLIT_WM_WMv, R::VM_SPLIT_WMp_BASE>(
          parent, chart, ctx, deductions::exact_vm_base_begin_k(parent),
          deductions::exact_vm_end_k(parent, ctx), [&](cand_pos_t k) {
        if (k >= deductions::exact_vm_binary_begin_k(parent)) {
          if (children_nonzero<Semiring, R::VM_SPLIT_WM_WMv>(parent, k, chart)) selected[4].push_back(k);
          if (children_nonzero<Semiring, R::VM_SPLIT_WM_WMp>(parent, k, chart)) selected[5].push_back(k);
        }
        if (deductions::exact_vm_base_split_allows(parent, k, ctx) &&
            children_nonzero<Semiring, R::VM_SPLIT_WMp_BASE>(parent, k, chart)) selected[6].push_back(k);
      });
    } else if (parent.nonterminal == grammar::NonTerminal::WI) {
      if (!deductions::exact_round_weakly_closed_family_active(ctx, parent.i, parent.j)) return;
      for_each_right_live_start<Semiring, R::WI_SPLIT_V, R::WI_SPLIT_WMB>(
          parent, chart, ctx, parent.i, parent.j - ctx.turn() - 1, [&](cand_pos_t k) {
        if (ctx.can_pair(k, parent.j) &&
            children_nonzero<Semiring, R::WI_SPLIT_V>(parent, k, chart)) selected[7].push_back(k);
        if (children_nonzero<Semiring, R::WI_SPLIT_WMB>(parent, k, chart)) selected[8].push_back(k);
      });
    } else if (parent.nonterminal == grammar::NonTerminal::WIP) {
      if (!deductions::exact_round_weakly_closed_family_active(ctx, parent.i, parent.j)) return;
      for_each_right_live_start<Semiring, R::WIP_SPLIT_V, R::WIP_SPLIT_WMB>(
          parent, chart, ctx, parent.i + 1, parent.j - ctx.turn() - 1, [&](cand_pos_t k) {
        if (children_nonzero<Semiring, R::WIP_SPLIT_V>(parent, k, chart)) selected[9].push_back(k);
        if (children_nonzero<Semiring, R::WIP_SPLIT_WMB>(parent, k, chart)) selected[10].push_back(k);
        if (deductions::exact_wip_start_split_allows(parent, k, ctx)) {
          if (children_nonzero<Semiring, R::WIP_BASEPAIR_V>(parent, k, chart)) selected[11].push_back(k);
          if (children_nonzero<Semiring, R::WIP_BASEPAIR_WMB>(parent, k, chart)) selected[12].push_back(k);
        }
      });
    } else if (parent.nonterminal == grammar::NonTerminal::VPR) {
      if (!deductions::exact_round_open_family_active(ctx, parent.i, parent.j) ||
          !(ctx.pair_round(parent.j) < parent.j)) return;
      const cand_pos_t boundary = std::max(ctx.border_B(parent.i, parent.j), ctx.border_bp(parent.i, parent.j));
      for_each_left_live_end<Semiring, R::VPR_SPLIT_VP_WIP>(parent, chart, ctx,
          std::max(cand_pos_t(parent.i + 1), cand_pos_t(boundary + 1)),
          parent.j - 1, [&](cand_pos_t k) {
        if (children_nonzero<Semiring, R::VPR_SPLIT_VP_WIP>(parent, k, chart)) selected[13].push_back(k);
        if (ctx.can_pair_right_span(k, parent.j) &&
            children_nonzero<Semiring, R::VPR_SPLIT_VP_BASEPAIR>(parent, k, chart)) selected[14].push_back(k);
      });
    } else if (parent.nonterminal == grammar::NonTerminal::WMBP) {
      if (!deductions::exact_wmbp_item_active(parent, ctx)) return;
      for_each_wmbp_live_start<Semiring>(parent, chart, ctx,
          parent.i + 1, parent.j - 1, [&](cand_pos_t k) {
        const auto child = deductions::exact_wmbp_vp_child(parent, k);
        if (!Semiring::is_zero(chart.get(child))) selected[15].push_back(k);
      });
    }
  }
};

// A provider-owned buffer survives until the next item. The executor copies
// this view, not its candidates; evaluation still follows complete generation.
struct ExactDeductionBufferView {
  const std::vector<deductions::Deduction> &buffer;
  auto begin() const { return buffer.begin(); }
  auto end() const { return buffer.end(); }
};

template <class Chart>
typename Chart::value_type runtime_child_total_for_deduction(
    const Chart &chart,
    const deductions::Deduction &deduction,
    std::size_t child_idx);

template <class Semiring, bool SkipZeroChildren = false>
struct RuntimeProductContribution {
  template <class Chart, class LocalWeight>
  typename Semiring::value_type operator()(
      const deductions::Deduction &deduction, const Chart &chart,
      LocalWeight local) const {
    auto contribution = Semiring::from_weight(local);
    for (std::size_t c = 0; c < deduction.child_count; ++c)
      contribution = Semiring::times(contribution,
          runtime_child_total_for_deduction(chart, deduction, c));
    return contribution;
  }
};

template <class Semiring, class Chart>
void materialize_storage_only_exact_carriers_for_deduction(const Chart &chart,
                                                           const deductions::Deduction &deduction,
                                                           double local);

template <class Semiring>
class StorageOnlyExactCarrierObserver {
 public:
  void operator()(const deductions::ItemKey &,
                  const deductions::Deduction &deduction,
                  typename Semiring::value_type,
                  const engine::Chart<typename Semiring::value_type> &chart,
                  double local) const {
    materialize_storage_only_exact_carriers_for_deduction<Semiring>(
        chart, deduction, local);
  }
};

template <class Semiring,
          class WeightModel,
          class DeductionProviderFn,
          class DeductionAllowsFn,
          class CandidateObserverFn,
          class BeforeItemObserverFn,
          class AfterItemObserverFn,
          class ContributionFn>
void update_w_final_exact_dp_item_with_provider_and_contribution(
    const deductions::ItemKey &item,
    engine::Chart<typename Semiring::value_type> &chart,
    WeightModel &weight,
    DeductionProviderFn &deduction_provider,
    DeductionAllowsFn &deduction_allows,
    CandidateObserverFn &candidate_observer,
    BeforeItemObserverFn &before_item_observer,
    AfterItemObserverFn &after_item_observer,
    ContributionFn &contribution_for) {
  using value_type = typename Semiring::value_type;
  before_item_observer(item, chart);
  const auto candidates = [&] {
    if constexpr (std::is_invocable_v<DeductionProviderFn &,
        const deductions::ItemKey &>)
      return deduction_provider(item);
    else
      return deduction_provider(item, std::as_const(chart));
  }();
  value_type total = Semiring::zero();
  for (const auto &deduction : candidates) {
    engine::require_legacy_span_schedule_support(deduction);
    if (!deduction_allows(deduction)) {
      continue;
    }
    if constexpr (std::is_same_v<std::decay_t<ContributionFn>,
                                RuntimeProductContribution<Semiring, true>>) {
      // Only the canonical numeric product opts in. Custom/beam products can
      // reconstruct private aggregates and need not equal these child values.
      std::array<value_type, 3> children;
      bool dead = false;
      for (std::size_t c = 0; c < deduction.child_count; ++c) {
        children[c] = runtime_child_total_for_deduction(chart, deduction, c);
        if (Semiring::is_zero(children[c])) { dead = true; break; }
      }
      if (dead) continue;
      const auto local = weight.local_weight_for(deduction);
      if (weights::is_zero(local)) continue;
      auto contribution = Semiring::from_weight(local);
      // Preserve local, child0, child1, child2 arithmetic order exactly.
      for (std::size_t c = 0; c < deduction.child_count; ++c)
        contribution = Semiring::times(contribution, children[c]);
      total = Semiring::plus(total, contribution);
      candidate_observer(item, deduction, contribution, chart, local);
    } else {
      const auto local = weight.local_weight_for(deduction);
      if (weights::is_zero(local)) continue;
      const value_type contribution = contribution_for(deduction, chart, local);
      total = Semiring::plus(total, contribution);
      candidate_observer(item, deduction, contribution, chart, local);
    }
  }
  chart.set(item, total);
  after_item_observer(item, chart);
}

template <class Semiring,
          class Context,
          class DeductionProviderFn,
          class DeductionAllowsFn,
          class CandidateObserverFn,
          class BeforeItemObserverFn,
          class AfterItemObserverFn,
          class ContributionFn>
engine::Chart<typename Semiring::value_type>
run_w_final_exact_dp_schedule_with_provider_item_observers_and_contribution(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx,
    DeductionProviderFn &&deduction_provider,
    DeductionAllowsFn &&deduction_allows,
    CandidateObserverFn &&candidate_observer,
    BeforeItemObserverFn &&before_item_observer,
    AfterItemObserverFn &&after_item_observer,
    ContributionFn &&contribution_for,
    cand_pos_t dense_length = 0) {
  engine::AllLocalWeightModel<Context> weight(&ctx);
  engine::Chart<typename Semiring::value_type> chart(Semiring::zero(), nullptr, dense_length);
  for (const auto &item : schedule) {
    update_w_final_exact_dp_item_with_provider_and_contribution<Semiring>(
        item,
        chart,
        weight,
        deduction_provider,
        deduction_allows,
        candidate_observer,
        before_item_observer,
        after_item_observer,
        contribution_for);
  }
  return chart;
}

template <class Semiring,
          class Context,
          class DeductionAllowsFn,
          class CandidateObserverFn,
          class BeforeItemObserverFn,
          class AfterItemObserverFn,
          class ContributionFn>
engine::Chart<typename Semiring::value_type>
run_w_final_exact_dp_schedule_with_item_observers_and_contribution(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx,
    DeductionAllowsFn &&deduction_allows,
    CandidateObserverFn &&candidate_observer,
    BeforeItemObserverFn &&before_item_observer,
    AfterItemObserverFn &&after_item_observer,
    ContributionFn &&contribution_for,
    cand_pos_t dense_length = 0,
    bool omit_unused_vm = false,
    bool filter_zero_candidates = false,
    bool reuse_candidate_buffer = false) {
  ExactNumericSplitScratch scratch;
  const auto generate = [&](const deductions::ItemKey &item, const auto &chart,
                        std::vector<deductions::Deduction> &out) {
    if (omit_unused_vm && item.nonterminal == grammar::NonTerminal::VM &&
        !deductions::exact_v_family_active(ctx, item.i, item.j))
      return;
    if (filter_zero_candidates && ExactNumericSplitScratch::supports_item(item.nonterminal)) {
      if constexpr (std::is_same_v<Semiring, semirings::LogSumExpSemiring> ||
                    std::is_same_v<Semiring, semirings::LogMaxSemiring>) {
        scratch.prepare<Semiring>(item, chart, ctx);
        const auto views = scratch.views();
        deductions::append_exact_basic_deductions_for_item(out, item, ctx, nullptr, &views);
        return;
      } else {
        throw std::invalid_argument("zero-candidate selection requires a numeric log semiring");
      }
    }
    deductions::append_exact_basic_deductions_for_item(out, item, ctx);
  };
  const auto run = [&](auto &&provider) {
    return run_w_final_exact_dp_schedule_with_provider_item_observers_and_contribution<Semiring>(
      schedule, ctx, provider,
      std::forward<DeductionAllowsFn>(deduction_allows),
      std::forward<CandidateObserverFn>(candidate_observer),
      std::forward<BeforeItemObserverFn>(before_item_observer),
      std::forward<AfterItemObserverFn>(after_item_observer),
      std::forward<ContributionFn>(contribution_for), dense_length);
  };
  if (reuse_candidate_buffer) {
    std::vector<deductions::Deduction> buffer;
    return run([&](const deductions::ItemKey &item, const auto &chart) {
      buffer.clear();
      generate(item, chart, buffer);
      return ExactDeductionBufferView{buffer};
    });
  }
  return run([&](const deductions::ItemKey &item, const auto &chart) {
    std::vector<deductions::Deduction> out;
    generate(item, chart, out);
    return out;
  });
}

template <class Semiring,
          class Context,
          class DeductionAllowsFn,
          class CandidateObserverFn,
          class ItemObserverFn,
          class ContributionFn>
engine::Chart<typename Semiring::value_type>
run_w_final_exact_dp_schedule_with_item_and_contribution(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx,
    DeductionAllowsFn &&deduction_allows,
    CandidateObserverFn &&candidate_observer,
    ItemObserverFn &&item_observer,
    ContributionFn &&contribution_for) {
  return run_w_final_exact_dp_schedule_with_item_observers_and_contribution<
      Semiring>(
      schedule,
      ctx,
      std::forward<DeductionAllowsFn>(deduction_allows),
      std::forward<CandidateObserverFn>(candidate_observer),
      std::forward<ItemObserverFn>(item_observer),
      [](const deductions::ItemKey &, const auto &) {},
      std::forward<ContributionFn>(contribution_for));
}

template <class Semiring,
          class Context,
          class DeductionAllowsFn,
          class CandidateObserverFn,
          class ContributionFn>
engine::Chart<typename Semiring::value_type>
run_w_final_exact_dp_schedule_with_contribution(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx,
    DeductionAllowsFn &&deduction_allows,
    CandidateObserverFn &&candidate_observer,
    ContributionFn &&contribution_for) {
  return run_w_final_exact_dp_schedule_with_item_and_contribution<Semiring>(
      schedule,
      ctx,
      std::forward<DeductionAllowsFn>(deduction_allows),
      std::forward<CandidateObserverFn>(candidate_observer),
      [](const deductions::ItemKey &, const auto &) {},
      std::forward<ContributionFn>(contribution_for));
}

template <class Semiring, class Chart, class LocalWeight>
typename Semiring::value_type runtime_deduction_contribution(
    const Chart &chart,
    const deductions::Deduction &deduction,
    LocalWeight local) {
  return RuntimeProductContribution<Semiring>{}(deduction, chart, local);
}

template <class Semiring,
          class Context,
          class DeductionAllowsFn,
          class CandidateObserverFn>
engine::Chart<typename Semiring::value_type> run_w_final_exact_dp_schedule(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx,
    DeductionAllowsFn &&deduction_allows,
    CandidateObserverFn &&candidate_observer,
    ExactExecutionOptions options = {}) {
  if constexpr (weights::LogContext<Context>::value &&
      (std::is_same_v<Semiring, semirings::LogSumExpSemiring> ||
       std::is_same_v<Semiring, semirings::LogMaxSemiring>)) {
    if (options.skip_zero_children) {
      return run_w_final_exact_dp_schedule_with_item_observers_and_contribution<Semiring>(
          schedule, ctx, std::forward<DeductionAllowsFn>(deduction_allows),
          std::forward<CandidateObserverFn>(candidate_observer),
          [](const auto &, const auto &) {}, [](const auto &, const auto &) {},
          RuntimeProductContribution<Semiring, true>{},
          options.dense_primary ? ctx.sequence_length() : 0, options.omit_unused_vm,
          options.filter_zero_candidates, options.reuse_candidate_buffer);
    }
  } else {
    if (options.skip_zero_children || options.omit_unused_vm || options.dense_primary ||
        options.filter_zero_candidates || options.reuse_candidate_buffer)
      throw std::invalid_argument("numeric exact execution requires a numeric log semiring");
  }
  return run_w_final_exact_dp_schedule_with_item_observers_and_contribution<Semiring>(
      schedule,
      ctx,
      std::forward<DeductionAllowsFn>(deduction_allows),
      std::forward<CandidateObserverFn>(candidate_observer),
      [](const auto &, const auto &) {}, [](const auto &, const auto &) {},
      RuntimeProductContribution<Semiring>{},
      options.dense_primary ? ctx.sequence_length() : 0, options.omit_unused_vm,
      options.filter_zero_candidates, options.reuse_candidate_buffer);
}

template <class Semiring = semirings::SumProductSemiring>
engine::Chart<typename Semiring::value_type> run_w_final_exact_inside_schedule(
    weights::RuntimeSnapshot &owner,
    const std::vector<deductions::ItemKey> &schedule) {
  WFinalExactInsideContext ctx(owner);
  return run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      [](const deductions::Deduction &) { return true; },
      [](const deductions::ItemKey &,
         const deductions::Deduction &,
         typename Semiring::value_type,
         const engine::Chart<typename Semiring::value_type> &,
         double) {});
}

template <class Semiring = semirings::MaxProductSemiring, class Context, class DeductionAllowsFn>
engine::ViterbiRunResult<double> run_w_final_exact_viterbi_schedule_with_context(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx,
    DeductionAllowsFn &&deduction_allows,
    ExactExecutionOptions options = {}) {
  engine::ChoiceMap<double> choices;
  StorageOnlyExactCarrierObserver<Semiring> carrier_observer;
  auto chart = run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      std::forward<DeductionAllowsFn>(deduction_allows),
      [&](const deductions::ItemKey &item,
          const deductions::Deduction &deduction,
          double contribution,
          const engine::Chart<double> &current_chart,
          auto local) {
        if (contribution != Semiring::zero()) {
          choices.observe(item,
                          deduction,
                          contribution,
                          [](double lhs, double rhs) { return lhs > rhs; });
        }
        if constexpr (!weights::LogContext<Context>::value)
          carrier_observer(item, deduction, contribution, current_chart, local);
      }, options);
  return {std::move(chart), std::move(choices)};
}

template <bool LogDomain = false>
engine::ViterbiTracebackResult<double> run_w_final_exact_viterbi_traceback(
    weights::RuntimeSnapshot &owner,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure) {
  WFinalInsideContext<LogDomain> ctx(owner, policy);
  using S = std::conditional_t<LogDomain, semirings::LogMaxSemiring, semirings::MaxProductSemiring>;
  if constexpr (!LogDomain) {
    storage::exact_store::reset();
  }
  auto run = run_w_final_exact_viterbi_schedule_with_context<S>(
      engine::build_pkfree_basic_exact_schedule(static_cast<cand_pos_t>(owner.length()), ctx),
      ctx, [](const auto &) { return true; },
      LogDomain ? ExactExecutionOptions::numeric() : ExactExecutionOptions{});
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
  return engine::build_viterbi_traceback_result(root, std::move(run), ctx);
}

template <bool LogDomain = false>
engine::ViterbiTracebackResult<double>
run_w_final_exact_output_constrained_viterbi_traceback(
    weights::RuntimeSnapshot &owner,
    const std::vector<engine::TracebackPair> &target_pairs,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure) {
  WFinalInsideContext<LogDomain> ctx(owner, policy);
  using S = std::conditional_t<LogDomain, semirings::LogMaxSemiring, semirings::MaxProductSemiring>;
  const auto generated_target = engine::generated_traceback_pairs(target_pairs, ctx);
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};

  const auto has_exact_generated_output = [&](const auto &result) {
    return engine::same_traceback_pair_set(
        result.generated_pairs, generated_target);
  };

  const auto run_with = [&](auto &&allows) {
    if constexpr (!LogDomain) {
      storage::exact_store::reset();
    }
    auto run = run_w_final_exact_viterbi_schedule_with_context<S>(
        schedule, ctx, std::forward<decltype(allows)>(allows),
        LogDomain ? ExactExecutionOptions::numeric() : ExactExecutionOptions{});
    return engine::build_viterbi_traceback_result(root, std::move(run), ctx);
  };

  auto result = run_with([&](const deductions::Deduction &deduction) {
    return engine::deduction_locally_matches_generated_target(
        deduction, ctx, generated_target);
  });
  if (has_exact_generated_output(result)) {
    return result;
  }

  result.run.chart.set(root, S::zero());
  result.trace = engine::DeductionTraceNode<double>{
      .item = root,
      .chart_value = S::zero(),
  };
  result.pairs.clear();
  result.lane_structures = engine::render_traceback_pairs(n, result.pairs);
  return result;
}

template <class Chart>
double chart_child_product(const Chart &chart, const deductions::Deduction &deduction) {
  double out = 1.0;
  for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
    out *= chart.get(deduction.children[idx].item);
  }
  return out;
}

template <class Chart>
typename Chart::value_type runtime_child_total_for_deduction(
    const Chart &chart,
    const deductions::Deduction &deduction,
    std::size_t child_idx) {
  return chart.get(engine::RuntimeDeductionChildItem{}(deduction, child_idx));
}

template <class Chart>
double chart_child_product_except(const Chart &chart,
                                  const deductions::Deduction &deduction,
                                  std::size_t skip_idx) {
  double out = 1.0;
  for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
    if (idx == skip_idx) {
      continue;
    }
    out *= chart.get(deduction.children[idx].item);
  }
  return out;
}

inline std::optional<std::size_t> find_child_index(const deductions::Deduction &deduction,
                                                   grammar::NonTerminal nt) {
  for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
    if (deduction.children[idx].item.nonterminal == nt) {
      return idx;
    }
  }
  return std::nullopt;
}

template <class Semiring>
void accumulate_exact_item_total(const storage::ExactItemKey &item, double candidate) {
  if (Semiring::is_zero(candidate)) {
    return;
  }
  storage::exact_store::set_total(
      item,
      Semiring::plus(storage::exact_store::total(item), candidate));
}

inline storage::ExactItemKey make_vp_case3_left_exact_item_key(
    const deductions::Deduction &deduction) {
  storage::ExactItemKey item;
  item.nt = storage::NonTerminal::VP_CASE3_LEFT;
  item.i = deduction.parent.i;
  item.j = deduction.parent.j;
  item.aux0 = deduction.split.p;
  item.aux1 = deduction.split.q;
  item.aux2 = deduction.split.k;
  item.aux3 = deduction.split.l;
  return item;
}

template <class Semiring>
void propagate_case3_left_items_to_parent(cand_pos_t parent_i,
                                          cand_pos_t parent_j,
                                          cand_pos_t child_i,
                                          cand_pos_t child_j,
                                          double wrapper_coeff) {
  if (wrapper_coeff == 0.0) {
    return;
  }
  storage::exact_store::for_each_item_span(
      storage::NonTerminal::VP_CASE3_LEFT,
      child_i,
      child_j,
      [&](const storage::ExactItemKey &child_item, pf_t child_total) {
        storage::ExactItemKey parent_item = child_item;
        parent_item.i = parent_i;
        parent_item.j = parent_j;
        accumulate_exact_item_total<Semiring>(parent_item, child_total * wrapper_coeff);
      });
}

inline std::optional<storage::ExactItemKey> make_wmb_bstep_exact_item_key(
    const deductions::Deduction &deduction) {
  const grammar::WmbpStepPrefixKind prefix_kind =
      grammar::classify_wmbp_step_prefix_kind(deduction.rule);
  if (prefix_kind == grammar::WmbpStepPrefixKind::Unknown) {
    return std::nullopt;
  }

  storage::ExactItemKey item;
  item.nt = storage::NonTerminal::WMB_BSTEP;
  item.i = deduction.parent.i;
  item.j = deduction.parent.j;
  item.aux0 = deduction.split.k;
  item.aux1 = static_cast<cand_pos_t>(prefix_kind);
  item.aux2 = -1;
  if (prefix_kind == grammar::WmbpStepPrefixKind::Empty) {
    item.aux0 = deduction.parent.i;
  } else if (prefix_kind == grammar::WmbpStepPrefixKind::PrefixWi) {
    item.aux2 = deduction.split.p;
  } else if (prefix_kind == grammar::WmbpStepPrefixKind::PrefixTail) {
    item.aux2 = deduction.split.l;
  }
  return item;
}

template <class Semiring>
void propagate_wmb_bstep_items_to_parent(const storage::ExactItemKey &parent_item,
                                         cand_pos_t child_i,
                                         cand_pos_t child_j,
                                         double wrapper_coeff) {
  if (wrapper_coeff == 0.0) {
    return;
  }
  storage::exact_store::for_each_item_span(
      storage::NonTerminal::VP_CASE3_LEFT,
      child_i,
      child_j,
      [&](const storage::ExactItemKey &, pf_t child_total) {
        accumulate_exact_item_total<Semiring>(parent_item, child_total * wrapper_coeff);
      });
}

inline storage::ExactItemKey make_wmb_nonempty_exact_item_key(
    const deductions::Deduction &deduction,
    const storage::ExactItemKey &child_item) {
  storage::ExactItemKey item;
  item.nt = storage::NonTerminal::WMB_NONEMPTY;
  item.i = deduction.parent.i;
  item.j = deduction.parent.j;
  item.aux0 = static_cast<cand_pos_t>(deduction.rule);
  item.aux1 = child_item.aux0;
  item.aux2 = child_item.aux1;
  item.aux3 = child_item.aux2;
  if (deduction.rule == grammar::RuleId::WMB_SPLIT_BE_WMBP_WI) {
    item.ip = deduction.split.k;
    item.jp = deduction.split.q;
  }
  return item;
}

inline storage::ExactItemKey relay_exact_item_to_parent(storage::NonTerminal parent_nt,
                                                         cand_pos_t parent_i,
                                                         cand_pos_t parent_j,
                                                         const storage::ExactItemKey &child_item) {
  storage::ExactItemKey item = child_item;
  item.nt = parent_nt;
  item.i = parent_i;
  item.j = parent_j;
  return item;
}

template <class Semiring, class Chart>
void materialize_storage_only_exact_carriers_for_deduction(const Chart &chart,
                                                           const deductions::Deduction &deduction,
                                                           double local) {
  switch (deduction.rule) {
    case grammar::RuleId::VP_WI_CASE3:
      if (deduction.parent.nonterminal == grammar::NonTerminal::VP) {
        if (deduction.split.l >= deduction.split.p && deduction.split.p >= 0) {
          const double contribution = local * chart_child_product(chart, deduction);
          accumulate_exact_item_total<Semiring>(
              make_vp_case3_left_exact_item_key(deduction),
              contribution);
        }
      }
      break;

    case grammar::RuleId::VP_STACK:
    case grammar::RuleId::VP_INTERNAL_LOOP:
      if (deduction.parent.nonterminal == grammar::NonTerminal::VP && deduction.child_count == 1) {
        const auto &child = deduction.children[0].item;
        propagate_case3_left_items_to_parent<Semiring>(
            deduction.parent.i, deduction.parent.j, child.i, child.j, local);
      }
      break;

    case grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP:
    case grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP:
    case grammar::RuleId::WMBP_DIRECT_VP:
    case grammar::RuleId::WMBP_SPLIT_BE_WI_VP: {
      const auto vp_idx = find_child_index(deduction, grammar::NonTerminal::VP);
      const auto parent_item = make_wmb_bstep_exact_item_key(deduction);
      if (!vp_idx.has_value() || !parent_item.has_value()) {
        break;
      }
      const auto &vp_child = deduction.children[*vp_idx].item;
      const double wrapper_coeff = local * chart_child_product_except(chart, deduction, *vp_idx);
      propagate_wmb_bstep_items_to_parent<Semiring>(
          *parent_item, vp_child.i, vp_child.j, wrapper_coeff);
      break;
    }

    case grammar::RuleId::WMB_DIRECT_WMBP:
    case grammar::RuleId::WMB_SPLIT_BE_WMBP_WI: {
      const auto wmbp_idx = find_child_index(deduction, grammar::NonTerminal::WMBP);
      if (!wmbp_idx.has_value()) {
        break;
      }
      const auto &wmbp_child = deduction.children[*wmbp_idx].item;
      const double wrapper_coeff = local * chart_child_product_except(chart, deduction, *wmbp_idx);
      storage::exact_store::for_each_item_span(
          storage::NonTerminal::WMB_BSTEP,
          wmbp_child.i,
          wmbp_child.j,
          [&](const storage::ExactItemKey &child_item, pf_t child_total) {
            accumulate_exact_item_total<Semiring>(
                make_wmb_nonempty_exact_item_key(deduction, child_item),
                child_total * wrapper_coeff);
          });
      break;
    }

    case grammar::RuleId::WM_START_WMB:
    case grammar::RuleId::WM_SPLIT_WMB: {
      const auto wmb_idx = find_child_index(deduction, grammar::NonTerminal::WMB);
      if (!wmb_idx.has_value()) {
        break;
      }
      const auto &wmb_child = deduction.children[*wmb_idx].item;
      const double wrapper_coeff = local * chart_child_product_except(chart, deduction, *wmb_idx);
      storage::exact_store::for_each_item_span(
          storage::NonTerminal::WMB_NONEMPTY,
          wmb_child.i,
          wmb_child.j,
          [&](const storage::ExactItemKey &child_item, pf_t child_total) {
            accumulate_exact_item_total<Semiring>(
                relay_exact_item_to_parent(
                    storage::NonTerminal::WM, deduction.parent.i, deduction.parent.j, child_item),
                child_total * wrapper_coeff);
          });
      break;
    }

    case grammar::RuleId::WMp_STEM_WMB: {
      const auto wmb_idx = find_child_index(deduction, grammar::NonTerminal::WMB);
      if (!wmb_idx.has_value()) {
        break;
      }
      const auto &wmb_child = deduction.children[*wmb_idx].item;
      const double wrapper_coeff = local * chart_child_product_except(chart, deduction, *wmb_idx);
      storage::exact_store::for_each_item_span(
          storage::NonTerminal::WMB_NONEMPTY,
          wmb_child.i,
          wmb_child.j,
          [&](const storage::ExactItemKey &child_item, pf_t child_total) {
            accumulate_exact_item_total<Semiring>(
                relay_exact_item_to_parent(
                    storage::NonTerminal::WMp, deduction.parent.i, deduction.parent.j, child_item),
                child_total * wrapper_coeff);
          });
      break;
    }

    case grammar::RuleId::WMp_EXTEND_UNPAIRED: {
      const auto wmp_idx = find_child_index(deduction, grammar::NonTerminal::WMp);
      if (!wmp_idx.has_value()) {
        break;
      }
      const auto &wmp_child = deduction.children[*wmp_idx].item;
      const double wrapper_coeff = local * chart_child_product_except(chart, deduction, *wmp_idx);
      storage::exact_store::for_each_item_span(
          storage::NonTerminal::WMp,
          wmp_child.i,
          wmp_child.j,
          [&](const storage::ExactItemKey &child_item, pf_t child_total) {
            accumulate_exact_item_total<Semiring>(
                relay_exact_item_to_parent(
                    storage::NonTerminal::WMp, deduction.parent.i, deduction.parent.j, child_item),
                child_total * wrapper_coeff);
          });
      break;
    }

    default:
      break;
  }
}

}  // namespace cparty::scfg2::replay

#endif  // CPARTY_SCFG2_REPLAY_W_FINAL_EXACT_INSIDE_HH

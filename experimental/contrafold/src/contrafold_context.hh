#ifndef PK_CONTRAFOLD_CONTEXT_HH
#define PK_CONTRAFOLD_CONTEXT_HH

#include "scfg2/engine/runtime_child_item.hh"
#include "contra_scores.hh"

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
class ContrafoldContext final : public engine::BasicExactWeightContext,
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

  explicit ContrafoldContext(weights::RuntimeSnapshot &owner,
      WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure)
      : owner_(owner), wmb_closure_policy_(policy),
        tree_(owner_.scaffold(), owner_.length()) {}

  double exterior_unpaired_weight() const override {
    return project_log_factor(pkfit::external_unpaired_score());
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
    if (i < 1 || j < 1 || j-i < 4) {
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
    return project_log_factor(length * pkfit::multi_unpaired_score());
  }
  double cp_factor(cand_pos_t length) const {
    return project_log_factor(length == 0 ? 0. : length * owner_.log_pk_parameters().cp);
  }
  double pup_factor(cand_pos_t length) const {
    return project_log_factor(length == 0 ? 0. : length * owner_.log_pk_parameters().pup);
  }
  double pk_factor(double weights::RuntimeSnapshot::LogPkParameters::*member) const {
#ifdef PK_FIXED_CF_BRANCH
    if(member == &weights::RuntimeSnapshot::LogPkParameters::b)
      return project_log_factor(pkfit::multi_paired_score());
#endif
    return project_log_factor(owner_.log_pk_parameters().*member);
  }
  double mlclosing_factor() const {
    return project_log_factor(pkfit::multi_base_score());
  }
  double hairpin_factor(cand_pos_t i, cand_pos_t j) const {
    return project_log_factor(pkfit::hairpin(owner_.sequence(),i,j));
  }
  double internal_log_factor(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
    return pkfit::internal(owner_.sequence(),i,j,k,l);
  }
  double internal_factor(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) const {
    return project_log_factor(internal_log_factor(i,j,k,l));
  }
  double stack_pk_factor(cand_pos_t i, cand_pos_t j) const {
    const double value=internal_log_factor(i,j,i+1,j-1);
    return std::isfinite(value) ? project_log_factor(owner_.pseudoknot_model().stack_multiplier*value) : impossible();
  }
  double internal_pk_factor(cand_pos_t i,cand_pos_t k,cand_pos_t l,cand_pos_t j) const {
    if(k==i+1 && l==j-1)return impossible();
    const double value=internal_log_factor(i,j,k,l);
    return std::isfinite(value) ? project_log_factor(owner_.pseudoknot_model().internal_multiplier*value) : impossible();
  }
  double stem_factor(cand_pos_t i,cand_pos_t j,bool multiloop,bool closing) const {
    return project_log_factor(pkfit::stem(owner_.sequence(),i,j,multiloop,closing));
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


} // namespace cparty::scfg2::replay
#endif

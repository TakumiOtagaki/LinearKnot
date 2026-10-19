#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

int main(int argc, char **argv) {
  if (argc < 4 || argc > 7) {
    std::cerr << "usage: " << argv[0]
              << " SEQUENCE SCAFFOLD WIDTH [inside|max] [PARAM_FILE]"
                 " [REPETITIONS]\n";
    return EXIT_FAILURE;
  }
  const std::string sequence = argv[1];
  const std::string scaffold = argv[2];
  const std::size_t width = std::stoull(argv[3]);
  const std::string mode = argc >= 5 ? argv[4] : "inside";
  const std::size_t repetitions = argc == 7 ? std::stoull(argv[6]) : 1;
  if (sequence.size() != scaffold.size() ||
      (mode != "inside" && mode != "max") || repetitions == 0) {
    std::cerr << "invalid sequence/scaffold length or mode\n";
    return EXIT_FAILURE;
  }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc >= 6 ? argv[5] : "params/rna_Turner04.par",
         1);

  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto widths =
      cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(width);
  cparty::scfg2::replay::WFinalOnlineNonterminalCounts visits{};
  cparty::scfg2::replay::WFinalOnlineNonterminalCounts positive{};
  double root = 0.0;
  std::size_t wmbp_parent_visits = 0;
  std::size_t vm_binary_join_probes = 0;
  std::size_t vm_base_prefix_cell_visits = 0;
  std::size_t vm_base_prefix_selected = 0;
  std::size_t wm_binary_join_probes = 0;
  std::size_t wm_extension_probes = 0;
  std::size_t wm_start_prefix_cell_visits = 0;
  std::size_t wm_start_prefix_selected = 0;
  std::size_t wi_binary_join_probes = 0;
  std::size_t wi_direct_probes = 0;
  std::size_t wi_extension_probes = 0;
  std::size_t wip_binary_join_probes = 0;
  std::size_t wip_extension_probes = 0;
  std::size_t wip_start_prefix_cell_visits = 0;
  std::size_t wip_start_prefix_selected = 0;
  std::size_t wmbw_binary_join_probes = 0;
  std::size_t vpl_prefix_cell_visits = 0;
  std::size_t vpl_prefix_selected = 0;
  std::size_t vpr_binary_join_probes = 0;
  std::size_t vpr_base_state_probes = 0;
  std::size_t vpr_base_state_selected = 0;
  std::size_t w_binary_join_probes = 0;
  std::size_t w_direct_probes = 0;
  std::size_t w_extension_probes = 0;
  std::size_t w_root_decomposition_probes = 0;
  std::size_t wmb_split_positive_sources = 0;
  std::size_t wmb_direct_positive_sources = 0;
  std::size_t wmb_empty_positive_sources = 0;
  std::size_t wmb_btail_join_probes = 0;
  std::size_t wmb_btail_direct_probes = 0;
  std::size_t wmb_btail_selected = 0;
  std::size_t wmb_be_btail_join_probes = 0;
  std::size_t wmb_retained_direct_probes = 0;
  std::size_t wmb_retained_split_selected = 0;
  std::size_t wmbp_shared_wmbp_left_probes = 0;
  std::size_t wmbp_shared_wmbw_left_probes = 0;
  std::size_t wmbp_shared_emitted_sources = 0;
  std::size_t wmbp_wi_history_insertions = 0;
  std::size_t wmbp_wi_be_wi_probes = 0;
  std::size_t wmbp_wi_empty_probes = 0;
  std::size_t wmbp_wi_prefix_candidates = 0;
  std::size_t wmbp_wi_prefix_vp_probes = 0;
  std::size_t wmbp_wi_direct_probes = 0;
  std::size_t wmbp_wi_emitted_sources = 0;
  std::size_t vp_loop_parent_probes = 0;
  std::size_t vp_stack_event_sources = 0;
  std::size_t vp_internal_event_sources = 0;
  std::size_t case3_active_expiry_probes = 0;
  std::size_t case3_tail_materialization_probes = 0;
  std::size_t case3_prepared_score_probes = 0;
  std::size_t case3_max_active_segments = 0;
  std::size_t v_hairpin_parent_probes = 0;
  std::size_t v_internal_parent_probes = 0;
  std::size_t v_vm_parent_probes = 0;
  std::size_t v_internal_sources = 0;
  std::size_t case3_heap_registered_segments = 0;
  std::size_t case3_heap_expired_segments = 0;
  std::size_t case3_max_heap_size = 0;
  double elapsed_sum_ms = 0.0;
  double elapsed_min_ms = std::numeric_limits<double>::infinity();
  for (std::size_t repetition = 0; repetition < repetitions; ++repetition) {
    const auto started = std::chrono::steady_clock::now();
    if (mode == "inside") {
      auto run = cparty::scfg2::replay::run_w_final_vp_wmbp_online<
          cparty::scfg2::semirings::SumProductSemiring>(
          ctx, static_cast<cand_pos_t>(sequence.size()), widths);
      visits = run.item_visits;
      positive = run.positive_items;
      wmbp_parent_visits = run.schedule.wmbp_parent_visits;
      vm_binary_join_probes = run.vm_binary_join_probes;
      vm_base_prefix_cell_visits = run.vm_base_prefix_cell_visits;
      vm_base_prefix_selected = run.vm_base_prefix_selected;
      wm_binary_join_probes = run.wm_binary_join_probes;
      wm_extension_probes = run.wm_extension_probes;
      wm_start_prefix_cell_visits = run.wm_start_prefix_cell_visits;
      wm_start_prefix_selected = run.wm_start_prefix_selected;
      wi_binary_join_probes = run.wi_binary_join_probes;
      wi_direct_probes = run.wi_direct_probes;
      wi_extension_probes = run.wi_extension_probes;
      wip_binary_join_probes = run.wip_binary_join_probes;
      wip_extension_probes = run.wip_extension_probes;
      wip_start_prefix_cell_visits = run.wip_start_prefix_cell_visits;
      wip_start_prefix_selected = run.wip_start_prefix_selected;
      wmbw_binary_join_probes = run.wmbw_binary_join_probes;
      vpl_prefix_cell_visits = run.vpl_prefix_cell_visits;
      vpl_prefix_selected = run.vpl_prefix_selected;
      vpr_binary_join_probes = run.vpr_binary_join_probes;
      vpr_base_state_probes = run.vpr_base_state_probes;
      vpr_base_state_selected = run.vpr_base_state_selected;
      w_binary_join_probes = run.w_binary_join_probes;
      w_direct_probes = run.w_direct_probes;
      w_extension_probes = run.w_extension_probes;
      w_root_decomposition_probes = run.w_root_decomposition_probes;
      wmb_split_positive_sources = run.wmb_split_positive_sources;
      wmb_direct_positive_sources = run.wmb_direct_positive_sources;
      wmb_empty_positive_sources = run.wmb_empty_positive_sources;
      wmb_btail_join_probes = run.wmb_btail_join_probes;
      wmb_btail_direct_probes = run.wmb_btail_direct_probes;
      wmb_btail_selected = run.wmb_btail_selected;
      wmb_be_btail_join_probes = run.wmb_be_btail_join_probes;
      wmb_retained_direct_probes = run.wmb_retained_direct_probes;
      wmb_retained_split_selected = run.wmb_retained_split_selected;
      wmbp_shared_wmbp_left_probes =
          run.wmbp_shared_wmbp_left_probes;
      wmbp_shared_wmbw_left_probes =
          run.wmbp_shared_wmbw_left_probes;
      wmbp_shared_emitted_sources = run.wmbp_shared_emitted_sources;
      wmbp_wi_history_insertions = run.wmbp_wi_history_insertions;
      wmbp_wi_be_wi_probes = run.wmbp_wi_be_wi_probes;
      wmbp_wi_empty_probes = run.wmbp_wi_empty_probes;
      wmbp_wi_prefix_candidates = run.wmbp_wi_prefix_candidates;
      wmbp_wi_prefix_vp_probes = run.wmbp_wi_prefix_vp_probes;
      wmbp_wi_direct_probes = run.wmbp_wi_direct_probes;
      wmbp_wi_emitted_sources = run.wmbp_wi_emitted_sources;
      vp_loop_parent_probes = run.vp_loop_parent_probes;
      vp_stack_event_sources = run.vp_stack_event_sources;
      vp_internal_event_sources = run.vp_internal_event_sources;
      case3_active_expiry_probes = run.case3_active_expiry_probes;
      case3_tail_materialization_probes =
          run.case3_tail_materialization_probes;
      case3_prepared_score_probes = run.case3_prepared_score_probes;
      case3_max_active_segments = run.case3_max_active_segments;
      v_hairpin_parent_probes = run.v_hairpin_parent_probes;
      v_internal_parent_probes = run.v_internal_parent_probes;
      v_vm_parent_probes = run.v_vm_parent_probes;
      v_internal_sources = run.v_internal_sources;
      case3_heap_registered_segments = run.case3_heap_registered_segments;
      case3_heap_expired_segments = run.case3_heap_expired_segments;
      case3_max_heap_size = run.case3_max_heap_size;
      root = run.chart.get(cparty::scfg2::deductions::ItemKey{
          cparty::scfg2::grammar::NonTerminal::W,
          1,
          static_cast<cand_pos_t>(sequence.size()),
          -1,
          -1});
    } else {
      auto run = cparty::scfg2::replay::run_w_final_vp_wmbp_online<
          cparty::scfg2::semirings::MaxProductSemiring>(
          ctx, static_cast<cand_pos_t>(sequence.size()), widths);
      visits = run.item_visits;
      positive = run.positive_items;
      wmbp_parent_visits = run.schedule.wmbp_parent_visits;
      vm_binary_join_probes = run.vm_binary_join_probes;
      vm_base_prefix_cell_visits = run.vm_base_prefix_cell_visits;
      vm_base_prefix_selected = run.vm_base_prefix_selected;
      wm_binary_join_probes = run.wm_binary_join_probes;
      wm_extension_probes = run.wm_extension_probes;
      wm_start_prefix_cell_visits = run.wm_start_prefix_cell_visits;
      wm_start_prefix_selected = run.wm_start_prefix_selected;
      wi_binary_join_probes = run.wi_binary_join_probes;
      wi_direct_probes = run.wi_direct_probes;
      wi_extension_probes = run.wi_extension_probes;
      wip_binary_join_probes = run.wip_binary_join_probes;
      wip_extension_probes = run.wip_extension_probes;
      wip_start_prefix_cell_visits = run.wip_start_prefix_cell_visits;
      wip_start_prefix_selected = run.wip_start_prefix_selected;
      wmbw_binary_join_probes = run.wmbw_binary_join_probes;
      vpl_prefix_cell_visits = run.vpl_prefix_cell_visits;
      vpl_prefix_selected = run.vpl_prefix_selected;
      vpr_binary_join_probes = run.vpr_binary_join_probes;
      vpr_base_state_probes = run.vpr_base_state_probes;
      vpr_base_state_selected = run.vpr_base_state_selected;
      w_binary_join_probes = run.w_binary_join_probes;
      w_direct_probes = run.w_direct_probes;
      w_extension_probes = run.w_extension_probes;
      w_root_decomposition_probes = run.w_root_decomposition_probes;
      wmb_split_positive_sources = run.wmb_split_positive_sources;
      wmb_direct_positive_sources = run.wmb_direct_positive_sources;
      wmb_empty_positive_sources = run.wmb_empty_positive_sources;
      wmb_btail_join_probes = run.wmb_btail_join_probes;
      wmb_btail_direct_probes = run.wmb_btail_direct_probes;
      wmb_btail_selected = run.wmb_btail_selected;
      wmb_be_btail_join_probes = run.wmb_be_btail_join_probes;
      wmb_retained_direct_probes = run.wmb_retained_direct_probes;
      wmb_retained_split_selected = run.wmb_retained_split_selected;
      wmbp_shared_wmbp_left_probes =
          run.wmbp_shared_wmbp_left_probes;
      wmbp_shared_wmbw_left_probes =
          run.wmbp_shared_wmbw_left_probes;
      wmbp_shared_emitted_sources = run.wmbp_shared_emitted_sources;
      wmbp_wi_history_insertions = run.wmbp_wi_history_insertions;
      wmbp_wi_be_wi_probes = run.wmbp_wi_be_wi_probes;
      wmbp_wi_empty_probes = run.wmbp_wi_empty_probes;
      wmbp_wi_prefix_candidates = run.wmbp_wi_prefix_candidates;
      wmbp_wi_prefix_vp_probes = run.wmbp_wi_prefix_vp_probes;
      wmbp_wi_direct_probes = run.wmbp_wi_direct_probes;
      wmbp_wi_emitted_sources = run.wmbp_wi_emitted_sources;
      vp_loop_parent_probes = run.vp_loop_parent_probes;
      vp_stack_event_sources = run.vp_stack_event_sources;
      vp_internal_event_sources = run.vp_internal_event_sources;
      case3_active_expiry_probes = run.case3_active_expiry_probes;
      case3_tail_materialization_probes =
          run.case3_tail_materialization_probes;
      case3_prepared_score_probes = run.case3_prepared_score_probes;
      case3_max_active_segments = run.case3_max_active_segments;
      v_hairpin_parent_probes = run.v_hairpin_parent_probes;
      v_internal_parent_probes = run.v_internal_parent_probes;
      v_vm_parent_probes = run.v_vm_parent_probes;
      v_internal_sources = run.v_internal_sources;
      case3_heap_registered_segments = run.case3_heap_registered_segments;
      case3_heap_expired_segments = run.case3_heap_expired_segments;
      case3_max_heap_size = run.case3_max_heap_size;
      root = run.chart.get(cparty::scfg2::deductions::ItemKey{
          cparty::scfg2::grammar::NonTerminal::W,
          1,
          static_cast<cand_pos_t>(sequence.size()),
          -1,
          -1});
    }
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - started)
            .count();
    elapsed_sum_ms += elapsed_ms;
    elapsed_min_ms = std::min(elapsed_min_ms, elapsed_ms);
  }
  const double elapsed_mean_ms = elapsed_sum_ms / repetitions;

  std::cout << "n\tmode\twidth\trepetitions\telapsed_mean_ms"
               "\telapsed_min_ms\troot\twmbp_parent_visits"
               "\tvm_binary_join_probes\tvm_base_prefix_cell_visits"
               "\tvm_base_prefix_selected"
               "\twm_binary_join_probes\twm_extension_probes"
               "\twm_start_prefix_cell_visits"
               "\twm_start_prefix_selected"
               "\twi_binary_join_probes\twi_direct_probes"
               "\twi_extension_probes"
               "\twip_binary_join_probes\twip_extension_probes"
               "\twip_start_prefix_cell_visits"
               "\twip_start_prefix_selected"
               "\twmbw_binary_join_probes"
               "\tvpl_prefix_cell_visits\tvpl_prefix_selected"
               "\tvpr_binary_join_probes\tvpr_base_state_probes"
               "\tvpr_base_state_selected"
               "\tw_binary_join_probes\tw_direct_probes"
               "\tw_extension_probes\tw_root_decomposition_probes"
               "\twmb_split_positive_sources\twmb_direct_positive_sources"
               "\twmb_empty_positive_sources"
               "\twmb_btail_join_probes\twmb_btail_direct_probes"
               "\twmb_btail_selected\twmb_be_btail_join_probes"
               "\twmb_retained_direct_probes"
               "\twmb_retained_split_selected"
               "\twmbp_shared_wmbp_left_probes"
               "\twmbp_shared_wmbw_left_probes"
               "\twmbp_shared_emitted_sources"
               "\twmbp_wi_history_insertions\twmbp_wi_be_wi_probes"
               "\twmbp_wi_empty_probes\twmbp_wi_prefix_candidates"
               "\twmbp_wi_prefix_vp_probes\twmbp_wi_direct_probes"
               "\twmbp_wi_emitted_sources\tvp_loop_parent_probes"
               "\tvp_stack_event_sources\tvp_internal_event_sources"
               "\tcase3_active_expiry_probes"
               "\tcase3_tail_materialization_probes"
               "\tcase3_prepared_score_probes"
               "\tcase3_max_active_segments"
               "\tv_hairpin_parent_probes\tv_internal_parent_probes"
               "\tv_vm_parent_probes\tv_internal_sources"
               "\tcase3_heap_registered_segments"
               "\tcase3_heap_expired_segments\tcase3_max_heap_size\n"
            << sequence.size() << '\t' << mode << '\t' << width << '\t'
            << repetitions << '\t' << elapsed_mean_ms << '\t'
            << elapsed_min_ms << '\t' << root << '\t' << wmbp_parent_visits
            << '\t' << vm_binary_join_probes << '\t'
            << vm_base_prefix_cell_visits
            << '\t' << vm_base_prefix_selected
            << '\t' << wm_binary_join_probes << '\t' << wm_extension_probes
            << '\t' << wm_start_prefix_cell_visits
            << '\t' << wm_start_prefix_selected
            << '\t' << wi_binary_join_probes << '\t' << wi_direct_probes
            << '\t' << wi_extension_probes
            << '\t' << wip_binary_join_probes
            << '\t' << wip_extension_probes
            << '\t' << wip_start_prefix_cell_visits
            << '\t' << wip_start_prefix_selected
            << '\t' << wmbw_binary_join_probes
            << '\t' << vpl_prefix_cell_visits << '\t' << vpl_prefix_selected
            << '\t' << vpr_binary_join_probes
            << '\t' << vpr_base_state_probes
            << '\t' << vpr_base_state_selected
            << '\t' << w_binary_join_probes << '\t' << w_direct_probes
            << '\t' << w_extension_probes
            << '\t' << w_root_decomposition_probes
            << '\t' << wmb_split_positive_sources
            << '\t' << wmb_direct_positive_sources
            << '\t' << wmb_empty_positive_sources
            << '\t' << wmb_btail_join_probes
            << '\t' << wmb_btail_direct_probes
            << '\t' << wmb_btail_selected
            << '\t' << wmb_be_btail_join_probes
            << '\t' << wmb_retained_direct_probes
            << '\t' << wmb_retained_split_selected
            << '\t' << wmbp_shared_wmbp_left_probes
            << '\t' << wmbp_shared_wmbw_left_probes
            << '\t' << wmbp_shared_emitted_sources
            << '\t' << wmbp_wi_history_insertions
            << '\t' << wmbp_wi_be_wi_probes
            << '\t' << wmbp_wi_empty_probes
            << '\t' << wmbp_wi_prefix_candidates
            << '\t' << wmbp_wi_prefix_vp_probes
            << '\t' << wmbp_wi_direct_probes
            << '\t' << wmbp_wi_emitted_sources
            << '\t' << vp_loop_parent_probes
            << '\t' << vp_stack_event_sources
            << '\t' << vp_internal_event_sources
            << '\t' << case3_active_expiry_probes
            << '\t' << case3_tail_materialization_probes
            << '\t' << case3_prepared_score_probes
            << '\t' << case3_max_active_segments
            << '\t' << v_hairpin_parent_probes
            << '\t' << v_internal_parent_probes
            << '\t' << v_vm_parent_probes
            << '\t' << v_internal_sources
            << '\t' << case3_heap_registered_segments
            << '\t' << case3_heap_expired_segments
            << '\t' << case3_max_heap_size
            << "\n\nnonterminal\tvisits\tpositive\n";
  for (std::size_t idx = 0;
       idx < cparty::scfg2::replay::w_final_online_nonterminal_count;
       ++idx) {
    const auto nt = static_cast<cparty::scfg2::grammar::NonTerminal>(idx);
    std::cout << cparty::scfg2::grammar::nonterminal_name(nt) << '\t'
              << visits[idx] << '\t' << positive[idx] << '\n';
  }
  return EXIT_SUCCESS;
}

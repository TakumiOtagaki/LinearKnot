#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build_dir="${SCFG2_BUILD_DIR:-$root_dir/build}"

targets=(
  scfg2_grammar_smoke_test
  scfg2_engine_smoke_test
  scfg2_traceback_lifecycle_test
  scfg2_mfe_output_test
  scfg2_compact_choices_test
  scfg2_agenda_source_transfer_test
  scfg2_runtime_snapshot_test
  scfg2_inside_pair_max_test
  scfg2_streaming_consumers_test
  scfg2_logsumexp_tape_test
  scfg2_factorized_outside_test
  scfg2_regenerated_outside_test
  scfg2_be_core_left_retained_events_test
  scfg2_be_wmbp_regenerated_sources_test
  scfg2_vpr_regenerated_base_cells_test
  scfg2_wmb_regenerated_btail_cells_test
  scfg2_retained_outside_test
  scfg2_vp_recursive_split_retained_events_test
  scfg2_vp_loop_retained_events_test
  scfg2_vm_wmp_base_prefix_test
  scfg2_wm_start_prefix_test
  scfg2_wip_start_prefix_test
  scfg2_wip_hairpin_boundary_test
  scfg2_vpl_start_prefix_test
  scfg2_wmbp_prefix_recompute_test
  scfg2_beam_analytic_outside_oracle_test
  scfg2_exact_outside_derivative_test
  scfg2_exact_regenerated_outside_test
  scfg2_beam_outside_retained_dag_test
  scfg2_outside_api_test
  scfg2_outside_random_consistency_test
  scfg2_v_transition_smoke_test
  scfg2_v_retained_agenda_test
  scfg2_rules_core_exact_weight_test
  scfg2_fixed_rule_expansion_test
  scfg2_internal_loop_weight_test
  scfg2_vp_internal_border_regression_test
  scfg2_wmbp_be_availability_test
  scfg2_exact_inside_oracle_test
  scfg2_wmb_closure_policy_test
  scfg2_constrained_small_n_grouping_test
  scfg2_pk_profiles_oracle_test
  scfg2_pk_profiles_beam_test
  scfg2_pk_sampling_test
  scfg2_lazy_kbest_test
  scfg2_subopt_test
  scfg2_mea_test
  scfg2_maximal_pair_index_test
  scfg2_split_membership_test
  scfg2_border_query_test
  scfg2_padding_predicate_test
  scfg2_be_empty_guard_test
  scfg2_schedule_invariant_test
  scfg2_static_range_minimum_test
  scfg2_bounded_ranked_bucket_test
  scfg2_wmbp_wi_prefix_store_test
  scfg2_vp_case3_nonleft_last_child_index_test
  scfg2_vp_case12_last_child_index_test
  scfg2_vp_case3_endpoint_prefix_frontier_test
  scfg2_case3_setup_width_test
  scfg2_prefix_lifetime_test
  scfg2_right_boundary_frontier_test
  scfg2_beam_scaffold_forcing_test
  scfg2_pilot_guided_inside_test
  scfg2_outer_ranked_outside_test
  scfg2_outer_ranking_api_test
  scfg2_exact_api_test
  scfg2_api_input_boundary_test
  scfg2_log_weight_test
  scfg2_beam_api_test
  scfg2_beam_memory_observation_test
  scfg2_hotspot_wavefront_test
  scfg2_hotspot_viterbi_api_test
  scfg2_empty_scaffold_pkfree_contract_test
  scfg2_random_consistency_test
)

tool_targets=(
  linearknot
  scfg2_hotspot_proposal_tool
  scfg2_fixed_scaffold_viterbi_tool
  scfg2_hotspot_viterbi_tool
)

cmake --build "$build_dir" -j8 --target "${targets[@]}" "${tool_targets[@]}"

for target in "${targets[@]}"; do
  echo "[scfg2-required] start $target"
  "$build_dir/$target"
  echo "[scfg2-required] ok $target"
done

echo "[scfg2-required] start scfg2_pair_accuracy"
python3 "$root_dir/test/tools/test_scfg2_pair_accuracy.py"
echo "[scfg2-required] ok scfg2_pair_accuracy"

echo "[scfg2-required] start make_scaffold"
python3 "$root_dir/test/tools/test_make_scaffold.py"
echo "[scfg2-required] ok make_scaffold"

echo "[scfg2-required] start linearknot_mea_cli"
python3 "$root_dir/test/tools/test_linearknot_mea.py"
echo "[scfg2-required] ok linearknot_mea_cli"

echo "[scfg2-required] start linearknot_subopt_cli"
python3 "$root_dir/test/tools/test_linearknot_subopt.py"
echo "[scfg2-required] ok linearknot_subopt_cli"

echo "[scfg2-required] start linearknot_access_cli"
python3 "$root_dir/test/tools/test_linearknot_access.py"
echo "[scfg2-required] ok linearknot_access_cli"

echo "[scfg2-required] start outer_ranking_cli"
python3 "$root_dir/test/tools/test_outer_ranking_cli.py"
echo "[scfg2-required] ok outer_ranking_cli"

echo "[scfg2-required] start scfg2_installed_api"
python3 "$root_dir/test/tools/test_scfg2_installed_api.py" "$build_dir"
echo "[scfg2-required] ok scfg2_installed_api"

echo "[scfg2-required] start scfg2_exact_outside_memory"
python3 "$root_dir/test/tools/test_scfg2_exact_outside_memory.py"
echo "[scfg2-required] ok scfg2_exact_outside_memory"

echo "[scfg2-required] start scfg2_pseudobase_audit"
python3 "$root_dir/test/tools/test_audit_pseudobase_benchmark.py"
echo "[scfg2-required] ok scfg2_pseudobase_audit"

echo "[scfg2-required] start scfg2_pseudobase_viterbi_runner"
python3 "$root_dir/test/tools/test_run_scfg2_pseudobase_viterbi.py"
echo "[scfg2-required] ok scfg2_pseudobase_viterbi_runner"

echo "[scfg2-required] start ipknot_pseudobase_runner"
python3 "$root_dir/test/tools/test_run_ipknot_pseudobase.py"
echo "[scfg2-required] ok ipknot_pseudobase_runner"

echo "[scfg2-required] start pseudobase_benchmark_summary"
python3 "$root_dir/test/tools/test_summarize_pseudobase_benchmark.py"
echo "[scfg2-required] ok pseudobase_benchmark_summary"

"$root_dir/test/tools/check_scfg2_scfg_dependency_gate.sh"

echo "[scfg2-required] required set ok"

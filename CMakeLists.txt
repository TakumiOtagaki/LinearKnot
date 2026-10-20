cmake_minimum_required(VERSION 3.8)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" LINEARKNOT_VERSION LIMIT_COUNT 1)
project(LinearKnot VERSION ${LINEARKNOT_VERSION} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Bundled ViennaRNA uses pre-C23 function declarations.
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
option(LINEARKNOT_ENABLE_LTO "Enable LTO in Release builds" ON)
set(LINEARKNOT_LTO_JOBS "" CACHE STRING
    "GNU Release LTO workers per executable link; empty keeps compiler default")
if(NOT "${LINEARKNOT_LTO_JOBS}" STREQUAL "" AND
   NOT "${LINEARKNOT_LTO_JOBS}" MATCHES "^[1-9][0-9]*$")
  message(FATAL_ERROR "LINEARKNOT_LTO_JOBS must be empty or a positive integer")
endif()
if(LINEARKNOT_ENABLE_LTO)
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " -flto")
  if(NOT "${LINEARKNOT_LTO_JOBS}" STREQUAL "")
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      message(FATAL_ERROR "LINEARKNOT_LTO_JOBS requires the GNU C++ compiler when LTO is enabled")
    endif()
    # Keep object compilation unchanged; only parallelize executable linking.
    string(APPEND CMAKE_EXE_LINKER_FLAGS_RELEASE " -flto=${LINEARKNOT_LTO_JOBS}")
  endif()
endif()

# Python wheels put native commands in the installer's scripts directory.
# Ordinary CMake/Conda installs continue to use PREFIX/bin.
if(SKBUILD)
  set(LINEARKNOT_RUNTIME_DESTINATION "${SKBUILD_SCRIPTS_DIR}")
else()
  set(LINEARKNOT_RUNTIME_DESTINATION bin)
endif()

include(cmake/BundledViennaRNA.cmake)
linearknot_add_bundled_vienna()

set(SCFG2_EXACT_RUNTIME_SOURCES
  src/scfg2/api/exact.cc
  src/scfg2/weights/runtime_snapshot.cc
  src/scfg2/storage/exact_store.cc
  src/rna_pair_table.cc
)

add_library(cparty_scfg2_exact STATIC ${SCFG2_EXACT_RUNTIME_SOURCES})
add_library(LinearKnot::scfg2_exact ALIAS cparty_scfg2_exact)
target_include_directories(
  cparty_scfg2_exact
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
         $<INSTALL_INTERFACE:include>
)
target_compile_features(cparty_scfg2_exact PUBLIC cxx_std_17)
target_link_libraries(cparty_scfg2_exact PUBLIC RNA)

set_target_properties(cparty_scfg2_exact PROPERTIES EXPORT_NAME scfg2_exact)
install(TARGETS cparty_scfg2_exact EXPORT LinearKnotTargets ARCHIVE DESTINATION lib)
install(DIRECTORY src/scfg2/ DESTINATION include/scfg2 FILES_MATCHING PATTERN "*.hh")

add_library(cparty_scfg2_beam STATIC src/scfg2/api/beam.cc src/scfg2/api/profiles.cc src/scfg2/api/sampling.cc src/scfg2/api/mea.cc src/scfg2/api/subopt.cc)
add_library(LinearKnot::scfg2_beam ALIAS cparty_scfg2_beam)
target_include_directories(
  cparty_scfg2_beam
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
         $<INSTALL_INTERFACE:include>
)
target_compile_features(cparty_scfg2_beam PUBLIC cxx_std_17)
target_link_libraries(cparty_scfg2_beam PUBLIC LinearKnot::scfg2_exact)
set_target_properties(cparty_scfg2_beam PROPERTIES EXPORT_NAME scfg2_beam)
install(TARGETS cparty_scfg2_beam EXPORT LinearKnotTargets ARCHIVE DESTINATION lib)

add_executable(linearknot
  src/scfg2/cli/main.cc
  src/scfg2/cli/fold.cc
  src/scfg2/cli/profile.cc
  src/scfg2/cli/access.cc
  src/scfg2/cli/sample.cc
  src/scfg2/cli/subopt.cc
  src/scfg2/cli/decode.cc
)
target_link_libraries(linearknot PRIVATE LinearKnot::scfg2_beam)
target_compile_definitions(linearknot PRIVATE LINEARKNOT_VERSION="${PROJECT_VERSION}")
install(TARGETS linearknot RUNTIME DESTINATION "${LINEARKNOT_RUNTIME_DESTINATION}" COMPONENT cli)

add_library(
  cparty_scfg2_hotspot_proposal STATIC
  src/scfg2/proposal/hotspot_wavefront.cc
)
add_library(
  LinearKnot::scfg2_hotspot_proposal ALIAS cparty_scfg2_hotspot_proposal
)
target_include_directories(
  cparty_scfg2_hotspot_proposal
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
         $<INSTALL_INTERFACE:include>
)
target_compile_features(cparty_scfg2_hotspot_proposal PUBLIC cxx_std_17)
set_target_properties(cparty_scfg2_hotspot_proposal PROPERTIES EXPORT_NAME scfg2_hotspot_proposal)
install(TARGETS cparty_scfg2_hotspot_proposal EXPORT LinearKnotTargets ARCHIVE DESTINATION lib)

add_library(
  cparty_scfg2_hotspot_viterbi STATIC
  src/scfg2/api/hotspot_viterbi.cc
)
add_library(
  LinearKnot::scfg2_hotspot_viterbi ALIAS cparty_scfg2_hotspot_viterbi
)
target_include_directories(
  cparty_scfg2_hotspot_viterbi
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
         $<INSTALL_INTERFACE:include>
)
target_compile_features(cparty_scfg2_hotspot_viterbi PUBLIC cxx_std_17)
target_link_libraries(
  cparty_scfg2_hotspot_viterbi
  PUBLIC LinearKnot::scfg2_beam LinearKnot::scfg2_hotspot_proposal
)
set_target_properties(cparty_scfg2_hotspot_viterbi PROPERTIES EXPORT_NAME scfg2_hotspot_viterbi)
install(TARGETS cparty_scfg2_hotspot_viterbi EXPORT LinearKnotTargets ARCHIVE DESTINATION lib)
install(EXPORT LinearKnotTargets NAMESPACE LinearKnot:: DESTINATION lib/cmake/LinearKnot)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/LinearKnotConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY ExactVersion
)
install(FILES cmake/LinearKnotConfig.cmake
              "${CMAKE_CURRENT_BINARY_DIR}/LinearKnotConfigVersion.cmake"
        DESTINATION lib/cmake/LinearKnot)
install(FILES LICENSE VERSION THIRD_PARTY_NOTICES.md DESTINATION share/LinearKnot)
install(FILES src/ViennaRNA/COPYING src/ViennaRNA/AUTHORS
        DESTINATION share/LinearKnot/licenses/ViennaRNA)

add_custom_target(
  scfg2_scfg_dependency_gate
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test/tools/check_scfg2_scfg_dependency_gate.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)

include(CTest)
if(BUILD_TESTING)
  add_executable(scfg2_exact_numeric_execution_test test/scfg2_exact_numeric_execution_test.cc)
  target_link_libraries(scfg2_exact_numeric_execution_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_exact_numeric_execution COMMAND scfg2_exact_numeric_execution_test)
  set_tests_properties(scfg2_exact_numeric_execution PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(scfg2_exact_live_support_test test/scfg2_exact_live_support_test.cc)
  target_link_libraries(scfg2_exact_live_support_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_exact_live_support COMMAND scfg2_exact_live_support_test)
  set_tests_properties(scfg2_exact_live_support PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(scfg2_exact_api_test test/scfg2_exact_api_test.cc)
  target_link_libraries(scfg2_exact_api_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_api_input_boundary_test test/scfg2_api_input_boundary_test.cc)
  target_link_libraries(scfg2_api_input_boundary_test PRIVATE LinearKnot::scfg2_beam)
  target_compile_definitions(scfg2_api_input_boundary_test PRIVATE
    SCFG2_TEST_PARAM_FILE="${CMAKE_CURRENT_SOURCE_DIR}/params/rna_Turner04.par")
  add_test(NAME scfg2_api_input_boundary COMMAND scfg2_api_input_boundary_test)

  add_executable(scfg2_log_weight_test test/scfg2_log_weight_test.cc)
  target_link_libraries(scfg2_log_weight_test PRIVATE LinearKnot::scfg2_beam)
  add_test(NAME scfg2_log_weight COMMAND scfg2_log_weight_test)
  set_tests_properties(scfg2_log_weight PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(scfg2_beam_api_test test/scfg2_beam_api_test.cc)
  target_link_libraries(scfg2_beam_api_test PRIVATE LinearKnot::scfg2_beam)

  add_executable(scfg2_beam_memory_tool test/scfg2_beam_memory_tool.cc)
  target_link_libraries(scfg2_beam_memory_tool PRIVATE LinearKnot::scfg2_beam)

  add_executable(scfg2_beam_memory_observation_test test/scfg2_beam_memory_observation_test.cc)
  target_link_libraries(scfg2_beam_memory_observation_test PRIVATE LinearKnot::scfg2_beam)

  add_executable(scfg2_mfe_output_test test/scfg2_mfe_output_test.cc)
  add_executable(scfg2_compact_choices_test test/scfg2_compact_choices_test.cc)
  target_link_libraries(scfg2_compact_choices_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_compact_choices COMMAND scfg2_compact_choices_test)
  add_executable(scfg2_agenda_source_transfer_test test/scfg2_agenda_source_transfer_test.cc)
  target_link_libraries(scfg2_agenda_source_transfer_test PRIVATE LinearKnot::scfg2_exact)
  target_compile_definitions(scfg2_agenda_source_transfer_test PRIVATE
    SCFG2_TEST_PARAM_FILE="${CMAKE_CURRENT_SOURCE_DIR}/params/rna_Turner04.par")
  add_test(NAME scfg2_agenda_source_transfer COMMAND scfg2_agenda_source_transfer_test)
  add_executable(scfg2_runtime_snapshot_test test/scfg2_runtime_snapshot_test.cc)
  target_link_libraries(scfg2_runtime_snapshot_test PRIVATE LinearKnot::scfg2_exact)
  target_compile_definitions(scfg2_runtime_snapshot_test PRIVATE
    SCFG2_TEST_PARAM_FILE="${CMAKE_CURRENT_SOURCE_DIR}/params/rna_Turner04.par")
  add_test(NAME scfg2_runtime_snapshot COMMAND scfg2_runtime_snapshot_test)
  target_link_libraries(scfg2_mfe_output_test PRIVATE LinearKnot::scfg2_beam)
  target_compile_definitions(scfg2_mfe_output_test PRIVATE
    SCFG2_TEST_PARAM_FILE="${CMAKE_CURRENT_SOURCE_DIR}/params/rna_Turner04.par")
  add_test(NAME scfg2_mfe_output COMMAND scfg2_mfe_output_test)

  add_executable(scfg2_inside_pair_max_test test/scfg2_inside_pair_max_test.cc)
  target_link_libraries(scfg2_inside_pair_max_test PRIVATE LinearKnot::scfg2_beam)
  target_compile_definitions(scfg2_inside_pair_max_test PRIVATE
    SCFG2_TEST_PARAM_FILE="${CMAKE_CURRENT_SOURCE_DIR}/params/rna_Turner04.par")
  add_test(NAME scfg2_inside_pair_max COMMAND scfg2_inside_pair_max_test)

  add_executable(scfg2_streaming_consumers_test test/scfg2_streaming_consumers_test.cc)
  target_link_libraries(scfg2_streaming_consumers_test PRIVATE LinearKnot::scfg2_beam)
  target_compile_definitions(scfg2_streaming_consumers_test PRIVATE
    SCFG2_TEST_PARAM_FILE="${CMAKE_CURRENT_SOURCE_DIR}/params/rna_Turner04.par")
  add_test(NAME scfg2_streaming_consumers COMMAND scfg2_streaming_consumers_test)

  add_executable(
    scfg2_hotspot_wavefront_test
    test/scfg2_hotspot_wavefront_test.cc
  )
  target_link_libraries(
    scfg2_hotspot_wavefront_test
    PRIVATE LinearKnot::scfg2_hotspot_proposal
  )

  add_executable(
    scfg2_hotspot_viterbi_api_test
    test/scfg2_hotspot_viterbi_api_test.cc
  )
  target_link_libraries(
    scfg2_hotspot_viterbi_api_test
    PRIVATE LinearKnot::scfg2_hotspot_viterbi
  )

  add_executable(
    scfg2_empty_scaffold_pkfree_contract_test
    test/scfg2_empty_scaffold_pkfree_contract_test.cc
  )
  target_link_libraries(
    scfg2_empty_scaffold_pkfree_contract_test
    PRIVATE LinearKnot::scfg2_beam
  )

  add_executable(
    scfg2_hotspot_viterbi_tool
    test/scfg2_hotspot_viterbi_tool.cc
  )
  target_link_libraries(
    scfg2_hotspot_viterbi_tool
    PRIVATE LinearKnot::scfg2_hotspot_viterbi
  )

  add_executable(
    scfg2_hotspot_proposal_tool
    test/scfg2_hotspot_proposal_tool.cc
  )
  target_link_libraries(
    scfg2_hotspot_proposal_tool
    PRIVATE LinearKnot::scfg2_hotspot_proposal
  )

  add_executable(
    scfg2_fixed_scaffold_viterbi_tool
    test/scfg2_fixed_scaffold_viterbi_tool.cc
  )
  target_link_libraries(
    scfg2_fixed_scaffold_viterbi_tool
    PRIVATE LinearKnot::scfg2_beam
  )

  add_executable(
    scfg2_beam_inside_saturation_tool
    test/scfg2_beam_inside_saturation_tool.cc
  )
  target_link_libraries(
    scfg2_beam_inside_saturation_tool
    PRIVATE LinearKnot::scfg2_beam
  )

  add_executable(
    scfg2_beam_outside_saturation_tool
    test/scfg2_beam_outside_saturation_tool.cc
  )
  target_link_libraries(
    scfg2_beam_outside_saturation_tool
    PRIVATE LinearKnot::scfg2_beam
  )

  add_executable(
    scfg2_random_consistency_test
    test/scfg2_random_consistency_test.cc
  )
  target_link_libraries(scfg2_random_consistency_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(
    scfg2_exact_inside_oracle_test
    test/scfg2_exact_inside_oracle_test.cc
  )
  target_link_libraries(scfg2_exact_inside_oracle_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_wip_hairpin_boundary_test test/scfg2_wip_hairpin_boundary_test.cc)
  target_link_libraries(scfg2_wip_hairpin_boundary_test PRIVATE LinearKnot::scfg2_beam)
  add_test(NAME scfg2_wip_hairpin_boundary COMMAND scfg2_wip_hairpin_boundary_test)
  set_tests_properties(scfg2_wip_hairpin_boundary PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_wmb_closure_policy_test test/scfg2_wmb_closure_policy_test.cc)
  target_link_libraries(scfg2_wmb_closure_policy_test PRIVATE LinearKnot::scfg2_beam)

  add_executable(scfg2_internal_loop_weight_test test/scfg2_internal_loop_weight_test.cc)
  target_link_libraries(scfg2_internal_loop_weight_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(
    scfg2_constrained_small_n_grouping_test
    test/scfg2_constrained_small_n_grouping_test.cc
  )
  target_link_libraries(
    scfg2_constrained_small_n_grouping_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(scfg2_tape_reachability_tool test/scfg2_tape_reachability_tool.cc)
  target_link_libraries(scfg2_tape_reachability_tool PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_beam_root_trace_tool test/scfg2_beam_root_trace_tool.cc)
  target_link_libraries(scfg2_beam_root_trace_tool PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_beam_scaffold_forcing_test test/scfg2_beam_scaffold_forcing_test.cc)
  target_link_libraries(scfg2_beam_scaffold_forcing_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_pilot_guided_inside_test test/scfg2_pilot_guided_inside_test.cc)
  target_link_libraries(scfg2_pilot_guided_inside_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_pilot_guided_inside COMMAND scfg2_pilot_guided_inside_test)
  set_tests_properties(scfg2_pilot_guided_inside PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_bounded_outer_alpha_lookup_test test/scfg2_bounded_outer_alpha_lookup_test.cc)
  target_link_libraries(scfg2_bounded_outer_alpha_lookup_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_bounded_outer_alpha_lookup COMMAND scfg2_bounded_outer_alpha_lookup_test)
  set_tests_properties(scfg2_bounded_outer_alpha_lookup PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_outer_max_pool_test test/scfg2_outer_max_pool_test.cc)
  target_link_libraries(scfg2_outer_max_pool_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_outer_max_pool COMMAND scfg2_outer_max_pool_test)
  set_tests_properties(scfg2_outer_max_pool PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_primary_lookup_index_test test/scfg2_primary_lookup_index_test.cc)
  target_link_libraries(scfg2_primary_lookup_index_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_primary_lookup_index COMMAND scfg2_primary_lookup_index_test)
  set_tests_properties(scfg2_primary_lookup_index PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_numeric_primary_lookup_test test/scfg2_numeric_primary_lookup_test.cc)
  target_link_libraries(scfg2_numeric_primary_lookup_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_numeric_primary_lookup COMMAND scfg2_numeric_primary_lookup_test)
  set_tests_properties(scfg2_numeric_primary_lookup PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_beam_numeric_lookup_test test/scfg2_beam_numeric_lookup_test.cc)
  target_link_libraries(scfg2_beam_numeric_lookup_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_beam_numeric_lookup COMMAND scfg2_beam_numeric_lookup_test)
  set_tests_properties(scfg2_beam_numeric_lookup PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_outer_ranked_outside_test test/scfg2_outer_ranked_outside_test.cc)
  target_link_libraries(scfg2_outer_ranked_outside_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_outer_ranked_outside COMMAND scfg2_outer_ranked_outside_test)
  set_tests_properties(scfg2_outer_ranked_outside PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_outer_ranking_api_test test/scfg2_outer_ranking_api_test.cc)
  target_link_libraries(scfg2_outer_ranking_api_test PRIVATE LinearKnot::scfg2_beam)
  add_test(NAME scfg2_outer_ranking_api COMMAND scfg2_outer_ranking_api_test)
  set_tests_properties(scfg2_outer_ranking_api PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_pilot_guided_inside_probe test/tools/scfg2_pilot_guided_inside_probe.cc)
  target_link_libraries(scfg2_pilot_guided_inside_probe PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_outer_cost_probe test/tools/scfg2_outer_cost_probe.cc)
  target_link_libraries(scfg2_outer_cost_probe PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_main_inside_cost_probe test/tools/scfg2_main_inside_cost_probe.cc)
  target_link_libraries(scfg2_main_inside_cost_probe PRIVATE LinearKnot::scfg2_exact)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_executable(scfg2_main_inside_profile_probe EXCLUDE_FROM_ALL test/tools/scfg2_main_inside_cost_probe.cc)
    target_compile_definitions(scfg2_main_inside_profile_probe PRIVATE SCFG2_COST_GPROF)
    target_compile_options(scfg2_main_inside_profile_probe PRIVATE -pg)
    target_link_libraries(scfg2_main_inside_profile_probe PRIVATE LinearKnot::scfg2_exact -pg)
  endif()

  add_executable(scfg2_border_query_test test/scfg2_border_query_test.cc)
  target_link_libraries(scfg2_border_query_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_padding_predicate_test test/scfg2_padding_predicate_test.cc)
  target_link_libraries(scfg2_padding_predicate_test PRIVATE LinearKnot::scfg2_beam)
  add_test(NAME scfg2_padding_predicate COMMAND scfg2_padding_predicate_test)
  set_tests_properties(scfg2_padding_predicate PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_be_empty_guard_test test/scfg2_be_empty_guard_test.cc)
  target_link_libraries(scfg2_be_empty_guard_test PRIVATE LinearKnot::scfg2_beam)
  add_test(NAME scfg2_be_empty_guard COMMAND scfg2_be_empty_guard_test)
  set_tests_properties(scfg2_be_empty_guard PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_schedule_invariant_test test/scfg2_schedule_invariant_test.cc)
  target_link_libraries(scfg2_schedule_invariant_test PRIVATE LinearKnot::scfg2_exact)
  # Dependency bugs must remain visible when assertions are disabled.
  target_compile_definitions(scfg2_schedule_invariant_test PRIVATE NDEBUG)
  add_test(NAME scfg2_schedule_invariant COMMAND scfg2_schedule_invariant_test)
  set_tests_properties(scfg2_schedule_invariant PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(scfg2_bounded_ranked_bucket_test test/scfg2_bounded_ranked_bucket_test.cc)
  target_link_libraries(scfg2_bounded_ranked_bucket_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_wmbp_wi_prefix_store_test test/scfg2_wmbp_wi_prefix_store_test.cc)
  target_link_libraries(scfg2_wmbp_wi_prefix_store_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_be_wmbp_regenerated_sources_test test/scfg2_be_wmbp_regenerated_sources_test.cc)
  target_link_libraries(scfg2_be_wmbp_regenerated_sources_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_vpr_regenerated_base_cells_test test/scfg2_vpr_regenerated_base_cells_test.cc)
  target_link_libraries(scfg2_vpr_regenerated_base_cells_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_wmb_regenerated_btail_cells_test test/scfg2_wmb_regenerated_btail_cells_test.cc)
  target_link_libraries(scfg2_wmb_regenerated_btail_cells_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_retained_outside_test test/scfg2_retained_outside_test.cc)
  target_link_libraries(scfg2_retained_outside_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_static_range_minimum_test test/scfg2_static_range_minimum_test.cc)
  target_link_libraries(scfg2_static_range_minimum_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_split_membership_test test/scfg2_split_membership_test.cc)
  target_link_libraries(scfg2_split_membership_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_maximal_pair_index_test test/scfg2_maximal_pair_index_test.cc)
  target_link_libraries(scfg2_maximal_pair_index_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(scfg2_pk_sampling_test test/scfg2_pk_sampling_test.cc)
  target_link_libraries(scfg2_pk_sampling_test PRIVATE LinearKnot::scfg2_beam)

  add_executable(scfg2_lazy_kbest_test test/scfg2_lazy_kbest_test.cc)
  target_link_libraries(scfg2_lazy_kbest_test PRIVATE LinearKnot::scfg2_beam)
  add_test(NAME scfg2_lazy_kbest COMMAND scfg2_lazy_kbest_test)

  add_executable(scfg2_subopt_test test/scfg2_subopt_test.cc)
  target_link_libraries(scfg2_subopt_test PRIVATE LinearKnot::scfg2_beam)
  add_test(NAME scfg2_subopt COMMAND scfg2_subopt_test)
  set_tests_properties(scfg2_subopt PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

  add_executable(
    scfg2_pk_profiles_beam_test
    test/scfg2_pk_profiles_beam_test.cc
  )
  target_link_libraries(scfg2_pk_profiles_beam_test PRIVATE LinearKnot::scfg2_beam)

  add_executable(
    scfg2_pk_profiles_oracle_test
    test/scfg2_pk_profiles_oracle_test.cc
  )
  target_link_libraries(scfg2_pk_profiles_oracle_test PRIVATE LinearKnot::scfg2_exact)

  add_executable(
    scfg2_rules_core_exact_weight_test
    test/scfg2_rules_core_exact_weight_test.cc
  )
  target_include_directories(scfg2_rules_core_exact_weight_test PRIVATE src)

  add_executable(
    scfg2_fixed_rule_expansion_test
    test/scfg2_fixed_rule_expansion_test.cc
  )
  target_include_directories(scfg2_fixed_rule_expansion_test PRIVATE src)

  add_executable(
    scfg2_vp_internal_border_regression_test
    test/scfg2_vp_internal_border_regression_test.cc
  )
  target_link_libraries(
    scfg2_vp_internal_border_regression_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(scfg2_grammar_smoke_test test/scfg2_grammar_smoke_test.cc)
  add_executable(scfg2_engine_smoke_test test/scfg2_engine_smoke_test.cc)
  add_executable(scfg2_traceback_lifecycle_test test/scfg2_traceback_lifecycle_test.cc)
  target_include_directories(scfg2_traceback_lifecycle_test PRIVATE src)
  add_executable(
    scfg2_logsumexp_tape_test
    test/scfg2_logsumexp_tape_test.cc
  )
  add_executable(scfg2_factorized_outside_test test/scfg2_factorized_outside_test.cc)
  target_include_directories(scfg2_factorized_outside_test PRIVATE src)
  add_executable(scfg2_regenerated_outside_test test/scfg2_regenerated_outside_test.cc)
  target_link_libraries(scfg2_regenerated_outside_test PRIVATE LinearKnot::scfg2_exact)
  add_executable(scfg2_wmbp_prefix_recompute_test test/scfg2_wmbp_prefix_recompute_test.cc)
  target_include_directories(scfg2_wmbp_prefix_recompute_test PRIVATE src)
  add_executable(scfg2_mea_test test/scfg2_mea_test.cc)
  target_link_libraries(scfg2_mea_test PRIVATE LinearKnot::scfg2_beam)
  add_executable(
    scfg2_exact_outside_derivative_test
    test/scfg2_exact_outside_derivative_test.cc
  )
  target_link_libraries(
    scfg2_exact_outside_derivative_test
    PRIVATE LinearKnot::scfg2_exact
  )
  add_executable(scfg2_exact_regenerated_outside_test test/scfg2_exact_regenerated_outside_test.cc)
  target_link_libraries(scfg2_exact_regenerated_outside_test PRIVATE LinearKnot::scfg2_beam)
  add_executable(
    scfg2_beam_outside_retained_dag_test
    test/scfg2_beam_outside_retained_dag_test.cc
  )
  target_link_libraries(
    scfg2_beam_outside_retained_dag_test
    PRIVATE LinearKnot::scfg2_beam
  )
  add_executable(scfg2_beam_analytic_outside_oracle_test test/scfg2_beam_analytic_outside_oracle_test.cc)
  target_link_libraries(scfg2_beam_analytic_outside_oracle_test PRIVATE LinearKnot::scfg2_beam)
  add_executable(scfg2_analytic_outside_memory_tool test/scfg2_analytic_outside_memory_tool.cc)
  target_link_libraries(scfg2_analytic_outside_memory_tool PRIVATE LinearKnot::scfg2_beam)
  add_executable(scfg2_outside_pair_audit_tool test/scfg2_outside_pair_audit_tool.cc)
  target_link_libraries(scfg2_outside_pair_audit_tool PRIVATE LinearKnot::scfg2_beam)
  add_executable(scfg2_wmbp_be_availability_test test/scfg2_wmbp_be_availability_test.cc)
  target_link_libraries(scfg2_wmbp_be_availability_test PRIVATE LinearKnot::scfg2_beam)
  add_executable(
    scfg2_outside_api_test
    test/scfg2_outside_api_test.cc
  )
  target_link_libraries(
    scfg2_outside_api_test
    PRIVATE LinearKnot::scfg2_beam
  )
  add_executable(
    scfg2_outside_random_consistency_test
    test/scfg2_outside_random_consistency_test.cc
  )
  target_link_libraries(
    scfg2_outside_random_consistency_test
    PRIVATE LinearKnot::scfg2_beam
  )
  add_executable(scfg2_v_transition_smoke_test test/scfg2_v_transition_smoke_test.cc)
  foreach(target_name
      scfg2_grammar_smoke_test
      scfg2_engine_smoke_test
      scfg2_logsumexp_tape_test
      scfg2_v_transition_smoke_test)
    target_include_directories(${target_name} PRIVATE src)
  endforeach()

  add_executable(
    scfg2_virtual_beam_candidate_test
    test/scfg2_virtual_beam_candidate_test.cc
  )
  target_link_libraries(
    scfg2_virtual_beam_candidate_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_virtual_beam_profile_tool
    test/scfg2_virtual_beam_profile_tool.cc
  )
  target_link_libraries(
    scfg2_virtual_beam_profile_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case3_tail_cache_test
    test/scfg2_vp_case3_tail_cache_test.cc
  )
  target_link_libraries(
    scfg2_vp_case3_tail_cache_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case3_endpoint_prefix_frontier_test
    test/scfg2_vp_case3_endpoint_prefix_frontier_test.cc
  )
  target_link_libraries(
    scfg2_vp_case3_endpoint_prefix_frontier_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case3_left_cache_test
    test/scfg2_vp_case3_left_cache_test.cc
  )
  target_link_libraries(
    scfg2_vp_case3_left_cache_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case3_left_top_width_test
    test/scfg2_vp_case3_left_top_width_test.cc
  )
  target_link_libraries(
    scfg2_vp_case3_left_top_width_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case3_tail_top_width_test
    test/scfg2_vp_case3_tail_top_width_test.cc
  )
  target_link_libraries(
    scfg2_vp_case3_tail_top_width_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_be_frame_reuse_profile_test
    test/scfg2_be_frame_reuse_profile_test.cc
  )
  target_link_libraries(
    scfg2_be_frame_reuse_profile_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmbp_be_vp_prefix_profile_test
    test/scfg2_wmbp_be_vp_prefix_profile_test.cc
  )
  target_link_libraries(
    scfg2_wmbp_be_vp_prefix_profile_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmbp_be_vp_prefix_top_width_test
    test/scfg2_wmbp_be_vp_prefix_top_width_test.cc
  )
  target_link_libraries(
    scfg2_wmbp_be_vp_prefix_top_width_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmbp_be_vp_split_index_profile_tool
    test/scfg2_wmbp_be_vp_split_index_profile_tool.cc
  )
  target_link_libraries(
    scfg2_wmbp_be_vp_split_index_profile_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmbp_vp_phase_barrier_test
    test/scfg2_wmbp_vp_phase_barrier_test.cc
  )
  target_link_libraries(
    scfg2_wmbp_vp_phase_barrier_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmbp_vp_child_beam_test
    test/scfg2_wmbp_vp_child_beam_test.cc
  )
  target_link_libraries(
    scfg2_wmbp_vp_child_beam_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmbp_item_beam_test
    test/scfg2_wmbp_item_beam_test.cc
  )
  target_link_libraries(
    scfg2_wmbp_item_beam_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_wmbp_composed_beam_test
    test/scfg2_vp_wmbp_composed_beam_test.cc
  )
  target_link_libraries(
    scfg2_vp_wmbp_composed_beam_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_online_beam_work_profile_tool
    test/scfg2_online_beam_work_profile_tool.cc
  )
  target_link_libraries(
    scfg2_online_beam_work_profile_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_logmax_width_profile_tool
    test/scfg2_logmax_width_profile_tool.cc
  )
  target_link_libraries(
    scfg2_logmax_width_profile_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_coupled_log_probability_regression_test
    test/scfg2_coupled_log_probability_regression_test.cc
  )
  target_link_libraries(
    scfg2_coupled_log_probability_regression_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_finite_setup_profile_tool
    test/scfg2_finite_setup_profile_tool.cc
  )
  target_link_libraries(
    scfg2_finite_setup_profile_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_approximation_quality_tool
    test/scfg2_approximation_quality_tool.cc
  )
  target_link_libraries(
    scfg2_approximation_quality_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_viterbi_width_divergence_tool
    test/scfg2_viterbi_width_divergence_tool.cc
  )
  target_link_libraries(
    scfg2_viterbi_width_divergence_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_semiring_dag_divergence_tool
    test/scfg2_semiring_dag_divergence_tool.cc
  )
  target_link_libraries(
    scfg2_semiring_dag_divergence_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmb_consumer_view_monotonicity_test
    test/scfg2_wmb_consumer_view_monotonicity_test.cc
  )
  target_link_libraries(
    scfg2_wmb_consumer_view_monotonicity_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_source_inventory_tool
    test/scfg2_vp_source_inventory_tool.cc
  )
  target_link_libraries(
    scfg2_vp_source_inventory_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_loop_retained_events_test
    test/scfg2_vp_loop_retained_events_test.cc
  )
  target_link_libraries(
    scfg2_vp_loop_retained_events_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case3_nonleft_last_child_index_test
    test/scfg2_vp_case3_nonleft_last_child_index_test.cc
  )
  target_link_libraries(
    scfg2_vp_case3_nonleft_last_child_index_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case12_last_child_index_test
    test/scfg2_vp_case12_last_child_index_test.cc
  )
  target_link_libraries(
    scfg2_vp_case12_last_child_index_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(scfg2_prefix_lifetime_test test/scfg2_prefix_lifetime_test.cc)
  target_link_libraries(scfg2_prefix_lifetime_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_prefix_lifetime COMMAND scfg2_prefix_lifetime_test)
  set_tests_properties(scfg2_prefix_lifetime PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(scfg2_right_boundary_frontier_test test/scfg2_right_boundary_frontier_test.cc)
  target_link_libraries(scfg2_right_boundary_frontier_test PRIVATE LinearKnot::scfg2_exact)
  add_test(NAME scfg2_right_boundary_frontier COMMAND scfg2_right_boundary_frontier_test)

  add_executable(
    scfg2_vp_recursive_split_retained_events_test
    test/scfg2_vp_recursive_split_retained_events_test.cc
  )
  target_link_libraries(
    scfg2_vp_recursive_split_retained_events_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_be_base_stack_retained_events_test
    test/scfg2_be_base_stack_retained_events_test.cc
  )
  target_link_libraries(
    scfg2_be_base_stack_retained_events_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_be_core_left_retained_events_test
    test/scfg2_be_core_left_retained_events_test.cc
  )
  target_link_libraries(
    scfg2_be_core_left_retained_events_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_be_right_tail_retained_core_index_test
    test/scfg2_be_right_tail_retained_core_index_test.cc
  )
  target_link_libraries(
    scfg2_be_right_tail_retained_core_index_test
    PRIVATE LinearKnot::scfg2_exact
  )


  add_executable(
    scfg2_v_retained_agenda_test
    test/scfg2_v_retained_agenda_test.cc
  )
  target_link_libraries(
    scfg2_v_retained_agenda_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_case3_setup_width_test
    test/scfg2_case3_setup_width_test.cc
  )
  target_link_libraries(
    scfg2_case3_setup_width_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmbp_shared_retained_join_test
    test/scfg2_wmbp_shared_retained_join_test.cc
  )
  target_link_libraries(
    scfg2_wmbp_shared_retained_join_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmbp_wi_timing_test
    test/scfg2_wmbp_wi_timing_test.cc
  )
  target_link_libraries(
    scfg2_wmbp_wi_timing_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_be_right_tail_top_width_test
    test/scfg2_be_right_tail_top_width_test.cc
  )
  target_link_libraries(
    scfg2_be_right_tail_top_width_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_be_right_tail_index_profile_tool
    test/scfg2_be_right_tail_index_profile_tool.cc
  )
  target_link_libraries(
    scfg2_be_right_tail_index_profile_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_be_left_right_frame_top_width_test
    test/scfg2_be_left_right_frame_top_width_test.cc
  )
  target_link_libraries(
    scfg2_be_left_right_frame_top_width_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case3_tail_schedule_tool
    test/scfg2_vp_case3_tail_schedule_tool.cc
  )
  target_link_libraries(
    scfg2_vp_case3_tail_schedule_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case3_tail_border_profile_tool
    test/scfg2_vp_case3_tail_border_profile_tool.cc
  )
  target_link_libraries(
    scfg2_vp_case3_tail_border_profile_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vp_case3_online_profile_tool
    test/scfg2_vp_case3_online_profile_tool.cc
  )
  target_link_libraries(
    scfg2_vp_case3_online_profile_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_right_boundary_dependency_tool
    test/scfg2_right_boundary_dependency_tool.cc
  )
  target_link_libraries(
    scfg2_right_boundary_dependency_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_right_boundary_dependency_test
    test/scfg2_right_boundary_dependency_test.cc
  )
  target_link_libraries(
    scfg2_right_boundary_dependency_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_right_boundary_schedule_test
    test/scfg2_right_boundary_schedule_test.cc
  )
  target_link_libraries(
    scfg2_right_boundary_schedule_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_right_boundary_item_frontier_test
    test/scfg2_right_boundary_item_frontier_test.cc
  )
  target_link_libraries(
    scfg2_right_boundary_item_frontier_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vm_wmp_base_prefix_test
    test/scfg2_vm_wmp_base_prefix_test.cc
  )
  target_link_libraries(
    scfg2_vm_wmp_base_prefix_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wm_start_prefix_test
    test/scfg2_wm_start_prefix_test.cc
  )
  target_link_libraries(
    scfg2_wm_start_prefix_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wip_start_prefix_test
    test/scfg2_wip_start_prefix_test.cc
  )
  target_link_libraries(
    scfg2_wip_start_prefix_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vpl_start_prefix_test
    test/scfg2_vpl_start_prefix_test.cc
  )
  target_link_libraries(
    scfg2_vpl_start_prefix_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_vpr_base_prefix_test
    test/scfg2_vpr_base_prefix_test.cc
  )
  target_link_libraries(
    scfg2_vpr_base_prefix_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_wmb_btail_factor_test
    test/scfg2_wmb_btail_factor_test.cc
  )
  target_link_libraries(
    scfg2_wmb_btail_factor_test
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(
    scfg2_restricted_partition_function_tool
    test/scfg2_restricted_partition_function_tool.cc
  )
  target_link_libraries(
    scfg2_restricted_partition_function_tool
    PRIVATE LinearKnot::scfg2_exact
  )

  add_executable(scfg2_exact_trace_tool test/scfg2_exact_trace_tool.cc)
  target_link_libraries(scfg2_exact_trace_tool PRIVATE LinearKnot::scfg2_exact)

  add_test(NAME scfg2_exact_api_test COMMAND scfg2_exact_api_test)
  add_test(NAME scfg2_beam_api_test COMMAND scfg2_beam_api_test)
  add_test(NAME scfg2_beam_memory_observation_test COMMAND scfg2_beam_memory_observation_test)
  add_test(NAME scfg2_factorized_outside COMMAND scfg2_factorized_outside_test)
  add_test(NAME scfg2_regenerated_outside COMMAND scfg2_regenerated_outside_test)
  add_test(NAME scfg2_wmbp_prefix_recompute COMMAND scfg2_wmbp_prefix_recompute_test)
  set_tests_properties(scfg2_regenerated_outside PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME scfg2_mea COMMAND scfg2_mea_test)
  set_tests_properties(scfg2_mea PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME scfg2_beam_analytic_outside_oracle COMMAND scfg2_beam_analytic_outside_oracle_test)
  add_test(
    NAME scfg2_hotspot_wavefront_test
    COMMAND scfg2_hotspot_wavefront_test
  )
  add_test(
    NAME scfg2_hotspot_viterbi_api_test
    COMMAND scfg2_hotspot_viterbi_api_test
  )
  add_test(
    NAME scfg2_empty_scaffold_pkfree_contract_test
    COMMAND scfg2_empty_scaffold_pkfree_contract_test
  )
  add_test(NAME scfg2_random_consistency COMMAND scfg2_random_consistency_test)
  add_test(NAME scfg2_exact_inside_oracle COMMAND scfg2_exact_inside_oracle_test)
  add_test(NAME scfg2_wmb_closure_policy COMMAND scfg2_wmb_closure_policy_test)
  set_tests_properties(scfg2_wmb_closure_policy PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME scfg2_internal_loop_weight COMMAND scfg2_internal_loop_weight_test)
  add_test(
    NAME scfg2_constrained_small_n_grouping
    COMMAND scfg2_constrained_small_n_grouping_test
  )
  add_test(NAME scfg2_rules_core_exact_weight COMMAND scfg2_rules_core_exact_weight_test)
  add_test(NAME scfg2_fixed_rule_expansion COMMAND scfg2_fixed_rule_expansion_test)
  add_test(NAME scfg2_pk_profiles_oracle COMMAND scfg2_pk_profiles_oracle_test)
  add_test(NAME scfg2_exact_regenerated_outside COMMAND scfg2_exact_regenerated_outside_test)
  set_tests_properties(scfg2_exact_regenerated_outside PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME scfg2_pk_profiles_beam COMMAND scfg2_pk_profiles_beam_test)
  add_test(NAME scfg2_pk_sampling COMMAND scfg2_pk_sampling_test)
  add_test(NAME scfg2_maximal_pair_index COMMAND scfg2_maximal_pair_index_test)
  add_test(NAME scfg2_split_membership COMMAND scfg2_split_membership_test)
  add_test(NAME scfg2_border_query COMMAND scfg2_border_query_test)
  add_test(NAME scfg2_bounded_ranked_bucket COMMAND scfg2_bounded_ranked_bucket_test)
  add_test(NAME scfg2_wmbp_wi_prefix_store COMMAND scfg2_wmbp_wi_prefix_store_test)
  add_test(NAME scfg2_be_wmbp_regenerated_sources COMMAND scfg2_be_wmbp_regenerated_sources_test)
  add_test(NAME scfg2_vpr_regenerated_base_cells COMMAND scfg2_vpr_regenerated_base_cells_test)
  add_test(NAME scfg2_wmb_regenerated_btail_cells COMMAND scfg2_wmb_regenerated_btail_cells_test)
  add_test(NAME scfg2_retained_outside COMMAND scfg2_retained_outside_test)
  add_test(NAME scfg2_static_range_minimum COMMAND scfg2_static_range_minimum_test)
  add_test(NAME scfg2_beam_scaffold_forcing COMMAND scfg2_beam_scaffold_forcing_test)
  add_test(
    NAME scfg2_vp_internal_border_regression
    COMMAND scfg2_vp_internal_border_regression_test
  )
  add_test(NAME scfg2_grammar_smoke COMMAND scfg2_grammar_smoke_test)
  add_test(NAME scfg2_engine_smoke COMMAND scfg2_engine_smoke_test)
  add_test(NAME scfg2_traceback_lifecycle COMMAND scfg2_traceback_lifecycle_test)
  add_test(NAME scfg2_v_transition_smoke COMMAND scfg2_v_transition_smoke_test)
  add_test(
    NAME scfg2_virtual_beam_candidate
    COMMAND scfg2_virtual_beam_candidate_test
  )
  add_test(
    NAME scfg2_vp_case3_tail_cache
    COMMAND scfg2_vp_case3_tail_cache_test
  )
  add_test(
    NAME scfg2_vp_case3_endpoint_prefix_frontier
    COMMAND scfg2_vp_case3_endpoint_prefix_frontier_test
  )
  add_test(
    NAME scfg2_vp_case3_left_cache
    COMMAND scfg2_vp_case3_left_cache_test
  )
  add_test(
    NAME scfg2_vp_case3_left_top_width
    COMMAND scfg2_vp_case3_left_top_width_test
  )
  add_test(
    NAME scfg2_vp_case3_tail_top_width
    COMMAND scfg2_vp_case3_tail_top_width_test
  )
  add_test(
    NAME scfg2_be_frame_reuse_profile
    COMMAND scfg2_be_frame_reuse_profile_test
  )
  add_test(
    NAME scfg2_wmbp_be_vp_prefix_profile
    COMMAND scfg2_wmbp_be_vp_prefix_profile_test
  )
  add_test(
    NAME scfg2_wmbp_be_vp_prefix_top_width
    COMMAND scfg2_wmbp_be_vp_prefix_top_width_test
  )
  add_test(
    NAME scfg2_wmbp_vp_phase_barrier
    COMMAND scfg2_wmbp_vp_phase_barrier_test
  )
  add_test(
    NAME scfg2_wmbp_vp_child_beam
    COMMAND scfg2_wmbp_vp_child_beam_test
  )
  add_test(
    NAME scfg2_wmbp_item_beam
    COMMAND scfg2_wmbp_item_beam_test
  )
  add_test(
    NAME scfg2_vp_wmbp_composed_beam
    COMMAND scfg2_vp_wmbp_composed_beam_test
  )
  add_test(
    NAME scfg2_coupled_log_probability_regression
    COMMAND scfg2_coupled_log_probability_regression_test
  )
  add_test(
    NAME scfg2_wmbp_shared_retained_join
    COMMAND scfg2_wmbp_shared_retained_join_test
  )
  add_test(
    NAME scfg2_wmbp_wi_timing
    COMMAND scfg2_wmbp_wi_timing_test
  )
  add_test(
    NAME scfg2_vp_loop_retained_events
    COMMAND scfg2_vp_loop_retained_events_test
  )
  add_test(
    NAME scfg2_vp_case3_nonleft_last_child_index
    COMMAND scfg2_vp_case3_nonleft_last_child_index_test
  )
  add_test(
    NAME scfg2_vp_case12_last_child_index
    COMMAND scfg2_vp_case12_last_child_index_test
  )
  add_test(
    NAME scfg2_vp_recursive_split_retained_events
    COMMAND scfg2_vp_recursive_split_retained_events_test
  )
  add_test(
    NAME scfg2_be_base_stack_retained_events
    COMMAND scfg2_be_base_stack_retained_events_test
  )
  add_test(
    NAME scfg2_be_core_left_retained_events
    COMMAND scfg2_be_core_left_retained_events_test
  )
  add_test(
    NAME scfg2_be_right_tail_retained_core_index
    COMMAND scfg2_be_right_tail_retained_core_index_test
  )
  add_test(
    NAME scfg2_v_retained_agenda
    COMMAND scfg2_v_retained_agenda_test
  )
  add_test(
    NAME scfg2_case3_setup_width
    COMMAND scfg2_case3_setup_width_test
  )
  add_test(
    NAME scfg2_be_right_tail_top_width
    COMMAND scfg2_be_right_tail_top_width_test
  )
  add_test(
    NAME scfg2_be_left_right_frame_top_width
    COMMAND scfg2_be_left_right_frame_top_width_test
  )
  add_test(
    NAME scfg2_right_boundary_dependency
    COMMAND scfg2_right_boundary_dependency_test
  )
  add_test(
    NAME scfg2_right_boundary_schedule
    COMMAND scfg2_right_boundary_schedule_test
  )
  add_test(
    NAME scfg2_right_boundary_item_frontier
    COMMAND scfg2_right_boundary_item_frontier_test
  )
  add_test(
    NAME scfg2_vm_wmp_base_prefix
    COMMAND scfg2_vm_wmp_base_prefix_test
  )
  add_test(
    NAME scfg2_wm_start_prefix
    COMMAND scfg2_wm_start_prefix_test
  )
  add_test(
    NAME scfg2_wip_start_prefix
    COMMAND scfg2_wip_start_prefix_test
  )
  add_test(
    NAME scfg2_vpl_start_prefix
    COMMAND scfg2_vpl_start_prefix_test
  )
  add_test(
    NAME scfg2_vpr_base_prefix
    COMMAND scfg2_vpr_base_prefix_test
  )
  add_test(
    NAME scfg2_wmb_btail_factor
    COMMAND scfg2_wmb_btail_factor_test
  )
  add_test(
    NAME scfg2_approximation_quality
    COMMAND scfg2_approximation_quality_tool
      "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU"
      "....(...)((....)..)((.(.((...))))(.(....)((...))))"
      "1,2,5,10,20,40"
      "${CMAKE_CURRENT_SOURCE_DIR}/params/rna_Turner04.par"
  )
  add_test(
    NAME scfg2_wmb_consumer_view_monotonicity
    COMMAND scfg2_wmb_consumer_view_monotonicity_test
  )
  add_test(
    NAME scfg2_scfg_dependency_gate
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test/tools/check_scfg2_scfg_dependency_gate.sh
  )
endif()

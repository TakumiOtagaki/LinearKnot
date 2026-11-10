add_library(linearknot_contrafold_scores STATIC src/contra_scores.cc)
target_include_directories(linearknot_contrafold_scores PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_include_directories(linearknot_contrafold_scores PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party/linearsampling")
target_compile_features(linearknot_contrafold_scores PUBLIC cxx_std_17)
add_executable(linearknot_contrafold_probe src/dp09_probe.cc)
target_compile_definitions(linearknot_contrafold_probe PRIVATE PK_FIXED_CF_BRANCH=1)
target_link_libraries(linearknot_contrafold_probe PRIVATE linearknot_contrafold_scores LinearKnot::scfg2_exact)
add_executable(linearknot_contrafold_score src/score_structure.cc)
target_link_libraries(linearknot_contrafold_score PRIVATE linearknot_contrafold_scores)
if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME contrafold_exact_backend
      COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test/check_backend.py"
      --probe "$<TARGET_FILE:linearknot_contrafold_probe>"
      --score "$<TARGET_FILE:linearknot_contrafold_score>")
    set_tests_properties(contrafold_exact_backend PROPERTIES TIMEOUT 180)
  endif()
endif()

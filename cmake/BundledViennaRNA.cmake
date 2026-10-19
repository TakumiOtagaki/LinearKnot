# Bundled ViennaRNA build boundary. Source origin, local changes and terms
# are documented in THIRD_PARTY_NOTICES.md and src/ViennaRNA/COPYING.
# Keep the source inventory explicit and its order stable for archive builds.
function(linearknot_add_bundled_vienna)
  include(CheckSymbolExists)
  check_symbol_exists(erand48 "stdlib.h" HAVE_ERAND48)
  check_symbol_exists(strdup "string.h" HAVE_STRDUP)
  configure_file("${PROJECT_SOURCE_DIR}/cmake/vienna_config.h.in"
                 "${PROJECT_BINARY_DIR}/generated/config.h")

  set(vienna_units
    constraints/constraints.c
    constraints/hard.c
    constraints/SHAPE.c
    constraints/soft.c
    datastructures/basic_datastructures.c
    datastructures/char_stream.c
    datastructures/lists.c
    io/file_formats.c
    io/io_utils.c
    landscape/move.c
    loops/external_pf.c
    loops/external.c
    loops/hairpin.c
    loops/internal.c
    loops/multibranch.c
    params/default.c
    params/params.c
    params/io.c
    utils/cpu.c
    utils/higher_order_functions.c
    utils/string_utils.c
    utils/structure_utils.c
    utils/utils.c
    alphabet.c
    boltzmann_sampling.c
    centroid.c
    cofold.c
    commands.c
    dp_matrices.c
    equilibrium_probs.c
    eval.c
    fold_compound.c
    fold.c
    gquad.c
    grammar.c
    MEA.c
    mfe.c
    mm.c
    model.c
    part_func.c
    ribo.c
    sequence.c
    subopt.c
    unstructured_domains.c
  )
  set(vienna_sources)
  foreach(unit IN LISTS vienna_units)
    list(APPEND vienna_sources "${PROJECT_SOURCE_DIR}/src/ViennaRNA/${unit}")
  endforeach()

  add_library(RNA STATIC ${vienna_sources})
  target_include_directories(RNA
    PRIVATE "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src"
    PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>
           $<INSTALL_INTERFACE:include/linearknot-vienna>
  )
  install(TARGETS RNA EXPORT LinearKnotTargets ARCHIVE DESTINATION lib)
  install(FILES "${PROJECT_BINARY_DIR}/generated/config.h"
          DESTINATION include/linearknot-vienna)
endfunction()

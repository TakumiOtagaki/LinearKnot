#!/usr/bin/env bash
set -euo pipefail

: "${SRC_DIR:?conda-build must set SRC_DIR}"
: "${PREFIX:?conda-build must set PREFIX}"

# CMAKE_ARGS is the space-separated set of platform/toolchain options
# supplied by conda's compiler activation. Keep those options before ours.
cmake ${CMAKE_ARGS:-} -S "$SRC_DIR" -B "$SRC_DIR/build/conda-package" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF \
  -DLINEARKNOT_ENABLE_LTO=OFF \
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
cmake --build "$SRC_DIR/build/conda-package" -j8
cmake --install "$SRC_DIR/build/conda-package"

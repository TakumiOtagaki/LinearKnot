# Test Directory Policy

Canonical automated test assets live under `test/`.

Rules:
- Put test inputs, fixtures, helper scripts, and baselines in `test/`.
- Register runnable tests from source-controlled CMake (`CMakeLists.txt`) so fresh builds regenerate the same CTest list.
- Do not rely on stale `build/` artifacts for test discovery or pass/fail decisions.
- Do not introduce an alternative top-level `tests/` directory; use `test/` as the single canonical location.

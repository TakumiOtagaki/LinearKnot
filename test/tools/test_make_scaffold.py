#!/usr/bin/env python3
"""Unit test for scripts/scaffold/make_scaffold.py (no LinearPartition run)."""

import importlib.util
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
spec = importlib.util.spec_from_file_location(
    "make_scaffold", os.path.join(ROOT, "scripts", "scaffold", "make_scaffold.py"))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

sequence = "G" * 22
mea = "((((....))))..((....))"
bpp = {(1, 12): 0.9, (2, 11): 0.9, (3, 10): 0.2, (4, 9): 0.9,
       (15, 22): 0.8, (16, 21): 0.8}

# threshold only: (3,10) is dropped, the rest kept.
assert module.scaffold(sequence, mea, bpp, 0.5, 1) == "((.(....).))..((....))"
# min helix 2: the isolated (4,9) is dropped; (1,12)-(2,11) and the
# second helix (15,22)-(16,21) survive.
assert module.scaffold(sequence, mea, bpp, 0.5, 2) == "((........))..((....))"
# min helix 3: nothing is long enough.
assert module.scaffold(sequence, mea, bpp, 0.5, 3) == "." * 22
print("test_make_scaffold=ok")

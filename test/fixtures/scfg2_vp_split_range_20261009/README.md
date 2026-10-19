# Single maximal G pair does not fix a VP split

`single_child25.txt` has G pairs (1,22) and (11,17). VP(2,25) contains just
the latter maximal pair, yet the full legacy WIP*VPR range has a positive
k=10 term: WIP(3,9)*VPR(10,24). The retired shortcut chose k=11 instead.

At the time of the split-range correction, patched legacy inside gave
7.652115935675325 kcal/mol, compared with 7.652115940547774 before the fix.
Its reflection also lost a positive right-oriented contribution. These
historical values precede the separate multiloop-neighbor correction.
Current patched-legacy energy goldens are in the
[inside-oracle test](../../scfg2_exact_inside_oracle_test.cc).
See the [VP deduction rules](../../../src/scfg2/deductions/generate_exact_basic.hh) for the
full split ranges and the [validation guide](../../../documents/validation.md)
for the maintained checks.

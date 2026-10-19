# Single maximal G pair does not fix a VP split

`single_child25.txt` has G pairs (1,22) and (11,17). VP(2,25) contains just
the latter maximal pair, yet the full legacy WIP*VPR range has a positive
k=10 term: WIP(3,9)*VPR(10,24). The retired shortcut chose k=11 instead.

Patched legacy inside ensemble energy: 7.652115935675325 kcal/mol.
Before correction: 7.652115940547774. Its reflection also loses a positive
right-oriented contribution before correction. See
`documents/eval_structure/scfg2_vp_split_range_correction_20261009.md` and
`results/vp_split_range_20261009/` for source provenance and comparisons.

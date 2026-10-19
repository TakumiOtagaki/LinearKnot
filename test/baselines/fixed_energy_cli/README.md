Baselines in this directory are TSV files used for CLI energy regression tests.

- `density2_mfe_random_20260306_len50_70.tsv` was generated from the
  `compare_mfe_energy_legacy_viterbi_predict_random.sh` strong test with:
  seed `20260306`, 200 cases each for lengths 50 and 70, using
  `worktree_legacy_debug/build/CParty` and `build/scfg_viterbi_predict_tool`.
  One legacy case (`rand50_0133`) failed with
  `backtracking impossible for VM[9, 19]` and is excluded, leaving 399 rows.

# E3NN Protocol Fix Postmortem - 2026-06-04

- Extracted the clean EFG split/centering machinery into `analysis/e3nn_protocol.py`.
- `equiv_t2_e3nn.py` now uses the shared blocked/purged frame split, stores explicit `g_tr`/`g_te`, centers targets by train rows only, normalizes ring invariants from train source rows only, and centers predictions with `center_mask=g_tr`.
- `equiv_t2_backbone_e3nn.py` now uses the same shared blocked/purged split, train-only target centering, train-source-only per-kind feature normalization, and forward-time train-mask prediction centering.
- The broad evidence collector now uses `pack["g_tr"]` instead of `~g_te` so purged rows are not treated as train rows.
- EFG remains the template path; it now imports the shared helper while preserving the same helper names.
- Cheap verification: `/usr/bin/python3 -m py_compile` passed for the shared helper, ring, broad, EFG, and broad evidence scripts.
- Cheap verification: import-only check passed for `e3nn_protocol`, `equiv_t2_e3nn`, `equiv_t2_backbone_e3nn`, and `equiv_t2_efg_e3nn` with `LD_LIBRARY_PATH` set per `analysis/ENV.md`.
- Structural/tiny synthetic assertions passed: train-only atom centering did not use held-out rows; feature normalization used train rows only; blocked split purged neighbouring frames and had zero train/test lag-1 adjacency.
- Structural assertions passed: ring/broad build paths call the shared `make_split_masks`, `centred_by_train_atom`, and `normalised_by_train_rows`; old random split calls and `g_tr = ~g_te` are gone from those paths.
- Structural assertions passed: ring/broad forward paths expose `center_mask` and training calls pass `center_mask=g_tr`.
- E3NN RE-RUN NOT RUN - held per lead; clean-vs-leaky verdict pending.

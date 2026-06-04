# LOAO Fix Postmortem - 2026-06-04
Scope: code-only fix; no DFT/extraction; no full stage-2 regeneration.
Changed `analysis/stage2_law_fits.py` so the main LOAO/between path uses shared `true_loao_*` atom-mean machinery.
LOAO/between now trains on training-atom means, transforms the held-out atom mean with training-atom feature/target stats, and scores physical held-out atom means.
Retired own-heldout-atom centering from main kernel/unified/stage22/stage23 LOAO paths; kept old modulation scorer only under explicit `mislabeled_within_loao_*` historical diagnostic.
Within-frameblock evaluator untouched.
Verification 1: `python3 -m py_compile src/rediscover/analysis/stage2_law_fits.py` passed.
Verification 2: `PYTHONPATH=src/rediscover/analysis python3 -c import stage2_law_fits` passed.
Verification 3: synthetic true-LOAO smoke passed: R2=1.000000000000; held atom 30 prediction component0=12.000000; training mean0=1.666667 vs held own mean0=2.000000.
Verification 4: stage22 LOAO smoke selected_alpha is NaN and bootstrap_replicates=0, confirming no row-level own-centering ridge path.
Verification 5: within-frameblock unchanged vs pre-patch baseline: R2=0.999684707655913; alpha=0.01; coeffs=[1.600794753259740,-0.183435105301205].
FULL RE-FIT NOT RUN — held per lead.

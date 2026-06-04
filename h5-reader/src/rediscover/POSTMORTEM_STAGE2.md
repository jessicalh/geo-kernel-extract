# Stage 2 Postmortem - 2026-06-04

Run dir: `/tmp/rediscover-runs/2026-06-04-stage2-fits`
Substrate: `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4`

Chunk 1 decompose: `allatom_fit_piece2.py` is now a wrapper; active code split into
`allatom_fit_common.py`, `allatom_fit_build3.py`, `allatom_fit_legacy.py`, plus
`stage2_law_fits.py` for the new fits.

Chunk 1 gate: GREEN. Reproduced saved Build-3 artifacts byte-for-byte:
`allatom_fit_score_table.csv`, `allatom_fit_score_long.csv`,
`partition_response_curves.csv`, `partition_response_curves_long.csv`,
`partition_favorable_partitions.csv`, `partition_casehunter_shortlist.csv`.

Basis/integrity: frozen `change_of_basis.get_C()`, `|C.T C-I|max=1.11e-16`; scored
five-component total-T2 only. Python read only emitted Build-4 sidecars and CaseHunter manifests.

Convention ledger: JB/BS/HM kept as current-loop diagnostics; ringchi excluded as separate
convention. Water EFG sign flipped into APBS/MOPAC `-Hessian` convention. Larsen ppm tensors
reported separately from the geometric H-bond shadow and not used in the unified D_ab sum.

Per-kernel verdicts, clean strata selected by CaseHunter atoms + input-side dominance:

| kernel | bucket | coeff | within R2 | LOAO R2 | path agreement | note |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| ring | form-recovered-scale-fitted | +0.693 [0.212,1.174] | 0.275 | 0.165 | 0.009 | real but thin: 5 atoms |
| charge | recovered-law | +9.305 [8.095,10.514] | 0.278 | 0.377 | 0.005 | fixed D_ab shape supported |
| McConnell | can't-make-it-work-for-now | -5.393 [-18.208,7.423] | -0.045 | -0.001 | 0.261 | coefficient not nonzero |
| field | recovered-law | -0.830 [-1.232,-0.427] | 0.033 | 0.034 | 0.087 | weak but nonzero, PySR agrees |
| H-bond | can't-make-it-work-for-now | +24.656 [-45.095,94.408] | -0.048 | -0.013 | 0.050 | geometric/Larsen paths do not rescue |

PySR path: ran from `analysis/venv`; each primary shadow reduced to `shadow * k` (or no useful
law for McConnell), agreeing with ridge/equivariant-Schur coefficients where the law holds.

Unified D_ab-sum: constrained linear fit on through-space atoms
(`max(charge,mc,field,hbond) dominance >=0.5`, ring dominance `<0.7`), 3,903 rows / 25 atoms /
26 terms. Held-out recovery: within-frameblock R2=0.432, LOAO R2=0.258.

Unified intensities vs literature scaling: recovery is real, calibration is not literature-clean.
Charge_total=+15.47 [8.11,22.83], mopac_field_backbone=-8.81 [-12.71,-4.91]; many mc/pq/disp
terms are huge or weakly identified. Treat disagreeing shadows as diagnostic, not averaged truth.

Extended navigable manifests written under `run_dir/equations/{ring,charge,mc,field,hbond}/`.

Disk: final `/tmp` free 158G; `/tmp/rediscover-runs` 6.7G (<15G). Build4 and build1 kept.

# Stage 2.3 Postmortem - 2026-06-04

> **Historical — not current truth (trued 2026-06-04).** This probability
> close record is preserved, but its positive 1P9J LOAO/between values were
> later retracted as between evidence. Use `NOW.md`,
> `POSTMORTEM_TRUE_LOAO_2026-06-04.md`, and
> `MATHS_AUDIT_CHECKLIST_2026-06-04.md` before quoting any transfer claim.

Run dir: `/tmp/rediscover-runs/2026-06-04-stage2_3-probability`
Substrate: `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4`; Build4 CSV/NPY only; frozen `get_C`, `|C.T C-I|max=1.11e-16`; five-component total-T2.
Nulls: within shuffles the DFT target across frames within atom; LOAO shuffles the DFT target across atoms; 1000 shuffles each row.
| recovery | statistical position vs null | determinability | lead-scale placement |
| --- | --- | --- | --- |
| unified_Dab_sum within | R2 0.432; pct 100.0; p 0.001; z 262.8 | attributable to mopac_field_backbone +0.198; mc_category_4 +0.169; mc_category_0 +0.054; charge_total +0.024 | indicative of field + McConnell category mixture, determinable |
| unified_Dab_sum LOAO | R2 0.258; pct 100.0; p 0.001; z 2.0 | mc_category_4 +0.358; mopac_field_backbone +0.338; mc_category_0 +0.209; charge_total +0.067; atom-axis attribution not determined | potentially indicative but INDETERMINATE -> needs more atoms / second-protein atom-axis test |
| charge within | R2 0.278; pct 100.0; p 0.001; z 451.2 | single q/r3 shadow; coefficient CI off zero | indicative of charge q/r3, determinable |
| charge LOAO | R2 0.377; pct 100.0; p 0.001; z 32.6 | single q/r3 shadow; coefficient CI off zero | indicative of charge q/r3, determinable |
| field within | R2 0.033; pct 100.0; p 0.001; z 16.8 | single MOPAC-Coulomb shadow; recovery is 0.03-class | not indicative (~null) |
| field LOAO | R2 0.034; pct 100.0; p 0.001; z 73.8 | single MOPAC-Coulomb shadow; recovery is 0.03-class | not indicative (~null) |
| ring within | R2 0.275; pct 100.0; p 0.001; z 155.1 | single JB ring shadow; thin atom support | indicative of ring current, determinable |
| ring LOAO | R2 0.165; pct 100.0; p 0.001; z 3.2 | single JB ring shadow; thin atom support | indicative of ring current, determinable |
| mc within | R2 -0.045; pct 0.0; p 1.000; z -66.6 | single McConnell shadow; CI spans zero | not indicative (~null) |
| mc LOAO | R2 -0.001; pct 76.4; p 0.237; z 0.5 | single McConnell shadow; CI spans zero | not indicative (~null) |
| hbond within | R2 -0.048; pct 0.0; p 1.000; z -25.0 | single geometric H-bond shadow; CI spans zero | not indicative (~null) |
| hbond LOAO | R2 -0.013; pct 12.4; p 0.876; z -0.9 | single geometric H-bond shadow; CI spans zero | not indicative (~null) |
Cutoff probability curve reports p/R2 at each input-side cut; lower p means farther above the shuffled-target null.
D_ab cut (ring<0.70) within p/R2: 0.30 0.001/0.446; 0.40 0.001/0.446; 0.50 0.001/0.432; 0.55 0.001/0.347; 0.60 0.001/0.598
D_ab cut (ring<0.70) LOAO p/R2: 0.30 0.001/0.199; 0.40 0.001/0.192; 0.50 0.001/0.258; 0.55 0.121/-0.065; 0.60 0.001/0.513
ring cap (D_ab>=0.50) within p/R2: 0.60 0.001/0.476; 0.65 0.001/0.462; 0.70 0.001/0.432; 0.80 0.001/0.479; 0.90 0.001/0.477
ring cap (D_ab>=0.50) LOAO p/R2: 0.60 0.001/0.254; 0.65 0.001/0.267; 0.70 0.001/0.258; 0.80 0.001/0.177; 0.90 0.001/0.206
~0.03-class reads as ~null here; ~0.1-0.2-class entries are placed by determinability, not by size alone.
Artifacts: `stage2_3_recovery_probability.csv`, `stage2_3_unified_drop_one.csv`, `stage2_3_cutoff_probability_curve.csv`, `stage2_3_run_audit.json`.
Disk: `/tmp/rediscover-runs` 6.7G before write (<15G); output drop-old=True; build4+build1 kept.

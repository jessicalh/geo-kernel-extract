# Stage 2.1 Postmortem - 2026-06-04

> **Historical run record — not current truth (trued 2026-06-04).** Preserve as
> provenance. Check `NOW.md`, the corrected top of `STATE.md`, and the June 4
> audit postmortems before quoting quantitative claims.

Run dir: `/tmp/rediscover-runs/2026-06-04-stage2_1-sweep`
Substrate: `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4`
Scope: Build4 CSV/NPY sidecars only; no extraction/DFT/per-source; frozen `get_C`, `|C.T C-I|max=1.11e-16`; five-component total-T2.
Happy-spot shapes are dominance / isolation / driver-modulation; `rises_to_clean` = clean POP, `falls_to_clean` = washout.
- charge: flat / flat / falls_to_clean; clean POP no; support thin_atoms.
- field: flat / flat / rises_to_clean; clean POP modulation; support thin_atoms.
- ring: falls_to_clean / falls_to_clean / flat; clean POP no; support thin_atoms;thin_rows.
- mc: falls_to_clean / flat / falls_to_clean; clean POP no; rescue no; support thin_atoms.
- hbond: flat / falls_to_clean / unavailable; clean POP no; rescue no; support thin_atoms.
- unified_Dab_sum: rises_to_clean / falls_to_clean / rises_to_clean; clean POP dominance,modulation; support thin_atoms.
McConnell/H-bond clean-end rescue criterion: primary coefficient CI off zero at the strict end of an input-side axis; no DFT fit is used to define the axis.
Frame-count ablation uses centered contiguous 20 ps-stride frame blocks; 50 frames is the 1 ns-at-20 ps ubiquitin proxy.
- charge: knee 200 fr; 50fr within 0.171, LOAO 0.398, N_eff 154.6 (rho 0.492); 660fr N_eff 713.1.
- field: knee 50 fr; 50fr within 0.261, LOAO 0.025, N_eff 107.0 (rho 0.026); 660fr N_eff 812.3.
- ring: knee 20 fr; 50fr within 0.213, LOAO 0.027, N_eff 71.4 (rho 0.496); 660fr N_eff 423.6.
- mc: knee no_positive_plateau; 50fr within -0.006, LOAO 0.002, N_eff 167.8 (rho 0.460); 660fr N_eff 1033.3.
- hbond: knee no_positive_plateau; 50fr within -0.149, LOAO -0.017, N_eff 92.8 (rho 0.359); 660fr N_eff 370.1.
- unified_Dab_sum: knee 50 fr; 50fr within 0.415, LOAO 0.104, N_eff 161.5 (rho 0.381); 660fr N_eff 515.8.
50-frame read: use per-kernel rows above; charge/unified positive at 50 means the cheap 1 ns run keeps the main signal, null kernels remain provisional.
Geometric-noise flag: build4 has no noise axis distinct from driver modulation; motion and geometric jitter are not separable here. If clean POP tracks modulation, a future C++ emit should split signal modulation from geometric jitter.
Artifacts: `stage2_1_happy_spot_curves.csv`, `stage2_1_happy_spot_shape_summary.csv`, `stage2_1_frame_count_ablation.csv`, plots_written=True.
Disk: `/tmp/rediscover-runs` 6.7G before write (<15G); build4+build1 kept; output drop-old=True.

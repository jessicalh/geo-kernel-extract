# BUILD 4 Postmortem - 2026-06-03

Run dir: `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4`

Scope completed in C++ only. No git commands were run, and no Python source was added.

Field dominance source is MOPAC-Coulomb, not FF14SB, Amber, or APBS. Per-frame MOPAC charges were not resident; the resident atom-axis `MopacChargeWelfordMean` slab was present, so Build 4 computes per-source `|q_mopac * rhat / r^2|` from that charge and resident geometry in memory.

H-bond exact Larsen per-partner rows were not resident. Build 4 uses resident DSSP/H-bond donor-N/acceptor-O partner geometry plus the H-bond dipolar kernel form for per-partner dominance, with self/same-residue/singularity/near-field guards.

New outputs:
- `per_atom_substrate_features_dominance.npy`
- `per_atom_substrate_dominance_bins.npy`
- `query_results/dominance_fractions_build4_rows.csv`
- `equations/field/cases_manifest.csv`
- `equations/hbond/cases_manifest.csv`
- `equations/charge_wide/cases_manifest.csv`

Legacy `equations/charge/cases_manifest.csv` was kept byte-identical for oracle parity; the widened C/CA charge habitat is banked as `charge_wide`.

Gates:
- Build target `h5reader_extract h5reader_rediscover_tests`: pass.
- Rediscover tests, including deliberate DFT selection-guard leak: 7/7 pass.
- Disk guard: final `/tmp` and `/shared` free 161G; Build 4 3.2G; total rediscover 3.4G.
- Oracle parity: 31 pre-existing CSV/NPY files checked, 0 changed, 0 missing.
- Counts: 846 atoms, 660 DFT frames, 558360 substrate rows.
- Dominance two-path: field max abs delta 5e-13; hbond max abs delta 5e-13.
- New dominance bin IDs: all 5 columns non-degenerate.
- CaseHunter manifests present for ring, charge, mc, field, hbond, charge_wide.

Notes:
- The manifest records the field source as `MOPAC-Coulomb` and charge array as `MopacChargeWelfordMean`.
- APBS remains a grid reference only and is not used for per-source field dominance.
- Build 1 substrate was left in place.

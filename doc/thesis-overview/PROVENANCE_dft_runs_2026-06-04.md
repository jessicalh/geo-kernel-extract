# Provenance of ORCA DFT Runs - 2026-06-04

This is a best-effort forensic provenance pass over the ORCA "orcas" using read-only evidence from local NVMe, the spinner, and scan1/scan2/scan3. I did not infer a DFT setting unless an ORCA input/output echo or generator source showed it. Where I infer lineage, the line is labelled `INFERENCE`.

## 1. STRUCTURE PREP

### 1.1 Main run families found

**1P9J 15 ns ORCA campaign.** The consolidated campaign is explicitly described as "ORCA r2SCAN/def2-SVP NMR shielding calculations on the 1P9J SH3 (BMRB 5801, Wingens 2003) 15 ns MD trajectory" and says jobs were dispatched by `orca-fleet` across batcave + scan1 + scan2 + scan3 (local: `/shared/2026Thesis/1p9j-orcas/README.md`, lines 1-7). The campaign layout records one completed job directory containing `{job_id}.pdb`, `{job_id}.xyz`, `{job_id}_{ts}_nmr.out`, and `{job_id}_meta.json` (local: `/shared/2026Thesis/1p9j-orcas/README.md`, lines 16-26). The same README gives the trajectory source and sampling: 15 ns, 1501 frames at 10 ps, every second frame -> 751 sampled at 20 ps (local: `/shared/2026Thesis/1p9j-orcas/README.md`, lines 101-104).

Evidence example: `/shared/2026Thesis/1p9j-orcas/jobs/batcave_local_15ns_optB_20260501T144807Z_1P9J_5801_f000004_t40.0/..._meta.json` records `protein_id: 1P9J_5801`, `machine: scan2`, `charge: -1`, `atom_count: 846`, `orca_exit_code: 0`, `frame_ps: 40.0`, source machine `scan2:dft_fleet_1P9J_15ns/results`, and source PDB root `/shared/2026Thesis/orca-fleet/tree_1P9J_5801_15ns/pdbs`.

**10-protein calibration DFT set.** The April handoff says 10 proteins were selected because they had successful DFT shielding calculations in the `dft_md` database and were "the smallest proteins in the fleet (335-712 protein atoms)" (local: `/shared/2026Thesis/fleet-ops/SESSION_HANDOFF_20260414.md`, lines 216-229). The collection manifest records campaign `orca-fleet metadynamics calibration`, `job_count: 260`, proteins `1B1V_4292` through `1I8X_4351`, per-machine counts, source DB, and source provenance paths (local: `/shared/2026Thesis/fleet_calibration_dft/collection_manifest.json`, first 40 lines). Example job metadata for `1B1V_4292_P00` records source PDB `/shared/2026Thesis/fleet_calibration-working/1B1V_4292/analysis_output/1B1V_4292_analysis_0ns/1B1V_4292_analysis_0ns.pdb`, `charge: 1`, `atom_count: 335`, and ORCA exit 0 (local: `/shared/2026Thesis/fleet_calibration_dft/1B1V_4292/1B1V_4292_analysis_0ns/1B1V_4292_P00_meta.json`).

**AlphaFold WT/ALA mutation calibration set.** The manifest under `orca-alphafold-and-mutants` records UniProt-like `protein_id` values, `has_wt_out`, `has_ala_out`, WT/ALA charges, and associated `.prmtop`, `.inpcrd`, `.xyz`, `_nmr.out`, `_amber.pdb`, and `_tleap.*` files (local: `/shared/2026Thesis/shielding-calcsets/data/orca-alphafold-and-mutants/manifest.json`, first entries). The same artifact family exists on scan3 under `/home/jessica3/orca_calibration/results/`, including `.out`, `.xyz`, `.prmtop`, and `_tleap.in` files (scan3: `/home/jessica3/orca_calibration/results/...`, bounded listing on 2026-06-04).

### 1.2 Size and "has rings" selection evidence

**Size criterion found.** Two independent local sources document size-based selection:

- 10-protein calibration: "smallest proteins in the fleet (335-712 protein atoms)" (local: `/shared/2026Thesis/fleet-ops/SESSION_HANDOFF_20260414.md`, lines 220-222).
- 200-candidate ORCA triage list: "200 smallest clean proteins with fes_min pose metadata" (local: `/shared/2026Thesis/fleet-ops/SESSION_HANDOFF_20260424.md`, lines 170-173). The candidate CSV has columns `protein_id,n_residues,n_atoms,fes_min_walker,fes_min_time_ps,...` and includes `1P9J_5801` at line 14 (local: `/shared/2026Thesis/fleet-ops/audit_20260423/orca_triage_candidates_200.csv`, lines 1 and 14).

**"Has rings" criterion not found as an explicit persisted predicate.** I found ring-relevant residue composition in `final_choices.csv` rows, e.g. `1B1V_4292` has residue-count JSON containing `PHE`, `TYR`, and `PRO`, while `1P9J_5801` contains `PHE`, `TYR`, `TRP`, `PRO`, and protonated histidines (local: `/shared/2026Thesis/fleet-ops/audit_20260423/final_choices.csv`, rows 6 and 160). However, that is not the same as a documented selector named `has_rings` or "HAS RINGS".

`INFERENCE`: if "has rings" was applied, the likely data basis was residue composition/aromatic-residue presence (`PHE`, `TYR`, `TRP`, `HIS`/`HSP`, possibly `PRO`) in audit/prep tables such as `final_choices.csv`. This is only an inference from residue-count JSON; I did not find a file that states the actual predicate.

### 1.3 Chain extraction, chemistry, protonation, charges

**Old MD/fleet chain and chemistry audit.** The fleet audit says the 685-protein dataset went through a 2026-03 prep pipeline with PubChem compendium, three-agent chemistry assessment, PROPKA protonation, and GROMACS `pdb2gmx` (local: `/shared/2026Thesis/fleet-ops/audit_20260423/AUDIT_STATE_20260423.md`, lines 3-8). The same audit states primary-source snapshots, structural audit, `final_choices.csv`, and recipe/execution agreement checks, including force field/water `charmm36-feb2026_cgenff-5.0` with `tip3p` for 685/685 (local: `/shared/2026Thesis/fleet-ops/audit_20260423/AUDIT_STATE_20260423.md`, lines 42-71).

`1B1V_4292` row evidence: extracted chain `A`, pH 6.0, T 283 K, phosphate 20 mM, PROPKA/disulfide note, protonation source `chemistry.json.protonation_decisions`, force field `charmm36-feb2026_cgenff-5.0`, water `tip3p`, `genion` concentration 0.0224, residue counts, and BMRB/GRO sequence agreement (local: `/shared/2026Thesis/fleet-ops/audit_20260423/final_choices.csv`, row 6).

`1P9J_5801` row evidence: extracted chain `A`, pH 6.3, T 298 K, phosphate 50 mM, `-conc 0.0612 -neutral`, protonation source `protonation.json`, force field `charmm36-feb2026_cgenff-5.0`, water `tip3p`, 3 disulfide pairs, residue counts including 3 `HSP`, and BMRB/GRO sequence agreement (local: `/shared/2026Thesis/fleet-ops/audit_20260423/final_choices.csv`, row 160).

**Charge source.** The ORCA fleet direct prep source states that MD PDBs come from GROMACS MD trajectories with correct protonation and that "Charges come from the GROMACS topology, stored in protein_charges.json" (local: `/shared/2026Thesis/orca-fleet/src/prep.rs`, lines 1-5). It loads `charge` fields from JSON (same file, lines 14-29), skips waters/ions while parsing PDBs (same file, lines 51-70), writes XYZ, warns if a hydrogen-count charge estimate differs, and uses the GROMACS formal charge (same file, lines 204-236). `1p9j-orcas/campaign/protein_charges.json` records `1P9J_5801: charge -1, n_atoms 846` and the 10 calibration proteins including `1B1V_4292: charge 1, n_atoms 335` (local: `/shared/2026Thesis/1p9j-orcas/campaign/protein_charges.json`, first 60 lines). A separate charge audit explains why GROMACS topology charges were preferred: 685 proteins had GROMACS charges and 307 mismatched rounded chemistry estimates (local: `/shared/2026Thesis/dft-pipeline/charge_audit.md`, lines 1-12).

### 1.4 PDB -> ORCA-ready XYZ preparation

There are two prep paths in the code.

**Direct MD/fleet path:** `orca-fleet/src/prep.rs` parses protein `ATOM` records, excludes solvent/ions (`HOH`, `SOL`, `WAT`, `NA`, `CL`, etc.), writes XYZ, uses the GROMACS topology charge, then invokes `orca_caller` with charge and multiplicity 1 (local: `/shared/2026Thesis/orca-fleet/src/prep.rs`, lines 51-70 and 204-252). `dft-pipeline/src/worker.rs` describes the same direct mode as "MD Boltzmann minimum pose: PDB already has correct protonation" and "PDB -> XYZ with formal charge from GROMACS topology -> orca_caller" (local: `/shared/2026Thesis/dft-pipeline/src/worker.rs`, lines 178-210).

**AMBER/tleap path:** `dft-pipeline/src/worker.rs` describes the original calibration path as `pdb4amber -> tleap -> amber_prep.py -> orca_caller` (local: `/shared/2026Thesis/dft-pipeline/src/worker.rs`, lines 222-258). The prep code runs `pdb4amber` with `--nohyd --dry` (local: `/shared/2026Thesis/dft-pipeline/src/prep.rs`, lines 6-25), builds AMBER topology/coordinates with `leaprc.protein.ff14SB` and `leaprc.water.fb3` (same file, lines 38-59), then calls `amber_prep.py` for charge extraction and XYZ output (same file, lines 85-95). `amber_prep.py` sums `NonbondedForce` charges from AMBER topology, rounds to net integer charge, writes a protonated PDB from AMBER topology + coordinates, converts first model to XYZ, and prints `CHARGE=` and `ATOMS=` (local: `/shared/2026Thesis/dft-pipeline/scripts/amber_prep.py`, lines 20-99).

The scan3 WT/ALA calibration artifacts match that AMBER path: `_tleap.in` contains `source leaprc.protein.ff14SB`, `source leaprc.water.fb3`, `loadpdb ..._amber.pdb`, and `saveamberparm ... .prmtop ... .inpcrd` (scan3: `/home/jessica3/orca_calibration/results/A0A101DZD3_WT_tleap.in`, first 5 lines; local mirror example: `/shared/2026Thesis/shielding-calcsets/data/orca-alphafold-and-mutants/A0A822ILJ2/A0A822ILJ2_WT_tleap.in`).

## 2. ORCA / DFT SETTINGS

### 2.1 Strongest evidence: actual ORCA output echoes

**1P9J sample output:** ORCA reports `Program Version 6.1.1`; the echoed input is:

```text
! r2SCAN def2-SVP def2/J NMR CPCM(Water)
%maxcore 12000
%pal nprocs 4 end
%scf MaxIter 300 end
%eprnmr
  TAU DOBSON
* xyzfile -1 1 ...xyz
```

Evidence paths: local `/shared/2026Thesis/1p9j-orcas/jobs/...f000004_t40.0/..._20260508_093046_nmr.out`, lines 54 and 223-230, terminated normally at line 51489; scan2 original `/home/jessica2/dft_fleet_1P9J_15ns/results/batcave_local_15ns_optB_20260501T144807Z_1P9J_5801_f000004_t40.0_20260508_093046_nmr.out`, same echoed settings.

**10-protein calibration sample output:** ORCA 6.1.1; same r2SCAN/def2-SVP/def2/J/NMR/CPCM(Water), `%maxcore 12000`, `%pal nprocs 6`, `%scf MaxIter 300`, `%eprnmr TAU DOBSON`, `* xyzfile 1 1 ...xyz`, terminated normally (local: `/shared/2026Thesis/fleet_calibration_dft/1B1V_4292/1B1V_4292_analysis_0ns/1B1V_4292_P00_20260415_004548_nmr.out`, lines 54 and 223-230, line 21045).

**AlphaFold WT/ALA calibration sample output:** ORCA 6.1.1; same r2SCAN/def2-SVP/def2/J/NMR/CPCM(Water), `%maxcore 8000`, `%pal nprocs 6`, `%scf MaxIter 300`, `%eprnmr TAU DOBSON`, `* xyzfile 2 1 ...xyz`, terminated normally (local: `/shared/2026Thesis/shielding-calcsets/data/orca-alphafold-and-mutants/A0A822ILJ2/A0A822ILJ2_WT_20260401_093917_nmr.out`, lines 54 and 223-230, line 33205; scan3 original `/home/jessica3/orca_calibration/results/A0A822ILJ2_WT_20260401_093917_nmr.out`, same echoed settings).

### 2.2 Generator source

The C generator writes the core input:

```text
! r2SCAN def2-SVP def2/J NMR%s
%maxcore %d
%pal nprocs %d end
%scf MaxIter %d end
%eprnmr
  TAU DOBSON
end
* xyzfile %d %d %s
```

Evidence: `/shared/2026Thesis/1p9j-orcas/campaign/orca_caller.c`, lines 302-311. The same source defines ORCA binary `/opt/orca/orca`, ORCA dir `/opt/orca`, `MAXCORE 12000`, `NPROCS 6`, and `SCF_MAXITER 300` (local: `/shared/2026Thesis/1p9j-orcas/campaign/orca_caller.c`, lines 38-44). The `--no-solvent` option omits `CPCM(Water)` (same file, lines 833-839); actual sampled outputs did not use `--no-solvent`.

**Important mismatch:** the campaign README says the settings include `SlowConv` and `NPROCS 4` "per campaign/orca_caller.c" (local: `/shared/2026Thesis/1p9j-orcas/README.md`, lines 101-105), but the current C source lines 302-311 do not include `SlowConv`, and sampled output echoes also do not include `SlowConv`. The output echo is the stronger evidence for actual run settings. Worker logs show `nprocs: 4` on scan1/scan2/scan3 and `nprocs: 6` on batcave, all with maxcore 12000 and SCF maxiter 300 (local: `/shared/2026Thesis/1p9j-orcas/worker_logs/{scan1,scan2,scan3,batcave}.log`, first ORCA headers).

### 2.3 Per-run variation found

- `nprocs`: observed 4 on scan1/scan2/scan3 1P9J worker logs and outputs; observed 6 on batcave and in the 1B1V calibration sample (local: `/shared/2026Thesis/1p9j-orcas/worker_logs/scan1.log`, line 14; `/shared/2026Thesis/1p9j-orcas/worker_logs/batcave.log`, line 14; `/shared/2026Thesis/fleet_calibration_dft/.../1B1V_4292_P00_20260415_004548_nmr.out`, line 225).
- `maxcore`: 12000 MB for 1P9J and 10-protein calibration samples; 8000 MB for the WT/ALA AlphaFold calibration sample (local paths in section 2.1).
- `charge`: embedded in `* xyzfile <charge> 1`; examples are `-1 1` for 1P9J, `1 1` for 1B1V, `2 1` for A0A822ILJ2 WT (local output paths in section 2.1).
- `SlowConv` / `%scf MaxIter 500`: found only in experimental 1DCJ inputs under `/shared/2026Thesis/dft-pipeline/orca_experiments/1DCJ_4819/exp{1,2,3}_*.inp`, with `SlowConv`, `%scf MaxIter 500`, and nprocs 1/2/6. I did not find those settings in the sampled main completed outputs.

### 2.4 Settings not found

I did not find explicit ORCA grid keywords in the main generator or sampled output echoes (`orca_caller.c`, lines 302-311; output echo lines cited above). I did not find explicit dispersion keywords such as D3/D4 in the main generator or sampled output echoes. Do not state a grid level or dispersion correction in the thesis methods from these files alone.

## 3. ALPHAFOLD / OF3 REFERENCE

### 3.1 RefDB / OF3 corpus sources

The RefDB download README says it is a fresh pull from `https://refdb.wishartlab.com/`, fetched 2026-04-24, with `RefDB_files.tar.gz`, `bmrpdbnew.txt`, 2429 `.str.corr` files, 2427 data rows, and 685 overlap with the thesis fleet (local: `/shared/2026Thesis/of3_breakout/refdb_download/README.md`, lines 1-38). Its intended use was to parse sequences/shifts, filter standard amino acids, submit to OpenFold3, and pair embeddings with re-referenced shifts (same file, lines 61-76).

The OF3 prep terrain notes say the directory maps real OpenFold3 outputs into GROMACS-prep terrain, receives about 685 OF3 outputs, and stores per-protein workspaces under `proteins/<id>/` (local: `/shared/2026Thesis/of3_prep_terrain/CLAUDE.md`, lines 3-8, 18-25, 34-36, 41-61).

The later coherent OF3 input tree says the 685-protein MD-prep working set was rebuilt from `of3-scratch`, with each `of3-scratch/bmr<id>/` subtree copied verbatim, and that per-protein contents include raw BMRB, RefDB, PDB, `bmr<id>/seed_42/` OF3 outputs, config JSONs, `provenance.json`, and `_merge_record.json` with SHA256 of every OF3 file (local: `/shared/2026Thesis/nmr-data/md-prep/data/of3-inputs/_BUILD_PROVENANCE.txt`, lines 11-14 and 30-52). It records 685/685 built and all 685 carrying 10 OF3 manifest-verified files (same file, lines 54-77). The spinner has the same backup lineage under `/mnt/expansion/nmr-data-backup-20260602-pre-coherent-tree/runnable-685/_BUILD_PROVENANCE.txt` and an archive provenance under `/mnt/expansion/nmr-data-archive-20260516/_ARCHIVE_PROVENANCE.txt`.

### 3.2 OF3 per-structure references found

**1B1V_4292 / BMRB 4292.** OF3 reference files exist under `/shared/2026Thesis/nmr-data/md-prep/data/of3-inputs/proteins/bmr4292/`. The query file records seed 42, query `bmr4292`, chain `A`, and sequence `ENFNGGCLAGYMRTADGRCKPTF` (local: `.../bmr4292/inference_query_set.json`, lines 1-15). The same tree contains `bmr4292/seed_42/bmr4292_seed_42_sample_1_model.cif`, raw BMRB files, RefDB `bmr4292.str.corr`, and primary PDB files under `pdb/refdb_primary_1B1V/` (local: bounded file listing under `.../proteins/bmr4292/`). `provenance.json` records OF3 version `0.4.2.dev25+g5ce86ee82`, latent path, latent SHA256, 168 template candidates, `of3_returncode: 0`, and `elapsed_sec: 137.3` (local: `.../bmr4292/provenance.json`, lines 8-26). `_merge_record.json` records `file_count: 10`, model CIF and latent output relative paths, merge script `/shared/2026Thesis/nmr-data/scripts/merge_of3.py`, source dir `/shared/2026Thesis/nmr-data/20260501_path_b/outputs/bmr4292`, and source run `of3-main` (local: `.../bmr4292/_merge_record.json`, lines 3, 6, 21, and 58-62).

Confidence for `1B1V_4292`: high that an OF3 reference exists for BMRB 4292; medium/low that the April ORCA `1B1V_4292_P00` geometry came directly from that OF3 model, because the ORCA job metadata points to `/shared/2026Thesis/fleet_calibration-working/.../1B1V_4292_analysis_0ns.pdb` rather than to the OF3 CIF (local: `/shared/2026Thesis/fleet_calibration_dft/.../1B1V_4292_P00_meta.json`).

**1P9J_5801 / BMRB 5801.** OF3 reference files exist under `/shared/2026Thesis/nmr-data/md-prep/data/of3-inputs/proteins/bmr5801/`. The query file records seed 42, query `bmr5801`, chain `A`, and sequence `VVSHFNDCPLSHDGYCLHDGVCMYIEALDKYACNCVVGYIGERCQYRDLKWWEL` (local: `.../bmr5801/inference_query_set.json`, lines 1-15). The same tree contains `bmr5801/seed_42/bmr5801_seed_42_sample_1_model.cif`, raw BMRB files, RefDB `bmr5801.str.corr`, and primary PDB files under `pdb/refdb_primary_1P9J/` (local: bounded file listing under `.../proteins/bmr5801/`). `provenance.json` records OF3 version `0.4.2.dev25+g5ce86ee82`, latent path, latent SHA256, 69 template candidates, `of3_returncode: 0`, and `elapsed_sec: 168.0` (local: `.../bmr5801/provenance.json`, lines 9-27). `_merge_record.json` records `file_count: 10`, model CIF and latent output relative paths, merge script `/shared/2026Thesis/nmr-data/scripts/merge_of3.py`, source dir `/shared/2026Thesis/nmr-data/20260501_path_b/outputs/bmr5801`, and source run `of3-main` (local: `.../bmr5801/_merge_record.json`, lines 3, 6, 21, and 58-62).

Confidence for `1P9J_5801`: high that an OF3 reference exists for BMRB 5801; medium/low that the 1P9J ORCA frames came directly from that OF3 model, because the ORCA campaign metadata points to a 15 ns MD trajectory source under `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/...` and source PDBs under `/shared/2026Thesis/orca-fleet/tree_1P9J_5801_15ns/pdbs`, not directly to the OF3 CIF (local: `/shared/2026Thesis/1p9j-orcas/README.md`, lines 101-104 and `/shared/2026Thesis/1p9j-orcas/jobs/..._meta.json`).

### 3.3 Newer OF3-to-GROMACS pipeline status

The newer MD-prep pipeline document says the goal is 646 vetted OF3 predicted structures -> GROMACS with AMBER ff14SB/TIP3P (local: `/shared/2026Thesis/nmr-data/md-prep/docs/PIPELINE.md`, lines 1-10). It says the inputs are `runnable-685/proteins/bmr<id>/bmr<id>/seed_42/bmr<id>_seed_42_sample_1_model.cif`, single chain A, canonical residues, heavy atoms only, no chemistry (same file, lines 62-65). It records chemistry/protonation ledger tags `final-2026-05-30`, `maths-v1`, `propka-v1`, `protonate-v1` (same file, lines 67-78). It also says the structure rewrite step was only 3/646 and the full 646 pass had not run (same file, lines 89-95).

Confidence: high for the OF3 source/reference corpus; low that this newer 646 pipeline produced the older March/April/May ORCA outputs, because the document itself says the structure-prep pass had not run at full scale.

### 3.4 AlphaFold WT/ALA accession status

The WT/ALA dataset is named `orca-alphafold-and-mutants` and has UniProt-like IDs and ORCA artifacts (local: `/shared/2026Thesis/shielding-calcsets/data/orca-alphafold-and-mutants/manifest.json`, first entries). I found no AlphaFold DB model accession file, no AFDB URL, and no per-protein AlphaFold source manifest in the local `orca-alphafold-and-mutants` manifest or the scan3 `orca_calibration/results` sample. Treat the UniProt-like IDs as protein identifiers only, not as verified AlphaFold accessions.

## 4. WHERE IT LIVES + CONFIDENCE

**Local NVMe - high confidence for consolidated artifacts.**

- `/shared/2026Thesis/1p9j-orcas/`: consolidated 1P9J campaign, exact local `.out`, `.xyz`, `.pdb`, metadata, worker logs, campaign manifest, source snapshots (local: `/shared/2026Thesis/1p9j-orcas/README.md`, layout and sources sections).
- `/shared/2026Thesis/fleet_calibration_dft/`: 260-job 10-protein calibration collection with per-job `.out`, `.xyz`, `.pdb`, `_meta.json`, and `collection_manifest.json` (local: `/shared/2026Thesis/fleet_calibration_dft/collection_manifest.json`).
- `/shared/2026Thesis/shielding-calcsets/data/orca-alphafold-and-mutants/`: WT/ALA ORCA outputs, `.xyz`, `.prmtop`, `.inpcrd`, `_tleap.in`, `_tleap.log`, and manifest (local: manifest path above).
- `/shared/2026Thesis/orca-fleet/`, `/shared/2026Thesis/dft-pipeline/`, `/shared/2026Thesis/1p9j-orcas/campaign/orca_caller.c`: prep and invocation source code (local code paths cited above).
- `/shared/2026Thesis/nmr-data/md-prep/data/of3-inputs/` and `/shared/2026Thesis/nmr-data/allrefdb/`: OF3/RefDB/BMRB/PDB reference trees (local paths cited above).

**Spinner - medium confidence as historical/backup copy.**

- `/mnt/expansion/nmr-data-backup-20260602-pre-coherent-tree/runnable-685/_BUILD_PROVENANCE.txt`: backup provenance for the 685 OF3 input tree.
- `/mnt/expansion/nmr-data-archive-20260516/_ARCHIVE_PROVENANCE.txt`: archive provenance for `please-last-run-of3`, `of3-working`, and related NMR data trees.
- `/mnt/expansion/orca_calibration_backup/`: older WT/ALA ORCA output logs such as `A0A7C5BHP9_WT_20260312_220746_nmr.out`; useful as antecedent storage, but the local `shielding-calcsets` and scan3 copies provided cleaner sampled setting evidence.
- `/mnt/expansion/nmr-shielding-release-cleanup-20260528/`: mostly analysis/speculative docs; no `.inp` files found at bounded maxdepth 4 during this pass.

**scan hosts - high confidence for original run locations, sampled only.**

- scan1: `/home/jessica1/dft_fleet`, `/home/jessica1/dft_md`, `/home/jessica1/dft_fleet_1P9J_15ns`, `/home/jessica1/orca_calibration` exist. A completed 1P9J output at `/home/jessica1/dft_fleet_1P9J_15ns/results/batcave_local_15ns_optB_20260501T144807Z_1P9J_5801_f000040_t400.0_20260509_001617_nmr.out` echoes ORCA 6.1.1 and the same r2SCAN/def2-SVP/NMR/CPCM(Water) settings.
- scan2: `/home/jessica2/dft_fleet`, `/home/jessica2/dft_md`, `/home/jessica2/dft_fleet_1P9J_15ns`, `/home/jessica2/orca_calibration` exist; `/home/jessica2/dft_md/slowconv_test` exists but no shallow `.inp`/`.out` was found there. A completed 1P9J output at `/home/jessica2/dft_fleet_1P9J_15ns/results/...f000004_t40.0_20260508_093046_nmr.out` echoes the same settings.
- scan3: `/home/jessica3/dft_fleet`, `/home/jessica3/dft_md`, `/home/jessica3/dft_fleet_1P9J_15ns`, `/home/jessica3/orca_calibration` exist. `/home/jessica3/orca_calibration/results/A0A822ILJ2_WT_20260401_093917_nmr.out` echoes ORCA 6.1.1 with `%maxcore 8000`; `/home/jessica3/orca_calibration/results/A0A101DZD3_WT_tleap.in` shows AMBER ff14SB/water.fb3 topology prep.

## GAPS

- I did not find a persisted `has_rings`, "HAS RINGS", or equivalent explicit selection predicate. I found size selection and ring-relevant residue-count evidence, but not the actual selector. Looked in local `of3_prep_terrain`, `of3_breakout`, `dft-pipeline`, `orca-runner`, `orca-fleet`, `1p9j-orcas`, `fleet_calibration_dft`, `nmr-720`, `nmr-data`, `shielding-calcsets`, `fleet-ops`, the spinner cleanup snapshot, and shallow scan directories.
- I did not prove direct lineage from the OF3 CIFs (`bmr4292_seed_42_sample_1_model.cif`, `bmr5801_seed_42_sample_1_model.cif`) to the older April/May ORCA geometries. The OF3 references exist; the ORCA metadata points to MD-derived PDBs and trajectory sources.
- I did not find AlphaFold DB accessions or AFDB URLs for the WT/ALA `orca-alphafold-and-mutants` proteins. The manifest gives UniProt-like IDs and ORCA artifacts, not source model accessions.
- I did not find explicit ORCA grid keywords in the main `.out` echoes or generator source. Do not report a grid setting from this evidence.
- I did not find explicit dispersion keywords in the main `.out` echoes or generator source. Do not report D3/D4/dispersion from this evidence.
- The 1P9J README has stale/conflicting settings (`SlowConv`, `NPROCS 4`) relative to the current C source and actual output echoes. The actual `.out` echoes and worker logs should be cited for settings.
- Remote scan reads were deliberately shallow and bounded. I sampled original run locations and outputs, but did not recursively inventory every scan result directory or read live `.orca_tmp_*` files, because the scans may still be running jobs.

# Fleet Manager State Dump

Last updated: 2026-04-10

## Machines

| Machine | Hostname | Role | GPU | SSH alias | User | Status |
|---------|----------|------|-----|-----------|------|--------|
| batcave | batcave | NAS, dev/test/ML, orchestration | RTX 5090 | local | jessica | UP |
| gotham | gotham | DGX Spark, E3NN training, OF3 | GB10 Blackwell | gotham | jessica | UP |
| scan1 | (VPN) | MD production, MOPAC extraction | RTX 5090 | scan1 | jessica1 | UP |
| scan2 | (VPN) | MD production, MOPAC extraction | RTX 5090 | scan2 | jessica2 | UP |
| scan3 | (VPN) | MD production, MOPAC extraction | RTX 5090 | scan3 | jessica3 | DOWN since ~2026-04-09 |

SSH config in `~/.ssh/config`. NEVER use hardcoded IPs. VPN via iPad, Apple may reset it.

## Current Working Machine

We are on **batcave**. `/shared` is batcave's local NVMe (1.8 TB). The spinner `/mnt/expansion` (5.5 TB) is also local to batcave. Gotham mounts batcave's `/shared` via NFS as `/shared-gotham`.

## The Canonical Data Tree

**`/shared/fleet/clean_fleet/results/`** is king. 685 proteins, all walker data, all poses. This is the single working tree.

### What's in results/{protein}/
```
{protein_id}/
  harvest_receipt.json     # metadata, source machine, pathology, FES minimum
  walker_0/                # COLVAR, HILLS, GRID.dat, md.xtc, md.log, md.gro, reference.pdb
  walker_1/
  walker_2/
  walker_3/
  walker_4/
  run_params/              # reprep_run.sh, prod.mdp, plumed.dat, prod.tpr, machine.txt
  initial-samples/         # pose extractor output (Boltzmann min + dihedral extremes)
    {protein}_01_boltzmann_minimum_{weight}.pdb
    {protein}_02_max_dev_{angle}_{weight}.pdb
    ...
    ensemble.json          # weights, temperature, chain, method description
    dihedrals.tsv          # full dihedral time series
    angle_stats.tsv        # per-angle variance stats
    selection_report.tsv   # why each pose was picked
```

### Counts
- 685 proteins harvested (of 687 in thesis tree, of 696 staged)
- 2 failures never completed MD: 1CKV_4431, 2HA1_7128
- 685 have initial-samples (pose extractor ran on all harvested)
- 0 flagged for pathology
- Backed up to `/mnt/expansion/fleet_results_685/` (checksum verified 2026-04-09)

## Scan Machine Trees

Each scan machine has all 696 proteins staged at:
```
/opt/fleet-reprep1/staging/{protein_id}/
  walker_0/ through walker_4/    # only populated for proteins that ran on THIS machine
  initial-samples/               # deployed 2026-04-09 for ALL 685 proteins on ALL machines
  {protein}_*.pdb/gro/top/...    # prep files
  full_detail.json
```

### Which proteins ran where
- scan1: 228 proteins (results file: `/opt/fleet-reprep1/scan1_reprep_results.txt`)
- scan2: 229 proteins (results file: `/opt/fleet-reprep1/scan2_reprep_results.txt`)
- scan3: 228 proteins (results file: `/opt/fleet-reprep1/scan3_reprep_results.txt`)
- Zero overlap between machines
- All 685 initial-samples deployed to all 3 machines (for MOPAC load balancing)

### GROMACS on scan machines
```
/opt/gromacs/2026.0/bin/gmx_mpi
```
Invoked via mpirun with PLUMED at `/opt/plumed/2.10/lib/libplumedKernel.so`.

## Other Fleet Locations (DO NOT DELETE, document only)

| Path | What | Size | Proteins |
|------|------|------|----------|
| `/shared/fleet/clean_fleet/thesis_tree/` | Older canonical tree, 312 have walker data, 373 don't | 48 GB | 687 |
| `/shared/fleet/clean_fleet/staging/` | Prep artifacts (topology, protonation, run configs) | 2.4 GB | 696 |
| `/shared/fleet/clean_fleet/production/` | Older prep artifacts (subset of staging) | 2.2 GB | 696 |
| `/shared/fleet/clean_fleet/[0-9]*/` | Source PDBs + full_detail.json + tag_alignment.json | 1.6 GB | 696 |
| `/shared/fleet/[0-9]*/` | Original fleet (AF models, NMR chains, equilibration) | ~9 GB | 701 |
| `/shared/fleet_v2/` | Intermediate generation (chain extraction, protonation) | 16 MB | 447 |
| `/shared/fleet_v3/` | Third gen (PROPKA, ionic strength v5, bones) | 1.9 GB | 701 |
| `/shared/final-working-protein-tree/` | Earliest gen, pre-triage, AlphaFold/OpenFold models | 7.7 GB | 1191 |
| `/shared/fleet_scratch/` | Triage logs, batch files, no protein dirs | 182 MB | 0 |
| `/shared/fleet_test/` | Test set | 15 MB | 5 |
| `/shared/2026Thesis/fleet/` | Empty skeleton (proteins/ dir, ops.jsonl) | 8 KB | 0 |
| `/mnt/expansion/fleet_backup_20260326/` | Old pre-harvest backup | 3.9 GB | — |
| `/mnt/expansion/fleet_results_685/` | Verified backup of results/ (made 2026-04-09) | 101 GB | 685 |

## Key Scripts and Tools

### harvest.py — `/shared/fleet/clean_fleet/harvest.py`
Pulls completed MD runs from scan machines into `results/`. Fixed 2026-04-08:
- Dynamic walker count detection (was hardcoded 5)
- Verifies ALL walkers finished (was only walker_0)
- Removed broken gmx trjconv pose extraction (use fes_sampler instead)
- Removed expensive md5 checksumming (rsync handles integrity)
- Added partial-harvest cleanup (clean_partial) for interrupted runs
- Added --partial to rsync for resumability
- Reads results files: `scan{1,2,3}_reprep_results.txt`

Usage:
```bash
python3 harvest.py --status
python3 harvest.py                    # harvest all
python3 harvest.py --machine scan1    # one machine
python3 harvest.py --protein 1UBQ_15410
```

### fes_sampler (C++, Boltzmann pose extraction)
Two copies:
- `/shared/fleet/clean_fleet/fes_sampler/` — original, referenced in CLAUDE.md
- `/shared/2026Thesis/fes-sampler/` — active development copy, built for x86-64 (batcave)
  - Binary: `/shared/2026Thesis/fes-sampler/build/fes_sampler`
  - Test tree: `/shared/2026Thesis/fes-sampler/test_tree/` (10 proteins, full layout)

### Pose extractor (the one that just ran)
Produced `initial-samples/` in all 685 results proteins. Uses dihedral variance selection:
pose 0 = Boltzmann minimum, remaining = thermal extremes of most variable backbone dihedrals.
Typically 5-6 poses per protein.

### Feature extractor (MOPAC-based, NOT YET DEPLOYED)
Lives in `/shared/2026Thesis/nmr-shielding/` (or will). Runs MOPAC on the PDB poses.
~2 hours per protein. Will be deployed to scan machines once ready.

### gmx26.sh — `/shared/fleet/clean_fleet/gmx26.sh`
GROMACS wrapper that strips ORCA from LD_LIBRARY_PATH. Used on scan machines.
Hardcodes `/opt/gromacs/2026.0/bin/gmx_mpi`.

## Pipeline Status

- [x] 687 proteins curated with source + prep + provenance
- [x] 685/687 completed 10ns 5-walker well-tempered metadynamics
- [x] All 685 harvested to batcave results/ (2026-04-09)
- [x] Backed up to spinner (checksum verified)
- [x] Pose extractor ran on all 685 (initial-samples/ with Boltzmann + dihedral poses)
- [x] Poses deployed to all 3 scan machines (2026-04-09)
- [ ] Feature extractor (MOPAC) — awaiting handoff from dev session, then deploy to fleet
- [ ] ORCA DFT shielding — ~400 smallest proteins
- [ ] OpenFold3 embeddings
- [ ] E3NN training
- [ ] Longer metadynamics runs (43ns feasible in 31 days on 3 machines if starting May 1)

## Known Issues

- scan3 VPN down since ~2026-04-09. Data verified locally before dropout.
- 1CKV_4431: failed on scan1, got through prep but never started production. No COLVAR.
- 2HA1_7128: failed, not investigated.
- thesis_tree is 373 proteins behind results/ (not merged, may not need to be).
- Extensive fleet sprawl across /shared — documented above, do not delete.

## Metadynamics Parameters (uniform, all proteins)

| Parameter | Value |
|-----------|-------|
| CVs | RMSD + Radius of gyration |
| Force field | CHARMM36m |
| Height | 1.2 kJ/mol |
| Pace | 500 steps (1 ps) |
| Bias factor | 12 |
| Walkers | 5 per protein (WALKERS_MPI) |
| Production | 10 ns per walker |
| Temperature | Per-protein from BMRB (277-323 K) |
| Trajectory stride | 5000 steps (10 ps), protein atoms only |

## Performance Reference

- 685 proteins: 48-874 ns/day, median 364 ns/day, mean 373 ns/day
- Smallest: 1HD6_4820, 9,614 atoms, 874 ns/day
- Largest: 2GD7_6985, 316,463 atoms, 48 ns/day
- 10ns run wall time: ~7 days across 3 machines
- 43ns run estimate: ~31 days across 3 machines (max feasible for May 1-June 1 window)

## Cost Planning

Scan machines are cloud instances on VPN. Potential surrender date: June 2026.
43 ns is the longest uniform run fitting in 31 days (May 1 - June 1).
40 ns gives 2.4 days buffer for restarts/failures.

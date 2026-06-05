# Methodology: 1P9J Deep-Trajectory Protein

This document records what was done and what is now in hand for the single
deep-trajectory protein, 1P9J. 1P9J is used here as the **within-protein
instrument**: it gives a dense atom-by-frame axis for asking whether geometric
kernels correlate with DFT tensor modulation within one protein. It is not, by
itself, a between-protein transferability instrument. Its between-axis is thin
by construction: one sequence, one fold, one experimental structure family, and
one MD trajectory. The between/transferability axis is therefore reserved for
the 720-WT static corpus.

The central target discipline is tensorial. The DFT NMR tensor is decomposed as

```text
sigma = T0 I + T1 + T2
T0 = Tr(sigma) / 3
T2 = symmetric traceless rank-2 tensor, stored as five components
```

The extractor emits **geometric kernels and features**, not calibrated
shielding predictors. Calibration can turn kernels into shielding estimates;
the methodology here asks whether the emitted objects correlate with the DFT
tensor targets under transparent cuts. It does not assume pointwise equality to
experimental chemical shifts.

## Outline

1. Structure and literature provenance: PDB 1P9J and Wingens et al. 2003.
2. BMRB entry 5801: measured chemical-shift data and sample conditions.
3. RefDB corrected entry: re-referenced shifts and correction metadata.
4. GROMACS MD: force field, water, ensemble, temperature, length, output stride,
   and `.mdp`/`.tpr` settings.
5. ORCA DFT: r2SCAN/def2-SVP NMR calculations with CPCM water.
6. Extractor products: H5 trajectory substrate, topology sidecars, and
   per-atom NPY resources, including available resources not yet evaluated.

## Pipeline

```text
[PDB 1P9J + Wingens 2003 provenance]
                  |
                  v
[BMRB 5801 measured shifts] ---> [RefDB re-referenced shifts]
                  |
                  v
[GROMACS 15 ns NPT trajectory, 846-atom protein]
                  |
                  v
[ORCA r2SCAN/def2-SVP NMR tensor calculations]
                  |
                  v
[nmr_extract trajectory.h5 + topology sidecars]
                  |
                  v
[PerAtomSubstrate: atom-frame tensor targets + geometric kernels]
```

Plain PNG version: `pipeline.png`.

## 1. Structure And Wingens 2003

1P9J is the NMR solution-structure ensemble for T1E, an EGF/TGF-alpha chimera.
The Wingens et al. paper is used here as provenance for the structure,
sequence, sample, restraints, and deposited NMR context. It is not a driver of
the kernel methodology: the present work does not lean on the paper for ring
current, DFT, or geometric-kernel claims.

| Item | Value | Source |
|---|---:|---|
| PDB ID | 1P9J | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/pdb/refdb_primary_1P9J/1P9J_rcsb_metadata.json` |
| BMRB ID linked to structure | 5801 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/pdb/_pdb_index.json` |
| Molecule | T1E, an EGF/TGF-alpha chimera | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/pdb/refdb_primary_1P9J/1P9J.pdb.gz` |
| Experimental method | Solution NMR | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/pdb/refdb_primary_1P9J/1P9J.pdb.gz` |
| Deposited ensemble | 36 models | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/pdb/refdb_primary_1P9J/1P9J_rcsb_metadata.json` |
| Protein length | 54 residues, chain A | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801.prot.fasta` |
| Protein atom count used by the MD/DFT campaign | 846 atoms | `/shared/2026Thesis/1p9j-orcas/campaign/protein_charges.json`; `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| Disulfides | Cys8-Cys22, Cys16-Cys33, Cys35-Cys44 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| Sample conditions in deposited NMR record | pH 6.3, 298 K | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |

| Wingens 2003 provenance item | Recorded value | Source |
|---|---:|---|
| Paper | Wingens et al., J. Biol. Chem. 278, 39114-39123, 2003 | `references-meta/wingens-walma-2003-t1e-egf-tgfa-chimera-summary.txt` |
| DOI | `10.1074/jbc.M305603200` | `references-meta/wingens-walma-2003-t1e-egf-tgfa-chimera-summary.txt` |
| PubMed ID | 12869572 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| NMR temperature | 298 K | `references-meta/wingens-walma-2003-t1e-egf-tgfa-chimera-summary.txt` |
| Spectrometers | 600 and 800 MHz Varian INOVA | `references-meta/wingens-walma-2003-t1e-egf-tgfa-chimera-summary.txt`; `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/pdb/refdb_primary_1P9J/1P9J.pdb.gz` |
| Structural restraints | 660 NOEs, 9 H-bond restraints, 98 TALOS phi/psi restraints, 3 disulfides | `references-meta/wingens-walma-2003-t1e-egf-tgfa-chimera-summary.txt` |
| Refinement software | X-PLOR 3.851 with CHARMM22 water refinement | `references-meta/wingens-walma-2003-t1e-egf-tgfa-chimera-summary.txt`; `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/pdb/refdb_primary_1P9J/1P9J.pdb.gz` |

The sequence used by the data record is:

```text
VVSHFNDCPLSHDGYCLHDGVCMYIEALDKYACNCVVGYIGERCQYRDLKWWEL
```

Source: `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801.prot.fasta`.

## 2. BMRB Measured Chemical Shifts

The measured chemical shifts are taken from BMRB entry 5801. The local BMRB
mirror records the entry fetch and checksum metadata, and the NMR-STAR entry
contains both the experimental chemical-shift loop and the sample/reference
metadata.

| BMRB item | Value | Source |
|---|---:|---|
| Entry ID | 5801 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| Submission/accession date | 2003-05-19 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| Release date | 2003-10-13 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| NMR-STAR version | 3.1.1.61 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| Local fetch date | 2026-05-14 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/_fetch_record.json` |
| Raw BMRB files fetched | 6 files, 249843 bytes total, MD5 check true | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/_fetch_record.json` |
| Sample | 0.8 mM T1E, U-15N labelled; 50 mM phosphate; 95% H2O / 5% D2O | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| Conditions | pH 6.3 +/- 0.2; 298 +/- 1 K | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| Chemical reference | 1H DSS methyl protons at 0.0 ppm, internal direct; 15N indirect DSS ratio 0.101329118 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |

| Measured BMRB shift inventory | Count | Source |
|---|---:|---|
| Raw atom chemical-shift rows parsed from the BMRB atom-shift loop | 402 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| 1H chemical shifts listed in BMRB entry statistics | 346 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| 15N chemical shifts listed in BMRB entry statistics | 56 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| Residues represented in the raw atom-shift rows | 54 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| T1 relaxation values listed in BMRB entry statistics | 52 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |
| Coupling constants listed in BMRB entry statistics | 51 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/bmrb/bmr5801_3.str` |

These measured shifts provide the experimental NMR reference frame for the
protein. In the present tensor-kernel work, they are not used as pointwise
targets for the per-frame geometric kernels. The measured shifts and RefDB
corrections instead define the deposited NMR data provenance around the 1P9J
structure.

## 3. RefDB Re-Referenced Shifts

The RefDB file is the local corrected/re-referenced counterpart of BMRB 5801.
It records that correction was performed using PDB structure `1P9JA` and gives
the shift-offset summary. RefDB's corrected loop is used as a curated chemical
shift resource alongside the raw BMRB entry.

| RefDB item | Value | Source |
|---|---:|---|
| Local corrected file | `bmr5801.str.corr` | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| Local extraction date | 2026-05-14 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/_fetch_record.json` |
| RefDB source archive | `https://refdb.wishartlab.com/RefDB_files.tar.gz` | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/_fetch_record.json` |
| BMRB accession in corrected file | 5801 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| PDB used for correction | 1P9JA | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| Corrected atom-shift rows parsed from RefDB loop | 354 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| Atom types in corrected loop | 298 H rows; 56 N rows | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| Residues represented in corrected loop | 54 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |

RefDB defines the corrected value as:

```text
Observed* = Observed + Offset correction
```

| RefDB correction summary | Value | Source |
|---|---:|---|
| 15N offset added to original shifts | -1.07 ppm | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| Mean corrected-minus-predicted difference, N | -1.07 ppm | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| Mean corrected-minus-predicted difference, HN | -0.17 ppm | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| Mean corrected-minus-predicted difference, HA | 0.02 ppm | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| RefDB N correlation | 0.522 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| RefDB HN correlation | 0.638 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| RefDB HA correlation | 0.831 | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| RefDB N RMSD after correction | 1.806 ppm | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| RefDB HN RMSD after correction | 0.323 ppm | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |
| RefDB HA RMSD after correction | 0.152 ppm | `/shared/2026Thesis/nmr-data/allrefdb/bmr5801/refdb/bmr5801.str.corr` |

The RefDB record also includes 52 T1 rows, 52 T1rho rows, and 52 heteronuclear
NOE rows. These are available NMR observables around the same protein system,
but the current per-atom tensor substrate is built from the MD/DFT trajectory
axis rather than fitting those relaxation records.

## 4. GROMACS Molecular Dynamics

The production trajectory used for the 1P9J deep-trajectory work is the
`batcave_local_15ns_optB_20260501T144807Z` GROMACS run. The MD was prepared
from the 1P9J protein structure using Amber ff14SB, explicit TIP3P water, ions,
equilibration, and a 15 ns production NPT simulation at 298 K and 1 bar.

The production ensemble is the usual constant-composition NPT target:

```text
pi(q,p,V) proportional to exp[-beta (K(p) + U(q,V) + P_ext V)]
N = fixed atom composition, T = 298 K, P_ext = 1 bar
```

The production length follows directly from the `.mdp`:

```text
t = nsteps * dt = 7,500,000 * 0.002 ps = 15,000 ps = 15 ns
```

| MD item | Value | Source |
|---|---:|---|
| Production directory | `batcave_local_15ns_optB_20260501T144807Z` | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/` |
| GROMACS version | 2026.0 | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/prep_run_record.json`; `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.run_record.json` |
| Protein force field | Amber ff14SB (`amber14sb.ff/forcefield.itp`) | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/topol.top`; `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/prep_run_record.json` |
| Water model | TIP3P topology (`amber14sb.ff/tip3p.itp`) | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/topol.top` |
| Solvent coordinate pack used during solvation | `spc216.gro` | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/prep_run_record.json` |
| Box | Dodecahedron, 2.00 nm padding | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/prep_run_record.json`; `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/dispatch_record.json` |
| Ions | 21 NA, 20 CL; neutralized at 0.0612 M requested concentration | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/topol.top`; `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/prep_run_record.json` |
| Total production system size | 52337 atoms | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/dispatch_record.json` |
| Production ensemble | NPT, 298 K, 1 bar | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| Thermostat | V-rescale, groups Protein and Non-Protein, tau_t 0.1 ps, ref_t 298/298 K | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| Barostat | C-rescale, isotropic, tau_p 5.0 ps, ref_p 1.0 bar | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| Integrator and timestep | `md`, dt 0.002 ps | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| Steps and length | 7,500,000 steps; 15 ns | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp`; `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/dispatch_record.json` |
| Constraints | H-bonds constrained with LINCS | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| Nonbonded scheme | Verlet; PME electrostatics; rcoulomb 1.0 nm; rvdw 1.0 nm; DispCorr EnerPres | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| `.tpr` generation | `gmx_mpi grompp -f production.mdp -c relax.gro -t relax.cpt -p topol.top -o production.tpr -maxwarn 1` | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/mdout.mdp` |
| Production completion | Exit OK; 2026-05-01 | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.run_record.json` |

| MD phase | Settings | Source |
|---|---|---|
| `pdb2gmx` | `-ff amber14sb -water tip3p -ter -his -ignh` | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/prep_run_record.json` |
| NVT equilibration | 100000 steps, dt 0.002 ps, V-rescale, 298 K, position restraints | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/prep_run_record.json`; NVT `.mdp` in the same prep directory |
| NPT equilibration | 250000 steps, dt 0.002 ps, V-rescale + C-rescale, 1 bar, position restraints | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/prep_run_record.json`; NPT `.mdp` in the same prep directory |
| Unrestrained relaxation | 1,000,000 steps, dt 0.002 ps, 2 ns | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/relax_2.0ns_20260501T142819Z/relax.mdp` |
| Production | 7,500,000 steps, dt 0.002 ps, 15 ns | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |

| Production output | Stride in `.mdp` | Physical stride | Source |
|---|---:|---:|---|
| TRR positions and velocities | `nstxout = 5000`, `nstvout = 5000` | 10 ps | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| XTC compressed positions | `nstxout-compressed = 2500` | 5 ps | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| Energies | `nstenergy = 500` | 1 ps | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| Log | `nstlog = 5000` | 10 ps | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| Forces | `nstfout = 0` | not written | `/shared/2026Thesis/nmr-shielding/tests/data/fleet_amber/1P9J_5801/prep_run_20260501T141627Z/batcave_local_15ns_optB_20260501T144807Z/production.mdp` |
| Extractor/ORCA sampling axis | every second 10 ps frame, including frame 0 | 751 frames at 20 ps over 15 ns | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json`; `/shared/2026Thesis/1p9j-orcas/README.md` |

## 5. ORCA DFT Settings

The ORCA campaign used the MD frames as protein-only DFT geometries. The local
campaign documentation and representative ORCA output agree on the calculation
level:

```text
! r2SCAN def2-SVP def2/J NMR CPCM(Water)
%maxcore 12000
%scf MaxIter 300 end
%eprnmr TAU DOBSON
* xyzfile -1 1 <frame>.xyz
```

| ORCA item | Value | Source |
|---|---:|---|
| Campaign root | `/shared/2026Thesis/1p9j-orcas/` | `doc/thesis-overview/PROVENANCE_dft_runs_2026-06-04.md`; `doc/thesis-overview/orcas-720/WRITEUP.md` |
| ORCA version | 6.1.1 | `doc/thesis-overview/PROVENANCE_dft_runs_2026-06-04.md`; representative `_nmr.out` under `/shared/2026Thesis/1p9j-orcas/jobs/` |
| Functional | r2SCAN | `doc/thesis-overview/PROVENANCE_dft_runs_2026-06-04.md`; `/shared/2026Thesis/1p9j-orcas/campaign/orca_caller.c` |
| Basis | def2-SVP | `doc/thesis-overview/PROVENANCE_dft_runs_2026-06-04.md`; `/shared/2026Thesis/1p9j-orcas/campaign/orca_caller.c` |
| Auxiliary basis / RI keyword | def2/J | `doc/thesis-overview/PROVENANCE_dft_runs_2026-06-04.md`; representative `_nmr.out` under `/shared/2026Thesis/1p9j-orcas/jobs/` |
| Solvation | CPCM(Water) | `doc/thesis-overview/PROVENANCE_dft_runs_2026-06-04.md`; representative `_nmr.out` under `/shared/2026Thesis/1p9j-orcas/jobs/` |
| NMR keyword | `NMR` with `%eprnmr TAU DOBSON` | `doc/thesis-overview/PROVENANCE_dft_runs_2026-06-04.md`; `/shared/2026Thesis/1p9j-orcas/campaign/orca_caller.c` |
| Charge and multiplicity | charge -1, multiplicity 1 | `/shared/2026Thesis/1p9j-orcas/campaign/protein_charges.json`; representative `_nmr.out` under `/shared/2026Thesis/1p9j-orcas/jobs/` |
| Protein atoms per ORCA geometry | 846 | `/shared/2026Thesis/1p9j-orcas/campaign/protein_charges.json` |
| SCF iteration cap | 300 | `/shared/2026Thesis/1p9j-orcas/campaign/orca_caller.c`; representative `_nmr.out` under `/shared/2026Thesis/1p9j-orcas/jobs/` |
| Per-job files | `.pdb`, `.xyz`, `_nmr.out`, `_meta.json` | `doc/thesis-overview/PROVENANCE_dft_runs_2026-06-04.md` |
| Targeted trajectory frames | 751 sampled frames at 20 ps, including frame 0 | `/shared/2026Thesis/1p9j-orcas/README.md`; `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| Current status-file inventory read for this write-up | 751 job rows; 730 completed job directories recorded | `/shared/2026Thesis/1p9j-orcas/status/jobs.csv`; `/shared/2026Thesis/1p9j-orcas/jobs/` |

No explicit grid or dispersion keyword is asserted here, because the cited
input generator and representative ORCA output do not record one. The tensor
targets used downstream come from the ORCA NMR tensor output; the extractor
then decomposes and aligns those tensor targets for comparison with geometric
kernels.

## 6. Extractor Products And Available Resources

The baseline extraction run is:

```text
/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/
```

It ran `nmr_extract` on the 1P9J trajectory with trajectory, MOPAC, vacuum
Coulomb, topology-sidecar, and per-frame substrate outputs. The run completed
successfully after 7 h 48 m and wrote a 751-frame, 846-atom H5 trajectory file.

| Extraction item | Value | Source |
|---|---:|---|
| Run label | `run_01_baseline_2026-05-21_sha7d8dbe9` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| Extractor | `nmr_extract` v0.2.0 | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/extraction_manifest.json` |
| Git SHA | `7d8dbe9` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| Configuration | `FullFatFrameExtraction (PerFrame + MOPAC + vacuum Coulomb)` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| Completion | success, 2026-05-22T07:43:00+0100 | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| Wall time | 7 h 48 m | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| H5 output | `trajectory.h5`, 751 frames, 846 atoms | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| H5 manifest size | 1.8 GB on disk | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| Actual H5 trajectory groups | 58 child groups under `/trajectory` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| AIMNet2 model attached | `/shared/2026Thesis/nmr-shielding/data/models/aimnet2_wb97m_0.jpt` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |

### Topology Sidecars

The topology sidecars establish the stable axes used by the H5 and NPY outputs.

| Sidecar | Shape / inventory | Source |
|---|---:|---|
| `atoms_category_info.npy` | 846 atom records, including element, residue, naming, role, formal charge, aromatic/exchangeable flags, and force-field atom type | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/atoms_category_info.npy`; `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/extraction_manifest.json` |
| `residues.npy` | 54 residue records | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/residues.npy`; `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/extraction_manifest.json` |
| `bonds.npy` | 862 bond records | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/bonds.npy`; `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/extraction_manifest.json` |
| `rings.npy` | 16 rings: 15 aromatic, 1 saturated | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/rings.npy`; `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/extraction_manifest.json` |
| `ring_membership.npy` | 96 ring-membership records | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/ring_membership.npy`; `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/extraction_manifest.json` |
| `extraction_manifest.json` | Axis contract and topology provenance | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/extraction_manifest.json` |
| `manifest.json` | Run configuration, wall time, and H5 summary | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/manifest.json` |
| `run.log`, `run.err`, `run.pid`, `io.mc` | Execution logs and run-time support files | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/` |

The sidecar alignment contract is explicit: atom NPYs share a common atom-row
order; residues, bonds, rings, and ring memberships each use their canonical
axes. Source: `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/extraction_manifest.json`.

### H5 Trajectory Inventory

The H5 file contains atom metadata and 58 trajectory child groups. The table
below lists the complete group inventory by function. Dataset shapes are
representative of the stored axes; full dataset names are in `trajectory.h5`.

| H5 group family | Groups / stored resources | Source |
|---|---|---|
| Atom metadata | `/atoms/element`, `/atoms/pdb_atom_name`, `/atoms/residue_index` with 846 entries | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| Frame and coordinate axis | `frames`, `positions` with positions `(846, 751, 3)` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| Ring-current and ring-path tensors | `bs_shielding_time_series`, `bs_welford`, `bs_t0_autocorrelation`, `hm_shielding_time_series`, `hm_welford`, `ringchi_shielding_time_series`, `piquad_shielding_time_series`, `disp_shielding_time_series`, `ring_neighbourhood_trajectory_stats`, `ring_pucker_time_series` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| McConnell / bond anisotropy resources | `mc_shielding_time_series`, `mc_welford`, `mopac_mc_shielding_time_series`, `mopac_bond_order_welford`, `bond_length_stats` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| Electrostatic field / EFG resources | `apbs_efield_time_series`, `apbs_efg_time_series`, `mopac_coulomb_shielding_time_series`, `mopac_vs_ff14sb_reconciliation`, `water_field_time_series`, `water_field_welford` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| Charge resources | `aimnet2_charge_time_series`, `aimnet2_charge_response_gradient_time_series`, `aimnet2_charge_response_gradient_welford`, `mopac_charge_welford`, `eeq_welford` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| AIMNet2 representation | `aimnet2_embedding_time_series` with embedding `(846, 751, 256)` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| Hydrogen-bond resources | `hbond_shielding_time_series`, `hbond_count_welford`, `larsen_hbond_1pHB_shielding_time_series`, `larsen_hbond_2pHB_shielding_time_series`, `larsen_hbond_1pHaB_shielding_time_series`, `larsen_hbond_2pHaB_shielding_time_series`, `larsen_hbond_count_time_series`, `larsen_hbond_water_term_time_series` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| Solvent exposure and hydration | `sasa_time_series`, `sasa_welford`, `hydration_shell_time_series`, `hydration_shell_welford`, `hydration_geometry_time_series`, `hydration_geometry_welford` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| Backbone/conformation conditioners | `dihedral_time_series`, `dihedral_bin_transition`, `dssp8_time_series`, `dssp8_transition`, `j_coupling_time_series`, `rmsd_tracking` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| GROMACS energy resources | `gromacs_energy_time_series`, `bonded_energy_time_series` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| Tripeptide/Larsen-style resources | `tripeptide_bb_shielding_time_series`, `tripeptide_bb_residual_vec_time_series`, `tripeptide_bb_method_tag_time_series`, `tripeptide_neighbor_shielding_time_series`, `tripeptide_neighbor_residual_vec_prev_time_series`, `tripeptide_neighbor_residual_vec_next_time_series` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |
| Selection/QC resources | `selections` | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5` |

Several H5 group names contain the legacy word `shielding`; in this methodology
they are treated according to their actual role. The ORCA output is a DFT NMR
tensor target. The classical time-series groups are geometric or semi-empirical
kernel carriers whose calibration and evaluation are separate.

### PerAtomSubstrate Build

The current per-atom aggregate substrate inventory cited here is:

```text
/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/
```

It emits one row per `(frame_slot, atom_index)` and uses the alignment contract:

```text
row_id == frame_slot * n_atoms + atom_index
sidecar row i == rows.csv row_id i
```

| PerAtomSubstrate item | Value | Source |
|---|---:|---|
| Relationship | `per_atom_substrate`, `per_atom_aggregate` | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json` |
| Atoms | 846 | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json` |
| DFT-aligned frames | 660 | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json` |
| Rows | 558360 | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json`; `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_rows.csv` |
| T2 component order | `[xy, yz, zz, xz, xx-yy]` | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_column_specs.json` |
| Tensor alignment check | max angle 0.000237 deg; max RMSD 0.000513 A; T2 frame-aligned; T1 frame verified | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json` |
| Normalization | raw lab frame | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json` |
| Ring cutoff | 8 A | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json`; `h5-reader/src/rediscover/PerAtomSubstrate.h` |
| Bond cutoff | 10 A | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json`; `h5-reader/src/rediscover/PerAtomSubstrate.h` |
| Charge cutoff | 6 A | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json`; `h5-reader/src/rediscover/PerAtomSubstrate.h` |
| McConnell near-field ratio | 0.5 | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json`; `h5-reader/src/rediscover/PerAtomSubstrate.h` |
| AIMNet2 embedding | separable, 256 dimensions | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json`; `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_aimnet2_embedding.npy` |

| Per-atom NPY / CSV output | Shape or rows | Role |
|---|---:|---|
| `per_atom_substrate_rows.csv` | 558360 rows, 71 columns | Row index, atom identity, frame slot, original frame index, time, residue/role labels, strata, and metadata |
| `per_atom_substrate_ring_identity.csv` | 96 rows, 10 columns | Ring membership identity sidecar |
| `per_atom_substrate_column_specs.json` | 641 column specifications | Array/column metadata, mechanisms, irreps, units, T2 order |
| `per_atom_substrate_manifest.json` | manifest | Relationship, cutoffs, DFT alignment, support rows, sidecar list |
| `per_atom_substrate_target_T0.npy` | `(558360,)`, float64 | DFT scalar target |
| `per_atom_substrate_target_T1.npy` | `(558360, 3)`, float64 | DFT antisymmetric diagnostic vector |
| `per_atom_substrate_target_T2.npy` | `(558360, 5)`, float64 | Primary DFT rank-2 tensor target |
| `per_atom_substrate_target_dia_T0.npy` | `(558360,)`, float64 | Diamagnetic T0 diagnostic |
| `per_atom_substrate_target_dia_T1.npy` | `(558360, 3)`, float64 | Diamagnetic T1 diagnostic |
| `per_atom_substrate_target_dia_T2.npy` | `(558360, 5)`, float64 | Diamagnetic T2 diagnostic |
| `per_atom_substrate_target_para_T0.npy` | `(558360,)`, float64 | Paramagnetic T0 diagnostic |
| `per_atom_substrate_target_para_T1.npy` | `(558360, 3)`, float64 | Paramagnetic T1 diagnostic |
| `per_atom_substrate_target_para_T2.npy` | `(558360, 5)`, float64 | Paramagnetic T2 diagnostic |
| `per_atom_substrate_features_classical.npy` | `(558360, 89)`, float64 | Classical geometric and charge features |
| `per_atom_substrate_features_ring_paths.npy` | `(558360, 226)`, float64 | Ring-current, ring EFG, pi-quadrupole, dispersion, and ring-geometry paths |
| `per_atom_substrate_features_method_paths.npy` | `(558360, 111)`, float64 | Method-comparison paths, electrostatic field/EFG, charge paths, and related method carriers |
| `per_atom_substrate_features_hbond_conditioning.npy` | `(558360, 73)`, float64 | Larsen/DSSP/H-bond conditioning and backup channels |
| `per_atom_substrate_features_conditioning.npy` | `(558360, 32)`, float64 | Secondary structure, dihedral, ring geometry, and general conditioners |
| `per_atom_substrate_features_dominance.npy` | `(558360, 10)`, float64 | Dominance and isolation scalars |
| `per_atom_substrate_partition_bins.npy` | `(558360, 25)`, int32 | Partition-bin IDs |
| `per_atom_substrate_dominance_bins.npy` | `(558360, 5)`, int32 | Dominance-bin IDs |
| `per_atom_substrate_driver_modulation_by_atom.npy` | `(846, 9)`, float64 | Per-atom driver modulation summaries |
| `per_atom_substrate_backbone_audit.npy` | `(558360, 14)`, float64 | Backbone audit and tensor/QC sidecar |
| `per_atom_substrate_aimnet2_embedding.npy` | `(558360, 256)`, float32 | AIMNet2 embedding resource |

Source for the NPY and CSV inventory:
`/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/`.

The manifest records full support for the currently emitted complete-row
channels: charge-complete rows, FF14SB charge, MOPAC Welford mean charge,
MOPAC/field resources, ring-path resources, SASA, hydration, water field,
H-bond count, H-bond kernel resources, dispersion, pi-quadrupole, HM, and
ringchi resources all have 558360 present rows where listed as present. Source:
`/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json`.

### Resources Now In Hand Beyond The Current Evaluated Cuts

The extractor output is deliberately richer than the cuts already evaluated.
These resources are available for later analysis without changing the core
provenance of the 1P9J trajectory:

| Available resource | What it contains | Source |
|---|---|---|
| Dia/para target splits | Separate T0/T1/T2 diamagnetic and paramagnetic target arrays; dia+para consistency was recorded in the rediscover state notes | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_target_dia_*.npy`; `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_target_para_*.npy`; `h5-reader/src/rediscover/STATE.md` |
| T1 diagnostic | Three-component antisymmetric diagnostic retained as a completeness/QC channel, not collapsed into a scalar shift claim | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_target_T1.npy`; `h5-reader/src/rediscover/STATE.md` |
| AIMNet2 embedding | 256-dimensional per-atom/per-frame embedding, separable sidecar, available for Stage 3 or T0/T2 comparative analysis | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_aimnet2_embedding.npy`; `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_manifest.json` |
| Ring-current path family | BS, Haigh-Mallion, ringchi, JB-style, per-type, total-B, count, ring EFG, pi-quadrupole, and dispersion channels | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_features_ring_paths.npy`; `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_column_specs.json`; `h5-reader/src/rediscover/STATE.md` |
| Field/EFG agreement family | MOPAC field/EFG, APBS field/EFG, AIMNet2 EFG/charge-response, water EFG/field resources | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5`; `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_features_method_paths.npy` |
| H-bond and Larsen family | Best/backup H-bond channels, donor/acceptor Larsen classes, DSSP backup, and hbond scalar support | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5`; `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_features_hbond_conditioning.npy`; `h5-reader/src/rediscover/STATE.md` |
| Conditioners | DSSP8, chi/omega/dihedral, pyramidalization, ring geometry, SASA, hydration, partition and dominance bins | `/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/trajectory.h5`; `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_features_conditioning.npy` |
| Driver modulation | Per-atom driver modulation summaries and dominance/isolation bins | `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_driver_modulation_by_atom.npy`; `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build4/per_atom_substrate_features_dominance.npy` |
| Stage 2 analysis scripts | Downstream analysis resources for convergence, overlays, essential dynamics, and variance comparisons | `/mnt/expansion/nmr-shielding-release-cleanup-20260528/learn/stage2/` |

The current state notes record the interpretation boundary for 1P9J: within
results stand as a positive use of the dense trajectory axis, while true
between-atom or between-protein claims are not assigned to this single-protein
instrument. Within 1P9J, the state notes record charge and ring-current
geometric kernels at about R2 0.28 and a unified within combine at about R2
0.43. These are within-protein correlation results, not proof of
between-protein transfer. Source: `h5-reader/src/rediscover/STATE.md`.

That framing is the strength of this system. The MD/DFT/extractor stack gives a
high-density, tensor-preserving within-protein substrate, and the emitted
resources retain the richer physics carriers needed for later work instead of
collapsing the problem to a scalar or a single mechanism.

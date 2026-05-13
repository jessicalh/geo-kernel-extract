# MD prep session handoff — two sample runs

**Date:** 2026-04-30
**Status:** Standalone. Self-contained. Scope is GROMACS / MDP only — what the prep session needs to know to set up and run MD for the two sample proteins.

H5 schema, extractor config, calculator decisions, and downstream analysis are handled by separate sessions and are explicitly **not** the prep session's concern.

---

## Section 1 — What the two sample runs are for

The 685-protein fleet was stopped after evaluation showed bad chain extractions in the structure preparation. The recovery path for the full-fleet redispatch will use OF3-generated structures direct from sequence.

**These two sample runs are not on OF3 structures.** They are on the same two BMRB-PDB-extracted structures that produced the first batch — known good for these specific proteins, with published S² / T1 / T2 NMR relaxation data and BMRB-assigned chemical shifts.

For these 15 ns × 1-replica runs the testable comparison downstream is **per-frame kernel T2 contributions vs sparse DFT anchors and per-frame agreement with BMRB shifts**. Published S² / T1 / T2 are the validation target for the calibration-µs scans (1Z9B 2 µs and 1P9J 5 µs already in flight), not for these 15 ns runs. The 15 ns sample runs are MD-setup validation plus a chemical-shift-level calculator check, not relaxation-observable validation.

Their purpose for the prep session: validate the MD setup (MDP, ions, pH, force field choices, run shape) before the full 685-fleet redispatch runs against OF3 structures with the same setup.

```
Prep session's gate     The sample runs land cleanly with the MDP and
                          ion / pH setup below. Trajectories produce
                          well-formed XTC / EDR / TRR / cpt outputs at
                          the expected size and wallclock.

What the prep session    Whether the H5 extraction downstream produces
does NOT need to gate    valid features, whether calculators land on
                          target — those are downstream sessions.
```

---

## Section 2 — Force field, water, ions, pH

```
Force field            ff14SB
Water model            TIP3P

Ions                   Neutralise the box; standard physiological
                         ionic strength (NaCl ~0.15 M or what the
                         prep_toolchain canonical method produces).
                         Whatever the toolchain produces for the OF3
                         path should be the same as what's used here
                         for consistency. The prep session is the
                         authority on the canonical ion / pH workflow;
                         this brief defers to that workflow.

pH / protonation       Per the prep_toolchain canonical protocol.
                         Titratable residues (HIS / ASP / GLU / CYS /
                         LYS / ARG / TYR) resolved according to
                         pKa-prediction + structural context per the
                         existing toolchain. The MD-side calculation
                         consumes whatever protonation state the
                         topology declares; correctness is upstream.
```

---

## Section 3 — Run shape and output cadence

```
Replica shape           1 × 15 ns single replica per protein
                         (chosen 2026-04-30 PM after analysis of
                         extraction-volume cost vs scientific value;
                         per-frame regression model does not require
                         ensemble convergence)
dt                      0.002 ps
nsteps                  7,500,000   (= 15 ns at dt = 0.002)

Output cadence

  nstxout-compressed    10000       (= 20 ps in MD output;
                                      stride 2 in extraction →
                                      40 ps cadence in analysis output;
                                      = 376 frames per 15 ns.
                                      DELIBERATE CHANGE from the prior
                                      production template's value of
                                      5000 (10 ps); 40 ps analysis
                                      cadence is what the per-frame
                                      calibration model is designed to
                                      consume. Picosecond-scale motions
                                      [crankshaft ψᵢ↔φᵢ₊₁, methyl rotation]
                                      are aliased — accepted limitation.)
  nstxout (TRR)         0           (XTC only is sufficient)
  nstvout               0           (no velocities written)
  nstfout               0           (no forces written)
  nstenergy             5000        (= 10 ps energy frames)
  nstcalcenergy         100         (= 0.2 ps energy resolution)
  nstcheckpoint         default     (15 min wallclock)
  nstlog                5000        (= 10 ps progress log)

Solvent in XTC         KEEP full system. Do NOT filter to protein-only
                         at write time. Downstream analyses (hydration
                         shell, water electric field) require all atoms
                         in the trajectory.

Constraints            Existing LINCS defaults from the toolchain.
PME                    Existing defaults from the toolchain.

Per-protein wallclock  ~25 minutes for ~25 K-atom proteins on a 5090
                         (linear extrapolation from the 10-protein
                         calibration fleet baseline at
                         `/shared/2026Thesis/fleet_calibration-backup/*/md.log`:
                         25 ns ran in 29-40 min wallclock at 890-1255
                         ns/day in mid-April, same MDP, same hardware
                         class).
                         For ~200 K-atom proteins (1Z9B-scale), expect
                         ~3-4× longer.
                         If a run substantially exceeds 60 min on a
                         small protein, something is wrong — flag it.

Per-protein output     md.tpr, md.xtc, md.edr, md.cpt, md.gro (plus
                         optional H5MD per Section 5 below).
```

---

## Section 4 — What NOT to change in the GROMACS setup

```
Replica shape           1 × 15 ns. Do NOT switch to 1 × 25 ns, 3 × 5 ns,
                         or any other shape. Deviation breaks
                         downstream-fleet homogeneity and the
                         validation matrix the sample runs anchor.

Cadence                 nstxout-compressed = 10000. Do NOT change to
                         5000, 1000, or use compressed-x-grps with
                         per-group mixed cadence. The per-group mixed-
                         cadence pattern in particular is fraught
                         (silent partial breakage at the index-group
                         level) and is explicitly off for these runs.

Force field             ff14SB / TIP3P. No switch to ff19SB / OPC or
                         any other combination.

energygrps              OFF. Adding `energygrps = Protein SOL Ion`
                         (or finer) forces the short-range non-bonded
                         calculation onto the CPU pathway — GROMACS
                         2026 docs explicitly: "not supported on
                         GPUs". Prohibitive at scale, no useful
                         downstream consumer.

Output filtering        Do NOT filter waters / ions out of the XTC
                         at write time. The full system stays in.
```

---

## Section 5 — H5MD output (optional)

GROMACS 2026.0 supports H5MD as an additive output format. The release notes:

> Trajectory data can now be written to the H5MD file format by gmx
> mdrun. Data written to this file format includes positions, velocities,
> forces, and the simulation box. Chemical bonds are also written to the
> connectivity group. For this release, only lossless output is supported.
> Note that support for this file format in GROMACS tools is experimental.
> Bugs should be expected.

This is **additive** — it does not replace XTC / EDR / TRR. The optional experiment: enable H5MD output on **only one** of the two sample runs (not both) via the appropriate mdrun flag (see `gmx mdrun -h` for the H5MD option in 2026.0). If the H5MD output corrupts or affects standard outputs, the second sample run is unaffected and the experiment can be reported back without losing one of the two sample anchors. If H5MD writes a clean parallel file alongside standard outputs, that establishes the option for the 685-fleet without burning both anchors on an experimental feature.

This experiment is genuinely optional. The sample runs succeed without it.

---

## Section 6 — GROMACS settings to verify in the production MDP template

These are settings that should already be in the existing `prod.mdp` template (they were correct in the prior batch) — listed here as a checklist to verify the upgrade to GROMACS 2026.0 didn't introduce drift. None of these should need changing; if any look wrong, flag it.

```
cutoff-scheme           Verlet
                          Required for modern GPU acceleration. The
                          older "group" scheme was deprecated; if the
                          template still has it, switch to Verlet.

integrator              md
                          Leap-frog Verlet. Standard for biomolecular
                          MD. Don't switch to md-vv unless there's a
                          specific reason.

constraints             h-bonds  (or all-bonds — match the existing
                          template; both work with dt = 0.002 ps)
constraint-algorithm    LINCS
lincs-iter              1
lincs-order             4
                          Standard. SHAKE alternative is unnecessary
                          for unrestrained free MD (Robustelli 2010
                          used SHAKE-on-all-bonds-and-angles for
                          CS-restrained MD; we're not CS-restrained).

tc-grps                 Protein Non-Protein  (or whatever the
                          toolchain canonical is)
ref-t                   298  (or project-standard temperature; match
                          equilibration runs)
tcoupl                  V-rescale  (or stochastic-velocity-rescale)
                          NOT Berendsen for production — Berendsen
                          gives wrong canonical fluctuations.

pcoupl                  Parrinello-Rahman  (or C-rescale)
                          NOT Berendsen for production — same reason.
                          Equilibration can use Berendsen but
                          production must switch.
ref-p                   1.0
tau-p                   2.0  (or template default for the chosen
                          barostat)
pcoupltype              isotropic  (typical for protein-in-water box)
compressibility         4.5e-5  (water at 298 K, 1 bar)

coulombtype             PME
rcoulomb                1.0       (or template default; match equilibration)
pme-order               4
fourier-spacing         0.16      (GROMACS default; matches the active
                                    prep_toolchain template; do NOT
                                    tighten to 0.12 unless the toolchain
                                    canonical changes)
ewald-rtol              1e-5

vdwtype                 Cut-off  (with modifier: Potential-shift,
                                   typically)
rvdw                    1.0       (or template default; match
                                   coulomb cutoff for consistency)

DispCorr                EnerPres
                          Recovers density / pressure under truncated
                          LJ. Standard for biomolecular MD; absent it,
                          density drifts. Worth verifying.

continuation            yes  (assuming production restarts from
                          npt.cpt; if running a fresh production
                          without an equilibration run, no, but
                          equilibration first is the standard path)

gen-vel                 no  (continuation runs); yes only if starting
                          fresh

pbc                     xyz  (cubic or octahedral box per the
                          existing toolchain; not 2D / no-pbc)

comm-mode               linear
nstcomm                 100  (= 0.2 ps)

verlet-buffer-tolerance 0.005  (auto-tuned default; should be
                          appropriate for the chosen rcoulomb / rvdw)
```

**Common pitfalls to watch for during the upgrade to GROMACS 2026:**

- **Equilibration template force-field convention.** Two templates currently exist in active `prep_run_*` directories on disk:
  - **AMBER ff14SB convention (REQUIRED for these runs):** `rcoulomb=1.0`, `rvdw=1.0`, no `vdw-modifier`, `DispCorr=EnerPres`, C-rescale at `tau_p=5.0`. This is what the prep_toolchain canonical produces (e.g. the recent `1Z9B_6577/prep_run_20260430T085928Z/npt.mdp`).
  - **CHARMM36m convention (NOT used here):** `rcoulomb=1.2`, `rvdw=1.2`, `vdw-modifier=force-switch`, `rvdw-switch=1.0`, `DispCorr=No`. Present in some older `prep_run_*` directories.

  These produce different equilibrations. The 685-fleet redispatch must be force-field-homogeneous; the equilibration template must emit the AMBER convention. If CHARMM-convention values appear in any of the equilibration .mdp files (em / nvt / npt / relax), **fail loudly and flag** rather than running through.

- The native `amber14sb_OL15.ff` (or equivalent) ff14SB port in GROMACS 2026 supersedes older translated topologies. Verify that `topol.top` for these sample runs comes from the same source as the calibration-µs scans (1Z9B 2 µs and 1P9J 5 µs); if they're different (e.g. acpype-translated vs native), that's quiet drift across the fleet that should be resolved before the 685 redispatch.
- `tc-grps = "Protein Non-Protein"` is the typical scheme. If the toolchain uses something more granular (e.g., `Protein SOL Ion`), that's fine; just make sure `ref-t`, `tau-t` arrays have matching length.
- Verlet buffer tolerance should match the accuracy expected; the auto-tuned default is usually fine but check the mdrun log for warnings about Verlet pair list scheme issues.
- PME grid changes with box size; if the box is different from prior runs (e.g., due to OF3 vs PDB-extraction structural differences in box-construction), grid auto-tuning will produce different `fourier_nx/ny/nz` values. Acceptable.

---

## Section 7 — GROMACS 2026.0 features to know about

The native AMBER ff14SB port is now in mainline GROMACS (was a separate amber-port previously). If the toolchain has been using a translated / amber-port version, switching to the native is appropriate as part of the upgrade — but the toolchain is the authority on which version's topology files match the prep workflow.

Other 2026.0 features generally NOT relevant to these sample runs:
- OPC / TIP4P-Ew / OPC3 water now native — we're staying with TIP3P (calibration anchor).
- NN/MM electrostatic embedding — out of scope; we're not doing QM/MM.
- FMM long-range electrostatics — we're using PME (standard).
- HIP for AMD GPUs — we're on NVIDIA / CUDA.
- Two new performance metrics (ms/step, Matom·steps/s) reported by mdrun — informative; useful for the wallclock report-back at the end of each run.

The H5MD additive output (Section 5 above) is the one feature this brief explicitly opts into experimenting with.

---

## Section 8 — What this brief explicitly DOES NOT touch

- **Enhanced sampling** of any kind (REMD, AMD, metadynamics, AWH, replica exchange) — unrestrained free MD only.
- **Restraints** of any kind (no NMR distance restraints, no CS restraints, no orientation restraints) — free MD only.
- **Pull / steered MD** — not used.
- **Free-energy calculations** — not used.

Standard relaxation-based protein MD in explicit solvent. Nothing exotic.

---

## Section 9 — Quick reference for execution

1. Confirm the two sample proteins (same two as the first batch — the BMRB-PDB-extracted ones with published S² / T1 / T2 data; the prep session knows which). **Use the same starting structure as the prior 10-protein-fleet calibration set** — the BMRB-PDB-extracted structure used 2026-04-14, NOT a fresh extraction from BMRB. Verify the input PDB hash matches the prior batch's; if they differ, flag before proceeding.
2. Use the existing `prod.mdp` template, with `nsteps = 7,500,000` (= 15 ns at dt = 0.002 ps). Other MDP values per Section 3 above.
3. Confirm ion / pH preparation per the canonical prep_toolchain protocol (Section 2 above defers to the toolchain). Verify equilibration .mdp files emit AMBER convention (Section 6 pitfalls).
4. Run mdrun.
5. After mdrun, the standard output files (md.tpr, md.xtc, md.edr, md.cpt, md.gro) plus optionally an H5MD file are handed off to the downstream extraction session.

**"Well-formed" output verification before handoff:**

```
gmx check -f md.xtc                  exits 0 (no corrupt frames)
gmx energy -f md.edr                  opens; expected per-frame entries
                                      present (Potential, Temperature,
                                      Pressure, Density, etc.)
gmx dump -s md.tpr | head             shows expected nsteps=7500000,
                                      dt=0.002, atom count
md.cpt                                final-step checkpoint; last step
                                      matches nsteps
md.log                                no fatal errors; no Verlet pair-list
                                      warnings; reasonable ns/day metric
```

The prep session reports back with:
- Per-protein wallclock (target ~25 min for ~25 K-atom proteins on a single 5090)
- Confirmation of the four checks above
- Any errors / surprises during prep, equilibration, or production mdrun
- Whether the H5MD option (if attempted on one of the two) worked cleanly
- Whether MPS multi-process scheduling was enabled (affects how the wallclock benchmark interprets)

The downstream extraction / H5 / calculator-output validation is a separate session; the prep session does not need to track those gates.

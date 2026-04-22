# Microsecond-MD Harvester and DFT Queue Design — 2026-04-22

**Status:** DECIDED (design) / NOT YET BUILT (implementation).
Captured 2026-04-22 pm during Session 0 literature-pass close-out.
Becomes au-courant ~2026-04-24 when the current 25 ns fleet run
completes.

**Not for the next session.** The next session is Session 1 drafting
of `spec/PHYSICS_FOUNDATIONS.md` parts 1 and 2. This harvester plan
is queued for after that.

---

## 1 — Context

- Current fleet MD run (25 ns × 10 calibration proteins on
  CHARMM36m) completes ~2026-04-24.
- Three rented RTX 5090 machines available for 41 days, until
  ~2026-06-02 (deadline after which the horses turn into mice —
  any unused compute is zero-valued).
- Batcave (local RTX 5090) is in use for ML work (equivariant
  training, etc.) — stays on ML during the rental window.
- Goal: extend MD on the same 10 calibration proteins from 25 ns
  to ~10 μs each (≈ 75–100 μs total wall time). Direct
  "what does 200× wall time buy us?" comparison against the
  existing 260-DFT calibration + 25 ns extraction pipeline.
- Secondary goal: acquire additional DFT snapshots sampling
  conformational substates not visited at 25 ns. No pre-committed
  DFT budget.

---

## 2 — Options considered

Before landing on the decided design, three framings were weighed.

**Frame A — 10 protein × ~10 μs extension.** Same proteins, longer
trajectories. Direct apples-to-apples against 50 ns pipeline on the
same proteins. Methyl S² converges (H13 μs regime). Scales with the
existing calibration; existing analysis H5 eats it unchanged.
*Chosen.*

**Frame B — 3-protein × ~25 μs ring-flip benchmarking** (BPTI, ubi,
GB3). Cleaner ring-flip rate claims but disconnected from the
260-DFT side. Rejected for thesis coherence.

**Frame C — hybrid 5-protein.** Compromise. Rejected; fewer proteins
at lower sampling gave neither benefit cleanly.

On the DFT side, three sub-options were weighed:

- **Uniform temporal DFT every ~55 ns** — rejected as
  oversamples populated basins and wastes budget on noise.
- **Cluster-centroid DFT** — rejected because requires trajectory
  completion for clustering (posthoc problem); MD and DFT need to
  run concurrently.
- **Daily scanner + AI-curator submission** — *chosen*. Separates
  mechanical pose extraction (trivial, deterministic) from DFT
  submission judgment (contextual, AI-session-driven).

---

## 3 — Decided: machine allocation

- **3 rented 5090s:** run MD on the 10 calibration proteins.
  Roughly 3–4 proteins per machine, rotating as individual runs
  finish. Target 10 μs per protein.
- **Batcave:** ML work during the rental window. After
  2026-06-02, batcave picks up any protein trajectories that
  didn't reach 10 μs on the rented machines.
- **No DFT budget hoarding.** Submit DFTs as pose candidates
  surface; unused budget after 2026-06-02 is zero-valued.

---

## 4 — Decided: pose harvester + AI-curator pattern

Two-layer design separating mechanical extraction from judgment.

### 4.1 — Layer 1: `nmr_extract --harvest-poses` (mechanical scanner)

New mode in `nmr_extract`. Runs daily (cron or user-triggered) over
each protein's trajectory.

Inputs per protein:

- `{protein_dir}/md.xtc` (growing)
- `{protein_dir}/md.tpr`
- `{protein_dir}/.harvest_checkpoint` (last trajectory frame
  scanned; absent on first run)
- `{protein_dir}/.harvest_state.json` (per-residue bin-occupancy
  counters, running-mean backbone structure for RMSD, last χ
  states per residue)

Behaviour:

- Read from checkpoint frame → current end of trajectory.
- Per-residue (φ, ψ) and χ₁ bin occupancy tracked in persistent
  state.
- Extract a pose when any of:
  1. **Novel substate:** a (φ, ψ) or χ₁ bin crosses a population
     threshold (tunable; default ≥ 50 frames) that it had not
     previously crossed. One pose per newly-populated bin per
     scan.
  2. **Large RMSD spike:** backbone RMSD vs running-mean > 3σ
     of running baseline.
  3. **Rotamer transition:** χ angle crosses rotamer boundary
     (±60°, ±180°) on a side-chain-bearing residue.
- Cap extractions at 20 poses per protein per scan run.

Outputs:

- `{protein_dir}/poses/pose_0001.pdb`, `pose_0002.pdb`, ...
  consecutive across all scans. Once written, never renumbered.
- `{protein_dir}/poses/pose_0001.json` metadata sidecar per pose:
  frame index, simulation time, trigger reason
  (`"new_phipsi_bin:R42"`, `"rmsd_spike"`,
  `"chi1_transition:R67"`), bin occupancy snapshot at time of
  extraction, RMSD value, any other cheap context.
- Updated `.harvest_checkpoint` and `.harvest_state.json`.

**No scoring, no prioritisation.** Mechanical only. Judgment is
Layer 2's job.

### 4.2 — Layer 2: AI-curator submission

Periodic Claude session. Walks `poses/` across the 10 proteins.
Reads metadata sidecars. Applies judgment for DFT submission:

- DFT budget remaining (user-tracked, not pre-committed).
- Coverage diversity across the 10 proteins — has protein X been
  sampled on multiple substates yet? Is protein Y over-represented?
- Transition / ring-flip captures of specific interest.
- Anything unexpected in trigger metadata.

Submits selected pose PDBs to the `orca-coord` fleet DFT queue
(existing infrastructure per `project_dft_pipeline_state` memory).
Logs which poses were submitted and why. Unsubmitted poses stay in
place — they can be revisited on subsequent curation rounds.

---

## 5 — Build plan

Engineering estimate: ~5 days total.

- **Day 1–2:** checkpoint mechanism, trajectory frame iteration,
  per-residue (φ, ψ) + χ₁ bin occupancy tracker.
- **Day 3:** trigger logic (novel bin, RMSD spike, χ rotamer
  transition).
- **Day 4:** pose PDB writer + metadata sidecar + checkpoint update.
- **Day 5:** dry-run on existing 25 ns calibration trajectories to
  sanity-check what it harvests. Calibrate thresholds. Verify
  that the scanner picks sensible frames vs the actual 260-DFT
  snapshots.

Builds on existing infrastructure:

- `nmr_extract --trajectory --analysis` machinery
- `orca-coord` DFT queue submission
- Analysis H5 schema (already carries dihedrals, RMSD, etc.)

Adds:

- `--harvest-poses` mode flag
- Scanner module tracking per-protein bin occupancy across runs
- Pose-writer utility (GROMACS / cifpp can handle PDB out)

Nice property: trajectories persist, so if we don't like what's
being harvested mid-run, the fix is a scanner code change — just
re-scan from an earlier checkpoint with new criteria. Low-stakes
iteration.

---

## 6 — Known unknowns (to resolve when this plan is picked up)

1. **Bin granularity.** 30° bins for (φ, ψ) and χ₁? Finer for χ₁
   on aromatic residues (ring-flip-adjacent)? Calibrate at dry-run.
2. **RMSD spike threshold.** 3σ starting point; may need per-
   protein calibration. Running-mean window length also a
   parameter.
3. **Pose cap per scan.** 20 per protein per day → up to 200/day
   across the fleet. If most get submitted (not filtered), 1800
   budget is consumed in 9 days. AI curator is the rate-limit.
4. **Scanner state serialisation format.** JSON for human-readable
   tracking; HDF5 if state grows large.
5. **Scan scheduling.** cron on batcave? watcher process? user-
   triggered? Likely cron for reliability.
6. **MD protocol consistency.** Extension runs must use same FF /
   solvent / ensemble as the 25 ns run. Confirm before launch.
7. **Substate-conditioned kernel-vs-DFT comparison analysis.**
   The downstream thesis analysis that actually uses this data is
   not yet scaffolded. Needs its own design pass when the harvester
   lands.

---

## 7 — Relationship to other specs

- **`spec/IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md` section 5.**
  The `chi_density_grid` and `phi_psi_density` `TrajectoryResult`
  fields are posthoc analyses; this harvester is the online
  companion for DFT acquisition during the trajectory run. The
  downstream substate-conditioned analysis strengthens the
  rollup's thesis-grade synergies (see PHYSICS_FOUNDATIONS 0.9).
- **`spec/PHYSICS_FOUNDATIONS.md` section 0.9.** Synergy #2
  (per-residue CSA tensor validation) and synergy #3 (per-atom-
  to-bulk χ reconstruction) both benefit from substate-resolved
  DFT data. This harvester is the acquisition mechanism for that
  data.

---

## 8 — This document's status

**Tentative. Decided. Not yet built.**

- Design settled in conversation 2026-04-22 pm after Session 0
  literature integration.
- Implementation queued for after ~2026-04-24 (current 25 ns fleet
  run completes).
- Next session is Session 1 drafting of PHYSICS_FOUNDATIONS parts
  1–2 — *not* this.
- Revisions welcome; date any amendments.

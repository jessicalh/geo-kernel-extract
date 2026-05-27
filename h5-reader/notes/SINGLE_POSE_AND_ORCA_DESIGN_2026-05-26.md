# Single-conformation reading + Orca ingestion — design (2026-05-26)

**Living design note.** Captures the architecture decided in the 2026-05-26
design conversation. Authoritative on these topics; where it disagrees with an
older note, this wins. Companion: memory `project_h5reader_formal_model_design`.

## 0. One-line framing

The reader is **basically a trajectory reader** (it stays that). A single pose
is not a new "mode" — it's a single-conformation *sibling* of the trajectory
under a shared base class. Orca shielding is a sparse per-frame overlay that
only some trajectories (and within them, some frames) carry; a standalone Orca
pose is the degenerate one-conformation end of that same spectrum.

## 1. Run shapes — what each extraction mode emits (verified on disk)

`shielding-calcsets/baselines/old-system-baseline/` has one run of each mode:

| Run | Mode | trajectory.h5 | per-atom NPYs | orca_*.npy |
|-----|------|---------------|---------------|------------|
| `trajectory_1p9j_*` | 5 `--trajectory` | **yes** (1.9 GB) | in `per_frame_npys/frame_NNNNNN/` | no |
| `orca-single_*` | 3 `--orca` | **no** | flat in run root | **yes** |
| `orca-mutant_*` | 4 `--mutant` | **no** | flat in run root | **yes** |
| `of3-protonate_*` | 1/2 `--pdb` | **no** | flat in run root | no |

- **Only the trajectory run emits an H5.** It is the dense layer: per-calculator
  time-series `(N, frames[, d])` + `frame_indices` + `frame_times` + Welford
  aggregates + an `/atoms` identity group. Inherently trajectory-shaped.
- **Single-pose runs emit no H5** — the 5 topology-sidecar NPYs +
  `extraction_manifest.json` + the ~90 `WriteAllFeatures` per-atom NPYs flat in
  the run root. Identical payload to one `frame_NNNNNN/` dir; the manifest is
  byte-shape-identical to the trajectory's.
- So a single pose is fully described by **sidecar + manifest → `QtProtein`**
  (loads today, unchanged) **+ one run-root NPY dir → one `QtConformationSnapshot`**.
  Nothing to animate ⇒ no H5 needed.

## 2. The conformation base class (not a dummy frame)

A shared Qt base (QObject-idiomatic — virtual frame-access interface,
`frameChanged`-style signals, proper parent/child ownership) with two concrete
subclasses:

- **H5-trajectory conformation** — the current `QtConformation`; *adds* the
  dense time-series + Welford on top of the base. N frames.
- **Single conformation** — one frame, no time-series; the snapshot *is* the
  run-root NPYs.

"Subclass, don't hack" (`feedback_vtk_subclassing`). The UI / inspector / the
result-group views bind to the base; playback over one frame is trivially short;
the time-series dock engages only for the H5 subclass. The current
`QtConformation` is lifted into base + H5-subclass.

**The only blocker today** is the entry point. `QtProteinLoader::Load(h5_path)`:
Steps 1-2+5 (`QtTopologySidecar::Load` → build `QtProtein`) need **no H5**;
Step 3 hard-builds `QtTrajectoryH5` and Step 6 builds `QtConformation`. Those
two steps are the entire trajectory coupling.

## 3. The dual-purpose snapshot loader (#4)

Write the per-frame loader **directory-agnostic**: it loads a per-atom NPY
directory — a trajectory `frame_NNNNNN/` *or* a single-pose run root — into one
`QtConformationSnapshot`. One mechanism serves both run shapes. A thin
single-conformation entry runs the sidecar→`QtProtein` steps and attaches one
snapshot, skipping the H5/QtConformation/playback layer. (The reader **never
writes H5** — charter — so the single conformation is constructed in memory, not
by synthesizing an H5 file.)

## 4. Orca — three tiers, increasing size

- **(A) Pre-extracted — `QtOrcaGroup` reads `orca_*.npy`.** Mutants and any
  single pose that already has the extracted tensors. `orca_total` /
  `orca_diamagnetic` / `orca_paramagnetic`, each N×9 SphericalTensor. This is a
  normal `#3` mirror result-group (3 `UnpackSphericalTensor` accessors). Only
  ever *populated* on a single-conformation open (never a trajectory frame).
- **(B) Single-conformation open** — §2 + §3. Unlocks actually viewing a DFT
  pose (with its Orca tensors via (A)) and the mutant WT–ALA delta pairs.
- **(C) Raw-`_nmr.out` trajectory overlay** — §5. The genuinely new capability:
  parse ORCA output directly + map onto the trajectory topology. The **3rd
  linked layer** (DFT-pose ↔ frame by stride index).

## 5. The raw-`_nmr.out` path — spec (the hard one; investigated + settled)

**Why raw.** For a trajectory, DFT is run on the per-stride PDBs as a separate
campaign (e.g. `data/trajectories/1p9j-calibration-with-dft/dft/jobs/{job_id}/`,
where `job_id = ..._f000000_t0.0` keys the frame index + time). Those jobs hold
the raw `_nmr.out` (+ `.pdb`/`.xyz`/`meta.json`), **not** extracted `orca_*.npy`,
and the team is **not** re-running the extractor. So the reader parses the
`_nmr.out` itself.

**Verified facts (from a real job + the library + `orca_caller.c`):**
- **Single-point, no optimization** — ORCA input `! r2SCAN def2-SVP def2/J NMR
  CPCM(Water)`, `* Single Point Calculation *`, zero geometry-opt markers, 846
  atoms == the 1P9J trajectory topology. ⇒ DFT geometry == frame geometry, so
  positions are exactly comparable.
- **`orca_caller.c` is copy-only** — counts atoms, copies the `.xyz` verbatim,
  points ORCA at it. No reorder / re-center / H-fiddle.
- **The `_nmr.out` echoes `CARTESIAN COORDINATES (ANGSTROEM)`** alongside the
  shielding section — one file gives tensors + coords + elements + count.

**The matchup verdict.** The original library matchup
(`OrcaShieldingResult.cpp:172-204`) is **already purely positional**: nucleus
`ai` → atom `ai` by index, verified by count-equality + per-atom **element**
equality, hard-failing ("every tensor would go to the wrong atom"). **No
position check** — it can do that because `OrcaRunLoader` builds its Protein
from the pose's own prmtop+xyz (`positions[ai]=xyz[ai]`), so order is correct by
construction.

The reader instead joins the `_nmr.out` to the **trajectory sidecar** (a
possibly different prep — these fleet PDBs came from `orca-fleet/.../pdbs/`). So
the reader **adds a position-to-a-hair + bijection gate** on top of
index+element+count:

> For nucleus *k*: element must match the sidecar atom, AND its echoed
> coordinate must match the frame position within ε ≈ 0.01 Å (cleanly between
> PDB's 1e-3 Å rounding and the ~0.9 Å minimum interatomic spacing). Plus exact
> count and a bijection (every nucleus → exactly one frame atom, none reused).
> Any violation → reject the whole overlay, loud, naming the offending atom.

This is **strictly stricter than the original** (it catches a same-element
reordering the original's element-only check would miss), so "at least as good
as the original model" is met by definition. Two flavors, both pure-positional,
zero protein-SNOBOL (index + coords + element-as-Z only): **(i) by index**
(primary; the original's assumption), **(ii) by coordinate** (spatial-hash
lookup; drop-in if the gate reveals an order difference). A global re-centering,
if ever present, is removed by a one-time Kabsch align (assert RMSD≈0) — still
positional.

**Parse + pack (replicate the library, deterministic, testable):**
- Parse = `OrcaShieldingResult::ParseOrcaNmrOutput`: scan `CHEMICAL SHIELDINGS
  (ppm)`, then per `Nucleus NNNX :` block → `Diamagnetic` 3×3, `Paramagnetic`
  3×3, `Total shielding` 3×3.
- Pack = `SphericalTensor::Decompose(Mat3)` → `PackFull9` (the same convention as
  `OrcaShieldingResult::WriteFeatures`; the reader's `UnpackSphericalTensor`
  reads the inverse layout `[T0, T1[3], T2[5]]`).
- **Regression test = the "at least as good" proof:** the reader's parse+pack
  must *byte-reproduce* the library's `orca_*.npy` on a pose that has both the
  `_nmr.out` and the NPY.

## 6. Sequencing (user-chosen: "finish mirror, then dual-purpose loader")

1. Finish the `#3` mirror: **`QtOrcaGroup` (A)** + the identity groups
   `RingContributions` / `RingGeometry`.
2. **`#4` dual-purpose snapshot loader** (§3) + the conformation base class (§2)
   + the single-conformation entry → unlocks **(B)**.
3. **Raw-`_nmr.out` overlay (C)** — the 3rd linked layer (§5).

## 7. Open items to confirm at build (don't assume)

- A pose with **both** `_nmr.out` and `orca_*.npy` for the §5 regression test
  (orca-single baseline has the NPYs; locate its source `_nmr.out`).
- `pos.npy` / H5-frame position **units** (Å vs nm) for the ε comparison.
- Whether the fleet-PDB atom order actually equals the sidecar order (picks
  flavor (i) vs (ii)); the position gate makes it safe either way.

## Library references (definitive; do not modify)

- `src/OrcaShieldingResult.cpp:66-145` (parser), `:172-204` (the positional
  matchup + element/count gate), `:235-254` (`WriteFeatures` pack).
- `src/OrcaRunLoader.cpp` (the 1:1 `positions[ai]=xyz[ai]` Protein build).
- `data/trajectories/1p9j-calibration-with-dft/dft/jobs/.../*_nmr.out` (real
  output), `dft/campaign/orca_caller.c` (copy-only driver).

# Run descriptor — `h5reader_run.toml` (format spec)

A small TOML file that sits at the head of a run directory and tells the reader
**where to get things from** — so the loader reads a documented map instead of
sniffing. This is the "no file discovery" rule made explicit (CLAUDE.md):
documented conventions only.

> **Status (2026-05-27): FORMAT SPEC ONLY — the descriptor PARSER is not yet
> implemented.** Today the loader opens a run by a bounded convention check
> (`QtProteinLoader::LoadRunPath`: a dir with `trajectory.h5` → load; else a
> single-pose dir; plus a file path straight to the `.h5`), and the strip chart
> locates the DFT campaign with `ReaderMainWindow::locateDftJobsDir` (the dataset
> root holds `extract/` and `dft/` as siblings). This note pins the format we
> will read once a descriptor parser lands (toml++ or a tightly-scoped parser
> with tests — never ad-hoc string splitting). It is written down now because
> the layouts are settled and the demo proved them.

All `[paths]` are relative to the descriptor file unless `relative_to` says
otherwise. Keep them relative so a directory can be moved or copied wholesale.

---

## Why (the ambiguities it removes)

The 1P9J calibration dataset (verified on disk 2026-05-27) is the motivating
case. Its real layout:

```
1p9j-calibration-with-dft/
├── topol.top                      GROMACS topology (parent-of-md convention)
├── md/                            production.{tpr,trr,xtc,edr,...}
├── extract/                       nmr_extract --trajectory output
│   ├── trajectory.h5              dense per-frame H5  (846 atoms, 751 frames)
│   ├── atoms_category_info.npy    ┐
│   ├── bonds.npy                  │ 5 topology sidecars + manifest
│   ├── residues.npy               │ (what QtTopologySidecar::Load expects)
│   ├── ring_membership.npy        │
│   ├── rings.npy                  │
│   ├── extraction_manifest.json   ┘
│   ├── npys/frame_NNNNNN/         per-frame kernel NPYs  (NOT "per_frame_npys")
│   └── pdbs/                      per-frame PDBs
└── dft/jobs/<run>_fNNNNNN_t<ps>/  per-frame ORCA shielding  (*_nmr.out, meta.json)
```

Two real mismatches the descriptor resolves cleanly (both currently handled by
convention, not sniffing): the H5 + sidecars live under `extract/`, not the run
root; and the per-frame dir is `npys/`, not the standard-nmr-extract
`per_frame_npys/`. The DFT campaign is a sibling `dft/`, keyed by the **original
TRR frame index** (0, 2, …, 1500) — the same key as the H5 `frame_indices` and
the `npys/frame_NNNNNN` dirs.

---

## Schema by case

### A. Standard nmr-extract `--trajectory` output

```toml
schema_version = 1
kind           = "h5_trajectory"
protein_id     = "1P9J_5801"

[topology]
dir = "."                       # the 5 sidecar NPYs + extraction_manifest.json

[trajectory]
h5                 = "trajectory.h5"
per_frame_npys_dir = "per_frame_npys"   # standard nmr-extract name (optional; absent ⇒ no per-frame detail)
per_frame_pattern  = "frame_%06d"       # keyed by ORIGINAL frame index

[dft]
enabled = false
```

### B. The 1P9J calibration dataset (H5 under `extract/`, DFT sibling)

```toml
schema_version = 1
kind           = "h5_trajectory"
protein_id     = "1P9J_5801"

[topology]
dir = "extract"                          # sidecars live here, not the root

[trajectory]
h5                 = "extract/trajectory.h5"
per_frame_npys_dir = "extract/npys"      # note: "npys", not "per_frame_npys"
per_frame_pattern  = "frame_%06d"

[dft]
enabled    = true
kind       = "orca"
jobs_dir   = "dft/jobs"
job_pattern = "*_f%06d_t*"               # original frame index, zero-padded to 6
out_primary_from = "meta.json:files.out_primary"   # the SUCCESSFUL .out (never glob *_nmr.out)
frame_key  = "original_index"            # == H5 frame_indices == npys/frame_NNNNNN
# campaign is partial: a missing job dir is an honest GAP, not an error.
```

### C. Single pose (`--orca` / `--mutant` / `--pdb`)

```toml
schema_version = 1
kind           = "single_pose"
protein_id     = "some_pose"

[topology]
dir = "."

[pose]
snapshot_dir = "."                       # flat per-atom NPY payload at the run root
frame_index  = 0
time_ps      = 0.0

[dft]
enabled = false
```

---

## Deliberately NOT in scope

- **A "no-H5 per-frame-dir trajectory" backing.** `nmr_extract --trajectory`
  **always emits a trajectory H5** (CLAUDE.md canonical 5-mode spec), so a
  genuinely H5-less trajectory is not a real production shape — it was only the
  transient state of this dataset mid-extraction (now complete). If a real
  H5-less multi-frame source ever appears, it would be a new `Conformation`
  subclass (`kind = "npy_trajectory"`) keyed off `[frames]` dir/pattern/stride —
  documented here so the schema has room, but not built (no live consumer; the
  project's no-speculative-substrate rule).
- **The descriptor parser itself.** Pending (see status note above).

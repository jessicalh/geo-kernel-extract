# Viewer Quick Start

## Binary location

The viewer binary is `build/nmr-viewer`, built alongside the library:

```
cmake --build build --target nmr-viewer
```

A second build tree at `build-ui/` is supported for standalone UI
work but the canonical launch path uses `build/`. The launch script
`ui/launch_viewer.sh` points at `build/nmr-viewer`.

## Prerequisites

`RuntimeEnvironment::Load()` reads `~/.nmr_tools.toml` at startup.
This file must exist and point to valid paths for mopac, reduce,
OpenBabel, mkdssp, xtb, the ff14SB parameter file, the AIMNet2
model, and (optionally) the tensorcs15 Postgres DSN + the Larsen
H-bond grid directory. If keys are missing the viewer still starts
but the corresponding calculators are silently skipped.

The `[logging]` section of the same TOML configures the
`OperationLog` UDP destination (canonical: `udp_host = "239.255.0.1"`,
`udp_port = 9998` — multicast, so the viewer's Log dock and
`ui/udp_listen.py` can co-receive every datagram).

## Launch script (always use this)

**Always launch through `ui/launch_viewer.sh`** — running
`build/nmr-viewer` directly fails on a cu13 `libnvrtc-builtins.so.13.0`
loader error. The wrapper sets `LD_LIBRARY_PATH` to the bundled cu13
runtime via `scripts/run_with_cuda_env.sh`. Env-var assignments from
inside the process don't reach the dynamic loader once dlopen has
cached its resolution path.

```
bash ui/launch_viewer.sh                              # default: 1ubq
bash ui/launch_viewer.sh --pdb myprotein.pdb          # bare PDB
bash ui/launch_viewer.sh --protonated-pdb pose.pdb    # skip reduce
bash ui/launch_viewer.sh --orca --root path/to/root   # ORCA DFT
bash ui/launch_viewer.sh --pdb pose.pdb --analysis-h5 traj.h5
```

## Command-line modes

The viewer uses the same `JobSpec` parser as `nmr_extract`. Four
interactive modes are available (Trajectory mode is CLI-only; use
`nmr_extract --trajectory` for ensemble extraction):

```
# Bare PDB — reduce protonates, ff14SB charges assigned.
# Reduce may re-evaluate HIS/LYS protonation; if loading an MD
# snapshot prefer --protonated-pdb to preserve the existing state.
nmr-viewer --pdb path/to/protein.pdb

# Pre-protonated PDB (from reduce, tleap, MD output, etc.).
# Trusts the existing protonation; skips reduce.
nmr-viewer --protonated-pdb path/to/protein.pdb

# ORCA DFT run — expects {root}.xyz, {root}.prmtop, optional {root}_nmr.out
nmr-viewer --orca --root path/to/A0A7C5FAR6_WT

# WT + ALA mutant pair — each root expands the same way
nmr-viewer --mutant --wt path/to/WT --ala path/to/ALA
```

### Trajectory companion (read-only)

`--analysis-h5 PATH` pairs the loaded protein with a `trajectory.h5`
companion (one HDF5 group per attached `TrajectoryResult` under
`/trajectory/<name>/`). The viewer reads atom identity, frame
metadata, Welford rollups, and frame-0 slabs from the file and
surfaces them in the **Time Series (H5)** dock tab per picked atom.
Sparse-set tolerant — only the groups present in the file render.

ComputeWorker performs an atom-count + per-index Element enum
identity check at load; if the H5 describes a different protein the
binding is refused and the Time Series tab shows "no H5 loaded".
The viewer never writes H5 or triggers extraction.

### Other flags

- `--config path/to/params.toml` — override calculator parameters
  (see `CALCULATOR_PARAMETER_API.md`).
- `--output DIR` — write feature NPY arrays to DIR on completion;
  optional in the viewer (export is also available via File menu and
  the REST `export_features` command).
- `--rest-port PORT` — REST server port, default 9147; `0` disables.

## REST interface

The viewer listens on TCP port 9147 for newline-delimited JSON
commands. Common ones:

```
echo '{"cmd":"status"}'                          | nc -q1 localhost 9147
echo '{"cmd":"atom_dump","atom":42}'             | nc -q1 localhost 9147
echo '{"cmd":"list_atoms"}'                      | nc -q1 localhost 9147
echo '{"cmd":"screenshot","path":"/tmp/out.png"}'| nc -q1 localhost 9147
echo '{"cmd":"reset_view"}'                      | nc -q1 localhost 9147
echo '{"cmd":"export_features","path":"/tmp/f"}' | nc -q1 localhost 9147
echo '{"cmd":"quit"}'                            | nc -q1 localhost 9147
```

`{"cmd":"quit"}` is preferred over `pkill nmr-viewer` — Qt drains
the write buffer, fires `aboutToQuit`, and `MainWindow::shutdown`
stops timers / workers / VTK before the process exits.

`atom_dump` returns the full typed inspector tree for one atom; when
`--analysis-h5` was supplied and the identity check passed, the
response includes an `h5` block carrying per-TR-group Welford
rollups + frame-0 values for that atom. The companion test client
`ui/tests/test_inspector.py` drives `atom_dump` over every atom to
report structural invariants (PASS / NOTE / WARN; never asserts,
per `feedback_log_overages_dont_assert`).

Full command list: `ui/src/REST_INTERFACE_SPEC.md`.

## Interaction

- **Rotate/zoom/pan**: standard VTK mouse controls (left-drag
  rotates, right-drag zooms, middle-drag pans).
- **Atom inspector**: double-click an atom to open the inspector
  dock. Shows identity, AMBER substrate (planar group, ring position,
  pseudoatom, prochiral stereo), charges, every per-calculator
  shielding contribution the live pipeline produces (BS, HM,
  McConnell, Pi-quadrupole, RingChi, Dispersion, HBond kernel-form +
  Larsen H-bond grid, Tripeptide BB + Neighbor, AIMNet2 EFG), ring
  neighbours with per-ring G and H decompositions, bond neighbours,
  vector fields, DSSP, ORCA DFT, MOPAC slice when present, AIMNet2
  polarisability + embedding, MolecularGraphResult graph topology.
- **Bond inspector**: shift+double-click → bond identity, MOPAC
  Wiberg order, endpoint charges, McConnell contribution.
- **Time Series (H5) tab**: per-picked-atom view of the trajectory
  companion. Rollup (mean ± std) alongside frame-0 slab, one
  section per TR group present in the file.
- **Sidebar**: per-calculator toggles, overlay controls, opacity
  sliders.
- **File > Export Features**: writes NPY arrays to a chosen directory.

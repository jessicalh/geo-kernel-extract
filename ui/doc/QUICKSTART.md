# Viewer Quick Start

## Binary location

The viewer binary is `nmr-viewer`. There are two build trees:

- `build/nmr-viewer` — built alongside the library (`cmake --build build --target nmr-viewer`)
- `build-ui/nmr-viewer` — standalone UI build

Both produce the same binary. The launch script uses `build-ui/`.

## Prerequisites

`RuntimeEnvironment::Load()` reads `~/.nmr_tools.toml` at startup.
This file must exist and point to valid paths for mopac, reduce,
OpenBabel, mkdssp, xtb, and the ff14SB parameter file. If it is
missing or incomplete the viewer will still start but some
calculators will fail silently.

## Command line

The viewer uses the same JobSpec parser as `nmr_extract`. All five
modes are available:

```
# Bare PDB — reduce protonates, ff14SB charges assigned
nmr-viewer --pdb path/to/protein.pdb

# Pre-protonated PDB (from reduce, tleap, GROMACS, etc.)
nmr-viewer --protonated-pdb path/to/protein.pdb

# ORCA DFT run — expects {root}.xyz, {root}.prmtop, optional {root}_nmr.out
nmr-viewer --orca --root path/to/A0A7C5FAR6_WT

# WT + ALA mutant pair — each root expands the same way
nmr-viewer --mutant --wt path/to/WT --ala path/to/ALA

# GROMACS ensemble — TPR topology + pose directory
nmr-viewer --fleet --tpr path/to/prod.tpr --poses path/to/poses
```

All modes accept `--config path/to/params.toml` to override
calculator parameters (see `CALCULATOR_PARAMETER_API.md`).

The `--output DIR` flag is accepted but optional in the viewer.
When given, feature arrays are written to that directory on
completion. Without it the viewer computes and displays but does
not write to disk. Use File > Export Features or the REST command
for on-demand export.

## Launch script

The script `ui/launch_viewer.sh` is the easiest way to start:

```
bash ui/launch_viewer.sh                              # 1ubq protonated (default)
bash ui/launch_viewer.sh --pdb myprotein.pdb          # bare PDB
bash ui/launch_viewer.sh --orca --root path/to/root   # ORCA DFT
```

It starts a UDP log listener (`ui/udp_listen.py`) in the
background, sets `DISPLAY=:1`, and passes all arguments through
to `nmr-viewer`. The log listener prints structured JSON log
messages to the terminal as the pipeline runs.

## REST interface

The viewer listens on TCP port 9147 for JSON commands. Send
with netcat:

```
echo '{"cmd":"status"}' | nc -q1 localhost 9147
echo '{"cmd":"screenshot","path":"/tmp/shot.png"}' | nc -q1 localhost 9147
echo '{"cmd":"reset_view"}' | nc -q1 localhost 9147
echo '{"cmd":"export_features","path":"/tmp/features"}' | nc -q1 localhost 9147
echo '{"cmd":"load_pdb","path":"/path/to/file.pdb"}' | nc -q1 localhost 9147
```

The full command list is in `ui/src/REST_INTERFACE_SPEC.md`.

## Interaction

- **Rotate/zoom/pan**: standard VTK mouse controls (left-drag
  rotates, right-drag zooms, middle-drag pans).
- **Atom inspector**: double-click an atom to open the inspector
  dock. Shows identity, charges, all 8 calculator shielding
  tensors, ring neighbours with per-ring G and H decompositions,
  bond neighbours, vector fields, DSSP, and ORCA DFT when present.
- **Sidebar**: per-calculator toggles, overlay controls, opacity
  sliders.
- **File > Export Features**: writes NPY arrays to a chosen
  directory.

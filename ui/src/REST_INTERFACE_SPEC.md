# Protein Tensor Viewer — REST Interface

## Transport

- **QTcpServer** on `localhost:9147` (tries 9147–9156 if port in use)
- Newline-delimited JSON: one object per line, one response per line
- Synchronous within the Qt event loop
- Starts automatically (disable with `--rest-port 0`)

## Response format

```json
{"ok": true, "result": { ... }}
{"ok": false, "error": "description"}
```

## Commands

### Session

#### `status`
```json
{"cmd": "status"}
```
Returns: protein name, atom/ring/residue/tier counts, compute state,
grid counts and T0 range, overlay visibility.

#### `load_pdb`
```json
{"cmd": "load_pdb", "path": "/abs/path/to/protein.pdb"}
```
Builds a JobSpec with mode=Pdb and triggers compute. Returns immediately
with `{"status": "loading"}` — poll `status` for completion.

#### `load_protein_dir`
**Deprecated.** Use `--orca --root` or `--mutant --wt/--ala` from the
command line instead. Retained for backwards compatibility — treats the
path as a bare PDB load.

#### `export_features`
```json
{"cmd": "export_features", "path": "/abs/path/to/output_dir"}
```
Writes all computed feature arrays (NPY) to the directory. Fleet mode
creates `frame_N/` subdirectories automatically.
Returns: `{"arrays": 44, "path": "...", "frames": N}` (frames only for fleet).

### Camera

#### `get_camera`
```json
{"cmd": "get_camera"}
```
Returns: position, focal_point, view_up (each `[x,y,z]`), distance, view_angle.

#### `set_camera`
```json
{"cmd": "set_camera", "position": [x,y,z], "focal_point": [x,y,z], "view_up": [0,1,0]}
```
All fields optional — omitted fields keep their current value.

#### `orbit`
```json
{"cmd": "orbit", "azimuth": 30, "elevation": 0}
```
Rotates camera around focal point. Angles in degrees.

#### `reset_view`
```json
{"cmd": "reset_view"}
```

#### `look_at_ring`
```json
{"cmd": "look_at_ring", "ring": 0, "view": "side", "distance": 15.0}
```
Frames camera on an aromatic ring by index. Views: `"side"` (sees butterfly
lobes), `"top"` (looks down normal), `"edge"`.
Returns: ring_index, ring_type, center, normal, view.

#### `look_at_atom`
```json
{"cmd": "look_at_atom", "atom": 42, "distance": 15.0}
```
Re-centers camera on atom, keeping current view direction.
Returns: atom_index, element, pdb_name, residue.

### Rendering

#### `set_render_mode`
```json
{"cmd": "set_render_mode", "mode": "ball_stick"}
```
Modes: `ball_stick`, `stick` (alias `liquorice`).

#### `set_glyph_scale`
```json
{"cmd": "set_glyph_scale", "scale": 0.5}
```

#### `set_opacity`
```json
{"cmd": "set_opacity", "value": 0.7}
```

### Overlays

#### `show_rings`
```json
{"cmd": "show_rings", "visible": true}
```

#### `show_bonds`
```json
{"cmd": "show_bonds", "visible": true}
```
Toggles peptide bond overlay.

#### `show_butterfly`
```json
{"cmd": "show_butterfly", "visible": true}
```
Toggles Biot-Savart magnetic field visualization.

#### `show_field_grid`
```json
{"cmd": "show_field_grid", "shielded": true, "deshielded": true}
```
Toggle shielded/deshielded isosurface regions independently, or use
`"visible": true` to toggle both.

#### `set_iso_threshold`
```json
{"cmd": "set_iso_threshold", "value": 0.5}
```
Isosurface contour level (0–1 maps to slider range).

#### `set_overlay`
```json
{"cmd": "set_overlay", "mode": "classical"}
```
Legacy — overlay modes have been removed. Returns a note directing
to per-calculator toggles.

#### `set_calculators`
```json
{"cmd": "set_calculators", "bs": true, "hm": false}
```
Legacy — per-calculator visualization toggles pending.

### Data

#### `list_rings`
```json
{"cmd": "list_rings"}
```
Returns array of ring objects: index, type, residue, center, normal,
radius, intensity.

#### `get_log`
```json
{"cmd": "get_log", "lines": 50}
```
Returns last N lines from the operations log panel. Alternative forms:
`{"first": 0, "last": 99}` for a specific range, or omit all for
the full log.

### Output

#### `screenshot`
```json
{"cmd": "screenshot", "path": "/abs/path/to/image.png"}
```
Saves current render. Returns: path, width, height.

### Introspection

These endpoints return per-atom and per-residue typed data as JSON,
suitable for test clients and scripted exploration. They surface the
same typed object model the Atom Inspector dock shows. No
test-specific shortcuts — adding a new section to `populateAtomInfo`
gets the same data into `atom_dump` for parity (the
Inspector / atom_dump / test_inspector triple is documented in
`ui/CLAUDE.md`).

#### `atom_dump`
```json
{"cmd": "atom_dump", "atom": 42}
```
Returns the full typed inspector tree for one atom: identity
(element, role, residue, terminal state, protonation variant),
AMBER substrate (planar group, polar H kind, prochiral stereo, ring
position, pseudoatom kind), charges, per-calculator shielding
contributions, ring neighbours, bond neighbours (covalent + McConnell
+ MOPAC), DSSP, ORCA DFT slice, and — when `--analysis-h5` was
supplied and the identity check passed — an `h5` block carrying the
trajectory companion data for this atom (Welford rollups + frame-0
slabs per TR group present in the file).

The `h5` block is absent (not null) when no companion H5 is loaded
or the identity check failed; consumers gate on `.get("h5")`.

#### `list_atoms`
```json
{"cmd": "list_atoms"}
{"cmd": "list_atoms", "filter": {"role": "AmideH", "residue_type": "HIS"}}
```
Returns a concise per-atom record (index, element, name, residue
type/index, role, ring/substrate flags). Optional `filter` narrows
by `role`, `element`, `residue_type`, `residue_index`, `in_ring`,
`planar_group`, `polar_h`, `is_amide_H`, `is_methyl`.

### Lifecycle

#### `quit`
```json
{"cmd": "quit"}
```
Polite shutdown. Returns `{"ok": true, "result": {"message": "shutting down"}}`
synchronously; the reply ships from the dispatch thread, then
`socket->disconnectFromHost` drains the write buffer before Qt
emits `disconnected`, at which point the server posts
`QCoreApplication::quit` via `QueuedConnection`. The `aboutToQuit`
signal triggers `MainWindow::shutdown` to stop timers / workers /
VTK in order. Use in preference to `pkill nmr-viewer` so widget
destructors and the worker-thread auto-delete chain run cleanly.

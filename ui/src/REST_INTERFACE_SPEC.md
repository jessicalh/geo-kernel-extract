# bs-viewer REST Interface Specification

## Purpose

Give external tools (Claude Code, scripts, tests) programmatic control of the
Qt/VTK molecular viewer. Modelled after ChimeraX MCP: load, configure,
screenshot, iterate. All existing visual features are exposed with full
parameterization — nothing is simplified for AI consumption.

## Transport

- **QTcpServer** on `localhost:9147` (port mnemonic: NMR)
- Newline-delimited JSON: one JSON object per line, one JSON response per line
- Synchronous request/response within the Qt event loop
- Server starts automatically with the viewer (opt-out via `--no-rest`)

## Response format

Every response:
```json
{"ok": true, "result": { ... }}
{"ok": false, "error": "description"}
```

## Commands

### Session

#### `load_pdb`
Load a PDB file and run full physics computation.
```json
{"cmd": "load_pdb", "path": "/abs/path/to/protein.pdb"}
```
Response includes atom count, ring count, protein name. Blocks until compute
finishes (may take seconds for large proteins).

#### `load_protein_dir`
Load a prepared protein directory (PDB + charges + APBS field).
```json
{"cmd": "load_protein_dir", "path": "/abs/path/to/prepared/P84477"}
```

#### `status`
```json
{"cmd": "status"}
```
Returns: protein name, atom count, ring count, residue count, compute state,
current overlay mode, current viz mode, which calculators are enabled.

---

### Camera

#### `set_camera`
Explicit camera placement.
```json
{"cmd": "set_camera", "position": [x,y,z], "focal_point": [x,y,z], "view_up": [0,1,0]}
```

#### `look_at_ring`
Frame camera on a specific aromatic ring. Positions camera along ring normal
at a comfortable distance, looking at ring center.
```json
{"cmd": "look_at_ring", "residue_number": 42, "ring_type": "PHE"}
```

#### `look_at_residue`
Frame on a residue's CA atom.
```json
{"cmd": "look_at_residue", "residue_number": 42}
```

#### `reset_view`
Reset camera to show full molecule.
```json
{"cmd": "reset_view"}
```

#### `orbit`
Rotate camera around the current focal point. Essential for verifying 3D
scenes from multiple angles — compensates for 2D screenshot limitation.
```json
{"cmd": "orbit", "azimuth": 30, "elevation": 0}
```
Angles in degrees. Azimuth rotates horizontally, elevation rotates vertically.

#### `turntable`
Take N screenshots at evenly spaced azimuth angles (full 360° orbit).
Saves to directory with numbered filenames.
```json
{"cmd": "turntable", "n_frames": 8, "output_dir": "/path/to/frames", "width": 1920, "height": 1080}
```
Returns list of saved file paths.

---

### Molecule Rendering

#### `set_render_mode`
```json
{"cmd": "set_render_mode", "mode": "ball_stick"}
```
Modes: `ball_stick`, `stick`, `vdw`, `backbone`

#### `color_atoms`
Color specific atoms or groups. Does not affect overlay colors.
```json
{"cmd": "color_atoms", "selection": "residue", "residue_number": 42, "color": [1.0, 0.5, 0.0]}
```
Selection types:
- `"all"` — every atom
- `"element"` + `"element": "C"` — by element
- `"residue"` + `"residue_number": 42` — by residue
- `"chain"` + `"chain_id": "A"` — by chain
- `"backbone"` / `"sidechain"` — structural role
- `"aromatic"` — atoms in aromatic rings
- `"range"` + `"start": 10, "end": 50` — residue range

Color is RGB float [0-1]. Optional `"opacity"` field (0-1).

#### `color_by_attribute`
Color all atoms by a computed attribute using the diverging colormap.
```json
{"cmd": "color_by_attribute", "attribute": "total_T0"}
```
Attributes: `total_T0`, `bs_T0`, `hm_T0`, `dft_T0`, `ml_T0`, `tier`,
`secondary_structure`, `bfactor`, `element`

---

### Overlay Control

#### `set_overlay`
Switch the active overlay mode.
```json
{"cmd": "set_overlay", "mode": "isosurface"}
```
Modes: `none`, `rings`, `bonds`, `tensor_glyph`, `ellipsoid`, `scalar_field`,
`arrows`, `isosurface`, `butterfly`, `comparison`

#### `set_glyph_style`
```json
{"cmd": "set_glyph_style", "style": "spherical_harmonic"}
```
Styles: `spherical_harmonic`, `ellipsoid`

#### `set_glyph_scale`
```json
{"cmd": "set_glyph_scale", "scale": 0.5}
```

#### `set_opacity`
```json
{"cmd": "set_opacity", "value": 0.7}
```

#### `set_iso_threshold`
Isosurface contour level.
```json
{"cmd": "set_iso_threshold", "value": 0.5}
```

#### `set_gaussian_radius`
Isosurface smoothing radius (Angstroms).
```json
{"cmd": "set_gaussian_radius", "value": 2.0}
```

#### `set_current_scale`
Ring current arrow/field scale factor.
```json
{"cmd": "set_current_scale", "value": 1.0}
```

#### `show_rings`
Toggle aromatic ring outlines.
```json
{"cmd": "show_rings", "visible": true}
```

#### `show_peptide_bonds`
```json
{"cmd": "show_peptide_bonds", "visible": true}
```

#### `show_butterfly`
Toggle Biot-Savart magnetic field streamlines.
```json
{"cmd": "show_butterfly", "visible": true}
```

---

### Physics Calculator Selection

#### `set_calculators`
Enable/disable individual physics contributions. Affects which calculators
contribute to total T0, tensor glyphs, isosurfaces, and scalar field coloring.
```json
{"cmd": "set_calculators", "bs": true, "hm": true, "mc": false, "ld": false, "ce": true, "pq": false, "rsa": false, "hb": false}
```

---

### Visualization Mode

#### `set_viz_mode`
What scalar value drives the overlay (color, isosurface, field).
```json
{"cmd": "set_viz_mode", "mode": "total_T0"}
```
Modes: `total_T0`, `per_calc`, `dft_delta`, `ml_pred`, `tier`

---

### Output

#### `screenshot`
Save current render to file.
```json
{"cmd": "screenshot", "path": "/abs/path/to/image.png", "width": 1920, "height": 1080}
```
Returns: `{"ok": true, "result": {"path": "/abs/path/to/image.png", "width": 1920, "height": 1080}}`

#### `get_atom_data`
Return per-atom computed data as JSON.
```json
{"cmd": "get_atom_data", "residue_number": 42}
```
Without residue_number, returns all atoms. Response includes per-calculator T0,
spherical tensor components, B-field, E-field, tier, DSSP, ML prediction.

#### `get_ring_data`
Return aromatic ring data.
```json
{"cmd": "get_ring_data"}
```
Response: center, normal, radius, intensity, residue, ring type for each ring.

#### `get_overlay_info`
Return current overlay state — what's visible, what parameters are set.
```json
{"cmd": "get_overlay_info"}
```

---

### Exploration Protocol

#### `cycle_overlays`
Step through overlay modes, taking a screenshot at each. For systematic
visual evaluation of all representations on the current protein.
```json
{"cmd": "cycle_overlays", "output_dir": "/path/to/cycle", "width": 1920, "height": 1080, "modes": ["none", "isosurface", "tensor_glyph", "scalar_field", "arrows", "butterfly"]}
```
Returns list of `{mode, path}` pairs.

#### `sweep_parameter`
Vary a single parameter across a range, taking screenshots. For finding
the right isosurface threshold, glyph scale, etc.
```json
{"cmd": "sweep_parameter", "parameter": "iso_threshold", "values": [0.1, 0.3, 0.5, 0.8, 1.0, 2.0], "output_dir": "/path/to/sweep", "width": 1920, "height": 1080}
```

---

## Implementation Notes

### Architecture

Add `RestServer` class (QTcpServer subclass) owned by MainWindow.
- Listens on port 9147
- Reads newline-delimited JSON from client socket
- Dispatches to MainWindow slots via queued invocations (thread-safe)
- Serializes responses back to client

### Thread Safety

All VTK operations must happen on the main thread. The RestServer's
`readyRead` handler should use `QMetaObject::invokeMethod(mainWindow,
..., Qt::QueuedConnection)` for anything that touches the renderer.
ComputeWorker already runs on a QThread with signal/slot communication.

### Screenshot Implementation

Use `vtkWindowToImageFilter` → `vtkPNGWriter` on the render window.
Already partially implemented in `MainWindow::saveScreenshot()` — factor
out the core logic and add width/height control via `SetSize()`.

### Startup

```
bs-viewer --rest-port 9147 --pdb data/prepared/P84477/protonated.pdb
```
Or in headless mode for CI:
```
bs-viewer --rest-port 9147 --offscreen --pdb ...
```
(Offscreen rendering requires VTK OSMesa build — note for future.)

### Dependencies

Only Qt6::Network (QTcpServer, QTcpSocket) on top of existing deps.
No external HTTP library needed — raw TCP with newline-delimited JSON
is simpler and sufficient for localhost control.

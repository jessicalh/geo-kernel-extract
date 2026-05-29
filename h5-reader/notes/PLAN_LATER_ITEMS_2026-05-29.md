# Plan: LATER items from Phase A-H cleanup (2026-05-29 handoff)

## Context

Phases A-H of the 5-TR surfacing landed (2026-05-29) plus the 4 NOW
items from adversarial review. The 10 LATER items below are polish,
defensive hardening, data extensions, and one new rendering surface
(the Mat3 ellipsoid glyph) that were deliberately deferred after the
"wrap through phases completely before going into long fix mode"
directive.

**Inputs to load at session start:**

- `h5-reader/notes/SCOPE_NEW_TRS_2026-05-29.md` — the original scope doc
- `h5-reader/notes/CLEANUP_AFTER_PHASE_AH_2026-05-29.md` — the punch list these items come from
- `/home/jessica/.claude/plans/validated-questing-metcalfe.md` — the Phase A-H plan that produced the current state
- The 5 producer-side TRs at `src/{DihedralAutocorrelation,IRedOrderParameter,KernelCoherence,KernelDynamics,ReorientationalDynamics}TrajectoryResult.{h,cpp}` (read-only — never modified)

**Baseline before starting:** 3/3 ctest targets green
(`h5reader_model_tests`, `h5reader_scene_math_tests`,
`h5reader_rest_smoke`). The new code is uncommitted; current `git
status` shows 15 modified + 12 untracked files under `h5-reader/`.

**Decision: commit or rebase before starting?** The Phase A-H + NOW
work is a single coherent feature landing. Recommend either:
- (a) Commit all current work as one "Surface 5 new TRs end-to-end"
  commit before starting this LATER plan, OR
- (b) Treat this LATER plan as continuation in the same uncommitted
  worktree and commit the whole thing together at the end.

The CLEANUP_AFTER doc lists every file; either path is clean.

## Critical files (by phase)

### Phase L-1: Low-risk housekeeping (5 items, ~30 lines total)

These are "land them as one commit, move on" items. No architectural
risk; each is a few lines.

- **L-1a (LATER-8): Drop `mutable` on `ownedPanels_`.**
  `h5-reader/src/app/DashboardDisplayController.h:~183`. The `mutable`
  qualifier is unused — only `takeOwnedPanels()` mutates it, and that
  method is non-const. Drop the keyword. Build, done.

- **L-1b (LATER-10): Status text for owned-panel signals.**
  `h5-reader/src/app/DashboardDisplayController.cpp` —
  `updateStatusText()` (~line 1380). Currently reports "No active
  strip display modes" when only owned-panel signals are active. Add
  a parallel `activeOwnedPanelCount_` member, increment it in the
  rebuild loop's panel-mode branch, and surface it in the status text:
  `"%1 strip + %2 panel signal(s)"`.

- **L-1c (LATER-3): REST `anchorToJson` BondVector arm.**
  `h5-reader/src/app/RestServer.cpp:~83`. The `anchorToJson()` helper
  has no `BondVectorAnchor` case — `/dashboard/signals` returns `{}`
  for the anchor field on any iRED/Reorient signal. Add a
  `bond_vector` JSON kind with `residue` and `kind` int fields,
  modelled on the existing `atom_tuple` arm. Then add a small REST
  test in `tests/rest/` that asserts the new shape.

- **L-1d (LATER-1): Extract `lookupBondVector` helper.**
  `h5-reader/src/app/SceneRevealOverlay.cpp` —
  `atomsForBinding` and `reveal()` both have a duplicated walk-iRED-
  then-walk-Reorient block. Pull it into a free function
  `lookupBondVector(h5, residue, kind) → optional<{tail, head}>` at
  the top of the file. Both call sites consume it. Saves ~20 lines.
  No behaviour change; reviewers will appreciate.

- **L-1e (LATER-2): HighFive narrowing + exception safety on the readers.**
  `h5-reader/src/io/QtTrajectoryH5.cpp` — every `Read*` reader
  function (ReadIRedOrderParameters, ReadKernelDynamics, ReadKernel-
  Coherence, ReadReorientationalDynamics, ReadDihedralAutocorrelation)
  performs unchecked HighFive reads. If the H5 dtype doesn't match
  the target vector type, HighFive throws and the whole trajectory
  load aborts. Wrap each reader body in `try { ... } catch
  (const HighFive::Exception& e) { WarnShapeMismatch(...); out.reset();
  }`. Same pattern for `std::exception`. This is a uniform safety net
  matching the existing `WarnGroupAbsent` / `WarnShapeMismatch` pattern.

  Also: the int32 identity columns currently match the producer's
  emission dtype. Add an explicit `getDimensions()` + dtype assertion
  for the integer arrays (residue_index, n_atom, h_atom, vector_kind,
  owning_atom, tail_atom, head_atom) so a future producer that emits
  int64 doesn't silently narrow.

### Phase L-2: Data extensions + chassis cleanup (~150 lines)

- **L-2a (LATER-6): Chi[0..3] dihedral autocorrelation.**
  Phase F landed phi/psi only. The producer emits `chi_acf` (R, 4, L)
  and `chi_corr_time_ps` (R, 4) for the 4 chi torsions per residue.
  Two design options:
  - Option A (8 separate descriptors): h5:dihedral_chi0_corr_time …
    chi3_corr_time + chi0_acf … chi3_acf. Clean, follows phi/psi
    pattern exactly. Total: 8 catalog entries + 8 dispatch branches.
  - Option B (2 composite descriptors with PerClassBlock + 4 channels):
    h5:dihedral_chi_corr_time and h5:dihedral_chi_acf, channels =
    {chi0, chi1, chi2, chi3}. Less catalog noise; per-channel dispatch
    inside the controller's denseH5Plan branch + panel builder.
  - **Recommend Option B.** Existing dihedral_time_series in the
    catalog uses PerClassBlock with phi/psi/omega/chi channels — same
    pattern. Files: `TrajectorySignalCatalog.cpp`,
    `DashboardDisplayController.cpp` denseH5Plan branch + panel
    builders, `QtPerResidueBuffers.h` extension for per-chi data,
    `QtTrajectoryH5.cpp` reader extension, `dashboard_model_tests.cpp`
    catalog presence row.

- **L-2b (LATER-7): Unify the three-bucket panel render order.**
  `h5-reader/src/app/StripStackWidget.cpp` currently paints in three
  sequential blocks (tracks_, spectrumTracks_, ownedPanels_). Phase B
  shortcut. For heterogeneous interleaving (a researcher wants
  SequenceBarPanel above a TemporalStripPanel for the same residue),
  the order needs to be controller-determined.

  Approach: refactor StripStackWidget to a single ordered
  `std::vector<std::unique_ptr<AbstractStripPanel>> panels_` member.
  setTracks/setSpectrumTracks become builders that produce
  `TemporalStripPanel`/`SpectrumStripPanel` instances and push them
  into panels_; setOwnedPanels appends owned panels. Need to move
  TemporalStripPanel + SpectrumStripPanel to take Track / SpectrumTrack
  by value (currently hold `const Track&` — invalid after the unified
  rebuild path replaces tracks_ but not the panel object).

  Files: `StripStackWidget.{h,cpp}`. Risk: the existing
  Track/SpectrumTrack API stays as setter input; only the internal
  storage shape changes. Tests should pick up any regression in
  `h5reader_rest_smoke`'s dashboard-path test.

### Phase L-3: New rendering surfaces (~250 lines)

These are the genuinely-new visual work. Each is a self-contained
panel addition.

- **L-3a (LATER-4): Mat3 tensor-glyph rendering.**
  Reorient's `bond_orientation_tensor` (per-vector 3x3) is loaded
  into the buffer but only rendered as `static.table` placeholder.
  Add a VTK ellipsoid actor in `SceneRevealOverlay` that, when a
  Mat3PerRow descriptor is the active reveal, paints a translucent
  ellipsoid attached to the bond. Eigendecomposition of the 3x3
  matrix gives the ellipsoid's principal axes and lengths. Use
  `vtkSphereSource` + `vtkTransform` to scale + rotate.

  Files: `SceneRevealOverlay.{h,cpp}` (new actor + lifecycle), maybe
  a small helper file `TensorGlyphMath.{h,cpp}` if the
  eigendecomposition warrants its own unit test (use Eigen). Phase B
  established the pattern of pulling math into a separately-testable
  header — follow that.

  Trigger: pick descriptor `h5:reorient_orientation_tensor` and add
  a `static.tensor.glyph` display mode to the catalog (extend
  `isPanelMode` in `DashboardDisplayController.cpp`, extend
  `SignalDisplayDialog.cpp`'s DisplayModeKind with a `TensorGlyph`
  arm — wait, that already exists). Actually `static.tensor` already
  maps to TensorGlyph kind in the dialog; just ensure the Reorient
  orientation descriptor has it in its `staticModes` list.

- **L-3b (LATER-5): FixedFreqBlock J(ω) panel.**
  Reorient's `spectral_density_j` is (V, 5) — J at the 5 KTB Larmor
  combination frequencies. Currently no dedicated panel. Build a
  small scatter+line plot panel that shows J(ω) at the 5 freqs as
  log-x markers connected by polylines, one trace per selected bond
  vector.

  Reuses PowerSpectrumPanel's painter machinery but with 5 discrete
  x positions instead of a dense F grid. Could either extend
  PowerSpectrumPanel to optionally render as scatter-with-markers
  when n_samples is small, OR add a new `FixedFreqPanel` subclass.
  **Recommend the latter** — keeps PowerSpectrumPanel focused on its
  dense PSD case.

  Files: `FixedFreqPanel.{h,cpp}` (new), catalog descriptor for
  `h5:reorient_spectral_density` with `static.fixed_freq` mode (need
  to add this to `isPanelMode` + dialog), controller builder + rebuild
  dispatch.

### Phase L-4: Feature add (~80 lines)

- **L-4 (LATER-9): SequenceBarPanel multi-channel overlay.**
  Today each SequenceBarPanel renders one scalar across residues. A
  twin-y overlay would let a researcher compare S² and τ_e on the
  same axis. Add an `addOverlay(scalar, color, label)` method to
  SequenceBarPanel that pushes a second SequenceBarRow vector. Paint
  renders the second series with offset bars (use the kind-offset
  machinery already in place — assign overlay rows kind=high so
  they slot to the right of the primary rows).

  Twin-y axis: scale each series independently if their units differ.
  Add a small per-series y-axis label on the right margin.

  Files: `SequenceBarPanel.{h,cpp}` only. Test: extend
  `dashboard_model_tests.cpp` with a `_data()` row covering overlay
  rendering (assertions on bar count + position).

## Phase order rationale

L-1 first because it's defensive and uncontroversial — sets the
floor for what "polished" means. L-2 next because the chi extension
and the panel-collection refactor exercise the chassis with real data.
L-3 last because tensor-glyph rendering is novel VTK work that
benefits from a stable chassis underneath. L-4 is optional polish
that doesn't block anything else.

After each phase: `ctest --output-on-failure` from
`h5-reader/build/linux-debug` must remain 100% green.

## Verification

```bash
cd /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-debug
cmake --build . -j
ctest --output-on-failure
# 3/3 must pass

# Manual end-to-end after L-3:
./h5reader --rest 9988 \
  /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft &
python ../../udp_listen.py    # UDP debug stream
# In app:
#   - Open Add Signal, confirm "Bar (sequence)", "Curve (lag)",
#     "Chord (coupling)" mode chips appear.
#   - Add iRED S² → SequenceBarPanel renders, survives playhead ticks.
#   - Add Reorient S² → 3 sub-slot bars per residue (NH/CaHa/CO).
#   - Add KernelDynamics PSD on one atom → PowerSpectrumPanel with
#     13 channel overlay.
#   - Add Reorient orientation tensor → ellipsoid glyph on the bond
#     in the 3D scene.
#   - Add Reorient J(ω) → FixedFreqPanel with 5 markers.
```

## Out of scope (named, not done)

These remain deferred even beyond the LATER list:

- Lagged cross-correlation matrix for KernelCoherence (producer
  documented `lagged_cross_correlation="deferred"`).
- Side-chain Lipari-Szabo (Reorient covers backbone NH/CaHa/CO only;
  the producer documents methyl/aromatic/CH as a deferred follow-up).
- Sonification synth (deferred to v2 per the original scope doc).
- Eigenvalue-spectrum global "Trajectory diagnostics" panel for iRED.
- Library / producer-side changes (`feedback_extractor_untouchable`).
- True heterogeneous panel order across all three buckets (L-2b only
  unifies the collection — controller-driven ordering may need
  additional sequencing logic if researchers want strip panels mixed
  with owned panels in arbitrary order; for now, all owned panels
  paint after all strip + spectrum panels).

## Suggested kickoff prompt for next session

> Resuming the h5-reader 5-TR surfacing project. Phases A-H + the 4
> NOW cleanup items landed in the previous session; the LATER plan is
> at `h5-reader/notes/PLAN_LATER_ITEMS_2026-05-29.md`. Load the
> qt6-cpp skill, read OBJECT_MODEL.md + PATTERNS.md in full per
> `feedback_session_depth`, then read the LATER plan. Start at Phase
> L-1 (5 housekeeping items) and progress through L-2, L-3, L-4. Each
> phase ends with `ctest --output-on-failure` green. Read the
> CLEANUP_AFTER_PHASE_AH_2026-05-29.md and SCOPE_NEW_TRS_2026-05-29.md
> docs as background. The original Phase A-H plan is at
> `/home/jessica/.claude/plans/validated-questing-metcalfe.md`.

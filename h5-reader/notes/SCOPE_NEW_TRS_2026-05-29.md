# Scope: surfacing the five new TRs in h5-reader

Date: 2026-05-29 (next-session brief)

Five TrajectoryResults landed producer-side this evening:

- `DihedralAutocorrelationTrajectoryResult`
- `IRedOrderParameterTrajectoryResult`
- `KernelCoherenceTrajectoryResult`
- `KernelDynamicsTrajectoryResult`
- `ReorientationalDynamicsTrajectoryResult`

All five are registered together in `src/RunConfiguration.cpp:249-253`. None
declares peer `Dependencies()`; attach order is free. Producer side is
complete — this doc is for the h5-reader side.

## Per-TR data shape

### 1. DihedralAutocorrelationTrajectoryResult
H5 group: `/trajectory/dihedral_autocorrelation/`

| dataset | shape | dtype | meaning |
|---|---|---|---|
| `phi_acf`, `psi_acf` | (R, L) | float64 | per-residue circular ACF over lag grid; C(0)=1; C(k) in [-1, 1] |
| `chi_acf` | (R, 4, L) | float64 | same, per-chi-index |
| `phi_corr_time_ps`, `psi_corr_time_ps` | (R,) | float64 | 1/e decorrelation time |
| `chi_corr_time_ps` | (R, 4) | float64 | same, per chi |
| `phi_defined`, `psi_defined` | (R,) | uint8 | NaN-distinguishing presence mask |
| `chi_defined` | (R, 4) | uint8 | same |
| `residue_index_per_atom` | (N,) | int32 | atom-axis broadcast |
| `lag_frames`, `lag_times_ps` | (L,) | uint64/float64 | grid metadata |

Native axis: **residue**. Channels per residue: `{phi, psi, chi[0], chi[1], chi[2], chi[3]}`.

### 2. IRedOrderParameterTrajectoryResult
H5 group: `/trajectory/ired_order_parameters/`

| dataset | shape | dtype | meaning |
|---|---|---|---|
| `s2_ired` | (M,) | float64 | per-N-H-vector iRED order parameter |
| `eigenvalues` | (M,) | float64 | descending; first 5 are tumbling modes |
| `residue_index` | (M,) | int32 | parent residue of each N-H vector |
| `n_atom`, `h_atom` | (M,) | int32 | atom indices defining each vector |

Native axis: **residue** (mapped via `residue_index`), one entry per non-Pro
residue. M ≈ R − Pro − terminus. **Static, not time-series.**

### 3. KernelCoherenceTrajectoryResult
H5 group: `/trajectory/kernel_coherence/`

| dataset | shape | dtype | meaning |
|---|---|---|---|
| `correlation_matrix` | (N, C, C) | float64 | per-atom 13×13 Pearson, symmetric, diag=1 |
| `channel_names`, `channel_units` | (C,) | string | C = 13 kernel reductions |

Native axis: **atom**, sub-keyed by channel pair `(channel_a, channel_b)`.
**Static, matrix-valued.**

### 4. KernelDynamicsTrajectoryResult
H5 group: `/trajectory/kernel_dynamics/`

| dataset | shape | dtype | meaning |
|---|---|---|---|
| `acf` | (N, C, L) | float64 | per-atom-per-channel autocorrelation over lag |
| `power_spectrum` | (N, C, F) | float64 | Parzen PSD, ≥ 0 |
| `decay_time_ps` | (N, C) | float64 | dt × sum(rho) to first zero crossing |
| `peak_freq_per_ps` | (N, C) | float64 | argmax frequency (excl. DC bin) |
| `spectral_centroid_per_ps` | (N, C) | float64 | sum(f S) / sum(S) |
| `channel_names`, `channel_units` | (C,) | string | C = 13 |
| `lag_frames`, `lag_times_ps` | (L,) | uint64/float64 | lag grid |
| `frequencies_per_ps` | (F,) | float64 | one-sided spectrum grid |

Native axis: **atom**, sub-keyed by **channel**. **Lag- and freq-domain.**

### 5. ReorientationalDynamicsTrajectoryResult
H5 group: `/trajectory/reorientational_dynamics/`

| dataset | shape | dtype | meaning |
|---|---|---|---|
| `bond_vector_autocorrelation` | (V, L) | float64 | per-vector internal C_I(k), body frame |
| `bond_vector_autocorrelation_lab` | (V, L) | float64 | per-vector lab-frame TCF |
| `order_parameter_S2` | (V,) | float64 | Henry–Szabo S² |
| `lipari_szabo_tau_e` | (V,) | float64 | ps, area method |
| `bond_orientation_tensor` | (V, 3, 3) | float64 | body-frame ⟨u⊗u⟩ |
| `vector_kind` | (V,) | uint8 | 1=NH, 2=CaHa, 3=CO |
| `owning_atom`, `tail_atom`, `head_atom` | (V,) | int32 | atom identity |
| `residue_index` | (V,) | int32 | parent residue |
| `spectral_density_j` | (V, 5) | float64 | J at [0, ω_N, ω_H−ω_N, ω_H, ω_H+ω_N] (NH only; NaN elsewhere) |
| `relaxation_R1`, `relaxation_R2` | (V,) | float64 | s⁻¹ (NH only) |
| `relaxation_NOE` | (V,) | float64 | dimensionless (NH only) |
| `relaxation_larmor_freqs_rad_per_s` | (5,) | float64 | the J sampling frequencies |

Native axis: **residue** (mapped via `residue_index`), with `vector_kind`
acting as a channel discriminator within the residue. V ≈ 3R for backbone
NH + CaHa + CO.

## Key simplification: analogous selection-keyed groupings per axis

h5-reader already groups "a set of things for a residue selection" — pick
a residue, see the signals keyed to that residue as one ordered block.
The same pattern is in place for rings via the recent
`AromaticRing`/`SaturatedRing` reveal work. The new TRs lean on the same
architectural shape extended to two more axes:

| Axis | Selection unit | New TRs that key here | Channel family |
|---|---|---|---|
| **Atom** (existing) | one atom | `KernelDynamics`, `KernelCoherence` | 13 kernel reductions per atom |
| **Residue** (existing) | one residue | `DihedralAutocorrelation` | `{phi, psi, chi0..3}` per residue |
| **Bond vector** (new — analogous shape) | one (residue, kind) bond vector | `IRedOrderParameter`, `ReorientationalDynamics` | scalars + curve + tensor per vector; the `vector_kind` enum `{NH, CaHa, CO}` is the sub-discriminator within the residue |

The bond-vector axis IS its own first-class grouping — it just gets
built using the same plumbing pattern as residue/ring, not a different
shape. `SignalAxis::BondVector` + `BondVectorAnchor` + a per-bond-vector
candidate enumerator in `NearbySignalModel` (analogous to how
`AromaticRing`/`SaturatedRing` candidates were added). Reveal maps a
bond-vector binding back to `[head_atom, tail_atom]` for highlighting,
the same way `RingMembership` maps back to the parent ring's atoms.

This keeps each axis honest:
- `IRedOrderParameter`'s S² is a property of the N-H *vector*, not of
  the residue. The residue is its parent; the vector is its identity.
- `ReorientationalDynamics`'s tensor + S² + tau_e are per-vector. A
  residue with no N-H (Pro) still has Cα-Hα and C=O vectors; the axis
  needs to express that.
- `DihedralAutocorrelation` is genuinely per-residue (phi/psi/chi are
  residue-scoped torsions), so it stays on the residue axis.

The infrastructure pattern is shared across axes: each axis has its own
typed buffer (`QtPerResidueScalar`, `QtPerBondVectorScalar`, etc.), its
own candidate enumerator in `NearbySignalModel`, its own anchor variant
in `SignalAnchor`, and its own reveal mapping in `SceneRevealOverlay`.
The cost per new axis is modest because the surrounding machinery
(picker model/view, descriptor catalog, strip rendering) is generic over
the axis.

## What new infrastructure is actually needed

After the residue-grouping simplification, the genuinely new pieces are:

1. **Curve display over lag** (and **over frequency**). Three of the five
   TRs emit `(*, L)` or `(*, F)` arrays whose natural display is a curve
   plotted vs. lag (or frequency), not vs. trajectory frame. The dashboard
   strip widget is the wrong host — it draws trajectory time series.
   Options:
   - Dock a sibling **"Static curve" panel** next to the strip area.
     Triggered when an active signal's display mode is `static.curve.lag`
     or `static.curve.freq`. Renders one curve per selected channel.
   - Or modal: right-click a residue/atom → "Show ACF" → opens a popup
     curve viewer. Less integrated but simpler v1.

2. **Matrix display** for `KernelCoherence`. One per atom, 13×13. Options:
   - Static heatmap panel with row/column = channel name.
   - Or projection to scalar via channel-pair anchor (`static.scalar` over
     a `RingContributionPair`-like channel-pair anchor, exposing one
     element of the matrix at a time). The matrix-panel route is more
     informative; the projection route reuses existing strip rails.

3. **Mat3 per residue/vector** for `bond_orientation_tensor`. Already have
   `SphericalTensor` displays at atom axis; need the analogous Mat3 typed
   buffer and a tensor table / ellipsoid glyph display mode at residue
   axis. Could borrow `static.tensor.glyph` from existing tensor modes.

4. **Lag-domain typed buffer**. The existing `QtAutocorrelation` is
   atom-axis for `bs_t0_autocorrelation`. Generalise to
   `QtPerRowCurve<axis>` (atom, residue, or "vector-block-row") so all
   three new lag arrays can use one buffer class. Same generalisation
   covers the (F,) frequency-grid PSD.

5. **Per-axis static scalar buffers**. Currently the residue/ring axes
   carry only time-series. New typed buffers:
   - `QtPerResidueScalar` — for `phi_corr_time_ps` (DihedralAutocorrelation).
   - `QtPerBondVectorScalar` — for `s2_ired`, `order_parameter_S2`,
     `lipari_szabo_tau_e`, R1, R2, NOE.
   Same shape as `QtScalarTimeSeries` minus the time axis. Single class
   pattern, two axis instantiations.

6. **`SignalAxis::BondVector` end-to-end**. New axis on the same plumbing
   as the existing residue/ring axes:
   - `BondVectorAnchor{vector_index}` variant arm in `SignalAnchor`.
   - `SignalAxis::BondVector` enum addition.
   - `QtBondVectorTable` typed buffer holding `(vector_kind, residue_index,
     owning_atom, tail_atom, head_atom)` per row, parsed once at H5 open.
   - `NearbySignalModel::CandidateKind::BondVector` enumerator + label
     helper (`"Lys17 N-H"`).
   - `SceneRevealOverlay` mapping: `BondVectorAnchor` → highlight
     `{tail_atom, head_atom}` (analogous to the existing bond-anchor reveal).
   - `AnchorMatchesAxis` widening so a vector-axis descriptor accepts
     either a `BondVectorAnchor` or a parent-`ResidueAnchor` (the
     residue-grouping ergonomic: pick a residue, get all its vectors as
     a sub-list, analogous to how a Ring axis accepts Aromatic/Saturated
     subkinds today).

7. **Per-channel sub-key plumbing**. The 13-channel kernel axis in
   `KernelDynamics`/`KernelCoherence` and the 3-channel bond-vector-kind
   axis sit *inside* their parent atom/residue selection. The existing
   `PerClassBlock` (residue × class) is the closest pattern. Audit
   whether it generalises to `(parent_axis × channel_axis)` cleanly or
   if `PerAtomPerChannel` / `PerResiduePerKind` deserve their own
   typedefs.

The new `SignalValueShape` enum additions:
- `Matrix` (for (C, C) per-row matrices)
- `Spectrum` (for (F,) frequency-grid curves)
- `CurveOverLag` (for (L,) lag-grid curves; could be unified with Spectrum
  as `CurveOverAxis<Lag|Freq>`)

## Per-TR cost estimate

| TR | Axis | New typed buffers | Catalog descriptors | New value shapes | New display modes | Estimated diff |
|---|---|---|---|---|---|---|
| DihedralAutocorrelation | residue | `QtPerRowCurve<Residue>` + `QtPerResidueScalar` + `QtPerResiduePerKindBlock` (phi/psi/chi₀..₃) | 2 (scalar reductions, ACF curves) | `CurveOverLag` | `static.curve.lag.animated` (in inspect dock), `static.bar.sequence` | ~120 lines |
| IRedOrderParameter | bond vector | `QtBondVectorTable` (identity) + `QtPerBondVectorScalar` | 1 (S² + eigenvalues table) | — (uses new BondVector axis) | `static.bar.sequence` (NMR-convention per-residue bars) | ~110 lines (incl. axis plumbing) |
| KernelCoherence | atom + channel-pair | `QtPerAtomMatrix` | 1 | `Matrix` | `static.chord.coupling` (chord diagram, NOT heatmap) | ~110 lines + `ChordCouplingPanel` subclass (~80 lines) |
| KernelDynamics | atom + channel | `QtPerAtomChannelCurve` (covers ACF + PSD) + `QtPerAtomChannelScalar` (scalar reductions) | 3 (acf, psd, scalar reductions) | `CurveOverLag`, `CurveOverFreq` | `static.spectrum.power` (line plot of power vs frequency; supports overlay of multiple channels), `static.curve.lag.animated`, `static.audio.sonify` (literal WAV of the 13-voice mix, v2), `strip.scalar` per channel | ~180 lines + `PowerSpectrumPanel` + `LagDecayPanel` subclasses (~120 lines combined; audio synth deferred) |
| ReorientationalDynamics | bond vector | `QtBondVectorTable` (reused) + `QtPerBondVectorScalar` (reused) + `QtPerRowCurve<BondVector>` + `QtPerBondVectorMat3` + `QtPerBondVectorFixedFreqBlock` (5-J + R1/R2/NOE) | 4-5 | `Mat3PerRow`, `FixedFreqBlock` | `static.bar.sequence` (S², τ_e, R1, R2, NOE per residue), `static.tensor.glyph` (body-frame ellipsoid), `static.curve.lag.animated`, `static.fixed_freq` | ~220 lines |

Total estimate: ~600 lines of new code + ~150 lines of test coverage +
**4 new panel subclasses** on the existing `StripStackWidget` chassis
(per-sequence bar, animated lag curve, chord-coupling, power-spectrum
line plot — no separate widget classes). Sonification synth deferred
to v2. Spread across roughly 7 new typed buffer classes (mostly small)
and 6-7 new value shapes / display mode IDs.

The four new panels share the existing stack widget's painter, layout,
reveal-button, hover, and selection plumbing; what each panel actually
adds is just its rendering code plus optional gesture handlers.

## Preferred visualisations (decided 2026-05-29)

The visualisation strategy is **accessibility-first**: each TR's natural
output gets mapped to a viz idiom the brain reads pre-cognitively, so
the data lands for non-specialists (advisor, committee, biologists)
without forcing them to learn a new visual grammar. Heatmaps are an
implicit literacy test the audience shouldn't have to pass; bar charts,
chord diagrams, equaliser displays, decaying curves, and actual sound
are all things people process without instruction.

The framing match: this work is "listening to proteins sing" — the
dynamics are auditory/temporal in spirit. The visualisation choices
honour that frame where the data shape permits.

### Per-TR preferred viz

- **KernelCoherence (13×13 per-atom Pearson matrix) → chord / arc diagram.**
  Channel labels around a circle; arcs between coupled pairs; thickness
  = `|r|`, hue = sign. Reads at a glance: "this atom couples BS strongly
  with HBond, weakly with everything else." Reuses Qt painting; no layout
  solver; ~80-line widget. Standard idiom in network and genomics viz —
  doesn't require fluency. **Heatmap is explicitly rejected** as the
  default: it surfaces values but not the coupling-structure story.

- **KernelDynamics (per-atom-per-channel ACF + PSD + scalar reductions)
  → power spectrum (line plot of power vs frequency) AND optional
  literal sonification.** Everyone loves a power spectrum because
  everyone *already reads* power spectra — audio, EE, physics,
  spectroscopy all share the idiom. The PSD on disk is a `(F,)` curve
  per (atom, channel); the natural display is exactly that: a line plot
  with frequency on x (linear or log), power on y (typically log).
  Selecting one (atom, channel) shows one curve; selecting multiple
  channels stacks them as colour-coded traces on the same axes with a
  legend (Qt Charts handles this cleanly). For a single atom the user
  can dump all 13 channels overlaid in one click to see which kernels
  carry the high-frequency content and which are slow.

  The ACF curves go to the inspect dock as animated decay traces (you
  watch the curve relax frame-by-frame as playback advances).
  Per-channel scalar reductions (decay_time_ps, peak_freq_per_ps,
  spectral_centroid_per_ps) are conventional strip-mode signals along
  the existing strip surface.

  The literal sonification path lands as a sibling display mode
  (`static.audio.sonify`): map each kernel's peak frequency to an
  audible tone via Qt's `QAudioSink` + a small sine-bank synth (~50
  lines), mix the 13 voices, export per atom or play directly. Coupling
  manifests as harmony; decoupling as dissonance. Defer the synth to a
  second pass — the power-spectrum line plot is the load-bearing v1.

- **ReorientationalDynamics + IRedOrderParameter (per-residue/vector
  S², τ_e, R1, R2, NOE) → per-residue bar chart along the sequence.**
  This is the standard NMR-relaxation viz idiom — every Lipari-Szabo
  paper uses it. S² as bars 0→1 with residue index on x, optionally
  shaded by `vector_kind` for ReorientationalDynamics's NH/CaHa/CO
  triple. R1/R2/NOE get the same treatment in adjacent panels.
  Universally legible to anyone who's read an NMR-dynamics paper;
  comparable to published baselines (Brüschweiler, Showalter, Bremi)
  without translation. The `bond_orientation_tensor` Mat3 per vector
  becomes a standard 3D ellipsoid glyph attached to the highlighted
  bond in the existing scene overlay.

- **DihedralAutocorrelation (per-residue corr_time_ps + ACF curves over
  lag) → per-residue bar chart + animated curve in inspect dock.**
  `phi_corr_time_ps` and `psi_corr_time_ps` are the bar-along-sequence
  story (the same Lipari-Szabo-paper idiom). The ACF curves themselves
  go to the inspect dock as decaying curves; the dock animates the
  curve relaxing as the playhead advances, so the user feels the
  decorrelation timescale directly rather than reading a number.

### Cross-cutting: all panels live in `StripStackWidget`

Decided 2026-05-29: the new viz idioms are **panel subclasses** on the
existing `StripStackWidget` chassis, not separate widgets. The recent
strip-widget rewrite already factored out `AbstractStripPanel`
polymorphism (`TemporalStripPanel` + `SpectrumStripPanel` are siblings
sharing the same painter, geometry, reveal-button, hover, selection);
the new viz types extend the same hierarchy:

```cpp
class AbstractStripPanel { ... };                                // exists
class TemporalStripPanel  : public AbstractStripPanel { ... };   // exists
class SpectrumStripPanel  : public AbstractStripPanel { ... };   // exists
class ChordCouplingPanel  : public AbstractStripPanel { ... };   // new
class PowerSpectrumPanel  : public AbstractStripPanel { ... };   // new
class LagDecayPanel       : public AbstractStripPanel { ... };   // new (animated)
class SequenceBarPanel    : public AbstractStripPanel { ... };   // new
```

Each new panel just implements `paint()`, `hasRevealBinding()`,
`revealBinding()`, and optionally `plotContains()` / `tooltipLine()`.
The stack handles layout, scrolling, painter setup, header rendering,
reveal-button gesture, hover, selection, surface ownership. No new
widget classes; no separate "inspect" dock; one event model.

**One stack, mixed panel types**: a researcher comparing S² (sequence
bar) against the BS power spectrum (line plot) against the dihedral ACF
(animated decay) for the same residue gets all three stacked in one
scrollable column, not scattered across three docks. The display mode
on each active signal selects its panel subclass; the stack composes
whatever heterogeneous set the user has built up.

**Small extensions to `AbstractStripPanel`** needed to accommodate the
new viz types:

- `virtual std::optional<double> preferredAspect() const` — chord wants
  square, sequence-bar wants wide; default `nullopt` (use full plot
  rect, as today). Stack letterboxes within the assigned rect when set.
- `virtual void mousePressInPlot(...)`, `mouseMoveInPlot(...)` —
  defaults to no-op. Temporal panels keep their existing
  drag-select-range gesture; chord panels add click-arc-to-highlight;
  sequence-bar adds click-bar-to-focus.
- `PaintContext` consumption is per-panel. Temporal + animated-lag
  panels honour `currentFrame` / `time`; static panels (chord,
  sequence-bar, power-spectrum) ignore it and don't repaint per frame.
  The stack's `update()` triggers are already conditional on what
  changed, so this is just a matter of each panel reading what it
  needs.

The cost saving is real: the four "new widgets" in the per-TR cost
table become four panel subclasses sharing one container. Widget
infrastructure cost drops to ~zero; only the per-panel painting code
remains.

### What this changes about the cost estimate

The chord widget replaces a more expensive heatmap+colormap widget at
roughly the same line count. The equaliser + sonification adds ~50
lines for the WAV path but is otherwise the same widget as the chord
(custom QPainter). The bar-along-sequence widget is the most-reused
piece (4 of the 5 TRs use it), so building it well early is the
highest-leverage move. The animated lag-curve widget is a single
reused class across DihedralAutocorrelation, KernelDynamics, and
ReorientationalDynamics.

Net: same ~700-line estimate, with the viz pieces shifted toward
idioms that don't require teaching the viewer a new visual grammar.

## Suggested task order for next session

1. **Design session first** (~30 min): nail the value-shape vocabulary
   (`Matrix`, `Spectrum`, `CurveOverLag`, possibly unified into
   `CurveOverAxis`) and the display mode IDs. No code yet. Choose between
   sibling-dock vs modal for static curve display.
2. **Generalise `QtAutocorrelation` → `QtPerRowCurve<axis>`**. Covers
   `bs_t0_autocorrelation` (existing) plus the new lag/freq arrays.
3. **Add `QtPerResidueScalar`**. Touches IRed, ReorientationalDynamics,
   DihedralAutocorrelation reductions.
4. **Land in payoff order**:
   - **KernelDynamics first** (the headliner — kernel-by-kernel ringing
     across the protein; all 13 channels become per-atom strips of their
     scalar reductions, plus the ACF/PSD curves as static-curve detail).
   - **ReorientationalDynamics second** (Lipari-Szabo proper — the
     decade-mainstream NMR-dynamics output; per-residue S² and τ_e are
     residue-axis scalars, R1/R2/NOE per residue, body-frame tensor
     glyph).
   - **DihedralAutocorrelation third** (per-residue dihedral
     decorrelation timescales as the "clicking" view; six channels per
     residue, scalars + curves).
   - **IRedOrderParameter fourth** (one number per residue, smallest
     surface; can land in a session quarter).
   - **KernelCoherence last** (matrix heatmap is the heaviest viz; can
     land projection-to-scalar first if the heatmap widget is deferred).
5. **`kDensePaths` ↔ `denseH5Plan` lockstep regression test**: add to
   `dashboard_model_tests` a test that loads the catalog, iterates every
   `DenseH5Trajectory` descriptor with Valid status, and asserts the
   descriptor's storagePath appears in `kDensePaths`. Cheap; catches
   future TR landings that forget step 9 of the checklist.
6. **REST + pytest** sanity per TR: one round-trip
   `GET /dashboard/signals` filter for each new descriptor id; assert it
   is listed and reports the expected `value_shape`. The screenshot
   coverage can wait for the static-curve widget to exist.

## Open questions to settle in the design session

Settled (recorded here, not relitigated):
- ~~Static curve display: dock or modal?~~ → **Sibling "Inspect" dock**
  (per the cross-cutting section above).
- ~~Matrix display: heatmap vs channel-pair projection?~~ → **Chord
  diagram** (per the preferred-viz section). Heatmap explicitly rejected.

Still open:
- **`PerClassBlock` extension vs new `Per{Atom,Residue,BondVector}Channel`
  typedefs?** Audit the existing block's shape first; the cleanest
  answer might be one generic `PerRowChannelBlock<RowAxis, ChannelAxis>`
  template.
- **Bond-vector axis sub-listing under residue.** Picking a residue
  should surface its bond vectors as a sub-list (analogous to ring
  candidates today). Decide whether this is a separate row in
  `NearbySignalModel` per vector, or a single residue row that expands
  into vectors on selection.
- **Channel naming convention for the 13 KernelDynamics channels.** They
  are emitted as strings (`channel_names`); the picker UI should group
  them in the canonical thesis-narrative order (BsT0, BsAbsT2, HmT0, …)
  rather than alphabetically. Decide whether the catalog enforces order
  or trusts the H5 string array.
- **Whether to expose the eigenvalue spectrum for IRed.** The (M,)
  eigenvalues array is a one-off per-trajectory diagnostic, not a
  per-residue signal — could go on a global "Trajectory diagnostics"
  panel rather than the dashboard.
- **Sonification scope.** EQ display lands in v1; the literal WAV
  synth (`static.audio.sonify`) is a second pass — decide whether to
  scope it for the same TR sweep or hold for a dedicated viz session
  after the visual displays are bedded in.

## Cross-references

- Producer source: `src/{DihedralAutocorrelation,IRedOrderParameter,KernelCoherence,KernelDynamics,ReorientationalDynamics}TrajectoryResult.{h,cpp}`
- Producer registration: `src/RunConfiguration.cpp:249-253`
- Master pipeline inventory: this file's preceding session conversation (5-lane survey of how a TR moves from src/ → SDK catalog → C++ generated catalog → QtTrajectoryH5 typed buffers → SignalCatalog → DashboardDisplayController sampler → strip)
- Test framework that should grow with the new descriptors:
  `h5-reader/tests/dashboard_model_tests.cpp`, `h5-reader/tests/rest/`
- Memory entries to consult: `project_streaming_observer_for_sigma_pred`
  (L-S deliverable framing — these TRs are exactly the producer-side L-S
  work that memory describes), `feedback_correlate_not_match`,
  `feedback_methods_accumulate` (no retirement; new methods land
  alongside).

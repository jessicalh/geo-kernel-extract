# Stabilized Local View Design

Date: 2026-05-29

Status: SUPERSEDED (trued 2026-06-04). The stabilisation WAS implemented since this note —
`TransformedConformation` (display-space transforms) + `CameraComposer` (camera locks); plane lock is
real and tested. This note remains the original design *rationale*; for the ACTUAL built state, the gap to
its three-atom-wedge design, and the cruft it carried, see `STABILISATION_FEATURE_EVAL_2026-06-04.md`. The
radius / local-isolation view it was meant to enable is still UNBUILT (the planned
`AtomFilter::WithinRadiusOf` path — see `UI_STATE_OVERVIEW_2026-06-04.md`).

## Purpose

The goal is a viva/explanatory view for shielding: choose exactly three atoms,
build a local reference frame around them, and let the rest of the protein move
relative to that frame. This is a display mode for explanation. It is not a
claim that the real molecule is planar, rigid, or transformed in the source
data.

The desired effect is stronger than the current camera plane lock. The current
lock can keep the selected atoms in sight and center the camera on them. The
stabilized view should make the selected local context feel steady enough that
the surrounding motion, shielding changes, and local geometry are easier to
explain.

## Concept

The view uses an artificial thick plane or wedge tied to exactly three selected
atoms. The plane has give: it stabilizes the view over a range of real molecular
motion without pretending the selected atoms are frozen or perfectly planar.

The selected atoms still move legally. Their residual movement is part of the
display, not an error to hide. The mode should try hard to keep the containing
wedge stable, with the first selected atom acting as the orientation/context
anchor, but it must not distort the chemistry to make the picture cleaner.

Changing the selection, or unselecting any of the three atoms, should stop the
mode. It only has a clear meaning for exactly three selected atoms.

## Geometry Contract

June 4 correction: the exact three-atom Kabsch path in this original design is
invalid for the current shared Kabsch implementation. That path rejects
rank-deficient/coplanar subsets, so a three-atom stabilized local view must use
a plane/wedge stabilizer instead of treating ordinary Kabsch as the ready
implementation route.

It is not possible to keep all three selected atoms exactly fixed across frames
without distorting geometry whenever their internal distances or angles change.
The implementation must choose an honest stabilization rule.

Acceptable first rules:

- Fix the centroid, keep one selected edge or first-atom direction stable, and
  let the third atom show its real residual motion in or near the plane.
- Use a dedicated plane/wedge stabilizer over the three selected atoms to
  minimize their apparent motion while leaving residual deformation visible.

The artificial plane should be thick enough to tolerate ordinary local motion.
If the selected atoms move outside the useful range, the UI should degrade
clearly rather than silently inventing a false stable plane.

## Implementation Guardrails

Do not rewrite raw conformation data or scientific values. `Conformation` and
`Conformation::atomPosition` remain the source of truth.

Implement the stabilized view as a display-space transform in the rendering
path, probably owned by `MoleculeScene`. While active, the molecule geometry,
selection overlay, measurement overlay, reveal overlay, picking, and camera
must all agree about the same display-space coordinates. A mixed state where
the molecule is transformed but picking or overlays still use raw coordinates
would be worse than no feature.

Prefer a small explicit helper such as `displayPosition(frame, atom)` or an
equivalent local transform object over a broad new architecture. This is a
view-mode problem, not a new data model.

Keep the existing playback path healthy. Playback already works, so the
stabilized view should be opt-in and should not disturb ordinary frame advance,
camera controls, or dashboard sampling.

## UI Shape

The mode can be exposed as "Stabilized view" or "Local frame". It should become
available only when exactly three atoms are selected. Any selection change
clears it.

The UI should make the active state obvious. It can be separate from the
current camera-only "Plane lock", because these are different promises:

- Camera plane lock: keep the three atoms in sight, center the camera, and look
  along the current selected plane.
- Stabilized local view: transform displayed positions into a forgiving local
  explanatory frame so the selected context stays visually steady.

## Test Requirements

Keep the current camera plane-lock smoke test. Add a separate stabilized-view
smoke before enabling this for regular use:

- Load the trajectory fixture.
- Select three atoms and activate the stabilized view.
- Advance through a frame window.
- Verify the selected local frame remains stable within a declared tolerance.
- Verify residual selected-atom motion is still measurable rather than forced
  to zero.
- Verify molecule geometry, selection/reveal overlays, measurement overlays,
  picking, and camera use the same display-space transform while active.
- Verify ordinary playback still renders with the feature inactive.

## Hold Point

This is recorded for the next implementation pass. For now, the conservative
camera plane lock is useful only as a centering/orientation aid. The stronger
stabilized explanatory view should not be attempted until the display-space
transform and tests can be done carefully.

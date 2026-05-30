// CameraComposer — owner of the per-frame absolute camera write.
//
// Given the current frame, the active CameraMode + OrientationPolicy,
// and any accumulated user-gesture deltas, produces the camera state
// and writes it via SetFocalPoint / SetPosition / SetViewUp +
// OrthogonalizeViewUp. NEVER delta-translates the camera. NEVER calls
// Render — the scheduler (MoleculeScene::requestRender) owns that.
//
// Stage 3 of spec/viewport_pipeline_2026-05-30.md. Replaces:
//   * MoleculeScene::CameraPlaneLock (struct + applyCameraPlaneLock)
//   * MoleculeScene::focusCameraOnReveal one-shot dihedral math
//   * The centroid-delta camera translation block in setFrame
//
// All math lives in FitTargetMath.h (pure functions); this class is
// the orchestrator that dispatches on CameraMode::Kind and assembles
// the absolute write each frame.
//
// User-gesture delta storage: kept per-composer (not per-mode-variant).
// Each setMode call resets the delta to zero — agent decision per the
// implementation prompt, §4-b: each lock acquisition is a fresh start.
// applyGesture is the slot the eventFilter calls; it mutates the delta
// and asks for a render via the scene's requestRender.
//
// Thread affinity: GUI thread only. ASSERT_THREAD at every entry point.

#pragma once

#include "CameraMode.h"
#include "OrientationPolicy.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/Conformation.h"
#include "../model/QtProtein.h"
#include "../model/Types.h"

#include <QObject>
#include <QPointer>

#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <cstddef>
#include <optional>
#include <vector>

namespace h5reader::app {

// User-gesture delta accumulated by the eventFilter. The composer
// composes this on top of the per-frame natural fit. Resets to zero
// every setMode call.
struct CameraGesture {
    enum class Kind {
        Azimuth,   // rotate around current ViewUp by deltaRadians
        Elevation, // rotate around right vector by deltaRadians
        Roll,      // rotate around view direction by deltaRadians
        Pan,       // translate focal point by (dxPx, dyPx) in screen space
        Dolly,     // multiplicative zoom: position -> focal + (position - focal) / factor
    };

    Kind   kind         = Kind::Azimuth;
    double deltaRadians = 0.0;   // Azimuth / Elevation / Roll
    double dxScreenPx   = 0.0;   // Pan
    double dyScreenPx   = 0.0;   // Pan
    double dollyFactor  = 1.0;   // Dolly; > 1 zooms in, < 1 zooms out
};

class MoleculeScene;

class CameraComposer final : public QObject {
    Q_OBJECT

public:
    // The renderer is shared with the main scene; the composer writes to
    // its active camera each frame. protein + conformation are non-owning
    // pointers — the composer is a child of MoleculeScene (Qt parent), so
    // lifetime tracks the scene which tracks the loader result.
    CameraComposer(vtkSmartPointer<vtkRenderer>  renderer,
                   const model::QtProtein*       protein,
                   model::Conformation*          conformation,
                   QObject*                      parent = nullptr);
    ~CameraComposer() override;

    // Switch the active mode + orientation policy. Bumps the reference
    // frame for any fit (Plane normal-sign continuity, Subset Kabsch
    // reference positions) so the camera state at the next write is
    // continuous with the user's current view. Resets accumulated user
    // deltas to zero. Idempotent on equal (mode, policy, currentFrame).
    // Emits modeChanged() at the end if anything moved.
    void setMode(CameraMode mode,
                 OrientationPolicy policy,
                 std::size_t currentFrame);

    // Absolute camera write for frame t. Returns false if the mode's
    // atom subset is degenerate at this frame (collinear, missing atoms,
    // zero-length axis); caller (MoleculeScene::setFrame) decides whether
    // to fall back to free camera. Writes happen only for non-Free modes;
    // Free mode is a no-op so the user's gesture-driven camera state
    // persists frame to frame.
    [[nodiscard]] bool write(std::size_t t);

    // Apply a user gesture (from CameraInputFilter). Updates the
    // accumulated delta state; the next write(t) re-composes the fit
    // and the delta. For CameraMode::Free, the delta is applied directly
    // to the camera state (no fit to compose with).
    void applyGesture(const CameraGesture& g);

    // Read-only inspectors for REST + tests.
    CameraMode mode() const { return mode_; }
    OrientationPolicy policy() const { return policy_; }

signals:
    // Emitted at the end of a successful setMode call.
    void modeChanged();

private:
    // Per-frame fit dispatch helpers, one per Kind. Each returns false
    // on degenerate input; the composer falls back to no write on that
    // frame (keeping the previous frame's camera state).
    bool writeAtom(std::size_t t);
    bool writeBond(std::size_t t);
    bool writeDihedral(std::size_t t);
    bool writePlane(std::size_t t);
    bool writeSubset(std::size_t t);
    bool writeFree(std::size_t t);  // strictly a no-op; recorded here for shape

    // Capture per-mode initial state from the current camera at setMode
    // time. Different modes need different fields (distance for Atom/Bond/
    // Plane; viewUp for Plane; subset reference positions for Subset).
    void captureInitialState(std::size_t referenceFrame);

    // Convenience: validate the active mode's atom indices against the
    // protein at construction time and per write. Returns false (with a
    // logged warning) if any index is out of range.
    bool validateAtomsForCurrentMode() const;

    // Resolve atom positions for the active mode at frame t. Returns
    // empty vector if any atom is out of range. Helper for the per-mode
    // writers.
    std::vector<model::Vec3> readAtomPositions(std::size_t t) const;

    // Compose the user gesture delta on top of an absolute camera state
    // and write it. focal/position/viewUp are the "natural" per-frame
    // fit output; this function rotates / pans / dollies them per the
    // accumulated delta and calls SetFocalPoint/SetPosition/SetViewUp.
    void writeCameraComposed(const model::Vec3& focal,
                             const model::Vec3& position,
                             const model::Vec3& viewUp);

    // Apply the accumulated delta (azimuth, elevation, roll, pan, dolly)
    // to a (focal, position, viewUp) triple, returning the modified one.
    // Pure transform; no VTK side effects.
    void applyDeltaToState(model::Vec3& focal,
                           model::Vec3& position,
                           model::Vec3& viewUp) const;

    vtkSmartPointer<vtkRenderer>   renderer_;
    const model::QtProtein*        protein_      = nullptr;
    QPointer<model::Conformation>  conformation_;

    CameraMode        mode_   = FreeMode();
    OrientationPolicy policy_ = DefaultPolicy();

    // Per-mode captured state at setMode time.
    double distance_ = 1.0;   // distance from focal to position; used by Atom/Bond/Plane

    // Plane mode: stored local view-up + normal-sign + last direction
    // (sign-continuity guard). Mirrors the original CameraPlaneLock
    // fields; lifted into the composer.
    model::Vec3 planeLocalViewUp_ = model::Vec3::Zero();
    double      planeNormalSign_  = 1.0;
    std::optional<model::Vec3> planeLastDirection_;

    // Dihedral mode: sign-continuity guard for the sight axis (Codex
    // finding #1). Mirrors planeLastDirection_ exactly. Stored direction
    // is the POST-flip axis used to write the camera last frame; the
    // next frame's axis is sign-flipped if it dots negative against this
    // reference, then this reference is updated. Reset to nullopt on
    // setMode transitions into/out of Dihedral so the first write picks
    // the natural axis direction.
    std::optional<model::Vec3> dihedralLastDirection_;

    // Atom mode: captured sight + up + cam-relative offset at lock
    // acquisition (Codex finding #2). The prior implementation derived
    // each frame's sight from the LIVE camera, which already contained
    // the composed user gesture — leading to drift on subsequent
    // setFrames because the gesture re-applied on top of itself. Now
    // each frame's natural pose comes from the captured reference, and
    // accumulated gestures compose exactly once per frame on top of it.
    // Sight is the unit camera-to-focal direction (normalised
    // (focal - position)) at lock acquisition; up is orthogonalised
    // against sight. atomReferenceCamRel_ holds (position - atom_t0)
    // so that each frame's camera position = atom_t + atomReferenceCamRel_
    // (preserving the captured zoom and angle).
    model::Vec3 atomReferenceSight_  = model::Vec3::Zero();
    model::Vec3 atomReferenceUp_     = model::Vec3::Zero();
    model::Vec3 atomReferenceCamRel_ = model::Vec3::Zero();

    // Bond mode: matched-shape capture for the bond lock (Codex finding
    // #2). bondReferenceMidpoint_ caches midpoint(a, b) at lock
    // acquisition so each frame's camera position composes against the
    // captured pose (camera = midpoint_t + bondReferenceCamRel_) rather
    // than the live camera that already contains accumulated gestures.
    model::Vec3 bondReferenceSight_   = model::Vec3::Zero();
    model::Vec3 bondReferenceUp_      = model::Vec3::Zero();
    model::Vec3 bondReferenceCamRel_  = model::Vec3::Zero();
    model::Vec3 bondReferenceMidpoint_ = model::Vec3::Zero();

    // Subset mode: reference positions captured at setMode time
    // (subsetAtoms_[i] -> reference position). The composer Kabsch-fits
    // the current subset to these per frame; the captured camera state
    // below is rotated by R^T so the camera follows the molecule's frame.
    std::vector<model::Vec3> subsetReference_;
    // Captured camera-to-centroid vector at lock acquisition; rotated by
    // R^T each frame and added to the current centroid to place the
    // camera at the same relative pose. Zero until captureInitialState
    // populates it for Subset mode.
    model::Vec3 subsetReferenceCamRel_  = model::Vec3::Zero();
    // Captured ViewUp at lock acquisition (orthogonalised against the
    // sight direction). Rotated by R^T each frame so the world-up of the
    // captured view stays attached to the molecule.
    model::Vec3 subsetReferenceUp_      = model::Vec3::Zero();
    // Centroid of subsetReference_ at lock acquisition; cached so the
    // per-frame translation is computed against the same anchor.
    model::Vec3 subsetReferenceCentroid_ = model::Vec3::Zero();

    // Accumulated user gesture delta. The deltas are in the camera's
    // local frame at the time the gesture fired; we re-apply them on top
    // of each per-frame fit. Reset to zero on setMode.
    double accumAzimuthRad_   = 0.0;
    double accumElevationRad_ = 0.0;
    double accumRollRad_      = 0.0;
    model::Vec3 accumPan_     = model::Vec3::Zero();  // world-space focal offset
    double accumDolly_        = 1.0;                  // multiplicative
};

}  // namespace h5reader::app

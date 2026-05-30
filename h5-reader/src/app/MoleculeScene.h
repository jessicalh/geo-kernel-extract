// MoleculeScene — VTK render pipeline for one protein + animated positions.
//
// Owns the vtkMolecule, its mapper, and the actor. Attaches the mapper
// via SetInputData (molecule is a vtkDataObject, not a pipeline source).
// On setFrame(t), updates every atom's position in place via
// vtkMolecule::SetAtomPosition + Modified() — the idiomatic VTK pattern
// for trajectory animation. Bond connectivity is static; lengths follow
// atom positions automatically because the mapper re-queries them.
//
// Two-layer renderer (per spec/viewport_pipeline_2026-05-30.md §5):
// the main renderer at layer 0 holds the molecule + ribbon/rings/field
// overlays; an overlay renderer at layer 1 holds the markers
// (MeasurementOverlay, SceneRevealOverlay). Markers paint after the
// main layer with depth reset between layers, so they remain visible
// regardless of depth occlusion — essential for the viewport harness's
// marker-blob analysis.
//
// Camera pipeline (per spec §2.3, §2.4, §2.5): per-frame absolute
// camera writes are owned by CameraComposer, dispatched on a typed
// CameraMode (Free / Atom / Bond / Dihedral / Plane / Subset). The
// public plane-lock API is now a thin shim over the composer; the REST
// surface and toolbar action continue to work without consumer
// changes. Centroid-delta camera translation is GONE.
//
// Thread affinity: GUI thread only. All VTK mutations must happen there.
// The frame-advance slot is connected to QtPlaybackController's
// frameChanged signal with the default (direct/auto) connection; since
// the controller lives on the same thread, that's safe.
//
// --- Overlay contract (applies to QtBackboneRibbonOverlay,
//     QtRingPolygonOverlay, QtFieldGridOverlay, QtBFieldStreamOverlay,
//     and any future overlay owned by MoleculeScene) ---------------
//
// 1. Build(protein, conformation) is called ONCE per Scene lifetime,
//    before any setFrame. Idempotent on the same (protein, conformation)
//    pointers — caller guards with a same-inputs check.
//
// 2. setFrame(int t) is called per frame change by MoleculeScene
//    AFTER molecule positions update and BEFORE the render is scheduled.
//    The overlay must update its VTK backing data (vtkPoints,
//    vtkImageData scalars, vtkStructuredGrid vectors, …) and call
//    Modified() as needed. The overlay MUST NOT call Render() itself;
//    MoleculeScene's render scheduler issues exactly one render per
//    setFrame via the Qt paint chain.
//
// 3. setVisible(bool) toggles actor visibility. An overlay that skips
//    expensive work while hidden (kernel eval, filter rerun) must NOT
//    run that work in setVisible; it waits for the next setFrame.
//    MoleculeScene's refreshCurrentFrame() is the path that re-invokes
//    setFrame on the current frame when a visibility flip needs data
//    populated; ReaderMainWindow's toggle callbacks call it.
//
// 4. Shared helpers (QtFrame::position, QtFrame::ringGeometry,
//    QtFrame::ringVertices, model::OrthoBasisFromNormal) live on the
//    model so every overlay reads the same interpretation of the H5.
//    Do NOT duplicate these in an overlay's anonymous namespace.
//
// 5. All VTK state mutations happen on the GUI thread. ASSERT_THREAD
//    at the top of each public method that mutates VTK objects.

#pragma once

#include "../model/Conformation.h"
#include "../model/QtProtein.h"
#include "../model/SignalDictionary.h"
#include "CameraMode.h"
#include "OrientationPolicy.h"
#include "PlaneFrameMath.h"

#include <QObject>
#include <QPointer>

#include <vtkActor.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkMolecule.h>
#include <vtkOpenGLMoleculeMapper.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <memory>
#include <optional>
#include <vector>

class QVTKOpenGLNativeWidget;

namespace h5reader::app {

class QtBackboneRibbonOverlay;
class QtRingPolygonOverlay;
class QtFieldGridOverlay;
class QtBFieldStreamOverlay;
class QtSelectionOverlay;
class MeasurementOverlay;
class SceneRevealOverlay;
class CameraComposer;

class MoleculeScene final : public QObject {
    Q_OBJECT

public:
    // Source-of-render hint set by requestRender, read by the EndEvent
    // observer (Stage 8 in spec/viewport_pipeline_2026-05-30.md) so each
    // logged render line carries the trigger that caused it.
    enum class RenderSource {
        Timer,       // QtPlaybackController::frameChanged -> setFrame
        CameraInput, // CameraInputFilter gesture
        Picker,      // double-click pick triggered an overlay update
        Rest,        // REST handler mutated state and needs a redraw
        Overlay,     // overlay visibility toggle from the UI
        Reveal,      // dashboard reveal binding
        Shutdown,    // explicit teardown render
        External,    // explicit consumer request that doesn't fit above
    };

    // Takes the VTK widget so the render scheduler can call
    // vtkWidget_->update() (the only render verb in app code per spec
    // §2.5). The widget owns the render window; the scene pulls the
    // window out and uses it for actor management and EndEvent
    // observation. parent is the Qt parent (ReaderMainWindow) so
    // destruction order is determined.
    explicit MoleculeScene(QVTKOpenGLNativeWidget* vtkWidget,
                           vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
                           QObject* parent = nullptr);
    ~MoleculeScene() override;

    // Build the vtkMolecule from the protein topology + frame-0 positions.
    // Must be called once after construction, before any setFrame() call.
    // Idempotent on the same protein/conformation — subsequent calls with
    // the same pointers no-op; different pointers rebuild from scratch.
    void Build(const model::QtProtein& protein,
               model::Conformation&    conformation);

    // Current renderer — for future overlay classes to attach actors.
    vtkRenderer* Renderer() const { return renderer_.Get(); }

    // Overlay-layer renderer — at SetLayer(1), shares the active camera
    // with the main renderer. Markers (MeasurementOverlay,
    // SceneRevealOverlay) paint here so they remain visible regardless
    // of depth occlusion. See spec §5.
    vtkRenderer* OverlayRenderer() const { return overlayRenderer_.Get(); }

    // Camera composer — owns the per-frame absolute camera write
    // dispatched on CameraMode + OrientationPolicy. Public so the REST
    // surface and CameraInputFilter can drive it. Lifetime tied to the
    // scene (Qt parent).
    CameraComposer* cameraComposer() const { return composer_; }

    // Reset camera to frame the molecule. Call after Build().
    void ResetCamera();

    // Overlays owned by the scene. MoleculeScene propagates setFrame
    // to each and issues a single Render() at the end. Nullable before
    // Build(); non-null after. Raw pointers because lifetime is tied to
    // the scene via QObject parent.
    QtBackboneRibbonOverlay* ribbonOverlay()     const { return ribbon_; }
    QtRingPolygonOverlay*    ringPolygonOverlay() const { return ringPolygons_; }
    QtFieldGridOverlay*      fieldGridOverlay()  const { return fieldGrid_; }
    QtBFieldStreamOverlay*   bfieldStreamOverlay() const { return bfieldStream_; }
    QtSelectionOverlay*      selectionOverlay()   const { return selection_; }
    MeasurementOverlay*      measurementOverlay() const { return measurement_; }
    SceneRevealOverlay*      revealOverlay()      const { return reveal_; }

public slots:
    // Update atom positions to frame t AND propagate to every overlay.
    // Early-returns when t == currentFrame to keep playback cheap; use
    // refreshCurrentFrame() to force a re-update on the same frame
    // (after toggling a visibility or changing an overlay parameter).
    void setFrame(int t);

    // Schedule one render via the Qt paint chain (widget->update()).
    // Coalesces multiple per-tick requests into one paint. The source
    // tag is recorded for the EndEvent observer to log.
    void requestRender(RenderSource src);

    // Back-compat shim — overlay toggle paths still call requestRender()
    // without a source. Forwards to requestRender(External).
    void requestRender();

    // Re-run setFrame() on the current frame, ignoring the early-return
    // guard. Use this when an overlay that skips expensive work while
    // hidden (QtFieldGridOverlay, future B-field overlay) is turned
    // back on — its kernel re-eval needs to run for the current frame.
    void refreshCurrentFrame();

    // Dashboard strip reveal: highlight and camera-focus the atom, residue, or
    // atom tuple represented by a strip binding without changing AtomSelection.
    void revealBinding(const model::SignalBinding& binding);
    void clearReveal();

    // ---- Plane lock compatibility shim ---------------------------------
    //
    // Public signatures unchanged for back-compat with the toolbar action
    // and REST surface (/plane-lock/*). The implementation now delegates
    // to the camera composer (CameraMode::Plane + OrientationPolicy::Default).
    // The internal CameraPlaneLock struct has been removed; the composer
    // owns the lock state.
    bool lockCameraToSelectionPlane(const std::vector<std::size_t>& atoms);
    void clearCameraPlaneLock();
    bool isCameraPlaneLocked() const;
    std::vector<std::size_t> cameraPlaneLockAtoms() const;

signals:
    void cameraPlaneLockChanged(bool active);

private:
    void focusCameraOnReveal(const model::SignalBinding& binding,
                             const std::vector<std::size_t>& atoms,
                             int frame);

    // Per-frame helper: push atom positions into vtkMolecule and
    // accumulate the bounds in one pass. Bounds are passed out to the
    // caller (setFrame) which uses them to set the renderer's clipping
    // range. See feedback_vtk_bounds_cache for why we compute bounds
    // ourselves rather than calling vtkActor::GetBounds().
    void PushAtomPositions(int t, double bounds[6]);

    QPointer<QVTKOpenGLNativeWidget>              vtkWidget_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkRenderer>                  overlayRenderer_;
    vtkSmartPointer<vtkMolecule>                  molecule_;
    vtkSmartPointer<vtkOpenGLMoleculeMapper>      mapper_;
    vtkSmartPointer<vtkActor>                     actor_;

    QtBackboneRibbonOverlay* ribbon_       = nullptr;   // QObject child
    QtRingPolygonOverlay*    ringPolygons_ = nullptr;   // QObject child
    QtFieldGridOverlay*      fieldGrid_    = nullptr;   // QObject child
    QtBFieldStreamOverlay*   bfieldStream_ = nullptr;   // QObject child
    QtSelectionOverlay*      selection_    = nullptr;   // QObject child (dormant; superseded by measurement_)
    MeasurementOverlay*      measurement_  = nullptr;   // QObject child
    SceneRevealOverlay*      reveal_       = nullptr;   // QObject child
    CameraComposer*          composer_     = nullptr;   // QObject child

    const model::QtProtein*       protein_      = nullptr;
    QPointer<model::Conformation> conformation_;
    int                           currentFrame_ = -1;
    std::optional<model::SignalBinding> activeRevealBinding_;

    // Render scheduler state. Single-pending-paint coalesce per
    // event-loop tick. lastRenderSource_ is set by requestRender and
    // read by the EndEvent observer in the constructor. Plain non-atomic
    // members because the scene is GUI-thread-only (ASSERT_THREAD enforced).
    bool         renderPending_    = false;
    RenderSource lastRenderSource_ = RenderSource::External;
};

}  // namespace h5reader::app

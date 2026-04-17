// MoleculeScene — VTK render pipeline for one protein + animated positions.
//
// Owns the vtkMolecule, its mapper, and the actor. Attaches the mapper
// via SetInputData (molecule is a vtkDataObject, not a pipeline source).
// On setFrame(t), updates every atom's position in place via
// vtkMolecule::SetAtomPosition + Modified() — the idiomatic VTK pattern
// for trajectory animation. Bond connectivity is static; lengths follow
// atom positions automatically because the mapper re-queries them.
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
//    AFTER molecule positions update and BEFORE renderWindow_->Render().
//    The overlay must update its VTK backing data (vtkPoints,
//    vtkImageData scalars, vtkStructuredGrid vectors, …) and call
//    Modified() as needed. The overlay MUST NOT call Render() itself
//    — MoleculeScene issues exactly one Render per setFrame.
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

#include "../model/QtConformation.h"
#include "../model/QtProtein.h"

#include <QObject>

#include <vtkActor.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkMolecule.h>
#include <vtkOpenGLMoleculeMapper.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <memory>

namespace h5reader::app {

class QtBackboneRibbonOverlay;
class QtRingPolygonOverlay;
class QtFieldGridOverlay;
class QtBFieldStreamOverlay;
class QtSelectionOverlay;

class MoleculeScene final : public QObject {
    Q_OBJECT

public:
    explicit MoleculeScene(vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
                           QObject* parent = nullptr);
    ~MoleculeScene() override;

    // Build the vtkMolecule from the protein topology + frame-0 positions.
    // Must be called once after construction, before any setFrame() call.
    // Idempotent on the same protein/conformation — subsequent calls with
    // the same pointers no-op; different pointers rebuild from scratch.
    void Build(const model::QtProtein&      protein,
               const model::QtConformation& conformation);

    // Current renderer — for future overlay classes to attach actors.
    vtkRenderer* Renderer() const { return renderer_.Get(); }

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

public slots:
    // Update atom positions to frame t AND propagate to every overlay.
    // Early-returns when t == currentFrame to keep playback cheap; use
    // refreshCurrentFrame() to force a re-update on the same frame
    // (after toggling a visibility or changing an overlay parameter).
    void setFrame(int t);

    // Render the current scene without running any overlay updates.
    // Intended for visibility / opacity / threshold toggles where
    // actor state changed but the underlying data is still valid.
    void requestRender();

    // Re-run setFrame() on the current frame, ignoring the early-return
    // guard. Use this when an overlay that skips expensive work while
    // hidden (QtFieldGridOverlay, future B-field overlay) is turned
    // back on — its kernel re-eval needs to run for the current frame.
    void refreshCurrentFrame();

private:
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkMolecule>                  molecule_;
    vtkSmartPointer<vtkOpenGLMoleculeMapper>      mapper_;
    vtkSmartPointer<vtkActor>                     actor_;

    QtBackboneRibbonOverlay* ribbon_       = nullptr;   // QObject child
    QtRingPolygonOverlay*    ringPolygons_ = nullptr;   // QObject child
    QtFieldGridOverlay*      fieldGrid_    = nullptr;   // QObject child
    QtBFieldStreamOverlay*   bfieldStream_ = nullptr;   // QObject child
    QtSelectionOverlay*      selection_    = nullptr;   // QObject child

    // Computes the centroid of the given frame's atom positions.
    // Used for the camera-follow feature so the molecule stays centred
    // in the view as it diffuses through the MD simulation box.
    model::Vec3 ComputeCentroid(size_t tIndex) const;

    const model::QtProtein*      protein_      = nullptr;
    const model::QtConformation* conformation_ = nullptr;
    int                          currentFrame_ = -1;

    // Last centroid the camera was pointed at. Initialised at Build()
    // time from frame 0; every setFrame() translates the camera by the
    // delta between the new centroid and this one, then updates it.
    // Keeps view direction and zoom stable while following translation.
    model::Vec3 lastCentroid_ = model::Vec3::Zero();
    bool        haveLastCentroid_ = false;
};

}  // namespace h5reader::app

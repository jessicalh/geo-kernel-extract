// QtOccupancyShellsOverlay — per-atom occupation-probability envelope shells.
//
// For the FOCUSED atom of the AtomSelection, gather its (already backbone-
// aligned) position across every trajectory frame, estimate the occupation
// density (OccupancyShellsMath), and draw two translucent nested isosurfaces:
// the 50% and 90% highest-density regions — the conformational "shadow" the
// still structure scrubs through. See notes/OCCUPANCY_SHELLS_SPEC_2026-06-13.md.
//
// MoleculeScene-owned; obeys the overlay contract (MoleculeScene.h §1-5) with
// ONE deliberate departure: the shells are a whole-trajectory AGGREGATE, so they
// are frame-INVARIANT — this overlay is NOT in the scene's setFrame fan. It
// rebuilds instead on (a) AtomSelection::focusChanged, (b)
// TransformedConformation::transformChanged (a fit-mode change moves every
// aligned position), and (c) becoming visible while dirty. All VTK mutation on
// the GUI thread; the overlay never calls Render() — the scene schedules it.
//
// Pipeline mirrors QtFieldGridOverlay: vtkImageData (point scalars) ->
// vtkTrivialProducer -> vtkContourFilter (x2, one per shell) -> translucent
// actors on the MAIN renderer (layer 0), so the shells depth-compose with the
// molecule. The math is hand-rolled in OccupancyShellsMath only where VTK has
// no fitting primitive (the anisotropic-covariance kernel); the contour itself
// is VTK's, because that IS the right tool.

#pragma once

#include "../model/Conformation.h"
#include "../model/QtProtein.h"
#include "OccupancyShellsMath.h"

#include <QObject>
#include <QPointer>

#include <vtkActor.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkImageData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>

#include <array>
#include <cstddef>

namespace h5reader::model {
class AtomSelection;
}

namespace h5reader::app {

class QtOccupancyShellsOverlay final : public QObject {
    Q_OBJECT

public:
    explicit QtOccupancyShellsOverlay(
        vtkSmartPointer<vtkRenderer>                  renderer,
        vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
        QObject*                                      parent = nullptr);
    ~QtOccupancyShellsOverlay() override;

    void Build(const model::QtProtein& protein,
               model::Conformation&    conformation);

    // Bind the selection whose FOCUS atom the shells follow. Wired by
    // ReaderMainWindow once the scene + selection both exist (the scene is
    // built before the selection, so this cannot happen in Build()).
    void setSelection(model::AtomSelection* selection);

    // World-space bbox of the current shells, for MoleculeScene to union into
    // the camera clipping range (the shells can extend well past the
    // current-frame molecule bounds). Returns false when no shells are shown.
    bool worldBounds(double out[6]) const;

public slots:
    // Focus atom changed -> rebuild for it (or mark dirty if hidden).
    void onFocusChanged(std::size_t atomIdx);
    // Selection emptied -> hide.
    void onSelectionCleared();
    // Alignment transform changed (fit-mode / smoothing window) -> every
    // aligned position moved, so the aggregate is stale; rebuild (or dirty).
    void onTransformChanged();

    // Master visibility. On enable, rebuilds for the current focus (the
    // dirty-state covers focus/transform changes that happened while hidden).
    void setVisible(bool on);

private:
    void rebuild();              // rebuild shells for the current focus atom
    void hideShells();
    void applyActorStyling();

    struct Shell {
        vtkSmartPointer<vtkContourFilter> contour;
        vtkSmartPointer<vtkActor>         actor;
        double                            fraction = 0.0;     // 0.5 / 0.9
        double                            opacity  = 0.3;
        std::array<double, 3>             color    = {0.55, 0.40, 0.80};
    };

    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    vtkSmartPointer<vtkImageData>                 imageData_;
    vtkSmartPointer<vtkDoubleArray>               scalars_;
    vtkSmartPointer<vtkTrivialProducer>           producer_;
    std::array<Shell, 2>                          shells_;   // [0]=50% inner, [1]=90% outer

    const model::QtProtein*        protein_ = nullptr;
    QPointer<model::Conformation>  conformation_;
    QPointer<model::AtomSelection> selection_;

    math::OccupancyConfig cfg_;
    bool   visible_   = false;   // off by default — user enables
    bool   dirty_     = false;   // focus/transform changed while hidden
    bool   hasShells_ = false;   // shells currently built + visible
    double worldBounds_[6] = {0, 0, 0, 0, 0, 0};
};

}  // namespace h5reader::app

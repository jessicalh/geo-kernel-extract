// QtRingPolygonOverlay — translucent filled polygon + normal arrow per
// aromatic ring. Ring geometry (center, normal, radius) comes from the
// H5's ring_geometry/data slab per frame; on frameChanged the overlay
// recomputes polygon vertices and re-renders.
//
// Adapts the pattern from ui/src/RingCurrentOverlay.cpp for trajectory
// animation. Colour is driven by the parent residue's amino-acid type
// via typed dispatch (AminoAcid enum), never by string matching.

#pragma once

#include "../model/QtConformation.h"
#include "../model/QtProtein.h"
#include "../model/QtRing.h"

#include <QObject>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <vector>

namespace h5reader::app {

class QtRingPolygonOverlay final : public QObject {
    Q_OBJECT

public:
    // renderer and renderWindow must already be set up by MoleculeScene.
    explicit QtRingPolygonOverlay(
        vtkSmartPointer<vtkRenderer>                  renderer,
        vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
        QObject* parent = nullptr);
    ~QtRingPolygonOverlay() override;

    // Build one polygon-and-arrow actor pair per ring. Must be called
    // before setFrame(). Idempotent on the same inputs.
    void Build(const model::QtProtein&      protein,
               const model::QtConformation& conformation);

public slots:
    // Update polygon vertices and arrow orientation from ring_geometry
    // at frame t. Called from QtPlaybackController::frameChanged by
    // MoleculeScene's signal relay.
    void setFrame(int t);

    // Show/hide every actor owned by this overlay.
    void setVisible(bool visible);

private:
    struct RingActor {
        // Two translucent polygons offset ±0.15 Å along the ring normal
        // for visual depth; update one vtkPoints each per frame.
        vtkSmartPointer<vtkPoints>   pointsAbove;
        vtkSmartPointer<vtkPoints>   pointsBelow;
        vtkSmartPointer<vtkPolyData> polyAbove;
        vtkSmartPointer<vtkPolyData> polyBelow;
        vtkSmartPointer<vtkActor>    actorAbove;
        vtkSmartPointer<vtkActor>    actorBelow;
        // Normal-direction arrow (thin cone via translate-scale)
        vtkSmartPointer<vtkActor>    arrowActor;
        // Polygon-side color, ring-type driven.
        unsigned char rgb[3] = {200, 200, 100};
    };

    void UpdateRingActor(RingActor& ra,
                         const model::RingGeometry& geo);

    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    const model::QtProtein*                       protein_      = nullptr;
    const model::QtConformation*                  conformation_ = nullptr;
    std::vector<RingActor>                        rings_;
    bool                                          visible_      = true;
};

}  // namespace h5reader::app

// QtBFieldStreamOverlay — B-field streamlines around each aromatic ring.
//
// The iconic NMR "butterfly" picture is a pair of complementary
// visualisations: scalar isosurfaces (T0 shielding, handled by
// QtFieldGridOverlay) AND vector streamlines of the B-field. This
// overlay renders the latter by tracing streamlines through a
// per-ring B-field grid, tubed via vtkTubeFilter and coloured by
// |B| with a Moreland diverging colormap.
//
// Per-frame update: on setFrame(t) we recompute the grid points
// (ring-local frame centred on the current frame's ring centre) and
// the B-field vectors via calculators::EvaluateBField, Modified() the
// structured grid, and vtkStreamTracer reruns on next Render.
//
// Ported from ui/src/ButterflyOverlay.cpp. Same stream tracer
// parameters; per-frame kernel evaluation added.

#pragma once

#include "../model/QtConformation.h"
#include "../model/QtProtein.h"

#include <QObject>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkStreamTracer.h>
#include <vtkStructuredGrid.h>
#include <vtkTubeFilter.h>

#include <vector>

namespace h5reader::app {

class QtBFieldStreamOverlay final : public QObject {
    Q_OBJECT

public:
    explicit QtBFieldStreamOverlay(
        vtkSmartPointer<vtkRenderer>                  renderer,
        vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
        QObject* parent = nullptr);
    ~QtBFieldStreamOverlay() override;

    void Build(const model::QtProtein&      protein,
               const model::QtConformation& conformation);

public slots:
    void setFrame(int t);
    void setVisible(bool visible);

private:
    struct RingStream {
        vtkSmartPointer<vtkPoints>         gridPoints;
        vtkSmartPointer<vtkDoubleArray>    vectors;      // 3-component B (Tesla)
        vtkSmartPointer<vtkDoubleArray>    magnitudes;   // |B| (Tesla)
        vtkSmartPointer<vtkStructuredGrid> grid;

        vtkSmartPointer<vtkPoints>         seedPoints;
        vtkSmartPointer<vtkCellArray>      seedVerts;
        vtkSmartPointer<vtkPolyData>       seedPoly;

        vtkSmartPointer<vtkStreamTracer>   tracer;
        vtkSmartPointer<vtkTubeFilter>     tubes;
        vtkSmartPointer<vtkActor>          actor;
    };

    void UpdateRing(size_t ri, int t);

    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    const model::QtProtein*                       protein_      = nullptr;
    const model::QtConformation*                  conformation_ = nullptr;
    std::vector<RingStream>                       rings_;
    bool                                          visible_      = false;
};

}  // namespace h5reader::app

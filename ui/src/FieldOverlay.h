#pragma once
#include "ComputeWorker.h"
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkColorTransferFunction.h>
#include <vector>

// Visualizes per-atom scalar fields as colored spheres and E-field as arrows
//
// Scalar fields: atom spheres colored by a diverging blue-white-red colormap
// E-field arrows: orange arrows at each atom pointing in the E-field direction

class FieldOverlay {
public:
    explicit FieldOverlay(vtkSmartPointer<vtkRenderer> renderer);
    ~FieldOverlay();

    // Show colored spheres at atom positions by scalar value
    void setScalarField(const std::vector<nmr::Vec3>& positions,
                        const std::vector<double>& values,
                        double sphereRadius = 0.25,
                        double opacity = 0.7);

    // Show arrow glyphs for E-field vectors
    void setArrows(const std::vector<nmr::Vec3>& positions,
                   const std::vector<nmr::Vec3>& vectors,
                   double scale = 2.0,
                   double opacity = 0.6);

    void clearScalar();
    void clearArrows();
    void clear();

    void setScalarVisible(bool visible);
    void setArrowsVisible(bool visible);

private:
    vtkSmartPointer<vtkRenderer> renderer_;
    vtkSmartPointer<vtkActor> scalarActor_;
    vtkSmartPointer<vtkActor> arrowActor_;
    vtkSmartPointer<vtkActor2D> scalarBarActor_;

};

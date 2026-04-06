#pragma once
#include "ComputeWorker.h"
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vector>

// Smooth isosurface overlay using Gaussian splatting + marching cubes.
// Produces translucent shielded (blue) / deshielded (red) surfaces
// that wrap around molecular regions instead of discrete atom spheres.

class IsosurfaceOverlay {
public:
    explicit IsosurfaceOverlay(vtkSmartPointer<vtkRenderer> renderer);
    ~IsosurfaceOverlay();

    void setData(const std::vector<nmr::Vec3>& positions,
                 const std::vector<double>& values,
                 double threshold = 0.5,
                 double gaussianRadius = 2.0,
                 double opacity = 0.4);

    void clear();
    void setVisible(bool visible);

private:
    vtkSmartPointer<vtkRenderer> renderer_;
    vtkSmartPointer<vtkActor> posActor_;   // deshielded (positive) surface
    vtkSmartPointer<vtkActor> negActor_;   // shielded (negative) surface
    vtkSmartPointer<vtkActor2D> barActor_;
};

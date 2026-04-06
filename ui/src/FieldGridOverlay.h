#pragma once
//
// FieldGridOverlay: renders isosurfaces from structured T0 field grids.
//
// Unlike the old IsosurfaceOverlay (which uses Gaussian splatting of point
// values), this renders physically accurate isosurfaces from the actual
// calculator field evaluated on a 3D grid.
//
// Two surfaces per ring:
//   Positive (red, translucent): T0 > +threshold (deshielded region)
//   Negative (blue, translucent): T0 < -threshold (shielded region)
//

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vector>
#include "ComputeWorker.h"  // ViewerFieldGrid

class FieldGridOverlay {
public:
    explicit FieldGridOverlay(vtkSmartPointer<vtkRenderer> renderer);
    ~FieldGridOverlay();

    // Build isosurfaces from all field grids.
    // mode: 0 = total T0, 1 = BS only
    void setData(const std::vector<ViewerFieldGrid>& grids,
                 double threshold, double opacity, int mode = 0);
    void clear();
    void setVisible(bool visible);
    void setThreshold(double threshold);

private:
    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkActor>> actors_;
    std::vector<ViewerFieldGrid> grids_;  // cached for threshold updates
    double opacity_ = 0.4;
    int mode_ = 0;
};

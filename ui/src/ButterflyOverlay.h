#pragma once
#include "ComputeWorker.h"
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vector>

// Visualizes Biot-Savart magnetic field streamlines around aromatic rings.
// Produces the classic "butterfly" pattern: field lines loop above and below
// the ring plane, diverging from the deshielding cone and converging through
// the shielding region above/below the ring.
//
// Streamlines are seeded on a circle at 1.5x ring radius in the ring plane
// and integrated in both directions using RK4-5.  Tubes are colored by
// field magnitude on a Moreland diverging (blue-white-red) colormap.

class ButterflyOverlay {
public:
    explicit ButterflyOverlay(vtkSmartPointer<vtkRenderer> renderer);
    ~ButterflyOverlay();

    void setData(const std::vector<ViewerButterflyData>& butterflies);
    void clear();
    void setVisible(bool visible);

private:
    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkActor>> actors_;
};

#pragma once
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vector>

namespace nmr {
class Protein;
class ProteinConformation;
}

// Visualizes aromatic rings as colored tube outlines in the 3D scene
// PHE=green, TYR=cyan, TRP=magenta, HIS=yellow

class RingCurrentOverlay {
public:
    RingCurrentOverlay(vtkSmartPointer<vtkRenderer> renderer,
                       const nmr::Protein& protein,
                       const nmr::ProteinConformation& conf);
    ~RingCurrentOverlay();

    void setVisible(bool visible);
    void setTubeRadius(double radius);

private:
    void buildActors(const nmr::Protein& protein,
                     const nmr::ProteinConformation& conf);

    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkActor>> actors_;
    double tubeRadius_ = 0.15;
};

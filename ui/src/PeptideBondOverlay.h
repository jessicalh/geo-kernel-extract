#pragma once
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vector>

namespace nmr {
class Protein;
class ProteinConformation;
}

// Visualizes bonds as colored tubes in the 3D scene
// Peptide bonds: red, Sidechain bonds: orange, Other: blue

class PeptideBondOverlay {
public:
    PeptideBondOverlay(vtkSmartPointer<vtkRenderer> renderer,
                       const nmr::Protein& protein,
                       const nmr::ProteinConformation& conf);
    ~PeptideBondOverlay();

    void setVisible(bool visible);
    void setShowSidechain(bool show);

private:
    void buildActors(const nmr::Protein& protein,
                     const nmr::ProteinConformation& conf);

    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkActor>> bbActors_;   // backbone / peptide
    std::vector<vtkSmartPointer<vtkActor>> scActors_;   // sidechain
    double tubeRadius_ = 0.06;
};

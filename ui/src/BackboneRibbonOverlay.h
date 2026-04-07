#pragma once
//
// BackboneRibbonOverlay: protein backbone ribbon colored by residue type.
//
// Uses vtkProteinRibbonFilter with arrays built from the library's Protein
// object model. No vtkPDBReader round-trip needed — we construct the
// required point data arrays directly from Protein/ProteinConformation.
//

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>

namespace nmr {
class Protein;
class ProteinConformation;
}

class BackboneRibbonOverlay {
public:
    explicit BackboneRibbonOverlay(vtkSmartPointer<vtkRenderer> renderer,
                                    const nmr::Protein& protein,
                                    const nmr::ProteinConformation& conf);
    ~BackboneRibbonOverlay();

    void setVisible(bool visible);
    bool isVisible() const;

private:
    vtkSmartPointer<vtkRenderer> renderer_;
    vtkSmartPointer<vtkActor> actor_;
};

// Colormap.h — perceptually uniform diverging colormap helper.
//
// Moreland cool-to-warm with CIELAB interpolation via
// vtkColorTransferFunction::SetColorSpaceToDiverging(). Symmetric
// about zero: cool blue (negative) → neutral grey → warm red (positive).
//
// Copied verbatim from ui/src/Colormap.h — the library-viewer reuse we
// need for the B-field streamline overlay. Header-only, no external
// dependencies beyond VTK.

#pragma once

#include <vtkColorTransferFunction.h>
#include <vtkSmartPointer.h>

#include <algorithm>
#include <cmath>

namespace h5reader::app {

inline vtkSmartPointer<vtkColorTransferFunction> MakeDivergingCTF(
    double minVal, double maxVal) {
    double absMax = std::max(std::abs(minVal), std::abs(maxVal));
    if (absMax < 1e-10) absMax = 1.0;

    auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->SetColorSpaceToDiverging();
    // Moreland cool-to-warm endpoints.
    ctf->AddRGBPoint(-absMax, 0.230, 0.299, 0.754);   // cool blue
    ctf->AddRGBPoint( 0.0,    0.865, 0.865, 0.865);   // neutral grey
    ctf->AddRGBPoint( absMax, 0.706, 0.016, 0.150);   // warm red
    ctf->Build();
    return ctf;
}

}  // namespace h5reader::app

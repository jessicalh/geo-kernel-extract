#pragma once
#include <vtkSmartPointer.h>
#include <vtkColorTransferFunction.h>
#include <cmath>
#include <algorithm>

// Perceptually uniform diverging colormap (Moreland cool-to-warm)
// Uses CIELAB interpolation via vtkColorTransferFunction::SetColorSpaceToDiverging()
// Symmetric about zero: blue (shielded) - white - red (deshielded)

inline vtkSmartPointer<vtkColorTransferFunction> makeDivergingCTF(double minVal, double maxVal) {
    double absMax = std::max(std::abs(minVal), std::abs(maxVal));
    if (absMax < 1e-10) absMax = 1.0;

    auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->SetColorSpaceToDiverging();
    // Moreland cool-to-warm endpoints
    ctf->AddRGBPoint(-absMax, 0.230, 0.299, 0.754);  // cool blue
    ctf->AddRGBPoint( 0.0,   0.865, 0.865, 0.865);   // neutral gray-white
    ctf->AddRGBPoint( absMax, 0.706, 0.016, 0.150);   // warm red
    ctf->Build();
    return ctf;
}

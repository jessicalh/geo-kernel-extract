#pragma once
#include "ComputeWorker.h"
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vector>

// Renders spherical harmonic surface glyphs for T2 tensor visualization
// f(theta,phi) = sum_m T2_m * Y_2^m(theta,phi)
// Radius = |f|, color: red = deshielded (f>0), blue = shielded (f<0)

class TensorGlyph {
public:
    explicit TensorGlyph(vtkSmartPointer<vtkRenderer> renderer);
    ~TensorGlyph();

    // Set tensor data: positions and spherical tensors
    void setData(const std::vector<nmr::Vec3>& positions,
                 const std::vector<nmr::SphericalTensor>& tensors,
                 double scale = 0.5,
                 double opacity = 0.7);

    void clear();
    void setVisible(bool visible);

private:
    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkActor>> actors_;

    // Generate a single SH glyph polydata at the origin
    vtkSmartPointer<vtkActor> createGlyph(
        const nmr::Vec3& center,
        const nmr::SphericalTensor& st,
        double scale,
        double opacity);

    static constexpr int THETA_RES = 20;
    static constexpr int PHI_RES = 40;
};

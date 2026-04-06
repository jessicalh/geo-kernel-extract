#pragma once
#include "ComputeWorker.h"
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vector>

// Batched ellipsoid tensor glyphs using VTK's vtkTensorGlyph pipeline.
// ONE actor for all glyphs — dramatically faster than TensorGlyph's per-atom actors.
// Ellipsoid axes show tensor anisotropy; color shows isotropic value.

class EllipsoidGlyph {
public:
    explicit EllipsoidGlyph(vtkSmartPointer<vtkRenderer> renderer);
    ~EllipsoidGlyph();

    void setData(const std::vector<nmr::Vec3>& positions,
                 const std::vector<nmr::Mat3>& tensors,
                 const std::vector<double>& isoValues,
                 double scale = 0.5,
                 double opacity = 0.7);

    void clear();
    void setVisible(bool visible);

private:
    vtkSmartPointer<vtkRenderer> renderer_;
    vtkSmartPointer<vtkActor> actor_;
    vtkSmartPointer<vtkActor2D> barActor_;
};

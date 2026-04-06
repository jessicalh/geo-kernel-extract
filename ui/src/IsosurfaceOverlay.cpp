#include "IsosurfaceOverlay.h"
#include "Colormap.h"
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkGaussianSplatter.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkTrivialProducer.h>
#include <vtkNew.h>
#include <cmath>
#include <algorithm>

IsosurfaceOverlay::IsosurfaceOverlay(vtkSmartPointer<vtkRenderer> renderer)
    : renderer_(renderer)
{
}

IsosurfaceOverlay::~IsosurfaceOverlay() {
    clear();
}

void IsosurfaceOverlay::clear() {
    if (posActor_) { renderer_->RemoveActor(posActor_); posActor_ = nullptr; }
    if (negActor_) { renderer_->RemoveActor(negActor_); negActor_ = nullptr; }
    if (barActor_) { renderer_->RemoveActor2D(barActor_); barActor_ = nullptr; }
}

void IsosurfaceOverlay::setVisible(bool visible) {
    int vis = visible ? 1 : 0;
    if (posActor_) posActor_->SetVisibility(vis);
    if (negActor_) negActor_->SetVisibility(vis);
    if (barActor_) barActor_->SetVisibility(vis);
}

void IsosurfaceOverlay::setData(
    const std::vector<nmr::Vec3>& positions,
    const std::vector<double>& values,
    double threshold,
    double gaussianRadius,
    double opacity)
{
    clear();
    if (positions.empty()) return;

    // Build point cloud with scalars
    vtkNew<vtkPoints> points;
    vtkNew<vtkFloatArray> scalars;
    scalars->SetName("Values");
    scalars->SetNumberOfComponents(1);

    double bounds[6] = {1e30, -1e30, 1e30, -1e30, 1e30, -1e30};
    double vmin = 1e30, vmax = -1e30;

    for (size_t i = 0; i < positions.size(); ++i) {
        points->InsertNextPoint(positions[i].data());
        float val = static_cast<float>(values[i]);
        scalars->InsertNextValue(val);
        vmin = std::min(vmin, static_cast<double>(val));
        vmax = std::max(vmax, static_cast<double>(val));
        for (int d = 0; d < 3; ++d) {
            bounds[2*d]   = std::min(bounds[2*d],   positions[i][d]);
            bounds[2*d+1] = std::max(bounds[2*d+1], positions[i][d]);
        }
    }

    // Pad bounds by 2x gaussian radius
    double pad = gaussianRadius * 2.0;
    for (int d = 0; d < 3; ++d) {
        bounds[2*d]   -= pad;
        bounds[2*d+1] += pad;
    }

    vtkNew<vtkPolyData> pd;
    pd->SetPoints(points);
    pd->GetPointData()->SetScalars(scalars);

    vtkNew<vtkTrivialProducer> producer;
    producer->SetOutput(pd);

    // Compute radius as fraction of model bounds diagonal
    double dx = bounds[1] - bounds[0];
    double dy = bounds[3] - bounds[2];
    double dz = bounds[5] - bounds[4];
    double diagonal = std::sqrt(dx*dx + dy*dy + dz*dz);
    double radiusFraction = (diagonal > 1e-6) ? gaussianRadius / diagonal : 0.1;

    // Gaussian splatter — creates a 3D scalar field
    vtkNew<vtkGaussianSplatter> splatter;
    splatter->SetInputConnection(producer->GetOutputPort());
    splatter->SetModelBounds(bounds);
    splatter->SetSampleDimensions(64, 64, 64);
    splatter->SetRadius(radiusFraction);
    splatter->SetScalarWarping(true);
    splatter->SetExponentFactor(-4.0);
    splatter->SetNullValue(0.0);

    double absMax = std::max(std::abs(vmin), std::abs(vmax));
    if (absMax < 1e-10) return;

    auto ctf = makeDivergingCTF(-absMax, absMax);

    // Helper: create contour actor at a given isosurface level
    auto makeContourActor = [&](double level, double r, double g, double b) -> vtkSmartPointer<vtkActor> {
        vtkNew<vtkContourFilter> contour;
        contour->SetInputConnection(splatter->GetOutputPort());
        contour->SetValue(0, level);

        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(contour->GetOutputPort());
        mapper->ScalarVisibilityOff();

        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(r, g, b);
        actor->GetProperty()->SetOpacity(opacity);
        actor->GetProperty()->SetInterpolationToPhong();
        actor->SetForceTranslucent(true);
        return actor;
    };

    // Positive (deshielded) isosurface — warm red
    if (vmax > threshold * 0.1) {
        posActor_ = makeContourActor(threshold, 0.706, 0.016, 0.150);
        renderer_->AddActor(posActor_);
    }

    // Negative (shielded) isosurface — cool blue
    if (vmin < -threshold * 0.1) {
        negActor_ = makeContourActor(-threshold, 0.230, 0.299, 0.754);
        renderer_->AddActor(negActor_);
    }

    // Color bar
    auto bar = vtkSmartPointer<vtkScalarBarActor>::New();
    bar->SetLookupTable(ctf);
    bar->SetNumberOfLabels(5);
    bar->SetWidth(0.08);
    bar->SetHeight(0.4);
    bar->SetPosition(0.02, 0.3);
    bar->GetLabelTextProperty()->SetFontSize(10);
    bar->GetLabelTextProperty()->SetColor(0.9, 0.9, 0.9);
    bar->GetTitleTextProperty()->SetColor(0.9, 0.9, 0.9);
    bar->SetTitle("ppm");
    barActor_ = bar;
    renderer_->AddActor2D(bar);
}

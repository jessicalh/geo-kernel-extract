#include "FieldGridOverlay.h"
#include "OperationLog.h"
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkTrivialProducer.h>
#include <vtkNew.h>
#include <cmath>

using namespace nmr;

FieldGridOverlay::FieldGridOverlay(vtkSmartPointer<vtkRenderer> renderer)
    : renderer_(renderer) {}

FieldGridOverlay::~FieldGridOverlay() { clear(); }

void FieldGridOverlay::clear() {
    for (auto& actor : shieldedActors_)
        renderer_->RemoveActor(actor);
    for (auto& actor : deshieldedActors_)
        renderer_->RemoveActor(actor);
    shieldedActors_.clear();
    deshieldedActors_.clear();
}

void FieldGridOverlay::setVisible(bool visible) {
    setShieldedVisible(visible);
    setDeshieldedVisible(visible);
}

void FieldGridOverlay::setShieldedVisible(bool visible) {
    for (auto& actor : shieldedActors_)
        actor->SetVisibility(visible ? 1 : 0);
}

void FieldGridOverlay::setDeshieldedVisible(bool visible) {
    for (auto& actor : deshieldedActors_)
        actor->SetVisibility(visible ? 1 : 0);
}

void FieldGridOverlay::setData(const std::vector<ViewerFieldGrid>& grids,
                                double threshold, double opacity, int mode) {
    clear();
    grids_ = grids;
    opacity_ = opacity;
    mode_ = mode;

    for (size_t gi = 0; gi < grids.size(); gi++) {
        const auto& grid = grids[gi];
        int nx = grid.dims[0], ny = grid.dims[1], nz = grid.dims[2];
        int nPoints = nx * ny * nz;

        OperationLog::Info(LogViewer, "FieldGridOverlay",
            "Grid " + std::to_string(gi) + " (" + grid.ringType +
            "): " + std::to_string(nx) + "x" + std::to_string(ny) +
            "x" + std::to_string(nz) + " = " + std::to_string(nPoints) + " points");

        if (nPoints == 0) continue;

        const std::vector<double>& values = (mode == 1) ? grid.bsT0 : grid.T0;
        if (values.empty()) continue;

        double vmin = 1e30, vmax = -1e30;
        for (double v : values) {
            if (v < vmin) vmin = v;
            if (v > vmax) vmax = v;
        }

        // Build VTK structured grid (vtkImageData)
        vtkNew<vtkImageData> imageData;
        imageData->SetDimensions(nx, ny, nz);
        imageData->SetOrigin(grid.origin[0], grid.origin[1], grid.origin[2]);
        imageData->SetSpacing(grid.spacing[0], grid.spacing[1], grid.spacing[2]);

        vtkNew<vtkFloatArray> scalars;
        scalars->SetName("T0");
        scalars->SetNumberOfTuples(nPoints);
        for (int i = 0; i < nPoints; i++)
            scalars->SetValue(i, static_cast<float>(values[i]));
        imageData->GetPointData()->SetScalars(scalars);

        vtkNew<vtkTrivialProducer> producer;
        producer->SetOutput(imageData);

        // Shielded isosurface (T0 < -threshold, above/below ring — sky blue)
        if (vmin < -threshold) {
            vtkNew<vtkContourFilter> contour;
            contour->SetInputConnection(producer->GetOutputPort());
            contour->SetValue(0, -threshold);
            contour->Update();

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(contour->GetOutputPort());
            mapper->ScalarVisibilityOff();

            auto actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->GetProperty()->SetColor(0.50, 0.70, 0.95);  // pastel sky blue
            actor->GetProperty()->SetOpacity(opacity);
            actor->GetProperty()->SetInterpolationToPhong();
            actor->SetForceTranslucent(true);
            renderer_->AddActor(actor);
            shieldedActors_.push_back(actor);
        }

        // Deshielded isosurface (T0 > +threshold, in ring plane — coral)
        if (vmax > threshold) {
            vtkNew<vtkContourFilter> contour;
            contour->SetInputConnection(producer->GetOutputPort());
            contour->SetValue(0, threshold);
            contour->Update();

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(contour->GetOutputPort());
            mapper->ScalarVisibilityOff();

            auto actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            actor->GetProperty()->SetColor(0.95, 0.55, 0.45);  // pastel coral
            actor->GetProperty()->SetOpacity(opacity);
            actor->GetProperty()->SetInterpolationToPhong();
            actor->SetForceTranslucent(true);
            renderer_->AddActor(actor);
            deshieldedActors_.push_back(actor);
        }
    }
}

void FieldGridOverlay::setThreshold(double threshold) {
    setData(grids_, threshold, opacity_, mode_);
}

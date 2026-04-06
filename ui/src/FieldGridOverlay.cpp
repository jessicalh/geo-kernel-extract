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
    for (auto& actor : actors_)
        renderer_->RemoveActor(actor);
    actors_.clear();
}

void FieldGridOverlay::setVisible(bool visible) {
    for (auto& actor : actors_)
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
        if (values.empty()) {
            OperationLog::Info(LogViewer, "FieldGridOverlay", "Empty values, skipping");
            continue;
        }

        double vmin = 1e30, vmax = -1e30;
        int nNonZero = 0;
        for (double v : values) {
            if (v < vmin) vmin = v;
            if (v > vmax) vmax = v;
            if (std::abs(v) > 1e-10) nNonZero++;
        }

        OperationLog::Info(LogViewer, "FieldGridOverlay",
            "Range: [" + std::to_string(vmin) + ", " + std::to_string(vmax) +
            "], nonzero=" + std::to_string(nNonZero) +
            ", threshold=" + std::to_string(threshold));

        // Build VTK structured grid (vtkImageData)
        OperationLog::Info(LogViewer, "FieldGridOverlay", "Creating vtkImageData...");
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

        OperationLog::Info(LogViewer, "FieldGridOverlay", "vtkImageData created, running contours...");

        // Use vtkTrivialProducer to properly connect imageData into the pipeline
        // (SetInputData on pipeline algorithms can crash -- use SetInputConnection)
        vtkNew<vtkTrivialProducer> producer;
        producer->SetOutput(imageData);

        // Positive isosurface (deshielded — warm red)
        if (vmax > threshold) {
            OperationLog::Info(LogViewer, "FieldGridOverlay", "Creating positive contour...");
            vtkNew<vtkContourFilter> contourPos;
            contourPos->SetInputConnection(producer->GetOutputPort());
            contourPos->SetValue(0, threshold);
            contourPos->Update();

            vtkNew<vtkPolyDataMapper> mapperPos;
            mapperPos->SetInputConnection(contourPos->GetOutputPort());
            mapperPos->ScalarVisibilityOff();

            auto actorPos = vtkSmartPointer<vtkActor>::New();
            actorPos->SetMapper(mapperPos);
            actorPos->GetProperty()->SetColor(0.706, 0.016, 0.150);  // Moreland warm red
            actorPos->GetProperty()->SetOpacity(opacity);
            actorPos->GetProperty()->SetInterpolationToPhong();
            actorPos->SetForceTranslucent(true);
            renderer_->AddActor(actorPos);
            actors_.push_back(actorPos);
            OperationLog::Info(LogViewer, "FieldGridOverlay", "Positive contour added");
        }

        // Negative isosurface (shielded — cool blue)
        if (vmin < -threshold) {
            OperationLog::Info(LogViewer, "FieldGridOverlay", "Creating negative contour...");
            vtkNew<vtkContourFilter> contourNeg;
            contourNeg->SetInputConnection(producer->GetOutputPort());
            contourNeg->SetValue(0, -threshold);
            contourNeg->Update();

            vtkNew<vtkPolyDataMapper> mapperNeg;
            mapperNeg->SetInputConnection(contourNeg->GetOutputPort());
            mapperNeg->ScalarVisibilityOff();

            auto actorNeg = vtkSmartPointer<vtkActor>::New();
            actorNeg->SetMapper(mapperNeg);
            actorNeg->GetProperty()->SetColor(0.230, 0.299, 0.754);  // Moreland cool blue
            actorNeg->GetProperty()->SetOpacity(opacity);
            actorNeg->GetProperty()->SetInterpolationToPhong();
            actorNeg->SetForceTranslucent(true);
            renderer_->AddActor(actorNeg);
            actors_.push_back(actorNeg);
            OperationLog::Info(LogViewer, "FieldGridOverlay", "Negative contour added");
        }
    }
}

void FieldGridOverlay::setThreshold(double threshold) {
    setData(grids_, threshold, opacity_, mode_);
}

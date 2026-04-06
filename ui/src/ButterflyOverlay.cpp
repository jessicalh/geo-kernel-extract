#include "ButterflyOverlay.h"
#include "Colormap.h"
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStreamTracer.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkColorTransferFunction.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkNew.h>
#include <cmath>
#include <algorithm>

ButterflyOverlay::ButterflyOverlay(vtkSmartPointer<vtkRenderer> renderer)
    : renderer_(renderer)
{
}

ButterflyOverlay::~ButterflyOverlay() {
    clear();
}

void ButterflyOverlay::clear() {
    for (auto& actor : actors_)
        renderer_->RemoveActor(actor);
    actors_.clear();
}

void ButterflyOverlay::setVisible(bool visible) {
    int vis = visible ? 1 : 0;
    for (auto& actor : actors_)
        actor->SetVisibility(vis);
}

void ButterflyOverlay::setData(const std::vector<ViewerButterflyData>& butterflies) {
    clear();

    using nmr::Vec3;

    for (const auto& bf : butterflies) {
        if (bf.positions.empty() || bf.fields.empty()) continue;
        if (bf.positions.size() != bf.fields.size()) continue;

        const int nx = bf.gridDims[0];
        const int ny = bf.gridDims[1];
        const int nz = bf.gridDims[2];
        const int nTotal = nx * ny * nz;
        if (nTotal <= 0 || static_cast<int>(bf.positions.size()) != nTotal) continue;

        // ---- Build vtkStructuredGrid from grid data ----

        vtkNew<vtkPoints> gridPoints;
        gridPoints->SetNumberOfPoints(nTotal);

        vtkNew<vtkFloatArray> fieldVectors;
        fieldVectors->SetName("B-field");
        fieldVectors->SetNumberOfComponents(3);
        fieldVectors->SetNumberOfTuples(nTotal);

        vtkNew<vtkFloatArray> fieldMagnitude;
        fieldMagnitude->SetName("B-magnitude");
        fieldMagnitude->SetNumberOfComponents(1);
        fieldMagnitude->SetNumberOfTuples(nTotal);

        double magMin = 1e30, magMax = -1e30;

        for (int i = 0; i < nTotal; ++i) {
            gridPoints->SetPoint(i,
                bf.positions[i].x(), bf.positions[i].y(), bf.positions[i].z());

            float bx = static_cast<float>(bf.fields[i].x());
            float by = static_cast<float>(bf.fields[i].y());
            float bz = static_cast<float>(bf.fields[i].z());
            fieldVectors->SetTuple3(i, bx, by, bz);

            double mag = bf.fields[i].norm();
            fieldMagnitude->SetValue(i, static_cast<float>(mag));
            magMin = std::min(magMin, mag);
            magMax = std::max(magMax, mag);
        }

        auto grid = vtkSmartPointer<vtkStructuredGrid>::New();
        grid->SetDimensions(nx, ny, nz);
        grid->SetPoints(gridPoints);
        grid->GetPointData()->SetVectors(fieldVectors);
        grid->GetPointData()->AddArray(fieldMagnitude);

        // ---- Create seed points on a circle at 1.5x ring radius ----

        Vec3 n = bf.ringNormal.normalized();

        // Build orthonormal basis (u, v) in the ring plane
        Vec3 arbitrary = (std::abs(n.x()) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
        Vec3 u = n.cross(arbitrary).normalized();
        Vec3 v = n.cross(u);

        constexpr int nSeeds = 12;
        double seedRadius = 1.5 * bf.ringRadius;

        vtkNew<vtkPoints> seedPoints;
        vtkNew<vtkCellArray> seedVerts;
        for (int i = 0; i < nSeeds; ++i) {
            double theta = 2.0 * M_PI * i / nSeeds;
            Vec3 p = bf.ringCenter + seedRadius * (std::cos(theta) * u + std::sin(theta) * v);
            vtkIdType pid = seedPoints->InsertNextPoint(p.x(), p.y(), p.z());
            seedVerts->InsertNextCell(1, &pid);
        }

        auto seedPoly = vtkSmartPointer<vtkPolyData>::New();
        seedPoly->SetPoints(seedPoints);
        seedPoly->SetVerts(seedVerts);

        // ---- Stream tracer ----

        vtkNew<vtkStreamTracer> tracer;
        tracer->SetInputData(grid);
        tracer->SetSourceData(seedPoly);
        tracer->SetIntegrationDirectionToBoth();
        tracer->SetMaximumPropagation(12.0);
        tracer->SetIntegratorTypeToRungeKutta45();
        tracer->SetInitialIntegrationStep(0.1);
        tracer->SetMinimumIntegrationStep(0.01);
        tracer->SetMaximumIntegrationStep(0.5);
        tracer->SetMaximumNumberOfSteps(2000);
        tracer->SetInputArrayToProcess(
            0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "B-field");

        // ---- Tube filter ----

        vtkNew<vtkTubeFilter> tubes;
        tubes->SetInputConnection(tracer->GetOutputPort());
        tubes->SetRadius(0.04);
        tubes->SetNumberOfSides(8);
        tubes->CappingOn();

        // ---- Color by field magnitude (blue-white-red diverging) ----

        double absMax = std::max(std::abs(magMin), std::abs(magMax));
        if (absMax < 1e-10) absMax = 1.0;

        auto ctf = makeDivergingCTF(0.0, absMax);

        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(tubes->GetOutputPort());
        mapper->SetLookupTable(ctf);
        mapper->SetScalarRange(0.0, absMax);
        mapper->SetScalarModeToUsePointFieldData();
        mapper->SelectColorArray("B-magnitude");

        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetOpacity(0.85);
        actor->GetProperty()->SetInterpolationToPhong();
        actor->SetForceTranslucent(true);

        renderer_->AddActor(actor);
        actors_.push_back(actor);
    }
}

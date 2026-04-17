#include "QtBFieldStreamOverlay.h"

#include "Colormap.h"

#include "../calculators/QtBiotSavartCalc.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/QtFrame.h"

#include <QElapsedTimer>
#include <QLoggingCategory>

#include <vtkDataObject.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <cmath>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cStream, "h5reader.overlay.bstream")

// Grid resolution per ring — 10³ = 1000 evaluation points. BS kernel
// at ~1 μs per eval → ~1 ms per ring per frame, well inside budget.
constexpr int kGridDim = 10;

// Grid extent around ring centre, in Ångström. Larger than the field-
// grid overlay because streamlines benefit from room to integrate.
constexpr double kGridExtentA = 6.0;

// Seed streamlines on a circle at 1.5× ring radius in the ring plane.
// 12 seeds give an evenly distributed flower pattern.
constexpr int kSeedCount = 12;
constexpr double kSeedRadiusScale = 1.5;

// Max streamline propagation length (Å) and integrator parameters.
// Matches ui/src/ButterflyOverlay.cpp's vtkStreamTracer settings.
constexpr double kMaxPropagation      = 12.0;
constexpr double kInitialStep         = 0.1;
constexpr double kMinStep             = 0.01;
constexpr double kMaxStep             = 0.5;
constexpr int    kMaxSteps            = 2000;

// Tube filter — visual thickness of the streamlines.
constexpr double kTubeRadius = 0.04;
constexpr int    kTubeSides  = 8;

}  // namespace

QtBFieldStreamOverlay::QtBFieldStreamOverlay(
    vtkSmartPointer<vtkRenderer>                  renderer,
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
    QObject* parent)
    : QObject(parent),
      renderer_(std::move(renderer)),
      renderWindow_(std::move(renderWindow))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtBFieldStreamOverlay"));
}

QtBFieldStreamOverlay::~QtBFieldStreamOverlay() {
    for (const auto& rs : rings_) {
        if (rs.actor) renderer_->RemoveActor(rs.actor);
    }
}

void QtBFieldStreamOverlay::Build(const model::QtProtein&      protein,
                                    const model::QtConformation& conformation) {
    ASSERT_THREAD(this);

    if (protein_ == &protein && conformation_ == &conformation && !rings_.empty())
        return;

    for (const auto& rs : rings_) {
        if (rs.actor) renderer_->RemoveActor(rs.actor);
    }
    rings_.clear();

    protein_      = &protein;
    conformation_ = &conformation;

    const int nTotal = kGridDim * kGridDim * kGridDim;
    const size_t n_rings = protein.ringCount();
    rings_.resize(n_rings);

    for (size_t ri = 0; ri < n_rings; ++ri) {
        auto& rs = rings_[ri];

        rs.gridPoints = vtkSmartPointer<vtkPoints>::New();
        rs.gridPoints->SetNumberOfPoints(nTotal);

        rs.vectors = vtkSmartPointer<vtkDoubleArray>::New();
        rs.vectors->SetName("B-field");
        rs.vectors->SetNumberOfComponents(3);
        rs.vectors->SetNumberOfTuples(nTotal);

        rs.magnitudes = vtkSmartPointer<vtkDoubleArray>::New();
        rs.magnitudes->SetName("B-magnitude");
        rs.magnitudes->SetNumberOfComponents(1);
        rs.magnitudes->SetNumberOfTuples(nTotal);

        rs.grid = vtkSmartPointer<vtkStructuredGrid>::New();
        rs.grid->SetDimensions(kGridDim, kGridDim, kGridDim);
        rs.grid->SetPoints(rs.gridPoints);
        rs.grid->GetPointData()->SetVectors(rs.vectors);
        rs.grid->GetPointData()->AddArray(rs.magnitudes);

        rs.seedPoints = vtkSmartPointer<vtkPoints>::New();
        rs.seedPoints->SetNumberOfPoints(kSeedCount);
        rs.seedVerts = vtkSmartPointer<vtkCellArray>::New();
        for (vtkIdType i = 0; i < kSeedCount; ++i) {
            rs.seedVerts->InsertNextCell(1, &i);
        }
        rs.seedPoly = vtkSmartPointer<vtkPolyData>::New();
        rs.seedPoly->SetPoints(rs.seedPoints);
        rs.seedPoly->SetVerts(rs.seedVerts);

        rs.tracer = vtkSmartPointer<vtkStreamTracer>::New();
        rs.tracer->SetInputData(rs.grid);
        rs.tracer->SetSourceData(rs.seedPoly);
        rs.tracer->SetIntegrationDirectionToBoth();
        rs.tracer->SetMaximumPropagation(kMaxPropagation);
        rs.tracer->SetIntegratorTypeToRungeKutta45();
        rs.tracer->SetInitialIntegrationStep(kInitialStep);
        rs.tracer->SetMinimumIntegrationStep(kMinStep);
        rs.tracer->SetMaximumIntegrationStep(kMaxStep);
        rs.tracer->SetMaximumNumberOfSteps(kMaxSteps);
        rs.tracer->SetInputArrayToProcess(
            0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "B-field");

        rs.tubes = vtkSmartPointer<vtkTubeFilter>::New();
        rs.tubes->SetInputConnection(rs.tracer->GetOutputPort());
        rs.tubes->SetRadius(kTubeRadius);
        rs.tubes->SetNumberOfSides(kTubeSides);
        rs.tubes->CappingOn();

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(rs.tubes->GetOutputPort());
        mapper->SetScalarModeToUsePointFieldData();
        mapper->SelectColorArray("B-magnitude");
        // Scalar range and lookup table set on first UpdateRing when
        // we know the per-frame magnitude range.

        rs.actor = vtkSmartPointer<vtkActor>::New();
        rs.actor->SetMapper(mapper);
        rs.actor->GetProperty()->SetOpacity(0.85);
        rs.actor->GetProperty()->SetInterpolationToPhong();
        rs.actor->SetForceTranslucent(true);
        rs.actor->SetVisibility(visible_ ? 1 : 0);

        renderer_->AddActor(rs.actor);
    }

    qCInfo(cStream).noquote()
        << "Built B-field streamline overlay |"
        << "rings=" << n_rings
        << "| grid=" << kGridDim << "³"
        << "| seeds=" << kSeedCount;
}

void QtBFieldStreamOverlay::UpdateRing(size_t ri, int t) {
    auto& rs = rings_[ri];
    if (!rs.grid) return;

    const auto& ring  = protein_->ring(ri);
    const auto& frame = conformation_->frame(static_cast<size_t>(t));
    const auto geo    = frame.ringGeometry(ri);
    if (geo.radius < 1e-6) {
        rs.actor->SetVisibility(0);
        return;
    }
    if (visible_) rs.actor->SetVisibility(1);

    const auto vertices = frame.ringVertices(ri);
    if (vertices.size() < 3) {
        rs.actor->SetVisibility(0);
        return;
    }

    // Ring-local orthonormal basis {u, v, n} via the shared helper.
    const auto basis = model::OrthoBasisFromNormal(geo.normal);
    const model::Vec3& n = basis.n;
    const model::Vec3& u = basis.u;
    const model::Vec3& v = basis.v;

    const double spacing = 2.0 * kGridExtentA / (kGridDim - 1);
    const double intensityNA = ring.LiteratureIntensity();
    const double lobeOffsetA = ring.JBLobeOffset();

    // Update grid points (ring-local → world) + B-field vectors.
    double magMin =  std::numeric_limits<double>::infinity();
    double magMax = -std::numeric_limits<double>::infinity();

    for (int iz = 0; iz < kGridDim; ++iz) {
        const double fz = -kGridExtentA + iz * spacing;
        for (int iy = 0; iy < kGridDim; ++iy) {
            const double fy = -kGridExtentA + iy * spacing;
            for (int ix = 0; ix < kGridDim; ++ix) {
                const double fx = -kGridExtentA + ix * spacing;

                const model::Vec3 p = geo.center
                                      + fx * u
                                      + fy * v
                                      + fz * n;
                const int idx = ix + iy * kGridDim + iz * kGridDim * kGridDim;
                rs.gridPoints->SetPoint(idx, p.x(), p.y(), p.z());

                const model::Vec3 B = calculators::EvaluateBField(
                    p, geo, vertices, lobeOffsetA, intensityNA);
                rs.vectors->SetTuple3(idx, B.x(), B.y(), B.z());

                const double mag = B.norm();
                rs.magnitudes->SetValue(idx, mag);
                if (std::isfinite(mag)) {
                    magMin = std::min(magMin, mag);
                    magMax = std::max(magMax, mag);
                }
            }
        }
    }
    rs.gridPoints->Modified();
    rs.vectors->Modified();
    rs.magnitudes->Modified();
    rs.grid->Modified();

    // Seed points: circle at kSeedRadiusScale × ring radius in the
    // ring plane (spanned by u, v).
    const double seedRadius = kSeedRadiusScale * geo.radius;
    for (int i = 0; i < kSeedCount; ++i) {
        const double theta = 2.0 * M_PI * i / kSeedCount;
        const model::Vec3 sp = geo.center
                               + seedRadius *
                                 (std::cos(theta) * u + std::sin(theta) * v);
        rs.seedPoints->SetPoint(i, sp.x(), sp.y(), sp.z());
    }
    rs.seedPoints->Modified();
    rs.seedPoly->Modified();

    // Update colour map range per frame so the field stays legible as
    // |B| varies (ring geometry shifts, atoms diffuse away).
    const double absMax = std::max(std::abs(magMin), std::abs(magMax));
    if (absMax > 1e-10) {
        auto ctf = MakeDivergingCTF(0.0, absMax);
        auto* mapper = dynamic_cast<vtkPolyDataMapper*>(rs.actor->GetMapper());
        if (mapper) {
            mapper->SetLookupTable(ctf);
            mapper->SetScalarRange(0.0, absMax);
        }
    }
}

void QtBFieldStreamOverlay::setFrame(int t) {
    ASSERT_THREAD(this);
    if (!protein_ || !conformation_) return;
    if (t < 0 || static_cast<size_t>(t) >= conformation_->frameCount()) return;
    if (!visible_) return;

    QElapsedTimer timer;
    timer.start();
    for (size_t ri = 0; ri < rings_.size(); ++ri) UpdateRing(ri, t);
    // DEBUG level — restored from INFO after the pipeline was verified
    // running end-to-end. Flip back to INFO if streamlines stop
    // appearing and we need to trace the pipeline again.
    qCDebug(cStream).noquote()
        << "frame" << t << "|" << rings_.size() << "rings | eval"
        << timer.elapsed() << "ms";
}

void QtBFieldStreamOverlay::setVisible(bool visible) {
    ASSERT_THREAD(this);
    visible_ = visible;
    for (auto& rs : rings_) {
        if (rs.actor) rs.actor->SetVisibility(visible ? 1 : 0);
    }
}

}  // namespace h5reader::app

#include "HeroshotButterflyOverlay.h"

#include "../calculators/QtBiotSavartCalc.h"
#include "../calculators/QtHaighMallionCalc.h"
#include "../calculators/QtPhysicalConstants.h"
#include "../model/Conformation.h"
#include "../model/ConformationGeometry.h"
#include "../model/QtProtein.h"
#include "../model/QtRing.h"

#include <vtkPointData.h>
#include <vtkProperty.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace h5reader::app {

namespace {

constexpr int kMinGridDim = 16;
constexpr int kMaxGridDim = 72;

void styleSurfaceActor(vtkActor* actor,
                       double r,
                       double g,
                       double b,
                       double opacity,
                       bool visible) {
    if (!actor)
        return;
    auto* prop = actor->GetProperty();
    prop->SetColor(r, g, b);
    prop->SetOpacity(std::clamp(opacity, 0.0, 1.0));
    prop->LightingOff();
    prop->SetAmbient(1.0);
    prop->SetDiffuse(0.0);
    prop->SetSpecular(0.0);
    actor->SetForceTranslucent(opacity < 1.0);
    actor->SetVisibility(visible ? 1 : 0);
}

double sampleT0(const model::Vec3& p,
                const model::RingGeometry& geo,
                const std::vector<model::Vec3>& vertices,
                double lobeOffsetA,
                double intensityNA,
                HeroshotButterflyOverlay::Mode mode) {
    switch (mode) {
        case HeroshotButterflyOverlay::Mode::BiotSavart: {
            const auto st = calculators::EvaluateShielding(
                p, geo, vertices, lobeOffsetA, intensityNA);
            return st.T0;
        }
        case HeroshotButterflyOverlay::Mode::HaighMallion: {
            const auto st = calculators::EvaluateShielding(
                p, geo, vertices, intensityNA);
            return st.T0;
        }
        case HeroshotButterflyOverlay::Mode::Sum: {
            const auto stBS = calculators::EvaluateShielding(
                p, geo, vertices, lobeOffsetA, intensityNA);
            const auto stHM = calculators::EvaluateShielding(
                p, geo, vertices, intensityNA);
            return stBS.T0 + stHM.T0;
        }
    }
    return 0.0;
}

}  // namespace

HeroshotButterflyOverlay::HeroshotButterflyOverlay(vtkSmartPointer<vtkRenderer> renderer)
    : renderer_(std::move(renderer)) {}

HeroshotButterflyOverlay::~HeroshotButterflyOverlay() {
    clear();
}

void HeroshotButterflyOverlay::clear() {
    if (renderer_) {
        for (const Pipeline& pipe : pipelines_) {
            if (pipe.actorShielded)
                renderer_->RemoveActor(pipe.actorShielded);
            if (pipe.actorDeshielded)
                renderer_->RemoveActor(pipe.actorDeshielded);
        }
    }
    pipelines_.clear();
    stats_ = Stats{};
}

bool HeroshotButterflyOverlay::show(const model::QtProtein& protein,
                                    model::Conformation& conformation,
                                    std::size_t ringIdx,
                                    std::size_t frame,
                                    const Style& style) {
    clear();
    if (!renderer_ || ringIdx >= protein.ringCount() || frame >= conformation.frameCount())
        return false;

    const std::vector<model::Vec3> vertices =
        model::RingVertices(conformation, ringIdx, frame);
    const model::RingGeometry geo = model::FitRingGeometry(vertices);
    if (vertices.size() < 3 || geo.radius < 1e-6 || geo.normal.norm() < 1e-12)
        return false;

    const model::QtRing& ring = protein.ring(ringIdx);
    const double intensityNA = ring.LiteratureIntensity();
    const double lobeOffsetA = ring.JohnsonBoveyLobeOffset();
    const int dim = std::clamp(style.gridDim, kMinGridDim, kMaxGridDim);
    const double extent = calculators::FIELD_GRID_EXTENT_A;
    const double spacing = 2.0 * extent / static_cast<double>(dim - 1);
    const double threshold = std::max(0.0, style.thresholdPpm);
    const double sigma = std::max(0.1, style.gaussianExtentA);
    const double peak = std::max(0.0, style.gaussianPeak);
    Stats stats;
    stats.minT0 = std::numeric_limits<double>::infinity();
    stats.maxT0 = -std::numeric_limits<double>::infinity();

    Pipeline pipe;
    pipe.scalars = vtkSmartPointer<vtkFloatArray>::New();
    pipe.scalars->SetName("heroshot_T0");
    pipe.scalars->SetNumberOfTuples(dim * dim * dim);

    pipe.imageData = vtkSmartPointer<vtkImageData>::New();
    pipe.imageData->SetDimensions(dim, dim, dim);
    pipe.imageData->SetSpacing(spacing, spacing, spacing);
    const double origin[3] = {
        geo.center.x() - extent,
        geo.center.y() - extent,
        geo.center.z() - extent,
    };
    pipe.imageData->SetOrigin(origin[0], origin[1], origin[2]);
    pipe.imageData->GetPointData()->SetScalars(pipe.scalars);

    for (int iz = 0; iz < dim; ++iz) {
        for (int iy = 0; iy < dim; ++iy) {
            for (int ix = 0; ix < dim; ++ix) {
                const model::Vec3 p(origin[0] + ix * spacing,
                                    origin[1] + iy * spacing,
                                    origin[2] + iz * spacing);
                double t0 = sampleT0(p, geo, vertices, lobeOffsetA,
                                     intensityNA, style.mode);
                const double rDist = (p - geo.center).norm();
                const double window = std::exp(-(rDist * rDist) /
                                               (2.0 * sigma * sigma));
                t0 *= peak * window;
                if (std::isfinite(t0)) {
                    stats.minT0 = std::min(stats.minT0, t0);
                    stats.maxT0 = std::max(stats.maxT0, t0);
                }
                const bool onBoundary =
                    (ix == 0 || ix == dim - 1 ||
                     iy == 0 || iy == dim - 1 ||
                     iz == 0 || iz == dim - 1);
                const int idx = ix + iy * dim + iz * dim * dim;
                pipe.scalars->SetValue(
                    idx,
                    onBoundary ? 0.0f
                               : (std::isfinite(t0) ? static_cast<float>(t0) : 0.0f));
            }
        }
    }
    pipe.scalars->Modified();
    pipe.imageData->Modified();

    pipe.producer = vtkSmartPointer<vtkTrivialProducer>::New();
    pipe.producer->SetOutput(pipe.imageData);

    pipe.contourShielded = vtkSmartPointer<vtkContourFilter>::New();
    pipe.contourShielded->SetInputConnection(pipe.producer->GetOutputPort());
    pipe.contourShielded->SetValue(0, -threshold);
    pipe.contourShielded->Update();
    stats.shieldedPoints = static_cast<std::size_t>(
        pipe.contourShielded->GetOutput()->GetNumberOfPoints());
    stats.shieldedCells = static_cast<std::size_t>(
        pipe.contourShielded->GetOutput()->GetNumberOfCells());

    pipe.contourDeshielded = vtkSmartPointer<vtkContourFilter>::New();
    pipe.contourDeshielded->SetInputConnection(pipe.producer->GetOutputPort());
    pipe.contourDeshielded->SetValue(0, threshold);
    pipe.contourDeshielded->Update();
    stats.deshieldedPoints = static_cast<std::size_t>(
        pipe.contourDeshielded->GetOutput()->GetNumberOfPoints());
    stats.deshieldedCells = static_cast<std::size_t>(
        pipe.contourDeshielded->GetOutput()->GetNumberOfCells());
    if (!std::isfinite(stats.minT0))
        stats.minT0 = 0.0;
    if (!std::isfinite(stats.maxT0))
        stats.maxT0 = 0.0;

    pipe.mapperShielded = vtkSmartPointer<vtkPolyDataMapper>::New();
    pipe.mapperShielded->SetInputConnection(pipe.contourShielded->GetOutputPort());
    pipe.mapperShielded->ScalarVisibilityOff();
    pipe.actorShielded = vtkSmartPointer<vtkActor>::New();
    pipe.actorShielded->SetMapper(pipe.mapperShielded);
    styleSurfaceActor(pipe.actorShielded, 0.50, 0.70, 0.95,
                      style.opacity, style.showShielded);

    pipe.mapperDeshielded = vtkSmartPointer<vtkPolyDataMapper>::New();
    pipe.mapperDeshielded->SetInputConnection(pipe.contourDeshielded->GetOutputPort());
    pipe.mapperDeshielded->ScalarVisibilityOff();
    pipe.actorDeshielded = vtkSmartPointer<vtkActor>::New();
    pipe.actorDeshielded->SetMapper(pipe.mapperDeshielded);
    styleSurfaceActor(pipe.actorDeshielded, 0.95, 0.55, 0.45,
                      style.opacity, style.showDeshielded);

    renderer_->AddActor(pipe.actorShielded);
    renderer_->AddActor(pipe.actorDeshielded);
    pipelines_.push_back(std::move(pipe));
    stats_ = stats;
    return true;
}

}  // namespace h5reader::app

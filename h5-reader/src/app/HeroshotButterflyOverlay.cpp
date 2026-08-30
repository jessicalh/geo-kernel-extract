#include "HeroshotButterflyOverlay.h"

#include "../calculators/QtBiotSavartCalc.h"
#include "../calculators/QtHaighMallionCalc.h"
#include "../calculators/QtPhysicalConstants.h"
#include "../model/Conformation.h"
#include "../model/ConformationGeometry.h"
#include "../model/QtProtein.h"
#include "../model/QtRing.h"
#include "../physics/CircularRingCurrent.h"

#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>

namespace h5reader::app {

namespace {

constexpr int kMinGridDim = 16;
constexpr int kMaxGridDim = 72;
constexpr double kPi = 3.141592653589793238462643383279502884;

struct SurfaceAppearance {
    std::array<double, 3> color = {1.0, 1.0, 1.0};
    double opacity = 1.0;
    bool visible = true;
};

struct RingSampleSource {
    std::size_t ring = 0;
    model::RingGeometry geometry;
    std::vector<model::Vec3> vertices;
    double planeRmsA = 0.0;
    double lobeOffsetA = 0.0;
    double intensityNA = 0.0;
    std::optional<physics::CircularRingParameters> circular;
};

void styleSurfaceActor(vtkActor* actor, const SurfaceAppearance& appearance) {
    if (!actor)
        return;
    auto* prop = actor->GetProperty();
    prop->SetColor(appearance.color[0], appearance.color[1], appearance.color[2]);
    prop->SetOpacity(std::clamp(appearance.opacity, 0.0, 1.0));
    prop->LightingOff();
    prop->SetAmbient(1.0);
    prop->SetDiffuse(0.0);
    prop->SetSpecular(0.0);
    actor->SetForceTranslucent(appearance.opacity < 1.0);
    actor->SetVisibility(appearance.visible ? 1 : 0);
}

vtkSmartPointer<vtkActor>
makeCircularSourceLoopActor(const RingSampleSource& source, const HeroshotButterflyOverlay::Style& style) {
    if (!source.circular)
        return nullptr;
    const double normalNorm = source.geometry.normal.norm();
    if (!(normalNorm > 1.0e-12))
        return nullptr;
    const model::Vec3 normal = source.geometry.normal / normalNorm;

    model::Vec3 u = source.vertices.front() - source.geometry.center;
    u -= normal * u.dot(normal);
    if (u.norm() <= 1.0e-12)
        u = normal.unitOrthogonal();
    else
        u.normalize();
    model::Vec3 v = normal.cross(u);
    if (v.norm() <= 1.0e-12)
        return nullptr;
    v.normalize();

    const int resolution = std::clamp(style.sourceLoopResolution, 32, 360);
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    for (const double side : {-1.0, 1.0}) {
        const vtkIdType first = points->GetNumberOfPoints();
        const model::Vec3 loopCenter = source.geometry.center + side * source.circular->lobeOffsetA * normal;
        for (int i = 0; i < resolution; ++i) {
            const double angle = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(resolution);
            const model::Vec3 point = loopCenter + source.circular->radiusA * (std::cos(angle) * u + std::sin(angle) * v);
            points->InsertNextPoint(point.x(), point.y(), point.z());
        }
        lines->InsertNextCell(static_cast<vtkIdType>(resolution) + 1);
        for (int i = 0; i < resolution; ++i)
            lines->InsertCellPoint(first + static_cast<vtkIdType>(i));
        lines->InsertCellPoint(first);
    }

    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(lines);

    auto producer = vtkSmartPointer<vtkTrivialProducer>::New();
    producer->SetOutput(polyData);
    auto tube = vtkSmartPointer<vtkTubeFilter>::New();
    tube->SetInputConnection(producer->GetOutputPort());
    tube->SetRadius(std::clamp(style.sourceLoopTubeRadiusA, 0.002, 0.12));
    tube->SetNumberOfSides(12);
    tube->CappingOff();

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tube->GetOutputPort());
    mapper->ScalarVisibilityOff();

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    auto* property = actor->GetProperty();
    property->SetColor(0.78, 0.52, 0.12);
    property->SetOpacity(std::clamp(style.sourceLoopOpacity, 0.0, 1.0));
    property->LightingOff();
    property->SetAmbient(1.0);
    property->SetDiffuse(0.0);
    property->SetSpecular(0.0);
    actor->SetForceTranslucent(style.sourceLoopOpacity < 1.0);
    return actor;
}

double sampleT0(const model::Vec3& p, const RingSampleSource& source, HeroshotButterflyOverlay::Mode mode) {
    switch (mode) {
    case HeroshotButterflyOverlay::Mode::BiotSavart: {
        const auto st =
            calculators::EvaluateShielding(p, source.geometry, source.vertices, source.lobeOffsetA, source.intensityNA);
        return st.T0;
    }
    case HeroshotButterflyOverlay::Mode::HaighMallion: {
        const auto st = calculators::EvaluateShielding(p, source.geometry, source.vertices, source.intensityNA);
        return st.T0;
    }
    case HeroshotButterflyOverlay::Mode::Sum: {
        const auto stBS =
            calculators::EvaluateShielding(p, source.geometry, source.vertices, source.lobeOffsetA, source.intensityNA);
        const auto stHM = calculators::EvaluateShielding(p, source.geometry, source.vertices, source.intensityNA);
        return stBS.T0 + stHM.T0;
    }
    case HeroshotButterflyOverlay::Mode::CircularCandidateA: {
        if (!source.circular)
            return 0.0;
        const physics::CircularRingPlane plane{source.geometry.center, source.geometry.normal, source.planeRmsA};
        const auto st = physics::EvaluateCircularShielding(p, plane, *source.circular);
        return st ? st->T0 : 0.0;
    }
    }
    return 0.0;
}

}  // namespace

HeroshotButterflyOverlay::HeroshotButterflyOverlay(vtkSmartPointer<vtkRenderer> renderer) : renderer_(std::move(renderer)) {}

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
            for (const auto& actor : pipe.sourceLoopActors)
                if (actor)
                    renderer_->RemoveActor(actor);
        }
    }
    pipelines_.clear();
    circularSources_.clear();
    stats_ = Stats{};
}

std::size_t HeroshotButterflyOverlay::sourceLoopActorCount() const {
    std::size_t count = 0;
    for (const Pipeline& pipeline : pipelines_)
        count += pipeline.sourceLoopActors.size();
    return count;
}

bool HeroshotButterflyOverlay::show(const model::QtProtein& protein,
                                    const model::Conformation& conformation,
                                    std::size_t ringIdx,
                                    std::size_t frame) {
    return show(protein, conformation, ringIdx, frame, Style{});
}

bool HeroshotButterflyOverlay::show(const model::QtProtein& protein,
                                    const model::Conformation& conformation,
                                    std::size_t ringIdx,
                                    std::size_t frame,
                                    const Style& style) {
    return show(protein, conformation, std::vector<std::size_t>{ringIdx}, frame, style);
}

bool HeroshotButterflyOverlay::show(const model::QtProtein& protein,
                                    const model::Conformation& conformation,
                                    const std::vector<std::size_t>& ringIndices,
                                    std::size_t frame,
                                    const Style& style) {
    clear();
    if (!renderer_ || ringIndices.empty() || frame >= conformation.frameCount())
        return false;
    if (style.showSourceLoops && style.mode != Mode::CircularCandidateA)
        return false;

    std::vector<RingSampleSource> sources;
    sources.reserve(ringIndices.size());
    model::Vec3 gridCenter = model::Vec3::Zero();
    for (const std::size_t ringIdx : ringIndices) {
        if (ringIdx >= protein.ringCount())
            return false;
        const model::QtRing& ring = protein.ring(ringIdx);
        RingSampleSource source;
        source.ring = ringIdx;
        source.vertices = model::RingVertices(conformation, ringIdx, frame);
        source.geometry = model::FitRingGeometry(source.vertices);
        if (source.vertices.size() < 3 || source.geometry.radius < 1e-6 || source.geometry.normal.norm() < 1e-12) {
            return false;
        }
        source.lobeOffsetA = ring.JohnsonBoveyLobeOffset();
        source.intensityNA = ring.LiteratureIntensity();
        if (style.mode == Mode::CircularCandidateA) {
            const model::Vec3 windingNormal = source.geometry.normal;
            const auto plane = physics::FitCircularRingPlane(source.vertices);
            if (!plane)
                return false;
            source.geometry.center = plane->center;
            source.geometry.normal = plane->normal;
            if (source.geometry.normal.dot(windingNormal) < 0.0)
                source.geometry.normal *= -1.0;
            source.planeRmsA = plane->planeRmsA;
            int protonationVariant = -1;
            if (ring.parentResidueIndex >= 0 && static_cast<std::size_t>(ring.parentResidueIndex) < protein.residueCount()) {
                protonationVariant = static_cast<int>(
                    protein.residue(static_cast<std::size_t>(ring.parentResidueIndex)).protonationVariantIndex);
            }
            source.circular = physics::CandidateACircularParameters(ring.TypeIndex(), protonationVariant);
            if (!source.circular)
                return false;
        }
        gridCenter += source.geometry.center;
        sources.push_back(std::move(source));
    }
    gridCenter /= static_cast<double>(sources.size());
    std::vector<CircularSource> computedCircularSources;
    if (style.mode == Mode::CircularCandidateA) {
        computedCircularSources.reserve(sources.size());
        for (const RingSampleSource& source : sources) {
            if (!source.circular)
                return false;
            const physics::CircularRingParameters& circular = *source.circular;
            computedCircularSources.push_back(CircularSource{
                source.ring,
                source.geometry.center,
                source.geometry.normal,
                source.planeRmsA,
                circular.radiusA,
                circular.lobeOffsetA,
                circular.currentNanoamperePerTesla,
            });
        }
    }

    const int dim = std::clamp(style.gridDim, kMinGridDim, kMaxGridDim);
    const double extent = calculators::FIELD_GRID_EXTENT_A;
    const double spacing = 2.0 * extent / static_cast<double>(dim - 1);
    const double threshold = std::max(0.0, style.thresholdPpm);
    const double sigma = std::max(0.0, style.gaussianExtentA);
    const double peak = std::max(0.0, style.gaussianPeak);
    Stats computedStats;
    computedStats.minT0 = std::numeric_limits<double>::infinity();
    computedStats.maxT0 = -std::numeric_limits<double>::infinity();

    Pipeline pipe;
    pipe.scalars = vtkSmartPointer<vtkFloatArray>::New();
    pipe.scalars->SetName("heroshot_T0");
    const vtkIdType voxelCount = static_cast<vtkIdType>(dim) * static_cast<vtkIdType>(dim) * static_cast<vtkIdType>(dim);
    pipe.scalars->SetNumberOfTuples(voxelCount);

    pipe.imageData = vtkSmartPointer<vtkImageData>::New();
    pipe.imageData->SetDimensions(dim, dim, dim);
    pipe.imageData->SetSpacing(spacing, spacing, spacing);
    const double origin[3] = {
        gridCenter.x() - extent,
        gridCenter.y() - extent,
        gridCenter.z() - extent,
    };
    pipe.imageData->SetOrigin(origin[0], origin[1], origin[2]);
    pipe.imageData->GetPointData()->SetScalars(pipe.scalars);

    for (int iz = 0; iz < dim; ++iz) {
        for (int iy = 0; iy < dim; ++iy) {
            for (int ix = 0; ix < dim; ++ix) {
                const model::Vec3 p(origin[0] + ix * spacing, origin[1] + iy * spacing, origin[2] + iz * spacing);
                double t0 = 0.0;
                for (const RingSampleSource& source : sources)
                    t0 += sampleT0(p, source, style.mode);
                const double rDist = (p - gridCenter).norm();
                const double window = sigma > 0.0 ? std::exp(-(rDist * rDist) / (2.0 * sigma * sigma)) : 1.0;
                t0 *= peak * window;
                if (std::isfinite(t0)) {
                    computedStats.minT0 = std::min(computedStats.minT0, t0);
                    computedStats.maxT0 = std::max(computedStats.maxT0, t0);
                }
                const bool onBoundary = (ix == 0 || ix == dim - 1 || iy == 0 || iy == dim - 1 || iz == 0 || iz == dim - 1);
                const int idx = ix + iy * dim + iz * dim * dim;
                pipe.scalars->SetValue(idx, onBoundary ? 0.0f : (std::isfinite(t0) ? static_cast<float>(t0) : 0.0f));
            }
        }
    }
    pipe.scalars->Modified();
    pipe.imageData->Modified();

    pipe.producer = vtkSmartPointer<vtkTrivialProducer>::New();
    pipe.producer->SetOutput(pipe.imageData);

    pipe.contourShielded = vtkSmartPointer<vtkContourFilter>::New();
    pipe.contourShielded->SetInputConnection(pipe.producer->GetOutputPort());
    const double shieldedContour = style.mode == Mode::CircularCandidateA ? threshold : -threshold;
    const double deshieldedContour = -shieldedContour;
    pipe.contourShielded->SetValue(0, shieldedContour);
    pipe.contourShielded->Update();
    computedStats.shieldedPoints = static_cast<std::size_t>(pipe.contourShielded->GetOutput()->GetNumberOfPoints());
    computedStats.shieldedCells = static_cast<std::size_t>(pipe.contourShielded->GetOutput()->GetNumberOfCells());

    pipe.contourDeshielded = vtkSmartPointer<vtkContourFilter>::New();
    pipe.contourDeshielded->SetInputConnection(pipe.producer->GetOutputPort());
    pipe.contourDeshielded->SetValue(0, deshieldedContour);
    pipe.contourDeshielded->Update();
    computedStats.deshieldedPoints = static_cast<std::size_t>(pipe.contourDeshielded->GetOutput()->GetNumberOfPoints());
    computedStats.deshieldedCells = static_cast<std::size_t>(pipe.contourDeshielded->GetOutput()->GetNumberOfCells());
    if (!std::isfinite(computedStats.minT0))
        computedStats.minT0 = 0.0;
    if (!std::isfinite(computedStats.maxT0))
        computedStats.maxT0 = 0.0;

    pipe.mapperShielded = vtkSmartPointer<vtkPolyDataMapper>::New();
    pipe.mapperShielded->SetInputConnection(pipe.contourShielded->GetOutputPort());
    pipe.mapperShielded->ScalarVisibilityOff();
    pipe.actorShielded = vtkSmartPointer<vtkActor>::New();
    pipe.actorShielded->SetMapper(pipe.mapperShielded);
    styleSurfaceActor(pipe.actorShielded, SurfaceAppearance{{0.50, 0.70, 0.95}, style.opacity, style.showShielded});

    pipe.mapperDeshielded = vtkSmartPointer<vtkPolyDataMapper>::New();
    pipe.mapperDeshielded->SetInputConnection(pipe.contourDeshielded->GetOutputPort());
    pipe.mapperDeshielded->ScalarVisibilityOff();
    pipe.actorDeshielded = vtkSmartPointer<vtkActor>::New();
    pipe.actorDeshielded->SetMapper(pipe.mapperDeshielded);
    styleSurfaceActor(pipe.actorDeshielded, SurfaceAppearance{{0.95, 0.55, 0.45}, style.opacity, style.showDeshielded});

    if (style.showSourceLoops) {
        for (const RingSampleSource& source : sources) {
            auto actor = makeCircularSourceLoopActor(source, style);
            if (!actor)
                return false;
            pipe.sourceLoopActors.push_back(actor);
        }
    }

    renderer_->AddActor(pipe.actorShielded);
    renderer_->AddActor(pipe.actorDeshielded);
    for (const auto& actor : pipe.sourceLoopActors)
        renderer_->AddActor(actor);
    pipelines_.push_back(std::move(pipe));
    circularSources_ = std::move(computedCircularSources);
    stats_ = computedStats;
    return true;
}

}  // namespace h5reader::app

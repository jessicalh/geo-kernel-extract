#include "AtomTrackOverlay.h"

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

namespace h5reader::app {

namespace {

void insertColor(vtkUnsignedCharArray* colors, const std::array<unsigned char, 3>& rgb) {
    const unsigned char tuple[3] = {rgb[0], rgb[1], rgb[2]};
    colors->InsertNextTypedTuple(tuple);
}

vtkSmartPointer<vtkActor> makePointActor(const std::vector<AtomTrackOverlay::Sample>& samples,
                                         const std::vector<std::array<unsigned char, 3>>& colorsIn,
                                         bool currentOnly,
                                         double pointSizePixels,
                                         double opacity) {
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto vertices = vtkSmartPointer<vtkCellArray>::New();
    auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("sample_intensity_color");

    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (currentOnly && !samples[i].current)
            continue;
        const vtkIdType id =
            points->InsertNextPoint(samples[i].position.x(),
                                    samples[i].position.y(),
                                    samples[i].position.z());
        vertices->InsertNextCell(1);
        vertices->InsertCellPoint(id);
        insertColor(colors, colorsIn[i]);
    }
    if (points->GetNumberOfPoints() == 0)
        return nullptr;

    auto poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(points);
    poly->SetVerts(vertices);
    poly->GetPointData()->SetScalars(colors);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(poly);
    mapper->ScalarVisibilityOn();
    mapper->SetColorModeToDirectScalars();
    mapper->SetScalarModeToUsePointData();

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    auto* prop = actor->GetProperty();
    prop->SetOpacity(std::clamp(opacity, 0.0, 1.0));
    prop->SetPointSize(static_cast<float>(std::clamp(pointSizePixels, 1.0, 72.0)));
    prop->SetAmbient(1.0);
    prop->SetDiffuse(0.0);
    prop->SetSpecular(0.0);
    prop->SetRepresentationToPoints();
    prop->RenderPointsAsSpheresOff();
    actor->SetForceTranslucent(opacity < 1.0);
    return actor;
}

vtkSmartPointer<vtkActor> makeSphereGlyphActor(
    const std::vector<AtomTrackOverlay::Sample>& samples,
    const std::vector<std::array<unsigned char, 3>>& colorsIn,
    bool currentOnly,
    double radiusA,
    double opacity) {
    auto sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetRadius(std::clamp(radiusA, 0.002, 0.20));
    sphere->SetThetaResolution(8);
    sphere->SetPhiResolution(6);
    sphere->Update();
    vtkPolyData* spherePoly = sphere->GetOutput();
    if (!spherePoly || !spherePoly->GetPoints() || !spherePoly->GetPolys())
        return nullptr;

    auto points = vtkSmartPointer<vtkPoints>::New();
    auto polys = vtkSmartPointer<vtkCellArray>::New();
    auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("sample_intensity_color");

    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (currentOnly && !samples[i].current)
            continue;
        const vtkIdType base = points->GetNumberOfPoints();
        for (vtkIdType p = 0; p < spherePoly->GetNumberOfPoints(); ++p) {
            double raw[3] = {};
            spherePoly->GetPoint(p, raw);
            points->InsertNextPoint(samples[i].position.x() + raw[0],
                                    samples[i].position.y() + raw[1],
                                    samples[i].position.z() + raw[2]);
            insertColor(colors, colorsIn[i]);
        }

        spherePoly->GetPolys()->InitTraversal();
        vtkIdType npts = 0;
        const vtkIdType* ids = nullptr;
        while (spherePoly->GetPolys()->GetNextCell(npts, ids)) {
            polys->InsertNextCell(static_cast<int>(npts));
            for (vtkIdType k = 0; k < npts; ++k)
                polys->InsertCellPoint(base + ids[k]);
        }
    }
    if (points->GetNumberOfPoints() == 0)
        return nullptr;

    auto poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(points);
    poly->SetPolys(polys);
    poly->GetPointData()->SetScalars(colors);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(poly);
    mapper->ScalarVisibilityOn();
    mapper->SetColorModeToDirectScalars();
    mapper->SetScalarModeToUsePointData();

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    auto* prop = actor->GetProperty();
    prop->SetOpacity(std::clamp(opacity, 0.0, 1.0));
    prop->LightingOff();
    prop->SetAmbient(1.0);
    prop->SetDiffuse(0.0);
    prop->SetSpecular(0.0);
    actor->SetForceTranslucent(opacity < 1.0);
    return actor;
}

vtkSmartPointer<vtkActor> makeLineActor(const std::vector<AtomTrackOverlay::Sample>& samples,
                                        const std::vector<std::array<unsigned char, 3>>& colorsIn,
                                        double lineWidthPixels,
                                        double opacity) {
    if (samples.size() < 2)
        return nullptr;

    auto points = vtkSmartPointer<vtkPoints>::New();
    for (const auto& sample : samples) {
        points->InsertNextPoint(sample.position.x(),
                                sample.position.y(),
                                sample.position.z());
    }

    auto lines = vtkSmartPointer<vtkCellArray>::New();
    auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("sample_edge_color");
    for (std::size_t i = 1; i < samples.size(); ++i) {
        lines->InsertNextCell(2);
        lines->InsertCellPoint(static_cast<vtkIdType>(i - 1));
        lines->InsertCellPoint(static_cast<vtkIdType>(i));
        const std::array<unsigned char, 3> c{
            static_cast<unsigned char>((static_cast<int>(colorsIn[i - 1][0]) +
                                        static_cast<int>(colorsIn[i][0])) / 2),
            static_cast<unsigned char>((static_cast<int>(colorsIn[i - 1][1]) +
                                        static_cast<int>(colorsIn[i][1])) / 2),
            static_cast<unsigned char>((static_cast<int>(colorsIn[i - 1][2]) +
                                        static_cast<int>(colorsIn[i][2])) / 2),
        };
        insertColor(colors, c);
    }

    auto poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(points);
    poly->SetLines(lines);
    poly->GetCellData()->SetScalars(colors);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(poly);
    mapper->ScalarVisibilityOn();
    mapper->SetColorModeToDirectScalars();
    mapper->SetScalarModeToUseCellData();

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    auto* prop = actor->GetProperty();
    prop->SetOpacity(std::clamp(opacity, 0.0, 1.0));
    prop->SetLineWidth(static_cast<float>(std::clamp(lineWidthPixels, 1.0, 12.0)));
    prop->SetAmbient(1.0);
    prop->SetDiffuse(0.0);
    prop->SetSpecular(0.0);
    actor->SetForceTranslucent(opacity < 1.0);
    return actor;
}

}  // namespace

AtomTrackOverlay::AtomTrackOverlay(vtkSmartPointer<vtkRenderer> sceneRenderer)
    : renderer_(std::move(sceneRenderer)) {}

AtomTrackOverlay::~AtomTrackOverlay() {
    clear();
}

void AtomTrackOverlay::addActor(const vtkSmartPointer<vtkActor>& actor) {
    if (!renderer_ || !actor)
        return;
    renderer_->AddActor(actor);
    actors_.push_back(actor);
}

std::array<unsigned char, 3> AtomTrackOverlay::colorFor(double value,
                                                        double scale,
                                                        ColorMode mode,
                                                        double gamma,
                                                        double minFraction) {
    if (!std::isfinite(value) || !std::isfinite(scale) || scale <= 1e-12)
        return {154, 158, 164};

    const double raw = (mode == ColorMode::Absolute) ? std::abs(value) : value;
    const double gammaSafe = std::clamp(gamma, 0.1, 4.0);
    const double minSafe = std::clamp(minFraction, 0.0, 0.6);
    double t = std::clamp(std::abs(raw) / scale, 0.0, 1.0);
    t = std::pow(t, gammaSafe);
    if (std::abs(raw) > 1e-14)
        t = std::max(t, minSafe);
    const std::array<double, 3> neutral{0.68, 0.70, 0.70};
    const std::array<double, 3> warm{1.00, 0.62, 0.14};
    const std::array<double, 3> cool{0.10, 0.68, 0.96};
    const std::array<double, 3> peak{1.00, 0.26, 0.16};
    const auto& end = mode == ColorMode::Absolute ? peak : (raw >= 0.0 ? warm : cool);

    std::array<unsigned char, 3> out{};
    for (int i = 0; i < 3; ++i)
        out[static_cast<std::size_t>(i)] = static_cast<unsigned char>(
            std::clamp(std::lround((neutral[i] + (end[i] - neutral[i]) * t) * 255.0),
                       0L, 255L));
    return out;
}

void AtomTrackOverlay::clear() {
    if (renderer_) {
        for (const auto& actor : actors_)
            if (actor) renderer_->RemoveActor(actor);
    }
    actors_.clear();
}

void AtomTrackOverlay::show(const std::vector<Sample>& samples) {
    show(samples, Style{});
}

void AtomTrackOverlay::show(const std::vector<Sample>& samples, const Style& style) {
    clear();
    if (!renderer_ || samples.empty())
        return;

    double scale = 0.0;
    if (std::isfinite(style.colorScale) && style.colorScale > 1e-12) {
        scale = style.colorScale;
    } else {
        for (const Sample& sample : samples) {
            if (std::isfinite(sample.intensity))
                scale = std::max(scale, std::abs(sample.intensity));
        }
    }
    if (scale <= 1e-12)
        scale = 1.0;

    std::vector<std::array<unsigned char, 3>> colors;
    colors.reserve(samples.size());
    for (const Sample& sample : samples)
        colors.push_back(colorFor(sample.intensity,
                                  scale,
                                  style.colorMode,
                                  style.colorGamma,
                                  style.minColorFraction));

    if (style.showLines)
        addActor(makeLineActor(samples, colors, style.lineWidthPixels, style.lineOpacity));

    if (style.showHalos)
        addActor(makePointActor(samples, colors, false,
                                style.pointSizePixels * style.haloScale,
                                style.haloOpacity));

    if (style.showPoints) {
        if (style.pointShape == PointShape::Sphere) {
            addActor(makeSphereGlyphActor(samples, colors, false,
                                          style.sphereRadiusA,
                                          style.pointOpacity));
        } else {
            addActor(makePointActor(samples, colors, false,
                                    style.pointSizePixels,
                                    style.pointOpacity));
        }
    }

    if (style.highlightCurrent) {
        if (style.pointShape == PointShape::Sphere) {
            addActor(makeSphereGlyphActor(samples, colors, true,
                                          style.sphereRadiusA * style.currentPointScale,
                                          std::min(1.0, style.pointOpacity + 0.04)));
        } else {
            addActor(makePointActor(samples, colors, true,
                                    style.pointSizePixels * style.currentPointScale,
                                    std::min(1.0, style.pointOpacity + 0.04)));
        }
    }
}

}  // namespace h5reader::app

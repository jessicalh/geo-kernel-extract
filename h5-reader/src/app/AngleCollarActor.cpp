#include "AngleCollarActor.h"

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>

#include <algorithm>
#include <cmath>
#include <utility>

namespace h5reader::app {

namespace {
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kMinNorm = 1e-9;

model::Vec3 normalizedOrZero(const model::Vec3& v) {
    const double n = v.norm();
    if (n <= kMinNorm)
        return model::Vec3::Zero();
    return v / n;
}

model::Vec3 projectPerpendicular(const model::Vec3& v, const model::Vec3& axisUnit) {
    return v - axisUnit * v.dot(axisUnit);
}

model::Vec3 rotateAroundAxis(const model::Vec3& v, const model::Vec3& axisUnit, double radians) {
    const double c = std::cos(radians);
    const double s = std::sin(radians);
    return v * c + axisUnit.cross(v) * s + axisUnit * (axisUnit.dot(v)) * (1.0 - c);
}

vtkSmartPointer<vtkActor> tubeActor(const std::vector<model::Vec3>& points,
                                    double radius,
                                    const std::array<double, 3>& color,
                                    double opacity) {
    auto vtkPointsObj = vtkSmartPointer<vtkPoints>::New();
    for (const auto& p : points)
        vtkPointsObj->InsertNextPoint(p.x(), p.y(), p.z());

    auto line = vtkSmartPointer<vtkCellArray>::New();
    const vtkIdType pointCount = static_cast<vtkIdType>(points.size());
    line->InsertNextCell(pointCount);
    for (vtkIdType i = 0; i < pointCount; ++i)
        line->InsertCellPoint(i);

    auto poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(vtkPointsObj);
    poly->SetLines(line);

    auto tube = vtkSmartPointer<vtkTubeFilter>::New();
    tube->SetInputData(poly);
    tube->SetRadius(std::max(0.001, radius));
    tube->SetNumberOfSides(18);
    tube->CappingOn();

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tube->GetOutputPort());

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color[0], color[1], color[2]);
    actor->GetProperty()->SetOpacity(std::clamp(opacity, 0.0, 1.0));
    actor->GetProperty()->SetAmbient(0.55);
    actor->GetProperty()->SetDiffuse(0.45);
    actor->GetProperty()->SetSpecular(0.05);
    return actor;
}

vtkSmartPointer<vtkActor> coneActor(const model::Vec3& neckCenter,
                                    const model::Vec3& rimCenter,
                                    const model::Vec3& referenceUnit,
                                    const model::Vec3& axisUnit,
                                    double neckRadius,
                                    double rimRadius,
                                    int segments,
                                    const std::array<double, 3>& color,
                                    double opacity) {
    const int n = std::max(24, segments);

    auto vtkPointsObj = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < n; ++i) {
        const double theta = (2.0 * kPi * static_cast<double>(i)) / static_cast<double>(n);
        const model::Vec3 dir = rotateAroundAxis(referenceUnit, axisUnit, theta);
        const model::Vec3 neck = neckCenter + dir * neckRadius;
        const model::Vec3 rim = rimCenter + dir * rimRadius;
        vtkPointsObj->InsertNextPoint(neck.x(), neck.y(), neck.z());
        vtkPointsObj->InsertNextPoint(rim.x(), rim.y(), rim.z());
    }

    auto polys = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 0; i < n; ++i) {
        const vtkIdType neck0 = static_cast<vtkIdType>(2 * i);
        const vtkIdType rim0 = neck0 + 1;
        const vtkIdType neck1 = static_cast<vtkIdType>(2 * ((i + 1) % n));
        const vtkIdType rim1 = neck1 + 1;
        polys->InsertNextCell(4);
        polys->InsertCellPoint(neck0);
        polys->InsertCellPoint(neck1);
        polys->InsertCellPoint(rim1);
        polys->InsertCellPoint(rim0);
    }

    auto poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(vtkPointsObj);
    poly->SetPolys(polys);

    auto normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals->SetInputData(poly);
    normals->ConsistencyOn();
    normals->SplittingOff();

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(normals->GetOutputPort());

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color[0], color[1], color[2]);
    actor->GetProperty()->SetOpacity(std::clamp(opacity, 0.0, 1.0));
    actor->GetProperty()->SetAmbient(0.50);
    actor->GetProperty()->SetDiffuse(0.45);
    actor->GetProperty()->SetSpecular(0.12);
    actor->GetProperty()->SetSpecularPower(18.0);
    actor->GetProperty()->BackfaceCullingOff();
    return actor;
}

std::vector<model::Vec3> arcPoints(const model::Vec3& center,
                                   const model::Vec3& referenceUnit,
                                   const model::Vec3& axisUnit,
                                   double radius,
                                   double angleRadians,
                                   int preferredSegments) {
    const int byAngle = static_cast<int>(std::ceil(std::abs(angleRadians) / (kPi / 48.0)));
    const int segments = std::max({8, byAngle, preferredSegments / 4});
    std::vector<model::Vec3> pts;
    pts.reserve(static_cast<std::size_t>(segments + 1));
    for (int i = 0; i <= segments; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(segments);
        const model::Vec3 dir = rotateAroundAxis(referenceUnit, axisUnit, angleRadians * t);
        pts.push_back(center + dir * radius);
    }
    return pts;
}

std::vector<model::Vec3> circlePoints(const model::Vec3& center,
                                      const model::Vec3& referenceUnit,
                                      const model::Vec3& axisUnit,
                                      double radius,
                                      int segments) {
    const int n = std::max(24, segments);
    std::vector<model::Vec3> pts;
    pts.reserve(static_cast<std::size_t>(n + 1));
    for (int i = 0; i <= n; ++i) {
        const double theta = (2.0 * kPi * static_cast<double>(i)) / static_cast<double>(n);
        pts.push_back(center + rotateAroundAxis(referenceUnit, axisUnit, theta) * radius);
    }
    return pts;
}

std::vector<model::Vec3> seamPoints(const model::Vec3& neckCenter,
                                    const model::Vec3& rimCenter,
                                    const model::Vec3& referenceUnit,
                                    const model::Vec3& axisUnit,
                                    double neckRadius,
                                    double rimRadius,
                                    double angleRadians) {
    const model::Vec3 dir = rotateAroundAxis(referenceUnit, axisUnit, angleRadians);
    return {neckCenter + dir * neckRadius, rimCenter + dir * rimRadius};
}
}  // namespace

AngleCollarActor::AngleCollarActor(vtkSmartPointer<vtkRenderer> sceneRenderer)
    : renderer_(std::move(sceneRenderer)) {}

AngleCollarActor::~AngleCollarActor() {
    clear();
}

void AngleCollarActor::clear() {
    if (renderer_) {
        for (const auto& actor : actors_)
            if (actor) renderer_->RemoveActor(actor);
    }
    actors_.clear();
}

void AngleCollarActor::show(const model::Vec3& axisStart,
                            const model::Vec3& axisEnd,
                            const model::Vec3& center,
                            const model::Vec3& referenceDirection,
                            const std::vector<Arc>& arcs,
                            const Style& style) {
    clear();
    if (!renderer_ || arcs.empty())
        return;

    const model::Vec3 axisUnit = normalizedOrZero(axisEnd - axisStart);
    if (axisUnit.squaredNorm() < kMinNorm)
        return;
    const model::Vec3 refUnit = normalizedOrZero(projectPerpendicular(referenceDirection, axisUnit));
    if (refUnit.squaredNorm() < kMinNorm)
        return;

    const double radius = std::max(0.05, style.radius);
    const double tube = std::max(0.001, style.tubeRadius);
    const double axisPad = std::max(0.0, style.axisPadding);
    const int segments = std::max(24, style.arcSegments);
    const double coneAxisSign = style.coneDirection < 0.0 ? -1.0 : 1.0;
    const model::Vec3 coneAxis = axisUnit * coneAxisSign;
    const double axisLength = (axisEnd - axisStart).norm();
    const double coneLength =
        style.coneLength > 0.0 ? style.coneLength : std::max(axisLength + axisPad, radius * 1.35);
    const double neckRadius =
        style.neckRadius > 0.0 ? style.neckRadius : std::max(0.10, radius * 0.34);
    const double rimRadius =
        style.rimRadius > 0.0 ? style.rimRadius : std::max(neckRadius * 1.8, radius * 1.24);
    const model::Vec3 neckCenter = center - coneAxis * std::min(0.18, coneLength * 0.18);
    const model::Vec3 rimCenter = center + coneAxis * coneLength;

    auto add = [&](vtkSmartPointer<vtkActor> actor) {
        renderer_->AddActor(actor);
        actors_.push_back(actor);
    };

    const std::array<double, 3> axisColor{{0.18, 0.20, 0.22}};
    const std::array<double, 3> collarColor{{0.78, 0.82, 0.86}};
    const std::array<double, 3> refColor{{0.24, 0.58, 0.86}};

    add(coneActor(neckCenter, rimCenter, refUnit, axisUnit, neckRadius, rimRadius,
                  segments, collarColor, style.coneOpacity));
    add(tubeActor({axisStart - axisUnit * axisPad, axisEnd + axisUnit * axisPad},
                  tube * 0.80, axisColor, 0.78));
    add(tubeActor(circlePoints(neckCenter, refUnit, axisUnit, neckRadius, segments),
                  tube * 0.55, collarColor, 0.46));
    add(tubeActor(circlePoints(rimCenter, refUnit, axisUnit, rimRadius, segments),
                  tube * 0.80, collarColor, 0.62));
    add(tubeActor(seamPoints(neckCenter, rimCenter, refUnit, axisUnit, neckRadius,
                             rimRadius, 0.0),
                  tube * 0.68, refColor, 0.88));

    for (const Arc& arc : arcs) {
        const double r = rimRadius * std::clamp(arc.radiusScale, 0.2, 3.0);
        const double neckR = neckRadius * std::clamp(arc.radiusScale, 0.45, 1.45);
        const model::Vec3 targetUnit = rotateAroundAxis(refUnit, axisUnit, arc.angleRadians);
        add(tubeActor(arcPoints(rimCenter, refUnit, axisUnit, r, arc.angleRadians, segments),
                      tube, arc.color, arc.opacity));
        add(tubeActor(seamPoints(neckCenter, rimCenter, refUnit, axisUnit, neckR, r,
                                 arc.angleRadians),
                      tube * 0.72, arc.color, std::min(1.0, arc.opacity + 0.12)));
        add(tubeActor({rimCenter + targetUnit * r, rimCenter + targetUnit * (r + tube * 7.5)},
                      tube * 0.90, arc.color, std::min(1.0, arc.opacity + 0.08)));
    }
}

}  // namespace h5reader::app

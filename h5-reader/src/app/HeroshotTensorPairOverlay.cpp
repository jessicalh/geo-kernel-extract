#include "HeroshotTensorPairOverlay.h"

#include <vtkColorTransferFunction.h>
#include <vtkDoubleArray.h>
#include <vtkLineSource.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkTrivialProducer.h>
#include <vtkTubeFilter.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <optional>
#include <utility>

namespace h5reader::app {

namespace {

vtkSmartPointer<vtkColorTransferFunction> signedTensorColors() {
    auto colors = vtkSmartPointer<vtkColorTransferFunction>::New();
    colors->AddRGBPoint(-1.0, 0.78, 0.16, 0.12);
    colors->AddRGBPoint(0.0, 0.84, 0.84, 0.82);
    colors->AddRGBPoint(1.0, 0.08, 0.35, 0.78);
    return colors;
}

vtkSmartPointer<vtkActor> makeTensorActor(const model::Vec3& center,
                                          const model::Mat3& input,
                                          const HeroshotTensorPairOverlay::Style& style,
                                          bool wireframe) {
    model::Mat3 tensor = 0.5 * (input + input.transpose());
    tensor -= (tensor.trace() / 3.0) * model::Mat3::Identity();
    if (!center.allFinite() || !tensor.allFinite())
        return nullptr;

    Eigen::SelfAdjointEigenSolver<model::Mat3> solver(tensor);
    if (solver.info() != Eigen::Success)
        return nullptr;
    const double maxResponse = solver.eigenvalues().cwiseAbs().maxCoeff();
    if (!(maxResponse > 1.0e-12) || !std::isfinite(maxResponse))
        return nullptr;

    auto sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetRadius(1.0);
    sphere->SetThetaResolution(style.thetaResolution);
    sphere->SetPhiResolution(style.phiResolution);
    sphere->Update();

    auto surface = vtkSmartPointer<vtkPolyData>::New();
    surface->DeepCopy(sphere->GetOutput());
    auto sign = vtkSmartPointer<vtkDoubleArray>::New();
    sign->SetName("tensor_response_sign");
    sign->SetNumberOfComponents(1);
    sign->SetNumberOfTuples(surface->GetNumberOfPoints());

    for (vtkIdType i = 0; i < surface->GetNumberOfPoints(); ++i) {
        double raw[3];
        surface->GetPoint(i, raw);
        model::Vec3 direction(raw[0], raw[1], raw[2]);
        const double norm = direction.norm();
        if (!(norm > 1.0e-12))
            return nullptr;
        direction /= norm;
        const double response = direction.dot(tensor * direction);
        const double radius = std::max(style.minimumRadiusA, style.scaleAperPpm * std::abs(response));
        const model::Vec3 point = center + radius * direction;
        surface->GetPoints()->SetPoint(i, point.x(), point.y(), point.z());
        sign->SetValue(i, std::clamp(response / maxResponse, -1.0, 1.0));
    }
    surface->GetPoints()->Modified();
    surface->GetPointData()->SetScalars(sign);

    auto producer = vtkSmartPointer<vtkTrivialProducer>::New();
    producer->SetOutput(surface);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    if (wireframe) {
        mapper->SetInputConnection(producer->GetOutputPort());
    } else {
        auto normals = vtkSmartPointer<vtkPolyDataNormals>::New();
        normals->SetInputConnection(producer->GetOutputPort());
        normals->ComputePointNormalsOn();
        normals->ComputeCellNormalsOff();
        normals->SplittingOff();
        normals->ConsistencyOn();
        mapper->SetInputConnection(normals->GetOutputPort());
    }
    mapper->SetLookupTable(signedTensorColors());
    mapper->SetColorModeToMapScalars();
    mapper->SetScalarModeToUsePointData();
    mapper->InterpolateScalarsBeforeMappingOn();
    mapper->SetScalarRange(-1.0, 1.0);

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    if (wireframe) {
        actor->GetProperty()->SetRepresentationToWireframe();
        actor->GetProperty()->SetLineWidth(static_cast<float>(style.referenceLineWidth));
        actor->GetProperty()->SetOpacity(style.referenceOpacity);
        actor->GetProperty()->LightingOff();
    } else {
        actor->GetProperty()->SetOpacity(style.candidateOpacity);
        actor->GetProperty()->SetAmbient(0.20);
        actor->GetProperty()->SetDiffuse(0.84);
        actor->GetProperty()->SetSpecular(0.32);
        actor->GetProperty()->SetSpecularPower(28.0);
        actor->GetProperty()->SetInterpolationToPhong();
    }
    return actor;
}

std::optional<model::Vec3> dominantDirector(const model::Mat3& input) {
    model::Mat3 tensor = 0.5 * (input + input.transpose());
    tensor -= (tensor.trace() / 3.0) * model::Mat3::Identity();
    if (!tensor.allFinite())
        return std::nullopt;

    Eigen::SelfAdjointEigenSolver<model::Mat3> solver(tensor);
    if (solver.info() != Eigen::Success)
        return std::nullopt;

    Eigen::Index index = 0;
    solver.eigenvalues().cwiseAbs().maxCoeff(&index);
    model::Vec3 direction = solver.eigenvectors().col(index);
    if (!direction.allFinite() || !(direction.norm() > 1.0e-12))
        return std::nullopt;
    return direction.normalized();
}

vtkSmartPointer<vtkActor> makeDirectorActor(const model::Vec3& center,
                                            const model::Vec3& direction,
                                            double halfLengthA,
                                            double radiusA,
                                            const model::Vec3& color,
                                            double opacity) {
    if (!center.allFinite() || !direction.allFinite() || !(direction.norm() > 1.0e-12) || !(halfLengthA > 0.0)
        || !(radiusA > 0.0)) {
        return nullptr;
    }

    const model::Vec3 axis = direction.normalized();
    const model::Vec3 first = center - halfLengthA * axis;
    const model::Vec3 second = center + halfLengthA * axis;

    auto line = vtkSmartPointer<vtkLineSource>::New();
    line->SetPoint1(first.x(), first.y(), first.z());
    line->SetPoint2(second.x(), second.y(), second.z());

    auto tube = vtkSmartPointer<vtkTubeFilter>::New();
    tube->SetInputConnection(line->GetOutputPort());
    tube->SetRadius(radiusA);
    tube->SetNumberOfSides(20);
    tube->CappingOn();

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tube->GetOutputPort());

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color.x(), color.y(), color.z());
    actor->GetProperty()->SetOpacity(opacity);
    actor->GetProperty()->SetAmbient(0.28);
    actor->GetProperty()->SetDiffuse(0.72);
    actor->GetProperty()->SetSpecular(0.15);
    actor->GetProperty()->SetSpecularPower(18.0);
    return actor;
}

}  // namespace

HeroshotTensorPairOverlay::HeroshotTensorPairOverlay(vtkSmartPointer<vtkRenderer> renderer) : renderer_(std::move(renderer)) {}

HeroshotTensorPairOverlay::~HeroshotTensorPairOverlay() {
    clear();
}

bool HeroshotTensorPairOverlay::show(const std::vector<Sample>& samples, const Style& style) {
    clear();
    if (!renderer_ || samples.empty() || !(style.scaleAperPpm > 0.0) || !(style.minimumRadiusA >= 0.0)
        || (!style.showCandidate && !style.showReference) || (!style.showSurfaces && !style.showDirectors)
        || !(style.candidateDirectorHalfLengthA > 0.0) || !(style.referenceDirectorHalfLengthA > 0.0)
        || !(style.referenceNormalHalfLengthA > 0.0) || !(style.directorRadiusA > 0.0)) {
        return false;
    }

    for (const Sample& sample : samples) {
        if (style.showSurfaces && style.showCandidate) {
            auto actor = makeTensorActor(sample.center, sample.candidate, style, false);
            if (!actor) {
                clear();
                return false;
            }
            renderer_->AddActor(actor);
            actors_.push_back(actor);
        }
        if (style.showSurfaces && style.showReference) {
            auto actor = makeTensorActor(sample.center, sample.reference, style, true);
            if (!actor) {
                clear();
                return false;
            }
            renderer_->AddActor(actor);
            actors_.push_back(actor);
        }
        if (style.showDirectors && style.showCandidate) {
            const auto direction = dominantDirector(sample.candidate);
            auto actor = direction ? makeDirectorActor(sample.center,
                                                       *direction,
                                                       style.candidateDirectorHalfLengthA,
                                                       1.45 * style.directorRadiusA,
                                                       model::Vec3(0.83, 0.48, 0.08),
                                                       0.96)
                                   : nullptr;
            if (!actor) {
                clear();
                return false;
            }
            renderer_->AddActor(actor);
            actors_.push_back(actor);
        }
        if (style.showDirectors && style.showReference) {
            const auto direction = dominantDirector(sample.reference);
            auto actor = direction ? makeDirectorActor(sample.center,
                                                       *direction,
                                                       style.referenceDirectorHalfLengthA,
                                                       0.72 * style.directorRadiusA,
                                                       model::Vec3(0.10, 0.35, 0.72),
                                                       1.0)
                                   : nullptr;
            if (!actor) {
                clear();
                return false;
            }
            renderer_->AddActor(actor);
            actors_.push_back(actor);
        }
        if (style.showReferenceNormal) {
            auto actor = makeDirectorActor(sample.center,
                                           sample.referenceNormal,
                                           style.referenceNormalHalfLengthA,
                                           0.45 * style.directorRadiusA,
                                           model::Vec3(0.25, 0.27, 0.29),
                                           0.72);
            if (!actor) {
                clear();
                return false;
            }
            renderer_->AddActor(actor);
            actors_.push_back(actor);
        }
    }
    sampleCount_ = samples.size();
    return true;
}

void HeroshotTensorPairOverlay::clear() {
    if (renderer_) {
        for (const auto& actor : actors_)
            if (actor)
                renderer_->RemoveActor(actor);
    }
    actors_.clear();
    sampleCount_ = 0;
}

}  // namespace h5reader::app

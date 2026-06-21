#include "CsaTensorOverlay.h"

#include <vtkCellData.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

namespace h5reader::app {

namespace {
constexpr double kBaseRadius = 0.8;   // Angstrom, the glyph's neutral size
constexpr double kAniso = 0.6;        // fractional elongation by deviation
constexpr double kAxisLength = 2.8;   // Angstrom, molecular-frame axis arms
}  // namespace

CsaTensorOverlay::CsaTensorOverlay(vtkSmartPointer<vtkRenderer> renderer, QObject* parent)
    : QObject(parent), renderer_(std::move(renderer)) {}

CsaTensorOverlay::~CsaTensorOverlay() {
    if (renderer_) {
        if (ellipsoidActor_) renderer_->RemoveActor(ellipsoidActor_);
        if (axesActor_) renderer_->RemoveActor(axesActor_);
    }
}

void CsaTensorOverlay::ensureActors() {
    if (actorsBuilt_) return;

    ellipsoidSource_ = vtkSmartPointer<vtkSphereSource>::New();
    ellipsoidSource_->SetRadius(1.0);
    ellipsoidSource_->SetThetaResolution(32);
    ellipsoidSource_->SetPhiResolution(32);
    ellipsoidTransform_ = vtkSmartPointer<vtkTransform>::New();
    ellipsoidFilter_ = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    ellipsoidFilter_->SetInputConnection(ellipsoidSource_->GetOutputPort());
    ellipsoidFilter_->SetTransform(ellipsoidTransform_);
    auto ellipMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    ellipMapper->SetInputConnection(ellipsoidFilter_->GetOutputPort());
    ellipsoidActor_ = vtkSmartPointer<vtkActor>::New();
    ellipsoidActor_->SetMapper(ellipMapper);
    ellipsoidActor_->GetProperty()->SetColor(0.85, 0.75, 0.25);
    ellipsoidActor_->GetProperty()->SetOpacity(0.55);
    ellipsoidActor_->SetVisibility(false);
    renderer_->AddActor(ellipsoidActor_);

    axesPoints_ = vtkSmartPointer<vtkPoints>::New();
    axesLines_ = vtkSmartPointer<vtkCellArray>::New();
    axesColors_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
    axesColors_->SetNumberOfComponents(3);
    axesData_ = vtkSmartPointer<vtkPolyData>::New();
    axesData_->SetPoints(axesPoints_);
    axesData_->SetLines(axesLines_);
    axesData_->GetCellData()->SetScalars(axesColors_);
    auto axesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    axesMapper->SetInputData(axesData_);
    axesMapper->SetScalarModeToUseCellData();
    axesActor_ = vtkSmartPointer<vtkActor>::New();
    axesActor_->SetMapper(axesMapper);
    axesActor_->GetProperty()->SetLineWidth(2.5);
    axesActor_->SetVisibility(false);
    renderer_->AddActor(axesActor_);

    actorsBuilt_ = true;
}

void CsaTensorOverlay::show(const model::Vec3& atomPos,
                            const model::CsaShape& shape,
                            const std::optional<model::Mat3>& molecularAxes) {
    if (!renderer_ || !shape.valid) {
        clear();
        return;
    }
    ensureActors();

    // Ellipsoid radii: neutral kBaseRadius elongated by each principal value's
    // signed deviation from isotropic, normalized so the glyph stays a visible,
    // positive ovaloid regardless of the absolute ppm scale.
    const double iso = shape.sigma_iso;
    const std::array<double, 3> dev{
        shape.principal_values[0] - iso,
        shape.principal_values[1] - iso,
        shape.principal_values[2] - iso,
    };
    double maxAbs = 0.0;
    for (double d : dev) maxAbs = std::max(maxAbs, std::abs(d));
    std::array<double, 3> radii{};
    for (std::size_t i = 0; i < 3; ++i) {
        const double frac = (maxAbs > 1e-9) ? dev[i] / maxAbs : 0.0;
        radii[i] = kBaseRadius * (1.0 + kAniso * frac);
    }

    // 4x4: column c = PAS director c scaled by radius c; column 3 = atomPos.
    auto m = vtkSmartPointer<vtkMatrix4x4>::New();
    m->Identity();
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col)
            m->SetElement(row, col,
                          shape.pas_axes(row, col) * radii[static_cast<std::size_t>(col)]);
        m->SetElement(row, 3, atomPos[row]);
    }
    ellipsoidTransform_->SetMatrix(m);
    ellipsoidActor_->SetVisibility(true);

    // Molecular-frame axes: x/y/z = R/G/B segments from the atom centre.
    axesPoints_->Reset();
    axesLines_->Reset();
    axesColors_->Reset();
    axesColors_->SetNumberOfComponents(3);
    hasAxes_ = molecularAxes.has_value();
    if (hasAxes_) {
        static const unsigned char kRgb[3][3] = {
            {220, 60, 60}, {60, 200, 60}, {70, 110, 230}};
        for (int axis = 0; axis < 3; ++axis) {
            const model::Vec3 dir = molecularAxes->col(axis);
            const vtkIdType p0 =
                axesPoints_->InsertNextPoint(atomPos[0], atomPos[1], atomPos[2]);
            const vtkIdType p1 = axesPoints_->InsertNextPoint(
                atomPos[0] + dir[0] * kAxisLength,
                atomPos[1] + dir[1] * kAxisLength,
                atomPos[2] + dir[2] * kAxisLength);
            axesLines_->InsertNextCell(2);
            axesLines_->InsertCellPoint(p0);
            axesLines_->InsertCellPoint(p1);
            axesColors_->InsertNextTypedTuple(kRgb[axis]);
        }
        axesData_->Modified();
        axesActor_->SetVisibility(true);
    } else {
        axesActor_->SetVisibility(false);
    }

    active_ = true;
}

void CsaTensorOverlay::clear() {
    if (ellipsoidActor_) ellipsoidActor_->SetVisibility(false);
    if (axesActor_) axesActor_->SetVisibility(false);
    active_ = false;
}

void CsaTensorOverlay::setVisible(bool on) {
    if (ellipsoidActor_) ellipsoidActor_->SetVisibility(on && active_);
    if (axesActor_) axesActor_->SetVisibility(on && active_ && hasAxes_);
}

}  // namespace h5reader::app

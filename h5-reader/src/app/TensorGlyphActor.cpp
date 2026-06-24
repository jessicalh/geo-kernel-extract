#include "TensorGlyphActor.h"

#include <vtkDoubleArray.h>
#include <vtkMatrix4x4.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

namespace h5reader::app {

namespace {
constexpr double kOvaloidRadius = 1.5;   // ovaloid extent along the most-anisotropic axis
constexpr double kOvaloidFloor = 0.015;  // tiny floor so the pinch does not degenerate
constexpr double kArrowReach = 1.30;     // arrow tip beyond the surface extent
constexpr double kArrowMinLen = 0.70;    // floor so a near-iso axis still shows an arrow
constexpr double kArrowInnerGap = 0.10;  // tail offset from centre (+/- arrows don't coincide)
constexpr double kArrowWidth = 1.25;     // radial scale of the principal-axis arrows
constexpr double kAxisRgb[3][3] = {      // distinct per-axis colours (principal index 0/1/2);
    {0.96, 0.66, 0.16},  // 0 amber   mirrored in the Atom Info colour key so the arrows
    {0.18, 0.74, 0.74},  // 1 teal    stay decodable without in-scene labels
    {0.74, 0.36, 0.86},  // 2 violet
};
}  // namespace

TensorGlyphActor::TensorGlyphActor(vtkSmartPointer<vtkRenderer> sceneRenderer)
    : renderer_(std::move(sceneRenderer)) {}

TensorGlyphActor::~TensorGlyphActor() {
    if (renderer_) {
        if (glyphActor_) renderer_->RemoveActor(glyphActor_);
        for (auto& a : arrowActors_)
            if (a) renderer_->RemoveActor(a);
    }
}

void TensorGlyphActor::ensureActors() {
    if (actorsBuilt_) return;

    ovaloidSource_ = vtkSmartPointer<vtkSphereSource>::New();
    ovaloidSource_->SetRadius(1.0);
    ovaloidSource_->SetThetaResolution(64);
    ovaloidSource_->SetPhiResolution(64);

    glyphLocal_ = vtkSmartPointer<vtkPolyData>::New();

    glyphTransform_ = vtkSmartPointer<vtkTransform>::New();
    glyphFilter_ = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    glyphFilter_->SetInputData(glyphLocal_);
    glyphFilter_->SetTransform(glyphTransform_);

    // Re-derive normals AFTER the radial deform so lighting conveys the 3-D form
    // (deforming the sphere invalidates the source normals -> flat shading).
    glyphNormals_ = vtkSmartPointer<vtkPolyDataNormals>::New();
    glyphNormals_->SetInputConnection(glyphFilter_->GetOutputPort());
    glyphNormals_->ComputePointNormalsOn();
    glyphNormals_->ComputeCellNormalsOff();
    glyphNormals_->SplittingOff();
    glyphNormals_->ConsistencyOn();

    // Diverging sign map: below-iso (-1) red, iso (0) pale, above-iso (+1) blue.
    glyphLut_ = vtkSmartPointer<vtkColorTransferFunction>::New();
    glyphLut_->AddRGBPoint(-1.0, 0.82, 0.10, 0.10);
    glyphLut_->AddRGBPoint(0.0, 0.92, 0.92, 0.90);
    glyphLut_->AddRGBPoint(1.0, 0.09, 0.28, 0.78);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyphNormals_->GetOutputPort());
    mapper->SetLookupTable(glyphLut_);
    mapper->SetColorModeToMapScalars();
    mapper->SetScalarModeToUsePointData();
    mapper->InterpolateScalarsBeforeMappingOn();
    mapper->SetScalarRange(-1.0, 1.0);

    glyphActor_ = vtkSmartPointer<vtkActor>::New();
    glyphActor_->SetMapper(mapper);
    glyphActor_->GetProperty()->SetAmbient(0.20);
    glyphActor_->GetProperty()->SetDiffuse(0.85);
    glyphActor_->GetProperty()->SetSpecular(0.40);
    glyphActor_->GetProperty()->SetSpecularPower(30);
    glyphActor_->GetProperty()->SetOpacity(0.50);           // translucent: arrows lead, surface = context
    glyphActor_->GetProperty()->SetInterpolationToPhong();  // smooth shading, like the isosurfaces
    glyphActor_->SetVisibility(false);
    renderer_->AddActor(glyphActor_);  // MAIN renderer (depth-peeled)

    // Principal-axis arrows: a shared arrow source + mapper, six opaque actors
    // (ring-overlay proportions), each placed by a per-show UserMatrix and
    // coloured by principal-value index. Same depth-peeled renderer as the surface.
    arrowSource_ = vtkSmartPointer<vtkArrowSource>::New();
    arrowSource_->SetTipResolution(18);
    arrowSource_->SetShaftResolution(18);
    arrowSource_->SetTipLength(0.32);
    arrowSource_->SetTipRadius(0.12);
    arrowSource_->SetShaftRadius(0.042);
    arrowMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
    arrowMapper_->SetInputConnection(arrowSource_->GetOutputPort());
    for (auto& a : arrowActors_) {
        a = vtkSmartPointer<vtkActor>::New();
        a->SetMapper(arrowMapper_);
        a->GetProperty()->SetAmbient(0.30);
        a->GetProperty()->SetDiffuse(0.80);
        a->GetProperty()->SetSpecular(0.30);
        a->GetProperty()->SetSpecularPower(25);
        a->SetVisibility(false);
        renderer_->AddActor(a);
    }

    actorsBuilt_ = true;
}

void TensorGlyphActor::hideAll() {
    if (glyphActor_) glyphActor_->SetVisibility(false);
    for (auto& a : arrowActors_)
        if (a) a->SetVisibility(false);
}

void TensorGlyphActor::show(const model::Vec3& center,
                            const std::array<double, 3>& principalValues,
                            const model::Mat3& pasAxes,
                            double iso) {
    if (!renderer_) {
        clear();
        return;
    }
    ensureActors();

    std::array<double, 3> dev{
        principalValues[0] - iso,
        principalValues[1] - iso,
        principalValues[2] - iso,
    };
    double maxAbs = 0.0;
    for (double d : dev) maxAbs = std::max(maxAbs, std::abs(d));
    if (maxAbs < 1e-9) maxAbs = 1.0;

    // Ovaloid surface: a unit sphere deformed radially by r ~ |dev(n)| (the
    // deviation-from-isotropic surface), pinching toward the dev(n)=0 cone so the
    // above-/below-iso lobes separate. Local axes are the principal axes 0/1/2
    // directly; orientation is the transform. The scope limits the surface-build
    // locals (npts / signScalar / m).
    std::array<double, 3> surfaceExtent{};
    {
        ovaloidSource_->Update();
        glyphLocal_->DeepCopy(ovaloidSource_->GetOutput());
        const vtkIdType npts = glyphLocal_->GetNumberOfPoints();
        auto signScalar = vtkSmartPointer<vtkDoubleArray>::New();
        signScalar->SetName("dev_sign");
        signScalar->SetNumberOfComponents(1);
        signScalar->SetNumberOfTuples(npts);
        for (vtkIdType i = 0; i < npts; ++i) {
            double p[3];
            glyphLocal_->GetPoint(i, p);  // unit sphere: p is the direction n
            const double len = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            double nx = 0.0, ny = 0.0, nz = 1.0;
            if (len > 1e-9) { nx = p[0] / len; ny = p[1] / len; nz = p[2] / len; }
            const double devN = nx * nx * dev[0] + ny * ny * dev[1] + nz * nz * dev[2];
            const double r = kOvaloidRadius * std::max(kOvaloidFloor, std::abs(devN) / maxAbs);
            glyphLocal_->GetPoints()->SetPoint(i, nx * r, ny * r, nz * r);
            signScalar->SetValue(i, std::clamp(devN / maxAbs, -1.0, 1.0));
        }
        glyphLocal_->GetPoints()->Modified();
        glyphLocal_->GetPointData()->SetScalars(signScalar);

        auto m = vtkSmartPointer<vtkMatrix4x4>::New();
        m->Identity();  // unit eigenvector columns: rotate principal frame -> lab
        for (int col = 0; col < 3; ++col) {
            const model::Vec3 axisVec = pasAxes.col(col);
            for (int row = 0; row < 3; ++row)
                m->SetElement(row, col, axisVec[row]);
        }
        for (int row = 0; row < 3; ++row)
            m->SetElement(row, 3, center[row]);
        glyphTransform_->SetMatrix(m);
        for (int i = 0; i < 3; ++i)
            surfaceExtent[static_cast<std::size_t>(i)] =
                kOvaloidRadius * std::abs(dev[static_cast<std::size_t>(i)]) / maxAbs;
    }
    glyphFilter_->Modified();
    glyphActor_->SetVisibility(true);

    // Principal-axis arrows -- the descriptive element. Each PAS axis is drawn
    // double-headed (+/- director), index-coloured (amber/teal/violet for
    // index 0/1/2), length tracking the surface extent (with a floor) so the
    // arrowheads clear the translucent surface. The values + colour key are
    // shown in the Atom Info panel, not labelled here.
    for (int axis = 0; axis < 3; ++axis) {
        const std::size_t ai = static_cast<std::size_t>(axis);
        const model::Vec3 d = pasAxes.col(axis).normalized();
        const model::Vec3 ref = (std::abs(d[0]) < 0.9) ? model::Vec3(1, 0, 0) : model::Vec3(0, 1, 0);
        const model::Vec3 e1 = ref.cross(d).normalized();
        const model::Vec3 e2 = d.cross(e1);
        const double tipDist = std::max(kArrowMinLen, surfaceExtent[ai] * kArrowReach);
        const double* col = kAxisRgb[ai];

        for (int s = 0; s < 2; ++s) {
            const std::size_t slot = ai * 2 + static_cast<std::size_t>(s);
            const model::Vec3 dir = d * (s == 0 ? 1.0 : -1.0);
            const model::Vec3 tail = center + dir * kArrowInnerGap;
            const double shaftLen = std::max(0.05, tipDist - kArrowInnerGap);
            auto um = vtkSmartPointer<vtkMatrix4x4>::New();
            um->Identity();
            for (int r = 0; r < 3; ++r) {
                um->SetElement(r, 0, dir[r] * shaftLen);
                um->SetElement(r, 1, e1[r] * kArrowWidth);
                um->SetElement(r, 2, e2[r] * kArrowWidth);
                um->SetElement(r, 3, tail[r]);
            }
            arrowActors_[slot]->SetUserMatrix(um);
            arrowActors_[slot]->GetProperty()->SetColor(col[0], col[1], col[2]);
            arrowActors_[slot]->SetVisibility(true);
        }
    }

    active_ = true;
}

void TensorGlyphActor::clear() {
    hideAll();
    active_ = false;
}

void TensorGlyphActor::setVisible(bool on) {
    const bool v = on && active_;
    if (glyphActor_) glyphActor_->SetVisibility(v);
    for (auto& a : arrowActors_)
        if (a) a->SetVisibility(v);
}

}  // namespace h5reader::app

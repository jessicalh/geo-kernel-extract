#include "CsaTensorOverlay.h"

#include <vtkDoubleArray.h>
#include <vtkMatrix4x4.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkTextProperty.h>

#include <QByteArray>
#include <QFile>
#include <QLatin1String>
#include <QString>

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

namespace h5reader::app {

namespace {
constexpr double kBaseRadius = 1.3;   // Angstrom, the most-anisotropic axis extent
constexpr double kFloorFrac = 0.14;   // min radius fraction (near-iso axis thickness)
constexpr double kSharpGamma = 3.0;   // Kindlmann edge-sharpness exponent
constexpr double kRoundMin = 0.15;    // clamp roundness away from a degenerate box
constexpr double kLabelGap = 0.75;    // Angstrom beyond the principal extent
constexpr double kOvaloidRadius = 1.5;   // ovaloid extent along the most-anisotropic axis
constexpr double kOvaloidFloor = 0.015;  // tiny floor so the pinch does not degenerate
constexpr const char* kGreekFont = "C:/Windows/Fonts/arialbd.ttf";  // has Greek sigma

// "sigma_11" etc. as UTF-8 (Greek sigma U+03C3). Source stays ASCII.
QByteArray axisLabelText(int axis) {
    const QChar sigma(char16_t(0x03C3));
    const char* idx = (axis == 0) ? "11" : (axis == 1) ? "22" : "33";
    return (QString(sigma) + QLatin1String(idx)).toUtf8();
}
}  // namespace

CsaTensorOverlay::CsaTensorOverlay(vtkSmartPointer<vtkRenderer> sceneRenderer,
                                   vtkSmartPointer<vtkRenderer> hudRenderer,
                                   QObject* parent)
    : QObject(parent),
      renderer_(std::move(sceneRenderer)),
      hudRenderer_(std::move(hudRenderer)) {}

CsaTensorOverlay::~CsaTensorOverlay() {
    if (renderer_ && glyphActor_) renderer_->RemoveActor(glyphActor_);
    if (hudRenderer_) {
        for (auto& t : axisLabels_)
            if (t) hudRenderer_->RemoveViewProp(t);
        if (readout_) hudRenderer_->RemoveViewProp(readout_);
    }
}

void CsaTensorOverlay::ensureActors() {
    if (actorsBuilt_) return;

    glyphSource_ = vtkSmartPointer<vtkSuperquadricSource>::New();
    glyphSource_->SetThetaResolution(48);
    glyphSource_->SetPhiResolution(48);
    glyphSource_->ToroidalOff();
    glyphSource_->SetSize(1.0);
    glyphSource_->SetScale(1.0, 1.0, 1.0);

    ovaloidSource_ = vtkSmartPointer<vtkSphereSource>::New();
    ovaloidSource_->SetRadius(1.0);
    ovaloidSource_->SetThetaResolution(64);
    ovaloidSource_->SetPhiResolution(64);

    glyphLocal_ = vtkSmartPointer<vtkPolyData>::New();

    glyphTransform_ = vtkSmartPointer<vtkTransform>::New();
    glyphFilter_ = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    glyphFilter_->SetInputData(glyphLocal_);
    glyphFilter_->SetTransform(glyphTransform_);

    // Re-derive normals AFTER the deform/scale so lighting conveys the 3-D form
    // (the ovaloid radial deform invalidates the source sphere's normals, which is
    // what made the lobes read flat).
    glyphNormals_ = vtkSmartPointer<vtkPolyDataNormals>::New();
    glyphNormals_->SetInputConnection(glyphFilter_->GetOutputPort());
    glyphNormals_->ComputePointNormalsOn();
    glyphNormals_->ComputeCellNormalsOff();
    glyphNormals_->SplittingOff();
    glyphNormals_->ConsistencyOn();

    // Diverging sign map: deshielded (-1) red, iso (0) pale, shielded (+1) blue.
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
    glyphActor_->GetProperty()->SetSpecular(0.40);  // specular form-shading -> reads 3-D
    glyphActor_->GetProperty()->SetSpecularPower(30);
    glyphActor_->GetProperty()->SetOpacity(0.62);           // translucent: seamless via depth peeling
    glyphActor_->GetProperty()->SetInterpolationToPhong();  // smooth shading, like the isosurfaces
    glyphActor_->SetVisibility(false);
    renderer_->AddActor(glyphActor_);  // MAIN renderer (depth-peeled)

    for (auto& t : axisLabels_) {
        t = vtkSmartPointer<vtkBillboardTextActor3D>::New();
        vtkTextProperty* tp = t->GetTextProperty();
        tp->SetFontSize(22);
        tp->SetBold(true);
        tp->SetJustificationToCentered();
        tp->SetVerticalJustificationToCentered();
        tp->SetColor(0.97, 0.97, 0.92);
        tp->SetBackgroundColor(0.05, 0.05, 0.06);
        tp->SetBackgroundOpacity(0.5);
        tp->SetFrame(true);
        tp->SetFrameColor(0.3, 0.3, 0.3);
        if (QFile::exists(QString::fromLatin1(kGreekFont))) {
            tp->SetFontFamily(VTK_FONT_FILE);  // Greek sigma in the labels
            tp->SetFontFile(kGreekFont);
        }
        t->SetVisibility(false);
        hudRenderer_->AddViewProp(t);  // overlay renderer (always readable)
    }

    readout_ = vtkSmartPointer<vtkCornerAnnotation>::New();
    readout_->SetLinearFontScaleFactor(2);
    readout_->SetNonlinearFontScaleFactor(1);
    readout_->SetMaximumFontSize(20);
    vtkTextProperty* rp = readout_->GetTextProperty();
    rp->SetColor(0.97, 0.97, 0.88);
    rp->SetBold(true);
    rp->SetBackgroundColor(0.05, 0.05, 0.08);
    rp->SetBackgroundOpacity(0.55);
    readout_->SetVisibility(false);
    hudRenderer_->AddViewProp(readout_);  // overlay renderer (always readable)

    actorsBuilt_ = true;
}

void CsaTensorOverlay::hideAll() {
    if (glyphActor_) glyphActor_->SetVisibility(false);
    for (auto& t : axisLabels_)
        if (t) t->SetVisibility(false);
    if (readout_) readout_->SetVisibility(false);
}

void CsaTensorOverlay::show(const model::Vec3& atomPos,
                            const model::CsaShape& shape,
                            const std::optional<model::Mat3>& /*molecularAxes*/) {
    if (!renderer_ || !hudRenderer_ || !shape.valid) {
        clear();
        return;
    }
    ensureActors();

    const double iso = shape.sigma_iso;
    std::array<double, 3> dev{
        shape.principal_values[0] - iso,
        shape.principal_values[1] - iso,
        shape.principal_values[2] - iso,
    };
    double maxAbs = 0.0;
    for (double d : dev) maxAbs = std::max(maxAbs, std::abs(d));
    if (maxAbs < 1e-9) maxAbs = 1.0;

    // Per-principal-axis surface extent along PAS axis i, for label placement.
    std::array<double, 3> surfaceExtent{};

    if (style_ == GlyphStyle::Ovaloid) {
        // TensorView ovaloid: a unit sphere deformed radially by r ~ |dev(n)| (the
        // anisotropy shielding surface). It pinches toward the dev(n) = 0 cone, so
        // the shielded (dev>0) and deshielded (dev<0) lobes separate -- sign by
        // SHAPE as well as colour, the NMR-community "more correct" form. Local
        // axes are the principal axes 0/1/2 directly; orientation is the transform.
        ovaloidSource_->Update();
        glyphLocal_->DeepCopy(ovaloidSource_->GetOutput());
        const vtkIdType npts = glyphLocal_->GetNumberOfPoints();
        auto signScalar = vtkSmartPointer<vtkDoubleArray>::New();
        signScalar->SetName("shielding_sign");
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
            const model::Vec3 axisVec = shape.pas_axes.col(col);
            for (int row = 0; row < 3; ++row)
                m->SetElement(row, col, axisVec[row]);
        }
        for (int row = 0; row < 3; ++row)
            m->SetElement(row, 3, atomPos[row]);
        glyphTransform_->SetMatrix(m);
        for (int i = 0; i < 3; ++i)
            surfaceExtent[static_cast<std::size_t>(i)] =
                kOvaloidRadius * std::abs(dev[static_cast<std::size_t>(i)]) / maxAbs;
    } else {
        // Superquadric (Kindlmann) / ellipsoid: a unit superquadric scaled by the
        // |dev| radii and oriented onto the PAS via the transform. The shape is
        // symmetric in size; the SIGN of each direction is carried by colour.
        std::array<double, 3> radii{};
        for (int i = 0; i < 3; ++i)
            radii[i] = kBaseRadius * (kFloorFrac + (1.0 - kFloorFrac) * std::abs(dev[i]) / maxAbs);

        // Kindlmann linear/planar anisotropy on |dev| (sorted hi >= mid >= lo)
        // drives the edge sharpness: linear -> cigar along the largest-|dev| axis;
        // planar -> disc whose normal is the smallest-|dev| axis; isotropic ->
        // sphere -- the disambiguation a plain ellipsoid cannot give.
        std::array<int, 3> ord{0, 1, 2};
        std::sort(ord.begin(), ord.end(),
                  [&](int a, int b) { return std::abs(dev[a]) > std::abs(dev[b]); });
        const int hi = ord[0], mid = ord[1], lo = ord[2];
        const double sum = std::abs(dev[hi]) + std::abs(dev[mid]) + std::abs(dev[lo]);
        const double cl = (sum > 1e-9) ? (std::abs(dev[hi]) - std::abs(dev[mid])) / sum : 0.0;
        const double cp = (sum > 1e-9) ? 2.0 * (std::abs(dev[mid]) - std::abs(dev[lo])) / sum : 0.0;

        auto clampR = [](double r) { return std::clamp(r, kRoundMin, 1.0); };
        int uniqueAxis, axisA, axisB;
        double thetaR, phiR;
        if (cl >= cp) {  // linear -> cigar along the unique (largest-|dev|) axis
            uniqueAxis = hi; axisA = mid; axisB = lo;
            phiR = clampR(std::pow(1.0 - cl, kSharpGamma));
            thetaR = clampR(std::pow(1.0 - cp, kSharpGamma));
        } else {  // planar -> disc whose normal is the smallest-|dev| axis
            uniqueAxis = lo; axisA = hi; axisB = mid;
            phiR = clampR(std::pow(1.0 - cp, kSharpGamma));
            thetaR = clampR(std::pow(1.0 - cl, kSharpGamma));
        }
        if (style_ == GlyphStyle::Ellipsoid) {  // roundness-1 superquadric = ellipsoid
            thetaR = 1.0;
            phiR = 1.0;
        }

        glyphSource_->SetThetaRoundness(thetaR);
        glyphSource_->SetPhiRoundness(phiR);
        glyphSource_->Update();

        // Deep-copy the unit superquadric and attach the per-vertex sign scalar
        // (signed directional deviation / maxAbs in [-1,1]). Local axes map
        // X->axisA, Y->uniqueAxis (the vtk superquadric symmetry axis), Z->axisB.
        glyphLocal_->DeepCopy(glyphSource_->GetOutput());
        const vtkIdType npts = glyphLocal_->GetNumberOfPoints();
        auto signScalar = vtkSmartPointer<vtkDoubleArray>::New();
        signScalar->SetName("shielding_sign");
        signScalar->SetNumberOfComponents(1);
        signScalar->SetNumberOfTuples(npts);
        for (vtkIdType i = 0; i < npts; ++i) {
            double p[3];
            glyphLocal_->GetPoint(i, p);
            const double len = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            double s = 0.0;
            if (len > 1e-9) {
                const double nx = p[0] / len, ny = p[1] / len, nz = p[2] / len;
                s = (nx * nx * dev[axisA] + ny * ny * dev[uniqueAxis] + nz * nz * dev[axisB]) / maxAbs;
            }
            signScalar->SetValue(i, std::clamp(s, -1.0, 1.0));
        }
        glyphLocal_->GetPointData()->SetScalars(signScalar);

        // Transform: local X/Y/Z -> (axisA, uniqueAxis, axisB) eigenvectors scaled
        // by their radii; translate to the atom.
        const std::array<int, 3> localToPrincipal{axisA, uniqueAxis, axisB};
        auto m = vtkSmartPointer<vtkMatrix4x4>::New();
        m->Identity();
        for (int col = 0; col < 3; ++col) {
            const int pax = localToPrincipal[static_cast<std::size_t>(col)];
            const model::Vec3 axisVec = shape.pas_axes.col(pax) * radii[static_cast<std::size_t>(pax)];
            for (int row = 0; row < 3; ++row)
                m->SetElement(row, col, axisVec[row]);
        }
        for (int row = 0; row < 3; ++row)
            m->SetElement(row, 3, atomPos[row]);
        glyphTransform_->SetMatrix(m);
        for (int i = 0; i < 3; ++i)
            surfaceExtent[static_cast<std::size_t>(i)] = radii[static_cast<std::size_t>(i)];
    }
    glyphFilter_->Modified();
    glyphActor_->SetVisibility(true);

    // sigma_11/22/33 labels at each principal extent, tinted to the surface sign.
    for (int axis = 0; axis < 3; ++axis) {
        const model::Vec3 d = shape.pas_axes.col(axis).normalized();
        const model::Vec3 tip = atomPos + d * (surfaceExtent[static_cast<std::size_t>(axis)] + kLabelGap);
        const QByteArray lbl = axisLabelText(axis);
        axisLabels_[static_cast<std::size_t>(axis)]->SetInput(lbl.constData());
        axisLabels_[static_cast<std::size_t>(axis)]->SetPosition(tip[0], tip[1], tip[2]);
        if (dev[axis] > 0.0)
            axisLabels_[static_cast<std::size_t>(axis)]->GetTextProperty()->SetColor(0.66, 0.78, 0.98);
        else
            axisLabels_[static_cast<std::size_t>(axis)]->GetTextProperty()->SetColor(0.98, 0.70, 0.66);
        axisLabels_[static_cast<std::size_t>(axis)]->SetVisibility(true);
    }

    // Corner readout -- the shape numbers (absolute shielding, NOT shift) + key.
    // ASCII words (the embedded VTK font lacks Greek), so this always renders.
    QString txt = QStringLiteral("CSA tensor  (DFT, viewer-derived)\n");
    txt += QStringLiteral("iso  = %1 ppm   (absolute shielding)\n").arg(shape.sigma_iso, 0, 'f', 1);
    txt += QStringLiteral("span = %1 ppm\n").arg(shape.span, 0, 'f', 1);
    txt += QStringLiteral("skew = %1\n").arg(shape.skew, 0, 'f', 2);
    txt += QStringLiteral("eta  = %1\n").arg(shape.eta, 0, 'f', 2);
    txt += QStringLiteral("blue = shielded    red = deshielded");
    const QByteArray readoutBytes = txt.toUtf8();
    readout_->SetText(vtkCornerAnnotation::UpperLeft, readoutBytes.constData());
    readout_->SetVisibility(true);

    active_ = true;
}

void CsaTensorOverlay::clear() {
    hideAll();
    active_ = false;
}

void CsaTensorOverlay::setVisible(bool on) {
    const bool v = on && active_;
    if (glyphActor_) glyphActor_->SetVisibility(v);
    for (auto& t : axisLabels_)
        if (t) t->SetVisibility(v);
    if (readout_) readout_->SetVisibility(v);
}

}  // namespace h5reader::app

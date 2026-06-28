#include "QtFieldGridOverlay.h"

#include "../calculators/QtBiotSavartCalc.h"
#include "../calculators/QtHaighMallionCalc.h"
#include "../calculators/QtPhysicalConstants.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/ConformationGeometry.h"
#include "../model/QtResidue.h"

#include <QElapsedTimer>
#include <QLoggingCategory>

#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cField, "h5reader.overlay.field")

using calculators::FIELD_GRID_DIM;
using calculators::FIELD_GRID_EXTENT_A;

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kMagicConeTan = 1.41421356237309504880; // tan(acos(1/sqrt(3)))
constexpr int kNullConeSegments = 96;

}  // namespace

QtFieldGridOverlay::QtFieldGridOverlay(
    vtkSmartPointer<vtkRenderer>                  renderer,
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
    QObject* parent)
    : QObject(parent),
      renderer_(std::move(renderer)),
      renderWindow_(std::move(renderWindow))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtFieldGridOverlay"));
}

QtFieldGridOverlay::~QtFieldGridOverlay() {
    for (const auto& rg : rings_) {
        if (rg.actorShielded)   renderer_->RemoveActor(rg.actorShielded);
        if (rg.actorDeshielded) renderer_->RemoveActor(rg.actorDeshielded);
        if (rg.actorNullCone)   renderer_->RemoveActor(rg.actorNullCone);
    }
}

void QtFieldGridOverlay::Build(const model::QtProtein& protein,
                                model::Conformation&    conformation) {
    ASSERT_THREAD(this);

    if (protein_ == &protein && conformation_ == &conformation && !rings_.empty())
        return;

    for (const auto& rg : rings_) {
        if (rg.actorShielded)   renderer_->RemoveActor(rg.actorShielded);
        if (rg.actorDeshielded) renderer_->RemoveActor(rg.actorDeshielded);
        if (rg.actorNullCone)   renderer_->RemoveActor(rg.actorNullCone);
    }
    rings_.clear();

    protein_      = &protein;
    conformation_ = &conformation;

    const size_t n_rings = protein.ringCount();
    rings_.resize(n_rings);

    const int dim = FIELD_GRID_DIM;
    const int nPoints = dim * dim * dim;
    const double spacing = 2.0 * FIELD_GRID_EXTENT_A / (dim - 1);

    for (size_t ri = 0; ri < n_rings; ++ri) {
        auto& rg = rings_[ri];

        rg.scalars = vtkSmartPointer<vtkFloatArray>::New();
        rg.scalars->SetName("T0");
        rg.scalars->SetNumberOfTuples(nPoints);
        // Initial fill with zero — RecomputeRingScalars overwrites on setFrame.
        for (int i = 0; i < nPoints; ++i) rg.scalars->SetValue(i, 0.0f);

        rg.imageData = vtkSmartPointer<vtkImageData>::New();
        rg.imageData->SetDimensions(dim, dim, dim);
        rg.imageData->SetSpacing(spacing, spacing, spacing);
        rg.imageData->SetOrigin(0.0, 0.0, 0.0);    // updated per frame
        rg.imageData->GetPointData()->SetScalars(rg.scalars);

        rg.producer = vtkSmartPointer<vtkTrivialProducer>::New();
        rg.producer->SetOutput(rg.imageData);

        rg.contourShielded = vtkSmartPointer<vtkContourFilter>::New();
        rg.contourShielded->SetInputConnection(rg.producer->GetOutputPort());
        rg.contourShielded->SetValue(0, -thresholdPpm_);

        rg.contourDeshielded = vtkSmartPointer<vtkContourFilter>::New();
        rg.contourDeshielded->SetInputConnection(rg.producer->GetOutputPort());
        rg.contourDeshielded->SetValue(0, +thresholdPpm_);

        auto mapperS = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapperS->SetInputConnection(rg.contourShielded->GetOutputPort());
        mapperS->ScalarVisibilityOff();
        rg.actorShielded = vtkSmartPointer<vtkActor>::New();
        rg.actorShielded->SetMapper(mapperS);

        auto mapperD = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapperD->SetInputConnection(rg.contourDeshielded->GetOutputPort());
        mapperD->ScalarVisibilityOff();
        rg.actorDeshielded = vtkSmartPointer<vtkActor>::New();
        rg.actorDeshielded->SetMapper(mapperD);

        rg.nullConePoints = vtkSmartPointer<vtkPoints>::New();
        rg.nullConePolys = vtkSmartPointer<vtkCellArray>::New();
        rg.nullConePoly = vtkSmartPointer<vtkPolyData>::New();
        rg.nullConePoly->SetPoints(rg.nullConePoints);
        rg.nullConePoly->SetPolys(rg.nullConePolys);

        auto mapperN = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapperN->SetInputData(rg.nullConePoly);
        mapperN->ScalarVisibilityOff();
        rg.actorNullCone = vtkSmartPointer<vtkActor>::New();
        rg.actorNullCone->SetMapper(mapperN);

        ApplyActorStyling(ri, rg);

        renderer_->AddActor(rg.actorShielded);
        renderer_->AddActor(rg.actorDeshielded);
        renderer_->AddActor(rg.actorNullCone);
    }

    qCInfo(cField).noquote()
        << "Built field grid overlay |"
        << "rings=" << n_rings
        << "| grid=" << dim << "³"
        << "| extent=" << FIELD_GRID_EXTENT_A << "Å";

    setFrame(0);
}

bool QtFieldGridOverlay::RingIsAllowed(size_t ringIdx) const {
    return !visibleRing_ || *visibleRing_ == ringIdx;
}

void QtFieldGridOverlay::RecomputeVisibleRings() {
    if (!visible_ || !protein_ || !conformation_) return;
    for (size_t ri = 0; ri < rings_.size(); ++ri) {
        if (RingIsAllowed(ri))
            RecomputeRingScalars(ri, currentFrame_);
    }
}

void QtFieldGridOverlay::UpdateNullConeGeometry(size_t ringIdx,
                                                RingGrid& rg,
                                                const model::RingGeometry& geo) {
    if (!rg.actorNullCone || !rg.nullConePoints || !rg.nullConePolys || !rg.nullConePoly)
        return;
    if (!nullConeVisible_ || !RingIsAllowed(ringIdx) || geo.radius < 1e-6) {
        rg.actorNullCone->SetVisibility(0);
        return;
    }

    const auto basis = model::OrthoBasisFromNormal(geo.normal);
    const model::Vec3& u = basis.u;
    const model::Vec3& v = basis.v;
    const model::Vec3& n = basis.n;

    const double z0 = std::clamp(geo.radius * 0.22, 0.18, 0.45);
    const double z1 = std::max(z0 + 0.05, nullConeLengthA_);

    rg.nullConePoints->Reset();
    rg.nullConePolys->Reset();

    for (int side = 0; side < 2; ++side) {
        const double sign = side == 0 ? 1.0 : -1.0;
        const vtkIdType base = rg.nullConePoints->GetNumberOfPoints();
        for (int i = 0; i < kNullConeSegments; ++i) {
            const double theta = 2.0 * kPi * static_cast<double>(i)
                               / static_cast<double>(kNullConeSegments);
            const model::Vec3 radial = u * std::cos(theta) + v * std::sin(theta);
            const model::Vec3 inner = geo.center + n * (sign * z0)
                                    + radial * (z0 * kMagicConeTan);
            const model::Vec3 outer = geo.center + n * (sign * z1)
                                    + radial * (z1 * kMagicConeTan);
            rg.nullConePoints->InsertNextPoint(inner.x(), inner.y(), inner.z());
            rg.nullConePoints->InsertNextPoint(outer.x(), outer.y(), outer.z());
        }
        for (int i = 0; i < kNullConeSegments; ++i) {
            const vtkIdType i0 = base + static_cast<vtkIdType>(2 * i);
            const vtkIdType o0 = i0 + 1;
            const vtkIdType i1 = base + static_cast<vtkIdType>(2 * ((i + 1) % kNullConeSegments));
            const vtkIdType o1 = i1 + 1;
            rg.nullConePolys->InsertNextCell(4);
            rg.nullConePolys->InsertCellPoint(i0);
            rg.nullConePolys->InsertCellPoint(i1);
            rg.nullConePolys->InsertCellPoint(o1);
            rg.nullConePolys->InsertCellPoint(o0);
        }
    }

    rg.nullConePoints->Modified();
    rg.nullConePolys->Modified();
    rg.nullConePoly->Modified();
    rg.actorNullCone->SetVisibility(1);
}

void QtFieldGridOverlay::UpdateVisibleNullCones() {
    if (!nullConeVisible_ || !protein_ || !conformation_) return;
    for (size_t ri = 0; ri < rings_.size(); ++ri) {
        const auto geo = model::RingGeometryAt(*conformation_, ri, static_cast<size_t>(currentFrame_));
        UpdateNullConeGeometry(ri, rings_[ri], geo);
    }
}

void QtFieldGridOverlay::ApplyActorStyling(size_t ringIdx, RingGrid& rg) {
    // Sky blue for shielded (T0 < -threshold): atoms above/below the ring
    // where ring-current shielding is strongest. Coral for deshielded
    // (T0 > +threshold): in-plane where the current anti-shields.
    // Matches the library viewer's FieldGridOverlay colour choices.
    const bool showRing = visible_ && RingIsAllowed(ringIdx);
    rg.actorShielded->GetProperty()->SetColor(0.50, 0.70, 0.95);
    rg.actorShielded->GetProperty()->SetOpacity(opacity_);
    rg.actorShielded->GetProperty()->SetInterpolationToPhong();
    rg.actorShielded->SetForceTranslucent(true);
    rg.actorShielded->SetVisibility(showRing && shieldedVisible_ ? 1 : 0);

    rg.actorDeshielded->GetProperty()->SetColor(0.95, 0.55, 0.45);
    rg.actorDeshielded->GetProperty()->SetOpacity(opacity_);
    rg.actorDeshielded->GetProperty()->SetInterpolationToPhong();
    rg.actorDeshielded->SetForceTranslucent(true);
    rg.actorDeshielded->SetVisibility(showRing && deshieldedVisible_ ? 1 : 0);

    rg.actorNullCone->GetProperty()->SetColor(0.72, 0.76, 0.82);
    rg.actorNullCone->GetProperty()->SetOpacity(nullConeOpacity_);
    rg.actorNullCone->GetProperty()->LightingOff();
    rg.actorNullCone->GetProperty()->BackfaceCullingOff();
    rg.actorNullCone->SetForceTranslucent(true);
    rg.actorNullCone->SetVisibility(nullConeVisible_ && RingIsAllowed(ringIdx) ? 1 : 0);
}

void QtFieldGridOverlay::RecomputeRingScalars(size_t ri, int t) {
    auto& rg = rings_[ri];
    if (!rg.imageData) return;

    const auto& ring = protein_->ring(ri);
    const auto geo = model::RingGeometryAt(*conformation_, ri, static_cast<size_t>(t));
    if (geo.radius < 1e-6) {
        // No valid geometry (e.g. ring_geometry slab missing or loader
        // flagged the layout invalid). Zero the scalars and let the
        // contour filter produce nothing; actors stay hidden.
        for (int i = 0; i < rg.scalars->GetNumberOfTuples(); ++i)
            rg.scalars->SetValue(i, 0.0f);
        rg.scalars->Modified();
        return;
    }

    const auto vertices = model::RingVertices(*conformation_, ri, static_cast<size_t>(t));
    if (vertices.size() < 3) {
        rg.actorShielded->SetVisibility(0);
        rg.actorDeshielded->SetVisibility(0);
        return;
    }

    const int dim = FIELD_GRID_DIM;
    const double extent = FIELD_GRID_EXTENT_A;
    const double spacing = 2.0 * extent / (dim - 1);

    const double origin[3] = {
        geo.center.x() - extent,
        geo.center.y() - extent,
        geo.center.z() - extent
    };
    rg.imageData->SetOrigin(origin[0], origin[1], origin[2]);
    rg.imageData->SetSpacing(spacing, spacing, spacing);

    const double intensityNA = ring.LiteratureIntensity();
    const double lobeOffsetA = ring.JohnsonBoveyLobeOffset();

    // Evaluate kernel at every grid point. PointInValidRange inside
    // each calculator handles "too close" / "inside ring" / "too far".
    for (int iz = 0; iz < dim; ++iz) {
        for (int iy = 0; iy < dim; ++iy) {
            for (int ix = 0; ix < dim; ++ix) {
                const model::Vec3 p(
                    origin[0] + ix * spacing,
                    origin[1] + iy * spacing,
                    origin[2] + iz * spacing);

                double t0 = 0.0;
                switch (mode_) {
                    case FieldGridMode::BiotSavart: {
                        const auto st = calculators::EvaluateShielding(
                            p, geo, vertices, lobeOffsetA, intensityNA);
                        t0 = st.T0;
                        break;
                    }
                    case FieldGridMode::HaighMallion: {
                        const auto st = calculators::EvaluateShielding(
                            p, geo, vertices, intensityNA);
                        t0 = st.T0;
                        break;
                    }
                    case FieldGridMode::Sum: {
                        const auto stBS = calculators::EvaluateShielding(
                            p, geo, vertices, lobeOffsetA, intensityNA);
                        const auto stHM = calculators::EvaluateShielding(
                            p, geo, vertices, intensityNA);
                        t0 = stBS.T0 + stHM.T0;
                        break;
                    }
                }

                // Physics-preserving Gaussian radial taper: keep the true
                // kernel value near the ring, smoothly kill the 1/r³ tail so
                // the lobe has a bounded, tunable REACH (extent σ). gaussianPeak_
                // is an amplitude gain (1.0 = real ppm). This is what bounds the
                // field so the isosurface closes; the boundary clamp below is a
                // belt-and-suspenders guard for very large σ.
                const double rDist  = (p - geo.center).norm();
                const double window = std::exp(
                    -(rDist * rDist) /
                    (2.0 * gaussianExtentA_ * gaussianExtentA_));
                t0 *= gaussianPeak_ * window;

                const int idx = ix + iy * dim + iz * dim * dim;
                // Guaranteed closure: force the outermost grid shell to 0
                // (sub-threshold) so vtkContourFilter always emits a CLOSED
                // isosurface inside the box. Without this, any lobe still above
                // |threshold| at the box face exits the volume as an open/torn
                // surface — the "sheared butterfly" artifact. With the box
                // (FIELD_GRID_EXTENT_A) sized so the field decays below the
                // default threshold before the boundary, this normally clamps
                // only already-sub-threshold voxels; it also makes lowering the
                // threshold at runtime degrade gracefully (a cap at the box)
                // instead of tearing.
                const bool onBoundary =
                    (ix == 0 || ix == dim - 1 ||
                     iy == 0 || iy == dim - 1 ||
                     iz == 0 || iz == dim - 1);
                rg.scalars->SetValue(idx,
                    onBoundary ? 0.0f
                               : (std::isfinite(t0) ? static_cast<float>(t0) : 0.0f));
            }
        }
    }
    rg.scalars->Modified();
    rg.imageData->Modified();
}

void QtFieldGridOverlay::setFrame(int t) {
    ASSERT_THREAD(this);
    if (!protein_ || !conformation_) return;
    if (t < 0 || static_cast<size_t>(t) >= conformation_->frameCount()) return;
    currentFrame_ = t;
    if (nullConeVisible_) {
        for (size_t ri = 0; ri < rings_.size(); ++ri) {
            const auto geo = model::RingGeometryAt(*conformation_, ri, static_cast<size_t>(t));
            UpdateNullConeGeometry(ri, rings_[ri], geo);
        }
    }
    if (!visible_) return;   // skip the expensive re-eval when off

    // Mirror QtBFieldStreamOverlay's two-stage timing — eval (kernel
    // grid fill) and pipe (forced contour-filter Update).
    QElapsedTimer evalT;
    evalT.start();
    for (size_t ri = 0; ri < rings_.size(); ++ri) {
        if (RingIsAllowed(ri))
            RecomputeRingScalars(ri, t);
    }
    const qint64 evalMs = evalT.elapsed();

    QElapsedTimer pipeT;
    pipeT.start();
    for (size_t ri = 0; ri < rings_.size(); ++ri) {
        if (!RingIsAllowed(ri)) continue;
        auto& rg = rings_[ri];
        if (visible_ && shieldedVisible_   && rg.contourShielded)
            rg.contourShielded->Update();
        if (visible_ && deshieldedVisible_ && rg.contourDeshielded)
            rg.contourDeshielded->Update();
    }
    const qint64 pipeMs = pipeT.elapsed();

    qCDebug(cField).noquote()
        << "frame" << t
        << "|" << rings_.size() << "rings"
        << "| eval=" << evalMs << "ms"
        << "| pipe=" << pipeMs << "ms"
        << "| total=" << (evalMs + pipeMs) << "ms";
}

void QtFieldGridOverlay::setMode(FieldGridMode mode) {
    ASSERT_THREAD(this);
    if (mode == mode_) return;
    mode_ = mode;
    // Force re-eval at the current visible state.
    if (visible_) {
        for (size_t ri = 0; ri < rings_.size(); ++ri) {
            if (RingIsAllowed(ri))
                RecomputeRingScalars(ri, currentFrame_);
        }
    }
}

void QtFieldGridOverlay::UpdateThresholds() {
    for (auto& rg : rings_) {
        if (rg.contourShielded)   rg.contourShielded->SetValue(0, -thresholdPpm_);
        if (rg.contourDeshielded) rg.contourDeshielded->SetValue(0, +thresholdPpm_);
    }
}

void QtFieldGridOverlay::setThresholdPpm(double threshold) {
    ASSERT_THREAD(this);
    thresholdPpm_ = std::max(0.0, threshold);
    UpdateThresholds();
}

void QtFieldGridOverlay::setGaussianExtent(double sigmaA) {
    ASSERT_THREAD(this);
    gaussianExtentA_ = std::max(0.1, sigmaA);   // guard div-by-zero in the window
    // The taper changes the scalar field itself, so re-evaluate the grid at the
    // current frame (cheap one-exp-per-point added on top of the kernel). The
    // contour filters rerun on the next render via Modified().
    RecomputeVisibleRings();
}

void QtFieldGridOverlay::setGaussianPeak(double amplitude) {
    ASSERT_THREAD(this);
    gaussianPeak_ = std::max(0.0, amplitude);
    RecomputeVisibleRings();
}

void QtFieldGridOverlay::setOpacity(double opacity) {
    ASSERT_THREAD(this);
    opacity_ = std::clamp(opacity, 0.0, 1.0);
    for (size_t ri = 0; ri < rings_.size(); ++ri)
        ApplyActorStyling(ri, rings_[ri]);
}

void QtFieldGridOverlay::setVisible(bool visible) {
    ASSERT_THREAD(this);
    const bool wasVisible = visible_;
    visible_ = visible;
    for (size_t ri = 0; ri < rings_.size(); ++ri)
        ApplyActorStyling(ri, rings_[ri]);
    // When turning on, re-eval at the current frame — we skipped while off.
    if (visible && !wasVisible && protein_ && conformation_) {
        // MoleculeScene will drive the next setFrame; nothing to do here.
    }
}

void QtFieldGridOverlay::setShieldedVisible(bool visible) {
    ASSERT_THREAD(this);
    shieldedVisible_ = visible;
    for (size_t ri = 0; ri < rings_.size(); ++ri)
        ApplyActorStyling(ri, rings_[ri]);
}

void QtFieldGridOverlay::setDeshieldedVisible(bool visible) {
    ASSERT_THREAD(this);
    deshieldedVisible_ = visible;
    for (size_t ri = 0; ri < rings_.size(); ++ri)
        ApplyActorStyling(ri, rings_[ri]);
}

void QtFieldGridOverlay::setNullConeVisible(bool visible) {
    ASSERT_THREAD(this);
    nullConeVisible_ = visible;
    for (size_t ri = 0; ri < rings_.size(); ++ri)
        ApplyActorStyling(ri, rings_[ri]);
    if (visible)
        UpdateVisibleNullCones();
}

void QtFieldGridOverlay::setNullConeOpacity(double opacity) {
    ASSERT_THREAD(this);
    nullConeOpacity_ = std::clamp(opacity, 0.0, 1.0);
    for (size_t ri = 0; ri < rings_.size(); ++ri)
        ApplyActorStyling(ri, rings_[ri]);
}

void QtFieldGridOverlay::setNullConeLength(double lengthA) {
    ASSERT_THREAD(this);
    nullConeLengthA_ = std::clamp(lengthA, 0.25, FIELD_GRID_EXTENT_A);
    UpdateVisibleNullCones();
}

void QtFieldGridOverlay::setVisibleRing(std::optional<size_t> ringIdx) {
    ASSERT_THREAD(this);
    if (ringIdx && *ringIdx >= rings_.size())
        ringIdx.reset();
    if (visibleRing_ == ringIdx)
        return;
    visibleRing_ = ringIdx;
    for (size_t ri = 0; ri < rings_.size(); ++ri)
        ApplyActorStyling(ri, rings_[ri]);
    RecomputeVisibleRings();
    UpdateVisibleNullCones();
}

}  // namespace h5reader::app

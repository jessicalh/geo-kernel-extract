#include "QtFieldGridOverlay.h"

#include "../calculators/QtBiotSavartCalc.h"
#include "../calculators/QtHaighMallionCalc.h"
#include "../calculators/QtPhysicalConstants.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/QtFrame.h"
#include "../model/QtResidue.h"

#include <QElapsedTimer>
#include <QLoggingCategory>

#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <limits>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cField, "h5reader.overlay.field")

using calculators::FIELD_GRID_DIM;
using calculators::FIELD_GRID_EXTENT_A;

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
    }
}

void QtFieldGridOverlay::Build(const model::QtProtein&      protein,
                                const model::QtConformation& conformation) {
    ASSERT_THREAD(this);

    if (protein_ == &protein && conformation_ == &conformation && !rings_.empty())
        return;

    for (const auto& rg : rings_) {
        if (rg.actorShielded)   renderer_->RemoveActor(rg.actorShielded);
        if (rg.actorDeshielded) renderer_->RemoveActor(rg.actorDeshielded);
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

        ApplyActorStyling(rg);

        renderer_->AddActor(rg.actorShielded);
        renderer_->AddActor(rg.actorDeshielded);
    }

    qCInfo(cField).noquote()
        << "Built field grid overlay |"
        << "rings=" << n_rings
        << "| grid=" << dim << "³"
        << "| extent=" << FIELD_GRID_EXTENT_A << "Å";

    setFrame(0);
}

void QtFieldGridOverlay::ApplyActorStyling(RingGrid& rg) {
    // Sky blue for shielded (T0 < -threshold): atoms above/below the ring
    // where ring-current shielding is strongest. Coral for deshielded
    // (T0 > +threshold): in-plane where the current anti-shields.
    // Matches the library viewer's FieldGridOverlay colour choices.
    rg.actorShielded->GetProperty()->SetColor(0.50, 0.70, 0.95);
    rg.actorShielded->GetProperty()->SetOpacity(opacity_);
    rg.actorShielded->GetProperty()->SetInterpolationToPhong();
    rg.actorShielded->SetForceTranslucent(true);
    rg.actorShielded->SetVisibility(visible_ && shieldedVisible_ ? 1 : 0);

    rg.actorDeshielded->GetProperty()->SetColor(0.95, 0.55, 0.45);
    rg.actorDeshielded->GetProperty()->SetOpacity(opacity_);
    rg.actorDeshielded->GetProperty()->SetInterpolationToPhong();
    rg.actorDeshielded->SetForceTranslucent(true);
    rg.actorDeshielded->SetVisibility(visible_ && deshieldedVisible_ ? 1 : 0);
}

void QtFieldGridOverlay::RecomputeRingScalars(size_t ri, int t) {
    auto& rg = rings_[ri];
    if (!rg.imageData) return;

    const auto& ring = protein_->ring(ri);
    const auto& frame = conformation_->frame(static_cast<size_t>(t));
    const auto geo = frame.ringGeometry(ri);
    if (geo.radius < 1e-6) {
        // No valid geometry (e.g. ring_geometry slab missing or loader
        // flagged the layout invalid). Zero the scalars and let the
        // contour filter produce nothing; actors stay hidden.
        for (int i = 0; i < rg.scalars->GetNumberOfTuples(); ++i)
            rg.scalars->SetValue(i, 0.0f);
        rg.scalars->Modified();
        return;
    }

    const auto vertices = frame.ringVertices(ri);
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
    const double lobeOffsetA = ring.JBLobeOffset();

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

                const int idx = ix + iy * dim + iz * dim * dim;
                rg.scalars->SetValue(idx,
                    std::isfinite(t0) ? static_cast<float>(t0) : 0.0f);
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
    if (!visible_) return;   // skip the expensive re-eval when off

    QElapsedTimer timer;
    timer.start();
    for (size_t ri = 0; ri < rings_.size(); ++ri) {
        RecomputeRingScalars(ri, t);
    }
    qCDebug(cField).noquote()
        << "frame" << t << "|" << rings_.size() << "rings |"
        << timer.elapsed() << "ms";
}

void QtFieldGridOverlay::setMode(FieldGridMode mode) {
    ASSERT_THREAD(this);
    if (mode == mode_) return;
    mode_ = mode;
    // Force re-eval at the current visible state.
    if (visible_) {
        for (size_t ri = 0; ri < rings_.size(); ++ri) {
            RecomputeRingScalars(ri, /*t*/ 0);
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

void QtFieldGridOverlay::setOpacity(double opacity) {
    ASSERT_THREAD(this);
    opacity_ = std::clamp(opacity, 0.0, 1.0);
    for (auto& rg : rings_) ApplyActorStyling(rg);
}

void QtFieldGridOverlay::setVisible(bool visible) {
    ASSERT_THREAD(this);
    const bool wasVisible = visible_;
    visible_ = visible;
    for (auto& rg : rings_) ApplyActorStyling(rg);
    // When turning on, re-eval at the current frame — we skipped while off.
    if (visible && !wasVisible && protein_ && conformation_) {
        // MoleculeScene will drive the next setFrame; nothing to do here.
    }
}

void QtFieldGridOverlay::setShieldedVisible(bool visible) {
    ASSERT_THREAD(this);
    shieldedVisible_ = visible;
    for (auto& rg : rings_) ApplyActorStyling(rg);
}

void QtFieldGridOverlay::setDeshieldedVisible(bool visible) {
    ASSERT_THREAD(this);
    deshieldedVisible_ = visible;
    for (auto& rg : rings_) ApplyActorStyling(rg);
}

}  // namespace h5reader::app

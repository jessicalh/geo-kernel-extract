#include "QtOccupancyShellsOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/AtomSelection.h"

#include <QLoggingCategory>

#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <vector>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cOcc, "h5reader.overlay.occupancy")
}  // namespace

QtOccupancyShellsOverlay::QtOccupancyShellsOverlay(
    vtkSmartPointer<vtkRenderer>                  renderer,
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
    QObject* parent)
    : QObject(parent),
      renderer_(std::move(renderer)),
      renderWindow_(std::move(renderWindow))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtOccupancyShellsOverlay"));

    // One hue (violet), two alphas: inner 50% reads denser, outer 90% faint.
    // cfg_.massFractions default {0.5, 0.9} — keep shells_ in the same order so
    // r.shells[s] maps to shells_[s].
    shells_[0].fraction = 0.5;
    shells_[0].opacity  = 0.32;
    shells_[0].color    = {0.50, 0.35, 0.80};
    shells_[1].fraction = 0.9;
    shells_[1].opacity  = 0.14;
    shells_[1].color    = {0.50, 0.35, 0.80};
}

QtOccupancyShellsOverlay::~QtOccupancyShellsOverlay() {
    for (const auto& s : shells_)
        if (s.actor) renderer_->RemoveActor(s.actor);
}

void QtOccupancyShellsOverlay::Build(const model::QtProtein& protein,
                                     model::Conformation&    conformation) {
    ASSERT_THREAD(this);

    if (protein_ == &protein && conformation_ == &conformation && imageData_)
        return;  // idempotent on the same inputs

    for (const auto& s : shells_)
        if (s.actor) renderer_->RemoveActor(s.actor);

    protein_      = &protein;
    conformation_ = &conformation;
    hasShells_    = false;
    dirty_        = false;

    // Density grid — reconfigured (dims/spacing/origin + scalar values) on each
    // rebuild, since each focus atom gets its own grid. Point scalars, double
    // precision so the contour level matches the HDR level computed on the same
    // values exactly. x-fastest order == vtkImageData == GridSpec::index, so the
    // fill is a single std::copy.
    scalars_ = vtkSmartPointer<vtkDoubleArray>::New();
    scalars_->SetName("occupancy");
    scalars_->SetNumberOfComponents(1);

    imageData_ = vtkSmartPointer<vtkImageData>::New();
    imageData_->GetPointData()->SetScalars(scalars_);

    producer_ = vtkSmartPointer<vtkTrivialProducer>::New();
    producer_->SetOutput(imageData_);

    for (auto& s : shells_) {
        s.contour = vtkSmartPointer<vtkContourFilter>::New();
        s.contour->SetInputConnection(producer_->GetOutputPort());
        s.contour->SetValue(0, 0.0);  // real level set per rebuild

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(s.contour->GetOutputPort());
        mapper->ScalarVisibilityOff();

        s.actor = vtkSmartPointer<vtkActor>::New();
        s.actor->SetMapper(mapper);
        renderer_->AddActor(s.actor);
    }

    applyActorStyling();
    hideShells();  // nothing focused yet

    qCInfo(cOcc).noquote()
        << "Built occupancy-shells overlay (idle until an atom is focused)";
}

void QtOccupancyShellsOverlay::applyActorStyling() {
    for (auto& s : shells_) {
        if (!s.actor) continue;
        auto* prop = s.actor->GetProperty();
        prop->SetColor(s.color[0], s.color[1], s.color[2]);
        prop->SetOpacity(s.opacity);
        prop->SetInterpolationToPhong();
        // Per qt6-cpp skill 3d-vtk.md §8: translucent surfaces route through the
        // depth-peeling pass deterministically with ForceTranslucent, and
        // backface culling stops a closed shell rendering its own far wall
        // through the near one.
        s.actor->SetForceTranslucent(true);
        prop->SetBackfaceCulling(true);
    }
}

void QtOccupancyShellsOverlay::setSelection(model::AtomSelection* selection) {
    selection_ = selection;
}

bool QtOccupancyShellsOverlay::worldBounds(double out[6]) const {
    if (!hasShells_) return false;
    for (int i = 0; i < 6; ++i) out[i] = worldBounds_[i];
    return true;
}

void QtOccupancyShellsOverlay::onFocusChanged(std::size_t /*atomIdx*/) {
    ASSERT_THREAD(this);
    if (!visible_) { dirty_ = true; return; }
    rebuild();
}

void QtOccupancyShellsOverlay::onSelectionCleared() {
    ASSERT_THREAD(this);
    hideShells();
}

void QtOccupancyShellsOverlay::onTransformChanged() {
    ASSERT_THREAD(this);
    if (!visible_) { dirty_ = true; return; }
    rebuild();  // every aligned position moved — the aggregate is stale
}

void QtOccupancyShellsOverlay::setVisible(bool on) {
    ASSERT_THREAD(this);
    visible_ = on;
    if (on) {
        dirty_ = false;
        rebuild();  // populate for the current focus (covers dirty-while-hidden)
    } else {
        hideShells();
    }
}

void QtOccupancyShellsOverlay::hideShells() {
    for (auto& s : shells_)
        if (s.actor) s.actor->SetVisibility(0);
    hasShells_ = false;
}

void QtOccupancyShellsOverlay::rebuild() {
    ASSERT_THREAD(this);
    if (!visible_ || !protein_ || !conformation_ || !imageData_ ||
        !selection_ || !selection_->hasFocus()) {
        hideShells();
        return;
    }
    const std::size_t idx = selection_->focus();
    if (idx >= protein_->atomCount()) {
        hideShells();
        return;
    }

    const std::size_t T = conformation_->frameCount();
    std::vector<model::Vec3> pos;
    pos.reserve(T);
    for (std::size_t t = 0; t < T; ++t)
        pos.push_back(conformation_->atomPosition(t, idx));

    const auto r = math::computeOccupancy(pos, cfg_);
    if (!r.computed) {
        qCInfo(cOcc).noquote()
            << "atom" << idx << "— no shells:" << QString::fromStdString(r.note);
        hideShells();
        return;
    }

    // Reconfigure the grid and fill the density (x-fastest == GridSpec::index).
    const auto& g = r.field.grid;
    imageData_->SetDimensions(g.dims[0], g.dims[1], g.dims[2]);
    imageData_->SetSpacing(g.spacing, g.spacing, g.spacing);
    imageData_->SetOrigin(g.origin.x(), g.origin.y(), g.origin.z());
    const vtkIdType n = static_cast<vtkIdType>(r.field.values.size());
    scalars_->SetNumberOfTuples(n);
    std::copy(r.field.values.begin(), r.field.values.end(), scalars_->GetPointer(0));
    imageData_->GetPointData()->SetScalars(scalars_);
    scalars_->Modified();
    imageData_->Modified();

    // Set each shell's iso-level from the HDR result and force the contour to
    // execute now (GUI thread) so a degenerate/empty surface is visible in the
    // log before render, not a silent blank.
    for (std::size_t s = 0; s < shells_.size(); ++s) {
        if (s < r.shells.size() && r.shells[s].valid) {
            shells_[s].contour->SetValue(0, r.shells[s].isoValue);
            shells_[s].contour->Update();
            shells_[s].actor->SetVisibility(1);
        } else {
            shells_[s].actor->SetVisibility(0);
        }
    }

    worldBounds_[0] = g.origin.x();
    worldBounds_[1] = g.origin.x() + g.spacing * (g.dims[0] - 1);
    worldBounds_[2] = g.origin.y();
    worldBounds_[3] = g.origin.y() + g.spacing * (g.dims[1] - 1);
    worldBounds_[4] = g.origin.z();
    worldBounds_[5] = g.origin.z() + g.spacing * (g.dims[2] - 1);
    hasShells_ = true;

    qCInfo(cOcc).noquote()
        << "shells | atom" << idx
        << "| frames=" << T
        << "| RMSF=" << QString::number(r.stats.rmsf, 'f', 2) << "A"
        << "| n_eff=" << QString::number(r.stats.nEff, 'f', 1)
        << "| 50% iso=" << r.shells[0].isoValue
        << "(mass " << QString::number(r.shells[0].includedMass, 'f', 3) << ")"
        << "| 90% iso=" << r.shells[1].isoValue
        << "(mass " << QString::number(r.shells[1].includedMass, 'f', 3) << ")"
        << (r.note.empty() ? QString() : QString(" | ") + QString::fromStdString(r.note));
}

}  // namespace h5reader::app

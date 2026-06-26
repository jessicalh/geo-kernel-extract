#include "QtAtomTrajectoryOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/AtomSelection.h"

#include <QElapsedTimer>
#include <QLoggingCategory>

#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cTrajectory, "h5reader.overlay.trajectory")

std::array<double, 3> shellColor(std::optional<double> delta, double scale) {
    const std::array<double, 3> base = {0.50, 0.35, 0.80};
    if (!delta || !std::isfinite(*delta) || !std::isfinite(scale) || scale <= 1e-12)
        return base;

    const std::array<double, 3> warm = {0.92, 0.58, 0.18};
    const std::array<double, 3> cool = {0.20, 0.66, 0.86};
    const auto& end = *delta >= 0.0 ? warm : cool;
    const double t = 0.45 * std::clamp(std::abs(*delta) / scale, 0.0, 1.0);

    std::array<double, 3> out{};
    for (int i = 0; i < 3; ++i)
        out[static_cast<std::size_t>(i)] = base[i] + (end[i] - base[i]) * t;
    return out;
}
}  // namespace

QtAtomTrajectoryOverlay::QtAtomTrajectoryOverlay(vtkSmartPointer<vtkRenderer> renderer,
                                                 QObject* parent)
    : QObject(parent),
      renderer_(std::move(renderer)) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtAtomTrajectoryOverlay"));

    cfg_.voxelTarget = 0.18;
    cfg_.maxDim = 72;
    cfg_.marginSigma = 3.0;
    cfg_.floorVoxelFactor = 1.2;
    cfg_.rigidRmsfFactor = 0.75;
    cfg_.minFrames = 8;
    cfg_.massFractions = {0.25, 0.35};

    scalars_ = vtkSmartPointer<vtkDoubleArray>::New();
    scalars_->SetName("trajectory_density");
    scalars_->SetNumberOfComponents(1);

    imageData_ = vtkSmartPointer<vtkImageData>::New();
    imageData_->GetPointData()->SetScalars(scalars_);

    producer_ = vtkSmartPointer<vtkTrivialProducer>::New();
    producer_->SetOutput(imageData_);

    shells_[0].fraction = 0.25;
    shells_[0].opacity = 0.26;
    shells_[1].fraction = 0.35;
    shells_[1].opacity = 0.035;

    for (auto& shell : shells_) {
        shell.contour = vtkSmartPointer<vtkContourFilter>::New();
        shell.contour->SetInputConnection(producer_->GetOutputPort());
        shell.contour->SetValue(0, 0.0);

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(shell.contour->GetOutputPort());
        mapper->ScalarVisibilityOff();

        shell.actor = vtkSmartPointer<vtkActor>::New();
        shell.actor->SetMapper(mapper);
        shell.actor->SetVisibility(0);
        renderer_->AddActor(shell.actor);
    }
    applyShellStyling(std::nullopt, 0.0);
}

QtAtomTrajectoryOverlay::~QtAtomTrajectoryOverlay() {
    if (!renderer_)
        return;
    for (const auto& shell : shells_) {
        if (shell.actor)
            renderer_->RemoveActor(shell.actor);
    }
}

void QtAtomTrajectoryOverlay::Build(const model::QtProtein& protein,
                                    model::Conformation& conformation) {
    ASSERT_THREAD(this);
    protein_ = &protein;
    conformation_ = &conformation;
    dirty_ = true;
    hideEnvelope();
}

void QtAtomTrajectoryOverlay::setSelection(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    selection_ = selection;
    dirty_ = true;
}

void QtAtomTrajectoryOverlay::setDftStore(model::DftShieldingStore* store) {
    ASSERT_THREAD(this);
    dftStore_ = store;
    dirty_ = true;
}

void QtAtomTrajectoryOverlay::setFrame(int frame) {
    ASSERT_THREAD(this);
    currentFrame_ = std::max(0, frame);
    if (!visible_)
        return;
    if (dirty_ || !hasEnvelope_)
        rebuild();
}

void QtAtomTrajectoryOverlay::onFocusChanged(std::size_t /*atomIdx*/) {
    ASSERT_THREAD(this);
    dirty_ = true;
    if (visible_)
        rebuild();
}

void QtAtomTrajectoryOverlay::onSelectionCleared() {
    ASSERT_THREAD(this);
    dirty_ = true;
    hideEnvelope();
}

void QtAtomTrajectoryOverlay::onTransformChanged() {
    ASSERT_THREAD(this);
    dirty_ = true;
    if (visible_)
        rebuild();
}

void QtAtomTrajectoryOverlay::setVisible(bool on) {
    ASSERT_THREAD(this);
    visible_ = on;
    if (on) {
        rebuild();
    } else {
        hideEnvelope();
    }
}

void QtAtomTrajectoryOverlay::hideEnvelope() {
    for (auto& shell : shells_) {
        if (shell.actor)
            shell.actor->SetVisibility(0);
    }
    hasEnvelope_ = false;
    windowStart_ = -1;
    windowEnd_ = -1;
}

void QtAtomTrajectoryOverlay::applyShellStyling(std::optional<double> trendDelta,
                                                double trendScale) {
    const auto color = shellColor(trendDelta, trendScale);
    for (auto& shell : shells_) {
        if (!shell.actor)
            continue;
        auto* prop = shell.actor->GetProperty();
        prop->SetColor(color[0], color[1], color[2]);
        prop->SetOpacity(shell.opacity);
        prop->SetInterpolationToPhong();
        prop->SetSpecular(0.12);
        prop->SetAmbient(0.22);
        prop->SetDiffuse(0.78);
        shell.actor->SetForceTranslucent(true);
        prop->SetBackfaceCulling(true);
    }
}

void QtAtomTrajectoryOverlay::clearScalarCacheForAtom(std::size_t atom) {
    if (cachedAtom_ && *cachedAtom_ == atom)
        return;
    cachedAtom_ = atom;
    orcaT0ByOriginal_.clear();
}

std::optional<double> QtAtomTrajectoryOverlay::sampleOrcaT0(std::size_t frame,
                                                            std::size_t atom) {
    if (!conformation_ || !dftStore_)
        return std::nullopt;
    clearScalarCacheForAtom(atom);
    const std::size_t original = conformation_->originalFrameIndex(frame);
    const auto cached = orcaT0ByOriginal_.find(original);
    if (cached != orcaT0ByOriginal_.end())
        return cached->second;
    if (!dftStore_->hasJob(original)) {
        orcaT0ByOriginal_.emplace(original, std::nullopt);
        return std::nullopt;
    }
    dftStore_->requestFrame(original);
    std::optional<double> value =
        dftStore_->sample(original, atom,
                          model::DftPart::Total,
                          model::DftScalar::IsotropicT0);
    if (value && !std::isfinite(*value))
        value.reset();
    orcaT0ByOriginal_.emplace(original, value);
    return value;
}

void QtAtomTrajectoryOverlay::rebuild() {
    ASSERT_THREAD(this);
    dirty_ = false;
    if (!visible_ || !protein_ || !conformation_ || !selection_ ||
        !selection_->hasFocus()) {
        hideEnvelope();
        return;
    }

    const std::size_t atom = selection_->focus();
    clearScalarCacheForAtom(atom);
    if (atom >= protein_->atomCount()) {
        hideEnvelope();
        return;
    }

    const int frameCount = static_cast<int>(conformation_->frameCount());
    if (frameCount < static_cast<int>(cfg_.minFrames)) {
        hideEnvelope();
        return;
    }

    QElapsedTimer timer;
    timer.start();

    const int center = std::clamp(currentFrame_, 0, frameCount - 1);
    const int start = 0;
    const int end = frameCount - 1;
    const int nFrames = end - start + 1;
    if (nFrames < static_cast<int>(cfg_.minFrames)) {
        hideEnvelope();
        return;
    }

    emit rebuildStarted(nFrames);

    std::vector<model::Vec3> positions;
    positions.reserve(static_cast<std::size_t>(nFrames));
    for (int frame = start; frame <= end; ++frame) {
        positions.push_back(
            conformation_->atomPosition(static_cast<std::size_t>(frame), atom));
    }

    const std::optional<double> currentSigma =
        sampleOrcaT0(static_cast<std::size_t>(center), atom);
    int dftSamples = currentSigma ? 1 : 0;

    const auto r = math::computeOccupancy(positions, cfg_);
    if (!r.computed) {
        qCInfo(cTrajectory).noquote()
            << "envelope | atom=" << atom
            << "| frames=" << start << ".." << end
            << "| skipped=" << QString::fromStdString(r.note);
        hideEnvelope();
        emit rebuildFinished(nFrames, dftSamples, static_cast<int>(timer.elapsed()));
        return;
    }

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

    applyShellStyling(std::nullopt, 0.0);
    for (std::size_t s = 0; s < shells_.size(); ++s) {
        if (s < r.shells.size() && r.shells[s].valid) {
            shells_[s].contour->SetValue(0, r.shells[s].isoValue);
            shells_[s].contour->Update();
            shells_[s].actor->SetVisibility(1);
        } else if (shells_[s].actor) {
            shells_[s].actor->SetVisibility(0);
        }
    }

    hasEnvelope_ = true;
    windowStart_ = start;
    windowEnd_ = end;

    if (dftStore_ && conformation_) {
        const std::size_t currentOriginal =
            conformation_->originalFrameIndex(static_cast<std::size_t>(center));
        if (dftStore_->hasJob(currentOriginal))
            dftStore_->requestFrame(currentOriginal);
    }

    qCInfo(cTrajectory).noquote()
        << "envelope | atom=" << atom
        << "| frames=" << start << ".." << end
        << "| dft_samples=" << dftSamples << "/" << nFrames
        << "| cache_entries=" << static_cast<int>(orcaT0ByOriginal_.size())
        << "| RMSF=" << QString::number(r.stats.rmsf, 'f', 3) << "A"
        << "| 25% iso=" << r.shells[0].isoValue
        << "| 35% iso=" << r.shells[1].isoValue
        << "| current_t0=" << (currentSigma ? QString::number(*currentSigma, 'f', 3)
                                           : QStringLiteral("n/a"))
        << "| load_ms=" << timer.elapsed();
    emit rebuildFinished(nFrames, dftSamples, static_cast<int>(timer.elapsed()));
}

}  // namespace h5reader::app

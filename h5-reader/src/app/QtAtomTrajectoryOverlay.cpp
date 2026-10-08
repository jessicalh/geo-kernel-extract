#include "QtAtomTrajectoryOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/AtomSelection.h"

#include <QElapsedTimer>
#include <QLoggingCategory>
#include <QThread>

#include <vtkCallbackCommand.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
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

// No model or live rendering objects cross this boundary. Results are read
// only after QThread::finished reaches the overlay's GUI-thread slot.
class QtAtomTrajectoryOverlay::EnvelopeWorker final : public QThread {
public:
    EnvelopeWorker(std::vector<model::Vec3> positions,
                   math::OccupancyConfig config, std::uint64_t generation,
                   std::size_t atom, QElapsedTimer timer)
        : generation(generation), atom(atom),
          frameCount(static_cast<int>(positions.size())), timer(timer),
          positions_(std::move(positions)), config_(std::move(config)) {
        CENSUS_REGISTER(this);
        setObjectName(QStringLiteral("TrajectoryEnvelopeWorker"));
    }

    const std::uint64_t generation;
    const std::size_t atom;
    const int frameCount;
    const QElapsedTimer timer;
    math::OccupancyResult result;
    std::array<vtkSmartPointer<vtkPolyData>, 2> geometry;
    QString error;

protected:
    void run() override {
        try {
            if (isInterruptionRequested())
                return;
            result = math::computeOccupancy(positions_, config_);
            if (!result.computed || isInterruptionRequested())
                return;

            auto scalars = vtkSmartPointer<vtkDoubleArray>::New();
            scalars->SetName("trajectory_density");
            scalars->SetNumberOfComponents(1);
            scalars->SetNumberOfTuples(static_cast<vtkIdType>(result.field.values.size()));
            std::copy(result.field.values.begin(), result.field.values.end(),
                      scalars->GetPointer(0));

            const auto& g = result.field.grid;
            auto image = vtkSmartPointer<vtkImageData>::New();
            image->SetDimensions(g.dims[0], g.dims[1], g.dims[2]);
            image->SetSpacing(g.spacing, g.spacing, g.spacing);
            image->SetOrigin(g.origin.x(), g.origin.y(), g.origin.z());
            image->GetPointData()->SetScalars(scalars);

            auto producer = vtkSmartPointer<vtkTrivialProducer>::New();
            producer->SetOutput(image);
            auto errors = vtkSmartPointer<vtkCallbackCommand>::New();
            errors->SetClientData(&error);
            errors->SetCallback([](vtkObject*, unsigned long, void* data, void* message) {
                *static_cast<QString*>(data) = message
                    ? QString::fromUtf8(static_cast<const char*>(message))
                    : QStringLiteral("VTK contour failed");
            });

            for (std::size_t s = 0; s < geometry.size(); ++s) {
                if (isInterruptionRequested())
                    return;
                if (s >= result.shells.size() || !result.shells[s].valid)
                    continue;
                auto contour = vtkSmartPointer<vtkContourFilter>::New();
                contour->AddObserver(vtkCommand::ErrorEvent, errors);
                contour->SetInputConnection(producer->GetOutputPort());
                contour->SetValue(0, result.shells[s].isoValue);
                contour->Update();
                if (!error.isEmpty() || isInterruptionRequested())
                    return;
                // Detach output from the worker's pipeline before GUI handoff.
                geometry[s] = vtkSmartPointer<vtkPolyData>::New();
                geometry[s]->DeepCopy(contour->GetOutput());
            }
        } catch (const std::exception& e) {
            error = QString::fromUtf8(e.what());
        } catch (...) {
            error = QStringLiteral("Unknown occupancy-envelope failure");
        }
    }

private:
    const std::vector<model::Vec3> positions_;
    const math::OccupancyConfig config_;
};

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

    shells_[0].fraction = 0.25;
    shells_[0].opacity = 0.26;
    shells_[1].fraction = 0.35;
    shells_[1].opacity = 0.035;

    for (auto& shell : shells_) {
        shell.producer = vtkSmartPointer<vtkTrivialProducer>::New();
        shell.producer->SetOutput(vtkSmartPointer<vtkPolyData>::New());

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(shell.producer->GetOutputPort());
        mapper->ScalarVisibilityOff();

        shell.actor = vtkSmartPointer<vtkActor>::New();
        shell.actor->SetMapper(mapper);
        shell.actor->SetVisibility(0);
        renderer_->AddActor(shell.actor);
    }
    applyShellStyling(std::nullopt, 0.0);
}

QtAtomTrajectoryOverlay::~QtAtomTrajectoryOverlay() {
    ASSERT_THREAD(this);
    if (worker_) {
        worker_->requestInterruption();
        // Fallback only: normal close/replacement waits for becameIdle().
        if (!worker_->isFinished())
            qCWarning(cTrajectory) << "Destroying busy trajectory overlay; waiting for worker";
        worker_->wait();
        worker_.reset();
    }
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
    cachedAtom_.reset();
    orcaT0ByOriginal_.clear();
    invalidateEnvelope();
}

void QtAtomTrajectoryOverlay::setSelection(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    selection_ = selection;
    invalidateEnvelope();
    if (visible_)
        rebuild();
}

void QtAtomTrajectoryOverlay::setDftStore(model::DftShieldingStore* store) {
    ASSERT_THREAD(this);
    if (dftStore_ == store)
        return;
    if (dftStore_)
        disconnect(dftStore_, nullptr, this, nullptr);
    dftStore_ = store;
    cachedAtom_.reset();
    orcaT0ByOriginal_.clear();
    if (dftStore_) {
        Q_ASSERT(dftStore_->thread() == thread());
        connect(dftStore_, &model::DftShieldingStore::frameReady,
                this, &QtAtomTrajectoryOverlay::onDftFrameReady);
    }
}

bool QtAtomTrajectoryOverlay::isBusy() const {
    ASSERT_THREAD(this);
    return worker_ != nullptr;
}

bool QtAtomTrajectoryOverlay::captureReady(QString* error) const {
    ASSERT_THREAD(this);
    if (!visible_ || !selection_ || !selection_->hasFocus())
        return true;
    if (worker_ || dirty_)
        return false;
    if (!error_.isEmpty()) {
        *error = QStringLiteral("Trajectory envelope: %1").arg(error_);
        return false;
    }
    return true;
}

void QtAtomTrajectoryOverlay::cancelPending() {
    ASSERT_THREAD(this);
    ++generation_;
    dirty_ = false;
    cancelling_ = isBusy();
    if (worker_)
        worker_->requestInterruption();
    hideEnvelope();
}

void QtAtomTrajectoryOverlay::setFrame(int frame) {
    ASSERT_THREAD(this);
    currentFrame_ = std::max(0, frame);
    if (!visible_ || cancelling_)
        return;
    // The envelope covers all frames; playback only changes the ORCA sample.
    if (conformation_ && selection_ && selection_->hasFocus() && conformation_->frameCount()) {
        const auto center = std::min(static_cast<std::size_t>(currentFrame_),
                                     conformation_->frameCount() - 1);
        sampleOrcaT0(center, selection_->focus());
    }
    rebuild();
}

void QtAtomTrajectoryOverlay::onFocusChanged(std::size_t /*atomIdx*/) {
    ASSERT_THREAD(this);
    invalidateEnvelope();
    if (visible_)
        rebuild();
}

void QtAtomTrajectoryOverlay::onSelectionCleared() {
    ASSERT_THREAD(this);
    invalidateEnvelope();
}

void QtAtomTrajectoryOverlay::onTransformChanged() {
    ASSERT_THREAD(this);
    invalidateEnvelope();
    if (visible_)
        rebuild();
}

void QtAtomTrajectoryOverlay::setVisible(bool on) {
    ASSERT_THREAD(this);
    visible_ = on;
    invalidateEnvelope();
    if (on)
        rebuild();
}

void QtAtomTrajectoryOverlay::invalidateEnvelope() {
    ASSERT_THREAD(this);
    ++generation_;
    error_.clear();
    dirty_ = !cancelling_;
    if (worker_)
        worker_->requestInterruption();
    hideEnvelope();
}

void QtAtomTrajectoryOverlay::hideEnvelope() {
    ASSERT_THREAD(this);
    for (auto& shell : shells_) {
        if (shell.actor)
            shell.actor->SetVisibility(0);
    }
}

void QtAtomTrajectoryOverlay::applyShellStyling(std::optional<double> trendDelta,
                                                double trendScale) {
    ASSERT_THREAD(this);
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
    ASSERT_THREAD(this);
    if (cachedAtom_ && *cachedAtom_ == atom)
        return;
    cachedAtom_ = atom;
    orcaT0ByOriginal_.clear();
}

std::optional<double> QtAtomTrajectoryOverlay::sampleOrcaT0(std::size_t frame,
                                                            std::size_t atom) {
    ASSERT_THREAD(this);
    if (!conformation_ || !dftStore_)
        return std::nullopt;
    clearScalarCacheForAtom(atom);
    const std::size_t original = conformation_->originalFrameIndex(frame);
    const auto cached = orcaT0ByOriginal_.find(original);
    if (cached != orcaT0ByOriginal_.end())
        return cached->second;
    if (!dftStore_->hasJob(original) || dftStore_->hasFailedFrame(original)) {
        orcaT0ByOriginal_.emplace(original, std::nullopt);
        return std::nullopt;
    }
    if (!dftStore_->frame(original)) {
        dftStore_->requestFrameAsync(original);
        return std::nullopt;  // Pending is not a cached missing sample.
    }
    std::optional<double> value =
        dftStore_->sample(original, atom,
                          model::DftPart::Total,
                          model::DftScalar::IsotropicT0);
    if (value && !std::isfinite(*value))
        value.reset();
    orcaT0ByOriginal_.emplace(original, value);
    return value;
}

void QtAtomTrajectoryOverlay::onDftFrameReady(std::size_t original) {
    ASSERT_THREAD(this);
    if (!dftStore_ || sender() != dftStore_ || !visible_ || cancelling_ || !conformation_ ||
        !selection_ || !selection_->hasFocus() || !conformation_->frameCount())
        return;
    const auto center = std::min(static_cast<std::size_t>(currentFrame_),
                                 conformation_->frameCount() - 1);
    if (conformation_->originalFrameIndex(center) != original)
        return;
    // Another observer may have replaced the store's one resident frame.
    if (!dftStore_->frame(original) && !dftStore_->hasFailedFrame(original))
        return;
    const auto atom = selection_->focus();
    const auto value = sampleOrcaT0(center, atom);
    qCInfo(cTrajectory).noquote()
        << "ORCA sample | atom=" << atom << "| original_frame=" << original
        << "| current_t0=" << (value ? QString::number(*value, 'f', 3)
                                     : QStringLiteral("n/a"));
}

void QtAtomTrajectoryOverlay::rebuild() {
    ASSERT_THREAD(this);
    if (cancelling_ || worker_ || !dirty_)
        return;
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
        error_ = QStringLiteral("not enough frames to calculate the envelope");
        hideEnvelope();
        return;
    }

    QElapsedTimer timer;
    timer.start();

    try {
        std::vector<model::Vec3> positions;
        positions.reserve(static_cast<std::size_t>(frameCount));
        for (int frame = 0; frame < frameCount; ++frame) {
            positions.push_back(
                conformation_->atomPosition(static_cast<std::size_t>(frame), atom));
        }
        worker_ = std::make_unique<EnvelopeWorker>(
            std::move(positions), cfg_, generation_, atom, timer);
    } catch (const std::exception& e) {
        error_ = QString::fromUtf8(e.what());
        qCWarning(cTrajectory) << "Envelope snapshot failed | atom=" << atom << e.what();
        hideEnvelope();
        return;
    }

    connect(worker_.get(), &QThread::finished,
            this, &QtAtomTrajectoryOverlay::finishRebuild, Qt::QueuedConnection);
    worker_->start();
    qCDebug(cTrajectory) << "Envelope worker started | atom=" << atom
                        << "| generation=" << generation_ << "| frames=" << frameCount;
    const auto center = static_cast<std::size_t>(std::clamp(currentFrame_, 0, frameCount - 1));
    sampleOrcaT0(center, atom);
    emit rebuildStarted(frameCount);
}

void QtAtomTrajectoryOverlay::finishRebuild() {
    ASSERT_THREAD(this);
    const auto& work = *worker_;
    const auto& r = work.result;
    const bool current = !cancelling_ && work.generation == generation_ && visible_ &&
        conformation_ && selection_ && selection_->hasFocus() &&
        selection_->focus() == work.atom &&
        conformation_->frameCount() == static_cast<std::size_t>(work.frameCount);
    const int loadMs = static_cast<int>(work.timer.elapsed());
    int dftSamples = 0;
    if (!work.error.isEmpty()) {
        qCWarning(cTrajectory).noquote()
            << "Envelope calculation failed | atom=" << work.atom
            << "| generation=" << work.generation << "| error=" << work.error;
    }
    if (current)
        error_ = work.error;
    if (current && work.error.isEmpty()) {
        const auto center = static_cast<std::size_t>(
            std::clamp(currentFrame_, 0, work.frameCount - 1));
        const auto currentSigma = sampleOrcaT0(center, work.atom);
        dftSamples = currentSigma ? 1 : 0;
        if (r.computed) {
            applyShellStyling(std::nullopt, 0.0);
            for (std::size_t s = 0; s < shells_.size(); ++s) {
                if (work.geometry[s]) {
                    shells_[s].producer->SetOutput(work.geometry[s]);
                    shells_[s].actor->SetVisibility(1);
                } else {
                    shells_[s].actor->SetVisibility(0);
                }
            }
            qCInfo(cTrajectory).noquote()
                << "envelope | atom=" << work.atom
                << "| frames=0 .." << work.frameCount - 1
                << "| dft_samples=" << dftSamples << "/" << work.frameCount
                << "| cache_entries=" << static_cast<int>(orcaT0ByOriginal_.size())
                << "| RMSF=" << QString::number(r.stats.rmsf, 'f', 3) << "A"
                << "| 25% iso=" << r.shells[0].isoValue
                << "| 35% iso=" << r.shells[1].isoValue
                << "| current_t0=" << (currentSigma ? QString::number(*currentSigma, 'f', 3)
                                                     : QStringLiteral("n/a"))
                << "| load_ms=" << loadMs;
        } else {
            error_ = QString::fromStdString(r.note);
            qCInfo(cTrajectory).noquote()
                << "envelope | atom=" << work.atom
                << "| frames=0 .." << work.frameCount - 1
                << "| skipped=" << QString::fromStdString(r.note);
            hideEnvelope();
        }
    } else if (!current) {
        qCDebug(cTrajectory) << "Envelope result discarded | generation=" << work.generation;
    } else {
        hideEnvelope();
    }

    // Balance every started notification, including discarded/failed work.
    // Keep busy during callbacks so lifecycle consumers wait for becameIdle.
    const QPointer<QtAtomTrajectoryOverlay> guard(this);
    emit rebuildFinished(work.frameCount, dftSamples, loadMs);
    if (!guard)
        return;
    worker_.reset();
    cancelling_ = false;
    rebuild();
    if (guard && !isBusy())
        emit becameIdle();
}

}  // namespace h5reader::app

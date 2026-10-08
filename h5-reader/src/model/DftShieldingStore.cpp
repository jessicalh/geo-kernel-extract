#include "DftShieldingStore.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/DftShieldingLoader.h"

#include <QLoggingCategory>
#include <QPointer>
#include <QThread>

#include <exception>
#include <utility>

namespace h5reader::model {

namespace {
Q_LOGGING_CATEGORY(cDft, "h5reader.dft")

std::shared_ptr<const DftShieldingFrame> loadAndValidate(
    const QString& metaJson, const QtProtein* protein) {
    try {
        return h5reader::io::DftShieldingLoader::LoadAndValidate(metaJson, protein);
    } catch (const std::exception& error) {
        qCWarning(cDft).noquote() << "DFT load failed | meta=" << metaJson
                                << "|" << error.what();
    } catch (...) {
        qCWarning(cDft).noquote() << "DFT load failed | meta=" << metaJson
                                << "| unknown exception";
    }
    return nullptr;
}
}  // namespace

DftShieldingStore::DftShieldingStore(const QtProtein* protein,
                                     const std::vector<h5reader::io::DftFrame>& frames,
                                     QObject* parent)
    : QObject(parent), protein_(protein) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("DftShieldingStore"));

    metaByOriginal_.reserve(frames.size());
    for (const auto& fr : frames) {
        if (fr.frame_index < 0 || fr.meta_json_abspath.isEmpty()) continue;
        metaByOriginal_.emplace(static_cast<std::size_t>(fr.frame_index),
                                fr.meta_json_abspath);
    }
    qCInfo(cDft).noquote()
        << "DFT store initialised from .LGS frames |"
        << "frames=" << metaByOriginal_.size();
}

DftShieldingStore::~DftShieldingStore() {
    ASSERT_THREAD(this);
    // Fallback for owners that destroy a busy store. Normal close/replacement
    // uses cancelPending()/becameIdle() while retaining the protein.
    if (asyncThread_)
        asyncThread_->wait();
}

bool DftShieldingStore::hasJob(std::size_t originalIndex) const {
    return metaByOriginal_.find(originalIndex) != metaByOriginal_.end();
}

bool DftShieldingStore::hasFailedFrame(std::size_t originalIndex) const {
    return hasJob(originalIndex) && resolvedAbsent_.count(originalIndex) != 0;
}

const DftShieldingFrame* DftShieldingStore::frame(std::size_t originalIndex) const {
    if (!residentOriginal_ || *residentOriginal_ != originalIndex)
        return nullptr;
    return residentFrame_.get();
}

void DftShieldingStore::requestFrame(std::size_t originalIndex) {
    ASSERT_THREAD(this);
    if (cancelling_)
        return;
    // Idempotent: a resident frame or known-absent frame just re-announces.
    // The parsed frame is a temporary source object. Strips keep the durable
    // sampled display values in their ChannelBuffers.
    if (residentOriginal_ && *residentOriginal_ == originalIndex) {
        emit frameReady(originalIndex);
        return;
    }

    const auto it = metaByOriginal_.find(originalIndex);
    publishFrame(originalIndex,
                 it != metaByOriginal_.end() && !resolvedAbsent_.count(originalIndex)
                     ? loadAndValidate(it->second, protein_) : nullptr);
}

void DftShieldingStore::requestFrameAsync(std::size_t originalIndex) {
    ASSERT_THREAD(this);
    if (cancelling_ || asyncActiveOriginal_ == originalIndex
        || !asyncPendingSet_.insert(originalIndex).second)
        return;
    asyncPending_.push_back(originalIndex);
    if (!asyncActiveOriginal_)
        startNextAsyncFrame();
}

bool DftShieldingStore::isBusy() const {
    ASSERT_THREAD(this);
    return asyncActiveOriginal_.has_value();
}

void DftShieldingStore::cancelPending() {
    ASSERT_THREAD(this);
    asyncPending_.clear();
    asyncPendingSet_.clear();
    cancelling_ = isBusy();
}

void DftShieldingStore::publishFrame(
    std::size_t originalIndex,
    std::shared_ptr<const DftShieldingFrame> frame) {
    residentFrame_ = std::move(frame);
    if (residentFrame_)
        residentOriginal_ = originalIndex;
    else {
        residentOriginal_.reset();
        resolvedAbsent_.insert(originalIndex);
    }
    emit frameReady(originalIndex);
}

void DftShieldingStore::startNextAsyncFrame() {
    ASSERT_THREAD(this);
    while (!asyncPending_.empty()) {
        const std::size_t originalIndex = asyncPending_.front();
        asyncPending_.pop_front();
        asyncPendingSet_.erase(originalIndex);
        asyncActiveOriginal_ = originalIndex;

        const QPointer<DftShieldingStore> self(this);
        const auto it = metaByOriginal_.find(originalIndex);
        if (residentOriginal_ == originalIndex) {
            emit frameReady(originalIndex);
        } else if (resolvedAbsent_.count(originalIndex) || it == metaByOriginal_.end()) {
            publishFrame(originalIndex, nullptr);
        } else {
            const QString metaJson = it->second;
            const QtProtein* protein = protein_;
            const auto result = std::make_shared<std::shared_ptr<const DftShieldingFrame>>();
            asyncThread_ = QThread::create([metaJson, protein, result] {
                *result = loadAndValidate(metaJson, protein);
            });
            asyncThread_->setParent(this);
            connect(asyncThread_, &QThread::finished, this,
                    [this, originalIndex, result] {
                        finishAsyncFrame(originalIndex, std::move(*result));
                    }, Qt::QueuedConnection);
            asyncThread_->start();
            return;
        }
        // Direct frameReady handlers may enqueue, cancel, or destroy the store.
        if (!self)
            return;
        asyncActiveOriginal_.reset();
    }
    cancelling_ = false;
    emit becameIdle();
}

void DftShieldingStore::finishAsyncFrame(
    std::size_t originalIndex,
    std::shared_ptr<const DftShieldingFrame> frame) {
    ASSERT_THREAD(this);
    asyncThread_->deleteLater();
    asyncThread_ = nullptr;
    const QPointer<DftShieldingStore> self(this);
    if (!cancelling_)
        publishFrame(originalIndex, std::move(frame));
    if (!self)
        return;
    asyncActiveOriginal_.reset();
    startNextAsyncFrame();
}

std::optional<double> DftShieldingStore::sample(std::size_t originalIndex, std::size_t atom,
                                                DftPart part, DftScalar scalar) const {
    const DftShieldingFrame* framePtr = frame(originalIndex);
    if (!framePtr)
        return std::nullopt;  // not resident, or resolved-absent
    const DftShieldingFrame& fr = *framePtr;
    if (atom >= fr.atoms.size())
        return std::nullopt;
    const DftAtomShielding& a    = fr.atoms[atom];
    const SphericalTensor&  tens = (part == DftPart::Total) ? a.total
                                   : (part == DftPart::Dia) ? a.dia
                                                            : a.para;
    return (scalar == DftScalar::IsotropicT0) ? tens.T0 : tens.T2Magnitude();
}

}  // namespace h5reader::model

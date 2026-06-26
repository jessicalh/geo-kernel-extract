#include "DftShieldingStore.h"

#include "../diagnostics/ErrorBus.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/CalcsetManifest.h"
#include "../io/DftShieldingLoader.h"
#include "QtProtein.h"

#include <QLoggingCategory>
#include <QMetaObject>

#include <utility>

namespace h5reader::model {

namespace {
Q_LOGGING_CATEGORY(cDft, "h5reader.dft")
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
    for (std::thread& thread : asyncThreads_) {
        if (thread.joinable())
            thread.join();
    }
}

bool DftShieldingStore::hasJob(std::size_t originalIndex) const {
    return metaByOriginal_.find(originalIndex) != metaByOriginal_.end();
}

const DftShieldingFrame* DftShieldingStore::frame(std::size_t originalIndex) const {
    if (!residentOriginal_ || *residentOriginal_ != originalIndex)
        return nullptr;
    return residentFrame_.get();
}

void DftShieldingStore::requestFrame(std::size_t originalIndex) {
    ASSERT_THREAD(this);
    // Idempotent: a resident frame or known-absent frame just re-announces.
    // The parsed frame is a temporary source object. Strips keep the durable
    // sampled display values in their ChannelBuffers.
    if (residentOriginal_ && *residentOriginal_ == originalIndex) {
        emit frameReady(originalIndex);
        return;
    }

    residentOriginal_.reset();
    residentFrame_.reset();

    if (resolvedAbsent_.find(originalIndex) != resolvedAbsent_.end()) {
        emit frameReady(originalIndex);
        return;
    }

    residentFrame_ = loadAndValidate(originalIndex);
    if (residentFrame_) {
        residentOriginal_ = originalIndex;
    } else {
        resolvedAbsent_.insert(originalIndex);
    }
    emit frameReady(originalIndex);
}

void DftShieldingStore::requestFrameAsync(std::size_t originalIndex) {
    ASSERT_THREAD(this);

    if (residentOriginal_ && *residentOriginal_ == originalIndex) {
        emit frameReady(originalIndex);
        return;
    }
    if (resolvedAbsent_.find(originalIndex) != resolvedAbsent_.end())
        return;
    if (!asyncInFlight_.empty()) {
        asyncPendingOriginal_ = originalIndex;
        return;
    }
    if (asyncInFlight_.find(originalIndex) != asyncInFlight_.end())
        return;

    const auto it = metaByOriginal_.find(originalIndex);
    if (it == metaByOriginal_.end()) {
        resolvedAbsent_.insert(originalIndex);
        return;
    }

    const QString metaJson = it->second;
    const QtProtein* protein = protein_;
    asyncInFlight_.insert(originalIndex);
    asyncThreads_.emplace_back([this, originalIndex, metaJson, protein]() {
        std::shared_ptr<const DftShieldingFrame> loaded =
            h5reader::io::DftShieldingLoader::LoadAndValidate(metaJson, protein);
        QMetaObject::invokeMethod(this,
                                  [this, originalIndex, loaded]() {
                                      finishAsyncFrame(originalIndex, loaded);
                                  },
                                  Qt::QueuedConnection);
    });
}

void DftShieldingStore::finishAsyncFrame(
    std::size_t originalIndex,
    std::shared_ptr<const DftShieldingFrame> frame) {
    ASSERT_THREAD(this);
    asyncInFlight_.erase(originalIndex);
    if (frame) {
        residentOriginal_ = originalIndex;
        residentFrame_ = std::move(frame);
    } else {
        if (residentOriginal_ && *residentOriginal_ == originalIndex) {
            residentOriginal_.reset();
            residentFrame_.reset();
        }
        resolvedAbsent_.insert(originalIndex);
    }
    const std::optional<std::size_t> pending = asyncPendingOriginal_;
    asyncPendingOriginal_.reset();
    if (pending && *pending != originalIndex)
        requestFrameAsync(*pending);
    emit frameReady(originalIndex);
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

std::shared_ptr<const DftShieldingFrame>
DftShieldingStore::loadAndValidate(std::size_t originalIndex) const {
    const auto it = metaByOriginal_.find(originalIndex);
    if (it == metaByOriginal_.end()) return nullptr;
    return h5reader::io::DftShieldingLoader::LoadAndValidate(it->second, protein_);
}

}  // namespace h5reader::model

#include "DftShieldingStore.h"

#include "../diagnostics/ErrorBus.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/DftShieldingLoader.h"
#include "QtProtein.h"

#include <QDir>
#include <QLoggingCategory>
#include <QStringList>

namespace h5reader::model {

namespace {
Q_LOGGING_CATEGORY(cDft, "h5reader.dft")

using Severity = h5reader::diagnostics::Severity;
void report(Severity sev, const QString& msg, const QString& ctx) {
    h5reader::diagnostics::ErrorBus::Report(sev, QStringLiteral("DftShieldingStore"), msg, ctx);
}

}  // namespace

DftShieldingStore::DftShieldingStore(const QtProtein* protein, const QString& jobsDir, QObject* parent)
    : QObject(parent), protein_(protein), jobsDir_(jobsDir) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("DftShieldingStore"));

    QDir dir(jobsDir_);
    if (!dir.exists()) {
        report(Severity::Warning, QStringLiteral("dft jobs dir does not exist"), jobsDir_);
        return;
    }

    // Build originalIndex -> dir from the documented job-dir naming convention
    // (..._fNNNNNN_t<ps>): pull the digits between the last "_f" and the last
    // "_t". This is reading a documented layout (like TrajectoryConformation
    // scanning frame_NNNNNN), not glob-guessing.
    const QStringList entries = dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
    for (const QString& name : entries) {
        const int t = name.lastIndexOf(QStringLiteral("_t"));
        if (t < 0) continue;
        const int f = name.lastIndexOf(QStringLiteral("_f"), t);
        if (f < 0 || t <= f + 2) continue;
        bool ok = false;
        const qulonglong orig = name.mid(f + 2, t - (f + 2)).toULongLong(&ok);
        if (!ok) continue;
        dirByOriginal_.emplace(static_cast<std::size_t>(orig), dir.absoluteFilePath(name));
    }
    qCInfo(cDft).noquote() << "scanned dft jobs |" << jobsDir_ << "| jobs=" << dirByOriginal_.size();
}

bool DftShieldingStore::hasJob(std::size_t originalIndex) const {
    return dirByOriginal_.find(originalIndex) != dirByOriginal_.end();
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

std::shared_ptr<const DftShieldingFrame> DftShieldingStore::loadAndValidate(std::size_t originalIndex) const {
    return h5reader::io::DftShieldingLoader::LoadAndValidate(originalIndex, jobsDir_, protein_);
}

}  // namespace h5reader::model

#include "DftShieldingStore.h"

#include "../diagnostics/ErrorBus.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/OrcaShieldingParser.h"
#include "QtProtein.h"

#include <QByteArray>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLoggingCategory>
#include <QStringList>

#include <cmath>
#include <sstream>
#include <string>

namespace h5reader::model {

namespace {
Q_LOGGING_CATEGORY(cDft, "h5reader.dft")

using Severity = h5reader::diagnostics::Severity;
void report(Severity sev, const QString& msg, const QString& ctx) {
    h5reader::diagnostics::ErrorBus::Report(sev, QStringLiteral("DftShieldingStore"), msg, ctx);
}

// ppm tolerance for the ORCA identity total == dia + para. ORCA prints the
// tensors to ~3-4 decimals, so the reconstructed sum agrees to well under this.
constexpr double kIdentityTolPpm = 0.1;
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
    const auto it = cache_.find(originalIndex);
    return (it != cache_.end()) ? it->second.get() : nullptr;  // null value => resolved-absent
}

void DftShieldingStore::requestFrame(std::size_t originalIndex) {
    ASSERT_THREAD(this);
    // Idempotent: a resolved frame (resident OR known-absent) just re-announces.
    // (cache_ is GUI-thread-only in v1; add a mutex here when a prefetch worker
    // lands, exactly as the Conformation snapshot cache does.)
    if (cache_.find(originalIndex) == cache_.end())
        cache_.emplace(originalIndex, loadAndValidate(originalIndex));  // may cache null (absent/invalid)
    emit frameReady(originalIndex);
}

std::optional<double> DftShieldingStore::sample(std::size_t originalIndex, std::size_t atom,
                                                DftPart part, DftScalar scalar) const {
    const auto it = cache_.find(originalIndex);
    if (it == cache_.end() || !it->second)
        return std::nullopt;  // not resident, or resolved-absent
    const DftShieldingFrame& fr = *it->second;
    if (atom >= fr.atoms.size())
        return std::nullopt;
    const DftAtomShielding& a    = fr.atoms[atom];
    const SphericalTensor&  tens = (part == DftPart::Total) ? a.total
                                   : (part == DftPart::Dia) ? a.dia
                                                            : a.para;
    return (scalar == DftScalar::IsotropicT0) ? tens.T0 : tens.T2Magnitude();
}

std::shared_ptr<const DftShieldingFrame> DftShieldingStore::loadAndValidate(std::size_t originalIndex) const {
    const auto dirIt = dirByOriginal_.find(originalIndex);
    if (dirIt == dirByOriginal_.end())
        return nullptr;  // no job for this frame — an honest gap, not an error
    const QString jobDir = dirIt->second;
    const QString jobId  = QFileInfo(jobDir).fileName();

    // meta.json -> files.out_primary names the SUCCESSFUL run (a frame may have
    // retry .out files; never glob *_nmr.out).
    const QString metaPath = QStringLiteral("%1/%2_meta.json").arg(jobDir, jobId);
    QFile metaFile(metaPath);
    if (!metaFile.open(QIODevice::ReadOnly)) {
        report(Severity::Warning, QStringLiteral("cannot open meta.json"), metaPath);
        return nullptr;
    }
    const QJsonObject filesObj = QJsonDocument::fromJson(metaFile.readAll())
                                     .object()
                                     .value(QStringLiteral("files"))
                                     .toObject();
    const QString outRel = filesObj.value(QStringLiteral("out_primary")).toString();
    if (outRel.isEmpty()) {
        report(Severity::Warning, QStringLiteral("meta.json has no files.out_primary"), metaPath);
        return nullptr;
    }

    const QString outPath = QStringLiteral("%1/%2").arg(jobDir, outRel);
    QFile outFile(outPath);
    if (!outFile.open(QIODevice::ReadOnly)) {
        report(Severity::Warning, QStringLiteral("cannot open ORCA .out"), outPath);
        return nullptr;
    }
    // The parser is std::istream-based (Qt-free, testable); read via QFile (Qt
    // I/O) and feed an istringstream — "the app wraps it with QFile".
    const QByteArray  bytes = outFile.readAll();
    std::istringstream ss(bytes.toStdString());
    DftShieldingFrame fr = io::ParseOrcaNmrShielding(ss);

    // ---- Strict validation over the permissive parser. ----
    const std::size_t expected = protein_ ? protein_->atomCount() : 0;
    if (!fr.valid || fr.atoms.size() != expected) {
        report(Severity::Warning,
               QStringLiteral("DFT atom count %1 != topology %2 (or empty section)")
                   .arg(fr.atoms.size())
                   .arg(expected),
               outPath);
        return nullptr;
    }
    for (std::size_t i = 0; i < fr.atoms.size(); ++i) {
        const DftAtomShielding& a = fr.atoms[i];
        if (a.element == Element::Unknown) {  // a parser hole (default-filled index gap)
            report(Severity::Warning,
                   QStringLiteral("unparsed atom at index %1 (parser hole)").arg(i), outPath);
            return nullptr;
        }
        // Decomposition is linear, so the T0 identity stands in for all components.
        if (std::abs(a.total.T0 - (a.dia.T0 + a.para.T0)) > kIdentityTolPpm) {
            report(Severity::Warning,
                   QStringLiteral("atom %1 fails total==dia+para (%2 vs %3 ppm)")
                       .arg(i)
                       .arg(a.total.T0, 0, 'f', 3)
                       .arg(a.dia.T0 + a.para.T0, 0, 'f', 3),
                   outPath);
            return nullptr;
        }
    }

    qCDebug(cDft).noquote() << "loaded DFT frame |" << originalIndex << "| atoms=" << fr.atoms.size()
                            << "| from" << outRel;
    return std::make_shared<const DftShieldingFrame>(std::move(fr));
}

}  // namespace h5reader::model

#include "DftShieldingLoader.h"

#include "OrcaShieldingParser.h"

#include "../diagnostics/ErrorBus.h"
#include "../model/QtProtein.h"

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

namespace h5reader::io {

namespace {
Q_LOGGING_CATEGORY(cDft, "h5reader.dft")

using Severity = h5reader::diagnostics::Severity;

void report(Severity sev, const QString& msg, const QString& ctx) {
    h5reader::diagnostics::ErrorBus::Report(sev, QStringLiteral("DftShieldingLoader"), msg, ctx);
}

QString resolveJobDir(std::size_t originalIndex, const QString& jobsDir) {
    QDir dir(jobsDir);
    if (!dir.exists())
        return {};

    const QStringList entries = dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
    for (const QString& name : entries) {
        const int t = name.lastIndexOf(QStringLiteral("_t"));
        if (t < 0)
            continue;
        const int f = name.lastIndexOf(QStringLiteral("_f"), t);
        if (f < 0 || t <= f + 2)
            continue;
        bool ok = false;
        const qulonglong orig = name.mid(f + 2, t - (f + 2)).toULongLong(&ok);
        if (!ok || static_cast<std::size_t>(orig) != originalIndex)
            continue;
        return dir.absoluteFilePath(name);
    }
    return {};
}

// ppm tolerance for the ORCA identity total == dia + para. ORCA prints the
// tensors to ~3-4 decimals, so the reconstructed sum agrees to well under this.
constexpr double kIdentityTolPpm = 0.1;
}  // namespace

std::shared_ptr<const h5reader::model::DftShieldingFrame>
DftShieldingLoader::LoadAndValidate(std::size_t originalIndex,
                                    const QString& jobsDir,
                                    const h5reader::model::QtProtein* protein) {
    const QString jobDir = resolveJobDir(originalIndex, jobsDir);
    if (jobDir.isEmpty())
        return nullptr;  // no job for this frame -- an honest gap, not an error
    const QString jobId = QFileInfo(jobDir).fileName();

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
    // I/O) and feed an istringstream -- "the app wraps it with QFile".
    const QByteArray bytes = outFile.readAll();
    std::istringstream ss(bytes.toStdString());
    h5reader::model::DftShieldingFrame fr = ParseOrcaNmrShielding(ss);

    // ---- Strict validation over the permissive parser. ----
    const std::size_t expected = protein ? protein->atomCount() : 0;
    if (!fr.valid || fr.atoms.size() != expected) {
        report(Severity::Warning,
               QStringLiteral("DFT atom count %1 != topology %2 (or empty section)")
                   .arg(fr.atoms.size())
                   .arg(expected),
               outPath);
        return nullptr;
    }
    for (std::size_t i = 0; i < fr.atoms.size(); ++i) {
        const h5reader::model::DftAtomShielding& a = fr.atoms[i];
        if (a.element == h5reader::model::Element::Unknown) {  // a parser hole (default-filled index gap)
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
    return std::make_shared<const h5reader::model::DftShieldingFrame>(std::move(fr));
}

}  // namespace h5reader::io

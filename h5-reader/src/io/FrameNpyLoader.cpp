// FrameNpyLoader — implementation. One directory of per-atom NPYs → one
// QtConformationSnapshot, directory-agnostic (a trajectory frame_NNNNNN/ dir or
// a single-pose run root).

#include "FrameNpyLoader.h"

#include "QtFieldCatalog.gen.h"
#include "QtNpyReader.h"

#include "../diagnostics/ErrorBus.h"
#include "../model/QtConformationSnapshot.h"
#include "../model/QtProtein.h"

#include <QDir>
#include <QFileInfo>
#include <QLoggingCategory>

#include <optional>
#include <string>
#include <utility>

Q_LOGGING_CATEGORY(cFrameLoader, "h5reader.frameloader")

namespace h5reader::io {

namespace {
using h5reader::diagnostics::ErrorBus;
using h5reader::diagnostics::Severity;
}  // namespace

std::shared_ptr<h5reader::model::QtConformationSnapshot>
FrameNpyLoader::LoadSnapshotDir(const QString& dir,
                                const h5reader::model::QtProtein* protein,
                                std::size_t frameIndex,
                                double timePs) {
    QDir d(dir);
    if (!d.exists()) {
        ErrorBus::Report(Severity::Error, QStringLiteral("FrameNpyLoader"),
                         QStringLiteral("per-frame NPY directory does not exist"), dir);
        return nullptr;
    }

    auto snap = std::make_shared<h5reader::model::QtConformationSnapshot>(protein, frameIndex, timePs);
    const std::size_t atomCount = protein ? protein->atomCount() : 0;

    const QStringList npys = d.entryList(QStringList{QStringLiteral("*.npy")}, QDir::Files, QDir::Name);
    int loaded = 0;
    int skipped = 0;

    for (const QString& fname : npys) {
        const std::string stem = QFileInfo(fname).completeBaseName().toStdString();
        const std::optional<FieldKind> kind = FindFieldByStem(stem);
        if (!kind) {
            ++skipped;  // not a catalog field
            continue;
        }
        const FieldSpec& spec = FieldSpecFor(*kind);

        // Skip the structured sidecar / topology files — those are decoded by
        // QtTopologySidecar into the QtProtein spine, not stored as numeric
        // calculator columns. (Absent in a trajectory frame_NNNNNN/ dir; present
        // alongside the calc NPYs in a single-pose run root.)
        if (spec.group == FieldGroup::Topology || *kind == FieldKind::AtomsCategoryInfo) {
            ++skipped;
            continue;
        }

        QtNpyReader::WidenedArray arr = QtNpyReader::ReadArrayWidened(d.filePath(fname));
        if (!arr.ok) {
            ++skipped;  // malformed — already logged at the QtNpyReader seam; column stays absent
            continue;
        }

        // Shape into (rows, cols). A 1-D Protein-axis array (K,) is the protein
        // row, 1 x K (e.g. mopac_global); every other 1-D array stays N x 1.
        int rows = static_cast<int>(arr.rows);
        int cols = static_cast<int>(arr.cols);
        if (arr.cols == 1 && spec.axis == NativeAxis::Protein) {
            cols = static_cast<int>(arr.rows);
            rows = 1;
        }

        // Writer is definitive: when the NPY column count disagrees with the
        // catalog's (non -1) value, trust the NPY and flag the drift — this is
        // how gromacs_energy's 43-vs-catalog-42 surfaces at load.
        if (spec.cols != -1 && cols != spec.cols) {
            ErrorBus::Report(Severity::Warning, QStringLiteral("FrameNpyLoader"),
                             QStringLiteral("%1: NPY cols=%2 disagrees with catalog cols=%3 — trusting the NPY")
                                 .arg(QString::fromStdString(stem))
                                 .arg(cols)
                                 .arg(spec.cols),
                             dir);
        }
        // Atom-axis arrays should span the topology; a mismatch means dir and
        // spine disagree. Load anyway, but flag it loud.
        if (spec.axis == NativeAxis::Atom && atomCount != 0 && static_cast<std::size_t>(rows) != atomCount) {
            ErrorBus::Report(Severity::Warning, QStringLiteral("FrameNpyLoader"),
                             QStringLiteral("%1: rows=%2 != atom count %3")
                                 .arg(QString::fromStdString(stem))
                                 .arg(rows)
                                 .arg(atomCount),
                             dir);
        }

        model::NpyColumn& col = snap->mutableColumn(*kind);
        col.present = true;
        col.rows = rows;
        col.cols = cols;
        col.data = std::move(arr.data);
        ++loaded;
    }

    if (loaded == 0) {
        ErrorBus::Report(Severity::Warning, QStringLiteral("FrameNpyLoader"),
                         QStringLiteral("no recognised calculator NPYs found (skipped %1)").arg(skipped), dir);
        return nullptr;
    }

    qCInfo(cFrameLoader).noquote() << "snapshot frame" << frameIndex << "loaded" << loaded << "arrays ("
                                   << skipped << "skipped) from" << dir;
    return snap;
}

}  // namespace h5reader::io

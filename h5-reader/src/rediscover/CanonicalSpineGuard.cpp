#include "CanonicalSpineGuard.h"

#include "CanonicalSpineInputs.h"

#include "../io/QtNpyReader.h"

#include <QDir>
#include <QFileInfo>
#include <QStringList>

namespace h5reader::rediscover {
namespace {

QString cleanAbs(const QString& path) {
    return QDir::cleanPath(QFileInfo(path).absoluteFilePath());
}

bool sameOrUnder(const QString& path, const QString& root) {
    const QString p = cleanAbs(path);
    const QString r = cleanAbs(root);
    return p == r || p.startsWith(r + QStringLiteral("/"));
}

QString stemPath(const QString& dir, const QString& stem) {
    return QDir(dir).filePath(stem + QStringLiteral(".npy"));
}

bool requireFile(const QString& path, QString* err_out) {
    if (QFileInfo::exists(path) && QFileInfo(path).isFile()) return true;
    if (err_out) *err_out = QStringLiteral("canonical sentinel file is missing: %1").arg(path);
    return false;
}

}  // namespace

bool ValidateCanonical720Root(const QString& root, QString* err_out) {
    const QString clean = cleanAbs(root);
    const QString expected = cleanAbs(canonical::k720Root);
    if (clean != expected) {
        if (err_out) {
            *err_out = QStringLiteral("non-canonical 720 root: %1 (expected exactly %2)")
                           .arg(clean, expected);
        }
        return false;
    }
    if (!QFileInfo(clean).isDir()) {
        if (err_out) *err_out = QStringLiteral("canonical 720 root does not exist: %1").arg(clean);
        return false;
    }
    return true;
}

bool ValidateCanonical720Pose(const QString& poseDir, QString* err_out) {
    const QString clean = cleanAbs(poseDir);
    if (!sameOrUnder(clean, canonical::k720Root)
        || !canonical::isUnderCanonical720(clean)) {
        if (err_out) {
            *err_out = QStringLiteral("non-canonical 720 pose: %1 (must be under %2)")
                           .arg(clean, cleanAbs(canonical::k720Root));
        }
        return false;
    }
    if (!QFileInfo(clean).isDir()) {
        if (err_out) *err_out = QStringLiteral("canonical 720 pose dir does not exist: %1").arg(clean);
        return false;
    }

    QDir d(clean);
    const QStringList npys = d.entryList(QStringList{QStringLiteral("*.npy")}, QDir::Files, QDir::Name);
    if (npys.size() != canonical::k720PoseNpyCount) {
        if (err_out) {
            *err_out = QStringLiteral("canonical 720 schema sentinel failed for %1: saw %2 npys, expected %3")
                           .arg(clean)
                           .arg(npys.size())
                           .arg(canonical::k720PoseNpyCount);
        }
        return false;
    }

    for (const QString& stem : canonical::kCurrentSchemaSentinelStems) {
        if (!requireFile(stemPath(clean, stem), err_out)) return false;
    }

    const QString ringPath = stemPath(clean, QStringLiteral("ring_contributions"));
    if (!requireFile(ringPath, err_out)) return false;
    const io::QtNpyReader::WidenedArray ring = io::QtNpyReader::ReadArrayWidened(ringPath);
    if (!ring.ok) {
        if (err_out) *err_out = ring.error;
        return false;
    }
    if (ring.cols != static_cast<std::size_t>(canonical::kRingContributionsCols)) {
        if (err_out) {
            *err_out = QStringLiteral("canonical 720 schema sentinel failed for %1: ring_contributions has %2 cols, expected %3")
                           .arg(clean)
                           .arg(ring.cols)
                           .arg(canonical::kRingContributionsCols);
        }
        return false;
    }
    return true;
}

bool ValidateCanonical1p9jRoot(const QString& rootOrLgsPath, QString* err_out) {
    QFileInfo fi(rootOrLgsPath);
    const QString root = fi.isFile() && fi.suffix().compare(QStringLiteral("lgs"), Qt::CaseInsensitive) == 0
                             ? fi.absoluteDir().absolutePath()
                             : fi.absoluteFilePath();
    const QString clean = cleanAbs(root);
    const QString expected = cleanAbs(canonical::k1p9jRoot);
    if (clean != expected || !canonical::isCanonical1p9j(clean)) {
        if (err_out) {
            *err_out = QStringLiteral("non-canonical 1P9J root: %1 (expected exactly %2)")
                           .arg(clean, expected);
        }
        return false;
    }
    QDir d(clean);
    if (!d.exists()) {
        if (err_out) *err_out = QStringLiteral("canonical 1P9J root does not exist: %1").arg(clean);
        return false;
    }
    for (const QString& name : {QStringLiteral("trajectory.h5"),
                                QStringLiteral("extraction_manifest.json"),
                                QStringLiteral("atoms_category_info.npy"),
                                QStringLiteral("residues.npy"),
                                QStringLiteral("bonds.npy"),
                                QStringLiteral("rings.npy"),
                                QStringLiteral("ring_membership.npy")}) {
        if (!requireFile(d.filePath(name), err_out)) return false;
    }
    if (!QFileInfo(d.filePath(QStringLiteral("npys"))).isDir()) {
        if (err_out) *err_out = QStringLiteral("canonical 1P9J root is missing npys/: %1").arg(clean);
        return false;
    }
    return true;
}

}  // namespace h5reader::rediscover

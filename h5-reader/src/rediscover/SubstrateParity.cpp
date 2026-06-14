#include "SubstrateParity.h"

#include "SphericalBasis.h"

#include "../io/QtFieldCatalog.gen.h"
#include "../io/QtNpyReader.h"

#include <QFile>
#include <QTextStream>

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <vector>

namespace h5reader::rediscover {
namespace {

struct RowKey {
    std::size_t row_id = 0;
    std::size_t atom = 0;
    std::size_t h5_row = 0;
    std::size_t frame_slot = 0;
    std::size_t original_index = 0;
    bool dft_present = false;
};

std::optional<std::size_t> parseSize(const QString& s) {
    bool ok = false;
    const qulonglong v = s.trimmed().toULongLong(&ok);
    if (!ok) return std::nullopt;
    return static_cast<std::size_t>(v);
}

std::optional<int> parseInt(const QString& s) {
    bool ok = false;
    const int v = s.trimmed().toInt(&ok);
    if (!ok) return std::nullopt;
    return v;
}

std::unordered_map<QString, int> headerMap(const QString& header) {
    const QStringList cols = header.split(QLatin1Char(','));
    std::unordered_map<QString, int> out;
    for (int i = 0; i < cols.size(); ++i) out[cols[i]] = i;
    return out;
}

std::optional<QString> requireColumn(const std::unordered_map<QString, int>& cols,
                                     const QString& name) {
    if (cols.find(name) != cols.end()) return std::nullopt;
    return QStringLiteral("missing required CSV column '%1'").arg(name);
}

std::vector<RowKey> readRows(const QString& path,
                             bool requireFrameSlot,
                             SubstrateParityStats* stats) {
    std::vector<RowKey> rows;
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
        stats->errors << QStringLiteral("could not open rows CSV: %1").arg(path);
        return rows;
    }
    QTextStream in(&f);
    if (in.atEnd()) {
        stats->errors << QStringLiteral("empty rows CSV: %1").arg(path);
        return rows;
    }
    const QString header = in.readLine();
    const auto cols = headerMap(header);
    const QStringList required = requireFrameSlot
        ? QStringList{QStringLiteral("row_id"), QStringLiteral("atom_index"),
                      QStringLiteral("h5_row"), QStringLiteral("frame_slot"),
                      QStringLiteral("original_frame_index"), QStringLiteral("dft_present")}
        : QStringList{QStringLiteral("row_id"), QStringLiteral("atom_index"),
                      QStringLiteral("h5_row"), QStringLiteral("original_index"),
                      QStringLiteral("dft_present")};
    for (const QString& name : required) {
        if (const auto err = requireColumn(cols, name)) {
            stats->errors << *err;
            return rows;
        }
    }

    auto col = [&](const QString& name) -> int { return cols.at(name); };
    const int rowIdCol = col(QStringLiteral("row_id"));
    const int atomCol = col(QStringLiteral("atom_index"));
    const int h5RowCol = col(QStringLiteral("h5_row"));
    const int frameSlotCol = requireFrameSlot ? col(QStringLiteral("frame_slot")) : -1;
    const int originalCol = requireFrameSlot ? col(QStringLiteral("original_frame_index"))
                                             : col(QStringLiteral("original_index"));
    const int dftCol = col(QStringLiteral("dft_present"));

    int lineNo = 1;
    while (!in.atEnd()) {
        ++lineNo;
        const QString line = in.readLine();
        if (line.trimmed().isEmpty()) continue;
        const QStringList parts = line.split(QLatin1Char(','), Qt::KeepEmptyParts);
        const int minCols = std::max({rowIdCol, atomCol, h5RowCol, originalCol, dftCol, frameSlotCol}) + 1;
        if (parts.size() < minCols) {
            stats->errors << QStringLiteral("%1:%2 has %3 columns, need at least %4")
                                 .arg(path)
                                 .arg(lineNo)
                                 .arg(parts.size())
                                 .arg(minCols);
            continue;
        }
        const auto rowId = parseSize(parts[rowIdCol]);
        const auto atom = parseSize(parts[atomCol]);
        const auto h5Row = parseSize(parts[h5RowCol]);
        const auto frameSlot = requireFrameSlot ? parseSize(parts[frameSlotCol])
                                                : std::optional<std::size_t>(rows.size());
        const auto orig = parseSize(parts[originalCol]);
        const auto dft = parseInt(parts[dftCol]);
        if (!rowId || !atom || !h5Row || !frameSlot || !orig || !dft) {
            stats->errors << QStringLiteral("%1:%2 has an unparsable row key").arg(path).arg(lineNo);
            continue;
        }
        rows.push_back(RowKey{*rowId, *atom, *h5Row, *frameSlot, *orig, *dft != 0});
    }
    stats->rows_seen = rows.size();
    return rows;
}

bool near(double a, double b, const SubstrateParityOptions& options) {
    if (std::isnan(a) && std::isnan(b)) return true;
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    const double diff = std::abs(a - b);
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return diff <= options.abs_tolerance || diff <= options.rel_tolerance * scale;
}

void compareValue(SubstrateParityStats* stats,
                  const QString& label,
                  double actual,
                  double expected,
                  const SubstrateParityOptions& options) {
    if (!near(actual, expected, options)) {
        stats->errors << QStringLiteral("%1 mismatch: actual=%2 expected=%3")
                             .arg(label)
                             .arg(actual, 0, 'g', 17)
                             .arg(expected, 0, 'g', 17);
    }
}

io::QtNpyReader::WidenedArray readNpy(const QString& path, SubstrateParityStats* stats) {
    const io::QtNpyReader::WidenedArray a = io::QtNpyReader::ReadArrayWidened(path);
    if (!a.ok) stats->errors << a.error;
    return a;
}

io::QtNpyReader::NumericArray readNumericNpy(const QString& path, SubstrateParityStats* stats) {
    const io::QtNpyReader::NumericArray a = io::QtNpyReader::ReadNumericArrayWidened(path);
    if (!a.ok) stats->errors << a.error;
    return a;
}

std::optional<double> npyAt(const io::QtNpyReader::WidenedArray& a,
                            std::size_t row,
                            std::size_t col,
                            const QString& label,
                            SubstrateParityStats* stats) {
    if (row >= a.rows || col >= a.cols) {
        stats->errors << QStringLiteral("%1 out of range row=%2 col=%3 shape=(%4,%5)")
                             .arg(label)
                             .arg(static_cast<qulonglong>(row))
                             .arg(static_cast<qulonglong>(col))
                             .arg(static_cast<qulonglong>(a.rows))
                             .arg(static_cast<qulonglong>(a.cols));
        return std::nullopt;
    }
    return a.data[row * a.cols + col];
}

std::optional<double> npyAt3(const io::QtNpyReader::NumericArray& a,
                             std::size_t row,
                             std::size_t i,
                             std::size_t j,
                             const QString& label,
                             SubstrateParityStats* stats) {
    if (a.shape.size() != 3 || row >= a.shape[0] || i >= a.shape[1] || j >= a.shape[2]) {
        QString shape;
        for (std::size_t k = 0; k < a.shape.size(); ++k) {
            if (k) shape += QLatin1Char(',');
            shape += QString::number(static_cast<qulonglong>(a.shape[k]));
        }
        stats->errors << QStringLiteral("%1 out of range row=%2 i=%3 j=%4 shape=(%5)")
                             .arg(label)
                             .arg(static_cast<qulonglong>(row))
                             .arg(static_cast<qulonglong>(i))
                             .arg(static_cast<qulonglong>(j))
                             .arg(shape);
        return std::nullopt;
    }
    return a.data[(row * a.shape[1] + i) * a.shape[2] + j];
}

void checkShape(const io::QtNpyReader::WidenedArray& a,
                const QString& label,
                std::size_t rows,
                std::size_t cols,
                SubstrateParityStats* stats) {
    if (!a.ok) return;
    if (a.rows != rows || a.cols != cols) {
        stats->errors << QStringLiteral("%1 shape mismatch: got (%2,%3), expected (%4,%5)")
                             .arg(label)
                             .arg(static_cast<qulonglong>(a.rows))
                             .arg(static_cast<qulonglong>(a.cols))
                             .arg(static_cast<qulonglong>(rows))
                             .arg(static_cast<qulonglong>(cols));
    }
}

void checkShape3(const io::QtNpyReader::NumericArray& a,
                 const QString& label,
                 std::size_t rows,
                 std::size_t d1,
                 std::size_t d2,
                 SubstrateParityStats* stats) {
    if (!a.ok) return;
    if (a.shape.size() != 3 || a.shape[0] != rows || a.shape[1] != d1 || a.shape[2] != d2) {
        QString shape;
        for (std::size_t k = 0; k < a.shape.size(); ++k) {
            if (k) shape += QLatin1Char(',');
            shape += QString::number(static_cast<qulonglong>(a.shape[k]));
        }
        stats->errors << QStringLiteral("%1 shape mismatch: got (%2), expected (%3,%4,%5)")
                             .arg(label)
                             .arg(shape)
                             .arg(static_cast<qulonglong>(rows))
                             .arg(static_cast<qulonglong>(d1))
                             .arg(static_cast<qulonglong>(d2));
    }
}

SubstrateParityStats auditTargets(const QString& rowsPath,
                                  bool perAtomContract,
                                  std::size_t atomCount,
                                  const QString& t0Path,
                                  const QString& t1Path,
                                  const QString& t2Path,
                                  const QString& rawPath,
                                  const ProducerTargetLookup& producerTarget,
                                  const SubstrateParityOptions& options) {
    SubstrateParityStats stats;
    const std::vector<RowKey> rows = readRows(rowsPath, perAtomContract, &stats);
    if (!stats.errors.isEmpty()) return stats;

    const io::QtNpyReader::WidenedArray t0 = t0Path.isEmpty()
        ? io::QtNpyReader::WidenedArray{}
        : readNpy(t0Path, &stats);
    const io::QtNpyReader::WidenedArray t1 = t1Path.isEmpty()
        ? io::QtNpyReader::WidenedArray{}
        : readNpy(t1Path, &stats);
    const io::QtNpyReader::WidenedArray t2 = readNpy(t2Path, &stats);
    const io::QtNpyReader::NumericArray raw = rawPath.isEmpty()
        ? io::QtNpyReader::NumericArray{}
        : readNumericNpy(rawPath, &stats);
    if (!stats.errors.isEmpty()) return stats;

    if (!t0Path.isEmpty()) checkShape(t0, t0Path, rows.size(), 1, &stats);
    if (!t1Path.isEmpty()) checkShape(t1, t1Path, rows.size(), 3, &stats);
    checkShape(t2, t2Path, rows.size(), 5, &stats);
    if (!rawPath.isEmpty()) checkShape3(raw, rawPath, rows.size(), 3, 3, &stats);
    if (!stats.errors.isEmpty()) return stats;

    const std::size_t limit = options.max_rows == 0
        ? rows.size()
        : std::min(rows.size(), options.max_rows);
    for (std::size_t i = 0; i < limit; ++i) {
        const RowKey& row = rows[i];
        if (row.row_id != i) {
            stats.errors << QStringLiteral("row %1 row_id=%2, expected sidecar row id %1")
                                .arg(static_cast<qulonglong>(i))
                                .arg(static_cast<qulonglong>(row.row_id));
        }
        if (perAtomContract) {
            const std::size_t expected = row.frame_slot * atomCount + row.atom;
            if (row.row_id != expected) {
                stats.errors << QStringLiteral("per_atom row contract mismatch at sidecar row %1: "
                                               "row_id=%2 expected frame_slot*n_atoms+atom_index=%3")
                                    .arg(static_cast<qulonglong>(i))
                                    .arg(static_cast<qulonglong>(row.row_id))
                                    .arg(static_cast<qulonglong>(expected));
            }
        }
        if (!row.dft_present) continue;

        QString reason;
        const std::optional<model::Mat3> expectedRaw = producerTarget(row.atom, row.h5_row, &reason);
        if (!expectedRaw) {
            stats.errors << QStringLiteral("producer target absent for atom=%1 h5_row=%2: %3")
                                .arg(static_cast<qulonglong>(row.atom))
                                .arg(static_cast<qulonglong>(row.h5_row))
                                .arg(reason);
            continue;
        }
        const model::SphericalTensor expected = DecomposeLibrary(*expectedRaw);

        if (!t0Path.isEmpty()) {
            const std::optional<double> actual = npyAt(t0, i, 0, t0Path, &stats);
            if (actual) compareValue(&stats, QStringLiteral("target_T0 row %1").arg(i),
                                     *actual, expected.T0, options);
            ++stats.target_T0_checked;
        }
        if (!t1Path.isEmpty()) {
            for (std::size_t c = 0; c < 3; ++c) {
                const std::optional<double> actual = npyAt(t1, i, c, t1Path, &stats);
                if (actual) compareValue(&stats, QStringLiteral("target_T1 row %1 comp %2").arg(i).arg(c),
                                         *actual, expected.T1[c], options);
                ++stats.target_T1_checked;
            }
        }
        for (std::size_t c = 0; c < 5; ++c) {
            const std::optional<double> actual = npyAt(t2, i, c, t2Path, &stats);
            if (actual) compareValue(&stats, QStringLiteral("target_T2 row %1 comp %2").arg(i).arg(c),
                                     *actual, expected.T2[c], options);
            ++stats.target_T2_checked;
        }
        if (!rawPath.isEmpty()) {
            for (std::size_t r = 0; r < 3; ++r)
                for (std::size_t c = 0; c < 3; ++c) {
                    const std::optional<double> actual = npyAt3(raw, i, r, c, rawPath, &stats);
                    if (actual) {
                        const double expectedComp = (*expectedRaw)(static_cast<int>(r),
                                                                   static_cast<int>(c));
                        compareValue(&stats,
                                     QStringLiteral("target_raw row %1 comp %2,%3").arg(i).arg(r).arg(c),
                                     *actual, expectedComp, options);
                    }
                    ++stats.target_raw_components_checked;
                }
        }
        ++stats.rows_checked;
    }
    return stats;
}

ProducerTargetLookup catalogTargetLookup(const Body& body) {
    return [&body](std::size_t atom, std::size_t frame, QString* reason_out)
        -> std::optional<model::Mat3> {
        model::Mat3 m = model::Mat3::Zero();
        for (int c = 0; c < 9; ++c) {
            QString reason;
            const std::optional<double> v =
                body.catalog.value(body, io::FieldKind::OrcaTotal, atom, frame, c, &reason);
            if (!v) {
                if (reason_out) *reason_out = reason;
                return std::nullopt;
            }
            m(c / 3, c % 3) = *v;
        }
        if (reason_out) reason_out->clear();
        return m;
    };
}

}  // namespace

SubstrateParityStats AuditPerAtomSubstrateSidecars(
    const QString& outDir,
    std::size_t atomCount,
    const ProducerTargetLookup& producerTarget,
    const SubstrateParityOptions& options) {
    return auditTargets(QStringLiteral("%1/per_atom_substrate_rows.csv").arg(outDir),
                        true,
                        atomCount,
                        QStringLiteral("%1/per_atom_substrate_target_T0.npy").arg(outDir),
                        QStringLiteral("%1/per_atom_substrate_target_T1.npy").arg(outDir),
                        QStringLiteral("%1/per_atom_substrate_target_T2.npy").arg(outDir),
                        QString(),
                        producerTarget,
                        options);
}

SubstrateParityStats AuditAllAtomEquivariantSidecars(
    const QString& outDir,
    const ProducerTargetLookup& producerTarget,
    const SubstrateParityOptions& options) {
    return auditTargets(QStringLiteral("%1/all_atom_equivariant_targets.csv").arg(outDir),
                        false,
                        0,
                        QStringLiteral("%1/all_atom_equivariant_target_sigma_iso.npy").arg(outDir),
                        QString(),
                        QStringLiteral("%1/all_atom_equivariant_target_T2.npy").arg(outDir),
                        QStringLiteral("%1/all_atom_equivariant_target_raw.npy").arg(outDir),
                        producerTarget,
                        options);
}

SubstrateParityStats AuditPerAtomSubstrateAgainstProducer(
    const Body& body,
    const QString& outDir,
    const SubstrateParityOptions& options) {
    const std::size_t atomCount = body.run.protein ? body.run.protein->atomCount() : 0;
    return AuditPerAtomSubstrateSidecars(outDir, atomCount, catalogTargetLookup(body), options);
}

SubstrateParityStats AuditAllAtomEquivariantAgainstProducer(
    const Body& body,
    const QString& outDir,
    const SubstrateParityOptions& options) {
    return AuditAllAtomEquivariantSidecars(outDir, catalogTargetLookup(body), options);
}

}  // namespace h5reader::rediscover

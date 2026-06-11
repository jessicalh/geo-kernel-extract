#include "SpineReachabilityProbe.h"

#include "Catalog.h"
#include "RunData.h"
#include "ScopedProducerCatalog.h"
#include "Verbs.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>

#include <algorithm>
#include <cstdint>
#include <map>
#include <numeric>

namespace h5reader::rediscover {

namespace {

using CountMap = std::map<QString, std::uint64_t>;

QJsonObject countsJson(const CountMap& counts) {
    QJsonObject out;
    for (const auto& item : counts) {
        out.insert(item.first, static_cast<double>(item.second));
    }
    return out;
}

bool writeJson(const QString& path, const QJsonObject& obj, QString* err_out) {
    QSaveFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot open %1 for write").arg(path);
        return false;
    }
    f.write(QJsonDocument(obj).toJson(QJsonDocument::Indented));
    if (!f.commit()) {
        if (err_out) *err_out = QStringLiteral("cannot commit %1").arg(path);
        return false;
    }
    return true;
}

std::vector<std::size_t> frameDomain(const RunData& run, const FieldAccessSpec& access) {
    if (!access.available) return {};
    if (access.provider == FieldProvider::SparseDftByOriginal) {
        return run.frameMap.dftRows();
    }
    std::vector<std::size_t> frames(access.frames);
    std::iota(frames.begin(), frames.end(), std::size_t{0});
    return frames;
}

QString sourceFindingForAxis(io::NativeAxis axis) {
    switch (axis) {
    case io::NativeAxis::Atom:
        return QStringLiteral("typed_atom_index");
    case io::NativeAxis::RingContributionPair:
        return QStringLiteral("native_pair_rows; spatial source finding via RingCenters cloud plus topology");
    case io::NativeAxis::Bond:
        return QStringLiteral("native_bond_rows; source finding via topology plus AllBondMidpoints/BondMidpoints clouds");
    case io::NativeAxis::MOPACBondNeighborPair:
    case io::NativeAxis::MOPACUniquePair:
        return QStringLiteral("native_mopac_pair_rows; source finding via topology/bond midpoint clouds when queried");
    case io::NativeAxis::Residue:
        return QStringLiteral("typed_residue_index");
    case io::NativeAxis::AromaticRing:
    case io::NativeAxis::SaturatedRing:
    case io::NativeAxis::Ring:
    case io::NativeAxis::RingMembership:
        return QStringLiteral("typed_topology_ring_index");
    case io::NativeAxis::Protein:
        return QStringLiteral("protein_scope");
    default:
        return QStringLiteral("unsupported");
    }
}

std::map<QString, double> loadFlatCoverageCounts(const QString& outDir, bool* compared) {
    std::map<QString, double> counts;
    if (compared) *compared = false;
    const QString path = QDir(outDir).filePath(QStringLiteral("catalog_coverage.json"));
    if (!QFileInfo::exists(path)) return counts;
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return counts;
    const QJsonDocument doc = QJsonDocument::fromJson(f.readAll());
    const QJsonArray fields = doc.object().value(QStringLiteral("fields")).toArray();
    for (const QJsonValue& v : fields) {
        const QJsonObject obj = v.toObject();
        const QString stem = obj.value(QStringLiteral("stem")).toString();
        double total = 0.0;
        const QJsonValue populated = obj.value(QStringLiteral("populated_counts"));
        if (populated.isArray()) {
            for (const QJsonValue& c : populated.toArray()) total += c.toDouble();
        } else if (populated.isDouble()) {
            total += populated.toDouble();
        }
        if (!stem.isEmpty()) counts[stem] = total;
    }
    if (compared) *compared = true;
    return counts;
}

}  // namespace

SpineProbeDatasetResult RunSpineReachabilityProbe(const Body& body,
                                                  const QString& datasetLabel,
                                                  const QString& outDir,
                                                  QString* err_out) {
    SpineProbeDatasetResult result;
    result.dataset_label = datasetLabel;
    QDir().mkpath(outDir);
    result.manifest_path =
        QDir(outDir).filePath(QStringLiteral("spine_reachability_%1.json").arg(datasetLabel));

    bool comparedFlat = false;
    const std::map<QString, double> flatCounts = loadFlatCoverageCounts(outDir, &comparedFlat);

    QJsonArray fields;
    for (const io::FieldSpec* spec : ScopedProducerCatalog()) {
        const FieldAccessSpec& access = body.catalog.fieldSpec(spec->kind);
        const std::vector<std::size_t> frames = frameDomain(body.run, access);
        CountMap absenceCounts;
        std::uint64_t presentCount = 0;
        std::uint64_t readCount = 0;
        std::uint64_t expectedCount = 0;
        bool failed = false;

        if (!access.available) {
            failed = true;
            absenceCounts[access.absence_reason.isEmpty() ? QStringLiteral("not-produced-in-dataset")
                                                          : access.absence_reason] = 1;
        } else if (access.structured) {
            expectedCount = static_cast<std::uint64_t>(access.native_rows * frames.size());
            for (std::size_t frame : frames) {
                for (std::size_t row = 0; row < access.native_rows; ++row) {
                    const FieldPresence p = verbs::present(body, spec->kind, row, frame, 0);
                    if (p.present) {
                        ++presentCount;
                        ++readCount;
                    } else {
                        ++absenceCounts[p.reason.isEmpty() ? QStringLiteral("not-present") : p.reason];
                    }
                }
            }
            failed = readCount != expectedCount;
        } else {
            expectedCount = static_cast<std::uint64_t>(access.native_rows)
                            * static_cast<std::uint64_t>(frames.size())
                            * static_cast<std::uint64_t>(access.components);
            for (std::size_t frame : frames) {
                for (std::size_t row = 0; row < access.native_rows; ++row) {
                    for (std::size_t component = 0; component < access.components; ++component) {
                        const FieldPresence p = verbs::present(body, spec->kind, row, frame, component);
                        if (!p.present) {
                            ++absenceCounts[p.reason.isEmpty() ? QStringLiteral("not-present") : p.reason];
                            continue;
                        }
                        ++presentCount;
                        const std::optional<double> v = verbs::value(body, spec->kind, row, frame, component);
                        if (v) {
                            ++readCount;
                        } else {
                            ++absenceCounts[QStringLiteral("read-returned-null")];
                        }
                    }
                }
            }
            failed = readCount != expectedCount;
        }

        const auto flatIt = flatCounts.find(access.stem);
        const double flatPopulated = flatIt == flatCounts.end() ? 0.0 : flatIt->second;
        if (comparedFlat && static_cast<double>(readCount) < flatPopulated) {
            failed = true;
            absenceCounts[QStringLiteral("flat-coverage-regression")] =
                static_cast<std::uint64_t>(flatPopulated - static_cast<double>(readCount));
        }

        QJsonObject f;
        f.insert(QStringLiteral("stem"), access.stem);
        f.insert(QStringLiteral("group"), FieldGroupName(spec->group));
        f.insert(QStringLiteral("axis"), NativeAxisName(spec->axis));
        f.insert(QStringLiteral("provider"), FieldProviderName(access.provider));
        f.insert(QStringLiteral("residence"), ArrayResidenceName(access.residence));
        f.insert(QStringLiteral("structured"), access.structured);
        f.insert(QStringLiteral("available"), access.available);
        f.insert(QStringLiteral("absence_reason"), access.absence_reason);
        f.insert(QStringLiteral("native_rows"), static_cast<double>(access.native_rows));
        f.insert(QStringLiteral("frames"), static_cast<double>(frames.size()));
        f.insert(QStringLiteral("frame_domain"),
                 access.provider == FieldProvider::SparseDftByOriginal
                     ? QStringLiteral("dft_frame_rows")
                     : QStringLiteral("resident_frames"));
        f.insert(QStringLiteral("components"), static_cast<double>(access.components));
        f.insert(QStringLiteral("expected_reads"), static_cast<double>(expectedCount));
        f.insert(QStringLiteral("present_reads"), static_cast<double>(presentCount));
        f.insert(QStringLiteral("read_values"), static_cast<double>(readCount));
        f.insert(QStringLiteral("absence_counts"), countsJson(absenceCounts));
        f.insert(QStringLiteral("source_finding"), sourceFindingForAxis(spec->axis));
        f.insert(QStringLiteral("flat_populated_count"), flatPopulated);
        f.insert(QStringLiteral("passed"), !failed);
        fields.push_back(f);

        ++result.field_count;
        if (failed) {
            ++result.failed_fields;
            result.failed_stems.push_back(access.stem);
        }
    }

    QJsonArray failed;
    for (const QString& stem : result.failed_stems) failed.push_back(stem);
    result.passed = result.failed_fields == 0;

    QJsonObject doc;
    doc.insert(QStringLiteral("schema_version"), 1);
    doc.insert(QStringLiteral("dataset_label"), datasetLabel);
    doc.insert(QStringLiteral("dataset_id"), body.run.manifest.dataset_id);
    doc.insert(QStringLiteral("field_count"), static_cast<double>(result.field_count));
    doc.insert(QStringLiteral("failed_fields"), static_cast<double>(result.failed_fields));
    doc.insert(QStringLiteral("passed"), result.passed);
    doc.insert(QStringLiteral("flat_coverage_compared"), comparedFlat);
    doc.insert(QStringLiteral("failed_stems"), failed);
    doc.insert(QStringLiteral("fields"), fields);

    if (!writeJson(result.manifest_path, doc, err_out)) {
        result.passed = false;
    }
    return result;
}

}  // namespace h5reader::rediscover

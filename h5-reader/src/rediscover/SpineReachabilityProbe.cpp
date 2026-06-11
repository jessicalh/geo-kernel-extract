#include "SpineReachabilityProbe.h"

#include "AnalysisBody.h"
#include "CanonicalSpineGuard.h"
#include "Catalog.h"
#include "ResidentIndexes.h"
#include "RunData.h"
#include "ScopedProducerCatalog.h"
#include "StaticRunData.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QTextStream>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <vector>

namespace h5reader::rediscover {
namespace {

QString qsv(std::string_view s) {
    return QString::fromUtf8(s.data(), static_cast<qsizetype>(s.size()));
}

QString fieldStem(const io::FieldSpec& spec) { return qsv(spec.stem); }

QString axisName(io::NativeAxis a) {
    switch (a) {
    case io::NativeAxis::Atom: return QStringLiteral("atom");
    case io::NativeAxis::RingContributionPair: return QStringLiteral("ring_contribution_pair");
    case io::NativeAxis::AromaticRing: return QStringLiteral("aromatic_ring");
    case io::NativeAxis::RediscoverSourceRow: return QStringLiteral("rediscover_source_row");
    case io::NativeAxis::RediscoverAggregatedRow: return QStringLiteral("rediscover_aggregated_row");
    case io::NativeAxis::RediscoverTargetRow: return QStringLiteral("rediscover_target_row");
    case io::NativeAxis::Protein: return QStringLiteral("protein");
    case io::NativeAxis::Bond: return QStringLiteral("bond");
    case io::NativeAxis::MOPACBondNeighborPair: return QStringLiteral("mopac_bond_neighbor_pair");
    case io::NativeAxis::MOPACUniquePair: return QStringLiteral("mopac_unique_pair");
    case io::NativeAxis::MutationMatchPair: return QStringLiteral("mutation_match_pair");
    case io::NativeAxis::Residue: return QStringLiteral("residue");
    case io::NativeAxis::SaturatedRing: return QStringLiteral("saturated_ring");
    case io::NativeAxis::Ring: return QStringLiteral("ring");
    case io::NativeAxis::RingMembership: return QStringLiteral("ring_membership");
    }
    return QStringLiteral("unknown");
}

bool structuredField(const io::FieldSpec& spec) {
    return spec.kind == io::FieldKind::AtomsCategoryInfo
           || spec.group == io::FieldGroup::Topology;
}

bool dftField(io::FieldKind kind) {
    return kind == io::FieldKind::OrcaTotal
           || kind == io::FieldKind::OrcaDiamagnetic
           || kind == io::FieldKind::OrcaParamagnetic;
}

struct FlatFieldCoverage {
    QString representation;
    std::vector<std::size_t> counts;
};

using CoverageMap = std::map<int, FlatFieldCoverage>;

bool loadCoverage(const QString& path, CoverageMap* out, QString* err_out) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) {
        if (err_out) *err_out = QStringLiteral("cannot open flat coverage: %1").arg(path);
        return false;
    }
    const QJsonDocument doc = QJsonDocument::fromJson(f.readAll());
    if (!doc.isObject()) {
        if (err_out) *err_out = QStringLiteral("flat coverage is not a JSON object: %1").arg(path);
        return false;
    }
    const QJsonArray fields = doc.object().value(QStringLiteral("fields")).toArray();
    if (fields.size() != static_cast<int>(ScopedProducerCatalog().size())) {
        if (err_out) {
            *err_out = QStringLiteral("flat coverage field count %1 != scoped catalog %2: %3")
                           .arg(fields.size())
                           .arg(ScopedProducerCatalog().size())
                           .arg(path);
        }
        return false;
    }
    for (const QJsonValue& v : fields) {
        const QJsonObject o = v.toObject();
        const QString stem = o.value(QStringLiteral("stem")).toString();
        const std::optional<io::FieldKind> kind = io::FindFieldByStem(stem.toStdString());
        if (!kind) {
            if (err_out) *err_out = QStringLiteral("flat coverage has unknown stem %1").arg(stem);
            return false;
        }
        FlatFieldCoverage c;
        c.representation = o.value(QStringLiteral("representation")).toString();
        for (const QJsonValue& cv : o.value(QStringLiteral("populated_counts")).toArray())
            c.counts.push_back(static_cast<std::size_t>(cv.toDouble()));
        (*out)[static_cast<int>(*kind)] = std::move(c);
    }
    return true;
}

struct FieldDatasetStats {
    QString provider;
    QString nativeAxis;
    std::size_t runs = 0;
    std::size_t framesVisited = 0;
    std::size_t frameExtentMin = std::numeric_limits<std::size_t>::max();
    std::size_t frameExtentMax = 0;
    std::size_t nativeRowsTotal = 0;
    std::size_t nativeRowsMin = std::numeric_limits<std::size_t>::max();
    std::size_t nativeRowsMax = 0;
    std::size_t componentsMax = 0;
    std::size_t valuesRead = 0;
    std::size_t failures = 0;
    std::vector<std::size_t> finiteCounts;
    QStringList errors;
    QStringList absenceReasons;
    QString flatRepresentation;
    bool flatCompared = false;
    bool flatOk = false;
};

void appendError(FieldDatasetStats* s, const QString& err) {
    ++s->failures;
    if (s->errors.size() < 8) s->errors.push_back(err);
}

QJsonArray countsJson(const std::vector<std::size_t>& counts) {
    QJsonArray a;
    for (std::size_t c : counts) a.push_back(static_cast<double>(c));
    return a;
}

QString extentString(const FieldDatasetStats& s) {
    const QString rows = s.nativeRowsMin == std::numeric_limits<std::size_t>::max()
                             ? QStringLiteral("0")
                             : QStringLiteral("%1..%2 total=%3")
                                   .arg(s.nativeRowsMin)
                                   .arg(s.nativeRowsMax)
                                   .arg(s.nativeRowsTotal);
    QString out = QStringLiteral("runs=%1 frames=%2..%3 visited=%4 native_rows=%5 comps=%6 values=%7")
        .arg(s.runs)
        .arg(s.frameExtentMin == std::numeric_limits<std::size_t>::max() ? 0 : s.frameExtentMin)
        .arg(s.frameExtentMax)
        .arg(s.framesVisited)
        .arg(rows)
        .arg(s.componentsMax)
        .arg(s.valuesRead);
    if (!s.absenceReasons.isEmpty())
        out += QStringLiteral("; absence=`%1`").arg(s.absenceReasons.front());
    return out;
}

QJsonObject statsJson(const FieldDatasetStats& s) {
    QJsonObject o;
    o.insert(QStringLiteral("provider"), s.provider);
    o.insert(QStringLiteral("native_axis"), s.nativeAxis);
    o.insert(QStringLiteral("runs"), static_cast<double>(s.runs));
    o.insert(QStringLiteral("frames_visited"), static_cast<double>(s.framesVisited));
    o.insert(QStringLiteral("frame_extent_min"),
             s.frameExtentMin == std::numeric_limits<std::size_t>::max()
                 ? 0.0
                 : static_cast<double>(s.frameExtentMin));
    o.insert(QStringLiteral("frame_extent_max"), static_cast<double>(s.frameExtentMax));
    o.insert(QStringLiteral("native_rows_total"), static_cast<double>(s.nativeRowsTotal));
    o.insert(QStringLiteral("native_rows_min"),
             s.nativeRowsMin == std::numeric_limits<std::size_t>::max()
                 ? 0.0
                 : static_cast<double>(s.nativeRowsMin));
    o.insert(QStringLiteral("native_rows_max"), static_cast<double>(s.nativeRowsMax));
    o.insert(QStringLiteral("components_max"), static_cast<double>(s.componentsMax));
    o.insert(QStringLiteral("values_read"), static_cast<double>(s.valuesRead));
    o.insert(QStringLiteral("finite_counts"), countsJson(s.finiteCounts));
    o.insert(QStringLiteral("failures"), static_cast<double>(s.failures));
    QJsonArray errors;
    for (const QString& e : s.errors) errors.push_back(e);
    o.insert(QStringLiteral("errors"), errors);
    QJsonArray absenceReasons;
    for (const QString& e : s.absenceReasons) absenceReasons.push_back(e);
    o.insert(QStringLiteral("absence_reasons"), absenceReasons);
    o.insert(QStringLiteral("flat_representation"), s.flatRepresentation);
    o.insert(QStringLiteral("flat_coverage_compared"), s.flatCompared);
    o.insert(QStringLiteral("flat_coverage_ok"), s.flatOk);
    return o;
}

void compareFlat(const io::FieldSpec& spec,
                 const CoverageMap& coverage,
                 FieldDatasetStats* s) {
    const auto it = coverage.find(static_cast<int>(spec.kind));
    if (it == coverage.end()) {
        appendError(s, QStringLiteral("%1: flat coverage missing").arg(fieldStem(spec)));
        s->flatCompared = false;
        s->flatOk = false;
        return;
    }
    s->flatCompared = true;
    s->flatRepresentation = it->second.representation;
    bool ok = true;
    if (s->finiteCounts.size() < it->second.counts.size())
        s->finiteCounts.resize(it->second.counts.size(), 0);
    for (std::size_t i = 0; i < it->second.counts.size(); ++i) {
        if (s->finiteCounts[i] < it->second.counts[i]) {
            appendError(s, QStringLiteral("%1: component %2 finite count %3 < flat count %4")
                               .arg(fieldStem(spec))
                               .arg(i)
                               .arg(s->finiteCounts[i])
                               .arg(it->second.counts[i]));
            ok = false;
        }
    }
    s->flatOk = ok;
}

bool scanField(const Body& body,
               const io::FieldSpec& spec,
               FieldDatasetStats* s) {
    QString providerReason;
    const FieldProvider provider = body.catalog.provider(body.run, spec.kind, &providerReason);
    s->provider = FieldProviderName(provider);
    s->nativeAxis = axisName(spec.axis);
    ++s->runs;

    if (provider == FieldProvider::DatasetAbsent) {
        s->frameExtentMin = std::min<std::size_t>(s->frameExtentMin, 0);
        s->frameExtentMax = std::max<std::size_t>(s->frameExtentMax, 0);
        s->nativeRowsMin = std::min<std::size_t>(s->nativeRowsMin, 0);
        s->nativeRowsMax = std::max<std::size_t>(s->nativeRowsMax, 0);
        if (spec.cols > 0) {
            s->componentsMax = std::max<std::size_t>(s->componentsMax,
                                                     static_cast<std::size_t>(spec.cols));
            if (s->finiteCounts.empty())
                s->finiteCounts.resize(static_cast<std::size_t>(spec.cols), 0);
        }
        QString reason;
        const bool p = body.catalog.present(body, spec.kind, 0, 0, 0, &reason);
        if (p) {
            appendError(s, QStringLiteral("%1 dataset-absent provider reported present")
                               .arg(fieldStem(spec)));
        } else {
            const QString recorded = reason.isEmpty() ? providerReason : reason;
            if (!recorded.isEmpty() && s->absenceReasons.size() < 8)
                s->absenceReasons.push_back(recorded);
        }
        return s->failures == 0;
    }

    if (provider == FieldProvider::Unsupported) {
        appendError(s, QStringLiteral("%1 unsupported in Catalog: %2")
                           .arg(fieldStem(spec), providerReason));
        return false;
    }

    if (structuredField(spec)) {
        const std::size_t rows = body.catalog.nativeRowCount(body, spec.kind, 0);
        s->frameExtentMin = std::min<std::size_t>(s->frameExtentMin, 1);
        s->frameExtentMax = std::max<std::size_t>(s->frameExtentMax, 1);
        ++s->framesVisited;
        s->nativeRowsMin = std::min(s->nativeRowsMin, rows);
        s->nativeRowsMax = std::max(s->nativeRowsMax, rows);
        s->nativeRowsTotal += rows;
        s->componentsMax = std::max<std::size_t>(s->componentsMax, 1);
        if (s->finiteCounts.empty()) s->finiteCounts.resize(1, 0);
        s->finiteCounts[0] += rows;
        s->valuesRead += rows;
        if (rows > 0) {
            QString reason;
            if (!body.catalog.present(body, spec.kind, 0, 0, -1, &reason))
                appendError(s, QStringLiteral("%1 structured row 0 absent: %2")
                                   .arg(fieldStem(spec), reason));
        }
        return s->failures == 0;
    }

    std::vector<std::size_t> frames;
    if (dftField(spec.kind)) {
        frames = body.run.frameMap.dftRows();
        if (body.run.frameMap.frameCount() > frames.size()) {
            for (std::size_t f = 0; f < body.run.frameMap.frameCount(); ++f) {
                if (std::find(frames.begin(), frames.end(), f) == frames.end()) {
                    QString reason;
                    const bool p = body.catalog.present(body, spec.kind, 0, f, 0, &reason);
                    if (p || reason != QStringLiteral("frame-gap"))
                        appendError(s, QStringLiteral("%1 DFT gap check failed at frame %2: present=%3 reason=%4")
                                           .arg(fieldStem(spec))
                                           .arg(f)
                                           .arg(p)
                                           .arg(reason));
                    break;
                }
            }
        }
    } else {
        const std::size_t n = body.catalog.frameCount(body, spec.kind);
        frames.reserve(n);
        for (std::size_t f = 0; f < n; ++f) frames.push_back(f);
    }

    s->frameExtentMin = std::min(s->frameExtentMin, frames.size());
    s->frameExtentMax = std::max(s->frameExtentMax, frames.size());
    for (std::size_t frame : frames) {
        ++s->framesVisited;
        const std::size_t rows = body.catalog.nativeRowCount(body, spec.kind, frame);
        const std::size_t comps = body.catalog.componentCount(body, spec.kind, frame);
        s->nativeRowsMin = std::min(s->nativeRowsMin, rows);
        s->nativeRowsMax = std::max(s->nativeRowsMax, rows);
        s->nativeRowsTotal += rows;
        s->componentsMax = std::max(s->componentsMax, comps);
        if (s->finiteCounts.size() < comps) s->finiteCounts.resize(comps, 0);
        for (std::size_t native = 0; native < rows; ++native) {
            for (std::size_t comp = 0; comp < comps; ++comp) {
                QString reason;
                const std::optional<double> v =
                    body.catalog.value(body, spec.kind, native, frame, static_cast<int>(comp), &reason);
                if (!v) {
                    appendError(s, QStringLiteral("%1 absent at frame=%2 native=%3 comp=%4: %5")
                                       .arg(fieldStem(spec))
                                       .arg(frame)
                                       .arg(native)
                                       .arg(comp)
                                       .arg(reason));
                    continue;
                }
                ++s->valuesRead;
                if (std::isfinite(*v)) ++s->finiteCounts[comp];
            }
        }
        if (rows > 0 && comps > 0) {
            QString reason;
            const bool p = body.catalog.present(body, spec.kind, 0, frame, static_cast<int>(comps), &reason);
            if (p || reason != QStringLiteral("component-out-of-range"))
                appendError(s, QStringLiteral("%1 component guard failed at frame=%2: present=%3 reason=%4")
                                   .arg(fieldStem(spec))
                                   .arg(frame)
                                   .arg(p)
                                   .arg(reason));
        }
    }
    if (frames.empty())
        appendError(s, QStringLiteral("%1 has zero readable frames").arg(fieldStem(spec)));
    return s->failures == 0;
}

bool scanRun(const RunData& run,
             std::map<int, FieldDatasetStats>* stats) {
    Catalog catalog(run);
    ResidentIndexes indexes = BuildResidentIndexes(run);
    const Body body{run, indexes, catalog};
    bool ok = true;
    for (const io::FieldSpec* spec : ScopedProducerCatalog()) {
        FieldDatasetStats& s = (*stats)[static_cast<int>(spec->kind)];
        if (!scanField(body, *spec, &s)) ok = false;
    }
    return ok;
}

bool writeBytes(const QString& path, const QByteArray& bytes, QString* err_out) {
    QSaveFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        if (err_out) *err_out = QStringLiteral("cannot write %1").arg(path);
        return false;
    }
    if (f.write(bytes) != bytes.size()) {
        if (err_out) *err_out = QStringLiteral("short write to %1").arg(path);
        return false;
    }
    if (!f.commit()) {
        if (err_out) *err_out = QStringLiteral("commit failed for %1").arg(path);
        return false;
    }
    return true;
}

bool writeReport(const QString& path,
                 const std::map<int, FieldDatasetStats>& staticStats,
                 const std::map<int, FieldDatasetStats>& trajStats,
                 bool pass,
                 QString* err_out) {
    QString text;
    QTextStream ts(&text);
    ts << "# SPINE wire-up report\n\n";
    ts << "Canonical inputs only. Probe result: " << (pass ? "PASS" : "FAIL") << "\n\n";
    ts << "## Repairs\n\n";
    ts << "- Wired `CanonicalSpineInputs.h` through the 720 and 1P9J load boundaries; non-canonical inputs hard-error before loading.\n";
    ts << "- Added loader-owned trajectory producer residency for scoped per-frame NPY fields, including variable native-row axes.\n";
    ts << "- Added FieldKind-native Catalog access with provider resolution, native row/frame/component reads, and absence reasons.\n";
    ts << "- Added full-detail reachability probe with flat coverage comparison.\n\n";
    ts << "- Made explicit `dataset_absent` a first-class provider: the probe verifies `present()` is false with a reason and still fails if flat coverage expects populated values.\n\n";
    ts << "- Removed the flat `RowDesignCatalogCoverage` generator, the row-design catalog-column bridge, and the three old row-design sidecar writers after the canonical probe passed.\n";
    ts << "- Deleted the canonical flat sidecar dump (`row_design_field_*` plus `catalog_sidecar_support.csv`) and the old `row_design_target_T2`, `row_design_ring_tensors`, and `row_design_aimnet2_embedding` files; retained only the small flat `catalog_coverage.json` floor snapshots for probe comparison.\n\n";
    ts << "## Per-field results\n\n";
    ts << "| field | native axis | 720 provider | 720 extent/result | 1P9J provider | 1P9J extent/result |\n";
    ts << "| --- | --- | --- | --- | --- | --- |\n";
    for (const io::FieldSpec* spec : ScopedProducerCatalog()) {
        const auto si = staticStats.find(static_cast<int>(spec->kind));
        const auto ti = trajStats.find(static_cast<int>(spec->kind));
        const FieldDatasetStats empty;
        const FieldDatasetStats& s = si == staticStats.end() ? empty : si->second;
        const FieldDatasetStats& t = ti == trajStats.end() ? empty : ti->second;
        const QString sResult = s.failures == 0 && s.flatOk && s.flatCompared ? QStringLiteral("PASS") : QStringLiteral("FAIL");
        const QString tResult = t.failures == 0 && t.flatOk && t.flatCompared ? QStringLiteral("PASS") : QStringLiteral("FAIL");
        ts << "| `" << fieldStem(*spec) << "` | `" << axisName(spec->axis) << "` | "
           << "`" << s.provider << "` | " << extentString(s) << "; flat=`"
           << s.flatRepresentation << "`; " << sResult << " | "
           << "`" << t.provider << "` | " << extentString(t) << "; flat=`"
           << t.flatRepresentation << "`; " << tResult << " |\n";
    }
    ts << "\n## Remaining\n\n";
    if (pass) {
        ts << "- None from the probe. The flat-path code cleanup, old sidecar-writer cleanup, and canonical flat dump cleanup are complete.\n";
    } else {
        ts << "- Probe failed; see `spine_reachability_manifest.json` for per-field errors. Flat path was not deleted.\n";
    }
    return writeBytes(path, text.toUtf8(), err_out);
}

}  // namespace

bool RunSpineReachabilityProbe(const SpineProbeConfig& cfg, QString* err_out) {
    if (!ValidateCanonical720Root(cfg.root720, err_out)) return false;
    if (!ValidateCanonical1p9jRoot(cfg.run1p9j, err_out)) return false;
    QDir().mkpath(cfg.outDir);

    CoverageMap flat720;
    CoverageMap flat1p9j;
    if (!loadCoverage(cfg.flat720Coverage, &flat720, err_out)) return false;
    if (!loadCoverage(cfg.flat1p9jCoverage, &flat1p9j, err_out)) return false;

    std::map<int, FieldDatasetStats> staticStats;
    std::map<int, FieldDatasetStats> trajStats;

    bool ok = true;
    QDir root(cfg.root720);
    const QStringList dirs = root.entryList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
    if (dirs.isEmpty()) {
        if (err_out) *err_out = QStringLiteral("canonical 720 root has no pose dirs: %1").arg(cfg.root720);
        return false;
    }
    for (const QString& name : dirs) {
        QString err;
        auto run = StaticRunData::Load(root.filePath(name), &err);
        if (!run) {
            if (err_out) *err_out = QStringLiteral("static probe load failed for %1: %2").arg(name, err);
            return false;
        }
        if (!scanRun(*run, &staticStats)) ok = false;
    }

    QString runErr;
    auto traj = RunLoader::Load(cfg.run1p9j, &runErr);
    if (!traj) {
        if (err_out) *err_out = QStringLiteral("1P9J probe load failed: %1").arg(runErr);
        return false;
    }
    if (!scanRun(*traj, &trajStats)) ok = false;

    for (const io::FieldSpec* spec : ScopedProducerCatalog()) {
        compareFlat(*spec, flat720, &staticStats[static_cast<int>(spec->kind)]);
        compareFlat(*spec, flat1p9j, &trajStats[static_cast<int>(spec->kind)]);
        const FieldDatasetStats& ss = staticStats[static_cast<int>(spec->kind)];
        const FieldDatasetStats& ts = trajStats[static_cast<int>(spec->kind)];
        if (ss.failures != 0 || ts.failures != 0 || !ss.flatCompared || !ts.flatCompared
            || !ss.flatOk || !ts.flatOk) {
            ok = false;
        }
    }

    QJsonObject doc;
    doc.insert(QStringLiteral("schema_version"), 1);
    doc.insert(QStringLiteral("pass"), ok);
    doc.insert(QStringLiteral("canonical_720_root"), QFileInfo(cfg.root720).absoluteFilePath());
    doc.insert(QStringLiteral("canonical_1p9j_root"), QFileInfo(cfg.run1p9j).absoluteFilePath());
    doc.insert(QStringLiteral("flat_720_coverage"), QFileInfo(cfg.flat720Coverage).absoluteFilePath());
    doc.insert(QStringLiteral("flat_1p9j_coverage"), QFileInfo(cfg.flat1p9jCoverage).absoluteFilePath());
    doc.insert(QStringLiteral("field_count"), static_cast<int>(ScopedProducerCatalog().size()));
    QJsonArray fields;
    for (const io::FieldSpec* spec : ScopedProducerCatalog()) {
        QJsonObject f;
        f.insert(QStringLiteral("stem"), fieldStem(*spec));
        f.insert(QStringLiteral("native_axis"), axisName(spec->axis));
        f.insert(QStringLiteral("catalog_cols"), spec->cols);
        f.insert(QStringLiteral("static_720"), statsJson(staticStats[static_cast<int>(spec->kind)]));
        f.insert(QStringLiteral("trajectory_1p9j"), statsJson(trajStats[static_cast<int>(spec->kind)]));
        fields.push_back(f);
    }
    doc.insert(QStringLiteral("fields"), fields);
    const QString manifestPath = QDir(cfg.outDir).filePath(QStringLiteral("spine_reachability_manifest.json"));
    if (!writeBytes(manifestPath, QJsonDocument(doc).toJson(QJsonDocument::Indented), err_out))
        return false;
    if (!cfg.reportPath.isEmpty()
        && !writeReport(cfg.reportPath, staticStats, trajStats, ok, err_out))
        return false;
    if (!ok && err_out)
        *err_out = QStringLiteral("spine reachability probe failed; see %1").arg(manifestPath);
    return ok;
}

}  // namespace h5reader::rediscover

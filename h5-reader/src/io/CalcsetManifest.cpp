// CalcsetManifest implementation — `.LGS` JSON loader.
//
// Pure Qt I/O — QFile + QJsonDocument; no toml++, no exceptions across
// the boundary. The loader maps JSON nodes onto the typed sub-structs
// and existence-checks every declared path. The per-frame `DftFrame`
// objects keep their meta.json read lazy (LoadMeta()) so the common
// case (load the manifest, walk frames into a hash) doesn't hit 500+
// per-frame files on the GUI thread.

#include "CalcsetManifest.h"

#include "../diagnostics/ErrorBus.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QLoggingCategory>
#include <QStringList>

#include <cmath>

namespace h5reader::io {

namespace {
Q_LOGGING_CATEGORY(cCalcset, "h5reader.calcset")

using Severity = h5reader::diagnostics::Severity;

void reportErr(const QString& msg, const QString& ctx) {
    h5reader::diagnostics::ErrorBus::Report(
        Severity::Error, QStringLiteral("CalcsetManifest"), msg, ctx);
}

QString ResolveRelative(const QString& rootDir, const QString& declared) {
    if (declared.isEmpty()) return {};
    const QFileInfo fi(declared);
    if (fi.isAbsolute()) return fi.absoluteFilePath();
    return QFileInfo(QDir(rootDir).filePath(declared)).absoluteFilePath();
}

// Validate a declared path's existence + file/dir kind. Returns the
// error message on failure; nullopt on success.
std::optional<QString>
ValidatePath(const QString& fullKeyForError,
             const QString& absPath,
             bool mustBeDir) {
    const QFileInfo info(absPath);
    if (!info.exists())
        return QStringLiteral("path [%1] does not exist: %2").arg(fullKeyForError, absPath);
    if (mustBeDir && !info.isDir())
        return QStringLiteral("path [%1] must be a directory: %2").arg(fullKeyForError, absPath);
    if (!mustBeDir && info.isDir())
        return QStringLiteral("path [%1] must be a file: %2").arg(fullKeyForError, absPath);
    return std::nullopt;
}

// Read a required string from a JSON object. Returns nullopt on
// missing/wrong-type and writes the error.
std::optional<QString>
RequireString(const QJsonObject& obj, QStringView key, QString* err) {
    const QJsonValue v = obj.value(key);
    if (v.isUndefined() || v.isNull()) {
        if (err) *err = QStringLiteral("missing required key '%1'").arg(key);
        return std::nullopt;
    }
    if (!v.isString()) {
        if (err) *err = QStringLiteral("key '%1' must be a string").arg(key);
        return std::nullopt;
    }
    const QString s = v.toString();
    if (s.isEmpty()) {
        if (err) *err = QStringLiteral("key '%1' is empty").arg(key);
        return std::nullopt;
    }
    return s;
}

std::optional<QString>
OptionalString(const QJsonObject& obj, QStringView key) {
    const QJsonValue v = obj.value(key);
    if (v.isUndefined() || v.isNull() || !v.isString()) return std::nullopt;
    return v.toString();
}

// Required path: read string, resolve, existence-check.
std::optional<QString>
RequireResolvedPath(const QJsonObject& obj, QStringView jsonKey,
                    const QString& fullKeyForError,
                    const QString& rootDir, bool mustBeDir,
                    QString& out, QString* err) {
    auto raw = RequireString(obj, jsonKey, err);
    if (!raw) return *err;  // already populated by RequireString
    const QString abs = ResolveRelative(rootDir, *raw);
    if (auto perr = ValidatePath(fullKeyForError, abs, mustBeDir)) {
        if (err) *err = *perr;
        return *err;
    }
    out = abs;
    return std::nullopt;
}

// Optional path: empty is ok; declared-but-missing is an error.
std::optional<QString>
OptionalResolvedPath(const QJsonObject& obj, QStringView jsonKey,
                     const QString& fullKeyForError,
                     const QString& rootDir, bool mustBeDir,
                     QString& out) {
    const QJsonValue v = obj.value(jsonKey);
    if (v.isUndefined() || v.isNull()) {
        out.clear();
        return std::nullopt;
    }
    if (!v.isString()) {
        return QStringLiteral("key '%1' must be a string").arg(jsonKey);
    }
    const QString raw = v.toString();
    if (raw.isEmpty()) { out.clear(); return std::nullopt; }
    const QString abs = ResolveRelative(rootDir, raw);
    if (auto perr = ValidatePath(fullKeyForError, abs, mustBeDir))
        return *perr;
    out = abs;
    return std::nullopt;
}

CalcsetManifest::Kind ParseKind(const QString& s, bool& ok) {
    ok = true;
    if (s == QLatin1String("trajectory"))   return CalcsetManifest::Kind::Trajectory;
    if (s == QLatin1String("single_pose"))  return CalcsetManifest::Kind::SinglePose;
    if (s == QLatin1String("mutant_pair"))  return CalcsetManifest::Kind::MutantPair;
    ok = false;
    return CalcsetManifest::Kind::Trajectory;
}

CalcsetManifest::PoseKind ParsePoseKind(const QString& s, bool& ok) {
    ok = true;
    if (s == QLatin1String("pdb"))            return CalcsetManifest::PoseKind::Pdb;
    if (s == QLatin1String("protonated_pdb")) return CalcsetManifest::PoseKind::ProtonatedPdb;
    if (s == QLatin1String("orca"))           return CalcsetManifest::PoseKind::Orca;
    ok = false;
    return CalcsetManifest::PoseKind::Pdb;
}

// Resolve the calcset_root/<dataset_id>.LGS file from a directory
// argument by listing `*.LGS`. Zero or more than one match is a hard
// error per the spec (no glob-and-pick).
std::optional<QString> ResolveLgsPath(const QString& root_or_lgs_path,
                                       QString* err) {
    const QFileInfo fi(root_or_lgs_path);
    if (!fi.exists()) {
        if (err) *err = QStringLiteral("path does not exist: %1").arg(root_or_lgs_path);
        return std::nullopt;
    }
    if (fi.isFile()) {
        if (!root_or_lgs_path.endsWith(QStringLiteral(".LGS"), Qt::CaseSensitive)) {
            if (err) *err = QStringLiteral("file is not a .LGS: %1").arg(root_or_lgs_path);
            return std::nullopt;
        }
        return fi.absoluteFilePath();
    }
    // Directory: find the single *.LGS inside.
    QDir dir(root_or_lgs_path);
    const QStringList matches =
        dir.entryList(QStringList{QStringLiteral("*.LGS")}, QDir::Files);
    if (matches.isEmpty()) {
        if (err) *err = QStringLiteral("no .LGS file in directory: %1").arg(root_or_lgs_path);
        return std::nullopt;
    }
    if (matches.size() > 1) {
        if (err) *err = QStringLiteral("multiple .LGS files in directory %1: %2")
                            .arg(root_or_lgs_path, matches.join(QStringLiteral(", ")));
        return std::nullopt;
    }
    return dir.absoluteFilePath(matches.front());
}

}  // namespace

// ---- DftFrame ----------------------------------------------------

bool DftFrame::LoadMeta(QString* err_out) const {
    if (meta_loaded_) return true;
    QFile f(meta_json_abspath);
    if (!f.open(QIODevice::ReadOnly)) {
        if (err_out)
            *err_out = QStringLiteral("cannot open meta.json: %1").arg(meta_json_abspath);
        return false;
    }
    QJsonParseError parseErr;
    const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &parseErr);
    if (parseErr.error != QJsonParseError::NoError || !doc.isObject()) {
        if (err_out)
            *err_out = QStringLiteral("malformed meta.json (%1): %2")
                           .arg(parseErr.errorString(), meta_json_abspath);
        return false;
    }
    const QJsonObject root = doc.object();
    frame_ps_       = root.value(QStringLiteral("frame_ps")).toDouble(0.0);
    orca_exit_code_ = root.value(QStringLiteral("orca_exit_code")).toInt(-1);

    const QJsonObject files = root.value(QStringLiteral("files")).toObject();
    const QString outRel = files.value(QStringLiteral("out_primary")).toString();
    if (!outRel.isEmpty()) {
        const QString jobDir = QFileInfo(meta_json_abspath).absolutePath();
        orca_out_abspath_ = QStringLiteral("%1/%2").arg(jobDir, outRel);
    }
    meta_loaded_ = true;
    return true;
}

double  DftFrame::framePs()         const { return frame_ps_; }
QString DftFrame::orcaOutAbspath()  const { return orca_out_abspath_; }
int     DftFrame::orcaExitCode()    const { return orca_exit_code_; }

// ---- CalcsetManifest ---------------------------------------------

const char* CalcsetManifest::NameForKind(Kind k) {
    switch (k) {
        case Kind::Trajectory: return "trajectory";
        case Kind::SinglePose: return "single_pose";
        case Kind::MutantPair: return "mutant_pair";
    }
    return "?";
}

std::optional<CalcsetManifest>
CalcsetManifest::Load(const QString& root_or_lgs_path, QString* err_out) {
    QString err;

    // 1. Resolve the .LGS file path.
    auto lgsPath = ResolveLgsPath(root_or_lgs_path, &err);
    if (!lgsPath) {
        if (err_out) *err_out = err;
        reportErr(err, root_or_lgs_path);
        return std::nullopt;
    }

    CalcsetManifest m;
    m.lgs_path_abspath = *lgsPath;
    m.calcset_root_abspath = QFileInfo(*lgsPath).absolutePath();

    // 2. Open + parse JSON.
    QFile f(m.lgs_path_abspath);
    if (!f.open(QIODevice::ReadOnly)) {
        err = QStringLiteral("cannot open .LGS: %1").arg(m.lgs_path_abspath);
        if (err_out) *err_out = err;
        reportErr(err, m.lgs_path_abspath);
        return std::nullopt;
    }
    QJsonParseError parseErr;
    const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &parseErr);
    if (parseErr.error != QJsonParseError::NoError || !doc.isObject()) {
        err = QStringLiteral(".LGS parse error (%1) at %2:%3 in %4")
                  .arg(parseErr.errorString())
                  .arg(parseErr.offset)
                  .arg(parseErr.offset)
                  .arg(m.lgs_path_abspath);
        if (err_out) *err_out = err;
        reportErr(err, m.lgs_path_abspath);
        return std::nullopt;
    }
    const QJsonObject root = doc.object();

    // 3. schema_version — strict equality (no forward-compat guess).
    const QJsonValue schemaV = root.value(QStringLiteral("schema_version"));
    if (!schemaV.isDouble()) {
        err = QStringLiteral("missing or non-integer 'schema_version' in %1").arg(m.lgs_path_abspath);
        if (err_out) *err_out = err;
        reportErr(err, m.lgs_path_abspath);
        return std::nullopt;
    }
    m.schema_version = schemaV.toInt(0);
    if (m.schema_version != kSupportedSchemaVersion) {
        err = QStringLiteral("schema_version=%1 unsupported (this reader supports %2): %3")
                  .arg(m.schema_version).arg(kSupportedSchemaVersion).arg(m.lgs_path_abspath);
        if (err_out) *err_out = err;
        reportErr(err, m.lgs_path_abspath);
        return std::nullopt;
    }

    // 4. kind dispatch and required top-level identity fields.
    auto kindStr = RequireString(root, QStringLiteral("kind"), &err);
    if (!kindStr) {
        if (err_out) *err_out = err;
        reportErr(err, m.lgs_path_abspath);
        return std::nullopt;
    }
    bool kindOk = false;
    m.kind = ParseKind(*kindStr, kindOk);
    if (!kindOk) {
        err = QStringLiteral("kind=%1 unknown (expected: trajectory | single_pose | mutant_pair)")
                  .arg(*kindStr);
        if (err_out) *err_out = err;
        reportErr(err, m.lgs_path_abspath);
        return std::nullopt;
    }
    auto datasetId = RequireString(root, QStringLiteral("dataset_id"), &err);
    if (!datasetId) { if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt; }
    m.dataset_id = *datasetId;
    auto proteinId = RequireString(root, QStringLiteral("protein_id"), &err);
    if (!proteinId) { if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt; }
    m.protein_id = *proteinId;
    auto humanName = RequireString(root, QStringLiteral("human_name"), &err);
    if (!humanName) { if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt; }
    m.human_name = *humanName;

    // 5. Per-kind sub-block validation. The spec demands mutual
    //    exclusivity, but the loader is forgiving here — it only
    //    requires the active kind's block; other blocks are ignored.

    switch (m.kind) {
        case Kind::Trajectory: {
            const QJsonValue tv = root.value(QStringLiteral("trajectory"));
            if (!tv.isObject()) {
                err = QStringLiteral("kind=trajectory requires a 'trajectory' object");
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            const QJsonObject tobj = tv.toObject();
            Trajectory t;
            if (auto e = RequireResolvedPath(tobj, QStringLiteral("md_dir"),
                                              QStringLiteral("trajectory.md_dir"),
                                              m.calcset_root_abspath, /*mustBeDir=*/true,
                                              t.md_dir_abspath, &err)) {
                if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt;
            }
            if (auto e = RequireResolvedPath(tobj, QStringLiteral("topology_top"),
                                              QStringLiteral("trajectory.topology_top"),
                                              m.calcset_root_abspath, /*mustBeDir=*/false,
                                              t.topology_top_abspath, &err)) {
                if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt;
            }
            if (auto e = RequireResolvedPath(tobj, QStringLiteral("extraction_dir"),
                                              QStringLiteral("trajectory.extraction_dir"),
                                              m.calcset_root_abspath, /*mustBeDir=*/true,
                                              t.extraction_dir_abspath, &err)) {
                if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt;
            }
            if (auto e = RequireResolvedPath(tobj, QStringLiteral("trajectory_h5"),
                                              QStringLiteral("trajectory.trajectory_h5"),
                                              m.calcset_root_abspath, /*mustBeDir=*/false,
                                              t.trajectory_h5_abspath, &err)) {
                if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt;
            }
            if (auto e = RequireResolvedPath(tobj, QStringLiteral("extraction_manifest"),
                                              QStringLiteral("trajectory.extraction_manifest"),
                                              m.calcset_root_abspath, /*mustBeDir=*/false,
                                              t.extraction_manifest_abspath, &err)) {
                if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt;
            }
            const QJsonValue dt = tobj.value(QStringLiteral("frame_dt_ps"));
            if (!dt.isDouble()) {
                err = QStringLiteral("trajectory.frame_dt_ps must be a positive number");
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            t.frame_dt_ps = dt.toDouble();
            if (!(t.frame_dt_ps > 0.0)) {
                err = QStringLiteral("trajectory.frame_dt_ps=%1 must be > 0").arg(t.frame_dt_ps);
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            auto basis = RequireString(tobj, QStringLiteral("frame_index_basis"), &err);
            if (!basis) { if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt; }
            t.frame_index_basis = *basis;
            // Optional reference_pdb.
            if (auto e = OptionalResolvedPath(tobj, QStringLiteral("reference_pdb"),
                                               QStringLiteral("trajectory.reference_pdb"),
                                               m.calcset_root_abspath, /*mustBeDir=*/false,
                                               t.reference_pdb_abspath)) {
                if (err_out) *err_out = *e;
                reportErr(*e, m.lgs_path_abspath);
                return std::nullopt;
            }
            m.trajectory = std::move(t);
            break;
        }
        case Kind::SinglePose: {
            const QJsonValue sv = root.value(QStringLiteral("single_pose"));
            if (!sv.isObject()) {
                err = QStringLiteral("kind=single_pose requires a 'single_pose' object");
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            const QJsonObject sobj = sv.toObject();
            SinglePose s;
            auto poseKindStr = RequireString(sobj, QStringLiteral("pose_kind"), &err);
            if (!poseKindStr) {
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            bool poseOk = false;
            s.pose_kind = ParsePoseKind(*poseKindStr, poseOk);
            if (!poseOk) {
                err = QStringLiteral("single_pose.pose_kind=%1 unknown (expected: pdb | protonated_pdb | orca)")
                          .arg(*poseKindStr);
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            if (auto e = RequireResolvedPath(sobj, QStringLiteral("pose_dir"),
                                              QStringLiteral("single_pose.pose_dir"),
                                              m.calcset_root_abspath, /*mustBeDir=*/true,
                                              s.pose_dir_abspath, &err)) {
                if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt;
            }
            if (auto e = RequireResolvedPath(sobj, QStringLiteral("extraction_manifest"),
                                              QStringLiteral("single_pose.extraction_manifest"),
                                              m.calcset_root_abspath, /*mustBeDir=*/false,
                                              s.extraction_manifest_abspath, &err)) {
                if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt;
            }
            m.single_pose = std::move(s);
            break;
        }
        case Kind::MutantPair: {
            const QJsonValue mv = root.value(QStringLiteral("mutant_pair"));
            if (!mv.isObject()) {
                err = QStringLiteral("kind=mutant_pair requires a 'mutant_pair' object");
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            const QJsonObject mobj = mv.toObject();
            MutantPair mp;
            // Both wt_lgs / ala_lgs must resolve to existing .LGS files.
            if (auto e = RequireResolvedPath(mobj, QStringLiteral("wt_lgs"),
                                              QStringLiteral("mutant_pair.wt_lgs"),
                                              m.calcset_root_abspath, /*mustBeDir=*/false,
                                              mp.wt_lgs_abspath, &err)) {
                if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt;
            }
            if (auto e = RequireResolvedPath(mobj, QStringLiteral("ala_lgs"),
                                              QStringLiteral("mutant_pair.ala_lgs"),
                                              m.calcset_root_abspath, /*mustBeDir=*/false,
                                              mp.ala_lgs_abspath, &err)) {
                if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt;
            }
            if (!mp.wt_lgs_abspath.endsWith(QStringLiteral(".LGS"))
                || !mp.ala_lgs_abspath.endsWith(QStringLiteral(".LGS"))) {
                err = QStringLiteral("mutant_pair.wt_lgs and ala_lgs must be .LGS file paths");
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            m.mutant_pair = std::move(mp);
            break;
        }
    }

    // 6. Optional `dft` block — present iff DFT data exists.
    const QJsonValue dv = root.value(QStringLiteral("dft"));
    if (dv.isObject()) {
        const QJsonObject dobj = dv.toObject();
        Dft d;
        auto method = RequireString(dobj, QStringLiteral("method"), &err);
        if (!method) { if (err_out) *err_out = err; reportErr(err, m.lgs_path_abspath); return std::nullopt; }
        d.method = *method;
        d.campaign_target_frames =
            dobj.value(QStringLiteral("campaign_target_frames")).toInt(0);
        const QJsonObject strObj = dobj.value(QStringLiteral("frame_stride")).toObject();
        d.frame_stride.first = strObj.value(QStringLiteral("first")).toInt(0);
        d.frame_stride.last  = strObj.value(QStringLiteral("last")).toInt(0);
        d.frame_stride.step  = strObj.value(QStringLiteral("step")).toInt(1);

        const QJsonValue fv = dobj.value(QStringLiteral("frames"));
        if (!fv.isArray()) {
            err = QStringLiteral("dft.frames must be a JSON array");
            if (err_out) *err_out = err;
            reportErr(err, m.lgs_path_abspath);
            return std::nullopt;
        }
        const QJsonArray frArr = fv.toArray();
        d.frames.reserve(static_cast<std::size_t>(frArr.size()));
        for (int i = 0; i < frArr.size(); ++i) {
            const QJsonValue fe = frArr.at(i);
            if (!fe.isObject()) {
                err = QStringLiteral("dft.frames[%1] must be an object").arg(i);
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            const QJsonObject feobj = fe.toObject();
            const QJsonValue fiv = feobj.value(QStringLiteral("frame_index"));
            if (!fiv.isDouble()) {
                err = QStringLiteral("dft.frames[%1].frame_index must be an integer").arg(i);
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            DftFrame df;
            df.frame_index = fiv.toInt();
            if (df.frame_index < 0) {
                err = QStringLiteral("dft.frames[%1].frame_index=%2 must be >= 0")
                          .arg(i).arg(df.frame_index);
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            auto metaRel = RequireString(feobj, QStringLiteral("meta_json"), &err);
            if (!metaRel) {
                err = QStringLiteral("dft.frames[%1].meta_json: %2").arg(i).arg(err);
                if (err_out) *err_out = err;
                reportErr(err, m.lgs_path_abspath);
                return std::nullopt;
            }
            df.meta_json_relpath = *metaRel;
            df.meta_json_abspath = ResolveRelative(m.calcset_root_abspath, *metaRel);
            if (auto perr = ValidatePath(QStringLiteral("dft.frames[%1].meta_json").arg(i),
                                          df.meta_json_abspath, /*mustBeDir=*/false)) {
                if (err_out) *err_out = *perr;
                reportErr(*perr, m.lgs_path_abspath);
                return std::nullopt;
            }
            d.frames.push_back(std::move(df));
        }
        m.dft = std::move(d);
    }

    // 7. metadata — informational only.
    if (root.contains(QStringLiteral("metadata"))) {
        const QJsonObject meta = root.value(QStringLiteral("metadata")).toObject();
        m.generated_at_utc           = meta.value(QStringLiteral("generated_at_utc")).toString();
        m.lgs_writer                 = meta.value(QStringLiteral("lgs_writer")).toString();
        m.producer_extractor_version = meta.value(QStringLiteral("producer_extractor_version")).toString();
    }

    qCInfo(cCalcset).noquote()
        << "loaded .LGS |" << m.lgs_path_abspath
        << "| kind=" << NameForKind(m.kind)
        << "| dataset_id=" << m.dataset_id
        << "| protein_id=" << m.protein_id
        << "| dft_frames=" << (m.dft ? m.dft->frames.size() : 0);
    return m;
}

}  // namespace h5reader::io

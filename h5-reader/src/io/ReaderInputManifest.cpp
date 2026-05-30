// ReaderInputManifest implementation — toml++ parser.
//
// toml++ does the TOML lexing/typing; this file is the mapping layer
// from typed table reads into our struct, plus path resolution and
// existence validation.

#include "ReaderInputManifest.h"

#include "../diagnostics/ErrorBus.h"

#include <QDir>
#include <QFileInfo>

#include <toml.hpp>

#include <cstdint>
#include <optional>
#include <sstream>
#include <string_view>

namespace h5reader::io {

namespace {

QString ResolveRelative(const QString& rootDir, const QString& declared) {
    if (declared.isEmpty()) return {};
    const QFileInfo fi(declared);
    if (fi.isAbsolute()) return fi.absoluteFilePath();
    return QFileInfo(QDir(rootDir).filePath(declared)).absoluteFilePath();
}

std::optional<QString> ReadString(const toml::table& t, std::string_view key) {
    if (auto v = t[key].value<std::string>()) {
        return QString::fromStdString(*v);
    }
    return std::nullopt;
}

std::optional<int> ReadInt(const toml::table& t, std::string_view key) {
    if (auto v = t[key].value<std::int64_t>()) {
        return static_cast<int>(*v);
    }
    return std::nullopt;
}

std::optional<QString> ValidatePath(const QString& key,
                                     const QString& absPath,
                                     bool mustBeDir) {
    const QFileInfo info(absPath);
    if (!info.exists())
        return QStringLiteral("path [%1] does not exist: %2").arg(key, absPath);
    if (mustBeDir && !info.isDir())
        return QStringLiteral("path [%1] must be a directory: %2").arg(key, absPath);
    if (!mustBeDir && info.isDir())
        return QStringLiteral("path [%1] must be a file: %2").arg(key, absPath);
    return std::nullopt;
}

// Read a path key from a toml::table, resolve, validate. Returns the
// error message on failure, nullopt on success (with `out` populated).
std::optional<QString> ReadValidatedPath(const toml::table& t,
                                          std::string_view key,
                                          const QString& fullKeyForError,
                                          const QString& rootDir,
                                          bool required,
                                          bool mustBeDir,
                                          QString& out) {
    auto raw = ReadString(t, key);
    if (!raw || raw->isEmpty()) {
        if (required)
            return QStringLiteral("missing required key [%1]").arg(fullKeyForError);
        out.clear();
        return std::nullopt;
    }
    const QString abs = ResolveRelative(rootDir, *raw);
    if (auto err = ValidatePath(fullKeyForError, abs, mustBeDir))
        return err;
    out = abs;
    return std::nullopt;
}

ReaderInputManifest::RunKind ParseRunKind(const QString& s, bool& ok) {
    ok = true;
    if (s == QLatin1String("trajectory"))  return ReaderInputManifest::RunKind::Trajectory;
    if (s == QLatin1String("single_pose")) return ReaderInputManifest::RunKind::SinglePose;
    if (s == QLatin1String("mutant_pair")) return ReaderInputManifest::RunKind::MutantPair;
    ok = false;
    return ReaderInputManifest::RunKind::Trajectory;
}

}  // namespace

const char* ReaderInputManifest::NameForRunKind(RunKind k) {
    switch (k) {
        case RunKind::Trajectory: return "trajectory";
        case RunKind::SinglePose: return "single_pose";
        case RunKind::MutantPair: return "mutant_pair";
    }
    return "?";
}

bool ReaderInputManifest::ExistsAt(const QString& dir) {
    return QFileInfo::exists(QDir(dir).filePath(QString::fromLatin1(kFileName)));
}

QString ReaderInputManifest::primaryPoseDir() const {
    return runKind == RunKind::MutantPair ? wtPoseDir : QString();
}

QString ReaderInputManifest::alternatePoseDir() const {
    return runKind == RunKind::MutantPair ? alaPoseDir : QString();
}

ReaderInputManifest ReaderInputManifest::Load(const QString& dir) {
    ReaderInputManifest m;
    m.rootDir = QDir(dir).absolutePath();
    m.manifestPath = QDir(m.rootDir).filePath(QString::fromLatin1(kFileName));

    if (!QFileInfo::exists(m.manifestPath)) {
        m.error = QStringLiteral("ReaderInputManifest: not found at %1").arg(m.manifestPath);
        // No ErrorBus report — absence is the caller's call to make
        // (existing fixtures legitimately have no manifest yet).
        return m;
    }

    toml::table root;
    try {
        root = toml::parse_file(m.manifestPath.toStdString());
    } catch (const toml::parse_error& e) {
        const std::string desc(e.description());
        m.error = QStringLiteral("ReaderInputManifest: TOML parse failed at %1: %2 (line %3, col %4)")
                      .arg(m.manifestPath, QString::fromStdString(desc))
                      .arg(e.source().begin.line)
                      .arg(e.source().begin.column);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                 QStringLiteral("ReaderInputManifest"),
                                                 m.error,
                                                 m.manifestPath);
        return m;
    } catch (const std::exception& e) {
        m.error = QStringLiteral("ReaderInputManifest: parse exception at %1: %2")
                      .arg(m.manifestPath, QString::fromUtf8(e.what()));
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                 QStringLiteral("ReaderInputManifest"),
                                                 m.error,
                                                 m.manifestPath);
        return m;
    }

    auto reportErr = [&](const QString& detail) {
        m.error = QStringLiteral("ReaderInputManifest at %1: %2").arg(m.manifestPath, detail);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                 QStringLiteral("ReaderInputManifest"),
                                                 m.error,
                                                 m.manifestPath);
    };

    // ----- top-level -----
    auto schemaVer = ReadInt(root, "schema_version");
    if (!schemaVer) {
        reportErr(QStringLiteral("missing required key schema_version"));
        return m;
    }
    m.schemaVersion = *schemaVer;
    if (m.schemaVersion != kSupportedSchemaVersion) {
        reportErr(QStringLiteral("schema_version=%1 unsupported (this reader supports %2)")
                      .arg(m.schemaVersion).arg(kSupportedSchemaVersion));
        return m;
    }

    auto runKindStr = ReadString(root, "run_kind");
    if (!runKindStr) {
        reportErr(QStringLiteral("missing required key run_kind"));
        return m;
    }
    bool kindOk = false;
    m.runKind = ParseRunKind(*runKindStr, kindOk);
    if (!kindOk) {
        reportErr(QStringLiteral("run_kind=%1 unknown (expected: trajectory | single_pose | mutant_pair)")
                      .arg(*runKindStr));
        return m;
    }

    auto proteinId = ReadString(root, "protein_id");
    if (!proteinId || proteinId->isEmpty()) {
        reportErr(QStringLiteral("missing required key protein_id"));
        return m;
    }
    m.proteinId = *proteinId;
    m.humanName = ReadString(root, "human_name").value_or(QString());

    // ----- per-run-kind table -----
    switch (m.runKind) {
        case RunKind::Trajectory: {
            auto* t = root["trajectory"].as_table();
            if (!t) {
                reportErr(QStringLiteral("run_kind=trajectory requires a [trajectory] table"));
                return m;
            }
            if (auto err = ReadValidatedPath(*t, "trajectory_h5",
                                              QStringLiteral("trajectory.trajectory_h5"),
                                              m.rootDir, true, false, m.trajectoryH5)) {
                reportErr(*err); return m;
            }
            if (auto err = ReadValidatedPath(*t, "topology_sidecar_dir",
                                              QStringLiteral("trajectory.topology_sidecar_dir"),
                                              m.rootDir, true, true, m.topologySidecarDir)) {
                reportErr(*err); return m;
            }
            if (auto err = ReadValidatedPath(*t, "per_frame_npys_dir",
                                              QStringLiteral("trajectory.per_frame_npys_dir"),
                                              m.rootDir, false, true, m.perFrameNpysDir)) {
                reportErr(*err); return m;
            }
            break;
        }
        case RunKind::SinglePose: {
            auto* t = root["single_pose"].as_table();
            if (!t) {
                reportErr(QStringLiteral("run_kind=single_pose requires a [single_pose] table"));
                return m;
            }
            if (auto err = ReadValidatedPath(*t, "pose_dir",
                                              QStringLiteral("single_pose.pose_dir"),
                                              m.rootDir, true, true, m.poseDir)) {
                reportErr(*err); return m;
            }
            break;
        }
        case RunKind::MutantPair: {
            auto* t = root["mutant_pair"].as_table();
            if (!t) {
                reportErr(QStringLiteral("run_kind=mutant_pair requires a [mutant_pair] table"));
                return m;
            }
            if (auto err = ReadValidatedPath(*t, "wt_pose_dir",
                                              QStringLiteral("mutant_pair.wt_pose_dir"),
                                              m.rootDir, true, true, m.wtPoseDir)) {
                reportErr(*err); return m;
            }
            if (auto err = ReadValidatedPath(*t, "ala_pose_dir",
                                              QStringLiteral("mutant_pair.ala_pose_dir"),
                                              m.rootDir, true, true, m.alaPoseDir)) {
                reportErr(*err); return m;
            }
            break;
        }
    }

    // ----- optional [dft] -----
    if (auto* d = root["dft"].as_table()) {
        if (auto err = ReadValidatedPath(*d, "jobs_dir",
                                          QStringLiteral("dft.jobs_dir"),
                                          m.rootDir, true, true, m.dftJobsDir)) {
            reportErr(*err); return m;
        }
    }

    // ----- optional top-level reference_pdb -----
    if (auto raw = ReadString(root, "reference_pdb"); raw && !raw->isEmpty()) {
        const QString abs = ResolveRelative(m.rootDir, *raw);
        if (auto err = ValidatePath(QStringLiteral("reference_pdb"), abs, false)) {
            reportErr(*err); return m;
        }
        m.referencePdb = abs;
    }

    // ----- optional [metadata] (informational) -----
    if (auto* meta = root["metadata"].as_table()) {
        // extracted_at: TOML datetime or string.
        if (auto dt = (*meta)["extracted_at"].value<toml::date_time>()) {
            std::ostringstream oss;
            oss << *dt;
            m.extractedAt = QString::fromStdString(oss.str());
        } else {
            m.extractedAt = ReadString(*meta, "extracted_at").value_or(QString());
        }
        m.extractVersion = ReadString(*meta, "extract_version").value_or(QString());
        m.fixtureKind    = ReadString(*meta, "fixture_kind").value_or(QString());
    }

    m.ok = true;
    return m;
}

}  // namespace h5reader::io

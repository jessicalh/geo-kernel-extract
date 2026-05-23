// QtManifest implementation — QJsonDocument-based parser for
// extraction_manifest.json. Cross-platform; no POSIX.
//
// Failure modes log via ErrorBus and return a QtManifest with
// ok=false + error populated. Partial parse (missing optional fields)
// logs Warn and continues with defaults; missing REQUIRED fields
// (schema_version, protein_id, axis_sizes.atom) fail the load.

#include "QtManifest.h"

#include "../diagnostics/ErrorBus.h"

#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>

namespace h5reader::io {

namespace {

bool ReadString(const QJsonObject& obj, const QString& key, QString& out) {
    if (!obj.contains(key))
        return false;
    const QJsonValue v = obj.value(key);
    if (!v.isString())
        return false;
    out = v.toString();
    return true;
}

bool ReadBool(const QJsonObject& obj, const QString& key, bool& out) {
    if (!obj.contains(key))
        return false;
    const QJsonValue v = obj.value(key);
    if (!v.isBool())
        return false;
    out = v.toBool();
    return true;
}

bool ReadSizeT(const QJsonObject& obj, const QString& key, std::size_t& out) {
    if (!obj.contains(key))
        return false;
    const QJsonValue v = obj.value(key);
    if (!v.isDouble())
        return false;
    const double d = v.toDouble();
    if (d < 0.0)
        return false;
    out = static_cast<std::size_t>(d);
    return true;
}

QString JoinAxisAlignmentValues(const QJsonObject& obj) {
    // Build a human-readable summary of the axis_alignment block. The
    // manifest emits one descriptive sentence per axis kind; we
    // concatenate them with newlines for logging / inspector display.
    QStringList lines;
    for (auto it = obj.constBegin(); it != obj.constEnd(); ++it) {
        if (it.value().isString()) {
            lines << QStringLiteral("[%1] %2").arg(it.key(), it.value().toString());
        }
    }
    return lines.join(QStringLiteral("\n"));
}

}  // namespace


QtManifest QtManifest::Load(const QString& path) {
    QtManifest m;

    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) {
        m.error = QStringLiteral("QtManifest: could not open %1").arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtManifest"),
                                                m.error,
                                                path);
        return m;
    }
    const QByteArray bytes = f.readAll();
    f.close();

    QJsonParseError perr{};
    const QJsonDocument doc = QJsonDocument::fromJson(bytes, &perr);
    if (perr.error != QJsonParseError::NoError || !doc.isObject()) {
        m.error = QStringLiteral("QtManifest: JSON parse failed at offset %1: %2").arg(perr.offset).arg(perr.errorString());
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtManifest"),
                                                m.error,
                                                path);
        return m;
    }
    const QJsonObject root = doc.object();

    // Required top-level
    if (!ReadString(root, QStringLiteral("schema_version"), m.schemaVersion)
        || !ReadString(root, QStringLiteral("protein_id"), m.proteinId)) {
        m.error = QStringLiteral("QtManifest: missing required keys schema_version / protein_id in %1").arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtManifest"),
                                                m.error,
                                                path);
        return m;
    }
    ReadString(root, QStringLiteral("extractor"), m.extractor);
    ReadString(root, QStringLiteral("extractor_version"), m.extractorVersion);
    ReadString(root, QStringLiteral("generated_at_utc"), m.generatedAtUtc);

    // Topology block (optional; missing → all-false stays default)
    if (root.contains(QStringLiteral("topology"))) {
        const QJsonObject t = root.value(QStringLiteral("topology")).toObject();
        ReadString(t, QStringLiteral("source"), m.topology.source);
        ReadBool(t, QStringLiteral("has_atom_semantic"), m.topology.hasAtomSemantic);
        ReadBool(t, QStringLiteral("has_ff_atom_types"), m.topology.hasFfAtomTypes);
        ReadBool(t, QStringLiteral("has_ff_mass"), m.topology.hasFfMass);
    }

    // Axis sizes (atom is required; rest optional but expected)
    if (!root.contains(QStringLiteral("axis_sizes"))) {
        m.error = QStringLiteral("QtManifest: missing axis_sizes block in %1").arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtManifest"),
                                                m.error,
                                                path);
        return m;
    }
    const QJsonObject ax = root.value(QStringLiteral("axis_sizes")).toObject();
    if (!ReadSizeT(ax, QStringLiteral("atom"), m.axisSizes.atom)) {
        m.error = QStringLiteral("QtManifest: missing axis_sizes.atom in %1").arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtManifest"),
                                                m.error,
                                                path);
        return m;
    }
    ReadSizeT(ax, QStringLiteral("residue"), m.axisSizes.residue);
    ReadSizeT(ax, QStringLiteral("bond"), m.axisSizes.bond);
    ReadSizeT(ax, QStringLiteral("aromatic_ring"), m.axisSizes.aromaticRing);
    ReadSizeT(ax, QStringLiteral("saturated_ring"), m.axisSizes.saturatedRing);
    ReadSizeT(ax, QStringLiteral("ring"), m.axisSizes.ring);
    ReadSizeT(ax, QStringLiteral("ring_membership"), m.axisSizes.ringMembership);

    // Axis-alignment block (optional, informational)
    if (root.contains(QStringLiteral("axis_alignment"))) {
        m.axisAlignmentSummary = JoinAxisAlignmentValues(root.value(QStringLiteral("axis_alignment")).toObject());
    }

    m.ok = true;
    return m;
}

}  // namespace h5reader::io

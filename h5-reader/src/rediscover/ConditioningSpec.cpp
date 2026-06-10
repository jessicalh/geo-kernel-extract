#include "ConditioningSpec.h"

#include <QCryptographicHash>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>

namespace h5reader::rediscover {

ConditioningSpec ConditioningSpec::Default() {
    ConditioningSpec s;
    s.root.insert(QStringLiteral("schema_version"), 1);
    s.root.insert(QStringLiteral("grain"), QStringLiteral("shielding_sigma_iso_T0_primary_T2_sidecar"));
    s.root.insert(QStringLiteral("target"), QStringLiteral("absolute_dft_shielding"));
    s.root.insert(QStringLiteral("angle_convention"), QStringLiteral("IUPAC_signed_radians"));
    s.root.insert(QStringLiteral("no_subsampling"), true);
    s.root.insert(QStringLiteral("spatial_cutoffs_A"), QJsonArray{4, 6, 8, 10});
    return s;
}

std::optional<ConditioningSpec> ConditioningSpec::Load(const QString& path, QString* err_out) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) {
        if (err_out) *err_out = QStringLiteral("cannot open conditioning spec: %1").arg(path);
        return std::nullopt;
    }
    QJsonParseError pe;
    const QJsonDocument doc = QJsonDocument::fromJson(f.readAll(), &pe);
    if (pe.error != QJsonParseError::NoError || !doc.isObject()) {
        if (err_out) *err_out = QStringLiteral("conditioning spec JSON parse failed: %1").arg(pe.errorString());
        return std::nullopt;
    }
    ConditioningSpec s;
    s.sourcePath = path;
    s.root = doc.object();
    return s;
}

QByteArray ConditioningSpec::canonicalJson() const {
    return QJsonDocument(root).toJson(QJsonDocument::Compact);
}

QString ConditioningSpec::configHash() const {
    return QString::fromLatin1(
        QCryptographicHash::hash(canonicalJson(), QCryptographicHash::Sha256).toHex());
}

}  // namespace h5reader::rediscover

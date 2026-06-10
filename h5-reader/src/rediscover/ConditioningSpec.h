#pragma once

#include <QJsonObject>
#include <QString>

#include <optional>

namespace h5reader::rediscover {

struct ConditioningSpec {
    QString sourcePath;
    QJsonObject root;
    bool fixture = false;

    static ConditioningSpec Default();
    static std::optional<ConditioningSpec> Load(const QString& path, QString* err_out);
    QByteArray canonicalJson() const;
    QString configHash() const;
};

}  // namespace h5reader::rediscover

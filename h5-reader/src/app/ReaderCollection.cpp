#include "ReaderCollection.h"

#include "io/CalcsetManifest.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonParseError>
#include <QLoggingCategory>
#include <QSet>

#include <utility>

namespace h5reader::app {
namespace {
Q_LOGGING_CATEGORY(collectionLog, "h5reader.collection")
}

std::optional<ReaderDocument> ReadReaderDocument(const QString& path, QString* error) {
    const auto lgsPath = io::ResolveLgsPath(path, error);
    if (!lgsPath)
        return std::nullopt;

    QFile file(*lgsPath);
    if (!file.open(QIODevice::ReadOnly)) {
        if (error)
            *error = QStringLiteral("Cannot open .LGS file %1: %2").arg(*lgsPath, file.errorString());
        return std::nullopt;
    }
    QJsonParseError parseError;
    const QJsonDocument json = QJsonDocument::fromJson(file.readAll(), &parseError);
    if (parseError.error != QJsonParseError::NoError || !json.isObject()) {
        if (error)
            *error = QStringLiteral("Invalid .LGS file %1: %2")
                         .arg(*lgsPath, parseError.error == QJsonParseError::NoError
                                            ? QStringLiteral("expected a JSON object")
                                            : parseError.errorString());
        return std::nullopt;
    }

    ReaderDocument document;
    document.lgsPath = *lgsPath;
    const QJsonObject root = json.object();
    if (root.value(QStringLiteral("kind")).toString() != QStringLiteral("collection"))
        return document;

    document.collection = ParseReaderCollection(root, QUrl::fromLocalFile(*lgsPath), error);
    if (!document.collection)
        return std::nullopt;
    return document;
}

std::optional<ReaderCollection> ParseReaderCollection(const QJsonObject& root, const QUrl& source, QString* error) {
    if (!source.isValid() || (!source.isLocalFile() && (source.scheme() != QStringLiteral("https") || source.host().isEmpty()))) {
        if (error)
            *error = QStringLiteral("A collection source must be a local file or an HTTPS URL.");
        return std::nullopt;
    }
    const QJsonArray rows = root.value(QStringLiteral("entries")).toArray();
    ReaderCollection collection;
    collection.title = root.value(QStringLiteral("title")).toString();
    collection.source = source;
    if (root.value(QStringLiteral("kind")).toString() != QStringLiteral("collection")
        || root.value(QStringLiteral("schema_version")).toInt() != 1 || collection.title.isEmpty() || rows.isEmpty()) {
        if (error)
            *error = QStringLiteral("Collection .LGS needs kind collection, schema_version 1, a title, and entries: %1")
                         .arg(source.toDisplayString());
        return std::nullopt;
    }

    QSet<QString> keys;
    for (int i = 0; i < rows.size(); ++i) {
        const QJsonObject row = rows[i].toObject();
        const QString relative = row.value(QStringLiteral("archive")).toString();
        const QString clean = QDir::cleanPath(relative);
        ReaderCollection::Entry entry;
        entry.key = row.value(QStringLiteral("key")).toString();
        entry.title = row.value(QStringLiteral("title")).toString();
        entry.group = row.value(QStringLiteral("group")).toString();
        entry.description = row.value(QStringLiteral("description")).toString();
        entry.organism = row.value(QStringLiteral("organism")).toString();
        entry.bmrb = row.value(QStringLiteral("bmrb")).toString();
        entry.pdb = row.value(QStringLiteral("pdb")).toString();
        for (const auto& keyword : row.value(QStringLiteral("keywords")).toArray())
            entry.keywords.append(keyword.toString());
        entry.frames = row.value(QStringLiteral("frames")).toInt();
        entry.entryPoint = row.value(QStringLiteral("entry_point")).toString();
        entry.archiveBytes = row.value(QStringLiteral("archive_bytes")).toInteger(-1);
        entry.expandedBytes = row.value(QStringLiteral("expanded_bytes")).toInteger(-1);
        if (entry.key.isEmpty() || entry.key == QStringLiteral(".") || entry.key == QStringLiteral("..")
            || entry.key.contains(QLatin1Char('/')) || entry.key.contains(QLatin1Char('\\'))
            || entry.key.contains(QLatin1Char(':')) || keys.contains(entry.key.toCaseFolded())
            || entry.title.isEmpty() || relative.isEmpty() || QDir::isAbsolutePath(relative)
            || relative.contains(QLatin1Char('\\')) || relative.contains(QLatin1Char(':'))
            || relative.contains(QLatin1Char('?')) || relative.contains(QLatin1Char('#'))
            || clean == QStringLiteral("..") || clean.startsWith(QStringLiteral("../"))
            || !clean.endsWith(QStringLiteral(".tar.xz"), Qt::CaseInsensitive)
            || entry.entryPoint != QStringLiteral("run.LGS") || entry.frames < 1
            || entry.archiveBytes < 1 || entry.expandedBytes < 1) {
            if (error)
                *error = QStringLiteral("Invalid collection entry %1 in %2").arg(i + 1).arg(source.toDisplayString());
            return std::nullopt;
        }
        if (source.isLocalFile()) {
            const QDir base(QFileInfo(source.toLocalFile()).absolutePath());
            entry.archiveUrl = QUrl::fromLocalFile(base.absoluteFilePath(clean));
        } else {
            QUrl relativeUrl;
            relativeUrl.setPath(clean);
            entry.archiveUrl = source.resolved(relativeUrl);
        }
        keys.insert(entry.key.toCaseFolded());
        collection.entries.append(std::move(entry));
    }
    qCInfo(collectionLog).noquote() << "Collection loaded" << source.toDisplayString() << "entries=" << collection.entries.size();
    return collection;
}

}  // namespace h5reader::app

#pragma once

#include <QString>
#include <QStringList>
#include <QUrl>
#include <QVector>

#include <optional>

class QJsonObject;

namespace h5reader::app {

struct ReaderCollection {
    struct Entry {
        QString key;
        QString title;
        QString group;
        QString description;
        QString organism;
        QString bmrb;
        QString pdb;
        QStringList keywords;
        QUrl archiveUrl;
        QString entryPoint;
        int frames = 0;
        qint64 archiveBytes = 0;
        qint64 expandedBytes = 0;
    };

    QString title;
    QUrl source;
    QVector<Entry> entries;
};

struct ReaderDocument {
    QString lgsPath;
    std::optional<ReaderCollection> collection;
};

std::optional<ReaderDocument> ReadReaderDocument(const QString& path, QString* error);
std::optional<ReaderCollection> ParseReaderCollection(const QJsonObject& root, const QUrl& source, QString* error);

}  // namespace h5reader::app

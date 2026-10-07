#include "ResumeFile.h"

#include <QDateTime>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLocale>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QSaveFile>
#include <filesystem>

namespace sciencefiles {
namespace {

bool strongEtag(const QByteArray &value) {
    if (value.size() < 2 || !value.startsWith('"') || !value.endsWith('"'))
        return false;
    for (qsizetype i = 1; i < value.size() - 1; ++i) {
        const auto c = static_cast<unsigned char>(value[i]);
        if (c < 0x21 || c == '"' || c == 0x7f)
            return false;
    }
    return true;
}

bool strongDate(const QByteArray &modified, const QByteArray &date) {
    if (modified.contains('\r') || modified.contains('\n') || date.contains('\r') || date.contains('\n'))
        return false;
    // HTTP dates use GMT, which Qt::RFC2822Date does not accept.
    const auto format = QStringLiteral("ddd, dd MMM yyyy HH:mm:ss t");
    const auto modifiedTime = QLocale::c().toDateTime(QString::fromLatin1(modified), format);
    const auto responseTime = QLocale::c().toDateTime(QString::fromLatin1(date), format);
    return modifiedTime.isValid() && responseTime.isValid() && modifiedTime.secsTo(responseTime) >= 60;
}

bool byteCount(const QByteArray &text, qint64 &value) {
    if (text.isEmpty())
        return false;
    for (char c : text) {
        if (c < '0' || c > '9')
            return false;
    }
    bool ok = false;
    value = text.toLongLong(&ok);
    return ok;
}

bool contentRange(const QByteArray &value, qint64 &first, qint64 &last, qint64 &total) {
    if (!value.startsWith("bytes "))
        return false;
    const auto dash = value.indexOf('-', 6);
    const auto slash = value.indexOf('/', dash + 1);
    return dash > 6 && slash > dash + 1 &&
           byteCount(value.mid(6, dash - 6), first) &&
           byteCount(value.mid(dash + 1, slash - dash - 1), last) &&
           byteCount(value.mid(slash + 1), total) && first <= last && last < total;
}

} // namespace

ResumeFile::ResumeFile(const QUrl &url, const QString &destination)
    : url_(url), destination_(destination), metadataPath_(destination + ".json"), part_(destination + ".part") {}

bool ResumeFile::fail(Download::Result result, const QString &message) {
    errorResult_ = result;
    error_ = message;
    return false;
}

bool ResumeFile::discard() {
    part_.close();
    // Remove metadata first: an interrupted cleanup must not leave a trusted prefix.
    for (const auto &path : {metadataPath_, part_.fileName()}) {
        if (QFileInfo::exists(path) && !QFile::remove(path))
            return fail(Download::Result::FileError, QStringLiteral("Cannot remove resume file: %1").arg(path));
    }
    etag_.clear();
    modified_.clear();
    date_.clear();
    offset_ = 0;
    length_ = -1;
    return true;
}

bool ResumeFile::prepare() {
    const QFileInfo destination(destination_);
    if (destination.exists() || destination.isSymLink())
        return fail(Download::Result::FileError, QStringLiteral("Resume destination already exists: %1").arg(destination_));
    for (const auto &path : {metadataPath_, part_.fileName()}) {
        const QFileInfo file(path);
        if (file.isSymLink() || (file.exists() && !file.isFile()))
            return fail(Download::Result::FileError, QStringLiteral("Not a regular resume file: %1").arg(path));
    }
    QFile metadata(metadataPath_);
    if (metadata.exists()) {
        if (!metadata.open(QIODevice::ReadOnly))
            return fail(Download::Result::FileError, metadata.errorString());
        const auto object = QJsonDocument::fromJson(metadata.readAll()).object();
        metadata.close();
        etag_ = object.value("etag").toString().toLatin1();
        modified_ = object.value("last_modified").toString().toLatin1();
        date_ = object.value("date").toString().toLatin1();
        length_ = object.value("length").toInteger(-1);
        offset_ = QFileInfo(part_).size();
        const bool validator = strongEtag(etag_) || (etag_.isEmpty() && strongDate(modified_, date_));
        if (object.value("url").toString() == url_.toString(QUrl::FullyEncoded) && validator &&
            part_.exists() && offset_ > 0 && offset_ <= length_) {
            prepared_ = true;
            return true;
        }
    }
    prepared_ = discard();
    return prepared_;
}

void ResumeFile::setRequestHeaders(QNetworkRequest &request) const {
    if (offset_ > 0) {
        request.setRawHeader("Range", "bytes=" + QByteArray::number(offset_) + '-');
        request.setRawHeader("If-Range", etag_.isEmpty() ? modified_ : etag_);
    }
}

bool ResumeFile::openResponse(const QNetworkReply &reply) {
    using Result = Download::Result;
    const int status = reply.attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
    if (status == 416) {
        if (!discard())
            return false;
        return fail(Result::HttpError, QStringLiteral("HTTP 416. Partial download discarded; retry from the beginning."));
    }
    if (status != 200 && status != 206)
        return fail(Result::HttpError, QStringLiteral("HTTP %1").arg(status));
    const QByteArray encoding = reply.rawHeader("Content-Encoding").trimmed().toLower();
    if (!encoding.isEmpty() && encoding != "identity")
        return fail(Result::HttpError, QStringLiteral("Resume requires an identity-encoded response."));

    const QByteArray etag = reply.rawHeader("ETag").trimmed();
    const QByteArray modified = reply.rawHeader("Last-Modified").trimmed();
    const QByteArray date = reply.rawHeader("Date").trimmed();
    qint64 contentLength = -1;
    const bool hasLength = reply.hasRawHeader("Content-Length");
    if (hasLength && !byteCount(reply.rawHeader("Content-Length").trimmed(), contentLength))
        return fail(Result::HttpError, QStringLiteral("Invalid Content-Length."));
    if (status == 206) {
        qint64 first = -1, last = -1, total = -1;
        if (offset_ == 0 || !contentRange(reply.rawHeader("Content-Range").trimmed(), first, last, total) ||
            first != offset_ || total != length_ || (hasLength && contentLength != last - first + 1))
            return fail(Result::HttpError, QStringLiteral("Content-Range does not match the saved partial download."));
        // If-Range permits omission of representation headers; reject any supplied mismatch.
        if ((!etag_.isEmpty() && !etag.isEmpty() && etag != etag_) ||
            (etag_.isEmpty() && !modified.isEmpty() && modified != modified_))
            return fail(Result::HttpError, QStringLiteral("The partial response has a different validator."));
        responseEnd_ = last + 1;
        if (!part_.open(QIODevice::WriteOnly | QIODevice::Append))
            return fail(Result::FileError, part_.errorString());
        return true;
    }

    if (!hasLength)
        return fail(Result::HttpError, QStringLiteral("Resume requires Content-Length on a full response."));
    // A 200 response replaces the representation, including its previous validator.
    if (!discard())
        return false;
    length_ = contentLength;
    responseEnd_ = length_;
    if (strongEtag(etag))
        etag_ = etag;
    else if (etag.isEmpty() && strongDate(modified, date)) {
        modified_ = modified;
        date_ = date;
    }
    if (!part_.open(QIODevice::WriteOnly | QIODevice::Truncate))
        return fail(Result::FileError, part_.errorString());
    if (!etag_.isEmpty() || !modified_.isEmpty()) {
        const QJsonObject object{{"url", url_.toString(QUrl::FullyEncoded)}, {"etag", QString::fromLatin1(etag_)},
                                 {"last_modified", QString::fromLatin1(modified_)},
                                 {"date", QString::fromLatin1(date_)}, {"length", length_}};
        const QByteArray bytes = QJsonDocument(object).toJson(QJsonDocument::Compact);
        QSaveFile metadata(metadataPath_);
        metadata.setDirectWriteFallback(false);
        if (!metadata.open(QIODevice::WriteOnly) || metadata.write(bytes) != bytes.size() || !metadata.commit())
            return fail(Result::FileError, metadata.errorString());
    }
    return true;
}

bool ResumeFile::write(const char *data, qint64 size) {
    if (size > responseEnd_ - part_.pos())
        return fail(Download::Result::HttpError, QStringLiteral("Response exceeds the declared byte range."));
    if (part_.write(data, size) != size)
        return fail(Download::Result::FileError, part_.errorString());
    return true;
}

bool ResumeFile::closeFile() {
    if (part_.isOpen()) {
        const bool flushed = part_.flush();
        const QString detail = part_.errorString();
        part_.close();
        if (!flushed) {
            discard();
            return fail(Download::Result::FileError, detail);
        }
    }
    return true;
}

bool ResumeFile::finish() {
    if (!closeFile())
        return false;
    if (prepared_ && etag_.isEmpty() && modified_.isEmpty())
        return discard();
    return true;
}

bool ResumeFile::commit() {
    if (part_.size() != length_)
        return fail(Download::Result::HttpError, QStringLiteral("Incomplete download; retry to continue."));
    if (!closeFile())
        return false;
    const QFileInfo destination(destination_);
    if (destination.exists() || destination.isSymLink())
        return fail(Download::Result::FileError, QStringLiteral("Destination already exists: %1").arg(destination_));
    // Unlike QFile::rename, filesystem::rename never falls back to a full file copy.
    std::error_code error;
    std::filesystem::rename(QFileInfo(part_).filesystemFilePath(), destination.filesystemFilePath(), error);
    if (error)
        return fail(Download::Result::FileError, QString::fromStdString(error.message()));
    if (QFileInfo::exists(metadataPath_) && !QFile::remove(metadataPath_))
        return fail(Download::Result::FileError, QStringLiteral("Saved download, but cannot remove metadata: %1").arg(metadataPath_));
    return true;
}

} // namespace sciencefiles

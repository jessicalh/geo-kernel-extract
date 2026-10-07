#pragma once

#include <sciencefiles/Download.h>
#include <QFile>

class QNetworkReply;
class QNetworkRequest;

namespace sciencefiles {

// Private, single-writer download files. All calls run in the download thread.
class ResumeFile final {
public:
    ResumeFile(const QUrl &url, const QString &destination);
    bool prepare();
    void setRequestHeaders(QNetworkRequest &request) const;
    bool openResponse(const QNetworkReply &reply);
    bool write(const char *data, qint64 size);
    bool commit();
    bool finish();

    qint64 offset() const { return offset_; }
    qint64 length() const { return length_; }
    Download::Result errorResult() const { return errorResult_; }
    QString errorString() const { return error_; }

private:
    bool closeFile();
    bool discard();
    bool fail(Download::Result result, const QString &message);

    QUrl url_;
    QString destination_;
    QString metadataPath_;
    QFile part_;
    QByteArray etag_;
    QByteArray modified_;
    QByteArray date_;
    qint64 offset_ = 0;
    qint64 length_ = -1;
    qint64 responseEnd_ = -1;
    bool prepared_ = false;
    Download::Result errorResult_ = Download::Result::FileError;
    QString error_;
};

} // namespace sciencefiles

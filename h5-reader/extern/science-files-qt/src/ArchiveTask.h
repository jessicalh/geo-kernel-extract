#pragma once

#include <sciencefiles/BundleCache.h>
#include <QThread>

namespace sciencefiles {

// Inputs are set before start(); the owner reads the outcome only after wait().
class ArchiveTask final : public QThread {
    Q_OBJECT
public:
    ArchiveTask(std::unique_ptr<QTemporaryDir> staging, const Bundle &bundle, const QString &destination, bool extract,
                BundleCache::Result result, const QString &detail);
    BundleCache::Result result() const { return result_; }
    QString detail() const { return detail_; }
    QString publishedPath() const { return publishedPath_; }

signals:
    void progress(qint64 completedBytes, qint64 totalBytes);
    void cleaning();

protected:
    void run() override;

private:
    struct ExtractionResult {
        BundleCache::Result result;
        QString detail;
    };
    ExtractionResult extractArchive();
    std::unique_ptr<QTemporaryDir> staging_;
    Bundle bundle_;
    QString destination_;
    bool extract_;
    BundleCache::Result result_;
    QString detail_;
    QString publishedPath_;
};

} // namespace sciencefiles

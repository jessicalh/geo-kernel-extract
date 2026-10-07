#pragma once

#include "ReaderCollection.h"

#include <sciencefiles/BundleCache.h>
#include <QDialog>
#include <QJsonObject>
#include <QNetworkAccessManager>
#include <QPointer>

#include <memory>

class QLabel;
class QLineEdit;
class QPushButton;
class QProgressBar;
class QTableWidget;
class QNetworkReply;

namespace h5reader::app {

class ReaderCollectionDialog final : public QDialog {
    Q_OBJECT

public:
    explicit ReaderCollectionDialog(ReaderCollection collection = {}, QWidget* parent = nullptr,
                                    const QString& bufferRoot = {}, const QString& installedRoot = {});
    bool setCollection(ReaderCollection collection, QString* error = nullptr);
    bool loadCatalog(const QUrl& source, QString* error = nullptr);
    bool openTrajectory(const QString& key, QString* error);
    bool clearTrajectory(const QString& key, QString* error);
    void setCurrentRun(const QString& path);
    QJsonObject state() const;
    bool isClearing() const;
    void cancel();
    void showLoadError(const QString& message);
    void runOpened(const QString& key);
    void runFailed(const QString& key, const QString& message);
    void cancelPendingOpen();
    void shutdown();
    bool isShutdown() const;
    bool isBusy() const;

signals:
    void openRequested(const QString& lgsPath, const QString& key);
    void shutdownFinished();

private:
    void refresh();
    void filterRows();
    void finishShutdown();
    void openSelected();
    const ReaderCollection::Entry* findEntry(const QString& key) const;
    QString installedPath(const ReaderCollection::Entry& entry) const;
    sciencefiles::Bundle bundleFor(const ReaderCollection::Entry& entry) const;
    void reject() override;

    ReaderCollection collection_;
    QNetworkAccessManager network_;
    QPointer<QNetworkReply> catalogReply_;
    QUrl retryCatalog_;
    std::unique_ptr<sciencefiles::BundleCache> buffer_;
    QString bufferRoot_;
    QString installedRoot_;
    QString currentRun_;
    QString lastError_;
    QString openingKey_;
    QString retryKey_;
    QString currentKey_;
    QString pendingLoadError_;
    bool closing_ = false;
    bool shutdownComplete_ = false;
    QLineEdit* search_ = nullptr;
    QTableWidget* table_ = nullptr;
    QLabel* detail_ = nullptr;
    QLabel* source_ = nullptr;
    QLabel* status_ = nullptr;
    QProgressBar* progress_ = nullptr;
    QPushButton* open_ = nullptr;
    QPushButton* cancel_ = nullptr;
};

}  // namespace h5reader::app

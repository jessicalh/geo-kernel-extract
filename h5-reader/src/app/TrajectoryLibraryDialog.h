#pragma once

#include <QDialog>
#include <QJsonObject>
#include <memory>
#include <sciencefiles/BundleCache.h>

class QLabel;
class QLineEdit;
class QProgressBar;
class QPushButton;
class QTableWidget;

namespace h5reader::app {

class TrajectoryLibraryDialog final : public QDialog {
    Q_OBJECT
public:
    explicit TrajectoryLibraryDialog(QWidget* parent = nullptr,
                                     const QString& catalogPath = {},
                                     const QString& cacheRoot = {},
                                     const QString& installedRoot = {});
    bool openTrajectory(const QString& key, QString* error);
    bool clearTrajectory(const QString& key, QString* error);
    bool isClearing() const;
    void setCurrentRun(const QString& path);
    void showLoadError(const QString& message);
    void shutdown();
    bool isShutdown() const;
    QJsonObject state() const;

public slots:
    void cancel();

signals:
    void openRequested(const QString& path);
    void shutdownFinished();

protected:
    void reject() override;

private:
    struct Entry {
        sciencefiles::Bundle bundle;
        QString title;
        QString description;
        int frames = 0;
        qint64 archiveBytes = 0;
    };
    bool readCatalog(const QString& path);
    void createCache(const QString& root);
    void refresh();
    void chooseCacheRoot();
    void confirmClear();
    const Entry* findEntry(const QString& key) const;
    const Entry* selectedEntry() const;
    QString installedPath(const Entry& entry) const;
    bool isCurrent(const QString& directory) const;

    QList<Entry> entries_;
    std::unique_ptr<sciencefiles::BundleCache> cache_;
    QString cacheRoot_;
    QString installedRoot_;
    QString currentRun_;
    QString openingKey_;
    QString pendingRoot_;
    QString lastError_;
    bool closing_ = false;
    QLineEdit* search_ = nullptr;
    QTableWidget* table_ = nullptr;
    QLabel* detail_ = nullptr;
    QLabel* status_ = nullptr;
    QLabel* location_ = nullptr;
    QProgressBar* progress_ = nullptr;
    QPushButton* open_ = nullptr;
    QPushButton* cancel_ = nullptr;
    QPushButton* clear_ = nullptr;
    QPushButton* locationButton_ = nullptr;
};

}  // namespace h5reader::app

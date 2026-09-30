#include "TrajectoryLibraryDialog.h"

#ifdef Q_OS_LINUX

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLabel>
#include <QLocale>
#include <QLineEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QRegularExpression>
#include <QSet>
#include <QTableWidget>

namespace h5reader::app {
namespace {
bool within(const QString& path, const QString& root) {
    return path == root || path.startsWith(root + QLatin1Char('/'));
}
bool relativePath(const QString& path) {
    return !path.isEmpty() && !QDir::isAbsolutePath(path) && !path.contains(QLatin1Char('\\'))
        && path != QStringLiteral(".") && path != QStringLiteral("..")
        && !path.startsWith(QStringLiteral("../")) && QDir::cleanPath(path) == path;
}
}

bool TrajectoryLibraryDialog::readLocalCatalog(const QString& path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        showLoadError(tr("Cannot read the local trajectory catalog: %1").arg(file.errorString()));
        return false;
    }
    QJsonParseError parse;
    const QJsonDocument document = QJsonDocument::fromJson(file.readAll(), &parse);
    const QJsonObject catalog = document.object();
    if (file.error() != QFileDevice::NoError || parse.error != QJsonParseError::NoError
        || catalog.value("schema_version").toInt() != 1
        || catalog.value("kind").toString() != QStringLiteral("h5reader-local-library")
        || !catalog.value("datasets").isArray() || catalog.value("datasets").toArray().isEmpty()) {
        showLoadError(tr("Invalid local trajectory catalog: %1").arg(path));
        return false;
    }
    static const QRegularExpression validKey(QStringLiteral("^[A-Za-z0-9][A-Za-z0-9_.-]{0,159}$"));
    QSet<QString> keys;
    QSet<QString> paths;
    QList<Entry> entries;
    for (const QJsonValue& value : catalog.value("datasets").toArray()) {
        const QJsonObject row = value.toObject();
        Entry entry;
        entry.bundle.key = row.value("key").toString();
        entry.bundle.entryPoint = QStringLiteral("run.LGS");
        entry.title = row.value("title").toString();
        entry.description = row.value("description").toString();
        entry.frames = row.value("frames").toInt();
        entry.dftFrames = row.value("dft_frames").toInt(-1);
        entry.sourceLgs = row.value("source_lgs").toString();
        if (!validKey.match(entry.bundle.key).hasMatch() || entry.title.isEmpty() || entry.frames < 1
            || entry.dftFrames < 0 || entry.dftFrames > entry.frames || !relativePath(entry.sourceLgs)
            || keys.contains(entry.bundle.key) || paths.contains(entry.sourceLgs)) {
            showLoadError(tr("Invalid or duplicate local trajectory: %1").arg(entry.bundle.key));
            return false;
        }
        keys.insert(entry.bundle.key);
        paths.insert(entry.sourceLgs);
        entries.append(entry);
    }
    entries_ = entries;
    search_->setPlaceholderText(tr("Search all %1 trajectories").arg(entries_.size()));
    table_->setRowCount(static_cast<int>(entries_.size()));
    for (int row = 0; row < entries_.size(); ++row) {
        const auto& entry = entries_[row];
        table_->setItem(row, 0, new QTableWidgetItem(entry.title));
        table_->setItem(row, 1, new QTableWidgetItem(QLocale().toString(entry.frames)));
        table_->setItem(row, 2, new QTableWidgetItem(QLocale().toString(entry.dftFrames)));
        table_->setItem(row, 3, new QTableWidgetItem);
        table_->item(row, 0)->setToolTip(entry.description);
    }
    return true;
}

QString TrajectoryLibraryDialog::localSourcePath(const Entry& entry, QString* error) const {
    const QString root = QFileInfo(sourceRoot_).canonicalFilePath();
    const QFileInfo file(QDir(sourceRoot_).filePath(entry.sourceLgs));
    const QString path = file.canonicalFilePath();
    if (root.isEmpty() || root == QStringLiteral("/") || !file.isFile() || path.isEmpty()) {
        if (error) *error = tr("This trajectory is unavailable. Check that the data drive is connected: %1").arg(entry.sourceLgs);
        return {};
    }
    if (!within(path, root)) {
        if (error) *error = tr("The trajectory path leaves the data drive.");
        return {};
    }
    return path;
}

void TrajectoryLibraryDialog::refreshLocal() {
    for (int row = 0; row < entries_.size(); ++row)
        table_->item(row, 3)->setText(tr("Data drive"));
    const Entry* entry = selectedEntry();
    open_->setEnabled(!closing_ && entry);
    open_->setText(tr("Open"));
    detail_->setText(entry ? entry->description : QString());
    progress_->hide();
}

} // namespace h5reader::app
#endif

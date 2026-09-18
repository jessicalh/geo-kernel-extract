#include "TrajectoryLibraryDialog.h"

#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QHeaderView>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QLoggingCategory>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QSet>
#include <QSettings>
#include <QStandardPaths>
#include <QStyle>
#include <QTableWidget>
#include <QVBoxLayout>
#include <utility>

namespace h5reader::app {
namespace {
Q_LOGGING_CATEGORY(libraryLog, "h5reader.library")
using Cache = sciencefiles::BundleCache;
const QString cacheSetting = QStringLiteral("viewer/library/cacheRoot");

QString defaultInstalledRoot() {
    return QDir(QCoreApplication::applicationDirPath()).filePath(QStringLiteral("datasets"));
}
}  // namespace

TrajectoryLibraryDialog::TrajectoryLibraryDialog(QWidget* parent,
                                                 const QString& catalogPath,
                                                 const QString& cacheRoot,
                                                 const QString& installedRoot)
    : QDialog(parent)
    , installedRoot_(installedRoot.isEmpty() ? defaultInstalledRoot() : installedRoot) {
    setObjectName(QStringLiteral("TrajectoryLibraryDialog"));
    setWindowTitle(tr("Published trajectories"));
    resize(800, 540);
    auto* layout = new QVBoxLayout(this);
    search_ = new QLineEdit(this);
    search_->setObjectName(QStringLiteral("trajectorySearch"));
    search_->setPlaceholderText(tr("Search trajectories"));
    search_->setClearButtonEnabled(true);
    layout->addWidget(search_);
    table_ = new QTableWidget(0, 4, this);
    table_->setObjectName(QStringLiteral("trajectoryCatalog"));
    table_->setHorizontalHeaderLabels({tr("Trajectory"), tr("Frames"), tr("Download"), tr("Local copy")});
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->setSelectionMode(QAbstractItemView::SingleSelection);
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_->verticalHeader()->hide();
    table_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    for (int column = 1; column < 4; ++column)
        table_->horizontalHeader()->setSectionResizeMode(column, QHeaderView::ResizeToContents);
    layout->addWidget(table_, 1);
    detail_ = new QLabel(this);
    detail_->setWordWrap(true);
    detail_->setMinimumHeight(40);
    layout->addWidget(detail_);
    progress_ = new QProgressBar(this);
    progress_->setRange(0, 1000);
    progress_->setTextVisible(false);
    layout->addWidget(progress_);
    status_ = new QLabel(this);
    status_->setWordWrap(true);
    status_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(status_);
    auto* buttons = new QHBoxLayout;
    open_ = new QPushButton(style()->standardIcon(QStyle::SP_DialogOpenButton), tr("Open"), this);
    cancel_ = new QPushButton(style()->standardIcon(QStyle::SP_DialogCancelButton), tr("Cancel"), this);
    clear_ = new QPushButton(style()->standardIcon(QStyle::SP_TrashIcon), tr("Clear download"), this);
    open_->setObjectName(QStringLiteral("openTrajectory"));
    cancel_->setObjectName(QStringLiteral("cancelDownload"));
    clear_->setObjectName(QStringLiteral("clearDownload"));
    for (auto* button : {open_, cancel_, clear_}) {
        button->setAutoDefault(false);
        buttons->addWidget(button);
    }
    buttons->addStretch();
    layout->addLayout(buttons);
    auto* cacheLine = new QHBoxLayout;
    location_ = new QLabel(this);
    location_->setWordWrap(true);
    locationButton_ = new QPushButton(style()->standardIcon(QStyle::SP_DirIcon), tr("Cache location..."), this);
    locationButton_->setObjectName(QStringLiteral("cacheLocation"));
    locationButton_->setAutoDefault(false);
    cacheLine->addWidget(location_, 1);
    cacheLine->addWidget(locationButton_);
    layout->addLayout(cacheLine);

    const QString catalog = catalogPath.isEmpty()
                                ? QDir(QCoreApplication::applicationDirPath()).filePath(QStringLiteral("trajectories.json"))
                                : catalogPath;
    readCatalog(catalog);
    QString root = cacheRoot;
    if (root.isEmpty())
        root = QSettings()
                   .value(cacheSetting,
                          QDir(QStandardPaths::writableLocation(QStandardPaths::CacheLocation))
                              .filePath(QStringLiteral("trajectories")))
                   .toString();
    createCache(root);
    connect(search_, &QLineEdit::textChanged, this, [this](const QString& query) {
        for (int row = 0; row < entries_.size(); ++row) {
            const auto& entry = entries_[row];
            table_->setRowHidden(
                row,
                !(entry.title + " " + entry.description + " " + entry.bundle.key).contains(query, Qt::CaseInsensitive));
        }
        refresh();
    });
    connect(table_, &QTableWidget::itemSelectionChanged, this, &TrajectoryLibraryDialog::refresh);
    connect(open_, &QPushButton::clicked, this, [this] {
        if (const auto* entry = selectedEntry()) {
            QString error;
            if (!openTrajectory(entry->bundle.key, &error))
                showLoadError(error);
        }
    });
    connect(table_, &QTableWidget::cellDoubleClicked, this, [this] { open_->click(); });
    connect(cancel_, &QPushButton::clicked, this, &TrajectoryLibraryDialog::cancel);
    connect(clear_, &QPushButton::clicked, this, &TrajectoryLibraryDialog::confirmClear);
    connect(locationButton_, &QPushButton::clicked, this, &TrajectoryLibraryDialog::chooseCacheRoot);
    if (!entries_.isEmpty())
        table_->selectRow(0);
    refresh();
}

bool TrajectoryLibraryDialog::readCatalog(const QString& path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        showLoadError(tr("Cannot read trajectory catalog: %1").arg(file.errorString()));
        return false;
    }
    QJsonParseError parse;
    const auto document = QJsonDocument::fromJson(file.readAll(), &parse);
    if (file.error() != QFileDevice::NoError || !document.isArray()) {
        showLoadError(tr("Invalid trajectory catalog: %1").arg(parse.errorString()));
        return false;
    }
    QSet<QString> keys;
    for (const auto& value : document.array()) {
        const auto row = value.toObject();
        Entry entry;
        entry.bundle = {row["key"].toString(),
                        QUrl(row["url"].toString()),
                        row["entry_point"].toString(),
                        row["expanded_bytes"].toInteger(-1)};
        entry.title = row["title"].toString();
        entry.description = row["description"].toString();
        entry.frames = row["frames"].toInt();
        entry.archiveBytes = row["archive_bytes"].toInteger(-1);
        if (entry.title.isEmpty() || entry.bundle.key.isEmpty() || entry.bundle.entryPoint.isEmpty()
            || entry.bundle.url.scheme() != "https" || entry.frames < 1 || entry.archiveBytes < 0
            || entry.bundle.expandedBytes < 0 || keys.contains(entry.bundle.key)) {
            entries_.clear();
            showLoadError(tr("Invalid or duplicate trajectory catalog entry: %1").arg(entry.bundle.key));
            return false;
        }
        keys.insert(entry.bundle.key);
        entries_.append(entry);
    }
    table_->setRowCount(static_cast<int>(entries_.size()));
    for (int row = 0; row < entries_.size(); ++row) {
        const auto& entry = entries_[row];
        table_->setItem(row, 0, new QTableWidgetItem(entry.title));
        table_->setItem(row, 1, new QTableWidgetItem(QString::number(entry.frames)));
        table_->setItem(row, 2, new QTableWidgetItem(QLocale().formattedDataSize(entry.archiveBytes)));
        table_->setItem(row, 3, new QTableWidgetItem);
    }
    qCInfo(libraryLog) << "Catalog loaded:" << path << "entries=" << entries_.size();
    return true;
}

void TrajectoryLibraryDialog::createCache(const QString& root) {
    cacheRoot_ = QDir(root).absolutePath();
    cache_ = std::make_unique<Cache>(cacheRoot_, this);
    location_->setText(tr("Download cache: %1").arg(QDir::toNativeSeparators(cacheRoot_)));
    connect(cache_.get(), &Cache::stateChanged, this, [this](Cache::State phase) {
        progress_->setRange(0, 0);
        progress_->setValue(0);
        if (phase == Cache::State::Downloading)
            status_->setText(tr("Downloading..."));
        if (phase == Cache::State::Extracting)
            status_->setText(tr("Unpacking..."));
        if (phase == Cache::State::Clearing)
            status_->setText(tr("Clearing downloaded data..."));
        if (phase == Cache::State::Cleaning)
            status_->setText(tr("Finishing and removing temporary files..."));
        refresh();
    });
    connect(cache_.get(), &Cache::progress, this, [this](qint64 bytes, qint64 total) {
        progress_->setRange(0, total > 0 ? 1000 : 0);
        if (total > 0)
            progress_->setValue(static_cast<int>(1000.0 * static_cast<double>(bytes) / static_cast<double>(total)));
    });
    connect(cache_.get(), &Cache::finished, this, [this](Cache::Result result, const QString& detail, const QString& path) {
        status_->setText(detail);
        lastError_ = result == Cache::Result::Ready || result == Cache::Result::Cleared || result == Cache::Result::Cancelled
                         ? QString()
                         : detail;
        const auto* entry = findEntry(openingKey_);
        openingKey_.clear();
        refresh();
        if (result == Cache::Result::Ready && entry && !closing_)
            emit openRequested(QDir(path).filePath(entry->bundle.entryPoint));
    });
    connect(
        cache_.get(),
        &Cache::shutdownFinished,
        this,
        [this] {
            if (!pendingRoot_.isEmpty() && !closing_) {
                const QString root = std::exchange(pendingRoot_, {});
                cache_.reset();
                createCache(root);
                QSettings().setValue(cacheSetting, root);
                status_->setText(tr("Cache location changed. Existing downloads were not moved or removed."));
                refresh();
            } else {
                emit shutdownFinished();
            }
        },
        Qt::QueuedConnection);
}

const TrajectoryLibraryDialog::Entry* TrajectoryLibraryDialog::findEntry(const QString& key) const {
    for (const auto& entry : entries_)
        if (entry.bundle.key == key)
            return &entry;
    return nullptr;
}

const TrajectoryLibraryDialog::Entry* TrajectoryLibraryDialog::selectedEntry() const {
    const int row = table_->currentRow();
    return row >= 0 && row < entries_.size() && !table_->isRowHidden(row) ? &entries_[row] : nullptr;
}

QString TrajectoryLibraryDialog::installedPath(const Entry& entry) const {
    const QString directory = QDir(installedRoot_).filePath(entry.bundle.key);
    return QFileInfo::exists(QDir(directory).filePath(entry.bundle.entryPoint)) ? directory : QString();
}

bool TrajectoryLibraryDialog::isCurrent(const QString& directory) const {
    if (directory.isEmpty() || currentRun_.isEmpty())
        return false;
    const QString relative = QDir(directory).relativeFilePath(currentRun_);
    return relative != ".." && !relative.startsWith("../") && !QDir::isAbsolutePath(relative);
}

void TrajectoryLibraryDialog::setCurrentRun(const QString& path) {
    currentRun_ = path;
    refresh();
}

bool TrajectoryLibraryDialog::openTrajectory(const QString& key, QString* error) {
    const auto* entry = findEntry(key);
    if (!entry || closing_ || !pendingRoot_.isEmpty() || cache_->busy()) {
        *error = tr("Unknown trajectory, or the library is busy.");
        return false;
    }
    lastError_.clear();
    const QString installed = installedPath(*entry);
    if (!installed.isEmpty()) {
        emit openRequested(QDir(installed).filePath(entry->bundle.entryPoint));
        return true;
    }
    openingKey_ = key;
    if (!cache_->fetch(entry->bundle)) {
        openingKey_.clear();
        *error = tr("The download cache is shutting down.");
        return false;
    }
    return true;
}

bool TrajectoryLibraryDialog::isClearing() const {
    return cache_->state() == Cache::State::Clearing;
}

bool TrajectoryLibraryDialog::clearTrajectory(const QString& key, QString* error) {
    const auto* entry = findEntry(key);
    if (!entry || closing_ || !pendingRoot_.isEmpty() || cache_->busy()) {
        *error = tr("Unknown trajectory, or the library is busy.");
        return false;
    }
    if (!installedPath(*entry).isEmpty()) {
        *error = tr("This example belongs to the installation and is not cleared with downloads.");
        return false;
    }
    if (isCurrent(QDir(cacheRoot_).filePath(key))) {
        *error = tr("Open a different trajectory before clearing this download.");
        return false;
    }
    lastError_.clear();
    if (!cache_->clear(key)) {
        *error = tr("The download cache is shutting down.");
        return false;
    }
    return true;
}

void TrajectoryLibraryDialog::refresh() {
    if (!cache_)
        return;
    const bool idle = !cache_->busy() && !closing_ && pendingRoot_.isEmpty();
    for (int row = 0; row < entries_.size(); ++row) {
        const auto& entry = entries_[row];
        const bool installed = !installedPath(entry).isEmpty();
        const bool downloaded = !cache_->cachedPath(entry.bundle).isEmpty();
        table_->item(row, 3)->setText(installed ? tr("Included") : downloaded ? tr("Downloaded") : tr("Not downloaded"));
    }
    const auto* entry = selectedEntry();
    const bool installed = entry && !installedPath(*entry).isEmpty();
    const bool downloaded = entry && !cache_->cachedPath(entry->bundle).isEmpty();
    open_->setEnabled(idle && entry);
    open_->setText(installed || downloaded ? tr("Open")
                   : lastError_.isEmpty()  ? tr("Download and open")
                                           : tr("Retry download"));
    clear_->setEnabled(idle && entry && downloaded && !installed && !isCurrent(QDir(cacheRoot_).filePath(entry->bundle.key)));
    cancel_->setEnabled(!closing_
                        && (cache_->state() == Cache::State::Downloading || cache_->state() == Cache::State::Extracting));
    locationButton_->setEnabled(idle);
    detail_->setText(entry ? entry->description : QString());
    progress_->setVisible(cache_->busy());
}

void TrajectoryLibraryDialog::confirmClear() {
    const auto* entry = selectedEntry();
    if (!entry)
        return;
    const QString key = entry->bundle.key;
    auto* question =
        new QMessageBox(QMessageBox::Question,
                        tr("Clear download"),
                        tr("Remove the downloaded copy of %1? Opening it again will require a new download.").arg(entry->title),
                        QMessageBox::Yes | QMessageBox::No,
                        this);
    question->setDefaultButton(QMessageBox::No);
    question->setAttribute(Qt::WA_DeleteOnClose);
    connect(question, &QMessageBox::finished, this, [this, key](int result) {
        if (result == QMessageBox::Yes) {
            QString error;
            if (!clearTrajectory(key, &error))
                showLoadError(error);
        }
    });
    question->open();
}

void TrajectoryLibraryDialog::chooseCacheRoot() {
    auto* chooser = new QFileDialog(this, tr("Choose download cache"), cacheRoot_);
    chooser->setFileMode(QFileDialog::Directory);
    chooser->setOption(QFileDialog::ShowDirsOnly);
    chooser->setAttribute(Qt::WA_DeleteOnClose);
    connect(chooser, &QFileDialog::fileSelected, this, [this](const QString& directory) {
        if (closing_ || cache_->busy() || !pendingRoot_.isEmpty())
            return;
        const QString root = QDir(directory).absolutePath();
        if (QDir(root) == QDir(cacheRoot_))
            return;
        pendingRoot_ = root;
        cache_->shutdown();
        refresh();
    });
    chooser->open();
}

void TrajectoryLibraryDialog::showLoadError(const QString& message) {
    lastError_ = message;
    status_->setText(message);
    qCWarning(libraryLog).noquote() << message;
}

void TrajectoryLibraryDialog::cancel() {
    openingKey_.clear();
    cache_->cancel();
}

void TrajectoryLibraryDialog::reject() {
    cancel();
    QDialog::reject();
}

void TrajectoryLibraryDialog::shutdown() {
    closing_ = true;
    pendingRoot_.clear();
    cache_->shutdown();
    refresh();
}

bool TrajectoryLibraryDialog::isShutdown() const {
    return cache_->isShutdown();
}

QJsonObject TrajectoryLibraryDialog::state() const {
    QJsonArray rows;
    for (const auto& entry : entries_) {
        const QString installed = installedPath(entry);
        const QString downloaded = cache_->cachedPath(entry.bundle);
        rows.append(QJsonObject{{"key", entry.bundle.key},
                                {"title", entry.title},
                                {"frames", entry.frames},
                                {"included", !installed.isEmpty()},
                                {"downloaded", !downloaded.isEmpty()},
                                {"directory", installed.isEmpty() ? downloaded : installed},
                                {"entry_point", entry.bundle.entryPoint},
                                {"archive_bytes", entry.archiveBytes},
                                {"expanded_bytes", entry.bundle.expandedBytes}});
    }
    return {{"entries", rows},
            {"cache_root", cacheRoot_},
            {"busy", cache_->busy()},
            {"phase", static_cast<int>(cache_->state())},
            {"shutdown", isShutdown()},
            {"error", lastError_},
            {"visible", isVisible()}};
}

}  // namespace h5reader::app

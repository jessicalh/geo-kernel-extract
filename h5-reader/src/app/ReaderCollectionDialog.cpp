#include "ReaderCollectionDialog.h"

#include <QAbstractItemView>
#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QDateTime>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QLoggingCategory>
#include <QMetaEnum>
#include <QNetworkReply>
#include <QProgressBar>
#include <QPushButton>
#include <QStandardPaths>
#include <QStyle>
#include <QTableWidget>
#include <QUrl>
#include <QVBoxLayout>

#include <utility>

namespace h5reader::app {
namespace {
Q_LOGGING_CATEGORY(collectionDialogLog, "h5reader.collection.dialog")
}

ReaderCollectionDialog::ReaderCollectionDialog(ReaderCollection collection, QWidget* parent,
                                               const QString& bufferRoot, const QString& installedRoot)
    : QDialog(parent) {
    setObjectName(QStringLiteral("ReaderCollectionDialog"));
    setWindowTitle(tr("Trajectories"));
    resize(1040, 560);

    auto* layout = new QVBoxLayout(this);
    source_ = new QLabel(this);
    source_->setWordWrap(true);
    source_->setTextFormat(Qt::PlainText);
    source_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(source_);
    search_ = new QLineEdit(this);
    search_->setObjectName(QStringLiteral("collectionSearch"));
    search_->setPlaceholderText(tr("Search BMRB, protein, organism, or PDB"));
    search_->setClearButtonEnabled(true);
    layout->addWidget(search_);

    table_ = new QTableWidget(0, 6, this);
    table_->setObjectName(QStringLiteral("collectionTable"));
    table_->setHorizontalHeaderLabels({tr("BMRB"), tr("Protein"), tr("Organism"), tr("PDB"), tr("Frames"), tr("Size")});
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->setSelectionMode(QAbstractItemView::SingleSelection);
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_->verticalHeader()->hide();
    table_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    table_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    table_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Stretch);
    table_->horizontalHeader()->setSectionResizeMode(3, QHeaderView::ResizeToContents);
    table_->horizontalHeader()->setSectionResizeMode(4, QHeaderView::ResizeToContents);
    table_->horizontalHeader()->setSectionResizeMode(5, QHeaderView::ResizeToContents);
    layout->addWidget(table_, 1);

    detail_ = new QLabel(this);
    detail_->setWordWrap(true);
    detail_->setTextFormat(Qt::PlainText);
    detail_->setMinimumHeight(40);
    layout->addWidget(detail_);
    progress_ = new QProgressBar(this);
    progress_->setObjectName(QStringLiteral("collectionProgress"));
    progress_->setRange(0, 1000);
    progress_->setTextVisible(true);
    layout->addWidget(progress_);
    status_ = new QLabel(this);
    status_->setObjectName(QStringLiteral("collectionStatus"));
    status_->setWordWrap(true);
    status_->setTextFormat(Qt::PlainText);
    status_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(status_);

    auto* buttons = new QHBoxLayout;
    open_ = new QPushButton(style()->standardIcon(QStyle::SP_DialogOpenButton), tr("Open"), this);
    open_->setObjectName(QStringLiteral("openCollectionRun"));
    open_->setAutoDefault(false);
    cancel_ = new QPushButton(style()->standardIcon(QStyle::SP_DialogCancelButton), tr("Cancel"), this);
    cancel_->setObjectName(QStringLiteral("cancelCollectionTransfer"));
    cancel_->setAutoDefault(false);
    auto* close = new QPushButton(tr("Close"), this);
    close->setAutoDefault(false);
    buttons->addWidget(open_);
    buttons->addWidget(cancel_);
    buttons->addWidget(close);
    buttons->addStretch();
    layout->addLayout(buttons);

    bufferRoot_ = bufferRoot.isEmpty()
        ? QDir(QStandardPaths::writableLocation(QStandardPaths::CacheLocation))
              .filePath(QStringLiteral("reader-buffer"))
        : bufferRoot;
    installedRoot_ = installedRoot.isEmpty()
        ? QDir(QCoreApplication::applicationDirPath()).filePath(QStringLiteral("datasets")) : installedRoot;
    buffer_ = std::make_unique<sciencefiles::BundleCache>(bufferRoot_, this);
    using Cache = sciencefiles::BundleCache;
    connect(buffer_.get(), &Cache::stateChanged, this, [this](Cache::State state) {
        progress_->setRange(0, 0);
        const auto* entry = findEntry(openingKey_);
        const QString title = entry ? entry->title : openingKey_;
        if (state == Cache::State::Downloading)
            status_->setText(tr("Fetching %1...").arg(title));
        else if (state == Cache::State::Extracting)
            status_->setText(tr("Unpacking %1...").arg(title));
        else if (state == Cache::State::Clearing)
            status_->setText(tr("Removing the previous local copy..."));
        else if (state == Cache::State::Cleaning)
            status_->setText(tr("Finishing and removing temporary files..."));
        refresh();
    });
    connect(buffer_.get(), &Cache::progress, this, [this](qint64 completed, qint64 total) {
        progress_->setRange(0, total > 0 ? 1000 : 0);
        if (total > 0)
            progress_->setValue(static_cast<int>(1000.0 * static_cast<double>(completed) / static_cast<double>(total)));
        progress_->setFormat(total > 0
            ? tr("%1 / %2").arg(QLocale().formattedDataSize(completed), QLocale().formattedDataSize(total))
            : QLocale().formattedDataSize(completed));
    });
    connect(buffer_.get(), &Cache::finished, this,
            [this](Cache::Result result, const QString& detail, const QString& directory) {
                const QString key = std::exchange(openingKey_, {});
                if (result == Cache::Result::Ready)
                    retryKey_.clear();
                status_->setText(result == Cache::Result::Ready
                                     ? tr("Ready.")
                                     : result == Cache::Result::Cleared && !pendingLoadError_.isEmpty()
                                         ? std::exchange(pendingLoadError_, {}) : detail);
                refresh();
                if (result == Cache::Result::Ready && !key.isEmpty() && !closing_)
                    emit openRequested(QDir(directory).filePath(QStringLiteral("run.LGS")), key);
                else if (result != Cache::Result::Cleared && result != Cache::Result::Cancelled
                         && result != Cache::Result::Ready)
                    showLoadError(detail);
            });
    connect(buffer_.get(), &Cache::shutdownFinished, this,
            &ReaderCollectionDialog::finishShutdown, Qt::QueuedConnection);

    connect(search_, &QLineEdit::textChanged, this, &ReaderCollectionDialog::filterRows);
    connect(table_, &QTableWidget::itemSelectionChanged, this, &ReaderCollectionDialog::refresh);
    connect(table_, &QTableWidget::cellDoubleClicked, this, &ReaderCollectionDialog::openSelected);
    connect(open_, &QPushButton::clicked, this, &ReaderCollectionDialog::openSelected);
    connect(cancel_, &QPushButton::clicked, this, &ReaderCollectionDialog::cancel);
    connect(close, &QPushButton::clicked, this, &QDialog::close);
    setCollection(std::move(collection));
    search_->setFocus();
    refresh();
}

void ReaderCollectionDialog::refresh() {
    const int row = table_->currentRow();
    const bool selected = row >= 0 && row < collection_.entries.size() && !table_->isRowHidden(row);
    open_->setEnabled((selected || !retryCatalog_.isEmpty()) && !isBusy() && !closing_);
    if (selected) {
        const auto& entry = collection_.entries[row];
        const auto bundle = bundleFor(entry);
        const bool ready = !installedPath(entry).isEmpty() || !buffer_->cachedPath(bundle).isEmpty();
        open_->setText(ready ? tr("Open") : entry.key == retryKey_ ? tr("Retry")
                       : entry.archiveUrl.isLocalFile() ? tr("Copy and open") : tr("Download and open"));
    } else if (!retryCatalog_.isEmpty()) {
        open_->setText(tr("Retry collection"));
    }
    cancel_->setEnabled(!closing_ && (catalogReply_ || buffer_->state() == sciencefiles::BundleCache::State::Downloading
                                     || buffer_->state() == sciencefiles::BundleCache::State::Extracting));
    progress_->setVisible(isBusy());
    detail_->setText(selected ? collection_.entries[row].description : QString());
}

void ReaderCollectionDialog::openSelected() {
    const int row = table_->currentRow();
    if (row < 0 || table_->isRowHidden(row)) {
        if (!retryCatalog_.isEmpty())
            loadCatalog(retryCatalog_);
        return;
    }
    QString error;
    if (!openTrajectory(collection_.entries[row].key, &error))
        showLoadError(error);
}

bool ReaderCollectionDialog::openTrajectory(const QString& key, QString* error) {
    const auto* entry = findEntry(key);
    if (!entry || isBusy() || closing_) {
        *error = tr("Unknown trajectory, or the collection is busy.");
        return false;
    }
    lastError_.clear();
    const QString installed = installedPath(*entry);
    if (!installed.isEmpty()) {
        emit openRequested(QDir(installed).filePath(entry->entryPoint), entry->key);
        return true;
    }
    const auto bundle = bundleFor(*entry);
    if (buffer_->cachedPath(bundle).isEmpty() && entry->archiveUrl.isLocalFile()) {
        const QFileInfo archive(entry->archiveUrl.toLocalFile());
        if (!archive.isFile() || archive.size() != entry->archiveBytes) {
            *error = tr("Package missing or incomplete: %1").arg(archive.filePath());
            return false;
        }
    }
    openingKey_ = bundle.key;
    retryKey_ = bundle.key;
    pendingLoadError_.clear();
    qCInfo(collectionDialogLog).noquote() << "Opening collection package" << entry->archiveUrl.toDisplayString();
    if (!buffer_->fetch(bundle)) {
        openingKey_.clear();
        *error = tr("The local buffer is busy or shutting down.");
        return false;
    }
    return true;
}

sciencefiles::Bundle ReaderCollectionDialog::bundleFor(const ReaderCollection::Entry& entry) const {
    return {entry.key, entry.archiveUrl, entry.entryPoint, entry.expandedBytes};
}

void ReaderCollectionDialog::runOpened(const QString& key) {
    const QString previous = std::exchange(currentKey_, key);
    if (!previous.isEmpty() && previous.compare(key, Qt::CaseInsensitive) != 0 && !buffer_->clear(previous))
        qCWarning(collectionDialogLog).noquote() << "Could not remove previous local copy" << previous;
}

void ReaderCollectionDialog::runFailed(const QString& key, const QString& message) {
    showLoadError(message);
    if (key.compare(currentKey_, Qt::CaseInsensitive) == 0)
        return;
    pendingLoadError_ = message;
    if (!buffer_->clear(key)) {
        pendingLoadError_.clear();
    }
}

void ReaderCollectionDialog::cancelPendingOpen() {
    if (openingKey_.isEmpty())
        return;
    openingKey_.clear();
    buffer_->cancel();
}

void ReaderCollectionDialog::shutdown() {
    closing_ = true;
    cancel();
    buffer_->shutdown();
    refresh();
}

bool ReaderCollectionDialog::isShutdown() const {
    return shutdownComplete_;
}

bool ReaderCollectionDialog::isBusy() const {
    return catalogReply_ || buffer_->busy();
}

void ReaderCollectionDialog::reject() {
    cancel();
    QDialog::reject();
}

void ReaderCollectionDialog::showLoadError(const QString& message) {
    status_->setText(message);
    lastError_ = message;
    qCWarning(collectionDialogLog).noquote() << message;
}

bool ReaderCollectionDialog::setCollection(ReaderCollection collection, QString* error) {
    if (isBusy() || closing_) {
        if (error)
            *error = tr("Finish or cancel the current transfer before changing collections.");
        return false;
    }
    table_->setRowCount(0);
    collection_ = std::move(collection);
    retryCatalog_.clear();
    retryKey_.clear();
    lastError_.clear();
    setWindowTitle(collection_.title.isEmpty() ? tr("Trajectories") : collection_.title);
    source_->setText(collection_.source.toDisplayString(QUrl::PreferLocalFile));
    table_->setRowCount(static_cast<int>(collection_.entries.size()));
    for (int row = 0; row < collection_.entries.size(); ++row) {
        const auto& entry = collection_.entries[row];
        const QStringList cells{entry.bmrb, entry.title, entry.organism, entry.pdb, QString::number(entry.frames),
                                QLocale().formattedDataSize(entry.archiveBytes)};
        for (int column = 0; column < cells.size(); ++column) {
            auto* item = new QTableWidgetItem(cells[column]);
            item->setToolTip(cells[column]);
            table_->setItem(row, column, item);
        }
    }
    search_->clear();
    status_->setText(tr("%1 trajectories").arg(collection_.entries.size()));
    filterRows();
    return true;
}

bool ReaderCollectionDialog::loadCatalog(const QUrl& source, QString* error) {
    const auto fail = [this, error](const QString& message) {
        if (error)
            *error = message;
        showLoadError(message);
        return false;
    };
    if (isBusy() || closing_)
        return fail(tr("Finish or cancel the current transfer before changing collections."));
    if (source.isLocalFile()) {
        QString detail;
        auto document = ReadReaderDocument(source.toLocalFile(), &detail);
        if (!document)
            return fail(detail);
        if (!document->collection)
            return fail(tr("This LGS describes one run, not a collection."));
        return setCollection(std::move(*document->collection), error);
    }
    if (!source.isValid() || source.scheme() != QStringLiteral("https") || source.host().isEmpty())
        return fail(tr("Choose a local collection LGS or an HTTPS collection URL."));

    lastError_.clear();
    QNetworkRequest request(source);
    request.setTransferTimeout(std::chrono::seconds(60));
    catalogReply_ = network_.get(request);
    retryCatalog_ = source;
    status_->setText(tr("Reading collection: %1").arg(source.toDisplayString()));
    progress_->setRange(0, 0);
    qCInfo(collectionDialogLog).noquote() << "Fetching collection" << source.toDisplayString();
    connect(catalogReply_, &QNetworkReply::finished, this, [this] {
        auto* reply = catalogReply_.data();
        catalogReply_.clear();
        const auto networkError = reply->error();
        const QString networkDetail = reply->errorString();
        const auto source = reply->url();
        const auto bytes = networkError == QNetworkReply::NoError ? reply->readAll() : QByteArray();
        reply->deleteLater();
        if (closing_ || networkError == QNetworkReply::OperationCanceledError) {
            retryCatalog_.clear();
            status_->setText(tr("Cancelled."));
        } else if (networkError != QNetworkReply::NoError) {
            showLoadError(tr("Could not read collection: %1").arg(networkDetail));
        } else {
            QJsonParseError parseError;
            const auto json = QJsonDocument::fromJson(bytes, &parseError);
            QString error;
            if (parseError.error != QJsonParseError::NoError || !json.isObject()) {
                showLoadError(tr("Invalid collection JSON: %1").arg(parseError.error == QJsonParseError::NoError
                    ? tr("expected an object") : parseError.errorString()));
            } else if (auto collection = ParseReaderCollection(json.object(), source, &error)) {
                setCollection(std::move(*collection));
            } else {
                showLoadError(error);
            }
        }
        refresh();
        finishShutdown();
    });
    refresh();
    return true;
}

void ReaderCollectionDialog::filterRows() {
    const auto words = search_->text().simplified().split(QLatin1Char(' '), Qt::SkipEmptyParts);
    int firstVisible = -1;
    for (int row = 0; row < collection_.entries.size(); ++row) {
        const auto& entry = collection_.entries[row];
        const QString text = QStringList{entry.key, entry.bmrb, entry.title, entry.group, entry.description,
                                        entry.organism, entry.pdb, entry.keywords.join(QLatin1Char(' '))}
                                 .join(QLatin1Char(' '));
        bool match = true;
        for (const auto& word : words)
            match = match && text.contains(word, Qt::CaseInsensitive);
        table_->setRowHidden(row, !match);
        if (match && firstVisible < 0)
            firstVisible = row;
    }
    if (firstVisible >= 0 && (table_->currentRow() < 0 || table_->isRowHidden(table_->currentRow())))
        table_->selectRow(firstVisible);
    else if (firstVisible < 0)
        table_->clearSelection();
    refresh();
}

const ReaderCollection::Entry* ReaderCollectionDialog::findEntry(const QString& key) const {
    for (const auto& entry : collection_.entries)
        if (entry.key == key)
            return &entry;
    return nullptr;
}

QString ReaderCollectionDialog::installedPath(const ReaderCollection::Entry& entry) const {
    const QString directory = QDir(installedRoot_).filePath(entry.key);
    return QFileInfo::exists(QDir(directory).filePath(entry.entryPoint)) ? directory : QString();
}

void ReaderCollectionDialog::setCurrentRun(const QString& path) {
    currentRun_ = path;
    refresh();
}

bool ReaderCollectionDialog::isClearing() const {
    return buffer_->state() == sciencefiles::BundleCache::State::Clearing;
}

bool ReaderCollectionDialog::clearTrajectory(const QString& key, QString* error) {
    const auto* entry = findEntry(key);
    if (!entry || isBusy() || closing_) {
        *error = tr("Unknown trajectory, or the collection is busy.");
        return false;
    }
    if (!installedPath(*entry).isEmpty()) {
        *error = tr("This example belongs to the installation.");
        return false;
    }
    const QString relative = QDir(QDir(bufferRoot_).filePath(key)).relativeFilePath(currentRun_);
    if (!currentRun_.isEmpty() && relative != QStringLiteral("..")
        && !relative.startsWith(QStringLiteral("../")) && !QDir::isAbsolutePath(relative)) {
        *error = tr("Open a different trajectory before clearing this local copy.");
        return false;
    }
    return buffer_->clear(key);
}

void ReaderCollectionDialog::cancel() {
    cancelPendingOpen();
    if (catalogReply_)
        catalogReply_->abort();
}

void ReaderCollectionDialog::finishShutdown() {
    if (!closing_ || shutdownComplete_ || catalogReply_ || !buffer_->isShutdown())
        return;
    shutdownComplete_ = true;
    emit shutdownFinished();
}

QJsonObject ReaderCollectionDialog::state() const {
    QJsonArray entries;
    for (const auto& entry : collection_.entries) {
        const QString installed = installedPath(entry);
        const QString cached = buffer_->cachedPath(bundleFor(entry));
        entries.append(QJsonObject{{"key", entry.key}, {"title", entry.title}, {"description", entry.description},
                                   {"organism", entry.organism}, {"bmrb", entry.bmrb}, {"pdb", entry.pdb},
                                   {"keywords", QJsonArray::fromStringList(entry.keywords)},
                                   {"frames", entry.frames}, {"archive_bytes", entry.archiveBytes},
                                   {"expanded_bytes", entry.expandedBytes}, {"url", entry.archiveUrl.toString()},
                                   {"included", !installed.isEmpty()}, {"downloaded", !cached.isEmpty()},
                                   {"directory", installed.isEmpty() ? cached : installed},
                                   {"entry_point", entry.entryPoint}});
    }
    return {{"entries", entries}, {"source", collection_.source.toString()}, {"busy", isBusy()},
            {"phase", catalogReply_ ? QStringLiteral("catalog")
                : QString::fromLatin1(QMetaEnum::fromType<sciencefiles::BundleCache::State>().valueToKey(static_cast<int>(buffer_->state())))},
            {"cache_root", bufferRoot_}, {"shutdown", isShutdown()}, {"error", lastError_},
            {"visible", isVisible()}};
}

}  // namespace h5reader::app

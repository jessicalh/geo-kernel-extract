#include "ArchiveTask.h"
#include "BundlePaths.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QLoggingCategory>
#include <QSet>
#include <QStorageInfo>
#include <archive.h>
#include <archive_entry.h>
#include <array>

namespace sciencefiles {
Q_LOGGING_CATEGORY(archiveLog, "science.archive")
Q_LOGGING_CATEGORY(archiveFilesLog, "science.archive.files", QtInfoMsg)

namespace {
constexpr size_t CopyBufferBytes = 64 * 1024;
constexpr qint64 ProgressIntervalBytes = 1024 * 1024;

struct Archive {
    archive *handle = archive_read_new();
    ~Archive() {
        if (handle && archive_read_free(handle) != ARCHIVE_OK)
            qCWarning(archiveLog) << "Archive cleanup failed after extraction was abandoned.";
    }
    QString error() const { return QString::fromUtf8(archive_error_string(handle)); }
};
} // namespace

ArchiveTask::ArchiveTask(std::unique_ptr<QTemporaryDir> staging, const Bundle &bundle, const QString &destination,
                         bool extract, BundleCache::Result result, const QString &detail)
    : staging_(std::move(staging)), bundle_(bundle), destination_(destination), extract_(extract), result_(result),
      detail_(detail) {}

void ArchiveTask::run() {
    if (extract_) {
        const auto extraction = extractArchive();
        result_ = extraction.result;
        detail_ = extraction.detail;
        if (isInterruptionRequested()) {
            result_ = BundleCache::Result::Cancelled;
            detail_ = tr("Cancelled. No incomplete bundle was kept.");
        } else if (result_ == BundleCache::Result::Ready) {
            if (!QDir().rename(staging_->filePath("content"), destination_)) {
                result_ = BundleCache::Result::FileError;
                detail_ = tr("Could not publish the completed bundle: %1").arg(destination_);
            } else {
                // This rename is the commit point. A later cancel cannot undo a complete bundle.
                publishedPath_ = destination_;
                detail_ = tr("Bundle ready (%1 expanded bytes).").arg(bundle_.expandedBytes);
            }
        }
    }
    if (extract_)
        emit cleaning();
    if (!staging_->remove()) {
        result_ = BundleCache::Result::FileError;
        detail_ += tr(" Temporary directory could not be removed: %1").arg(staging_->path());
    }
    qCInfo(archiveLog) << "Disk operation finished:" << bundle_.key << "code=" << static_cast<int>(result_) << detail_;
}

ArchiveTask::ExtractionResult ArchiveTask::extractArchive() {
    using Result = BundleCache::Result;
    const QStorageInfo volume(staging_->path());
    if (!volume.isValid() || !volume.isReady() || volume.bytesAvailable() < bundle_.expandedBytes)
        return {Result::FileError, tr("Not enough available space for %1 expanded bytes on %2.")
                                       .arg(bundle_.expandedBytes)
                                       .arg(volume.rootPath())};
    const QString archivePath = staging_->filePath("download.archive");
    Archive input;
    if (!input.handle)
        return {Result::ArchiveError, tr("Cannot allocate the archive reader.")};
    // Read past tar end markers so the xz checksum/footer is checked in this same pass.
    // Require native xz support: ARCHIVE_WARN would select an external decompressor.
    if (archive_read_support_format_zip_seekable(input.handle) != ARCHIVE_OK ||
        archive_read_support_format_tar(input.handle) != ARCHIVE_OK ||
        archive_read_support_filter_xz(input.handle) != ARCHIVE_OK ||
        archive_read_set_options(input.handle, "hdrcharset=UTF-8,tar:read_concatenated_archives") != ARCHIVE_OK)
        return {Result::ArchiveError, tr("Cannot configure the archive reader: %1").arg(input.error())};
#ifdef Q_OS_WIN
    const int opened = archive_read_open_filename_w(
        input.handle, reinterpret_cast<const wchar_t *>(archivePath.utf16()), CopyBufferBytes);
#else
    const int opened =
        archive_read_open_filename(input.handle, QFile::encodeName(archivePath).constData(), CopyBufferBytes);
#endif
    if (opened != ARCHIVE_OK)
        return {Result::ArchiveError, tr("Cannot open the downloaded archive: %1").arg(input.error())};
    const QString content = staging_->filePath("content");
    if (!QDir().mkdir(content))
        return {Result::FileError, tr("Cannot create the extraction directory: %1").arg(content)};

    qint64 expanded = 0;
    qint64 declared = 0;
    qint64 lastProgress = 0;
    QSet<QString> names;
    QSet<QString> directories{content};
    for (;;) {
        if (isInterruptionRequested())
            return {Result::Cancelled, {}};
        archive_entry *entry = nullptr;
        const int next = archive_read_next_header(input.handle, &entry);
        if (next == ARCHIVE_EOF)
            break;
        if (next != ARCHIVE_OK)
            return {Result::ArchiveError, tr("Cannot read an archive entry: %1").arg(input.error())};
        const QByteArray bytes(archive_entry_pathname_utf8(entry));
        QString name = QString::fromUtf8(bytes);
        const auto type = archive_entry_filetype(entry);
        const bool directory = type == AE_IFDIR;
        if (directory && name.endsWith('/'))
            name.chop(1);
        if (QString::fromUtf8(bytes).toUtf8() != bytes || !safeRelativePath(name))
            return {Result::ArchiveError, tr("Archive contains an unsafe or non-UTF-8 path: %1").arg(name)};
        if (archive_entry_is_encrypted(entry) || archive_entry_hardlink(entry) || archive_entry_hardlink_w(entry) ||
            (type != AE_IFREG && !directory))
            return {Result::ArchiveError, tr("Archive contains an encrypted or non-regular entry: %1").arg(name)};
#ifdef Q_OS_WIN
        const QString identity = name.toCaseFolded();
#else
        const QString identity = name;
#endif
        if (names.contains(identity))
            return {Result::ArchiveError, tr("Archive contains a duplicate path: %1").arg(name)};
        names.insert(identity);
        const qint64 size = archive_entry_size(entry);
        if (!archive_entry_size_is_set(entry) || size < 0 || size > bundle_.expandedBytes - declared)
            return {Result::ArchiveError, tr("Archive size is invalid or exceeds the catalogue's byte count.")};
        declared += size;
        const QString destination = QDir(content).filePath(name);
        const QString parent = directory ? destination : QFileInfo(destination).absolutePath();
        if (!directories.contains(parent)) {
            if (!QDir().mkpath(parent))
                return {Result::FileError, tr("Cannot create archive directory: %1").arg(parent)};
            directories.insert(parent);
        }
        if (directory) {
            if (size != 0)
                return {Result::ArchiveError, tr("Archive directory declares file data: %1").arg(name)};
            char byte;
            if (archive_read_data(input.handle, &byte, 1) != 0)
                return {Result::ArchiveError,
                        tr("Archive directory integrity check failed for %1: %2").arg(name, input.error())};
        } else {
            QFile output(destination);
            if (!output.open(QIODevice::WriteOnly | QIODevice::NewOnly)) {
                const auto result = QFileInfo::exists(destination) ? Result::ArchiveError : Result::FileError;
                return {result, tr("Cannot create %1: %2").arg(name, output.errorString())};
            }
            std::array<char, CopyBufferBytes> buffer;
            qint64 fileBytes = 0;
            for (;;) {
                if (isInterruptionRequested())
                    return {Result::Cancelled, {}};
                const auto count = archive_read_data(input.handle, buffer.data(), buffer.size());
                if (count < 0)
                    return {Result::ArchiveError,
                            tr("Archive decompression failed for %1: %2").arg(name, input.error())};
                if (count == 0)
                    break;
                if (count > size - fileBytes || count > bundle_.expandedBytes - expanded)
                    return {Result::ArchiveError, tr("Archive output exceeds the declared byte count.")};
                if (output.write(buffer.data(), count) != count)
                    return {Result::FileError, tr("Cannot write %1: %2").arg(name, output.errorString())};
                expanded += count;
                fileBytes += count;
                if (expanded - lastProgress >= ProgressIntervalBytes) {
                    emit progress(expanded, bundle_.expandedBytes);
                    lastProgress = expanded;
                }
            }
            if (fileBytes != size)
                return {Result::ArchiveError, tr("Archive byte count does not match for %1.").arg(name)};
            if (!output.flush())
                return {Result::FileError, tr("Cannot flush %1: %2").arg(name, output.errorString())};
            output.close();
            if (output.error() != QFileDevice::NoError)
                return {Result::FileError, tr("Cannot close %1: %2").arg(name, output.errorString())};
            qCDebug(archiveFilesLog) << "Extracted:" << name << "bytes=" << fileBytes;
        }
    }
    if (expanded != bundle_.expandedBytes || declared != expanded)
        return {Result::ArchiveError, tr("Archive expanded byte count does not match the catalogue.")};
    if (!QFileInfo(QDir(content).filePath(bundle_.entryPoint)).isFile())
        return {Result::ArchiveError,
                tr("Archive does not contain the expected entry file: %1").arg(bundle_.entryPoint)};
    if (archive_read_close(input.handle) != ARCHIVE_OK)
        return {Result::ArchiveError, tr("Cannot close the archive: %1").arg(input.error())};
    const int freed = archive_read_free(input.handle);
    input.handle = nullptr;
    if (freed != ARCHIVE_OK)
        return {Result::ArchiveError, tr("Cannot release the archive reader (code %1).").arg(freed)};
    emit progress(expanded, bundle_.expandedBytes);
    return {Result::Ready, {}};
}

} // namespace sciencefiles

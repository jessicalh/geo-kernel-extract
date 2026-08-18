#pragma once

#include "../model/ScientificAlignment.h"

#include <QJsonObject>
#include <QString>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace h5reader::io {

struct QtLoadResult;

struct ReaderAlignmentExportRequest {
    QString outputRoot;
    QString datasetId;
    QString runId;
    QString humanName;
    QString lgsPath;
    QString calcsetRoot;
    QString extractionDirectory;
    QString trajectoryH5;
    QString extractionManifest;
    QString readerVersion;
    qint64 processId = 0;

    std::size_t residueCount = 0;
    std::size_t bondCount = 0;
    std::size_t aromaticRingCount = 0;
    std::size_t saturatedRingCount = 0;
    std::size_t ringCount = 0;
    std::size_t ringMembershipCount = 0;

    model::ScientificPositionTable positions;
    std::vector<std::uint64_t> originalFrameIndices;
    std::vector<double> frameTimesPicoseconds;
    std::vector<std::size_t> primaryFitAtoms;
    std::vector<std::size_t> caFitAtoms;
};

struct ReaderAlignmentExportPreparation {
    bool ok = false;
    QString error;
    ReaderAlignmentExportRequest request;
};

class ReaderAlignmentExportControl final {
public:
    void requestCancel() noexcept {
        cancelRequested_.store(true, std::memory_order_release);
    }
    bool cancelRequested() const noexcept {
        return cancelRequested_.load(std::memory_order_acquire);
    }

private:
    std::atomic_bool cancelRequested_{false};
};

struct CompletedAlignmentValidation {
    bool ok = false;
    QString error;
    QJsonObject manifest;
};

struct ReaderAlignmentExportResult {
    bool ok = false;
    bool cancelled = false;
    bool alreadyComplete = false;
    QString error;
    QString memberId;
    QString finalDirectory;
    QJsonObject manifest;
    model::ScientificAlignmentResult primaryAlignment;
    model::ScientificAlignmentResult caAlignment;
};

class ReaderAlignmentExporter final {
public:
    // Capture the small, immutable part of the loaded Reader state needed by
    // the exporter. This must run on the Reader's GUI thread; the resulting
    // request owns its data and is safe to pass to QtConcurrent.
    static ReaderAlignmentExportPreparation Prepare(
        const QtLoadResult& loaded,
        const QString& outputRoot);

    // File validation, fitting, hashing, and publication. This overload does
    // not touch Reader or QObject state and may run on a Qt worker thread.
    static ReaderAlignmentExportResult Export(
        const ReaderAlignmentExportRequest& request,
        const ReaderAlignmentExportControl* control = nullptr);

    static CompletedAlignmentValidation ValidateCompletedMember(
        const QString& memberDirectory,
        const ReaderAlignmentExportControl* control = nullptr);

    static model::ScientificAlignmentResult LoadCompletedPrimaryAlignment(
        const QString& memberDirectory,
        const ReaderAlignmentExportControl* control = nullptr);

    static bool RebuildMembersTable(const QString& outputRoot,
                                    QString* error = nullptr,
                                    const ReaderAlignmentExportControl* control = nullptr);
};

}  // namespace h5reader::io

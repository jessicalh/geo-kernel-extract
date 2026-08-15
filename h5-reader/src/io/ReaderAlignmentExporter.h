#pragma once

#include "../model/ScientificAlignment.h"

#include <QJsonObject>
#include <QString>

namespace h5reader::io {

struct QtLoadResult;

struct CompletedAlignmentValidation {
    bool ok = false;
    QString error;
    QJsonObject manifest;
};

struct ReaderAlignmentExportResult {
    bool ok = false;
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
    static ReaderAlignmentExportResult Export(const QtLoadResult& loaded,
                                               const QString& outputRoot);

    static CompletedAlignmentValidation ValidateCompletedMember(
        const QString& memberDirectory);

    static model::ScientificAlignmentResult LoadCompletedPrimaryAlignment(
        const QString& memberDirectory);

    static bool RebuildMembersTable(const QString& outputRoot,
                                    QString* error = nullptr);
};

}  // namespace h5reader::io

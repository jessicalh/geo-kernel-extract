#include "RowDesignSink.h"

#include <QDir>
#include <QIODevice>

namespace h5reader::rediscover {

namespace {

QString csvEscape(QString v) {
    if (v.contains('"')) v.replace(QStringLiteral("\""), QStringLiteral("\"\""));
    if (v.contains(',') || v.contains('\n') || v.contains('"'))
        return QStringLiteral("\"%1\"").arg(v);
    return v;
}

}  // namespace

RowDesignSink::RowDesignSink(const QString& outDir, std::size_t embeddingDims)
    : outDir_(outDir), embeddingDims_(embeddingDims) {
    QDir().mkpath(outDir_);
    rowsPath_ = QDir(outDir_).filePath(QStringLiteral("row_design_rows.csv"));

    rowsFile_ = std::make_unique<QSaveFile>(rowsPath_);
    if (!rowsFile_->open(QIODevice::WriteOnly | QIODevice::Text)) return;
    rowsOut_ = std::make_unique<QTextStream>(rowsFile_.get());
    *rowsOut_ << RowDesignHeader().join(QStringLiteral(",")) << '\n';
    ok_ = true;
}

RowDesignSink::~RowDesignSink() = default;

bool RowDesignSink::WriteRow(const RowDesignRow& row) {
    if (!ok_ || row.values.size() != RowDesignSchema().size()) return false;
    for (std::size_t i = 0; i < row.values.size(); ++i) {
        if (i) *rowsOut_ << ',';
        *rowsOut_ << csvEscape(row.values[i]);
    }
    *rowsOut_ << '\n';
    ++rows_;
    return true;
}

bool RowDesignSink::Commit() {
    if (!ok_ || committed_) return false;
    rowsOut_->flush();
    if (!rowsFile_->commit()) return false;
    committed_ = true;
    return committed_;
}

}  // namespace h5reader::rediscover

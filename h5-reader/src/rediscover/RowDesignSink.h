#pragma once

#include "RowDesign.h"

#include <QSaveFile>
#include <QString>
#include <QStringList>
#include <QTextStream>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace h5reader::rediscover {

class RowDesignSink {
public:
    explicit RowDesignSink(const QString& outDir, std::size_t embeddingDims = 256);
    ~RowDesignSink();

    bool Ok() const { return ok_; }
    bool WriteRow(const RowDesignRow& row);
    bool Commit();

    std::size_t rowsWritten() const { return rows_; }
    const QStringList& sidecarFiles() const { return sidecarFiles_; }

private:
    QString outDir_;
    QString rowsPath_;
    QStringList sidecarFiles_;
    std::unique_ptr<QSaveFile> rowsFile_;
    std::unique_ptr<QTextStream> rowsOut_;
    bool ok_ = false;
    bool committed_ = false;
    std::size_t rows_ = 0;
    std::size_t embeddingDims_ = 256;
    std::vector<double> targetT2_;
    std::vector<double> ringTensors_;
    std::vector<float> embeddings_;
};

}  // namespace h5reader::rediscover

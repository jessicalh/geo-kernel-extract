#include "RowDesignSink.h"

#include <QByteArray>
#include <QDir>
#include <QIODevice>

#include <cmath>
#include <limits>

namespace h5reader::rediscover {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const float kNaNF32 = std::numeric_limits<float>::quiet_NaN();

QString csvEscape(QString v) {
    if (v.contains('"')) v.replace(QStringLiteral("\""), QStringLiteral("\"\""));
    if (v.contains(',') || v.contains('\n') || v.contains('"'))
        return QStringLiteral("\"%1\"").arg(v);
    return v;
}

template <typename T>
bool writeNpy(const QString& path, const std::vector<std::size_t>& shape,
              const std::vector<T>& data, const QByteArray& descr) {
    std::size_t n = 1;
    for (std::size_t dim : shape) n *= dim;
    if (n != data.size()) return false;

    QByteArray header;
    header += "{'descr': '";
    header += descr;
    header += "', 'fortran_order': False, 'shape': (";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        if (i) header += ", ";
        header += QByteArray::number(static_cast<qulonglong>(shape[i]));
    }
    if (shape.size() == 1) header += ",";
    header += "), }";

    constexpr int kPreambleBytes = 10;
    const int newlineBytes = 1;
    int pad = 16 - ((kPreambleBytes + header.size() + newlineBytes) % 16);
    if (pad == 16) pad = 0;
    header += QByteArray(pad, ' ');
    header += '\n';
    if (header.size() > 65535) return false;

    QByteArray prefix;
    prefix.append("\x93NUMPY", 6);
    prefix.append(char(1));
    prefix.append(char(0));
    const quint16 headerLen = static_cast<quint16>(header.size());
    prefix.append(char(headerLen & 0xff));
    prefix.append(char((headerLen >> 8) & 0xff));

    QSaveFile f(path);
    if (!f.open(QIODevice::WriteOnly)) return false;
    if (f.write(prefix) != prefix.size()) return false;
    if (f.write(header) != header.size()) return false;
    const qsizetype payloadBytes = static_cast<qsizetype>(data.size() * sizeof(T));
    if (payloadBytes > 0
        && f.write(reinterpret_cast<const char*>(data.data()), payloadBytes) != payloadBytes)
        return false;
    return f.commit();
}

}  // namespace

RowDesignSink::RowDesignSink(const QString& outDir, std::size_t embeddingDims)
    : outDir_(outDir), embeddingDims_(embeddingDims) {
    QDir().mkpath(outDir_);
    rowsPath_ = QDir(outDir_).filePath(QStringLiteral("row_design_rows.csv"));
    sidecarFiles_ = {
        QStringLiteral("row_design_target_T2.npy"),
        QStringLiteral("row_design_ring_tensors.npy"),
        QStringLiteral("row_design_aimnet2_embedding.npy"),
    };

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

    for (double v : row.targetT2) targetT2_.push_back(row.targetT2Present ? v : kNaN);
    for (double v : row.ringTensors) ringTensors_.push_back(row.ringTensorsPresent ? v : kNaN);
    const bool embOk = row.embedding && row.embeddingDims == embeddingDims_;
    for (std::size_t i = 0; i < embeddingDims_; ++i)
        embeddings_.push_back(embOk ? row.embedding[i] : kNaNF32);
    ++rows_;
    return true;
}

bool RowDesignSink::Commit() {
    if (!ok_ || committed_) return false;
    rowsOut_->flush();
    if (!rowsFile_->commit()) return false;
    const bool t2Ok = writeNpy<double>(QDir(outDir_).filePath(QStringLiteral("row_design_target_T2.npy")),
                                       {rows_, 5}, targetT2_, QByteArrayLiteral("<f8"));
    const bool ringOk = writeNpy<double>(QDir(outDir_).filePath(QStringLiteral("row_design_ring_tensors.npy")),
                                         {rows_, 162}, ringTensors_, QByteArrayLiteral("<f8"));
    const bool embOk = writeNpy<float>(QDir(outDir_).filePath(QStringLiteral("row_design_aimnet2_embedding.npy")),
                                       {rows_, embeddingDims_}, embeddings_, QByteArrayLiteral("<f4"));
    committed_ = t2Ok && ringOk && embOk;
    return committed_;
}

}  // namespace h5reader::rediscover

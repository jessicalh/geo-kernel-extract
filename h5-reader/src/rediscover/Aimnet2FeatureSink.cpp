#include "Aimnet2FeatureSink.h"

#include <QByteArray>
#include <QDir>
#include <QFileInfo>
#include <QIODevice>
#include <QLoggingCategory>

#include <cmath>
#include <limits>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cAimnet2Sink, "h5reader.rediscover.aimnet2_sink")

QString num(double v) { return QString::number(v, 'g', 12); }

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const float kNaNF32 = std::numeric_limits<float>::quiet_NaN();

void appendVec3(std::vector<double>& dst, const Vec3& v, bool present) {
    dst.push_back(present ? v.x() : kNaN);
    dst.push_back(present ? v.y() : kNaN);
    dst.push_back(present ? v.z() : kNaN);
}

void appendT2(std::vector<double>& dst, const std::array<double, 5>& t2, bool present) {
    for (double v : t2) dst.push_back(present ? v : kNaN);
}

template <typename T>
bool writeNpy(const QString& path, const std::vector<std::size_t>& shape,
              const std::vector<T>& data, const QByteArray& descr) {
    std::size_t n = 1;
    for (const std::size_t dim : shape) n *= dim;
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

bool finiteT2(const std::array<double, 5>& t2) {
    for (const double v : t2)
        if (!std::isfinite(v)) return false;
    return true;
}

const char* kHeader =
    "row_id,atom_index,residue_index,residue_number,amino_acid_ord,element_ord,"
    "atom_name,backbone_frame_class,h5_row,original_index,time_ps,"
    "frame_z_x,frame_z_y,frame_z_z,frame_x_x,frame_x_y,frame_x_z,"
    "frame_y_x,frame_y_y,frame_y_z,frame_variant,frame_valid,frame_anchor_atom_index,"
    "dft_present,dft_sigma_iso,"
    "dft_total_raw_00,dft_total_raw_01,dft_total_raw_02,"
    "dft_total_raw_10,dft_total_raw_11,dft_total_raw_12,"
    "dft_total_raw_20,dft_total_raw_21,dft_total_raw_22,"
    "dft_total_T2_0,dft_total_T2_1,dft_total_T2_2,dft_total_T2_3,dft_total_T2_4,"
    "dft_local_frame_valid,"
    "aimnet2_charge_present,aimnet2_charge,charge_source,"
    "aimnet2_charge_response_gradient_present,"
    "aimnet2_charge_response_gradient_scalar,"
    "aimnet2_charge_response_gradient_local_x,"
    "aimnet2_charge_response_gradient_local_y,"
    "aimnet2_charge_response_gradient_local_z,"
    "aimnet2_charge_response_gradient_lab_x,"
    "aimnet2_charge_response_gradient_lab_y,"
    "aimnet2_charge_response_gradient_lab_z,"
    "aimnet2_charge_response_gradient_label,"
    "aimnet2_embedding_present,aimnet2_embedding_dim,aimnet2_embedding_role";

}  // namespace

bool FiniteAimnetVec3(const Vec3& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

Aimnet2FeatureSink::Aimnet2FeatureSink(const QString& outDir, const QString& caseName,
                                       std::size_t embeddingDims)
    : caseName_(caseName), embeddingDims_(embeddingDims) {
    QDir().mkpath(outDir);
    aggregatedPath_ = QStringLiteral("%1/%2_aggregated.csv").arg(outDir, caseName);
    sidecarFiles_ = {
        QStringLiteral("%1_charge_response_gradient_local.npy").arg(caseName),
        QStringLiteral("%1_charge_response_gradient_lab.npy").arg(caseName),
        QStringLiteral("%1_charge_response_gradient_scalar.npy").arg(caseName),
        QStringLiteral("%1_embedding.npy").arg(caseName),
        QStringLiteral("%1_target_local_T2.npy").arg(caseName),
        QStringLiteral("%1_target_lab_T2.npy").arg(caseName),
    };

    aggregatedFile_ = std::make_unique<QSaveFile>(aggregatedPath_);
    if (!aggregatedFile_->open(QIODevice::WriteOnly | QIODevice::Text)) {
        qCWarning(cAimnet2Sink).noquote() << "cannot open AIMNet2 aggregated CSV"
                                          << aggregatedPath_;
        return;
    }
    aggregatedOut_ = std::make_unique<QTextStream>(aggregatedFile_.get());
    *aggregatedOut_ << kHeader << '\n';
    ok_ = true;
}

Aimnet2FeatureSink::~Aimnet2FeatureSink() = default;

void Aimnet2FeatureSink::Write(const Aimnet2FeatureRow& row) {
    if (!ok_) return;
    const int64_t rowId = nextRowId_++;
    const bool crgPresent =
        row.crg_present && std::isfinite(row.crg_scalar)
        && FiniteAimnetVec3(row.crg_vector_lab)
        && (!row.frame_valid || FiniteAimnetVec3(row.crg_vector_local));
    const bool localCrgPresent = crgPresent && row.frame_valid && FiniteAimnetVec3(row.crg_vector_local);
    const bool embeddingPresent =
        row.embedding_present && row.embedding && row.embedding_dims == embeddingDims_;
    const bool targetLocalPresent =
        row.dft_present && row.dft_local_frame_valid && finiteT2(row.dft_target_local_T2);

    QTextStream& out = *aggregatedOut_;
    out << rowId << ',' << row.atom_index << ',' << row.residue_index << ','
        << row.residue_number << ',' << row.amino_acid << ',' << row.element << ','
        << row.atom_name << ',' << row.frame_variant << ',' << row.h5_row << ','
        << row.original_index << ',' << num(row.time_ps) << ','
        << num(row.frame_z.x()) << ',' << num(row.frame_z.y()) << ',' << num(row.frame_z.z())
        << ',' << num(row.frame_x.x()) << ',' << num(row.frame_x.y()) << ','
        << num(row.frame_x.z()) << ',' << num(row.frame_y.x()) << ','
        << num(row.frame_y.y()) << ',' << num(row.frame_y.z()) << ','
        << row.frame_variant << ',' << (row.frame_valid ? 1 : 0) << ','
        << row.frame_anchor_atom_index << ','
        << (row.dft_present ? 1 : 0) << ',' << num(row.dft_present ? row.dft_total_decomp.T0 : kNaN);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            out << ',' << num(row.dft_present ? row.dft_total_raw(i, j) : kNaN);
    for (double v : row.dft_total_decomp.T2) out << ',' << num(row.dft_present ? v : kNaN);
    out << ',' << (row.dft_local_frame_valid ? 1 : 0) << ','
        << (row.aimnet2_charge_present ? 1 : 0) << ','
        << num(row.aimnet2_charge_present ? row.aimnet2_charge : kNaN)
        << ",aimnet2,"
        << (crgPresent ? 1 : 0) << ','
        << num(crgPresent ? row.crg_scalar : kNaN) << ','
        << num(localCrgPresent ? row.crg_vector_local.x() : kNaN) << ','
        << num(localCrgPresent ? row.crg_vector_local.y() : kNaN) << ','
        << num(localCrgPresent ? row.crg_vector_local.z() : kNaN) << ','
        << num(crgPresent ? row.crg_vector_lab.x() : kNaN) << ','
        << num(crgPresent ? row.crg_vector_lab.y() : kNaN) << ','
        << num(crgPresent ? row.crg_vector_lab.z() : kNaN)
        << ",aimnet2_charge_response_gradient_not_polarizability,"
        << (embeddingPresent ? 1 : 0) << ','
        << static_cast<qulonglong>(embeddingPresent ? row.embedding_dims : 0)
        << ",learnable_ceiling_feature_not_physical_law\n";

    appendVec3(crgVectorLocal_, row.crg_vector_local, localCrgPresent);
    appendVec3(crgVectorLab_, row.crg_vector_lab, crgPresent);
    crgScalar_.push_back(crgPresent ? row.crg_scalar : kNaN);
    if (embeddingPresent) {
        for (std::size_t i = 0; i < embeddingDims_; ++i) embedding_.push_back(row.embedding[i]);
    } else {
        for (std::size_t i = 0; i < embeddingDims_; ++i) embedding_.push_back(kNaNF32);
    }
    appendT2(targetLocalT2_, row.dft_target_local_T2, targetLocalPresent);
    appendT2(targetLabT2_, row.dft_total_decomp.T2, row.dft_present);
    ++rows_;
}

bool Aimnet2FeatureSink::Commit() {
    if (!ok_) return false;
    if (aggregatedOut_) aggregatedOut_->flush();
    const bool csvOk = aggregatedFile_ && aggregatedFile_->commit();
    const QString outDir = QFileInfo(aggregatedPath_).dir().absolutePath();
    bool npyOk = true;
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[0]),
                                {rows_, 3}, crgVectorLocal_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[1]),
                                {rows_, 3}, crgVectorLab_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[2]),
                                {rows_}, crgScalar_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<float>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[3]),
                               {rows_, embeddingDims_}, embedding_, QByteArray("<f4"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[4]),
                                {rows_, 5}, targetLocalT2_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[5]),
                                {rows_, 5}, targetLabT2_, QByteArray("<f8"));
    if (!csvOk) qCWarning(cAimnet2Sink).noquote() << "AIMNet2 aggregated commit failed" << aggregatedPath_;
    if (!npyOk) qCWarning(cAimnet2Sink).noquote() << "AIMNet2 sidecar NPY commit failed" << caseName_;
    qCInfo(cAimnet2Sink).noquote() << "committed" << caseName_ << "| rows=" << rows_
                                   << "| embedding_dims=" << embeddingDims_;
    return csvOk && npyOk;
}

}  // namespace h5reader::rediscover

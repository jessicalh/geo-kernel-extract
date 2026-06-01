#include "BuckinghamEfieldSink.h"

#include <QByteArray>
#include <QDir>
#include <QFileInfo>
#include <QIODevice>
#include <QLoggingCategory>

#include <cmath>
#include <limits>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cBuckinghamSink, "h5reader.rediscover.buckingham_sink")

QString num(double v) { return QString::number(v, 'g', 12); }

const double kNaN = std::numeric_limits<double>::quiet_NaN();

void appendVec3(std::vector<double>& dst, const Vec3& v, bool present) {
    dst.push_back(present ? v.x() : kNaN);
    dst.push_back(present ? v.y() : kNaN);
    dst.push_back(present ? v.z() : kNaN);
}

void appendArray(std::vector<double>& dst, const std::array<double, 3>& a, bool present) {
    for (double v : a) dst.push_back(present ? v : kNaN);
}

void appendArray(std::vector<double>& dst, const std::array<double, 5>& a, bool present) {
    for (double v : a) dst.push_back(present ? v : kNaN);
}

bool writeNpyF64(const QString& path, const std::vector<std::size_t>& shape,
                 const std::vector<double>& data) {
    std::size_t n = 1;
    for (const std::size_t dim : shape) n *= dim;
    if (n != data.size()) return false;

    QByteArray header;
    header += "{'descr': '<f8', 'fortran_order': False, 'shape': (";
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
    const qsizetype payloadBytes = static_cast<qsizetype>(data.size() * sizeof(double));
    if (payloadBytes > 0
        && f.write(reinterpret_cast<const char*>(data.data()), payloadBytes) != payloadBytes)
        return false;
    return f.commit();
}

const char* kHeader =
    "row_id,atom_index,residue_index,residue_number,amino_acid_ord,element_ord,"
    "atom_name,backbone_frame_class,h5_row,original_index,time_ps,"
    "frame_z_x,frame_z_y,frame_z_z,frame_x_x,frame_x_y,frame_x_z,"
    "frame_y_x,frame_y_y,frame_y_z,frame_variant,frame_valid,frame_anchor_atom_index,"
    "dft_present,apbs_efield_present,efield_units,"
    "efield_local_x,efield_local_y,efield_local_z,E_proj,E_mag,"
    "dft_sigma_iso,dft_total_raw_00,dft_total_raw_01,dft_total_raw_02,"
    "dft_total_raw_10,dft_total_raw_11,dft_total_raw_12,"
    "dft_total_raw_20,dft_total_raw_21,dft_total_raw_22,"
    "dft_total_T1_0,dft_total_T1_1,dft_total_T1_2,"
    "dft_total_T2_0,dft_total_T2_1,dft_total_T2_2,dft_total_T2_3,dft_total_T2_4,"
    "target_T1_status";

}  // namespace

bool FiniteVec3(const Vec3& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

BuckinghamEfieldSink::BuckinghamEfieldSink(const QString& outDir, const QString& caseName)
    : caseName_(caseName) {
    QDir().mkpath(outDir);
    aggregatedPath_ = QStringLiteral("%1/%2_aggregated.csv").arg(outDir, caseName);
    sidecarFiles_ = {
        QStringLiteral("%1_feature_field_local.npy").arg(caseName),
        QStringLiteral("%1_target_T1_unverified.npy").arg(caseName),
        QStringLiteral("%1_target_T2.npy").arg(caseName),
    };

    aggregatedFile_ = std::make_unique<QSaveFile>(aggregatedPath_);
    if (!aggregatedFile_->open(QIODevice::WriteOnly | QIODevice::Text)) {
        qCWarning(cBuckinghamSink).noquote() << "cannot open buckingham aggregated CSV"
                                             << aggregatedPath_;
        return;
    }
    aggregatedOut_ = std::make_unique<QTextStream>(aggregatedFile_.get());
    *aggregatedOut_ << kHeader << '\n';
    ok_ = true;
}

BuckinghamEfieldSink::~BuckinghamEfieldSink() = default;

void BuckinghamEfieldSink::Write(const BuckinghamEfieldRow& row) {
    if (!ok_) return;
    const int64_t rowId = nextRowId_++;
    const bool featurePresent =
        row.apbs_efield_present && row.frame_valid && FiniteVec3(row.efield_local)
        && std::isfinite(row.e_proj) && std::isfinite(row.e_mag);

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
        << (row.dft_present ? 1 : 0) << ',' << (row.apbs_efield_present ? 1 : 0)
        << ',' << row.efield_units << ','
        << num(featurePresent ? row.efield_local.x() : kNaN) << ','
        << num(featurePresent ? row.efield_local.y() : kNaN) << ','
        << num(featurePresent ? row.efield_local.z() : kNaN) << ','
        << num(featurePresent ? row.e_proj : kNaN) << ','
        << num(row.apbs_efield_present ? row.e_mag : kNaN) << ','
        << num(row.dft_present ? row.dft_total_decomp.T0 : kNaN);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            out << ',' << num(row.dft_present ? row.dft_total_raw(i, j) : kNaN);
    for (double v : row.dft_total_decomp.T1) out << ',' << num(row.dft_present ? v : kNaN);
    for (double v : row.dft_total_decomp.T2) out << ',' << num(row.dft_present ? v : kNaN);
    out << ",unverified_emit_only\n";

    appendVec3(featureFieldLocal_, row.efield_local, featurePresent);
    appendArray(targetT1_, row.dft_total_decomp.T1, row.dft_present);
    appendArray(targetT2_, row.dft_total_decomp.T2, row.dft_present);
    ++rows_;
}

bool BuckinghamEfieldSink::Commit() {
    if (!ok_) return false;
    if (aggregatedOut_) aggregatedOut_->flush();
    const bool csvOk = aggregatedFile_ && aggregatedFile_->commit();
    const QString outDir = QFileInfo(aggregatedPath_).dir().absolutePath();
    bool npyOk = true;
    npyOk = npyOk
            && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[0]),
                           {rows_, 3}, featureFieldLocal_);
    npyOk = npyOk
            && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[1]),
                           {rows_, 3}, targetT1_);
    npyOk = npyOk
            && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[2]),
                           {rows_, 5}, targetT2_);
    if (!csvOk)
        qCWarning(cBuckinghamSink).noquote() << "buckingham aggregated commit failed"
                                             << aggregatedPath_;
    if (!npyOk)
        qCWarning(cBuckinghamSink).noquote() << "buckingham sidecar NPY commit failed"
                                             << caseName_;
    qCInfo(cBuckinghamSink).noquote() << "committed" << caseName_ << "| rows=" << rows_;
    return csvOk && npyOk;
}

}  // namespace h5reader::rediscover

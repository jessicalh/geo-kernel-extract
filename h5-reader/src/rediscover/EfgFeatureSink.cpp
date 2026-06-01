#include "EfgFeatureSink.h"

#include "SphericalBasis.h"

#include <QByteArray>
#include <QDir>
#include <QFileInfo>
#include <QIODevice>
#include <QLoggingCategory>

#include <cmath>
#include <limits>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cEfgSink, "h5reader.rediscover.efg_sink")

QString num(double v) { return QString::number(v, 'g', 12); }

const double kNaN = std::numeric_limits<double>::quiet_NaN();

void appendT2(std::vector<double>& dst, const std::array<double, 5>& t2, bool present) {
    for (double v : t2) dst.push_back(present ? v : kNaN);
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

}  // namespace

std::array<double, 5> DecomposeEfgLibraryT2(const Mat3& efg) {
    return DecomposeLibrary(efg).T2;
}

Mat3 ReconstructLibraryT2(const std::array<double, 5>& t2) {
    const double kSqrt2 = std::sqrt(2.0);
    const double kSqrtThreeHalves = std::sqrt(3.0 / 2.0);

    const double Sxy = t2[0] / kSqrt2;
    const double Syz = t2[1] / kSqrt2;
    const double Szz = t2[2] / kSqrtThreeHalves;
    const double Sxz = t2[3] / kSqrt2;
    const double SxxMinusSyy = kSqrt2 * t2[4];
    const double Sxx = 0.5 * (-Szz + SxxMinusSyy);
    const double Syy = 0.5 * (-Szz - SxxMinusSyy);

    Mat3 m = Mat3::Zero();
    m(0, 0) = Sxx;
    m(1, 1) = Syy;
    m(2, 2) = Szz;
    m(0, 1) = m(1, 0) = Sxy;
    m(1, 2) = m(2, 1) = Syz;
    m(0, 2) = m(2, 0) = Sxz;
    return m;
}

double T2Magnitude(const std::array<double, 5>& t2) {
    double s = 0.0;
    for (double v : t2) s += v * v;
    return std::sqrt(s);
}

bool FiniteT2(const std::array<double, 5>& t2) {
    for (double v : t2)
        if (!std::isfinite(v)) return false;
    return true;
}

EfgFeatureSink::EfgFeatureSink(const QString& outDir, const QString& caseName)
    : caseName_(caseName) {
    QDir().mkpath(outDir);
    aggregatedPath_ = QStringLiteral("%1/%2_aggregated.csv").arg(outDir, caseName);
    sidecarFiles_ = {
        QStringLiteral("%1_feature_T2.npy").arg(caseName),
        QStringLiteral("%1_target_T2.npy").arg(caseName),
        QStringLiteral("%1_feature_lab_T2.npy").arg(caseName),
        QStringLiteral("%1_target_lab_T2.npy").arg(caseName),
    };

    aggregatedFile_ = std::make_unique<QSaveFile>(aggregatedPath_);
    if (!aggregatedFile_->open(QIODevice::WriteOnly | QIODevice::Text)) {
        qCWarning(cEfgSink).noquote() << "cannot open efg aggregated CSV" << aggregatedPath_;
        return;
    }
    aggregatedOut_ = std::make_unique<QTextStream>(aggregatedFile_.get());
    *aggregatedOut_
        << "row_id,atom_index,residue_index,residue_number,amino_acid_ord,"
           "element_ord,atom_name,backbone_frame_class,h5_row,original_index,time_ps,"
           "frame_z_x,frame_z_y,frame_z_z,frame_x_x,frame_x_y,frame_x_z,"
           "frame_y_x,frame_y_y,frame_y_z,frame_variant,frame_valid,"
           "frame_anchor_atom_index,dft_present,apbs_efg_present,"
           "efg_T2_magnitude,efg_T2_lab_magnitude,efg_units\n";
    ok_ = true;
}

EfgFeatureSink::~EfgFeatureSink() = default;

void EfgFeatureSink::Write(const EfgFeatureRow& row) {
    if (!ok_) return;
    const int64_t rowId = nextRowId_++;
    const bool localFeaturePresent =
        row.apbs_efg_present && row.frame_valid && FiniteT2(row.efg_feature_T2);
    const bool localTargetPresent =
        row.dft_present && row.frame_valid && FiniteT2(row.dft_target_T2);
    const double mag = localFeaturePresent ? T2Magnitude(row.efg_feature_T2) : kNaN;
    const double labMag = row.apbs_efg_present ? T2Magnitude(row.efg_feature_lab_T2) : kNaN;
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
        << (row.dft_present ? 1 : 0) << ',' << (row.apbs_efg_present ? 1 : 0) << ','
        << num(mag) << ',' << num(labMag) << ',' << row.efg_units << '\n';

    appendT2(featureT2_, row.efg_feature_T2, localFeaturePresent);
    appendT2(targetT2_, row.dft_target_T2, localTargetPresent);
    appendT2(featureLabT2_, row.efg_feature_lab_T2, row.apbs_efg_present);
    appendT2(targetLabT2_, row.dft_target_lab_T2, row.dft_present);
    ++rows_;
}

bool EfgFeatureSink::Commit() {
    if (!ok_) return false;
    if (aggregatedOut_) aggregatedOut_->flush();
    const bool csvOk = aggregatedFile_ && aggregatedFile_->commit();
    const QString outDir = QFileInfo(aggregatedPath_).dir().absolutePath();
    bool npyOk = true;
    npyOk = npyOk
            && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[0]),
                           {rows_, 5}, featureT2_);
    npyOk = npyOk
            && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[1]),
                           {rows_, 5}, targetT2_);
    npyOk = npyOk
            && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[2]),
                           {rows_, 5}, featureLabT2_);
    npyOk = npyOk
            && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[3]),
                           {rows_, 5}, targetLabT2_);
    if (!csvOk) qCWarning(cEfgSink).noquote() << "efg aggregated commit failed" << aggregatedPath_;
    if (!npyOk) qCWarning(cEfgSink).noquote() << "efg sidecar NPY commit failed" << caseName_;
    qCInfo(cEfgSink).noquote() << "committed" << caseName_ << "| rows=" << rows_;
    return csvOk && npyOk;
}

}  // namespace h5reader::rediscover

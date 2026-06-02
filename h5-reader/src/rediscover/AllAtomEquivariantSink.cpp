#include "AllAtomEquivariantSink.h"

#include <QByteArray>
#include <QDir>
#include <QFileInfo>
#include <QIODevice>
#include <QLoggingCategory>

#include <cmath>
#include <limits>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cAllAtomEquivariantSink, "h5reader.rediscover.all_atom_equivariant_sink")

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const float kNaNF32 = std::numeric_limits<float>::quiet_NaN();

QString num(double v) { return QString::number(v, 'g', 12); }

bool finiteVec3(const Vec3& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

bool finiteT2(const std::array<double, 5>& t2) {
    for (double v : t2)
        if (!std::isfinite(v)) return false;
    return true;
}

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

double t2Magnitude(const std::array<double, 5>& t2) {
    double s = 0.0;
    for (double v : t2) s += v * v;
    return std::sqrt(s);
}

const char* kTargetsHeader =
    "row_id,atom_index,residue_index,residue_number,amino_acid_ord,element_ord,"
    "atom_name,h5_row,original_index,time_ps,frame_name,frame_valid,"
    "dft_present,dft_sigma_iso,"
    "dft_total_raw_00,dft_total_raw_01,dft_total_raw_02,"
    "dft_total_raw_10,dft_total_raw_11,dft_total_raw_12,"
    "dft_total_raw_20,dft_total_raw_21,dft_total_raw_22,"
    "dft_total_T2_0,dft_total_T2_1,dft_total_T2_2,dft_total_T2_3,dft_total_T2_4,"
    "apbs_efield_present,apbs_efield_x,apbs_efield_y,apbs_efield_z,apbs_efield_mag,"
    "apbs_efg_present,apbs_efg_T2_0,apbs_efg_T2_1,apbs_efg_T2_2,"
    "apbs_efg_T2_3,apbs_efg_T2_4,apbs_efg_T2_mag,"
    "aimnet2_charge_present,aimnet2_charge,"
    "aimnet2_charge_response_gradient_present,aimnet2_charge_response_gradient_scalar,"
    "aimnet2_charge_response_gradient_x,aimnet2_charge_response_gradient_y,"
    "aimnet2_charge_response_gradient_z,aimnet2_charge_response_gradient_mag,"
    "aimnet2_embedding_present,aimnet2_embedding_dim";

const char* kSourcesHeader =
    "row_id,target_atom_index,target_residue_index,h5_row,original_index,time_ps,"
    "mechanism,source_kind,category,category_ord,source_id,"
    "source_atom_index,source_residue_index,source_residue_number,source_element_ord,"
    "source_atom_name,"
    "disp_x,disp_y,disp_z,r,inv_r3,cos_theta,dipolar,"
    "orientation_a_x,orientation_a_y,orientation_a_z,"
    "orientation_b_x,orientation_b_y,orientation_b_z,"
    "source_value,source_value_2,source_units,source_present,"
    "value_vec_x,value_vec_y,value_vec_z,value_vec_mag,"
    "tensor_T2_0,tensor_T2_1,tensor_T2_2,tensor_T2_3,tensor_T2_4,tensor_present,"
    "ring_index,ring_type_index,ring_aromaticity,ring_size,ring_nitrogen,"
    "ring_fused,ring_intensity,ring_jb_offset,"
    "bond_index,bond_category,bond_order,bond_atom_a,bond_atom_b,bond_elem_a,bond_elem_b,"
    "bond_length,source_is_self_or_bonded,"
    "charge_source,source_q_e,q_over_r3,"
    "aimnet2_embedding_present,aimnet2_embedding_dim";

}  // namespace

AllAtomEquivariantSink::AllAtomEquivariantSink(const QString& outDir, const QString& caseName,
                                               std::size_t embeddingDims)
    : caseName_(caseName), embeddingDims_(embeddingDims) {
    QDir().mkpath(outDir);
    sourcesPath_ = QStringLiteral("%1/%2_sources.csv").arg(outDir, caseName);
    targetsPath_ = QStringLiteral("%1/%2_targets.csv").arg(outDir, caseName);
    sidecarFiles_ = {
        QStringLiteral("%1_target_T2.npy").arg(caseName),
        QStringLiteral("%1_target_sigma_iso.npy").arg(caseName),
        QStringLiteral("%1_target_raw.npy").arg(caseName),
        QStringLiteral("%1_apbs_efield.npy").arg(caseName),
        QStringLiteral("%1_apbs_efg_T2.npy").arg(caseName),
        QStringLiteral("%1_aimnet2_charge.npy").arg(caseName),
        QStringLiteral("%1_aimnet2_charge_response_gradient.npy").arg(caseName),
        QStringLiteral("%1_aimnet2_charge_response_gradient_scalar.npy").arg(caseName),
        QStringLiteral("%1_aimnet2_embedding.npy").arg(caseName),
    };

    sourcesFile_ = std::make_unique<QSaveFile>(sourcesPath_);
    targetsFile_ = std::make_unique<QSaveFile>(targetsPath_);
    if (!sourcesFile_->open(QIODevice::WriteOnly | QIODevice::Text)) {
        qCWarning(cAllAtomEquivariantSink).noquote() << "cannot open all-atom sources CSV"
                                                     << sourcesPath_;
        return;
    }
    if (!targetsFile_->open(QIODevice::WriteOnly | QIODevice::Text)) {
        qCWarning(cAllAtomEquivariantSink).noquote() << "cannot open all-atom targets CSV"
                                                     << targetsPath_;
        return;
    }
    sourcesOut_ = std::make_unique<QTextStream>(sourcesFile_.get());
    targetsOut_ = std::make_unique<QTextStream>(targetsFile_.get());
    *sourcesOut_ << kSourcesHeader << '\n';
    *targetsOut_ << kTargetsHeader << '\n';
    ok_ = true;
}

AllAtomEquivariantSink::~AllAtomEquivariantSink() = default;

int64_t AllAtomEquivariantSink::WriteTarget(const AllAtomEquivariantTargetRecord& row) {
    if (!ok_) return -1;
    const int64_t rowId = nextRowId_++;
    const bool dftPresent = row.target.present && finiteT2(row.target.total_decomp.T2);
    const bool efieldPresent = row.apbs_efield_present && finiteVec3(row.apbs_efield_lab);
    const bool efgPresent = row.apbs_efg_present && finiteT2(row.apbs_efg_T2);
    const bool crgPresent =
        row.aimnet2_crg_present && std::isfinite(row.aimnet2_crg_scalar)
        && finiteVec3(row.aimnet2_crg_lab);
    const bool embeddingPresent =
        row.aimnet2_embedding_present && row.aimnet2_embedding
        && row.aimnet2_embedding_dims == embeddingDims_;

    QTextStream& out = *targetsOut_;
    out << rowId << ',' << row.atom_index << ',' << row.residue_index << ','
        << row.residue_number << ',' << row.amino_acid << ',' << row.element << ','
        << row.atom_name << ',' << row.h5_row << ',' << row.original_index << ','
        << num(row.time_ps)
        << ",molecular_lab_h5_orca_aligned,1,"
        << (row.target.present ? 1 : 0) << ','
        << num(row.target.present ? row.target.total_decomp.T0 : kNaN);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            out << ',' << num(row.target.present ? row.target.total_raw(i, j) : kNaN);
    for (double v : row.target.total_decomp.T2) out << ',' << num(dftPresent ? v : kNaN);
    out << ',' << (efieldPresent ? 1 : 0) << ','
        << num(efieldPresent ? row.apbs_efield_lab.x() : kNaN) << ','
        << num(efieldPresent ? row.apbs_efield_lab.y() : kNaN) << ','
        << num(efieldPresent ? row.apbs_efield_lab.z() : kNaN) << ','
        << num(efieldPresent ? row.apbs_efield_lab.norm() : kNaN)
        << ',' << (efgPresent ? 1 : 0);
    for (double v : row.apbs_efg_T2) out << ',' << num(efgPresent ? v : kNaN);
    out << ',' << num(efgPresent ? t2Magnitude(row.apbs_efg_T2) : kNaN)
        << ',' << (row.aimnet2_charge_present ? 1 : 0) << ','
        << num(row.aimnet2_charge_present ? row.aimnet2_charge : kNaN)
        << ',' << (crgPresent ? 1 : 0) << ','
        << num(crgPresent ? row.aimnet2_crg_scalar : kNaN) << ','
        << num(crgPresent ? row.aimnet2_crg_lab.x() : kNaN) << ','
        << num(crgPresent ? row.aimnet2_crg_lab.y() : kNaN) << ','
        << num(crgPresent ? row.aimnet2_crg_lab.z() : kNaN) << ','
        << num(crgPresent ? row.aimnet2_crg_lab.norm() : kNaN)
        << ',' << (embeddingPresent ? 1 : 0) << ','
        << static_cast<qulonglong>(embeddingPresent ? row.aimnet2_embedding_dims : 0)
        << '\n';

    appendT2(targetT2_, row.target.total_decomp.T2, dftPresent);
    targetSigmaIso_.push_back(row.target.present ? row.target.total_decomp.T0 : kNaN);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            targetRaw_.push_back(row.target.present ? row.target.total_raw(i, j) : kNaN);
    appendVec3(apbsEfield_, row.apbs_efield_lab, efieldPresent);
    appendT2(apbsEfgT2_, row.apbs_efg_T2, efgPresent);
    aimnet2Charge_.push_back(row.aimnet2_charge_present ? row.aimnet2_charge : kNaN);
    appendVec3(aimnet2Crg_, row.aimnet2_crg_lab, crgPresent);
    aimnet2CrgScalar_.push_back(crgPresent ? row.aimnet2_crg_scalar : kNaN);
    if (embeddingPresent) {
        for (std::size_t i = 0; i < embeddingDims_; ++i)
            aimnet2Embedding_.push_back(row.aimnet2_embedding[i]);
    } else {
        for (std::size_t i = 0; i < embeddingDims_; ++i) aimnet2Embedding_.push_back(kNaNF32);
    }
    ++targetRows_;
    return rowId;
}

void AllAtomEquivariantSink::WriteSource(int64_t rowId,
                                         const AllAtomEquivariantSourceRecord& row) {
    if (!ok_ || rowId < 0) return;
    const bool tensorPresent = row.tensor_present && finiteT2(row.tensor_T2);
    QTextStream& out = *sourcesOut_;
    out << rowId << ',' << row.target_atom_index << ',' << row.target_residue_index << ','
        << row.h5_row << ',' << row.original_index << ',' << num(row.time_ps) << ','
        << row.mechanism << ',' << row.source_kind << ',' << row.category << ','
        << row.category_ord << ',' << row.source_id << ','
        << row.source_atom_index << ',' << row.source_residue_index << ','
        << row.source_residue_number << ',' << row.source_element << ','
        << row.source_atom_name << ','
        << num(row.disp.x()) << ',' << num(row.disp.y()) << ',' << num(row.disp.z()) << ','
        << num(row.r) << ',' << num(row.inv_r3) << ',' << num(row.cos_theta) << ','
        << num(row.dipolar) << ','
        << num(row.orientation_a.x()) << ',' << num(row.orientation_a.y()) << ','
        << num(row.orientation_a.z()) << ','
        << num(row.orientation_b.x()) << ',' << num(row.orientation_b.y()) << ','
        << num(row.orientation_b.z()) << ','
        << num(row.source_value) << ',' << num(row.source_value_2) << ','
        << row.source_units << ',' << (row.source_present ? 1 : 0) << ','
        << num(row.value_vec.x()) << ',' << num(row.value_vec.y()) << ','
        << num(row.value_vec.z()) << ',' << num(row.value_vec_mag);
    for (double v : row.tensor_T2) out << ',' << num(tensorPresent ? v : kNaN);
    out << ',' << (tensorPresent ? 1 : 0) << ','
        << row.ring_index << ',' << row.ring_type_index << ',' << row.ring_aromaticity << ','
        << row.ring_size << ',' << row.ring_nitrogen << ',' << (row.ring_fused ? 1 : 0)
        << ',' << num(row.ring_intensity) << ',' << num(row.ring_jb_offset) << ','
        << row.bond_index << ',' << row.bond_category << ',' << row.bond_order << ','
        << row.bond_atom_a << ',' << row.bond_atom_b << ',' << row.bond_elem_a << ','
        << row.bond_elem_b << ',' << num(row.bond_length) << ','
        << (row.source_is_self_or_bonded ? 1 : 0) << ','
        << row.charge_source << ',' << num(row.source_q_e) << ',' << num(row.q_over_r3)
        << ',' << (row.aimnet2_embedding_present ? 1 : 0) << ','
        << static_cast<qulonglong>(row.aimnet2_embedding_present ? row.aimnet2_embedding_dims : 0)
        << '\n';
    ++sourceRows_;
}

bool AllAtomEquivariantSink::Commit() {
    if (!ok_) return false;
    if (sourcesOut_) sourcesOut_->flush();
    if (targetsOut_) targetsOut_->flush();
    const bool sourcesOk = sourcesFile_ && sourcesFile_->commit();
    const bool targetsOk = targetsFile_ && targetsFile_->commit();
    const QString outDir = QFileInfo(targetsPath_).dir().absolutePath();
    bool npyOk = true;
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[0]),
                                {targetRows_, 5}, targetT2_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[1]),
                                {targetRows_}, targetSigmaIso_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[2]),
                                {targetRows_, 3, 3}, targetRaw_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[3]),
                                {targetRows_, 3}, apbsEfield_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[4]),
                                {targetRows_, 5}, apbsEfgT2_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[5]),
                                {targetRows_}, aimnet2Charge_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[6]),
                                {targetRows_, 3}, aimnet2Crg_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<double>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[7]),
                                {targetRows_}, aimnet2CrgScalar_, QByteArray("<f8"));
    npyOk = npyOk
            && writeNpy<float>(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[8]),
                               {targetRows_, embeddingDims_}, aimnet2Embedding_,
                               QByteArray("<f4"));

    if (!sourcesOk)
        qCWarning(cAllAtomEquivariantSink).noquote() << "all-atom sources commit failed"
                                                     << sourcesPath_;
    if (!targetsOk)
        qCWarning(cAllAtomEquivariantSink).noquote() << "all-atom targets commit failed"
                                                     << targetsPath_;
    if (!npyOk)
        qCWarning(cAllAtomEquivariantSink).noquote() << "all-atom sidecar NPY commit failed"
                                                     << caseName_;
    qCInfo(cAllAtomEquivariantSink).noquote() << "committed" << caseName_
                                              << "| target_rows=" << targetRows_
                                              << "| source_rows=" << sourceRows_;
    return sourcesOk && targetsOk && npyOk;
}

}  // namespace h5reader::rediscover

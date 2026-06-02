#include "BroadBackboneSink.h"

#include "SphericalBasis.h"

#include "../diagnostics/ErrorBus.h"

#include <QByteArray>
#include <QDir>
#include <QFileInfo>
#include <QIODevice>
#include <QLoggingCategory>

#include <array>
#include <limits>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cBroadSink, "h5reader.rediscover.broad_sink")

using Severity = h5reader::diagnostics::Severity;

void report(Severity sev, const QString& msg, const QString& ctx) {
    h5reader::diagnostics::ErrorBus::Report(sev, QStringLiteral("BroadBackboneSink"), msg, ctx);
}

QString num(double v) { return QString::number(v, 'g', 9); }

double kernelComponent(const BroadKernelT2& k, int i) {
    if (!k.present) return std::numeric_limits<double>::quiet_NaN();
    return k.T2[static_cast<std::size_t>(i)];
}

void writeKernelT2(QTextStream& out, const BroadKernelT2& k) {
    out << ',' << (k.present ? 1 : 0);
    for (int i = 0; i < 5; ++i) out << ',' << num(kernelComponent(k, i));
}

void appendKernelT2(std::vector<double>& dst, const BroadKernelT2& k) {
    for (int i = 0; i < 5; ++i) dst.push_back(kernelComponent(k, i));
}

double presentValue(bool present, double value) {
    return present ? value : std::numeric_limits<double>::quiet_NaN();
}

void writeMcLitKernel(QTextStream& out, const SphericalTensor& k, bool present) {
    out << ',' << (present ? 1 : 0) << ',' << num(presentValue(present, k.T0));
    for (double v : k.T2) out << ',' << num(presentValue(present, v));
}

SphericalTensor sumMcLitKernel(const NeighborhoodRecord& rec, bool* present) {
    SphericalTensor out;
    bool any = false;
    for (const SourceSlot& s : rec.sources) {
        if (s.kind != SourceKind::Bond || !s.bond_mc_lit_kernel_present) continue;
        any = true;
        out.T0 += s.bond_mc_lit_kernel.T0;
        for (int i = 0; i < 3; ++i)
            out.T1[static_cast<std::size_t>(i)] +=
                s.bond_mc_lit_kernel.T1[static_cast<std::size_t>(i)];
        for (int i = 0; i < 5; ++i)
            out.T2[static_cast<std::size_t>(i)] +=
                s.bond_mc_lit_kernel.T2[static_cast<std::size_t>(i)];
    }
    if (present) *present = any;
    return out;
}

SphericalTensor sumValidMcLitKernel(const NeighborhoodRecord& rec, bool* present) {
    SphericalTensor out;
    bool any = false;
    for (const SourceSlot& s : rec.sources) {
        if (s.kind != SourceKind::Bond || !s.bond_mc_lit_kernel_present) continue;
        if (s.mc_source_is_self_or_bonded) continue;
        any = true;
        out.T0 += s.bond_mc_lit_kernel.T0;
        for (int i = 0; i < 3; ++i)
            out.T1[static_cast<std::size_t>(i)] +=
                s.bond_mc_lit_kernel.T1[static_cast<std::size_t>(i)];
        for (int i = 0; i < 5; ++i)
            out.T2[static_cast<std::size_t>(i)] +=
                s.bond_mc_lit_kernel.T2[static_cast<std::size_t>(i)];
    }
    if (present) *present = any;
    return out;
}

// Same little-endian f64 NPY writer as RecordSink (kept local — additive, no
// coupling). Returns false on any shape/IO mismatch.
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

// ── Source-row schema header (per (atom,frame,source); NO target columns) ───
const char* kSourceHeader =
    "row_id,atom_index,residue_number,element_ord,atom_name,frame_variant,"
    "mechanism,source_kind,"
    "disp_local_x,disp_local_y,disp_local_z,r,cos_theta,dipolar,"
    "ring_index,ring_type_index,ring_intensity,ring_nitrogen,is_self_or_bonded,"
    "bond_index,bond_category,bond_atom_a,bond_atom_b,"
    "source_atom_index,source_residue_number,source_element_ord,source_q_e,"
    "source_normal_local_x,source_normal_local_y,source_normal_local_z,"
    "bond_axis_local_x,bond_axis_local_y,bond_axis_local_z,"
    "mc_lit_kernel_present,mc_lit_T0,"
    "mc_lit_T2_local_0,mc_lit_T2_local_1,mc_lit_T2_local_2,"
    "mc_lit_T2_local_3,mc_lit_T2_local_4,"
    "mc_source_is_self_or_bonded";

// ── Aggregated-row schema header (per (atom,frame); the target lives ONCE) ──
const char* kAggregatedHeader =
    "row_id,atom_index,residue_index,residue_number,amino_acid_ord,element_ord,"
    "atom_name,backbone_frame_class,h5_row,original_index,time_ps,"
    "frame_z_x,frame_z_y,frame_z_z,frame_x_x,frame_x_y,frame_x_z,"
    "frame_y_x,frame_y_y,frame_y_z,frame_variant,frame_valid,frame_anchor_atom_index,"
    "ring_n,ring_sum_dipolar,ring_cutoff_A,"
    "bond_n,bond_sum_dipolar,bond_cutoff_A,"
    "mc_lit_kernel_present,mc_lit_T0,"
    "mc_lit_T2_local_0,mc_lit_T2_local_1,mc_lit_T2_local_2,"
    "mc_lit_T2_local_3,mc_lit_T2_local_4,"
    "charge_n,charge_source,charge_cutoff_A,"
    "field_local_x,field_local_y,field_local_z,field_z,field_mag,"
    "mu_local_x,mu_local_y,mu_local_z,"
    "literature_kernel_present,"
    "literature_kernel_T2_0,literature_kernel_T2_1,literature_kernel_T2_2,"
    "literature_kernel_T2_3,literature_kernel_T2_4,"
    "ring_literature_kernel_present,"
    "ring_literature_kernel_T2_0,ring_literature_kernel_T2_1,"
    "ring_literature_kernel_T2_2,ring_literature_kernel_T2_3,"
    "ring_literature_kernel_T2_4,"
    "bond_literature_kernel_present,"
    "bond_literature_kernel_T2_0,bond_literature_kernel_T2_1,"
    "bond_literature_kernel_T2_2,bond_literature_kernel_T2_3,"
    "bond_literature_kernel_T2_4,"
    "charge_literature_kernel_present,"
    "charge_literature_kernel_T2_0,charge_literature_kernel_T2_1,"
    "charge_literature_kernel_T2_2,charge_literature_kernel_T2_3,"
    "charge_literature_kernel_T2_4,"
    "dft_present,dft_sigma_iso,"
    "dft_total_raw_00,dft_total_raw_01,dft_total_raw_02,"
    "dft_total_raw_10,dft_total_raw_11,dft_total_raw_12,"
    "dft_total_raw_20,dft_total_raw_21,dft_total_raw_22,"
    "dft_total_T1_0,dft_total_T1_1,dft_total_T1_2,"
    "dft_total_T2_0,dft_total_T2_1,dft_total_T2_2,dft_total_T2_3,dft_total_T2_4,"
    "dft_local_frame_valid,"
    "bond_n_valid,bond_sum_dipolar_valid,mc_source_near_field_ratio,"
    "mc_lit_valid_kernel_present,mc_lit_T0_valid,"
    "mc_lit_T2_local_valid_0,mc_lit_T2_local_valid_1,mc_lit_T2_local_valid_2,"
    "mc_lit_T2_local_valid_3,mc_lit_T2_local_valid_4";

}  // namespace

BroadBackboneSink::BroadBackboneSink(const QString& outDir, const QString& caseName)
    : caseName_(caseName) {
    QDir().mkpath(outDir);
    sourcesPath_ = QStringLiteral("%1/%2_sources.csv").arg(outDir, caseName);
    aggregatedPath_ = QStringLiteral("%1/%2_aggregated.csv").arg(outDir, caseName);
    sidecarFiles_ = {
        QStringLiteral("%1_aggregated_target_T2.npy").arg(caseName),
        QStringLiteral("%1_aggregated_target_local_T2.npy").arg(caseName),
        QStringLiteral("%1_aggregated_field_local.npy").arg(caseName),
        QStringLiteral("%1_aggregated_literature_kernel_T2.npy").arg(caseName),
        QStringLiteral("%1_aggregated_ring_literature_kernel_T2.npy").arg(caseName),
        QStringLiteral("%1_aggregated_bond_literature_kernel_T2.npy").arg(caseName),
        QStringLiteral("%1_aggregated_charge_literature_kernel_T2.npy").arg(caseName),
    };

    sourcesFile_ = std::make_unique<QSaveFile>(sourcesPath_);
    aggregatedFile_ = std::make_unique<QSaveFile>(aggregatedPath_);
    if (!sourcesFile_->open(QIODevice::WriteOnly | QIODevice::Text)) {
        report(Severity::Error, QStringLiteral("cannot open broad sources CSV"), sourcesPath_);
        return;
    }
    if (!aggregatedFile_->open(QIODevice::WriteOnly | QIODevice::Text)) {
        report(Severity::Error, QStringLiteral("cannot open broad aggregated CSV"), aggregatedPath_);
        return;
    }
    sourcesOut_ = std::make_unique<QTextStream>(sourcesFile_.get());
    aggregatedOut_ = std::make_unique<QTextStream>(aggregatedFile_.get());
    *sourcesOut_ << kSourceHeader << '\n';
    *aggregatedOut_ << kAggregatedHeader << '\n';
    ok_ = true;
}

BroadBackboneSink::~BroadBackboneSink() = default;

void BroadBackboneSink::writeSourceRow(const NeighborhoodRecord& rec, const SourceSlot& s,
                                       int64_t rowId) {
    QTextStream& out = *sourcesOut_;
    // mechanism string mirrors source_kind for join-friendly grouping.
    const char* mechanism = s.kind == SourceKind::Ring    ? "ring"
                            : s.kind == SourceKind::Bond  ? "bond"
                                                          : "charge";
    const Vec3 sourceNormalLocal = s.kind == SourceKind::Ring ? s.source_normal_local
                                                              : Vec3::Zero();
    const Vec3 bondAxisLocal = s.kind == SourceKind::Bond ? s.bond_axis_local
                                                          : Vec3::Zero();
    // Identity columns are the TARGET atom's, minimal (full identity + target
    // live ONCE on the aggregated row; this row joins back on row_id).
    // frame_variant encodes the backbone atom class (N/CA/C/O/HA frame).
    out << rowId << ',' << rec.atom_index << ',' << rec.residue_number << ','
        << rec.element << ',' << rec.atom_name << ',' << rec.frame_variant << ','
        << mechanism << ',' << static_cast<int>(s.kind) << ','
        << num(s.disp_local.x()) << ',' << num(s.disp_local.y()) << ',' << num(s.disp_local.z())
        << ',' << num(s.r) << ',' << num(s.cos_theta) << ',' << num(s.dipolar) << ','
        << s.ring_index << ',' << s.ring_type_index << ',' << num(s.ring_intensity) << ','
        << s.ring_nitrogen << ',' << (s.is_self_or_bonded ? 1 : 0) << ','
        << s.bond_index << ',' << s.bond_category << ',' << s.bond_atom_a << ','
        << s.bond_atom_b << ',' << s.source_atom_index << ',' << s.source_residue_number << ','
        << s.source_element << ',' << num(s.source_q_e) << ','
        << num(sourceNormalLocal.x()) << ',' << num(sourceNormalLocal.y()) << ','
        << num(sourceNormalLocal.z()) << ','
        << num(bondAxisLocal.x()) << ',' << num(bondAxisLocal.y()) << ','
        << num(bondAxisLocal.z());
    writeMcLitKernel(out, s.bond_mc_lit_kernel,
                     s.kind == SourceKind::Bond && s.bond_mc_lit_kernel_present);
    out << ',' << (s.kind == SourceKind::Bond && s.mc_source_is_self_or_bonded ? 1 : 0)
        << '\n';
    ++sourceRows_;
}

void BroadBackboneSink::writeAggregatedRow(const NeighborhoodRecord& rec,
                                           const BroadAggregate& agg, int64_t rowId) {
    QTextStream& out = *aggregatedOut_;
    out << rowId << ',' << rec.atom_index << ',' << rec.residue_index << ','
        << rec.residue_number << ',' << rec.amino_acid << ',' << rec.element << ','
        << rec.atom_name << ',' << rec.frame_variant << ',' << rec.h5_row << ','
        << rec.original_index << ',' << num(rec.time_ps) << ','
        << num(rec.frame_z.x()) << ',' << num(rec.frame_z.y()) << ',' << num(rec.frame_z.z()) << ','
        << num(rec.frame_x.x()) << ',' << num(rec.frame_x.y()) << ',' << num(rec.frame_x.z()) << ','
        << num(rec.frame_y.x()) << ',' << num(rec.frame_y.y()) << ',' << num(rec.frame_y.z()) << ','
        << rec.frame_variant << ',' << (rec.frame_valid ? 1 : 0) << ','
        << rec.frame_anchor_atom_index << ','
        << agg.ring_n << ',' << num(agg.ring_sum_dipolar) << ',' << num(agg.ring_cutoff_A) << ','
        << agg.bond_n << ',' << num(agg.bond_sum_dipolar) << ',' << num(agg.bond_cutoff_A);
    bool mcLitPresent = false;
    const SphericalTensor mcLit = sumMcLitKernel(rec, &mcLitPresent);
    writeMcLitKernel(out, mcLit, mcLitPresent);
    out << ','
        << agg.charge_n << ',' << agg.charge_source << ',' << num(agg.charge_cutoff_A) << ','
        << num(agg.charge_field_local.x()) << ',' << num(agg.charge_field_local.y()) << ','
        << num(agg.charge_field_local.z()) << ',' << num(agg.charge_field_z) << ','
        << num(agg.charge_field_mag) << ','
        << num(agg.charge_mu_local.x()) << ',' << num(agg.charge_mu_local.y()) << ','
        << num(agg.charge_mu_local.z());
    writeKernelT2(out, agg.literature_kernel);
    writeKernelT2(out, agg.ring_literature_kernel);
    writeKernelT2(out, agg.bond_literature_kernel);
    writeKernelT2(out, agg.charge_literature_kernel);
    out << ','
        << (rec.target.present ? 1 : 0) << ',' << num(rec.target.total_decomp.T0);
    auto m9 = [&](const Mat3& m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out << ',' << num(m(i, j));
    };
    m9(rec.target.total_raw);
    for (double v : rec.target.total_decomp.T1) out << ',' << num(v);
    for (double v : rec.target.total_decomp.T2) out << ',' << num(v);
    out << ',' << (rec.target.local_frame_valid ? 1 : 0) << ','
        << agg.bond_n_valid << ',' << num(agg.bond_sum_dipolar_valid) << ','
        << num(agg.mc_near_field_ratio);
    bool mcLitValidPresent = false;
    const SphericalTensor mcLitValid = sumValidMcLitKernel(rec, &mcLitValidPresent);
    writeMcLitKernel(out, mcLitValid, mcLitValidPresent);
    out << '\n';

    // NPY payloads, in aggregated-row order (target ONCE; field once).
    for (double v : rec.target.total_decomp.T2) aggTargetT2_.push_back(v);
    if (rec.target.present && rec.target.local_frame_valid) {
        const std::array<double, 5> lt2 = DecomposeLibrary(rec.target.total_local).T2;
        for (double v : lt2) aggTargetLocalT2_.push_back(v);
    } else {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (int i = 0; i < 5; ++i) aggTargetLocalT2_.push_back(nan);
    }
    aggFieldLocal_.push_back(agg.charge_field_local.x());
    aggFieldLocal_.push_back(agg.charge_field_local.y());
    aggFieldLocal_.push_back(agg.charge_field_local.z());
    appendKernelT2(aggLiteratureKernelT2_, agg.literature_kernel);
    appendKernelT2(aggRingLiteratureKernelT2_, agg.ring_literature_kernel);
    appendKernelT2(aggBondLiteratureKernelT2_, agg.bond_literature_kernel);
    appendKernelT2(aggChargeLiteratureKernelT2_, agg.charge_literature_kernel);
    ++aggRows_;
}

void BroadBackboneSink::Write(const NeighborhoodRecord& rec, const BroadAggregate& agg) {
    if (!ok_) return;
    const int64_t rowId = nextRowId_++;
    // Source rows FIRST (one per source, carrying only row_id + source fields).
    for (const SourceSlot& s : rec.sources) writeSourceRow(rec, s, rowId);
    // The single aggregated row + the once-per-(atom,frame) target NPY payloads.
    writeAggregatedRow(rec, agg, rowId);
}

bool BroadBackboneSink::Commit() {
    if (!ok_) return false;
    if (sourcesOut_) sourcesOut_->flush();
    if (aggregatedOut_) aggregatedOut_->flush();
    const bool a = sourcesFile_ && sourcesFile_->commit();
    const bool b = aggregatedFile_ && aggregatedFile_->commit();
    const QString outDir = QFileInfo(aggregatedPath_).dir().absolutePath();
    bool sidecarsOk = true;
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[0]),
                                {aggRows_, 5}, aggTargetT2_);
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[1]),
                                {aggRows_, 5}, aggTargetLocalT2_);
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[2]),
                                {aggRows_, 3}, aggFieldLocal_);
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[3]),
                                {aggRows_, 5}, aggLiteratureKernelT2_);
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[4]),
                                {aggRows_, 5}, aggRingLiteratureKernelT2_);
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[5]),
                                {aggRows_, 5}, aggBondLiteratureKernelT2_);
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[6]),
                                {aggRows_, 5}, aggChargeLiteratureKernelT2_);
    if (!a) report(Severity::Error, QStringLiteral("broad sources commit failed"), sourcesPath_);
    if (!b) report(Severity::Error, QStringLiteral("broad aggregated commit failed"), aggregatedPath_);
    if (!sidecarsOk)
        report(Severity::Error, QStringLiteral("broad sidecar NPY commit failed"), caseName_);
    qCInfo(cBroadSink).noquote() << "committed" << caseName_
                                 << "| source_rows=" << sourceRows_
                                 << "| agg_rows=" << aggRows_;
    return a && b && sidecarsOk;
}

}  // namespace h5reader::rediscover

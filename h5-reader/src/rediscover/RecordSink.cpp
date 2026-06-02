#include "RecordSink.h"

#include "SphericalBasis.h"

#include "../diagnostics/ErrorBus.h"

#include <QDir>
#include <QFileInfo>
#include <QIODevice>
#include <QLoggingCategory>
#include <QSaveFile>

#include <array>
#include <limits>
#include <numeric>

namespace h5reader::rediscover {

namespace {
Q_LOGGING_CATEGORY(cSink, "h5reader.rediscover.sink")

using Severity = h5reader::diagnostics::Severity;

void report(Severity sev, const QString& msg, const QString& ctx) {
    h5reader::diagnostics::ErrorBus::Report(sev, QStringLiteral("RediscoverRecordSink"), msg, ctx);
}

// Number formatting: fixed-but-honest, 9 significant digits. NaN/Inf are
// written verbatim ("nan"/"inf") so a gap is visible, never silently zeroed.
QString num(double v) { return QString::number(v, 'g', 9); }

void writeHeader(QTextStream& out, const std::vector<FeatureColumn>& cols) {
    for (std::size_t i = 0; i < cols.size(); ++i) {
        if (i) out << ',';
        out << cols[i].name;
    }
    out << '\n';
}

void appendT2(std::vector<double>& dst, const std::array<double, 5>& t2) {
    for (double v : t2) dst.push_back(v);
}

void appendLocalTargetT2(std::vector<double>& dst, const NeighborhoodRecord& rec) {
    if (!rec.target.present || !rec.target.local_frame_valid) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (int i = 0; i < 5; ++i) dst.push_back(nan);
        return;
    }
    appendT2(dst, DecomposeLibrary(rec.target.total_local).T2);
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

QString sidecarName(const FeatureSchema& schema, const QString& rowKind, const QString& payload) {
    return QStringLiteral("rediscover_%1_%2_%3.npy").arg(schema.caseName, rowKind, payload);
}

// The identity + frame + local-frame columns shared by BOTH row kinds.
void writeIdentity(QTextStream& out, const NeighborhoodRecord& r) {
    out << r.atom_index << ',' << r.residue_index << ',' << r.residue_number << ','
        << r.amino_acid << ',' << r.element << ',' << r.atom_name << ',' << r.stratum << ','
        << r.h5_row << ',' << r.original_index << ',' << num(r.time_ps) << ','
        << num(r.frame_z.x()) << ',' << num(r.frame_z.y()) << ',' << num(r.frame_z.z()) << ','
        << num(r.frame_x.x()) << ',' << num(r.frame_x.y()) << ',' << num(r.frame_x.z()) << ','
        << num(r.frame_y.x()) << ',' << num(r.frame_y.y()) << ',' << num(r.frame_y.z()) << ','
        << r.frame_variant << ',' << (r.frame_valid ? 1 : 0) << ','
        << r.frame_anchor_atom_index;
}

// The DFT-target columns shared by BOTH row kinds. raw 3x3 (total/dia/para,
// 9 each, lab frame), library-basis decomposition of total (T0 σ_iso + T1[3] +
// T2[5]), and the total tensor in the local frame (9). present flag last.
void writeTarget(QTextStream& out, const DftTarget& t) {
    auto m9 = [&](const Mat3& m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out << ',' << num(m(i, j));
    };
    out << ',' << (t.present ? 1 : 0);
    m9(t.total_raw);
    m9(t.dia_raw);
    m9(t.para_raw);
    out << ',' << num(t.total_decomp.T0);  // σ_iso
    for (double v : t.total_decomp.T1) out << ',' << num(v);
    for (double v : t.total_decomp.T2) out << ',' << num(v);
    out << ',' << num(t.dia_decomp.T0) << ',' << num(t.para_decomp.T0);
    m9(t.total_local);
    out << ',' << (t.local_frame_valid ? 1 : 0);
}

void writeBareKernel(QTextStream& out, const NeighborhoodRecord& r) {
    out << ',' << (r.bare_kernel_present ? 1 : 0) << ',' << num(r.bare_kernel.T0);
    for (double v : r.bare_kernel.T1) out << ',' << num(v);
    for (double v : r.bare_kernel.T2) out << ',' << num(v);
}

double presentValue(bool present, double value) {
    return present ? value : std::numeric_limits<double>::quiet_NaN();
}

void writeMcLitKernel(QTextStream& out, const SphericalTensor& k, bool present) {
    out << ',' << (present ? 1 : 0) << ',' << num(presentValue(present, k.T0));
    for (double v : k.T2) out << ',' << num(presentValue(present, v));
}

SphericalTensor sumMcLitKernel(const NeighborhoodRecord& r, bool* present) {
    SphericalTensor out;
    bool any = false;
    for (const SourceSlot& s : r.sources) {
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

void writeRingSource(QTextStream& out, const SourceSlot& s) {
    out << ',' << static_cast<int>(s.kind) << ','
        << num(s.disp_local.x()) << ',' << num(s.disp_local.y()) << ',' << num(s.disp_local.z()) << ','
        << num(s.r) << ',' << num(s.cos_theta) << ',' << num(s.dipolar) << ','
        << num(s.ring_z) << ',' << num(s.ring_rho) << ',' << num(s.ring_in_plane_angle) << ','
        << s.ring_index << ',' << (s.is_self_or_bonded ? 1 : 0) << ','
        << s.ring_type_index << ',' << num(s.ring_intensity) << ',' << s.ring_nitrogen << ','
        << num(s.ring_jb_offset) << ',' << s.ring_aromaticity << ',' << s.ring_size << ','
        << (s.ring_fused ? 1 : 0) << ','
        << num(s.source_normal_local.x()) << ',' << num(s.source_normal_local.y()) << ','
        << num(s.source_normal_local.z()) << ','
        << (s.ring_jb_kernel_present ? 1 : 0) << ',' << num(s.ring_jb_unit_kernel.T0);
    for (double v : s.ring_jb_unit_kernel.T2) out << ',' << num(v);
    out << ',' << num(s.ring_jb_kernel.T0);
    for (double v : s.ring_jb_kernel.T2) out << ',' << num(v);
}

void writeBondSource(QTextStream& out, const SourceSlot& s, bool includeMcLitKernel) {
    out << ',' << static_cast<int>(s.kind) << ','
        << num(s.disp_local.x()) << ',' << num(s.disp_local.y()) << ',' << num(s.disp_local.z()) << ','
        << num(s.r) << ',' << num(s.cos_theta) << ',' << num(s.dipolar) << ','
        << s.bond_category << ',' << s.bond_order << ','
        << s.bond_elem_a << ',' << s.bond_elem_b << ',' << s.bond_index << ','
        << s.bond_atom_a << ',' << s.bond_atom_b << ','
        << num(s.bond_axis_local.x()) << ',' << num(s.bond_axis_local.y()) << ','
        << num(s.bond_axis_local.z());
    if (includeMcLitKernel)
        writeMcLitKernel(out, s.bond_mc_lit_kernel, s.bond_mc_lit_kernel_present);
}

void writeChargeSource(QTextStream& out, const SourceSlot& s) {
    out << ',' << static_cast<int>(s.kind) << ',' << s.source_charge_source << ','
        << s.source_atom_index << ',' << s.source_residue_index << ','
        << s.source_residue_number << ',' << s.source_amino_acid << ','
        << s.source_element << ',' << s.source_atom_name << ','
        << num(s.source_q_e) << ','
        << num(s.disp_local.x()) << ',' << num(s.disp_local.y()) << ','
        << num(s.disp_local.z()) << ',' << num(s.r);
}

}  // namespace

RecordSink::RecordSink(const QString& outDir, const FeatureSchema& schema) : schema_(schema) {
    QDir().mkpath(outDir);
    sourcesPath_ = QStringLiteral("%1/%2_sources.csv").arg(outDir, schema.caseName);
    aggregatedPath_ = QStringLiteral("%1/%2_aggregated.csv").arg(outDir, schema.caseName);
    sidecarFiles_ = {
        sidecarName(schema, QStringLiteral("sources"), QStringLiteral("target_T2")),
        sidecarName(schema, QStringLiteral("sources"), QStringLiteral("target_local_T2")),
        sidecarName(schema, QStringLiteral("aggregated"), QStringLiteral("target_T2")),
        sidecarName(schema, QStringLiteral("aggregated"), QStringLiteral("target_local_T2")),
    };
    if (schema_.includeBareKernel) {
        sidecarFiles_.insert(2, sidecarName(schema, QStringLiteral("sources"),
                                           QStringLiteral("bare_kernel_T2")));
        sidecarFiles_.append(sidecarName(schema, QStringLiteral("aggregated"),
                                        QStringLiteral("bare_kernel_T2")));
    }

    sourcesFile_ = std::make_unique<QSaveFile>(sourcesPath_);
    aggregatedFile_ = std::make_unique<QSaveFile>(aggregatedPath_);
    if (!sourcesFile_->open(QIODevice::WriteOnly | QIODevice::Text)) {
        report(Severity::Error, QStringLiteral("cannot open sources CSV"), sourcesPath_);
        return;
    }
    if (!aggregatedFile_->open(QIODevice::WriteOnly | QIODevice::Text)) {
        report(Severity::Error, QStringLiteral("cannot open aggregated CSV"), aggregatedPath_);
        return;
    }
    sourcesOut_ = std::make_unique<QTextStream>(sourcesFile_.get());
    aggregatedOut_ = std::make_unique<QTextStream>(aggregatedFile_.get());
    writeHeader(*sourcesOut_, schema_.sourceColumns);
    writeHeader(*aggregatedOut_, schema_.aggregatedColumns);
    ok_ = true;
}

RecordSink::~RecordSink() = default;

void RecordSink::WriteSourceRows(const NeighborhoodRecord& rec) {
    if (!ok_) return;
    QTextStream& out = *sourcesOut_;
    for (const SourceSlot& s : rec.sources) {
        writeIdentity(out, rec);
        if (schema_.sourceSchemaKind == SourceSchemaKind::Ring)
            writeRingSource(out, s);
        else if (schema_.sourceSchemaKind == SourceSchemaKind::Bond)
            writeBondSource(out, s, schema_.includeMcLitKernel);
        else if (schema_.sourceSchemaKind == SourceSchemaKind::Charge)
            writeChargeSource(out, s);
        if (schema_.includeBareKernel) writeBareKernel(out, rec);
        writeTarget(out, rec.target);
        out << '\n';
        appendT2(sourceTargetT2_, rec.target.total_decomp.T2);
        appendLocalTargetT2(sourceTargetLocalT2_, rec);
        if (schema_.includeBareKernel) appendT2(sourceBareKernelT2_, rec.bare_kernel.T2);
        ++sourceRows_;
    }
}

void RecordSink::WriteAggregatedRow(const NeighborhoodRecord& rec, double sumAll, double sumValid,
                                    int nSourcesValid, const std::vector<double>& perTypeSums,
                                    double cutoff_A) {
    if (!ok_) return;
    QTextStream& out = *aggregatedOut_;
    writeIdentity(out, rec);
    // n_sources (all), n_sources_valid (producer-valid: self/bonded excluded for
    // the ring case), sum_dipolar_all, sum_dipolar_producer_valid, per-type sums
    // (on the producer-valid set), then the recorded source-discovery cutoff
    // (NaN when the sources come from the producer's H5, whose cutoff we don't set).
    out << ',' << static_cast<int>(rec.sources.size()) << ',' << nSourcesValid
        << ',' << num(sumAll) << ',' << num(sumValid);
    for (double v : perTypeSums) out << ',' << num(v);
    out << ',' << num(cutoff_A);
    if (schema_.includeMcLitKernel) {
        bool mcLitPresent = false;
        const SphericalTensor mcLit = sumMcLitKernel(rec, &mcLitPresent);
        writeMcLitKernel(out, mcLit, mcLitPresent);
    }
    if (schema_.includeBareKernel) writeBareKernel(out, rec);
    writeTarget(out, rec.target);
    out << '\n';
    appendT2(aggregatedTargetT2_, rec.target.total_decomp.T2);
    appendLocalTargetT2(aggregatedTargetLocalT2_, rec);
    if (schema_.includeBareKernel) appendT2(aggregatedBareKernelT2_, rec.bare_kernel.T2);
    ++aggRows_;
}

void RecordSink::WriteChargeDipoleAggregatedRow(const NeighborhoodRecord& rec, const Vec3& mu,
                                                double cutoff_A, const QString& chargeSource,
                                                bool excludeResidue) {
    if (!ok_) return;
    QTextStream& out = *aggregatedOut_;
    writeIdentity(out, rec);
    out << ',' << static_cast<int>(rec.sources.size()) << ',' << chargeSource << ','
        << (excludeResidue ? 1 : 0) << ',' << num(cutoff_A) << ','
        << num(mu.x()) << ',' << num(mu.y()) << ',' << num(mu.z()) << ','
        << num(mu.norm());
    if (schema_.includeBareKernel) writeBareKernel(out, rec);
    writeTarget(out, rec.target);
    out << '\n';
    appendT2(aggregatedTargetT2_, rec.target.total_decomp.T2);
    appendLocalTargetT2(aggregatedTargetLocalT2_, rec);
    if (schema_.includeBareKernel) appendT2(aggregatedBareKernelT2_, rec.bare_kernel.T2);
    ++aggRows_;
}

bool RecordSink::Commit() {
    if (!ok_) return false;
    if (sourcesOut_) sourcesOut_->flush();
    if (aggregatedOut_) aggregatedOut_->flush();
    const bool a = sourcesFile_ && sourcesFile_->commit();
    const bool b = aggregatedFile_ && aggregatedFile_->commit();
    const QString outDir = QFileInfo(sourcesPath_).dir().absolutePath();
    bool sidecarsOk = true;
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[0]),
                                {sourceRows_, 5}, sourceTargetT2_);
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[1]),
                                {sourceRows_, 5}, sourceTargetLocalT2_);
    std::size_t aggOffset = 2;
    if (schema_.includeBareKernel) {
        sidecarsOk = sidecarsOk
                     && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[2]),
                                    {sourceRows_, 5}, sourceBareKernelT2_);
        aggOffset = 3;
    }
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[aggOffset]),
                                {aggRows_, 5}, aggregatedTargetT2_);
    sidecarsOk = sidecarsOk
                 && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[aggOffset + 1]),
                                {aggRows_, 5}, aggregatedTargetLocalT2_);
    if (schema_.includeBareKernel) {
        sidecarsOk = sidecarsOk
                     && writeNpyF64(QStringLiteral("%1/%2").arg(outDir, sidecarFiles_[aggOffset + 2]),
                                    {aggRows_, 5}, aggregatedBareKernelT2_);
    }
    if (!a) report(Severity::Error, QStringLiteral("commit failed"), sourcesPath_);
    if (!b) report(Severity::Error, QStringLiteral("commit failed"), aggregatedPath_);
    if (!sidecarsOk)
        report(Severity::Error, QStringLiteral("sidecar NPY commit failed"), schema_.caseName);
    qCInfo(cSink).noquote() << "committed" << schema_.caseName << "| source_rows=" << sourceRows_
                            << "| agg_rows=" << aggRows_;
    return a && b && sidecarsOk;
}

}  // namespace h5reader::rediscover

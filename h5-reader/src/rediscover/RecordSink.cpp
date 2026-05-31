#include "RecordSink.h"

#include "../diagnostics/ErrorBus.h"

#include <QDir>
#include <QIODevice>
#include <QLoggingCategory>

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

}  // namespace

RecordSink::RecordSink(const QString& outDir, const FeatureSchema& schema) : schema_(schema) {
    QDir().mkpath(outDir);
    sourcesPath_ = QStringLiteral("%1/%2_sources.csv").arg(outDir, schema.caseName);
    aggregatedPath_ = QStringLiteral("%1/%2_aggregated.csv").arg(outDir, schema.caseName);

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
        // Per-source geometry + identity.
        out << ',' << static_cast<int>(s.kind) << ','
            << num(s.disp_local.x()) << ',' << num(s.disp_local.y()) << ',' << num(s.disp_local.z()) << ','
            << num(s.r) << ',' << num(s.cos_theta) << ',' << num(s.dipolar) << ','
            << num(s.ring_z) << ',' << num(s.ring_rho) << ',' << num(s.ring_in_plane_angle) << ','
            << s.ring_index << ',' << (s.is_self_or_bonded ? 1 : 0) << ','
            << s.ring_type_index << ',' << num(s.ring_intensity) << ',' << s.ring_nitrogen << ','
            << num(s.ring_jb_offset) << ',' << s.ring_aromaticity << ',' << s.ring_size << ','
            << (s.ring_fused ? 1 : 0) << ','
            << s.bond_category << ',' << s.bond_order << ','
            << s.bond_elem_a << ',' << s.bond_elem_b << ',' << s.bond_index << ','
            << s.bond_atom_a << ',' << s.bond_atom_b << ','
            << num(s.bond_axis_local.x()) << ',' << num(s.bond_axis_local.y()) << ','
            << num(s.bond_axis_local.z()) << ','
            << num(s.source_normal_local.x()) << ',' << num(s.source_normal_local.y()) << ','
            << num(s.source_normal_local.z());
        writeBareKernel(out, rec);
        writeTarget(out, rec.target);
        out << '\n';
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
    writeBareKernel(out, rec);
    writeTarget(out, rec.target);
    out << '\n';
    ++aggRows_;
}

bool RecordSink::Commit() {
    if (!ok_) return false;
    if (sourcesOut_) sourcesOut_->flush();
    if (aggregatedOut_) aggregatedOut_->flush();
    const bool a = sourcesFile_ && sourcesFile_->commit();
    const bool b = aggregatedFile_ && aggregatedFile_->commit();
    if (!a) report(Severity::Error, QStringLiteral("commit failed"), sourcesPath_);
    if (!b) report(Severity::Error, QStringLiteral("commit failed"), aggregatedPath_);
    qCInfo(cSink).noquote() << "committed" << schema_.caseName << "| source_rows=" << sourceRows_
                            << "| agg_rows=" << aggRows_;
    return a && b;
}

}  // namespace h5reader::rediscover

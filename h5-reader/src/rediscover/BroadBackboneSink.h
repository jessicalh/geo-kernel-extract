// BroadBackboneSink — the heterogeneous two-kind carrier for the broad-backbone
// relationship, applying the target-repeat FIX (BROAD_BACKBONE_NEXT.md carrier
// fix; the charge_dipole per-source CSV bloated to 828 MB because the ~50-column
// DFT target repeated on EVERY source row).
//
// The fix: the DFT target + shared identity live ONCE per (atom, frame), on the
// AGGREGATED row + the target NPYs. The per-source rows carry only the source
// fields + an (atom, frame) `row_id` to join back. This is a NEW sink — it does
// NOT touch RecordSink, so the ring/mc composed path + the procedural cells keep
// their exact byte-parity oracle.
//
// Files written under <outDir>:
//   broad_backbone_aggregated.csv — one row per (atom, frame): row_id + identity
//     + frame basis + per-mechanism summed scalar features + the local Coulomb
//     FIELD (E in the local frame + axis-projected E_z) + μ + the full DFT
//     target (raw 3×3, T0=σ_iso, T1, T2, local). The well-posed scalar/ridge
//     input AND the once-per-(atom,frame) home of the target.
//   broad_backbone_sources.csv — one row per (atom, frame, source): row_id +
//     mechanism + source geometry (disp_local, r, cosθ, dipolar) + per-mechanism
//     identity (ring / bond / charge). NO target columns.
//   NPY payloads (one ArraySpec each in python/nmr_extract/_catalog.py):
//     broad_backbone_aggregated_target_T2.npy        (agg_rows, 5)
//     broad_backbone_aggregated_target_local_T2.npy  (agg_rows, 5)
//     broad_backbone_aggregated_field_local.npy      (agg_rows, 3)  Coulomb E
//   CSV-only mc_lit_* columns are documented locally in
//     src/rediscover/MC_LIT_CSV_CONTRACT.json
//
// Plain data class, no QObject. Streams rows (QSaveFile atomic-on-commit).

#pragma once

#include "RediscoverTypes.h"  // NeighborhoodRecord, SourceSlot, Vec3 (via model/Types.h)

#include <QSaveFile>
#include <QString>
#include <QStringList>
#include <QTextStream>

#include <cmath>
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace h5reader::rediscover {

struct BroadKernelT2 {
    bool present = false;
    std::array<double, 5> T2 = {};
};

// The per-mechanism aggregated scalar features broad-backbone emits, plus the
// field vector. One BroadAggregate per (atom, frame); the reducer builds it.
struct BroadAggregate {
    // Ring mechanism (KD ring-centres): Σ (3cos²θ−1)/r³ over discovered rings.
    double ring_sum_dipolar = 0.0;
    int    ring_n = 0;
    double ring_cutoff_A = 0.0;

    // Bond mechanism (KD bond-midpoints): Σ (3cos²θ−1)/r³ over aniso bonds.
    double bond_sum_dipolar = 0.0;
    int    bond_n = 0;
    double bond_sum_dipolar_valid = 0.0;
    int    bond_n_valid = 0;
    double bond_cutoff_A = 0.0;
    double mc_near_field_ratio = 0.0;

    // Charge mechanism (KD charge-sites, FF14SB): the local Coulomb FIELD, not μ.
    // E = Σ q_i (r_atom − r_i) / |r|³  (V·e/Å² style; units cancel in a corr fit,
    // recorded in the manifest). Expressed in the target local frame. field_z is
    // the axis-projected scalar feature (the Buckingham r=+0.46 signal). μ is
    // KEPT for the tensor story (charge_dipole found μ null vs σ_iso).
    Vec3   charge_field_local = Vec3::Zero();  // local-frame E
    double charge_field_z = 0.0;               // E·ẑ (axis projection)
    double charge_field_mag = 0.0;             // |E|
    Vec3   charge_mu_local = Vec3::Zero();     // Σ q_i (r_i − r_atom), local frame
    int    charge_n = 0;
    double charge_cutoff_A = 0.0;
    QString charge_source;

    // Fixed-coefficient T2 kernels for the de-circularising gate. Ring and bond
    // are read from the producer H5 kernels; charge is the broad FF14SB Coulomb
    // EFG built in this C++ reducer from the already-attached charge sources.
    BroadKernelT2 literature_kernel;
    BroadKernelT2 ring_literature_kernel;
    BroadKernelT2 bond_literature_kernel;
    BroadKernelT2 charge_literature_kernel;
};

inline BroadAggregate ReduceBroadBackboneSources(const std::vector<SourceSlot>& sources,
                                                 double ring_cutoff_A,
                                                 double bond_cutoff_A,
                                                 double charge_cutoff_A,
                                                 const QString& charge_source,
                                                 double mc_near_field_ratio) {
    BroadAggregate agg;
    agg.ring_cutoff_A = ring_cutoff_A;
    agg.bond_cutoff_A = bond_cutoff_A;
    agg.charge_cutoff_A = charge_cutoff_A;
    agg.charge_source = charge_source;
    agg.mc_near_field_ratio = mc_near_field_ratio;

    Vec3 field = Vec3::Zero();  // Σ q_i (r_atom − r_i)/r³, in the LOCAL frame
    Vec3 mu = Vec3::Zero();     // Σ q_i (r_i − r_atom), in the LOCAL frame
    for (const SourceSlot& s : sources) {
        switch (s.kind) {
        case SourceKind::Ring:
            if (std::isfinite(s.dipolar)) {
                agg.ring_sum_dipolar += s.dipolar;
                ++agg.ring_n;
            }
            break;
        case SourceKind::Bond:
            if (std::isfinite(s.dipolar)) {
                agg.bond_sum_dipolar += s.dipolar;
                ++agg.bond_n;
                if (!s.mc_source_is_self_or_bonded) {
                    agg.bond_sum_dipolar_valid += s.dipolar;
                    ++agg.bond_n_valid;
                }
            }
            break;
        case SourceKind::Charge:
            // disp_local = (r_i − r_atom) in the local frame; r = |disp|.
            // FIELD E = Σ q_i (r_atom − r_i)/r³ = −Σ q_i disp_local / r³.
            // μ = Σ q_i (r_i − r_atom) = Σ q_i disp_local.
            if (s.r > 1e-9 && std::isfinite(s.source_q_e)) {
                const double r3 = s.r * s.r * s.r;
                field += (-s.source_q_e / r3) * s.disp_local;
                mu += s.source_q_e * s.disp_local;
                ++agg.charge_n;
            }
            break;
        }
    }
    agg.charge_field_local = field;
    agg.charge_field_z = field.z();  // axis-projected (the Buckingham signal)
    agg.charge_field_mag = field.norm();
    agg.charge_mu_local = mu;
    return agg;
}

class BroadBackboneSink {
public:
    BroadBackboneSink(const QString& outDir, const QString& caseName);
    ~BroadBackboneSink();

    bool Ok() const { return ok_; }

    // One (atom, frame) record. The aggregated row + the target NPYs carry the
    // target ONCE; the per-source rows carry source fields + this row_id only.
    void Write(const NeighborhoodRecord& rec, const BroadAggregate& agg);

    bool Commit();

    std::size_t aggregatedRowsWritten() const { return aggRows_; }
    std::size_t sourceRowsWritten() const { return sourceRows_; }
    const QStringList& sidecarFiles() const { return sidecarFiles_; }

private:
    void writeSourceRow(const NeighborhoodRecord& rec, const SourceSlot& s, int64_t rowId);
    void writeAggregatedRow(const NeighborhoodRecord& rec, const BroadAggregate& agg,
                            int64_t rowId);

    QString caseName_;
    QString sourcesPath_;
    QString aggregatedPath_;
    QStringList sidecarFiles_;
    std::unique_ptr<QSaveFile> sourcesFile_;
    std::unique_ptr<QSaveFile> aggregatedFile_;
    std::unique_ptr<QTextStream> sourcesOut_;
    std::unique_ptr<QTextStream> aggregatedOut_;
    bool ok_ = false;
    std::size_t aggRows_ = 0;
    std::size_t sourceRows_ = 0;
    int64_t nextRowId_ = 0;

    std::vector<double> aggTargetT2_;       // (aggRows, 5)
    std::vector<double> aggTargetLocalT2_;  // (aggRows, 5)
    std::vector<double> aggFieldLocal_;     // (aggRows, 3)
    std::vector<double> aggLiteratureKernelT2_;       // (aggRows, 5)
    std::vector<double> aggRingLiteratureKernelT2_;   // (aggRows, 5)
    std::vector<double> aggBondLiteratureKernelT2_;   // (aggRows, 5)
    std::vector<double> aggChargeLiteratureKernelT2_; // (aggRows, 5)
};

}  // namespace h5reader::rediscover

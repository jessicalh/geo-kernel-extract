// QtBondVectorBuffers.h — bond-vector-axis typed buffers.
//
// Bond vectors are a first-class axis sibling to atom / residue / ring
// (per the residue-grouping ergonomic — see DashboardSignal.h's
// BondVectorAnchor + AnchorMatchesAxis widening). Two producer TRs key
// to this axis: IRedOrderParameterTrajectoryResult (M N-H vectors) and
// ReorientationalDynamicsTrajectoryResult (V backbone vectors covering
// N-H, Cα-Hα, C=O). Each TR has its OWN identity table — IRed enumerates
// only N-H so its row count differs from Reorient's. The semantic
// (residue, kind) anchor is invariant; the per-TR sampler resolves it
// to the table's row at lookup time.
//
// vector_kind enum matches the producer convention exactly:
//   1 = NH, 2 = CaHa, 3 = CO. 0 = unspecified / any.
//
// Layout: row-major per-vector identity arrays + an optional scalar
// payload per row (s2, τ_e, R1, R2, NOE, eigenvalue, etc.). Static
// (per-trajectory), not per-frame.

#pragma once

#include <QString>

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace h5reader::model {

// Per-vector identity table. Owned by each TR-specific buffer that
// keys to the bond-vector axis. Row count differs per TR.
struct QtBondVectorTable {
    std::size_t n_vectors = 0;
    std::vector<std::uint8_t> kind;          // 1=NH, 2=CaHa, 3=CO
    std::vector<std::int32_t> residue_index; // parent residue index per row
    std::vector<std::int32_t> owning_atom;   // primary atom for highlighting (H, HA, O)
    std::vector<std::int32_t> tail_atom;     // start of the vector (N, CA, C)
    std::vector<std::int32_t> head_atom;     // end of the vector   (H, HA, O)

    // Resolve a (residue, kind) anchor to a row index. Returns nullopt
    // when no matching row exists. `wantKind == 0` matches the first row
    // at `residue`, irrespective of kind.
    std::optional<std::size_t> rowFor(std::size_t residue, std::uint8_t wantKind) const {
        for (std::size_t i = 0; i < n_vectors; ++i) {
            if (residue_index[i] != static_cast<std::int32_t>(residue))
                continue;
            if (wantKind == 0 || wantKind == kind[i])
                return i;
        }
        return std::nullopt;
    }
};

// One scalar per bond vector. Used for IRed S², Reorient S² / τ_e /
// R1 / R2 / NOE / etc. The per-row values are read out by a sampler
// against an owning identity table.
struct QtPerBondVectorScalar {
    std::size_t n_vectors = 0;
    std::vector<double> values;  // (n_vectors,)
    QString units;
    QString result_name;
    QString dataset_name;

    double at(std::size_t row) const {
        return (row < n_vectors && row < values.size()) ? values[row] : 0.0;
    }
};

// Per-bond-vector 1D curve over an external grid (Reorient body / lab
// TCFs over lag). Single-curve-per-row (no channel axis), in contrast
// to KernelDynamics's (atom, channel, sample) shape.
struct QtPerBondVectorCurve {
    std::size_t n_vectors = 0;
    std::size_t n_samples = 0;
    std::vector<double> data;          // (n_vectors * n_samples,) row-major
    std::vector<double> axis_values;   // (n_samples,) — lag times
    QString axis_unit;                 // "ps"
    QString units;                     // "dimensionless" for second-rank TCF
    QString result_name;

    double at(std::size_t row, std::size_t sample) const {
        if (row >= n_vectors || sample >= n_samples || data.empty())
            return 0.0;
        return data[row * n_samples + sample];
    }
};

// Per-bond-vector 3x3 dense tensor (Reorient body-frame <u⊗u>).
// Symmetric; stored as row-major (V, 9).
struct QtPerBondVectorMat3 {
    std::size_t n_vectors = 0;
    std::vector<double> data;          // (n_vectors * 9,) row-major
    QString result_name;

    // Read the 3x3 row into a flat 9-element array.
    std::array<double, 9> at(std::size_t row) const {
        std::array<double, 9> out{};
        if (row >= n_vectors || data.size() < (row + 1) * 9)
            return out;
        for (int i = 0; i < 9; ++i) out[i] = data[row * 9 + i];
        return out;
    }
};

// Per-bond-vector K-vector sampled at K externally-fixed frequencies
// (Reorient J(ω) at the 5 KTB Larmor combinations). Frequencies are
// descriptor metadata, not per-row.
struct QtPerBondVectorFixedFreqBlock {
    std::size_t n_vectors = 0;
    std::size_t n_freqs = 0;
    std::vector<double> data;          // (n_vectors * n_freqs,) row-major
    std::vector<double> freq_values;   // (n_freqs,) rad/s
    QString units;                     // "s" for J(ω)
    QString result_name;

    double at(std::size_t row, std::size_t freqIdx) const {
        if (row >= n_vectors || freqIdx >= n_freqs || data.empty())
            return 0.0;
        return data[row * n_freqs + freqIdx];
    }
};

// IRed order parameters: M N-H vectors + per-vector S² + per-trajectory
// eigenvalue spectrum. Composite — owns its own identity table plus
// the S² scalar payload plus the descending eigenvalue list. The
// eigenvalues are a per-trajectory diagnostic (the 5 tumbling modes
// followed by M−5 internal modes); they're not per-vector signals
// despite sharing the M axis.
struct QtIRedOrderParameters {
    QtBondVectorTable identity;            // identity.n_vectors = M, all kind=NH
    std::vector<double> s2_ired;           // (M,) per-vector S²
    std::vector<double> eigenvalues;       // (M,) descending; first 5 are tumbling
    double separability_gap = 0.0;         // λ5 / λ6; NaN if M ≤ 5
    std::size_t n_frames = 0;
    QString reference;                     // "none_isotropic"
    QString frame;                         // "lab"
    QString vector_set;                    // "amide_N_H"
    QString result_name;
};

// ReorientationalDynamics composite: V backbone bond vectors (NH + CaHa
// + CO per residue, gated by atom presence). Identity shared by every
// per-vector signal; tau_m_ps is a single per-trajectory diagnostic
// (overall tumbling), with tau_m_converged flagging whether the run
// was long enough (~5 tau_m) for the rate-equation outputs to be
// trustworthy.
struct QtReorientationalDynamics {
    QtBondVectorTable identity;                       // n_vectors = V
    QtPerBondVectorCurve  acf_internal;               // body-frame TCF
    QtPerBondVectorCurve  acf_lab;                    // lab-frame TCF
    QtPerBondVectorScalar s2;                         // Henry-Szabo S²
    QtPerBondVectorScalar tau_e;                      // Lipari-Szabo τ_e (ps)
    QtPerBondVectorMat3   orientation_tensor;        // body-frame <u⊗u>
    QtPerBondVectorFixedFreqBlock spectral_density_J; // J at 5 KTB freqs
    QtPerBondVectorScalar r1;                         // (V,) s⁻¹, NH only
    QtPerBondVectorScalar r2;                         // (V,) s⁻¹, NH only
    QtPerBondVectorScalar noe;                        // (V,) NH only
    double tau_m_ps = 0.0;
    double trajectory_length_over_tau_m = 0.0;
    bool tau_m_converged = false;
    double relaxation_field_tesla = 0.0;
    QString result_name;
};

}  // namespace h5reader::model

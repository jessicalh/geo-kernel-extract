// QtPerResidueBuffers.h — residue-axis per-frame typed buffers.
//
// These TRs use (R, T, K) shape instead of (N, T, K). The reader's
// access pattern: walk residues (54 for 1P9J), not atoms. Mainly used
// by the dihedral / DSSP / ring-pucker docks.
//
// DSSP carries a complex shape: ss8 [R, T] uint8 + four hbond
// datasets [R, T, 2] (top-2 acceptor/donor partners + energies).
// Dihedral carries phi/psi/omega/omega_deviation + chi (R, T, 4) +
// rama_region + static per-residue markers (chi_exists, is_proline,
// is_glycine, is_pre_proline, omega_is_xpro, residue_terminal_state).

#pragma once

#include "Types.h"

#include <QString>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace h5reader::model {

// Common metadata for residue-axis TRs (mirrors the atom-axis shape).
struct QtPerResidueFrameMeta {
    std::vector<uint64_t> frame_indices;
    std::vector<double> frame_times;
    std::vector<uint8_t> source_attached;

    std::size_t frameCount() const { return frame_times.size(); }
    bool sourceAttachedAt(std::size_t t) const {
        return source_attached.empty() || (t < source_attached.size() && source_attached[t] != 0);
    }
};


// ──────────────────────────────────────────────────────────────────
// QtPerResidueScalarTimeSeries — (R, T) per-residue scalar TS.
// Generic shape used by ring_pucker_time_series and potential future
// per-residue scalars.
// ──────────────────────────────────────────────────────────────────

struct QtPerResidueScalarTimeSeries {
    std::size_t n_residues = 0;
    std::size_t n_frames = 0;

    std::vector<double> data;  // (R*T,) row-major
    QtPerResidueFrameMeta meta;

    QString units;
    QString result_name;
    QString dataset_name;

    double at(std::size_t resIdx, std::size_t t) const {
        if (resIdx >= n_residues || t >= n_frames)
            return 0.0;
        return data[resIdx * n_frames + t];
    }
    bool sourceAttachedAt(std::size_t t) const { return meta.sourceAttachedAt(t); }
};


// ──────────────────────────────────────────────────────────────────
// QtDihedralTimeSeries — /trajectory/dihedral_time_series/ shape.
//
// Per-residue per-frame:
//   - phi, psi, omega, omega_deviation  (R, T) float64 radians
//   - chi[k=0..3]                       (R, T, 4) float64 radians
//   - rama_region                       (R, T) uint8 enum
//     {0:unassigned, 1:alphaR, 2:beta, 3:alphaL, 4:PPII, 5:other}
//
// Per-residue static (per-residue boolean / enum, no T axis):
//   - chi_exists                        (R, 4) uint8 (per chi slot)
//   - is_proline, is_glycine, is_pre_proline, omega_is_xpro (R,) uint8
//   - residue_terminal_state            (R,) uint8 (TerminalState ordinal)
//
// NOTE on chain_id_per_residue (object/string array in the H5): we
// IGNORE that string column at load — residues.npy.chain_id is the
// canonical residue addressing (via QtResidue.address). Reading the
// in-TR string would be Snobol-redundant.
// ──────────────────────────────────────────────────────────────────

enum class RamaRegion : uint8_t {
    Unassigned = 0,
    AlphaR = 1,
    Beta = 2,
    AlphaL = 3,
    PPII = 4,
    Other = 5,
};

struct QtDihedralTimeSeries {
    std::size_t n_residues = 0;
    std::size_t n_frames = 0;

    // Per-residue per-frame
    std::vector<double> phi;              // (R*T,)
    std::vector<double> psi;              // (R*T,)
    std::vector<double> omega;            // (R*T,)
    std::vector<double> omega_deviation;  // (R*T,)
    std::vector<double> chi;              // (R*T*4,)
    std::vector<uint8_t> rama_region;     // (R*T,) → RamaRegion

    // Per-residue static
    std::vector<uint8_t> chi_exists;              // (R*4,)
    std::vector<uint8_t> is_proline;              // (R,)
    std::vector<uint8_t> is_glycine;              // (R,)
    std::vector<uint8_t> is_pre_proline;          // (R,)
    std::vector<uint8_t> omega_is_xpro;           // (R,)
    std::vector<uint8_t> residue_terminal_state;  // (R,) → TerminalState ordinal

    QtPerResidueFrameMeta meta;

    QString units;
    QString result_name;
    QString angle_convention;  // "IUPAC signed dihedral atan2(y,x)..."

    double phiAt(std::size_t r, std::size_t t) const { return (r < n_residues && t < n_frames) ? phi[r * n_frames + t] : 0.0; }
    double psiAt(std::size_t r, std::size_t t) const { return (r < n_residues && t < n_frames) ? psi[r * n_frames + t] : 0.0; }
    double omegaAt(std::size_t r, std::size_t t) const {
        return (r < n_residues && t < n_frames) ? omega[r * n_frames + t] : 0.0;
    }
    double chiAt(std::size_t r, std::size_t t, int k) const {
        if (r >= n_residues || t >= n_frames || k < 0 || k > 3)
            return 0.0;
        return chi[(r * n_frames + static_cast<std::size_t>(t)) * 4 + static_cast<std::size_t>(k)];
    }
    RamaRegion ramaAt(std::size_t r, std::size_t t) const {
        if (r >= n_residues || t >= n_frames)
            return RamaRegion::Unassigned;
        return static_cast<RamaRegion>(rama_region[r * n_frames + t]);
    }
    bool sourceAttachedAt(std::size_t t) const { return meta.sourceAttachedAt(t); }
};


// ──────────────────────────────────────────────────────────────────
// QtDssp8TimeSeries — /trajectory/dssp8_time_series/ shape.
//
// Per-residue per-frame:
//   - ss8_code            (R, T) uint8 → DsspCode enum (255 = absent)
//   - hbond_acceptor_partner / energy   (R, T, 2) int32 / float64
//   - hbond_donor_partner / energy      (R, T, 2) int32 / float64
//
// Per-atom cross-reference: residue_index_per_atom (N,) int32 — lets
// the inspector look up which residue an atom belongs to even when
// the user is browsing per-atom data (it duplicates QtAtom.residueIndex
// but is provided by the TR for self-contained access).
//
// Per-frame mask: source_attached_per_frame (T,) uint8 — DSSP is
// conditional-attach (skip_dssp = true frames have ss8_code = 255 +
// hbond_partner = -1 + hbond_energy = NaN).
// ──────────────────────────────────────────────────────────────────

struct QtDssp8TimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_residues = 0;
    std::size_t n_frames = 0;

    std::vector<uint8_t> ss8_code;                // (R*T,)
    std::vector<int32_t> hbond_acceptor_partner;  // (R*T*2,)
    std::vector<double> hbond_acceptor_energy;    // (R*T*2,) kcal/mol
    std::vector<int32_t> hbond_donor_partner;     // (R*T*2,)
    std::vector<double> hbond_donor_energy;       // (R*T*2,) kcal/mol
    std::vector<int32_t> residue_index_per_atom;  // (N,)

    QtPerResidueFrameMeta meta;

    QString result_name;
    QString ss8_legend;       // display only
    QString hbond_threshold;  // narrative description

    DsspCode codeAt(std::size_t resIdx, std::size_t t) const {
        if (resIdx >= n_residues || t >= n_frames)
            return DsspCode::Unknown;
        return static_cast<DsspCode>(ss8_code[resIdx * n_frames + t]);
    }
    bool sourceAttachedAt(std::size_t t) const { return meta.sourceAttachedAt(t); }
};


// ──────────────────────────────────────────────────────────────────
// QtJCouplingTimeSeries — per-residue NMR scalar couplings.
// Source: /trajectory/j_coupling_time_series/.
//
// 9 J-coupling channels, each (R, T) float64 (Hz). Each has an
// existence mask (R,) uint8 indicating which residues actually
// have the geometry to compute it (e.g., J_HN_Halpha exists only
// where both atoms are present; chi-related Js depend on chi
// definitions per AminoAcid).
//
// Cross-reference: residue_index_per_atom (N,) int32 — same shape
// as DSSP's; lets inspector code map an atom index to its residue
// for J display.
//
// Always-attached (J couplings are geometry-derived, no DFT gate).
// ──────────────────────────────────────────────────────────────────

struct QtJCouplingTimeSeries {
    std::size_t n_residues = 0;
    std::size_t n_frames = 0;

    // 9 J channels (R*T,) Hz
    std::vector<double> J_Cprime_Cgamma;
    std::vector<double> J_HN_Cbeta;
    std::vector<double> J_HN_Cprime;
    std::vector<double> J_HN_Halpha;
    std::vector<double> J_HN_Halpha_Vogeli;
    std::vector<double> J_Halpha_Cprime;
    std::vector<double> J_Halpha_Hbeta2;
    std::vector<double> J_Halpha_Hbeta3;
    std::vector<double> J_N_Cgamma;

    // Per-residue existence masks (R,)
    std::vector<uint8_t> J_Cprime_Cgamma_exists;
    std::vector<uint8_t> J_HN_Cbeta_exists;
    std::vector<uint8_t> J_HN_Cprime_exists;
    std::vector<uint8_t> J_HN_Halpha_exists;
    std::vector<uint8_t> J_Halpha_Cprime_exists;
    std::vector<uint8_t> J_Halpha_Hbeta_exists;
    std::vector<uint8_t> J_N_Cgamma_exists;
    std::vector<uint8_t> J_chi1_exists;

    // Per-atom cross-reference
    std::vector<int32_t> residue_index_per_atom;

    QtPerResidueFrameMeta meta;

    QString result_name;
    QString units;

    double at(const std::vector<double>& channel, std::size_t r, std::size_t t) const {
        if (r >= n_residues || t >= n_frames)
            return 0.0;
        return channel[r * n_frames + t];
    }
};



// ──────────────────────────────────────────────────────────────────
// QtDihedralAutocorrelation — /trajectory/dihedral_autocorrelation/
//
// Per-residue circular ACF of phi/psi/chi torsions + 1/e decorrelation
// times. v1 surfaces phi + psi only; chi[0..3] are present in the H5
// for future expansion. Static (no time axis on the curves).
// ──────────────────────────────────────────────────────────────────

struct QtPerResidueCurve {
    std::size_t n_residues = 0;
    std::size_t n_samples = 0;
    std::vector<double> data;          // (R * n_samples,) row-major
    std::vector<double> axis_values;   // (n_samples,) lag times
    QString axis_unit;
    QString units;
    QString result_name;

    double at(std::size_t residue, std::size_t sample) const {
        if (residue >= n_residues || sample >= n_samples || data.empty())
            return 0.0;
        return data[residue * n_samples + sample];
    }
};

struct QtPerResidueScalar {
    std::size_t n_residues = 0;
    std::vector<double> values;        // (R,)
    std::vector<uint8_t> defined;      // (R,) — 1 if structurally defined
    QString units;
    QString result_name;

    double at(std::size_t residue) const {
        return (residue < n_residues && residue < values.size())
                   ? values[residue] : 0.0;
    }
    bool isDefined(std::size_t residue) const {
        return (residue < defined.size()) ? (defined[residue] != 0) : true;
    }
};

struct QtDihedralAutocorrelation {
    std::size_t n_residues = 0;
    std::size_t n_lags = 0;
    QtPerResidueCurve  phi_acf;        // (R, L)
    QtPerResidueCurve  psi_acf;        // (R, L)
    QtPerResidueScalar phi_corr_time;  // (R,) ps
    QtPerResidueScalar psi_corr_time;  // (R,) ps
    // Chi[0..3] composite payload (Option B per 2026-05-29 planning).
    // chi_acf flat layout: (residue, channel ∈ [0..3], sample). One
    // PerClassBlock + 4-channel descriptor covers all four chi torsions
    // rather than landing eight separate descriptors. The lag axis is
    // shared with phi_acf / psi_acf — chi_acf_axis is a copy for
    // accessor symmetry.
    std::vector<double>  chi_acf;       // (R * 4 * L,) row-major
    std::vector<double>  chi_acf_axis;  // (L,) ps; same content as phi_acf.axis_values
    std::vector<double>  chi_corr_time; // (R * 4,) ps
    std::vector<uint8_t> chi_defined;   // (R * 4,) 1=defined, 0=N/A
    double sample_interval_ps = 0.0;
    QString result_name;

    double chiAcfAt(std::size_t residue, std::size_t chi, std::size_t sample) const {
        if (residue >= n_residues || chi >= 4 || sample >= n_lags || chi_acf.empty())
            return 0.0;
        return chi_acf[(residue * 4 + chi) * n_lags + sample];
    }
    double chiCorrTimeAt(std::size_t residue, std::size_t chi) const {
        if (residue >= n_residues || chi >= 4)
            return 0.0;
        const std::size_t idx = residue * 4 + chi;
        return idx < chi_corr_time.size() ? chi_corr_time[idx] : 0.0;
    }
    bool chiIsDefined(std::size_t residue, std::size_t chi) const {
        if (residue >= n_residues || chi >= 4)
            return false;
        const std::size_t idx = residue * 4 + chi;
        return idx < chi_defined.size() ? (chi_defined[idx] != 0) : false;
    }
};


}  // namespace h5reader::model

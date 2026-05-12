#pragma once
//
// LarsenHBondGrid: HDF5-backed loader for Larsen 2015 ProCS15 H-bond
// DFT grids. Provides a unified
//   `QueryNearest(donor_class, acceptor_class, LarsenHBondGeometry)`
// lookup over 6 (donor × acceptor) grids:
//
//   AmideHydrogen × {BackboneCarbonyl, Hydroxyl, Carboxylate}   (Δσ_HB)
//   AlphaHydrogen × {BackboneCarbonyl, Hydroxyl, Carboxylate}   (Δσ_HαB)
//
// Storage: 6 cubic-spline-pre-computed dense regular grids stored as
// HDF5 files under `data/larsen_hbond_grids/<STEM>_dense.h5`, produced
// by `scripts/larsen_hbond_grid_parse/pre_compute_dense_grids.py` and
// `convert_dense_to_h5.py`. Per-grid axes (r, θ, ρ) are regular; spline
// pre-compute already absorbed the small scattered-data noise from the
// raw DFT scans.
//
// Runtime lookup is trilinear (the cubic curvature is baked into the
// dense grid). ρ wraps periodically at ±180°. Out-of-range (r outside
// the r-axis bounds, θ outside [90°, 180°]) returns is_hit=false; the
// calculator side decides whether to clamp, fall back, or skip.
// A small FP tolerance (1e-9) absorbs computed-value noise at the
// axis bounds (e.g. θ = 180.0 + 1e-12 from cos/acos round-trip).
//
// ─────────────────────────────────────────────────────────────────────
// Geometry contract — what θ and ρ mean
// ─────────────────────────────────────────────────────────────────────
//
// The four atoms used to define the geometry are:
//
//   donor_H       — the donor hydrogen (Hα for ALA donor; HN for NMA
//                   donor).
//   acceptor_O    — the H-bond acceptor oxygen.
//   acceptor_C    — the heavy atom bonded to acceptor_O. For carbonyl
//                   acceptors this is the C=O carbon; for hydroxyl
//                   acceptors this is the sp3 C bonded to OH.
//   acceptor_third — the second neighbour used to disambiguate the
//                   dihedral. Per-acceptor-class mapping below.
//
// Geometry definitions:
//
//   r       = |donor_H − acceptor_O|                            (Å)
//   theta   = angle(donor_H — acceptor_O — acceptor_C),
//             vertex at acceptor_O. 90° to 180°.                (deg)
//   rho     = dihedral(donor_H, acceptor_O, acceptor_C,
//                       acceptor_third),
//             standard IUPAC convention (atan2 of cross products),
//             −180° to +180°.                                    (deg)
//
// IMPORTANT: Larsen 2015's filename ρ has the opposite sign from the
// IUPAC convention this loader uses. The parser keys the grid on the
// IUPAC-sign ρ (from the actual atomic positions in the Gaussian
// logs), so callers should compute IUPAC-sign ρ from the protein
// atoms — do NOT negate the Larsen filename value.
//
// ─────────────────────────────────────────────────────────────────────
// Per-class atom-role mapping
// ─────────────────────────────────────────────────────────────────────
//
//   Donor class         donor_H     donor_anchor   donor_third
//   ───────────────     ─────────   ────────────   ──────────────────
//   AlphaHydrogen       Hα(i)        Cα(i)          N(i)                (own amide N)
//   AmideHydrogen       HN(i)        N(i)           C'(i−1)             (PRECEDING residue's carbonyl C)
//
//   Acceptor class       acceptor_O   acceptor_C    acceptor_third
//   ─────────────────    ──────────   ───────────   ──────────────
//   BackboneCarbonyl     O(j)         C'(j)         N(j+1)
//   SidechainCarbonyl    Asn OD1 /    Asn CG /      Asn ND2 /
//                        Gln OE1      Gln CD        Gln NE2
//   HydroxylOxygen       SER OG /     SER CB /      SER HG /
//                        THR OG1 /    THR CB /      THR HG1 /
//                        TYR OH       TYR CZ        TYR HH
//   CarboxylateOxygen    closer of    ASP CG /      OTHER O on the
//                        Asp OD1/2,   GLU CD,       carboxylate C
//                        Glu OE1/2,   C-term C      (the symmetric
//                        or C-term O                 carboxylate O)
//
// The canonical donor frame (used to store / return tensors) is
// defined entirely by the three donor atoms, NOT by the protein
// secondary structure or by acceptor atoms. Tensor rotation at runtime
// is the caller's responsibility (see ComputeLarsenDonorFrame).
//
// ─────────────────────────────────────────────────────────────────────
// Canonical donor frame
// ─────────────────────────────────────────────────────────────────────
//
//   origin   = donor_H
//   z-axis   = normalize(donor_H − donor_anchor)        (anchor → H)
//   x-axis   = component of (donor_H − donor_third) orthogonal to z,
//              normalized   (i.e. in the (donor_anchor, donor_third,
//                            donor_H) plane, pointing toward donor_third
//                            relative to z)
//   y-axis   = cross(z, x)
//
// The free function `ComputeLarsenDonorFrame(h_pos, anchor_pos,
// third_pos)` builds this rotation matrix; the parser and the future
// calculator both use it to keep the tensor basis consistent.
//
// Per-archive readout atoms (in the canonical frame above):
//   ALA donor: donor {N, CA, CB, C, HA, HN}
//   NMA donor: donor {N, CA, C, HA, HN}      (no CB — NMA has no sidechain)
//   NMA acceptor (Backbone): acceptor {N, HN, HA}  (Δσ_2° contributions,
//                                                  representing residue i+1)
//   HOMe/Acetate acceptor: no acceptor readouts (Larsen 2015 does not
//                                                define 2° terms here).
//
// Per-atom-type contribution dispatch (Larsen 2015 Table 2) is done by
// the calculator (`LarsenHBondShieldingResult`, future session), not
// here — this layer just returns whichever readouts the archive has.
//
// SidechainCarbonyl acceptor (ASN ODE1, GLN OE1) is approximated by the
// BackboneCarbonyl grid; documented limitation, no separate grid in
// Larsen's pipeline.
//
// Reference subtraction: parser-side. Grids store Δσ (free-monomer σ
// subtracted via r-max-edge proxy; see parser README for the
// approximation). No further subtraction at runtime.
//
// Validity mask: each archive carries a uint8 mask (Nr × Nθ × Nρ)
// marking which dense-grid cells correspond to nominal bins where
// real DFT data existed (1) vs nominal bins that were filled by
// nearest-neighbour imputation at parse time (0). Larsen 2015 noted
// some close-approach DFT calcs failed; the pre-compute fills them
// for spline continuity. The dense mask is a nearest-nominal lookup
// (approximate — a dense cell whose nearest nominal is valid but
// whose cubic-spline stencil includes imputed nominals is still
// marked valid). Per-query, the `any_corner_imputed` field on
// `LarsenHBondRecord` is true iff any of the 8 trilinear corner
// cells is an imputed bin under this nearest-nominal definition.
//
//
// Tensor rotation contract (CRITICAL for the calculator):
//
//   Grid tensors are stored such that
//     σ_canonical_stored = R_log · σ_log · R_logᵀ
//   where R_log = ComputeLarsenDonorFrame applied to the log-frame
//   donor atoms. At runtime, the calculator computes
//     R_protein = ComputeLarsenDonorFrame(H_p, anchor_p, third_p)
//   from protein-side atom positions and recovers the protein-lab-
//   frame tensor via
//     σ_lab = R_proteinᵀ · σ_canonical · R_protein.
//   Use `RotateTensorToProteinLabFrame(σ_canonical, R_protein)` so
//   this multiplication direction is hidden at the call site.
//

#include "Types.h"

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace nmr {


// Donor class — which side of the H-bond is the donor.
// Encodes the typed substrate role of the donor atom in the protein.
enum class HBondDonorClass : std::uint8_t {
    AmideHydrogen   = 0,  // Backbone H (N-H). Maps to NMA donor grids.
    AlphaHydrogen   = 1,  // Backbone Hα (Cα-H). Maps to ALA donor grids.
};


// Acceptor class — what the acceptor O's chemical context is.
// Maps to a specific Larsen grid archive (with SidechainCarbonyl
// approximated by BackboneCarbonyl per the design doc).
enum class HBondAcceptorClass : std::uint8_t {
    BackboneCarbonyl    = 0,  // Backbone C=O (any residue's backbone O).
    SidechainCarbonyl   = 1,  // ASN ODE1, GLN OE1. Approximated by NMA grid.
    HydroxylOxygen      = 2,  // SER OG, THR OG1, TYR OH.
    CarboxylateOxygen   = 3,  // ASP OD1/OD2, GLU OE1/OE2, C-terminal O.
};


// Typed H-bond geometry. Encapsulates the three scan parameters
// `r`, `theta`, `rho` with explicit units in the field names so a
// caller can't accidentally swap degrees↔radians at the API boundary.
// All values are computed from the four atom positions per the
// contract documented in the file header.
struct LarsenHBondGeometry {
    double r_angstrom = 0.0;   // |donor_H − acceptor_O|, Å.
    double theta_deg  = 0.0;   // angle(donor_H, acceptor_O, acceptor_C), [90, 180].
    double rho_deg    = 0.0;   // dihedral(donor_H, acceptor_O, acceptor_C,
                               // acceptor_third), IUPAC sign, [−180, 180].
};


// Compute (r, θ, ρ) from the four atom positions per the geometry
// contract. Free function so both the parser (Python side, conceptually)
// and the calculator (C++ runtime, future) compute geometry the same
// way. ρ uses standard IUPAC dihedral convention (atan2 of cross
// products), matching what the parser keys grid bins on.
LarsenHBondGeometry ComputeLarsenHBondGeometry(
    const Vec3& donor_H_pos,
    const Vec3& acceptor_O_pos,
    const Vec3& acceptor_C_pos,
    const Vec3& acceptor_third_pos);


// Compute the canonical-donor-frame rotation matrix R such that
//     v_canonical = R * (v_log − donor_H_pos)
// (the canonical frame is centered at donor_H; this function returns
// only the rotation, not the translation). See the file header for
// the canonical-frame definition.
//
// Caller plugs in the three protein atom positions according to the
// per-class mapping table in the file header (Hα/Cα/N for ALA donor;
// HN/N/C'(i−1) for NMA donor — note the PRECEDING residue's carbonyl C
// for the amide-H case, NOT this residue's CA).
//
// Returns Mat3::Identity if any of the three reference vectors is
// degenerate (donor_H coincident with donor_anchor or donor_third, or
// donor_third on the donor_anchor → donor_H line). Logs a warning via
// OperationLog::Warn. Real protein geometry should never trigger this
// — it's a guard against misconfigured callers.
Mat3 ComputeLarsenDonorFrame(
    const Vec3& donor_H_pos,
    const Vec3& donor_anchor_pos,
    const Vec3& donor_third_pos);


// Transform a tensor returned by QueryNearest from the CANONICAL DONOR
// FRAME (in which grid tensors are stored) into the protein lab frame.
//
// `R_protein` is the rotation matrix returned by ComputeLarsenDonorFrame
// applied to the protein-side atom positions. The grid tensors satisfy
//   σ_canonical = R_log · σ_log · R_logᵀ
// at parse time, so the corresponding lab-frame tensor is
//   σ_lab = R_proteinᵀ · σ_canonical · R_protein.
// Use this helper to avoid hand-rolling the multiplication and getting
// the direction wrong.
Mat3 RotateTensorToProteinLabFrame(
    const Mat3& sigma_canonical,
    const Mat3& R_protein);


// Returned tensors per readout. Each Mat3 is the rotation-invariant
// shielding contribution in the CANONICAL DONOR FRAME (consumer
// rotates to protein lab frame at calculator runtime).
//
// The `has_*` flags indicate which readouts the archive populates.
// Currently:
//   - donor_N, donor_CA, donor_C, donor_HA, donor_HN: always populated
//     (ALA + NMA donor archives).
//   - donor_CB: populated only for ALA donor (NMA has no Cβ).
//   - acceptor_*: populated only for BackboneCarbonyl acceptor (NMA
//     acceptor archives, where the i+1 mapping is defined).
struct LarsenHBondRecord {
    Mat3 donor_N  = Mat3::Zero();
    Mat3 donor_CA = Mat3::Zero();
    Mat3 donor_CB = Mat3::Zero();
    Mat3 donor_C  = Mat3::Zero();
    Mat3 donor_HA = Mat3::Zero();
    Mat3 donor_HN = Mat3::Zero();

    Mat3 acceptor_N  = Mat3::Zero();
    Mat3 acceptor_HN = Mat3::Zero();
    Mat3 acceptor_HA = Mat3::Zero();

    bool has_donor_CB         = false;
    bool has_acceptor_readouts = false;

    // Snapped / interpolated query coordinates, after FP tolerance
    // clamp on r/θ and periodic wrap on ρ to [−180°, +180°). Use these
    // (not the raw input) when logging diagnostics so the canonical
    // grid coordinates are what's recorded.
    double r_angstrom = 0.0;
    double theta_deg  = 0.0;
    double rho_deg    = 0.0;

    // True iff any of the 8 trilinear corner cells was an imputed bin
    // (nearest-neighbour fill of a Larsen-failed DFT grid point).
    // The calculator may want to log + downweight or skip these.
    bool any_corner_imputed = false;

    bool is_hit = false;
    bool IsHit() const { return is_hit; }
};


// One archive's pre-computed dense grid. The axes are regular; the
// tensor arrays have shape (Nr, Nθ, Nρ, 3, 3) stored row-major.
// Internal storage detail; not exposed publicly.
struct LarsenHBondDenseGrid {
    std::vector<double> r_axis;
    std::vector<double> theta_axis;
    std::vector<double> rho_axis;
    int Nr = 0, Ntheta = 0, Nrho = 0;

    // Donor readouts (always present). Each std::vector has
    // Nr * Ntheta * Nrho * 9 doubles, row-major (i_r, i_θ, i_ρ, row, col).
    std::vector<float> donor_N;
    std::vector<float> donor_CA;
    std::vector<float> donor_CB;   // empty for NMA donor archives
    std::vector<float> donor_C;
    std::vector<float> donor_HA;
    std::vector<float> donor_HN;

    // Acceptor readouts (NMA acceptor archives only; empty otherwise).
    std::vector<float> acceptor_N;
    std::vector<float> acceptor_HN;
    std::vector<float> acceptor_HA;

    // Validity mask: Nr × Nθ × Nρ uint8, 1 = real DFT, 0 = imputed by
    // nearest-neighbour fill at parse time. Empty if the archive did
    // not emit a mask (older grids).
    std::vector<std::uint8_t> validity_mask;

    bool has_donor_CB         = false;
    bool has_acceptor_readouts = false;
    bool has_validity_mask    = false;
};


class LarsenHBondGrid {
public:
    // Ctor loads all 6 dense.h5 files from `data_dir`. Throws
    // std::runtime_error if the directory is missing or any expected
    // file is malformed.
    explicit LarsenHBondGrid(const std::string& data_dir);
    ~LarsenHBondGrid();

    LarsenHBondGrid(const LarsenHBondGrid&) = delete;
    LarsenHBondGrid& operator=(const LarsenHBondGrid&) = delete;

    // Query the grid at the given (donor_class, acceptor_class, geom)
    // tuple. Returns a record with is_hit=true if (r, θ) are within the
    // grid bounds (with ±1e-9 FP tolerance at the bounds). is_hit=false
    // if outside; tensors are then zeroed (the caller decides whether
    // to clamp, fall back, or skip).
    //
    // ρ is wrapped periodically to [-180, 180) before lookup. The
    // record's rho_deg field carries the wrapped value (NOT the raw
    // input) for diagnostic clarity.
    //
    // Trilinear interpolation across the 8 corner cells. The cubic
    // smoothness is baked into the dense grid via the Python pre-
    // compute step; trilinear is sufficient at the dense-grid
    // resolution (5×/2×/3× denser than the original DFT scan).
    //
    // The record's `any_corner_imputed` flag is set if any of the 8
    // corner cells came from nearest-neighbour fill at parse time.
    LarsenHBondRecord QueryNearest(
        HBondDonorClass    donor_class,
        HBondAcceptorClass acceptor_class,
        const LarsenHBondGeometry& geom) const;

    // Health probe.
    bool IsLoaded() const { return loaded_; }

    // The data directory passed to the ctor (for logging).
    const std::string& DataDir() const { return data_dir_; }

private:
    std::string data_dir_;
    bool loaded_ = false;

    // 6 grids in canonical (donor, acceptor) order. Indexed by
    // ArchiveIndex(donor_class, acceptor_class).
    std::array<LarsenHBondDenseGrid, 6> grids_;

    // Map (donor, acceptor) → grids_ index. SidechainCarbonyl folds
    // into BackboneCarbonyl per Larsen's grid set (no separate grid).
    static int ArchiveIndex(HBondDonorClass donor,
                            HBondAcceptorClass acceptor);

    // File name for a (donor, acceptor) pair.
    static const char* ArchiveStem(int idx);
};


}  // namespace nmr

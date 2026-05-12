#pragma once
//
// LarsenHBondGrid: HDF5-backed loader for Larsen 2015 ProCS15 H-bond
// DFT grids. Provides a unified `QueryNearest(donor_class, acceptor_class,
// r, theta, rho)` lookup over 6 (donor × acceptor) grids:
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
//
// Tensors are stored and returned in the CANONICAL DONOR FRAME defined
// by the parser:
//   - origin at donor H (Hα for ALA donor, amide H for NMA donor);
//   - z-axis along donor_anchor → donor_H (Cα→Hα or N→H);
//   - x-axis in the donor_anchor / donor_third / donor_H plane;
//   - y-axis = z × x.
// The consumer at calculator runtime must compute the same canonical
// frame from the protein's atom positions and rotate the returned
// tensor into the protein lab frame.
//
// Per-archive readout atoms:
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

    // Snapped / interpolated query coordinates (for diagnostic logging).
    double r     = 0.0;
    double theta = 0.0;
    double rho   = 0.0;

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

    bool has_donor_CB         = false;
    bool has_acceptor_readouts = false;
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

    // Query the grid at (r, theta, rho) for the given donor + acceptor
    // class. r in Å, theta in degrees, rho in degrees.
    //
    // Returns a record with is_hit=true if (r, theta) are within the
    // grid bounds. is_hit=false if outside; tensors are then zeroed
    // (the caller decides whether to clamp, fall back, or skip).
    //
    // rho is wrapped periodically to [-180, 180) before lookup.
    //
    // Trilinear interpolation across the 8 corner cells. The cubic
    // smoothness is baked into the dense grid via the Python pre-
    // compute step; trilinear is sufficient at the dense-grid
    // resolution (5×/2×/3× denser than the original DFT scan).
    LarsenHBondRecord QueryNearest(
        HBondDonorClass    donor_class,
        HBondAcceptorClass acceptor_class,
        double r, double theta, double rho) const;

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

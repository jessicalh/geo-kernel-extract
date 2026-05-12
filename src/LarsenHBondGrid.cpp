#include "LarsenHBondGrid.h"

#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace nmr {

namespace {

// Archive order: 0=NMANMA, 1=NMACOH, 2=NMACOO, 3=ALANMA, 4=ALACOH, 5=ALACOO.
constexpr int kArchiveNMANMA = 0;
constexpr int kArchiveNMACOH = 1;
constexpr int kArchiveNMACOO = 2;
constexpr int kArchiveALANMA = 3;
constexpr int kArchiveALACOH = 4;
constexpr int kArchiveALACOO = 5;

const char* const kArchiveStems[6] = {
    "NMANMA", "NMACOH", "NMACOO", "ALANMA", "ALACOH", "ALACOO",
};

// Axis-bound tolerance for FP round-trip noise (1e-9 of an axis step).
constexpr double kAxisBoundTolerance = 1e-9;


// Read a 5D float32 (Nr × Ntheta × Nrho × 3 × 3) dataset into a flat
// std::vector<float> (row-major). Returns empty vector if the dataset
// doesn't exist (some archives don't have CB or acceptor readouts).
// Validates ALL five dimensions against expected_Nr / Ntheta / Nrho;
// throws on mismatch (catches regenerated-archive shape drift).
std::vector<float> ReadFlatTensorOptional(HighFive::File& f,
                                          const std::string& name,
                                          int expected_Nr,
                                          int expected_Ntheta,
                                          int expected_Nrho) {
    if (!f.exist(name)) return {};
    auto ds = f.getDataSet(name);
    auto dims = ds.getDimensions();
    if (dims.size() != 5 ||
        static_cast<int>(dims[0]) != expected_Nr ||
        static_cast<int>(dims[1]) != expected_Ntheta ||
        static_cast<int>(dims[2]) != expected_Nrho ||
        dims[3] != 3 || dims[4] != 3) {
        throw std::runtime_error(
            "LarsenHBondGrid: dataset " + name +
            " has unexpected shape (expected (" +
            std::to_string(expected_Nr) + ", " +
            std::to_string(expected_Ntheta) + ", " +
            std::to_string(expected_Nrho) +
            ", 3, 3))");
    }
    std::size_t total = dims[0] * dims[1] * dims[2] * 9;
    std::vector<float> flat(total);
    ds.read(flat.data());
    // Reject NaN/Inf in stored tensors. Parser/pre-compute should never
    // emit these; the assertion catches drift in the upstream pipeline.
    for (float v : flat) {
        if (!std::isfinite(v)) {
            throw std::runtime_error(
                "LarsenHBondGrid: non-finite value in dataset " + name);
        }
    }
    return flat;
}


std::vector<std::uint8_t> ReadValidityMaskOptional(HighFive::File& f,
                                                   int Nr, int Ntheta,
                                                   int Nrho) {
    const std::string name = "validity_mask";
    if (!f.exist(name)) return {};
    auto ds = f.getDataSet(name);
    auto dims = ds.getDimensions();
    if (dims.size() != 3 ||
        static_cast<int>(dims[0]) != Nr ||
        static_cast<int>(dims[1]) != Ntheta ||
        static_cast<int>(dims[2]) != Nrho) {
        throw std::runtime_error(
            "LarsenHBondGrid: validity_mask shape mismatch");
    }
    std::vector<std::uint8_t> flat(
        static_cast<std::size_t>(Nr) * Ntheta * Nrho);
    ds.read(flat.data());
    return flat;
}


std::vector<double> ReadAxis(HighFive::File& f, const std::string& name) {
    auto ds = f.getDataSet(name);
    auto dims = ds.getDimensions();
    if (dims.size() != 1) {
        throw std::runtime_error(
            "LarsenHBondGrid: axis " + name + " has wrong rank " +
            std::to_string(dims.size()));
    }
    std::vector<double> out(dims[0]);
    ds.read(out.data());
    return out;
}


// Schema validation: each archive identity (ALA vs NMA donor, NMA vs
// HOMe/COO acceptor) has a specific MANDATORY set of readout
// datasets. Optional reads previously allowed silent zeroing if a
// mandatory dataset was missing (codex M3 finding). After loading,
// check the archive's expected set against what was read and throw on
// any mismatch.
//
// Mandatory rules per archive:
//   All archives: donor_N, donor_CA, donor_C, donor_HA, donor_HN.
//   ALA donor archives (ALANMA/ALACOH/ALACOO): donor_CB present.
//   NMA donor archives (NMANMA/NMACOH/NMACOO): donor_CB ABSENT
//     (NMA has no Cβ — presence would be a parser regression).
//   NMA acceptor archives (NMANMA, ALANMA): acceptor_N, _HN, _HA all
//     present.
//   HOMe / COO acceptor archives: acceptor_* all ABSENT (Larsen 2015
//     does not define 2° terms for those acceptor classes).
void ValidateSchema(const fs::path& h5_path,
                    const LarsenHBondDenseGrid& g,
                    const std::string& archive_stem) {
    auto fail = [&](const std::string& msg) {
        throw std::runtime_error(
            "LarsenHBondGrid: schema validation failed for " +
            h5_path.string() + " (" + archive_stem + "): " + msg);
    };
    if (g.donor_N.empty()) fail("missing mandatory donor_N");
    if (g.donor_CA.empty()) fail("missing mandatory donor_CA");
    if (g.donor_C.empty())  fail("missing mandatory donor_C");
    if (g.donor_HA.empty()) fail("missing mandatory donor_HA");
    if (g.donor_HN.empty()) fail("missing mandatory donor_HN");

    const bool is_ala_donor = archive_stem.substr(0, 3) == "ALA";
    const bool is_nma_donor = archive_stem.substr(0, 3) == "NMA";
    if (!is_ala_donor && !is_nma_donor) {
        fail("archive stem does not start with ALA or NMA");
    }
    if (is_ala_donor && g.donor_CB.empty()) {
        fail("ALA donor archive missing mandatory donor_CB");
    }
    if (is_nma_donor && !g.donor_CB.empty()) {
        fail("NMA donor archive has unexpected donor_CB (NMA has no Cβ)");
    }

    const bool is_nma_acceptor = (archive_stem == "NMANMA"
                                  || archive_stem == "ALANMA");
    if (is_nma_acceptor) {
        if (g.acceptor_N.empty())  fail("NMA acceptor missing acceptor_N");
        if (g.acceptor_HN.empty()) fail("NMA acceptor missing acceptor_HN");
        if (g.acceptor_HA.empty()) fail("NMA acceptor missing acceptor_HA");
    } else {
        if (!g.acceptor_N.empty() ||
            !g.acceptor_HN.empty() ||
            !g.acceptor_HA.empty()) {
            fail("non-NMA acceptor archive has unexpected acceptor_* data "
                 "(only NMA acceptor grids have a defined i+1 mapping)");
        }
    }
}


void LoadOne(const fs::path& h5_path, LarsenHBondDenseGrid& g,
             const std::string& archive_stem) {
    if (!fs::exists(h5_path)) {
        throw std::runtime_error(
            "LarsenHBondGrid: file not found: " + h5_path.string());
    }
    HighFive::File f(h5_path.string(), HighFive::File::ReadOnly);

    g.r_axis = ReadAxis(f, "r_axis");
    g.theta_axis = ReadAxis(f, "theta_axis");
    g.rho_axis = ReadAxis(f, "rho_axis");
    g.Nr = static_cast<int>(g.r_axis.size());
    g.Ntheta = static_cast<int>(g.theta_axis.size());
    g.Nrho = static_cast<int>(g.rho_axis.size());

    g.donor_N  = ReadFlatTensorOptional(f, "donor_N",  g.Nr, g.Ntheta, g.Nrho);
    g.donor_CA = ReadFlatTensorOptional(f, "donor_CA", g.Nr, g.Ntheta, g.Nrho);
    g.donor_CB = ReadFlatTensorOptional(f, "donor_CB", g.Nr, g.Ntheta, g.Nrho);
    g.donor_C  = ReadFlatTensorOptional(f, "donor_C",  g.Nr, g.Ntheta, g.Nrho);
    g.donor_HA = ReadFlatTensorOptional(f, "donor_HA", g.Nr, g.Ntheta, g.Nrho);
    g.donor_HN = ReadFlatTensorOptional(f, "donor_HN", g.Nr, g.Ntheta, g.Nrho);

    g.acceptor_N  = ReadFlatTensorOptional(f, "acceptor_N",  g.Nr, g.Ntheta, g.Nrho);
    g.acceptor_HN = ReadFlatTensorOptional(f, "acceptor_HN", g.Nr, g.Ntheta, g.Nrho);
    g.acceptor_HA = ReadFlatTensorOptional(f, "acceptor_HA", g.Nr, g.Ntheta, g.Nrho);

    g.validity_mask = ReadValidityMaskOptional(f, g.Nr, g.Ntheta, g.Nrho);

    g.has_donor_CB         = !g.donor_CB.empty();
    g.has_acceptor_readouts = !g.acceptor_N.empty();
    g.has_validity_mask    = !g.validity_mask.empty();

    ValidateSchema(h5_path, g, archive_stem);
}


// Trilinear interpolation of a flat 5D tensor array at (ir + fr, ith + fth, irho + frho)
// where 0 ≤ f* ≤ 1. Indices are clamped at boundaries (caller ensures
// they're in range; periodic ρ wrap is handled by the caller via the
// next_irho parameter).
Mat3 TrilinearMat3(const std::vector<float>& flat,
                   int Nth, int Nrho,
                   int ir, int ith, int irho,
                   int ir_next, int ith_next, int irho_next,
                   double fr, double fth, double frho) {
    auto idx = [&](int i_r, int i_th, int i_rho) -> std::size_t {
        return static_cast<std::size_t>(i_r) * Nth * Nrho * 9
             + static_cast<std::size_t>(i_th) * Nrho * 9
             + static_cast<std::size_t>(i_rho) * 9;
    };
    Mat3 out = Mat3::Zero();
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            double c000 = flat[idx(ir,      ith,      irho)      + row * 3 + col];
            double c001 = flat[idx(ir,      ith,      irho_next) + row * 3 + col];
            double c010 = flat[idx(ir,      ith_next, irho)      + row * 3 + col];
            double c011 = flat[idx(ir,      ith_next, irho_next) + row * 3 + col];
            double c100 = flat[idx(ir_next, ith,      irho)      + row * 3 + col];
            double c101 = flat[idx(ir_next, ith,      irho_next) + row * 3 + col];
            double c110 = flat[idx(ir_next, ith_next, irho)      + row * 3 + col];
            double c111 = flat[idx(ir_next, ith_next, irho_next) + row * 3 + col];

            double c00 = c000 * (1.0 - frho) + c001 * frho;
            double c01 = c010 * (1.0 - frho) + c011 * frho;
            double c10 = c100 * (1.0 - frho) + c101 * frho;
            double c11 = c110 * (1.0 - frho) + c111 * frho;

            double c0 = c00 * (1.0 - fth) + c01 * fth;
            double c1 = c10 * (1.0 - fth) + c11 * fth;

            out(row, col) = c0 * (1.0 - fr) + c1 * fr;
        }
    }
    return out;
}


// Check if any of the 8 corner cells (ir,ith,irho), (ir+1,...), ..., (ir+1,ith+1,irho+1)
// is an imputed bin (validity_mask == 0). Returns false if the grid has
// no mask.
bool AnyCornerImputed(const LarsenHBondDenseGrid& g,
                      int ir, int ith, int irho,
                      int ir_next, int ith_next, int irho_next) {
    if (!g.has_validity_mask) return false;
    auto idx = [&](int i_r, int i_th, int i_rho) -> std::size_t {
        return static_cast<std::size_t>(i_r) * g.Ntheta * g.Nrho
             + static_cast<std::size_t>(i_th) * g.Nrho
             + static_cast<std::size_t>(i_rho);
    };
    int rs[2] = {ir, ir_next};
    int ths[2] = {ith, ith_next};
    int rhos[2] = {irho, irho_next};
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 2; ++c)
                if (g.validity_mask[idx(rs[a], ths[b], rhos[c])] == 0)
                    return true;
    return false;
}


// For a regular ascending axis, return (idx, frac) such that
// axis[idx] + frac * (axis[idx+1] - axis[idx]) == value, with frac in
// [0, 1) when value is in range. Returns idx=-1 when strictly out-of-range
// (after a small FP tolerance clamp at the bounds).
struct AxisLookup {
    int idx = -1;
    int idx_next = -1;
    double frac = 0.0;
};

AxisLookup LookupAxis(const std::vector<double>& axis, double value,
                      bool periodic = false) {
    AxisLookup out;
    int n = static_cast<int>(axis.size());
    if (n < 2) return out;

    if (periodic) {
        // Wrap value to [axis[0], axis[0] + period).
        double step = axis[1] - axis[0];
        double period = axis.back() - axis[0] + step;  // full periodic period
        // Sanity check: for the H-bond ρ axis the period must be 360°.
        // If a future grid breaks this invariant the wrap is wrong; bail.
        if (std::abs(period - 360.0) > 1e-6) {
            return out;  // idx remains -1 — caller bails on miss
        }
        double v = value;
        // Bring v into the half-open interval starting at axis[0].
        v = axis[0] + std::fmod(v - axis[0], period);
        if (v < axis[0]) v += period;

        double f = (v - axis[0]) / step;
        int i = static_cast<int>(std::floor(f));
        i = std::min(i, n - 1);
        i = std::max(i, 0);
        int i_next = i + 1;
        double frac = f - i;
        if (i_next >= n) {
            i_next = 0;        // wrap
        }
        out.idx = i;
        out.idx_next = i_next;
        out.frac = frac;
        return out;
    }

    double step = axis[1] - axis[0];
    double tol = std::abs(step) * kAxisBoundTolerance;
    if (value < axis.front() - tol || value > axis.back() + tol) {
        return out;  // out-of-range
    }
    // Clamp tiny FP overshoot to axis bounds.
    double v = std::clamp(value, axis.front(), axis.back());
    double f = (v - axis.front()) / step;
    int i = static_cast<int>(std::floor(f));
    i = std::min(i, n - 2);
    i = std::max(i, 0);
    out.idx = i;
    out.idx_next = i + 1;
    out.frac = f - i;
    if (out.frac < 0) out.frac = 0;
    if (out.frac > 1) out.frac = 1;
    return out;
}


// Wrap ρ to canonical [-180, 180) range, matching what the loader
// queries the axis on.
double WrapRho(double rho_deg) {
    double w = std::fmod(rho_deg + 180.0, 360.0);
    if (w < 0.0) w += 360.0;
    return w - 180.0;
}


}  // namespace


// -----------------------------------------------------------------------------
// Free functions: ComputeLarsenHBondGeometry, ComputeLarsenDonorFrame
// -----------------------------------------------------------------------------

LarsenHBondGeometry ComputeLarsenHBondGeometry(
    const Vec3& donor_H_pos,
    const Vec3& acceptor_O_pos,
    const Vec3& acceptor_C_pos,
    const Vec3& acceptor_third_pos) {

    LarsenHBondGeometry geom;
    Vec3 H_to_O = acceptor_O_pos - donor_H_pos;
    geom.r_angstrom = H_to_O.norm();

    // theta at acceptor O: angle between (O→H) and (O→C).
    Vec3 O_to_H = donor_H_pos - acceptor_O_pos;
    Vec3 O_to_C = acceptor_C_pos - acceptor_O_pos;
    double cos_theta = O_to_H.dot(O_to_C) / (O_to_H.norm() * O_to_C.norm());
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);
    geom.theta_deg = std::acos(cos_theta) * (180.0 / M_PI);

    // rho dihedral: H - O - C - third. Standard IUPAC convention via
    // atan2(cross(n1, b2_norm) · n2, n1 · n2).
    Vec3 b1 = acceptor_O_pos     - donor_H_pos;
    Vec3 b2 = acceptor_C_pos     - acceptor_O_pos;
    Vec3 b3 = acceptor_third_pos - acceptor_C_pos;
    Vec3 n1 = b1.cross(b2);
    Vec3 n2 = b2.cross(b3);
    Vec3 m1 = n1.cross(b2.normalized());
    double x = n1.dot(n2);
    double y = m1.dot(n2);
    geom.rho_deg = std::atan2(y, x) * (180.0 / M_PI);

    return geom;
}


Mat3 ComputeLarsenDonorFrame(
    const Vec3& donor_H_pos,
    const Vec3& donor_anchor_pos,
    const Vec3& donor_third_pos) {

    constexpr double kTinyVec = 1e-9;

    // z = normalize(donor_H − donor_anchor). Bail on coincident atoms.
    Vec3 z_raw = donor_H_pos - donor_anchor_pos;
    if (z_raw.norm() < kTinyVec) {
        OperationLog::Warn("ComputeLarsenDonorFrame",
            "donor_H and donor_anchor coincide; returning identity rotation");
        return Mat3::Identity();
    }
    Vec3 z = z_raw.normalized();

    // x = component of (donor_H − donor_third) orthogonal to z, normalized.
    // Bail if third is coincident with H or lies on the anchor→H line
    // (the orthogonal component is then zero).
    Vec3 v3 = donor_H_pos - donor_third_pos;
    if (v3.norm() < kTinyVec) {
        OperationLog::Warn("ComputeLarsenDonorFrame",
            "donor_H and donor_third coincide; returning identity rotation");
        return Mat3::Identity();
    }
    Vec3 x_raw = v3 - (v3.dot(z)) * z;
    if (x_raw.norm() < kTinyVec) {
        OperationLog::Warn("ComputeLarsenDonorFrame",
            "donor_third on the anchor→H line; returning identity rotation");
        return Mat3::Identity();
    }
    Vec3 x = x_raw.normalized();
    Vec3 y = z.cross(x);

    // Rotation matrix: rows are canonical basis vectors expressed in log frame.
    // R @ v_log = v_canonical, so R has rows [x; y; z].
    Mat3 R;
    R.row(0) = x;
    R.row(1) = y;
    R.row(2) = z;
    return R;
}


Mat3 RotateTensorToProteinLabFrame(
    const Mat3& sigma_canonical,
    const Mat3& R_protein) {
    return R_protein.transpose() * sigma_canonical * R_protein;
}


// -----------------------------------------------------------------------------
// LarsenHBondGrid implementation
// -----------------------------------------------------------------------------

int LarsenHBondGrid::ArchiveIndex(HBondDonorClass donor,
                                  HBondAcceptorClass acceptor) {
    // SidechainCarbonyl approximated by BackboneCarbonyl (NMA acceptor
    // grid). Documented limitation per the design doc.
    bool nma_acceptor =
        acceptor == HBondAcceptorClass::BackboneCarbonyl ||
        acceptor == HBondAcceptorClass::SidechainCarbonyl;

    if (donor == HBondDonorClass::AmideHydrogen) {
        if (nma_acceptor)                                        return kArchiveNMANMA;
        if (acceptor == HBondAcceptorClass::HydroxylOxygen)      return kArchiveNMACOH;
        if (acceptor == HBondAcceptorClass::CarboxylateOxygen)   return kArchiveNMACOO;
    } else {  // AlphaHydrogen
        if (nma_acceptor)                                        return kArchiveALANMA;
        if (acceptor == HBondAcceptorClass::HydroxylOxygen)      return kArchiveALACOH;
        if (acceptor == HBondAcceptorClass::CarboxylateOxygen)   return kArchiveALACOO;
    }
    return -1;
}


const char* LarsenHBondGrid::ArchiveStem(int idx) {
    if (idx < 0 || idx >= 6) return "<invalid>";
    return kArchiveStems[idx];
}


LarsenHBondGrid::LarsenHBondGrid(const std::string& data_dir)
    : data_dir_(data_dir) {
    fs::path dir(data_dir);
    if (!fs::exists(dir) || !fs::is_directory(dir)) {
        throw std::runtime_error(
            "LarsenHBondGrid: data_dir not found: " + data_dir);
    }
    for (int i = 0; i < 6; ++i) {
        fs::path h5 = dir / (std::string(kArchiveStems[i]) + "_dense.h5");
        LoadOne(h5, grids_[i], kArchiveStems[i]);
    }
    loaded_ = true;
    OperationLog::Info(LogCalcOther, "LarsenHBondGrid",
        "loaded 6 dense H-bond grids from " + data_dir);
}


LarsenHBondGrid::~LarsenHBondGrid() = default;


LarsenHBondRecord LarsenHBondGrid::QueryNearest(
    HBondDonorClass donor_class,
    HBondAcceptorClass acceptor_class,
    const LarsenHBondGeometry& geom) const {

    LarsenHBondRecord rec;
    rec.r_angstrom = geom.r_angstrom;
    rec.theta_deg  = geom.theta_deg;
    rec.rho_deg    = WrapRho(geom.rho_deg);  // canonical-form ρ for diagnostics

    // NaN/Inf guard: malformed hydrogens or coincident atoms upstream
    // can send NaN through ComputeLarsenHBondGeometry's normalizations.
    // LookupAxis's std::floor / std::clamp paths are undefined on NaN;
    // return a clean non-hit here instead.
    if (!std::isfinite(geom.r_angstrom) ||
        !std::isfinite(geom.theta_deg)  ||
        !std::isfinite(geom.rho_deg)) {
        OperationLog::Warn("LarsenHBondGrid::QueryNearest",
            "non-finite query geometry (r="
            + std::to_string(geom.r_angstrom) + " theta="
            + std::to_string(geom.theta_deg)  + " rho="
            + std::to_string(geom.rho_deg) + ")");
        return rec;
    }

    int idx = ArchiveIndex(donor_class, acceptor_class);
    if (idx < 0) {
        OperationLog::Warn("LarsenHBondGrid::QueryNearest",
            "no archive maps to this (donor, acceptor) pair");
        return rec;
    }
    const LarsenHBondDenseGrid& g = grids_[idx];

    AxisLookup lr = LookupAxis(g.r_axis, geom.r_angstrom, /*periodic=*/false);
    if (lr.idx < 0) return rec;  // r out of range
    AxisLookup lth = LookupAxis(g.theta_axis, geom.theta_deg, /*periodic=*/false);
    if (lth.idx < 0) return rec;  // θ out of range
    AxisLookup lrho = LookupAxis(g.rho_axis, rec.rho_deg, /*periodic=*/true);
    if (lrho.idx < 0) return rec;  // ρ axis malformed (e.g. period != 360°)

    auto interp = [&](const std::vector<float>& flat) -> Mat3 {
        if (flat.empty()) return Mat3::Zero();
        return TrilinearMat3(flat, g.Ntheta, g.Nrho,
                             lr.idx, lth.idx, lrho.idx,
                             lr.idx_next, lth.idx_next, lrho.idx_next,
                             lr.frac, lth.frac, lrho.frac);
    };

    rec.donor_N  = interp(g.donor_N);
    rec.donor_CA = interp(g.donor_CA);
    rec.donor_CB = interp(g.donor_CB);
    rec.donor_C  = interp(g.donor_C);
    rec.donor_HA = interp(g.donor_HA);
    rec.donor_HN = interp(g.donor_HN);

    rec.has_donor_CB         = g.has_donor_CB;
    rec.has_acceptor_readouts = g.has_acceptor_readouts;

    if (g.has_acceptor_readouts) {
        rec.acceptor_N  = interp(g.acceptor_N);
        rec.acceptor_HN = interp(g.acceptor_HN);
        rec.acceptor_HA = interp(g.acceptor_HA);
    }

    // r/θ remain the (clamped, in-range) query values; ρ is the
    // canonical wrapped form set above. (Earlier code recomputed
    // r/θ from axis * frac which round-tripped the FP-clamped
    // input — equivalent but confusing; we now leave them as the
    // canonical query coords.)

    rec.any_corner_imputed = AnyCornerImputed(
        g, lr.idx, lth.idx, lrho.idx, lr.idx_next, lth.idx_next, lrho.idx_next);

    rec.is_hit = true;
    return rec;
}


}  // namespace nmr

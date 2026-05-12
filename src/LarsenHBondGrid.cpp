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


// Read a 3D float32 (Nr × Ntheta × Nrho × 3 × 3) dataset into a flat
// std::vector<float> (row-major). Returns empty vector if the dataset
// doesn't exist (some archives don't have CB or acceptor readouts).
std::vector<float> ReadFlatTensorOptional(HighFive::File& f,
                                          const std::string& name) {
    if (!f.exist(name)) return {};
    auto ds = f.getDataSet(name);
    auto dims = ds.getDimensions();
    if (dims.size() != 5 || dims[3] != 3 || dims[4] != 3) {
        throw std::runtime_error(
            "LarsenHBondGrid: dataset " + name +
            " has unexpected shape (expected (Nr, Nθ, Nρ, 3, 3))");
    }
    std::size_t total = dims[0] * dims[1] * dims[2] * 9;
    std::vector<float> flat(total);
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


void LoadOne(const fs::path& h5_path, LarsenHBondDenseGrid& g) {
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

    g.donor_N  = ReadFlatTensorOptional(f, "donor_N");
    g.donor_CA = ReadFlatTensorOptional(f, "donor_CA");
    g.donor_CB = ReadFlatTensorOptional(f, "donor_CB");
    g.donor_C  = ReadFlatTensorOptional(f, "donor_C");
    g.donor_HA = ReadFlatTensorOptional(f, "donor_HA");
    g.donor_HN = ReadFlatTensorOptional(f, "donor_HN");

    g.acceptor_N  = ReadFlatTensorOptional(f, "acceptor_N");
    g.acceptor_HN = ReadFlatTensorOptional(f, "acceptor_HN");
    g.acceptor_HA = ReadFlatTensorOptional(f, "acceptor_HA");

    g.has_donor_CB         = !g.donor_CB.empty();
    g.has_acceptor_readouts = !g.acceptor_N.empty();
}


// Trilinear interpolation of a flat 5D tensor array at (ir + fr, ith + fth, irho + frho)
// where 0 ≤ f* ≤ 1. Indices are clamped at boundaries (caller ensures
// they're in range; periodic ρ wrap is handled by the caller via the
// next_irho parameter).
//
// flat[ir, ith, irho, row, col] is at index
//   ir * (Nth * Nrho * 9) + ith * (Nrho * 9) + irho * 9 + row * 3 + col.
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


// For a regular ascending axis, return (idx, frac) such that
// axis[idx] + frac * (axis[idx+1] - axis[idx]) == value, with frac in
// [0, 1) when value is in range. Returns idx=-1 when out-of-range.
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
        double period = axis.back() - axis[0]
                        + (axis[1] - axis[0]);  // wrap step
        double v = value;
        // Bring v into the half-open interval starting at axis[0].
        v = axis[0] + std::fmod(v - axis[0], period);
        if (v < axis[0]) v += period;

        // Find idx such that axis[idx] ≤ v < axis[idx+1] (or last bin wraps).
        // Uniform-step assumption: step = axis[1] - axis[0].
        double step = axis[1] - axis[0];
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

    if (value < axis.front() || value > axis.back()) {
        return out;  // out-of-range
    }
    // Find bin via uniform-step shortcut.
    double step = axis[1] - axis[0];
    double f = (value - axis[0]) / step;
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


}  // namespace


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
        LoadOne(h5, grids_[i]);
    }
    loaded_ = true;
    OperationLog::Info(LogCalcOther, "LarsenHBondGrid",
        "loaded 6 dense H-bond grids from " + data_dir);
}


LarsenHBondGrid::~LarsenHBondGrid() = default;


LarsenHBondRecord LarsenHBondGrid::QueryNearest(
    HBondDonorClass donor_class,
    HBondAcceptorClass acceptor_class,
    double r, double theta, double rho) const {

    LarsenHBondRecord rec;
    rec.r = r;
    rec.theta = theta;
    rec.rho = rho;

    int idx = ArchiveIndex(donor_class, acceptor_class);
    if (idx < 0) {
        OperationLog::Warn("LarsenHBondGrid::QueryNearest",
            "no archive maps to this (donor, acceptor) pair");
        return rec;
    }
    const LarsenHBondDenseGrid& g = grids_[idx];

    AxisLookup lr = LookupAxis(g.r_axis, r, /*periodic=*/false);
    if (lr.idx < 0) return rec;  // r out of range
    AxisLookup lth = LookupAxis(g.theta_axis, theta, /*periodic=*/false);
    if (lth.idx < 0) return rec;  // θ out of range
    AxisLookup lrho = LookupAxis(g.rho_axis, rho, /*periodic=*/true);

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

    // Update snapped coords to the interpolation cell base for diagnostics.
    rec.r     = g.r_axis[lr.idx]   + lr.frac   * (g.r_axis[lr.idx_next]   - g.r_axis[lr.idx]);
    rec.theta = g.theta_axis[lth.idx] + lth.frac * (g.theta_axis[lth.idx_next] - g.theta_axis[lth.idx]);
    // ρ snapped: leave as input (already wrapped to grid range internally).
    rec.rho = rho;

    rec.is_hit = true;
    return rec;
}


}  // namespace nmr

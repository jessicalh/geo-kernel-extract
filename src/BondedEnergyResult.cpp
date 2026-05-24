#include "BondedEnergyResult.h"
#include "ConformationAtom.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "Types.h"

#include <cmath>
#include <cassert>

namespace nmr {

// ── Energy functions (CHARMM36m functional forms) ───────────────
// All positions in Angstroms, parameters in GROMACS native units
// (nm, kJ/mol). We convert distances A→nm before evaluating.

static constexpr double A_TO_NM = 0.1;

// Harmonic bond: E = ½k(r - r0)²
static double EvalBond(const Vec3& p0, const Vec3& p1,
                       double r0_nm, double k) {
    double r_nm = (p1 - p0).norm() * A_TO_NM;
    double dr = r_nm - r0_nm;
    return 0.5 * k * dr * dr;
}

// Harmonic angle: E = ½k(θ - θ0)²
static double EvalAngle(const Vec3& p0, const Vec3& p1, const Vec3& p2,
                        double theta0, double k) {
    Vec3 v1 = (p0 - p1).normalized();
    Vec3 v2 = (p2 - p1).normalized();
    double cos_theta = std::max(-1.0, std::min(1.0, v1.dot(v2)));
    double theta = std::acos(cos_theta);
    double dtheta = theta - theta0;
    return 0.5 * k * dtheta * dtheta;
}

// Urey-Bradley: E = ½k_ub(r13 - r13_0)²
static double EvalUB(const Vec3& p0, const Vec3& p2,
                     double r13_0_nm, double k_ub) {
    double r13_nm = (p2 - p0).norm() * A_TO_NM;
    double dr = r13_nm - r13_0_nm;
    return 0.5 * k_ub * dr * dr;
}

// Dihedral angle from four positions (radians).
static double DihedralAngle(const Vec3& p0, const Vec3& p1,
                            const Vec3& p2, const Vec3& p3) {
    Vec3 b1 = p1 - p0;
    Vec3 b2 = p2 - p1;
    Vec3 b3 = p3 - p2;
    Vec3 n1 = b1.cross(b2);
    Vec3 n2 = b2.cross(b3);
    double n1n = n1.norm();
    double n2n = n2.norm();
    if (n1n < 1e-10 || n2n < 1e-10) return 0.0;
    n1 /= n1n;
    n2 /= n2n;
    double cos_phi = std::max(-1.0, std::min(1.0, n1.dot(n2)));
    Vec3 m = n1.cross(b2.normalized());
    double sin_phi = m.dot(n2);
    return std::atan2(sin_phi, cos_phi);
}

// Periodic proper dihedral: E = k(1 + cos(n*φ - φ0))
static double EvalProperDih(const Vec3& p0, const Vec3& p1,
                            const Vec3& p2, const Vec3& p3,
                            double phi0, double k, int mult) {
    double phi = DihedralAngle(p0, p1, p2, p3);
    return k * (1.0 + std::cos(mult * phi - phi0));
}

// Harmonic improper dihedral: E = ½k(φ - φ0)²
static double EvalImproperDih(const Vec3& p0, const Vec3& p1,
                              const Vec3& p2, const Vec3& p3,
                              double phi0, double k) {
    double phi = DihedralAngle(p0, p1, p2, p3);
    // Wrap difference to [-π, π]
    double dphi = phi - phi0;
    while (dphi > M_PI)  dphi -= 2.0 * M_PI;
    while (dphi < -M_PI) dphi += 2.0 * M_PI;
    return 0.5 * k * dphi * dphi;
}

// CMAP: 2D bilinear interpolation on (φ, ψ) grid.
// Grid covers [-π, π] × [-π, π] with grid_spacing points.
static double EvalCMAP(const Vec3& p0, const Vec3& p1,
                       const Vec3& p2, const Vec3& p3, const Vec3& p4,
                       const std::vector<double>& grid, int spacing) {
    if (grid.empty() || spacing < 2) return 0.0;

    double phi = DihedralAngle(p0, p1, p2, p3);
    double psi = DihedralAngle(p1, p2, p3, p4);

    // Map angle [-π, π] to grid index [0, spacing-1]
    double dx = 2.0 * M_PI / spacing;
    double fi = (phi + M_PI) / dx;
    double fj = (psi + M_PI) / dx;

    // Bilinear interpolation (sufficient — CMAP grids are smooth)
    int i0 = static_cast<int>(fi) % spacing;
    int j0 = static_cast<int>(fj) % spacing;
    int i1 = (i0 + 1) % spacing;
    int j1 = (j0 + 1) % spacing;
    double wi = fi - std::floor(fi);
    double wj = fj - std::floor(fj);

    double v00 = grid[i0 * spacing + j0];
    double v10 = grid[i1 * spacing + j0];
    double v01 = grid[i0 * spacing + j1];
    double v11 = grid[i1 * spacing + j1];

    return (1-wi)*(1-wj)*v00 + wi*(1-wj)*v10
         + (1-wi)*wj*v01     + wi*wj*v11;
}


// ── Compute ────────────────────────────────────────────────────

std::unique_ptr<BondedEnergyResult> BondedEnergyResult::Compute(
        ProteinConformation& conf,
        const BondedParameters& params) {

    const size_t N = conf.AtomCount();
    auto result = std::make_unique<BondedEnergyResult>();
    result->bond_energy_.resize(N, 0.0);
    result->angle_energy_.resize(N, 0.0);
    result->ub_energy_.resize(N, 0.0);
    result->proper_energy_.resize(N, 0.0);
    result->improper_energy_.resize(N, 0.0);
    result->cmap_energy_.resize(N, 0.0);
    result->total_bonded_.resize(N, 0.0);

    size_t count_bond = 0, count_angle = 0, count_ub = 0;
    size_t count_proper = 0, count_improper = 0, count_cmap = 0;

    for (const auto& ix : params.interactions) {
        // All atom indices must be in range
        bool valid = true;
        for (int k = 0; k < ix.n_atoms; ++k) {
            if (ix.atoms[k] >= N) { valid = false; break; }
        }
        if (!valid) continue;

        double energy = 0.0;
        std::vector<double>* target = nullptr;

        switch (ix.type) {
            case BondedInteraction::Bond: {
                Vec3 p0 = conf.PositionAt(ix.atoms[0]);
                Vec3 p1 = conf.PositionAt(ix.atoms[1]);
                energy = EvalBond(p0, p1, ix.p[0], ix.p[1]);
                target = &result->bond_energy_;
                ++count_bond;
                break;
            }
            case BondedInteraction::Angle: {
                Vec3 p0 = conf.PositionAt(ix.atoms[0]);
                Vec3 p1 = conf.PositionAt(ix.atoms[1]);
                Vec3 p2 = conf.PositionAt(ix.atoms[2]);
                energy = EvalAngle(p0, p1, p2, ix.p[0], ix.p[1]);
                target = &result->angle_energy_;
                ++count_angle;
                break;
            }
            case BondedInteraction::UreyBradley: {
                Vec3 p0 = conf.PositionAt(ix.atoms[0]);
                Vec3 p2 = conf.PositionAt(ix.atoms[2]);
                energy = EvalUB(p0, p2, ix.p[0], ix.p[1]);
                target = &result->ub_energy_;
                ++count_ub;
                break;
            }
            case BondedInteraction::ProperDih: {
                Vec3 p0 = conf.PositionAt(ix.atoms[0]);
                Vec3 p1 = conf.PositionAt(ix.atoms[1]);
                Vec3 p2 = conf.PositionAt(ix.atoms[2]);
                Vec3 p3 = conf.PositionAt(ix.atoms[3]);
                energy = EvalProperDih(p0, p1, p2, p3,
                                       ix.p[0], ix.p[1],
                                       static_cast<int>(ix.p[2]));
                target = &result->proper_energy_;
                ++count_proper;
                break;
            }
            case BondedInteraction::ImproperDih: {
                Vec3 p0 = conf.PositionAt(ix.atoms[0]);
                Vec3 p1 = conf.PositionAt(ix.atoms[1]);
                Vec3 p2 = conf.PositionAt(ix.atoms[2]);
                Vec3 p3 = conf.PositionAt(ix.atoms[3]);
                energy = EvalImproperDih(p0, p1, p2, p3, ix.p[0], ix.p[1]);
                target = &result->improper_energy_;
                ++count_improper;
                break;
            }
            case BondedInteraction::CMAP: {
                int cmap_idx = static_cast<int>(ix.p[0]);
                if (cmap_idx < 0 ||
                    cmap_idx >= static_cast<int>(params.cmap_grids.size()))
                    continue;
                Vec3 p0 = conf.PositionAt(ix.atoms[0]);
                Vec3 p1 = conf.PositionAt(ix.atoms[1]);
                Vec3 p2 = conf.PositionAt(ix.atoms[2]);
                Vec3 p3 = conf.PositionAt(ix.atoms[3]);
                Vec3 p4 = conf.PositionAt(ix.atoms[4]);
                energy = EvalCMAP(p0, p1, p2, p3, p4,
                                  params.cmap_grids[cmap_idx],
                                  params.cmap_grid_spacing);
                target = &result->cmap_energy_;
                ++count_cmap;
                break;
            }
        }

        if (!target) continue;

        // Split energy evenly among participating atoms
        double share = energy / ix.n_atoms;
        for (int k = 0; k < ix.n_atoms; ++k)
            (*target)[ix.atoms[k]] += share;
    }

    // Total bonded = sum of all types
    for (size_t i = 0; i < N; ++i) {
        result->total_bonded_[i] =
            result->bond_energy_[i] + result->angle_energy_[i] +
            result->ub_energy_[i] + result->proper_energy_[i] +
            result->improper_energy_[i] + result->cmap_energy_[i];
    }

    OperationLog::Info(LogCalcOther, "BondedEnergyResult",
        "interactions: bond=" + std::to_string(count_bond) +
        " angle=" + std::to_string(count_angle) +
        " UB=" + std::to_string(count_ub) +
        " proper=" + std::to_string(count_proper) +
        " improper=" + std::to_string(count_improper) +
        " CMAP=" + std::to_string(count_cmap));

    return result;
}


// ── WriteFeatures ──────────────────────────────────────────────

int BondedEnergyResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    const size_t N = conf.AtomCount();
    assert(bond_energy_.size() == N);

    // Write (N, 7): bond, angle, UB, proper, improper, CMAP, total
    std::vector<double> data(N * 7);
    for (size_t i = 0; i < N; ++i) {
        data[i * 7 + 0] = bond_energy_[i];
        data[i * 7 + 1] = angle_energy_[i];
        data[i * 7 + 2] = ub_energy_[i];
        data[i * 7 + 3] = proper_energy_[i];
        data[i * 7 + 4] = improper_energy_[i];
        data[i * 7 + 5] = cmap_energy_[i];
        data[i * 7 + 6] = total_bonded_[i];
    }
    NpyWriter::WriteFloat64(output_dir + "/bonded_energy.npy",
                            data.data(), N, 7);
    return 1;
}

}  // namespace nmr

#include "ApbsFieldResult.h"
#include "Protein.h"
#include "PhysicalConstants.h"
#include "RuntimeEnvironment.h"
#include "NpyWriter.h"
#include "OperationLog.h"

// C bridge for APBS — isolates FETK's Vec3/Mat3 from Eigen's.
extern "C" {
#include "apbs_bridge.h"
}

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

namespace nmr {

// ============================================================================
// Grid interpolation utilities (from old ApbsSolver — no APBS dependency)
// ============================================================================

struct GridCache {
    Vec3 origin = Vec3::Zero();
    Vec3 spacing = Vec3::Zero();
    int dims[3] = {0, 0, 0};
    std::vector<double> data;
    bool valid = false;

    double Interpolate(const Vec3& point) const {
        if (!valid) return 0.0;

        Vec3 frac;
        for (int d = 0; d < 3; ++d)
            frac(d) = (point(d) - origin(d)) / spacing(d);

        // floor() not static_cast<int>(): truncation toward zero gives wrong
        // grid cell for negative fractional coordinates.
        int ix = static_cast<int>(std::floor(frac(0)));
        int iy = static_cast<int>(std::floor(frac(1)));
        int iz = static_cast<int>(std::floor(frac(2)));

        if (ix < 0 || ix >= dims[0]-1 || iy < 0 || iy >= dims[1]-1 ||
            iz < 0 || iz >= dims[2]-1)
            return 0.0;

        double fx = frac(0) - ix;
        double fy = frac(1) - iy;
        double fz = frac(2) - iz;

        auto idx = [&](int x, int y, int z) -> int {
            return x + y * dims[0] + z * dims[0] * dims[1];
        };

        return data[idx(ix,iy,iz)]     * (1-fx)*(1-fy)*(1-fz)
             + data[idx(ix+1,iy,iz)]   * fx*(1-fy)*(1-fz)
             + data[idx(ix,iy+1,iz)]   * (1-fx)*fy*(1-fz)
             + data[idx(ix+1,iy+1,iz)] * fx*fy*(1-fz)
             + data[idx(ix,iy,iz+1)]   * (1-fx)*(1-fy)*fz
             + data[idx(ix+1,iy,iz+1)] * fx*(1-fy)*fz
             + data[idx(ix,iy+1,iz+1)] * (1-fx)*fy*fz
             + data[idx(ix+1,iy+1,iz+1)] * fx*fy*fz;
    }
};

static Vec3 ElectricFieldFromGrid(const GridCache& grid, const Vec3& point) {
    Vec3 E;
    for (int d = 0; d < 3; ++d) {
        Vec3 plus = point, minus = point;
        plus(d)  += grid.spacing(d);
        minus(d) -= grid.spacing(d);
        // E = -grad(phi)
        E(d) = -(grid.Interpolate(plus) - grid.Interpolate(minus))
               / (2.0 * grid.spacing(d));
    }
    return E;
}

static Mat3 FieldGradientFromGrid(const GridCache& grid, const Vec3& point) {
    Mat3 EFG;
    for (int j = 0; j < 3; ++j) {
        Vec3 plus = point, minus = point;
        plus(j)  += grid.spacing(j);
        minus(j) -= grid.spacing(j);
        Vec3 Eplus  = ElectricFieldFromGrid(grid, plus);
        Vec3 Eminus = ElectricFieldFromGrid(grid, minus);
        for (int i = 0; i < 3; ++i)
            EFG(i, j) = (Eplus(i) - Eminus(i)) / (2.0 * grid.spacing(j));
    }

    // Traceless projection: remove the self-potential Laplacian.
    //
    // The APBS potential includes each atom's own Coulomb field, whose
    // Laplacian is -(q/epsilon)delta(r-r_i) — a delta function that the grid
    // discretizes into a large finite trace. The EFG from all EXTERNAL
    // sources (other atoms + solvent reaction field) satisfies Laplace's
    // equation and IS traceless. Subtracting trace/3 from the diagonal
    // removes exactly the self-interaction artifact.
    double trace = EFG.trace();
    EFG -= (trace / 3.0) * Mat3::Identity();

    return EFG;
}


// ============================================================================
// APBS solve path: calls the C bridge, extracts E-field and EFG per atom
// ============================================================================

static bool ComputeViaApbs(ProteinConformation& conf) {
    const Protein& protein = conf.ProteinRef();
    const size_t n_atoms = conf.AtomCount();

    // Separate x, y, z arrays for the C bridge
    std::vector<double> xArr(n_atoms), yArr(n_atoms), zArr(n_atoms);
    std::vector<double> charges(n_atoms), radii(n_atoms);

    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 pos = conf.PositionAt(i);
        xArr[i] = pos.x();
        yArr[i] = pos.y();
        zArr[i] = pos.z();
        charges[i] = conf.AtomAt(i).partial_charge;

        double r = conf.AtomAt(i).vdw_radius;
        if (r < 0.5) {
            // VdW radius not set or too small — use element defaults.
            // These approximate ff14SB typical values.
            switch (protein.AtomAt(i).element) {
                case Element::C: r = 1.70; break;
                case Element::N: r = 1.63; break;
                case Element::O: r = 1.48; break;
                case Element::S: r = 1.78; break;
                case Element::H: r = 1.00; break;
                default:         r = 1.50; break;
            }
        }
        radii[i] = r;
    }

    // Grid sizing: matches APBS mg-auto convention.
    // Fine grid: extent + 40A (20A padding per side), minimum 40A per axis.
    // Coarse grid: fine + 30A (for boundary condition accuracy).
    // Grid dimension: 161 per axis, targeting ~0.3-0.5A spacing.
    Vec3 lo(xArr[0], yArr[0], zArr[0]);
    Vec3 hi = lo;
    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 p(xArr[i], yArr[i], zArr[i]);
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    Vec3 extent = hi - lo;

    double fine_dims[3], coarse_dims[3];
    for (int d = 0; d < 3; ++d) {
        fine_dims[d]   = std::max(extent(d) + 40.0, 40.0);
        coarse_dims[d] = fine_dims[d] + 30.0;
    }

    int grid_dim = 161;

    // Standard PB parameters
    double pdie = 4.0;              // protein interior dielectric
    double sdie = 78.54;            // solvent dielectric (water, 25C)
    double temperature = 298.15;    // Kelvin
    double ionic_strength = 0.15;   // molar (physiological)

    OperationLog::Log(OperationLog::Level::Info, LogAPBS,
        "ApbsFieldResult",
        "APBS solve: " + std::to_string(n_atoms) + " atoms, " +
        "grid=" + std::to_string(grid_dim) + "^3, " +
        "fine=[" + std::to_string(fine_dims[0]) + "," +
        std::to_string(fine_dims[1]) + "," +
        std::to_string(fine_dims[2]) + "]A, " +
        "pdie=" + std::to_string(pdie) + " sdie=" + std::to_string(sdie));

    // Call the C bridge
    ApbsGridResult gridResult;
    int rc = apbs_solve(
        static_cast<int>(n_atoms),
        xArr.data(), yArr.data(), zArr.data(),
        charges.data(), radii.data(),
        pdie, sdie,
        temperature,
        ionic_strength,
        grid_dim, grid_dim, grid_dim,
        fine_dims[0], fine_dims[1], fine_dims[2],
        coarse_dims[0], coarse_dims[1], coarse_dims[2],
        &gridResult
    );

    if (rc != APBS_BRIDGE_OK) {
        std::string msg = "APBS solve failed: " + std::string(gridResult.error_msg);
        OperationLog::Warn("ApbsFieldResult::Compute", msg);
        apbs_free_grid(&gridResult);
        return false;
    }

    // Cache the grid for field/gradient extraction
    GridCache grid;
    grid.origin = Vec3(gridResult.origin[0], gridResult.origin[1], gridResult.origin[2]);
    grid.spacing = Vec3(gridResult.spacing[0], gridResult.spacing[1], gridResult.spacing[2]);
    grid.dims[0] = gridResult.dims[0];
    grid.dims[1] = gridResult.dims[1];
    grid.dims[2] = gridResult.dims[2];
    grid.data.assign(gridResult.data, gridResult.data + gridResult.n_points);
    grid.valid = true;

    OperationLog::Log(OperationLog::Level::Info, LogAPBS,
        "ApbsFieldResult",
        "Grid cached: " + std::to_string(gridResult.n_points) + " points, " +
        "spacing=" + std::to_string(grid.spacing(0)) + "A");

    apbs_free_grid(&gridResult);

    // Extract per-atom E-field and EFG from the potential grid
    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 pos = conf.PositionAt(i);

        Vec3 E = ElectricFieldFromGrid(grid, pos);
        Mat3 EFG = FieldGradientFromGrid(grid, pos);

        // Sanitise
        bool has_nan = false;
        for (int d = 0; d < 3; ++d) {
            if (std::isnan(E(d)) || std::isinf(E(d))) { has_nan = true; break; }
        }
        if (has_nan) {
            E = Vec3::Zero();
            EFG = Mat3::Zero();
        } else {
            double E_mag = E.norm();
            if (E_mag > APBS_SANITY_LIMIT) {
                E *= APBS_SANITY_LIMIT / E_mag;
            }
        }

        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                if (std::isnan(EFG(a,b)) || std::isinf(EFG(a,b)))
                    EFG(a,b) = 0.0;

        // Convert from APBS native units kT/(e*A) to V/A for E-field,
        // kT/(e*A^2) to V/A^2 for EFG. This makes APBS fields directly
        // comparable to CoulombResult (also in V/A, V/A^2).
        E *= KT_OVER_E_298K;
        EFG *= KT_OVER_E_298K;

        // Store on ConformationAtom: both Mat3 AND SphericalTensor
        auto& ca = conf.MutableAtomAt(i);
        ca.apbs_efield = E;
        ca.apbs_efg = EFG;
        ca.apbs_efg_spherical = SphericalTensor::Decompose(EFG);
    }

    OperationLog::Log(OperationLog::Level::Info, LogAPBS,
        "ApbsFieldResult::Compute",
        "APBS complete: E-field and EFG extracted for " +
        std::to_string(n_atoms) + " atoms");

    return true;
}


// ============================================================================
// Main compute: try APBS, fall back to vacuum Coulomb on failure
// ============================================================================

std::unique_ptr<ApbsFieldResult> ApbsFieldResult::Compute(
        ProteinConformation& conf) {

    if (!conf.HasResult<ChargeAssignmentResult>()) {
        return std::make_unique<ApbsFieldResult>();
    }

    OperationLog::Scope scope("ApbsFieldResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    auto result = std::make_unique<ApbsFieldResult>();
    result->conf_ = &conf;

    // Try APBS Poisson-Boltzmann solve first
    bool apbs_ok = ComputeViaApbs(conf);

    if (!apbs_ok) {
        OperationLog::Error("ApbsFieldResult::Compute",
            "APBS failed. No fallback — solvated fields require a working PB solver.");
        return nullptr;
    }

    return result;
}


Vec3 ApbsFieldResult::ElectricFieldAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).apbs_efield;
}

Mat3 ApbsFieldResult::FieldGradientAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).apbs_efg;
}

SphericalTensor ApbsFieldResult::FieldGradientSphericalAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).apbs_efg_spherical;
}


// ============================================================================
// WriteFeatures: apbs_E (N,3), apbs_efg (N,9)
// ============================================================================

static void PackST(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int ApbsFieldResult::WriteFeatures(const ProteinConformation& conf,
                                    const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    // apbs_E: (N, 3) — solvated Poisson-Boltzmann E-field in V/A
    {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const Vec3& E = conf.AtomAt(i).apbs_efield;
            data[i*3+0] = E.x(); data[i*3+1] = E.y(); data[i*3+2] = E.z();
        }
        NpyWriter::WriteFloat64(output_dir + "/apbs_E.npy", data.data(), N, 3);
    }

    // apbs_efg: (N, 9) — solvated EFG as SphericalTensor [T0,T1[3],T2[5]]
    {
        std::vector<double> data(N * 9);
        for (size_t i = 0; i < N; ++i)
            PackST(conf.AtomAt(i).apbs_efg_spherical, &data[i*9]);
        NpyWriter::WriteFloat64(output_dir + "/apbs_efg.npy", data.data(), N, 9);
    }

    return 2;
}

}  // namespace nmr

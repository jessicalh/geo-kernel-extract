#include "ApbsFieldResult.h"
#include "Protein.h"
#include "CalculatorConfig.h"
#include "PhysicalConstants.h"
#include "RuntimeEnvironment.h"
#include "NpyWriter.h"
#include "OperationLog.h"

// C bridge for APBS — isolates APBS/FETK headers from Eigen headers.
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
// Grid interpolation utilities (no APBS dependency).
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

        // OOB -> 0. The default padding is intended to keep atom stencils
        // interior; this remains the guard for any out-of-grid interpolation.
        if (ix < 0 || ix >= dims[0]-1 || iy < 0 || iy >= dims[1]-1 ||
            iz < 0 || iz >= dims[2]-1)
            return 0.0;

        double fx = frac(0) - ix;
        double fy = frac(1) - iy;
        double fz = frac(2) - iz;

        auto idx = [&](int x, int y, int z) -> int {
            return x + y * dims[0] + z * dims[0] * dims[1];
        };

        // trilinear interpolation (8 corner weights)
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

// E = -grad(phi) by central difference.
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
    // EFG(i,j) = dE_i/dr_j  (j = differentiation axis, i = field component).
    for (int j = 0; j < 3; ++j) {
        Vec3 plus = point, minus = point;
        plus(j)  += grid.spacing(j);
        minus(j) -= grid.spacing(j);
        Vec3 Eplus  = ElectricFieldFromGrid(grid, plus);
        Vec3 Eminus = ElectricFieldFromGrid(grid, minus);
        for (int i = 0; i < 3; ++i)
            EFG(i, j) = (Eplus(i) - Eminus(i)) / (2.0 * grid.spacing(j));
    }

    // Symmetrize. The independent finite-difference construction of
    // dE_i/dr_j and dE_j/dr_i can leave a tiny antisymmetric residue
    // (interpolation noise, not physics). Since the emit is T2-only,
    // Decompose would otherwise carry that residue as a spurious T1
    // pseudovector; explicit symmetrization pins T1 = 0 by construction.
    EFG = 0.5 * (EFG + EFG.transpose());

    // Traceless projection: remove the self-potential Laplacian.
    //
    // The APBS potential includes each atom's own Coulomb field, whose
    // Laplacian is -(q/epsilon)delta(r-r_i) — a delta function that the grid
    // discretizes into a large finite trace. The EFG from all EXTERNAL
    // sources (other atoms + solvent reaction field) satisfies Laplace's
    // equation and IS traceless. Subtracting trace/3 from the diagonal
    // removes that trace artifact.
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
    const auto& charge_result = conf.Result<ChargeAssignmentResult>();
    const int non_authoritative_radii =
        charge_result.ChargeTable().NonAuthoritativePbRadiusCount();
    if (non_authoritative_radii > 0) {
        OperationLog::Warn("ApbsFieldResult::Compute",
            "APBS ran with " + std::to_string(non_authoritative_radii) +
            " atoms on the flat 1.5 A placeholder PB radius "
            "(kCompatibilityPlaceholderPbRadiusAngstrom; real per-element "
            "ff14SB->PB radii are not yet wired). The PB dielectric boundary "
            "is therefore placeholder-quality: this frame's apbs_E / apbs_efg "
            "are fully populated and finite but NOT physically validated.");
    }

    // No atoms: the bbox below seeds from x_coords[0], and there is nothing to
    // solve. Guard placed after dependency retrieval so that path is unchanged.
    if (n_atoms == 0) return false;

    // Separate x, y, z arrays for the C bridge
    std::vector<double> x_coords(n_atoms), y_coords(n_atoms), z_coords(n_atoms);
    std::vector<double> charges(n_atoms), radii(n_atoms);

    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 pos = conf.PositionAt(i);
        x_coords[i] = pos.x();
        y_coords[i] = pos.y();
        z_coords[i] = pos.z();
        charges[i] = conf.AtomAt(i).partial_charge;

        double r = conf.AtomAt(i).pb_radius;
        if (r <= 0.0) {
            OperationLog::Error("ApbsFieldResult::Compute",
                "missing PB radius for atom " + std::to_string(i) +
                " (" + protein.AtomAt(i).pdb_atom_name + ")");
            return false;
        }
        radii[i] = r;
    }

    // Grid sizing follows the configured APBS padding conventions.
    // fine_dims:   extent + configured fine padding (default 40 A total),
    //              with configured minimum dimension per axis.
    // coarse_dims: fine + configured coarse padding (default 30 A).
    // grid_dim:    configured points per axis (default 161).
    int grid_dim = static_cast<int>(CalculatorConfig::Get("apbs_grid_dim"));

    Vec3 bbox_min(x_coords[0], y_coords[0], z_coords[0]);
    Vec3 bbox_max = bbox_min;
    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 p(x_coords[i], y_coords[i], z_coords[i]);
        bbox_min = bbox_min.cwiseMin(p);
        bbox_max = bbox_max.cwiseMax(p);
    }
    Vec3 extent = bbox_max - bbox_min;

    const double fine_padding   = CalculatorConfig::Get("apbs_fine_padding_A");
    const double fine_min_dim    = CalculatorConfig::Get("apbs_fine_min_dim_A");
    const double coarse_padding  = CalculatorConfig::Get("apbs_coarse_padding_A");
    double fine_dims[3], coarse_dims[3];
    for (int d = 0; d < 3; ++d) {
        fine_dims[d]   = std::max(extent(d) + fine_padding, fine_min_dim);
        coarse_dims[d] = fine_dims[d] + coarse_padding;
    }

    // Standard PB parameters
    double pdie = CalculatorConfig::Get("apbs_protein_dielectric");   // protein interior dielectric
    double sdie = CalculatorConfig::Get("apbs_solvent_dielectric");   // solvent dielectric (water, 25C)
    double temperature = CalculatorConfig::Get("apbs_temperature_K"); // Kelvin
    double ionic_strength = CalculatorConfig::Get("apbs_ionic_strength_M"); // molar (physiological)

    OperationLog::Log(OperationLog::Level::Info, LogAPBS,
        "ApbsFieldResult",
        "APBS solve: " + std::to_string(n_atoms) + " atoms, " +
        "grid=" + std::to_string(grid_dim) + "^3, " +
        "fine=[" + std::to_string(fine_dims[0]) + "," +
        std::to_string(fine_dims[1]) + "," +
        std::to_string(fine_dims[2]) + "]A, " +
        "pdie=" + std::to_string(pdie) + " sdie=" + std::to_string(sdie));

    // Call the C bridge
    ApbsGridResult apbs_grid;
    int rc = apbs_solve(
        static_cast<int>(n_atoms),
        x_coords.data(), y_coords.data(), z_coords.data(),
        charges.data(), radii.data(),
        pdie, sdie,
        temperature,
        ionic_strength,
        grid_dim, grid_dim, grid_dim,
        fine_dims[0], fine_dims[1], fine_dims[2],
        coarse_dims[0], coarse_dims[1], coarse_dims[2],
        &apbs_grid
    );

    if (rc != APBS_BRIDGE_OK) {
        std::string msg = "APBS solve failed: " + std::string(apbs_grid.error_msg);
        OperationLog::Warn("ApbsFieldResult::Compute", msg);
        apbs_free_grid(&apbs_grid);
        return false;
    }

    // Cache the grid for field/gradient extraction
    GridCache grid;
    grid.origin = Vec3(apbs_grid.origin[0], apbs_grid.origin[1], apbs_grid.origin[2]);
    grid.spacing = Vec3(apbs_grid.spacing[0], apbs_grid.spacing[1], apbs_grid.spacing[2]);
    grid.dims[0] = apbs_grid.dims[0];
    grid.dims[1] = apbs_grid.dims[1];
    grid.dims[2] = apbs_grid.dims[2];
    grid.data.assign(apbs_grid.data, apbs_grid.data + apbs_grid.n_points);
    grid.valid = true;

    OperationLog::Log(OperationLog::Level::Info, LogAPBS,
        "ApbsFieldResult",
        "Grid cached: " + std::to_string(apbs_grid.n_points) + " points, " +
        "spacing=" + std::to_string(grid.spacing(0)) + "A");

    apbs_free_grid(&apbs_grid);

    // Extract per-atom E-field and EFG from the potential grid
    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 pos = conf.PositionAt(i);

        Vec3 E = ElectricFieldFromGrid(grid, pos);
        Mat3 EFG = FieldGradientFromGrid(grid, pos);

        // finite-value guard (zero non-finite E + EFG)
        bool has_nan = false;
        for (int d = 0; d < 3; ++d) {
            if (std::isnan(E(d)) || std::isinf(E(d))) { has_nan = true; break; }
        }
        if (has_nan) {
            E = Vec3::Zero();
            EFG = Mat3::Zero();
        } else {
            // field magnitude cap: a finite-value guard on E only. EFG is
            // left to the traceless projection and its own finite-value guard;
            // coupling the rescale would invent a physical relationship that
            // is not in the code or the PB model.
            //
            // Intentionally NOT the shared efield_magnitude_sanity_clamp
            // config key (used by CoulombResult / WaterFieldResult /
            // MopacCoulombResult): that key is in V/A, but here E is still in
            // APBS-native kT/(e*A) units — the KT_OVER_E_298K conversion is
            // below. The same number means different physical magnitudes in
            // the two unit systems, so this guard keeps its own constant.
            // (Whether APBS should instead clamp post-conversion in V/A to
            // truly match the siblings is a physics/blessing question, not a
            // config-plumbing one.)
            double E_mag = E.norm();
            if (E_mag > APBS_SANITY_LIMIT) {
                E *= APBS_SANITY_LIMIT / E_mag;
            }
        }

        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                if (std::isnan(EFG(a,b)) || std::isinf(EFG(a,b)))
                    EFG(a,b) = 0.0;

        // kT/(e*A) -> V/A ; kT/(e*A^2) -> V/A^2  (comparable to CoulombResult)
        E *= KT_OVER_E_298K;
        EFG *= KT_OVER_E_298K;

        // store Mat3 + T2
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
// Main compute: APBS solve only. No fallback — nullptr on failure.
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

    // Run the APBS Poisson-Boltzmann solve.
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
// Feature export: apbs_E (N,3), apbs_efg (N,5) — T2 only (T0+T1 structural
// zeros). T2 components emitted inline below.
// ============================================================================

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

    // apbs_efg: (N, 5) — T2 only. The EFG is symmetrized + traceless-
    // projected above, so T0 and T1 are structural zeros.
    {
        std::vector<double> data(N * 5);
        for (size_t i = 0; i < N; ++i) {
            conf.AtomAt(i).apbs_efg_spherical.PackT2(&data[i*5]);
        }
        NpyWriter::WriteFloat64(output_dir + "/apbs_efg.npy", data.data(), N, 5);
    }

    return 2;
}

}  // namespace nmr

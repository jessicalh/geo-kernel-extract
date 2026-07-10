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
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

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

static Vec3 SanitizeApbsVector(Vec3 value,
                               std::size_t atom_index,
                               const char* field_name,
                               std::uint8_t mask_bit,
                               std::uint8_t& sanitizer_mask) {
    for (int d = 0; d < 3; ++d) {
        if (!std::isfinite(value(d))) {
            sanitizer_mask |= mask_bit;
            OperationLog::Warn("ApbsFieldResult::Compute",
                std::string("non-finite ") + field_name +
                " at atom " + std::to_string(atom_index) +
                "; applying established APBS whole-vector zero sanitizer");
            return Vec3::Zero();
        }
    }
    return value;
}

static Mat3 SanitizeApbsMatrix(Mat3 value,
                               std::size_t atom_index,
                               const char* field_name,
                               std::uint8_t mask_bit,
                               std::uint8_t& sanitizer_mask) {
    int replaced = 0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            if (!std::isfinite(value(a, b))) {
                value(a, b) = 0.0;
                ++replaced;
            }
        }
    }
    if (replaced != 0) {
        sanitizer_mask |= mask_bit;
        OperationLog::Warn("ApbsFieldResult::Compute",
            std::string("non-finite ") + field_name +
            " entries at atom " + std::to_string(atom_index) +
            "; replaced " + std::to_string(replaced) +
            " entries with zero under established APBS matrix sanitizer");
    }
    return value;
}

static bool SameApbsGridGeometry(const ApbsGridResult& total,
                                 const ApbsGridResult& reference,
                                 std::string& reason) {
    for (int d = 0; d < 3; ++d) {
        if (total.dims[d] <= 1 || reference.dims[d] <= 1) {
            reason = "invalid returned grid dimension on axis " +
                     std::to_string(d);
            return false;
        }
        if (total.dims[d] != reference.dims[d]) {
            reason = "dimension mismatch on axis " + std::to_string(d) +
                     ": total=" + std::to_string(total.dims[d]) +
                     " reference=" + std::to_string(reference.dims[d]);
            return false;
        }
        if (!std::isfinite(total.origin[d]) ||
            !std::isfinite(reference.origin[d]) ||
            std::abs(total.origin[d] - reference.origin[d]) > 1e-12) {
            reason = "origin mismatch/nonfinite on axis " + std::to_string(d);
            return false;
        }
        if (!std::isfinite(total.spacing[d]) ||
            !std::isfinite(reference.spacing[d]) ||
            std::abs(total.spacing[d] - reference.spacing[d]) > 1e-12) {
            reason = "spacing mismatch/nonfinite on axis " + std::to_string(d);
            return false;
        }
    }
    if (total.n_points != reference.n_points) {
        reason = "point-count mismatch: total=" +
                 std::to_string(total.n_points) + " reference=" +
                 std::to_string(reference.n_points);
        return false;
    }
    if (total.n_points <= 0 || !total.data || !reference.data) {
        reason = "successful solve returned empty grid data";
        return false;
    }
    const std::uint64_t expected_points =
        static_cast<std::uint64_t>(total.dims[0]) *
        static_cast<std::uint64_t>(total.dims[1]) *
        static_cast<std::uint64_t>(total.dims[2]);
    if (expected_points != static_cast<std::uint64_t>(total.n_points)) {
        reason = "returned point count does not equal dimensions product";
        return false;
    }
    return true;
}

static GridCache CacheGrid(const ApbsGridResult& source) {
    GridCache grid;
    grid.origin = Vec3(source.origin[0], source.origin[1], source.origin[2]);
    grid.spacing = Vec3(source.spacing[0], source.spacing[1], source.spacing[2]);
    grid.dims[0] = source.dims[0];
    grid.dims[1] = source.dims[1];
    grid.dims[2] = source.dims[2];
    grid.data.assign(source.data, source.data + source.n_points);
    grid.valid = true;
    return grid;
}

bool ApbsFieldResult::ComputeViaApbs(ProteinConformation& conf,
                                     ApbsFieldResult& result) {
    const Protein& protein = conf.ProteinRef();
    const size_t n_atoms = conf.AtomCount();
    const auto& charge_result = conf.Result<ChargeAssignmentResult>();
    const ForceFieldChargeTable& charge_table = charge_result.ChargeTable();
    const int matched_radii = charge_table.MatchedPbRadiusCount();
    const int derived_radii = charge_table.DerivedMbondi2PbRadiusCount();
    const int placeholder_radii = charge_table.PlaceholderPbRadiusCount();
    const int missing_radii = charge_table.MissingCount();
    if (placeholder_radii + missing_radii > 0) {
        OperationLog::Error("ApbsFieldResult::Compute",
            "PB radius provenance blocks APBS: matched=" +
            std::to_string(matched_radii) +
            " derived_mbondi2=" + std::to_string(derived_radii) +
            " placeholder=" + std::to_string(placeholder_radii) +
            " missing=" + std::to_string(missing_radii));
        return false;
    }
    if (derived_radii > 0) {
        OperationLog::Info(LogAPBS, "ApbsFieldResult::Compute",
            "PB radius provenance: matched=" +
            std::to_string(matched_radii) +
            " derived_mbondi2=" + std::to_string(derived_radii) +
            " placeholder=0 missing=0");
    }

    // No atoms: the bbox below seeds from x_coords[0], and there is nothing to
    // solve. Guard placed after dependency retrieval so that path is unchanged.
    if (n_atoms == 0) return false;

    // Separate x, y, z arrays for the C bridge
    std::vector<double> x_coords(n_atoms), y_coords(n_atoms), z_coords(n_atoms);
    std::vector<double> charges(n_atoms), radii(n_atoms);

    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 pos = conf.PositionAt(i);
        if (!pos.allFinite()) {
            OperationLog::Error("ApbsFieldResult::Compute",
                "non-finite position at atom " + std::to_string(i));
            return false;
        }
        x_coords[i] = pos.x();
        y_coords[i] = pos.y();
        z_coords[i] = pos.z();
        charges[i] = conf.AtomAt(i).partial_charge;
        if (!std::isfinite(charges[i])) {
            OperationLog::Error("ApbsFieldResult::Compute",
                "non-finite partial charge at atom " + std::to_string(i));
            return false;
        }

        double r = conf.AtomAt(i).pb_radius;
        if (!std::isfinite(r) || r <= 0.0) {
            OperationLog::Error("ApbsFieldResult::Compute",
                "invalid PB radius for atom " + std::to_string(i) +
                " (" + protein.AtomAt(i).pdb_atom_name + ")");
            return false;
        }
        radii[i] = r;
    }

    // One honest manual grid: bbox extent + configured total padding,
    // floored by the configured minimum length on each axis.
    const int grid_dim =
        static_cast<int>(CalculatorConfig::Get("apbs_grid_dim"));
    if (grid_dim <= 1) {
        OperationLog::Error("ApbsFieldResult::Compute",
            "apbs_grid_dim must be greater than one");
        return false;
    }

    Vec3 bbox_min(x_coords[0], y_coords[0], z_coords[0]);
    Vec3 bbox_max = bbox_min;
    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 p(x_coords[i], y_coords[i], z_coords[i]);
        bbox_min = bbox_min.cwiseMin(p);
        bbox_max = bbox_max.cwiseMax(p);
    }
    Vec3 extent = bbox_max - bbox_min;

    const double manual_padding =
        CalculatorConfig::Get("apbs_manual_grid_padding_A");
    const double manual_min_dim =
        CalculatorConfig::Get("apbs_manual_grid_min_dim_A");
    if (!std::isfinite(manual_padding) || manual_padding < 0.0 ||
        !std::isfinite(manual_min_dim) || manual_min_dim <= 0.0) {
        OperationLog::Error("ApbsFieldResult::Compute",
            "invalid APBS manual-grid padding or minimum dimension");
        return false;
    }
    std::array<double, 3> grid_lengths{};
    for (int d = 0; d < 3; ++d) {
        grid_lengths[static_cast<std::size_t>(d)] =
            std::max(extent(d) + manual_padding, manual_min_dim);
    }

    // Standard PB parameters
    const double pdie = CalculatorConfig::Get("apbs_protein_dielectric");
    const double sdie = CalculatorConfig::Get("apbs_solvent_dielectric");
    const double temperature = CalculatorConfig::Get("apbs_temperature_K");
    const double ionic_strength =
        CalculatorConfig::Get("apbs_ionic_strength_M");
    const double thermal_voltage = ApbsThermalVoltage(temperature);
    const double efield_clamp =
        CalculatorConfig::Get("efield_magnitude_sanity_clamp");
    if (!std::isfinite(pdie) || pdie <= 0.0 ||
        !std::isfinite(sdie) || sdie <= 0.0 ||
        !std::isfinite(temperature) || temperature <= 0.0 ||
        !std::isfinite(ionic_strength) || ionic_strength < 0.0 ||
        !std::isfinite(thermal_voltage) || thermal_voltage <= 0.0 ||
        !std::isfinite(efield_clamp) || efield_clamp < 0.0) {
        OperationLog::Error("ApbsFieldResult::Compute",
            "invalid APBS dielectric, temperature, ionic strength, or clamp");
        return false;
    }

    OperationLog::Log(OperationLog::Level::Info, LogAPBS,
        "ApbsFieldResult",
        "APBS solve: " + std::to_string(n_atoms) + " atoms, " +
        "grid=" + std::to_string(grid_dim) + "^3, " +
        "manual_lengths=[" + std::to_string(grid_lengths[0]) + "," +
        std::to_string(grid_lengths[1]) + "," +
        std::to_string(grid_lengths[2]) + "]A, " +
        "pdie=" + std::to_string(pdie) + " sdie=" + std::to_string(sdie));

    std::array<double, 2> ion_concentrations = {
        ionic_strength, ionic_strength};
    std::array<double, 2> ion_radii = {0.95, 1.81};
    std::array<double, 2> ion_charges = {1.0, -1.0};

    const ApbsSolveParams total_params = {
        grid_dim, grid_dim, grid_dim,
        grid_lengths[0], grid_lengths[1], grid_lengths[2],
        pdie, sdie, temperature,
        2,
        ion_concentrations.data(), ion_radii.data(), ion_charges.data()};

    ApbsSolveParams reference_params = total_params;
    reference_params.pdie = 1.0;
    reference_params.sdie = 1.0;
    reference_params.mobile_ion_count = 0;
    reference_params.mobile_ion_conc_M = nullptr;
    reference_params.mobile_ion_radius_A = nullptr;
    reference_params.mobile_ion_charge_e = nullptr;

    ApbsGridResult total_grid{};
    int rc = apbs_solve(
        static_cast<int>(n_atoms),
        x_coords.data(), y_coords.data(), z_coords.data(),
        charges.data(), radii.data(),
        &total_params, &total_grid);

    if (rc != APBS_BRIDGE_OK) {
        const std::string msg = "APBS total-PB solve failed: " +
            std::string(total_grid.error_msg);
        OperationLog::Warn("ApbsFieldResult::Compute", msg);
        apbs_free_grid(&total_grid);
        return false;
    }

    ApbsGridResult reference_grid{};
    rc = apbs_solve(
        static_cast<int>(n_atoms),
        x_coords.data(), y_coords.data(), z_coords.data(),
        charges.data(), radii.data(),
        &reference_params, &reference_grid);

    if (rc != APBS_BRIDGE_OK) {
        const std::string msg =
            "APBS homogeneous-vacuum reference solve failed: " +
            std::string(reference_grid.error_msg);
        OperationLog::Warn("ApbsFieldResult::Compute", msg);
        apbs_free_grid(&total_grid);
        apbs_free_grid(&reference_grid);
        return false;
    }

    std::string geometry_mismatch;
    if (!SameApbsGridGeometry(total_grid, reference_grid,
                              geometry_mismatch)) {
        OperationLog::Error("ApbsFieldResult::Compute",
            "APBS total/reference grid geometry mismatch: " +
            geometry_mismatch);
        apbs_free_grid(&total_grid);
        apbs_free_grid(&reference_grid);
        return false;
    }

    GridCache total_cache = CacheGrid(total_grid);
    GridCache reaction_cache = total_cache;
    for (int k = 0; k < total_grid.n_points; ++k) {
        reaction_cache.data[static_cast<std::size_t>(k)] =
            total_grid.data[k] - reference_grid.data[k];
    }

    result.grid_dims_ = {
        static_cast<std::uint64_t>(total_grid.dims[0]),
        static_cast<std::uint64_t>(total_grid.dims[1]),
        static_cast<std::uint64_t>(total_grid.dims[2])};
    result.grid_lengths_A_ = grid_lengths;
    result.grid_origin_A_ = {
        total_grid.origin[0], total_grid.origin[1], total_grid.origin[2]};
    result.grid_spacing_A_ = {
        total_grid.spacing[0], total_grid.spacing[1], total_grid.spacing[2]};
    result.manual_grid_padding_A_ = manual_padding;
    result.manual_grid_min_dim_A_ = manual_min_dim;
    result.temperature_K_ = temperature;
    result.thermal_voltage_V_ = thermal_voltage;
    result.efield_clamp_threshold_ = efield_clamp;

    OperationLog::Log(OperationLog::Level::Info, LogAPBS,
        "ApbsFieldResult",
        "matched total/reference manual grid: dims=[" +
        std::to_string(total_grid.dims[0]) + "," +
        std::to_string(total_grid.dims[1]) + "," +
        std::to_string(total_grid.dims[2]) + "] origin_A=[" +
        std::to_string(total_grid.origin[0]) + "," +
        std::to_string(total_grid.origin[1]) + "," +
        std::to_string(total_grid.origin[2]) + "] spacing_A=[" +
        std::to_string(total_grid.spacing[0]) + "," +
        std::to_string(total_grid.spacing[1]) + "," +
        std::to_string(total_grid.spacing[2]) + "]");

    apbs_free_grid(&total_grid);
    apbs_free_grid(&reference_grid);

    // Compute every per-atom value here. WriteFeatures is pure read-back.
    std::vector<Vec3> reaction_E(n_atoms);
    std::vector<Mat3> reaction_EFG(n_atoms);
    std::vector<double> reaction_phi(n_atoms);
    std::vector<std::uint8_t> clamp_mask(n_atoms, 0);
    std::vector<double> clamp_scale(n_atoms, 1.0);
    std::vector<std::uint8_t> nonfinite_sanitizer_mask(n_atoms, 0);
    std::vector<Vec3> total_E(n_atoms);
    std::vector<Mat3> total_EFG(n_atoms);

    for (size_t i = 0; i < n_atoms; ++i) {
        const Vec3 pos = conf.PositionAt(i);

        Vec3 E = ElectricFieldFromGrid(reaction_cache, pos) * thermal_voltage;
        // FieldGradientFromGrid returns the field gradient dE_i/dr_j = -d2(phi).
        // The Coulomb-family EFG (MOPAC/FF/AIMNet) is the potential Hessian +d2(phi)
        // = -dE/dr, the standard nuclear-EFG convention.
        Mat3 EFG =
            -FieldGradientFromGrid(reaction_cache, pos) * thermal_voltage;
        Vec3 E_total =
            ElectricFieldFromGrid(total_cache, pos) * thermal_voltage;
        Mat3 EFG_total =
            -FieldGradientFromGrid(total_cache, pos) * thermal_voltage;

        E = SanitizeApbsVector(
            E, i, "canonical reaction E-field", 1u << 0,
            nonfinite_sanitizer_mask[i]);
        EFG = SanitizeApbsMatrix(
            EFG, i, "canonical reaction EFG", 1u << 1,
            nonfinite_sanitizer_mask[i]);
        E_total = SanitizeApbsVector(
            E_total, i, "raw total-PB diagnostic E-field", 1u << 2,
            nonfinite_sanitizer_mask[i]);
        EFG_total = SanitizeApbsMatrix(
            EFG_total, i, "raw total-PB diagnostic EFG", 1u << 3,
            nonfinite_sanitizer_mask[i]);

        const double phi =
            reaction_cache.Interpolate(pos) * thermal_voltage;
        if (!std::isfinite(phi)) {
            OperationLog::Error("ApbsFieldResult::Compute",
                "non-finite reaction potential at atom " +
                std::to_string(i) +
                "; refusing to attach APBS result (no silent zero sanitizer)");
            return false;
        }

        const double original_norm = E.norm();
        if (original_norm > efield_clamp) {
            clamp_scale[i] = efield_clamp / original_norm;
            E *= clamp_scale[i];
            clamp_mask[i] = 1;
        }

        reaction_E[i] = E;
        reaction_EFG[i] = EFG;
        reaction_phi[i] = phi;
        total_E[i] = E_total;
        total_EFG[i] = EFG_total;
    }

    for (size_t i = 0; i < n_atoms; ++i) {
        auto& ca = conf.MutableAtomAt(i);
        ca.apbs_efield = reaction_E[i];
        ca.apbs_efg = reaction_EFG[i];
        ca.apbs_efg_spherical = SphericalTensor::Decompose(reaction_EFG[i]);
        ca.apbs_phi = reaction_phi[i];
        ca.apbs_efield_clamp_mask = clamp_mask[i];
        ca.apbs_efield_clamp_scale = clamp_scale[i];
        ca.apbs_nonfinite_sanitizer_mask = nonfinite_sanitizer_mask[i];
        ca.apbs_efield_total_diagnostic = total_E[i];
        ca.apbs_efg_total_diagnostic = total_EFG[i];
        ca.apbs_efg_total_diagnostic_spherical =
            SphericalTensor::Decompose(total_EFG[i]);
    }

    OperationLog::Log(OperationLog::Level::Info, LogAPBS,
        "ApbsFieldResult::Compute",
        "APBS complete: reaction potential/E/EFG and total diagnostics for " +
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
    bool apbs_ok = ApbsFieldResult::ComputeViaApbs(conf, *result);

    if (!apbs_ok) {
        OperationLog::Error("ApbsFieldResult::Compute",
            "APBS failed. No fallback — reaction fields require a working PB solver.");
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
// Feature export is pure read-back from ConformationAtom fields populated by
// ComputeViaApbs. No APBS quantity is re-derived here.
// ============================================================================

int ApbsFieldResult::WriteFeatures(const ProteinConformation& conf,
                                    const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    OperationLog::Info(LogAPBS, "ApbsFieldResult::WriteFeatures",
        "reaction-field schema: grid_mode=single_manual dims=[" +
        std::to_string(grid_dims_[0]) + "," +
        std::to_string(grid_dims_[1]) + "," +
        std::to_string(grid_dims_[2]) + "] lengths_A=[" +
        std::to_string(grid_lengths_A_[0]) + "," +
        std::to_string(grid_lengths_A_[1]) + "," +
        std::to_string(grid_lengths_A_[2]) + "] origin_A=[" +
        std::to_string(grid_origin_A_[0]) + "," +
        std::to_string(grid_origin_A_[1]) + "," +
        std::to_string(grid_origin_A_[2]) + "] spacing_A=[" +
        std::to_string(grid_spacing_A_[0]) + "," +
        std::to_string(grid_spacing_A_[1]) + "," +
        std::to_string(grid_spacing_A_[2]) + "] padding_A=" +
        std::to_string(manual_grid_padding_A_) + " min_dim_A=" +
        std::to_string(manual_grid_min_dim_A_) + " temperature_K=" +
        std::to_string(temperature_K_) + " thermal_voltage_V=" +
        std::to_string(thermal_voltage_V_) + " efield_clamp_V_per_A=" +
        std::to_string(efield_clamp_threshold_));

    // apbs_E: (N, 3) — canonical reaction Poisson-Boltzmann E-field in V/A
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

    // apbs_phi: (N,) — canonical reaction potential in V.
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).apbs_phi;
        NpyWriter::WriteFloat64(output_dir + "/apbs_phi.npy", data.data(), N);
    }

    // Clamp audit for canonical reaction E only.
    {
        std::vector<std::uint8_t> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).apbs_efield_clamp_mask;
        NpyWriter::WriteUInt8(output_dir + "/apbs_E_clamp_mask.npy",
                              data.data(), N);
    }
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).apbs_efield_clamp_scale;
        NpyWriter::WriteFloat64(output_dir + "/apbs_E_clamp_scale.npy",
                                data.data(), N);
    }

    // Bitwise audit of finite sanitization performed during Compute:
    // bit 0 reaction E, bit 1 reaction EFG, bit 2 total E, bit 3 total EFG.
    {
        std::vector<std::uint8_t> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).apbs_nonfinite_sanitizer_mask;
        NpyWriter::WriteUInt8(
            output_dir + "/apbs_nonfinite_sanitizer_mask.npy",
            data.data(), N);
    }

    // Raw total-PB derivatives are diagnostics: finite-sanitized, unclamped.
    {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const Vec3& E = conf.AtomAt(i).apbs_efield_total_diagnostic;
            data[i*3+0] = E.x();
            data[i*3+1] = E.y();
            data[i*3+2] = E.z();
        }
        NpyWriter::WriteFloat64(
            output_dir + "/apbs_E_total_diagnostic.npy",
            data.data(), N, 3);
    }
    {
        std::vector<double> data(N * 5);
        for (size_t i = 0; i < N; ++i) {
            conf.AtomAt(i).apbs_efg_total_diagnostic_spherical.PackT2(
                &data[i*5]);
        }
        NpyWriter::WriteFloat64(
            output_dir + "/apbs_efg_total_diagnostic.npy",
            data.data(), N, 5);
    }

    return 8;
}

}  // namespace nmr

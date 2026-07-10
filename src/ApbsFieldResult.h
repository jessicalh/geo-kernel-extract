#pragma once
//
// ApbsFieldResult: APBS reaction potential, E-field, and EFG at each atom.
//
// Dependencies: ChargeAssignmentResult (charges and PB radii).
//
// Canonical values are total linearized-PB minus a homogeneous-vacuum
// reference solve over the identical atom list and grid. Raw total-PB
// derivatives are retained only as explicitly named diagnostics.
//
// Linearised Poisson-Boltzmann solve via the APBS C bridge (apbs_bridge.h).
// The bridge takes in-memory arrays of positions, charges, and radii and
// returns a potential grid in kT/e. The total and reference grids are
// subtracted in APBS-native units before the canonical potential, E-field,
// and EFG are extracted.
//
// No fallback. If APBS fails, Compute() returns nullptr. Reaction fields
// require a working PB solver — substituting vacuum Coulomb would silently
// produce a different physical quantity in different units.
//
// E-field from grid: E = -grad(phi) via central differences.
// EFG is the sign-aligned potential Hessian, symmetrized and traceless-
// projected before storage.
//
// Units stored on ConformationAtom: V for potential, V/A for E-field,
// V/A^2 for EFG. The conversion uses the configured solve temperature.
//

#include "ConformationResult.h"
#include "ChargeAssignmentResult.h"
#include "PhysicalConstants.h"
#include "ProteinConformation.h"

#include <array>
#include <cstdint>
#include <limits>

namespace nmr {

inline double ApbsThermalVoltage(double temperature_K) {
    return KT_OVER_E_298K * temperature_K / 298.15;
}

class ApbsFieldResult : public ConformationResult {
public:
    std::string Name() const override { return "ApbsFieldResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return { std::type_index(typeid(ChargeAssignmentResult)) };
    }

    // Compute APBS E-field and EFG at each atom from charges.
    static std::unique_ptr<ApbsFieldResult> Compute(ProteinConformation& conf);

    // Query methods
    Vec3 ElectricFieldAt(size_t atom_index) const;
    Mat3 FieldGradientAt(size_t atom_index) const;
    SphericalTensor FieldGradientSphericalAt(size_t atom_index) const;

    const std::array<std::uint64_t, 3>& GridDims() const {
        return grid_dims_;
    }
    const std::array<double, 3>& GridLengthsA() const {
        return grid_lengths_A_;
    }
    const std::array<double, 3>& GridOriginA() const {
        return grid_origin_A_;
    }
    const std::array<double, 3>& GridSpacingA() const {
        return grid_spacing_A_;
    }
    double ManualGridPaddingA() const { return manual_grid_padding_A_; }
    double ManualGridMinDimA() const { return manual_grid_min_dim_A_; }
    double TemperatureK() const { return temperature_K_; }
    double ThermalVoltageV() const { return thermal_voltage_V_; }
    double EfieldClampThreshold() const { return efield_clamp_threshold_; }

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    static bool ComputeViaApbs(ProteinConformation& conf,
                               ApbsFieldResult& result);

    const ProteinConformation* conf_ = nullptr;
    std::array<std::uint64_t, 3> grid_dims_ = {0, 0, 0};
    std::array<double, 3> grid_lengths_A_ = {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    std::array<double, 3> grid_origin_A_ = grid_lengths_A_;
    std::array<double, 3> grid_spacing_A_ = grid_lengths_A_;
    double manual_grid_padding_A_ = std::numeric_limits<double>::quiet_NaN();
    double manual_grid_min_dim_A_ = std::numeric_limits<double>::quiet_NaN();
    double temperature_K_ = std::numeric_limits<double>::quiet_NaN();
    double thermal_voltage_V_ = std::numeric_limits<double>::quiet_NaN();
    double efield_clamp_threshold_ = std::numeric_limits<double>::quiet_NaN();
};

}  // namespace nmr

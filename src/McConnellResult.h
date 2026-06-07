#pragma once
//
// McConnellResult: clean neighbour magnetic-anisotropy tensor kernels.
//
// For each source bond within the configured cutoff (default 10 A), computes
// A = D(r) Qhat, where D(r) = (3 n n^T - I)/r^3 and
// Qhat = u u^T - I/3. Accumulates seven source categories in two channels:
// fixed source strength and MOPAC Wiberg bond-order strength.
//
// Emit surface: mc_<category>_<fixed|bo>.npy, one packed (N,9) array per
// category/channel in SphericalTensor::PackFull9 order. Units Angstrom^-3.
//

#include "ConformationResult.h"
#include "Types.h"

#include <nlohmann/json_fwd.hpp>

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

inline constexpr const char* kMcConnellPackFull9IrrepLayout =
    "0e,1e_x,1e_y,1e_z,2e_m-2..+2";

// NOTE: vestigial. The runtime cutoff comes from CalculatorConfig
// "mcconnell_bond_anisotropy_cutoff"; this constant is referenced only
// in comments/docs, never in the math path.
constexpr double MCCONNELL_CUTOFF_A = 10.0;


struct McConnellPairKernel {
    Mat3 response = Mat3::Zero();         // D(r) Qhat, Angstrom^-3
    Mat3 dipolar = Mat3::Zero();          // D(r), symmetric traceless
    Mat3 source_shape = Mat3::Zero();     // Qhat, unit susceptibility shape
    double pcs_scalar = 0.0;              // trace(D Qhat)/3 = n^T Qhat n/r^3
    double distance = 0.0;
    Vec3 direction = Vec3::Zero();        // source midpoint -> target
};

class McConnellResult : public ConformationResult {
public:
    std::string Name() const override { return "McConnellResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute McConnell tensors for all atoms.
    static std::unique_ptr<McConnellResult> Compute(
        ProteinConformation& conf);

    // Legacy query methods (return fixed-channel PCS scalar projections).
    double CategoryScalarSum(size_t atom_index, BondCategory cat) const;
    double NearestCOScalarContribution(size_t atom_index) const;

    // Grid sampling: evaluate McConnell kernel at an arbitrary 3D point.
    SphericalTensor SampleKernelAt(Vec3 point) const;

    // Test/support hook for the contract checks. This is the same pair
    // kernel used by Compute(); no SphericalTensor constants are duplicated.
    static McConnellPairKernel ComputePairKernel(
        const Vec3& atom_pos,
        const Vec3& source_center,
        const Vec3& source_axis);
    static McConnellPairKernel ComputePeptideCORhombicPairKernel(
        const ProteinConformation& conf,
        size_t bond_index,
        const Vec3& atom_pos);

    // Single source for extraction_manifest.json::feature_metadata.mcconnell.
    // TopologySidecar owns the top-level manifest file; McConnellResult owns
    // this calculator-specific block because it owns the emitted stems.
    static nlohmann::ordered_json FeatureMetadata(bool include_xh_sources);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr

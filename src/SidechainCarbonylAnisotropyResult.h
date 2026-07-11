#pragma once
//
// SidechainCarbonylAnisotropyResult
//
// Auditable side-chain C=O source inventory and local frames, paired with
// the canonical McConnell SidechainCO fixed/bond-order tensors.  The tensor
// values remain the established axial McConnell response; the raw local
// planes are emitted alongside them so no orientation information is lost.
//

#include "ConformationResult.h"
#include "SemanticEnums.h"
#include "Types.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <typeindex>
#include <vector>

namespace nmr {

class Protein;
class ProteinConformation;

namespace sidechain_carbonyl_anisotropy_detail {

// Stable emitted codes for the two typed oxygen predicates consumed here.
enum class OxygenSemanticClass : std::int32_t {
    Unknown = 0,
    SidechainAmide = 1,
    SidechainCarboxylate = 2,
};

struct SourceBond {
    std::size_t bond_index = SIZE_MAX;
    std::size_t carbon_atom = SIZE_MAX;
    std::size_t oxygen_atom = SIZE_MAX;
    std::size_t residue_index = SIZE_MAX;
    BondCategory bond_category = BondCategory::Unknown;
    PlanarGroupKind planar_group_kind = PlanarGroupKind::None;
    OxygenSemanticClass oxygen_semantic_class =
        OxygenSemanticClass::Unknown;
    std::size_t plane_reference_atom = SIZE_MAX;
    bool source_valid = false;
};

struct SourceFrame {
    Vec3 origin = Vec3::Constant(
        std::numeric_limits<double>::quiet_NaN());
    Vec3 x_axis = Vec3::Constant(
        std::numeric_limits<double>::quiet_NaN());
    Vec3 y_axis = Vec3::Constant(
        std::numeric_limits<double>::quiet_NaN());
    Vec3 z_axis = Vec3::Constant(
        std::numeric_limits<double>::quiet_NaN());
    double bond_length = std::numeric_limits<double>::quiet_NaN();
    double orthogonality_error = std::numeric_limits<double>::quiet_NaN();
    double normal_norm = std::numeric_limits<double>::quiet_NaN();
    bool frame_valid = false;
};

// Exact production classification and frame construction used by Compute().
// Tests call this named per-file surface directly rather than duplicating the
// semantic predicates or canonical axis gauge.
SourceBond ClassifySourceBond(const Protein& protein,
                              std::size_t bond_index);
SourceFrame BuildSourceFrame(const ProteinConformation& conf,
                             const SourceBond& source);

}  // namespace sidechain_carbonyl_anisotropy_detail


class SidechainCarbonylAnisotropyResult : public ConformationResult {
public:
    std::string Name() const override {
        return "SidechainCarbonylAnisotropyResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<SidechainCarbonylAnisotropyResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    std::vector<sidechain_carbonyl_anisotropy_detail::SourceBond> sources_;
    std::vector<sidechain_carbonyl_anisotropy_detail::SourceFrame> frames_;
    std::vector<SphericalTensor> fixed_tensors_;
    std::vector<SphericalTensor> bond_order_tensors_;
    std::vector<std::array<double, 4>> scalar_audit_;
    bool has_mopac_bond_orders_ = false;
};

}  // namespace nmr

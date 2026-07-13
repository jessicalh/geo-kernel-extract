#pragma once
//
// LocalBackboneGeometryResult: residue-local backbone valence geometry for
// one conformation.  The numerical kernels are shared with
// LocalBackboneGeometryTrajectoryResult; both surfaces therefore retain the
// established ideal-C-beta calculation and its global-Cartesian residual.
//

#include "ConformationResult.h"
#include "Types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;

namespace local_backbone_geometry {

struct Measurements {
    double tau_n_ca_c;
    double angle_n_ca_cb;
    double angle_cb_ca_c;
    double angle_cprev_n_ca;
    double angle_ca_c_nnext;
    double cb_deviation;
    Vec3 cb_residual_vector;
};

double BondAngleRadians(const Vec3& endpoint_a,
                        const Vec3& vertex,
                        const Vec3& endpoint_b);
Vec3 IdealCbPosition(const Vec3& n_pos,
                     const Vec3& ca_pos,
                     const Vec3& c_pos);
Vec3 CbDeviationVector(const Vec3& n_pos,
                       const Vec3& ca_pos,
                       const Vec3& c_pos,
                       const Vec3& observed_cb);
double CbDeviation(const Vec3& n_pos,
                   const Vec3& ca_pos,
                   const Vec3& c_pos,
                   const Vec3& observed_cb);

// Compute all seven measurements for one typed residue.  Backbone adjacency
// is resolved only through Protein::BackbonePredecessor/Successor.
Measurements MeasureResidue(const ProteinConformation& conf,
                            std::size_t residue_index);

}  // namespace local_backbone_geometry


class LocalBackboneGeometryResult : public ConformationResult {
public:
    std::string Name() const override { return "LocalBackboneGeometryResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<LocalBackboneGeometryResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    std::vector<double> tau_n_ca_c_;
    std::vector<double> angle_n_ca_cb_;
    std::vector<double> angle_cb_ca_c_;
    std::vector<double> angle_cprev_n_ca_;
    std::vector<double> angle_ca_c_nnext_;
    std::vector<double> cb_deviation_;
    std::vector<Vec3> cb_residual_vector_;

    std::vector<std::uint8_t> tau_n_ca_c_valid_;
    std::vector<std::uint8_t> angle_n_ca_cb_valid_;
    std::vector<std::uint8_t> angle_cb_ca_c_valid_;
    std::vector<std::uint8_t> angle_cprev_n_ca_valid_;
    std::vector<std::uint8_t> angle_ca_c_nnext_valid_;
    std::vector<std::uint8_t> cb_deviation_valid_;
    std::vector<std::uint8_t> cb_residual_vector_valid_;
};

}  // namespace nmr

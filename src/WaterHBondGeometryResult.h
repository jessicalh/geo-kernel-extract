#pragma once
//
// Explicit protein-water H-bond candidate geometry.  Waters remain in the
// per-frame SolventEnvironment; this result stores only raw pair rows and
// protein-atom summaries for NPY emission.

#include "ConformationResult.h"
#include "SolventEnvironment.h"
#include "Types.h"

#include <array>
#include <cstdint>
#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;
struct AtomSemanticTable;

namespace water_hbond_geometry_detail {

struct GeometryEvaluation {
    double donor_heavy_acceptor_distance_A = 0.0;
    double h_acceptor_distance_A = 0.0;
    double angle_deg = 0.0;
    bool passes_geometry = false;
};

// Exact production geometry used in both donation directions.
GeometryEvaluation EvaluateGeometry(const Vec3& donor_heavy,
                                    const Vec3& donor_hydrogen,
                                    const Vec3& acceptor);

// Typed acceptor projection shared with the emitted protein_role column:
// 0 none, 1 backbone carbonyl, 2 sidechain amide carbonyl,
// 3 carboxylate, 4 hydroxyl/oxide oxygen, 5 unprotonated ring N,
// 6 neutral other N/O/S.
std::uint8_t ProteinAcceptorClass(const AtomSemanticTable& semantic);

}  // namespace water_hbond_geometry_detail

class WaterHBondGeometryResult : public ConformationResult {
public:
    std::string Name() const override { return "WaterHBondGeometryResult"; }
    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<WaterHBondGeometryResult> Compute(
        ProteinConformation& conf, const SolventEnvironment& solvent);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    struct Candidate {
        size_t protein_atom = SIZE_MAX;
        size_t protein_residue = SIZE_MAX;
        size_t water_index = SIZE_MAX;
        int mode = 0;
        int protein_role = 0;
        double water_O_distance_A = 0.0;
        water_hbond_geometry_detail::GeometryEvaluation geometry;
        Vec3 water_O = Vec3::Zero();
        Vec3 water_H = Vec3::Zero();
    };

    std::vector<Candidate> candidates_;
    std::vector<std::array<std::int32_t, 6>> counts_;
    std::vector<std::array<double, 8>> nearest_;
};

}  // namespace nmr

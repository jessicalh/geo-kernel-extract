#pragma once
//
// DispersionResult: London dispersion (van der Waals) tensors from ring vertices.
//
// For each aromatic ring within range of each atom, computes the
// anisotropic dispersion kernel by summing over ring vertex atoms:
//
//   K_ab = sum_vertices S(r) * (3 d_a d_b / r^8 - delta_ab / r^6)
//   scalar = sum_vertices S(r) / r^6
//
// where d = r_atom - r_vertex, r = |d|, S(r) is a smooth switching
// function (CHARMM form, Brooks et al. 1983) that tapers to zero at
// DISP_VERTEX_R_CUT = 5.0 A.
//
// Properties of the tensor kernel per vertex:
//   - Symmetric: K_ab = K_ba
//   - Traceless: Tr(K) = S(r) * (3r^2/r^8 - 3/r^6) = 0
//   - Pure T2 per vertex; sum over vertices is also traceless
//
// The isotropic 1/r^6 scalar (attractive VdW) is stored separately.
// The tensor captures the anisotropic part — the angular dependence
// of dispersion from the discrete vertex positions around the ring.
//
// Filters:
//   - DipolarNearFieldFilter at ring level (source_extent = ring diameter):
//     same physics as all other ring calculators.
//   - Through-bond vertex exclusion: atoms bonded to any ring vertex are
//     excluded because the 1/r^6 kernel models through-space interaction,
//     not through-bond coupling. Uses protein bond connectivity, not
//     a distance heuristic.
//   - Smooth switching function: C^1 taper from 4.3A to 5.0A per vertex.
//     No hard cutoff discontinuity. Essential for MD frame continuity.
//
// Unit C6 = 1 — parameter is learnable per ring type.
// Units: tensor in Angstrom^-6, scalar in Angstrom^-6.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

class DispersionResult : public ConformationResult {
public:
    std::string Name() const override { return "DispersionResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute dispersion tensors for all atoms.
    static std::unique_ptr<DispersionResult> Compute(
        ProteinConformation& conf);

    // Grid sampling: evaluate dispersion kernel at arbitrary 3D point.
    SphericalTensor SampleShieldingAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr

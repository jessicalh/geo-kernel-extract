#pragma once
//
// DispersionResult: London dispersion (van der Waals) tensors from ring vertices.
//
// For each aromatic ring within range of each atom, computes the
// anisotropic dispersion kernel by summing over ring vertex atoms:
//
//   K_ab  = sum_vertices S(r) * (3 d_a d_b / r^8 - delta_ab / r^6)
//   scalar = sum_vertices S(r) / r^6
//
// where d = r_atom - r_vertex, r = |d|, and S(r) is a smooth switching
// function that tapers the kernel to zero at the vertex cutoff. The
// tensor captures the anisotropic part; the isotropic 1/r^6 scalar
// (attractive VdW) is stored separately. Implementation rationale —
// the kernel, its tracelessness, the switching-function derivation, and
// the ring/through-bond filter policy — lives in DispersionResult.cpp.
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
    SphericalTensor SampleKernelAt(Vec3 point) const;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr

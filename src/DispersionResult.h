#pragma once
//
// DispersionResult: unit-coefficient aromatic R^-6 proximity scalar
// from ring vertices.
//
// For each aromatic ring within range of each atom, computes the clean
// switched dispersion scalar by summing over ring vertex atoms:
//
//   scalar = sum_vertices S(r) / r^6
//
// where d = r_atom - r_vertex, r = |d|, and S(r) is a smooth switching
// function that tapers the kernel to zero at the vertex cutoff. The
// value is not D3/D4 dispersion energy and carries no fitted C6/C8 or
// damping coefficients. The former unit-C6 rank-2 tensor output is
// intentionally absent.
// Units: scalar in Angstrom^-6.
//

#include "ConformationResult.h"
#include "Ring.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

namespace dispersion_detail {

struct VertexResult {
    double scalar = 0.0;
    double distance = 0.0;
    bool valid = false;
};

VertexResult ComputeVertex(double distance);

// Exact production per-vertex evaluation for a finite ring.  In particular,
// it intentionally has no center-distance exclusion.
std::vector<VertexResult> ComputeRingVertices(
    const Vec3& atom_pos,
    const RingGeometry& geom);

}  // namespace dispersion_detail

class DispersionResult : public ConformationResult {
public:
    std::string Name() const override { return "DispersionResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute dispersion scalar/contact values for all atoms.
    static std::unique_ptr<DispersionResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr

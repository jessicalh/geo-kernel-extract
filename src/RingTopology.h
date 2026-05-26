#pragma once
// Typed ring storage derived from AtomSemanticTable ring-position
// slots. Pro pyrrolidine rings are stored separately from aromatic
// rings because the public Protein ring accessors expose aromatics.

#include "Types.h"
#include "Ring.h"
#include "SemanticEnums.h"
#include <vector>
#include <memory>

namespace nmr {

class Residue;
class CovalentTopology;

class RingTopology {
public:
    // Builds typed Ring objects by grouping residue atoms through
    // AtomSemanticTable RingPosition slots and applying the canonical
    // cyclic walk for each ring kind.
    //
    // Aromatic rings (PHE benzene, TYR phenol, HIS imidazole variants,
    // TRP benzene, TRP pyrrole, TRP indole perimeter) land in
    // aromatic_. Pro pyrrolidine lands in saturated_. The TRP fused-
    // partner index is set as a post-pass.
    //
    // FATAL+abort on:
    //   - Substrate ring slot whose label set is incomplete for the
    //     chemistry (e.g. Phe missing one of Ipso/Ortho/Meta/Para).
    //   - Substrate disambiguation failure (two BridgeFusion atoms
    //     where bond-graph adjacency cannot distinguish them).
    //
    // Returns an empty topology when atom_semantic is empty.
    //
    // Bonds disambiguate the two BridgeFusion atoms in the TRP-5 and
    // TRP-6 walks. Other walks use labels only, with Locant separating
    // HIS delta and epsilon heteroatoms.
    static std::unique_ptr<RingTopology> ConstructFromSubstrate(
        const std::vector<Residue>& residues,
        const std::vector<AtomSemanticTable>& atom_semantic,
        const CovalentTopology& bonds);

    size_t AromaticCount() const { return aromatic_.size(); }
    const Ring& AromaticAt(size_t i) const { return *aromatic_[i]; }
    const std::vector<std::unique_ptr<Ring>>& Aromatic() const {
        return aromatic_;
    }

    size_t SaturatedCount() const { return saturated_.size(); }
    const Ring& SaturatedAt(size_t i) const { return *saturated_[i]; }
    const std::vector<std::unique_ptr<Ring>>& Saturated() const {
        return saturated_;
    }

    // Public so LegacyAmberTopology can represent the empty-substrate
    // path with a non-null RingTopology.
    RingTopology() = default;

private:

    // The first three atoms of the walk determine the ring-normal sign
    // via Ring::ComputeGeometry's cross-product orientation fix
    // (Ring.cpp). Reversing a walk flips the normal and the ring-current
    // shielding sign at probe atoms.
    static std::vector<size_t> CanonicalCyclicWalk(
        RingSystemKind kind,
        const std::vector<size_t>& candidates,
        const std::vector<AtomSemanticTable>& atom_semantic,
        const CovalentTopology& bonds);

    std::vector<std::unique_ptr<Ring>> aromatic_;
    std::vector<std::unique_ptr<Ring>> saturated_;
};

}  // namespace nmr

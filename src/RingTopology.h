#pragma once
//
// RingTopology: typed aromatic + saturated ring storage owned by
// LegacyAmberTopology. The companion to CovalentTopology.
//
// Bundle C / Slice B (2026-05-07): rings move from Protein::rings_
// onto LegacyAmberTopology, parallel to how bonds live on
// CovalentTopology. The motivation is the legacy-amber implementation
// brief's authority map: LegacyAmberTopology owns the calculator
// interpretation; rings are calculator interpretation derived from
// the atom-semantic substrate. Putting ring storage on Protein was
// pre-existing technical debt that the substrate-driven construction
// path makes intolerable to compound. See
// `spec/plan/legacy-amber-implementation-brief-2026-04-29.md` and
// memory `feedback_object_model_scope_discipline`.
//
// Construction is substrate-driven. ConstructFromSubstrate reads
// per-atom RingPosition slots from the AtomSemanticTable substrate
// and produces typed Ring objects in canonical cyclic walk order.
// Pro pyrrolidine rings go to saturated_; everything else aromatic.
// FATAL+abort on substrate gaps; no string fallbacks (the substrate
// IS the source of truth post-Bundle-B).
//
// Calculator-facing API surface stays unchanged at the consumer:
// `protein.RingCount()` / `protein.RingAt(i)` (aromatic-only) keep
// meaning. The change is internal — Protein's accessors delegate
// through LegacyAmberTopology to RingTopology, parallel to how
// `protein.BondAt(i)` delegates through to CovalentTopology.
//

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
    // Substrate-driven construction. Walks each residue's atoms,
    // groups them by RingSystemKind via the AtomSemanticTable's
    // RingPosition slots, orders each ring's atoms via the explicit
    // CanonicalCyclicWalk dispatch, and constructs a typed Ring
    // object via the existing CreateRing factory.
    //
    // Aromatic rings (PHE benzene, TYR phenol, HIS imidazole variants,
    // TRP benzene, TRP pyrrole, TRP indole perimeter) land in
    // aromatic_. Pro pyrrolidine lands in saturated_. The TRP fused-
    // partner index is set as a post-pass.
    //
    // FATAL+abort on:
    //   - HIS residue with unresolved variant_index.
    //   - Substrate ring slot whose label set is incomplete for the
    //     chemistry (e.g. Phe missing one of Ipso/Ortho/Meta/Para).
    //   - Substrate disambiguation failure (two BridgeFusion atoms
    //     where bond-graph adjacency cannot distinguish them).
    //
    // Returns an empty topology when atom_semantic is empty (the
    // legitimate stub-fixture signal — calculator-physics tests with
    // raw element/position fixtures that never had PDB names). This
    // mirrors LegacyAmberTopology::HasAtomSemantic() behavior and
    // produces RingCount() == 0 at the Protein API.
    //
    // Bonds are needed for two specific disambiguations: TRP-5's two
    // BridgeFusion atoms (one bonded to NE1, one bonded to CG) and
    // TRP-6's two BridgeFusion atoms (one bonded to the Ortho1-labelled
    // atom, one bonded to the Ortho2-labelled atom). PHE / TYR /
    // HIS / Pro / TRP-9 walks are label-only (with Locant
    // disambiguating HIS's two Heteroatom_* atoms).
    static std::unique_ptr<RingTopology> ConstructFromSubstrate(
        const std::vector<Residue>& residues,
        const std::vector<AtomSemanticTable>& atom_semantic,
        const CovalentTopology& bonds);

    // ── Aromatic rings (PHE/TYR/HIS-variants/TRP-{benzene,pyrrole,9}) ──
    size_t AromaticCount() const { return aromatic_.size(); }
    const Ring& AromaticAt(size_t i) const { return *aromatic_[i]; }
    const std::vector<std::unique_ptr<Ring>>& Aromatic() const {
        return aromatic_;
    }

    // ── Saturated rings (Pro pyrrolidine; future: other saturated) ──
    size_t SaturatedCount() const { return saturated_.size(); }
    const Ring& SaturatedAt(size_t i) const { return *saturated_[i]; }
    const std::vector<std::unique_ptr<Ring>>& Saturated() const {
        return saturated_;
    }

    // Default constructor produces an empty topology. Public so
    // LegacyAmberTopology can synthesise an empty RingTopology in
    // its stub-fixture path (the constructor accepts null `rings`
    // and replaces it with `make_unique<RingTopology>()`). The
    // intended construction path is the static factory
    // ConstructFromSubstrate; an empty default-constructed
    // RingTopology is the legitimate stub-fixture signal.
    RingTopology() = default;

private:

    // CanonicalCyclicWalk: produce atom_indices for one ring instance
    // in the canonical cyclic walk order.
    //
    // The first three atoms of the walk determine the ring-normal sign
    // via Ring::ComputeGeometry's cross-product orientation fix
    // (Ring.cpp:32-37). The full walk determines serialized vertex
    // order in fileformat H5 emission. Bless bit-identity on the eight
    // standard ring types is the gate.
    //
    // Per-RingSystemKind dispatch:
    //   Benzene_Phe / Benzene_Tyr  -> WalkSixRingByLabel
    //   Imidazole_His              -> WalkHisImidazole
    //   Indole_Trp_5               -> WalkIndolePyrrole
    //   Indole_Trp_6               -> WalkIndoleBenzene
    //   Indole_Trp_9               -> WalkIndolePerimeter
    //   Pyrrolidine_Pro            -> WalkProPyrrolidine
    //
    // Where the substrate labels alone are symmetric, disambiguation
    // is explicit:
    //   - HIS two Heteroatom_* atoms: by Locant (Delta < Epsilon).
    //   - TRP-5 two BridgeFusion atoms: by bond-graph adjacency (the
    //     "first" is bonded to the Heteroatom_NH-labelled NE1; the
    //     "second" is bonded to the Ipso-labelled CG).
    //   - TRP-6 two BridgeFusion atoms: by bond-graph adjacency (the
    //     "first" is bonded to the Ortho1-labelled CE3; the "second"
    //     is bonded to the Ortho2-labelled CZ2).
    //
    // Each per-kind walker has a doc block citing today's
    // AminoAcidType.cpp atom-name order it reproduces and the
    // disambiguation rule used.
    static std::vector<size_t> CanonicalCyclicWalk(
        RingSystemKind kind,
        const std::vector<size_t>& candidates,
        const std::vector<AtomSemanticTable>& atom_semantic,
        const CovalentTopology& bonds);

    std::vector<std::unique_ptr<Ring>> aromatic_;
    std::vector<std::unique_ptr<Ring>> saturated_;
};

}  // namespace nmr

// MolecularFrameSelect -- the typed selector that fills a MolFrameSpec for an
// atom by dispatching on the viewer's ALREADY-TYPED chemistry (BackboneRole,
// PlanarGroupKind, ring membership), instead of AnalysisAtom.cpp's
// configureMolecularFrame IUPAC-name-string matching. Same frames, but chosen
// off the enriched model the loader already built -- protein-agnostic, no
// string fragility, and it reuses enrichment the viewer already paid for.
//
// This is the thin adapter layer that DOES couple to QtProtein/QtResidue; the
// frame GEOMETRY it feeds (MolecularFrame.h) stays pure and shared-able. It
// resolves anchor ATOM INDICES only -- per-frame axes come from
// BuildMolecularFrameAxes with a position lookup.
//
// Returns std::nullopt for atoms with no defined molecular frame; those fall
// back to the frame-free PAS / invariant view (CsaShape.h).

#pragma once

#include "MolecularFrame.h"
#include "QtProtein.h"  // QtProtein, QtAtom, QtResidue, QtRing, topology
#include "QtRing.h"

#include <cstddef>
#include <cstdint>
#include <optional>

namespace h5reader::model {

// First ring (by index) whose canonical atom list contains `atomIndex`, or
// nullopt. Linear over rings -- a protein has a handful, called once per frame
// selection, not per trajectory step.
inline std::optional<std::size_t> RingContainingAtom(const QtProtein& protein,
                                                     std::int32_t atomIndex) {
    for (std::size_t r = 0; r < protein.ringCount(); ++r) {
        const QtRing& ring = protein.ring(r);
        for (std::int32_t a : ring.atomIndices)
            if (a == atomIndex) return r;
    }
    return std::nullopt;
}

inline std::optional<MolFrameSpec> SelectMolecularFrameSpec(const QtProtein& protein,
                                                            std::size_t atomIndex) {
    if (atomIndex >= protein.atomCount()) return std::nullopt;
    const QtAtom& atom = protein.atom(atomIndex);
    if (atom.residueIndex < 0
        || static_cast<std::size_t>(atom.residueIndex) >= protein.residueCount())
        return std::nullopt;
    const QtResidue& res = protein.residue(static_cast<std::size_t>(atom.residueIndex));

    // ---- Backbone carbonyl C / O -> peptide-plane frame ----
    // origin C, x toward O, plane toward the next residue's amide N (the peptide
    // plane) falling back to this residue's CA at the C-terminus.
    if (atom.backboneRole == BackboneRole::CarbonylCarbon
        || atom.backboneRole == BackboneRole::CarbonylOxygen) {
        if (res.HasC() && res.HasO() && res.HasCA()) {
            MolFrameSpec spec;
            spec.kind = MolecularFrameKind::BackboneCarbonyl;
            spec.origin = res.C;
            spec.xAnchor = res.O;
            std::int32_t plane = res.CA;
            if (res.nextResidueIndex >= 0
                && static_cast<std::size_t>(res.nextResidueIndex) < protein.residueCount()) {
                const QtResidue& nextRes =
                    protein.residue(static_cast<std::size_t>(res.nextResidueIndex));
                if (nextRes.HasN()) plane = nextRes.N;
            }
            spec.planeAnchor = plane;
            return spec;
        }
    }

    // ---- Backbone amide N -> origin N, x toward amide H (else CA), plane
    // toward the previous residue's carbonyl C (else CA at the N-terminus). ----
    if (atom.backboneRole == BackboneRole::Nitrogen) {
        if (res.HasN() && res.HasCA()) {
            MolFrameSpec spec;
            spec.kind = MolecularFrameKind::BackboneAmideN;
            spec.origin = res.N;
            spec.xAnchor = (res.H != QtResidue::NONE) ? res.H : res.CA;
            std::int32_t plane = res.CA;
            if (res.prevResidueIndex >= 0
                && static_cast<std::size_t>(res.prevResidueIndex) < protein.residueCount()) {
                const QtResidue& prevRes =
                    protein.residue(static_cast<std::size_t>(res.prevResidueIndex));
                if (prevRes.HasC()) plane = prevRes.C;
            }
            spec.planeAnchor = plane;
            return spec;
        }
    }

    // ---- Aromatic-ring atom -> ring-local frame: x radial to the heavy atom
    // (an aromatic H borrows its heavy parent), z along the ring normal. ----
    if (atom.aromatic && atom.IsInAnyRing()) {
        const std::int32_t heavy =
            (atom.element == Element::H && atom.parentAtomIndex >= 0)
                ? atom.parentAtomIndex
                : static_cast<std::int32_t>(atomIndex);
        if (const auto ringIdx = RingContainingAtom(protein, heavy)) {
            MolFrameSpec spec;
            spec.kind = MolecularFrameKind::AromaticRingLocal;
            spec.ring = static_cast<std::int32_t>(*ringIdx);
            spec.heavy = heavy;
            return spec;
        }
    }

    // TODO(next): sidechain groups via PlanarGroupKind == Carboxylate /
    // Guanidinium / SidechainAmide (gather the group's C/O/N members from
    // res.atomIndices by element), and Met SD (element S + locant). Same typed
    // pattern; same BuildMolecularFrameAxes kinds (already implemented).
    return std::nullopt;
}

}  // namespace h5reader::model

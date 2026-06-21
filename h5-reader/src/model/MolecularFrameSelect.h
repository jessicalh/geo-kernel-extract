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

    // ---- Aromatic-ring atom, OR a hydrogen attached to one -> ring-local
    // frame (z = ring normal, x radial to the heavy atom). An aromatic ring H
    // is NOT itself flagged aromatic / in-ring, so resolve to its heavy parent
    // and test THAT -- the stats pass frames ring H's too, and we match it
    // (otherwise ~56 aromatic H's silently drop out of the framed set). ----
    {
        const std::int32_t heavy =
            (atom.element == Element::H && atom.parentAtomIndex >= 0)
                ? atom.parentAtomIndex
                : static_cast<std::int32_t>(atomIndex);
        if (heavy >= 0 && static_cast<std::size_t>(heavy) < protein.atomCount()) {
            const QtAtom& heavyAtom = protein.atom(static_cast<std::size_t>(heavy));
            if (heavyAtom.aromatic && heavyAtom.IsInAnyRing()) {
                if (const auto ringIdx = RingContainingAtom(protein, heavy)) {
                    MolFrameSpec spec;
                    spec.kind = MolecularFrameKind::AromaticRingLocal;
                    spec.ring = static_cast<std::int32_t>(*ringIdx);
                    spec.heavy = heavy;
                    return spec;
                }
            }
        }
    }

    // ---- Sidechain planar groups: gather the group's atoms from the residue
    // by shared PlanarGroupKind + element (the planar moiety tags all its
    // members), then anchor the typed frame. ----

    if (atom.planarGroup == PlanarGroupKind::Carboxylate) {
        // C + two O (ASP CG/OD1/OD2, GLU CD/OE1/OE2, C-terminus).
        std::int32_t c = -1, o1 = -1, o2 = -1;
        for (std::int32_t ai : res.atomIndices) {
            const QtAtom& a = protein.atom(static_cast<std::size_t>(ai));
            if (a.planarGroup != PlanarGroupKind::Carboxylate) continue;
            if (a.element == Element::C) { if (c < 0) c = ai; }
            else if (a.element == Element::O) { if (o1 < 0) o1 = ai; else if (o2 < 0) o2 = ai; }
        }
        if (c >= 0 && o1 >= 0 && o2 >= 0) {
            MolFrameSpec spec;
            spec.kind = MolecularFrameKind::SidechainCarboxylate;
            spec.origin = c;
            spec.xAnchor = o1;
            spec.secondAnchor = o2;
            return spec;
        }
    }

    if (atom.planarGroup == PlanarGroupKind::Guanidinium) {
        // CZ + the guanidinium N (ARG): origin CZ, x toward one N, plane another.
        std::int32_t cz = -1, n1 = -1, n2 = -1;
        for (std::int32_t ai : res.atomIndices) {
            const QtAtom& a = protein.atom(static_cast<std::size_t>(ai));
            if (a.planarGroup != PlanarGroupKind::Guanidinium) continue;
            if (a.element == Element::C) { if (cz < 0) cz = ai; }
            else if (a.element == Element::N) { if (n1 < 0) n1 = ai; else if (n2 < 0) n2 = ai; }
        }
        if (cz >= 0 && n1 >= 0 && n2 >= 0) {
            MolFrameSpec spec;
            spec.kind = MolecularFrameKind::SidechainGuanidinium;
            spec.origin = cz;
            spec.xAnchor = n1;
            spec.planeAnchor = n2;
            return spec;
        }
    }

    if (atom.planarGroup == PlanarGroupKind::SidechainAmide) {
        // C + O + N (ASN CG/OD1/ND2, GLN CD/OE1/NE2).
        std::int32_t c = -1, o = -1, n = -1;
        for (std::int32_t ai : res.atomIndices) {
            const QtAtom& a = protein.atom(static_cast<std::size_t>(ai));
            if (a.planarGroup != PlanarGroupKind::SidechainAmide) continue;
            if (a.element == Element::C) { if (c < 0) c = ai; }
            else if (a.element == Element::O) { if (o < 0) o = ai; }
            else if (a.element == Element::N) { if (n < 0) n = ai; }
        }
        if (c >= 0 && o >= 0 && n >= 0) {
            MolFrameSpec spec;
            spec.kind = MolecularFrameKind::SidechainCarboxamide;
            spec.origin = c;
            spec.xAnchor = o;
            spec.planeAnchor = n;
            return spec;
        }
    }

    // Methionine thioether sulphur (SD): origin SD, x toward CE, plane toward
    // CG. Cys SG has no epsilon carbon, so requiring CG + CE excludes it.
    if (atom.element == Element::S) {
        std::int32_t cg = -1, ce = -1;
        for (std::int32_t ai : res.atomIndices) {
            const QtAtom& a = protein.atom(static_cast<std::size_t>(ai));
            if (a.element != Element::C) continue;
            if (a.locant == Locant::Gamma) cg = ai;
            else if (a.locant == Locant::Epsilon) ce = ai;
        }
        if (cg >= 0 && ce >= 0) {
            MolFrameSpec spec;
            spec.kind = MolecularFrameKind::MetSd;
            spec.origin = static_cast<std::int32_t>(atomIndex);
            spec.xAnchor = ce;
            spec.planeAnchor = cg;
            return spec;
        }
    }

    return std::nullopt;
}

}  // namespace h5reader::model

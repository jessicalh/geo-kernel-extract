#include "CameraAnchorHelper.h"

#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"

#include <cstddef>

namespace h5reader::app {

namespace {

// Map a residue's typed backbone-atom-index cache (built at load from
// typed BackboneRole + Locant; no string scan — feedback_identity_from_
// chemistry_not_position) into the typed std::size_t shape the
// CameraMode constructors expect. NONE entries return nullopt.
std::optional<std::size_t> ToOptIndex(int32_t v) {
    if (v == model::QtResidue::NONE) return std::nullopt;
    if (v < 0) return std::nullopt;
    return static_cast<std::size_t>(v);
}

}  // namespace

FocusAnchorResult DeriveFocusAnchor(const model::QtProtein& protein,
                                     std::size_t focusAtomIdx,
                                     FocusAnchorKind kind) {
    FocusAnchorResult out;

    if (focusAtomIdx >= protein.atomCount()) {
        out.outcome = FocusAnchorOutcome::AtomIndexOutOfRange;
        return out;
    }
    const auto& atom = protein.atom(focusAtomIdx);
    if (atom.residueIndex < 0
        || static_cast<std::size_t>(atom.residueIndex) >= protein.residueCount()) {
        out.outcome = FocusAnchorOutcome::AtomHasNoResidue;
        return out;
    }
    const auto& residue =
        protein.residue(static_cast<std::size_t>(atom.residueIndex));

    const auto n  = ToOptIndex(residue.N);
    const auto ca = ToOptIndex(residue.CA);
    const auto c  = ToOptIndex(residue.C);
    if (!n || !ca || !c) {
        out.outcome = FocusAnchorOutcome::MissingBackboneAtoms;
        return out;
    }

    if (kind == FocusAnchorKind::Plane) {
        // 3-atom plane lock on backbone (N, CA, C). Default policy →
        // PerpendicularToPlane (sight along plane normal). The focal
        // point is the centroid of the three backbone atoms — the focus
        // atom sits within ~1-2 Å of that centroid for backbone foci and
        // ~3-5 Å for sidechain foci on this residue's local frame.
        out.mode = PlaneMode(*n, *ca, *c);
        out.policy = DefaultPolicy();
        out.outcome = FocusAnchorOutcome::Ok;
        return out;
    }

    // Dihedral kinds need the flanking residue's atoms for the full
    // torsion. Both phi and psi straddle the current residue's backbone.
    if (kind == FocusAnchorKind::DihedralPhi) {
        // Phi = (prev_C, N, CA, C). Sight axis = (N, CA).
        if (residue.prevResidueIndex < 0
            || static_cast<std::size_t>(residue.prevResidueIndex)
                   >= protein.residueCount()) {
            out.outcome = FocusAnchorOutcome::MissingDihedralNeighbor;
            return out;
        }
        const auto& prev = protein.residue(
            static_cast<std::size_t>(residue.prevResidueIndex));
        const auto prevC = ToOptIndex(prev.C);
        if (!prevC) {
            out.outcome = FocusAnchorOutcome::MissingDihedralNeighbor;
            return out;
        }
        out.mode = DihedralMode(*prevC, *n, *ca, *c);
        out.policy = DownAxisPolicy(*n, *ca);
        out.outcome = FocusAnchorOutcome::Ok;
        return out;
    }

    // DihedralPsi = (N, CA, C, next_N). Sight axis = (CA, C).
    if (residue.nextResidueIndex < 0
        || static_cast<std::size_t>(residue.nextResidueIndex)
               >= protein.residueCount()) {
        out.outcome = FocusAnchorOutcome::MissingDihedralNeighbor;
        return out;
    }
    const auto& next = protein.residue(
        static_cast<std::size_t>(residue.nextResidueIndex));
    const auto nextN = ToOptIndex(next.N);
    if (!nextN) {
        out.outcome = FocusAnchorOutcome::MissingDihedralNeighbor;
        return out;
    }
    out.mode = DihedralMode(*n, *ca, *c, *nextN);
    out.policy = DownAxisPolicy(*ca, *c);
    out.outcome = FocusAnchorOutcome::Ok;
    return out;
}

}  // namespace h5reader::app

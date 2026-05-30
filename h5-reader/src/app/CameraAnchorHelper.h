// CameraAnchorHelper — derive a typed CameraMode from a focus atom.
//
// The user-facing scenario: "I want to look at THIS atom and see its
// immediate neighborhood coherently as time advances." A 3-atom plane
// lock on the atom's residue backbone (N / CA / C) holds the local
// chemistry steady; the focus atom (which usually IS one of N / CA / C
// or sits near them) appears at the plane's centroid each frame, so the
// 1.5 px-floor of plane lock applies to the user's actual focus.
//
// Per spec/viewport_pipeline_2026-05-30.md §3.1: PlaneMode lifts to the
// composer's CameraMode::Plane variant. The math is unchanged; this
// helper is a typed-identity-driven shortcut that reaches into
// QtProtein/QtResidue for the backbone atom indices instead of asking
// the user to type three numeric atom IDs.
//
// Dihedral variant: derive the phi or psi torsion atoms for the focus
// residue. Phi = (prev_C, N, CA, C); psi = (N, CA, C, next_N). The
// caller picks which torsion via the Kind argument. Terminal residues
// missing the required neighbors return AnchorOutcome::MissingAtoms;
// the REST handler turns that into 422 Unprocessable Entity.
//
// Pure-function shape — no QObject, no slots, no I/O. Tested from
// scene_math_tests against synthetic typed-residue fixtures.

#pragma once

#include "CameraMode.h"
#include "OrientationPolicy.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {
class QtProtein;
struct QtResidue;
}  // namespace h5reader::model

namespace h5reader::app {

// Outcome of a focus-atom anchor derivation. The CameraMode + policy
// are valid only on Outcome::Ok; the others carry the failure reason
// so the REST handler can pick the right HTTP status.
enum class FocusAnchorOutcome {
    Ok,
    AtomIndexOutOfRange,    // atom >= protein.atomCount()
    AtomHasNoResidue,        // atom.residueIndex == -1
    MissingBackboneAtoms,    // residue lacks N / CA / C
    MissingDihedralNeighbor, // dihedral kind, terminal residue missing prev/next
};

// What kind of anchor to derive from the focus atom's residue.
enum class FocusAnchorKind {
    // 3-atom plane lock on the residue's backbone (N, CA, C). Recommended
    // for the "focus atom + local neighborhood coherent" use case: sub-
    // pixel floor on the plane lock applies to the residue around the
    // focus atom.
    Plane,
    // 4-atom dihedral lock on the residue's phi torsion (prev_C, N, CA, C).
    // Sights down the (N, CA) axis (DownAxisPolicy). Used by reveal
    // bindings that want a Newman-projection view of a specific torsion.
    DihedralPhi,
    // 4-atom dihedral lock on the residue's psi torsion (N, CA, C, next_N).
    // Sights down the (CA, C) axis. Symmetric counterpart to DihedralPhi.
    DihedralPsi,
};

struct FocusAnchorResult {
    FocusAnchorOutcome outcome = FocusAnchorOutcome::Ok;
    CameraMode         mode;
    OrientationPolicy  policy;
};

// Derive a CameraMode + OrientationPolicy from a focus atom and a
// requested anchor kind. Returns an Outcome describing success or the
// specific failure shape so the caller can format an appropriate error.
FocusAnchorResult DeriveFocusAnchor(const model::QtProtein& protein,
                                     std::size_t focusAtomIdx,
                                     FocusAnchorKind kind);

}  // namespace h5reader::app

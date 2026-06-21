// NewmanProjection — resolve a residue's backbone phi/psi torsion and project
// the substituent bonds into the 2D plane perpendicular to the sight axis, for
// a live Newman diagram. Pure geometry over (QtProtein, Conformation): no Qt
// widget / VTK dependency beyond QString labels, so it is unit-checkable and
// REST-verifiable on its own. Torsions are rigid-transform invariant, so this
// runs on the base Conformation (the displayed TransformedConformation gives
// the same angles).
//
//   phi_i = C(i-1) - N(i) - CA(i) - C(i)    sight down N -> CA   (front N, back CA)
//   psi_i = N(i)   - CA(i) - C(i)  - N(i+1) sight down CA -> C   (front CA, back C)
//
// Termini / missing neighbour C(i-1) or N(i+1) -> valid == false with a reason.

#pragma once

#include <QString>

#include <cstddef>
#include <vector>

namespace h5reader::model {
class QtProtein;
class Conformation;
}

namespace h5reader::app {

enum class NewmanKind { Phi, Psi };

// One substituent bond drawn in the projection disc. angleRad is measured in
// the plane perpendicular to the sight axis: the front reference bond sits at
// 0 and the back reference spoke at the signed torsion angle. The panel maps
// front spokes from the hub and back spokes from the rim.
struct NewmanSpoke {
    double  angleRad  = 0.0;
    QString label;              // short atom label, e.g. "C(-1)", "H", "CB", "O"
    bool    front     = false;  // substituent of the front (near) atom vs back (far)
    bool    reference = false;  // one of the two torsion-defining bonds
};

struct NewmanProjection {
    bool        valid = false;
    QString     invalidReason;       // terminus / missing atom / degenerate
    NewmanKind  kind  = NewmanKind::Phi;
    double      torsionDeg = 0.0;    // (-180, 180], Blondel-Karplus sign
    std::size_t residueIndex = 0;
    QString     residueLabel;        // 3-letter code projection
    QString     frontLabel;          // near atom: "N" (phi) / "CA" (psi)
    QString     backLabel;           // far atom:  "CA" (phi) / "C"  (psi)
    std::vector<NewmanSpoke> spokes;
};

NewmanProjection ComputeNewmanProjection(const model::QtProtein& protein,
                                         const model::Conformation& conf,
                                         std::size_t residueIndex, std::size_t frame,
                                         NewmanKind kind);

}  // namespace h5reader::app

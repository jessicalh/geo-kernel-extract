// QtAtomNames — projection layer for atom-name strings, separate from
// the typed-identity QtAtom struct.
//
// The three atom-name strings (AMBER ff14SB, IUPAC, BMRB) in
// atoms_category_info.npy are GENUINE projection data: AMBER is the
// canonical force-field name needed for the inspector display, IUPAC
// is the systematic name for cross-referencing publications, BMRB is
// the experimental-shift-database name for ML matching against
// measured chemical shifts.
//
// They are NOT chemistry inputs. Per the no-strings discipline
// (notes/H5_READER_REWRITE_DESIGN_2026-05-23.md §2), they live HERE,
// not on QtAtom. Consumers that want a label call
// `QtProtein::atomNames(i)` — an explicit projection access. Code
// that asks chemistry questions on QtAtom stays on typed enums.
//
// The two NamingProvenance fields record whether the IUPAC and BMRB
// names came from atom_nom.tbl (Match) or fell back to the AMBER name
// (MissLogged) — data-quality flag for the ML-matching consumer.
//
// Storage strategy (design §11.C — eager arrays): QtProtein holds a
// std::vector<QtAtomNames> parallel to its atoms array. Memory cost
// is trivial (~25 KB for a 1000-atom protein); call-site clarity
// wins over lazy materialisation.

#pragma once

#include "Types.h"  // NamingProvenance

#include <QString>

namespace h5reader::model {

struct QtAtomNames {
    QString amber;  // AMBER ff14SB canonical name (S8 from sidecar)
    QString iupac;  // atom_nom.tbl IUPAC projection
    QString bmrb;   // atom_nom.tbl BMRB projection

    NamingProvenance iupacProvenance = NamingProvenance::Match;
    NamingProvenance bmrbProvenance = NamingProvenance::Match;
};

}  // namespace h5reader::model

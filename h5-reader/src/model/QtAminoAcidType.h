// QtAminoAcidType — static registry of amino-acid chemistry, mirroring
// the library's `nmr::AminoAcidType` + `GetAminoAcidType()` shape.
//
// Reader-side single authority for the predicates the substrate
// doesn't carry directly per residue (is_aromatic, is_titratable,
// has_amide_H, chi_angle_count) plus the canonical 3-letter and
// 1-letter naming projections derived from the AminoAcid enum.
//
// Storing these as a static lookup table (one row per AminoAcid)
// rather than as fields on QtResidue is the no-strings discipline:
// the 3-letter code is a `const char*` literal once, in this table,
// not a QString stored per residue.
//
// Library mirror: src/AminoAcidType.h (the full library version
// includes variant tables, atom lists, ring lists, chi atom maps —
// the reader doesn't need any of that because the substrate is
// already populated per atom from the sidecar; the helper subset
// here covers only the per-residue predicates).

#pragma once

#include "Types.h"  // AminoAcid

namespace h5reader::model {

struct QtAminoAcidType {
    AminoAcid type = AminoAcid::Unknown;
    const char* three_letter_code = "UNK";
    char one_letter_code = 'X';
    bool is_aromatic = false;
    bool is_titratable = false;
    bool has_amide_H = true;  // false only for PRO
    int chi_angle_count = 0;
};

// Returns the static row for an AminoAcid. Unknown returns the UNK
// sentinel. Reference stable for the program lifetime.
const QtAminoAcidType& GetQtAminoAcidType(AminoAcid aa);

}  // namespace h5reader::model

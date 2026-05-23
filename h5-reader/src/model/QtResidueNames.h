// QtResidueNames — derivation helpers for residue-name display
// surfaces, replacing the previously-stored string fields.
//
// The library's `Residue` has no name-string field; it carries the
// typed `AminoAcid` enum and computes the 3-letter code on demand via
// `GetAminoAcidType(type).three_letter_code`. The reader does the
// same here.
//
// Three free functions cover the projection surface:
//   - AmberResidue3LetterFor(aa, variantIdx) — variant-aware (HIS →
//     HID/HIE/HIP, ASP → ASH, CYS → CYX/CYM)
//   - IupacResidue3LetterFor(aa) — canonical 3-letter (variant-blind)
//   - ResidueOneLetterFor(aa) — single-letter code
//
// They return `const char*` literals from the static
// QtAminoAcidType table or hardcoded variant-name tables. Zero string
// allocation per call. The display-side caller (inspector dock,
// tooltip) decides whether to wrap in QString.

#pragma once

#include "Types.h"  // AminoAcid

namespace h5reader::model {

// Variant-aware AMBER 3-letter code. For most residues, returns the
// canonical 3-letter. For HIS/ASP/GLU/LYS/CYS/TYR with non-default
// variantIdx, returns the variant name (HID, HIE, HIP, ASH, GLH,
// LYN, CYX, CYM, TYM).
//
// variantIdx semantics: -1 = default protonation (returns the
// canonical 3-letter); 0..N = variant index per the library's
// AminoAcidType::variants table for this AA.
const char* AmberResidue3LetterFor(AminoAcid aa, int variantIdx);

// IUPAC canonical 3-letter code, variant-blind. Same as
// `GetQtAminoAcidType(aa).three_letter_code` — convenience name.
const char* IupacResidue3LetterFor(AminoAcid aa);

// One-letter code from QtAminoAcidType registry.
char ResidueOneLetterFor(AminoAcid aa);

}  // namespace h5reader::model

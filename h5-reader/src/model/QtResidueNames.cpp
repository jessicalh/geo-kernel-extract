// QtResidueNames — derivation helpers for residue display names.
//
// AmberResidue3LetterFor() handles variant-aware naming for the six
// titratable residues (HIS/ASP/GLU/LYS/CYS/TYR). The variant index
// semantics follow the library's AminoAcidType::variants ordering;
// values outside the defined range fall back to the canonical IUPAC
// 3-letter code.
//
// The variant tables here are minimal — only the named variants
// commonly present in AMBER ff14SB protein systems. Extensions
// (e.g. ASH-with-sp2-O variants) would extend the variant_idx range
// and add cases here in lockstep with the library's
// AminoAcidType::variants[type].

#include "QtResidueNames.h"

#include "QtAminoAcidType.h"

namespace h5reader::model {

const char* AmberResidue3LetterFor(AminoAcid aa, int variantIdx) {
    // -1 == default protonation: return the IUPAC canonical 3-letter.
    if (variantIdx < 0)
        return IupacResidue3LetterFor(aa);

    switch (aa) {
    case AminoAcid::HIS:
        // ff14SB variants: 0=HID (delta), 1=HIE (epsilon), 2=HIP (doubly protonated)
        switch (variantIdx) {
        case 0:
            return "HID";
        case 1:
            return "HIE";
        case 2:
            return "HIP";
        default:
            return "HIS";
        }
    case AminoAcid::ASP:
        // 0=ASH (protonated carboxyl)
        return (variantIdx == 0) ? "ASH" : "ASP";
    case AminoAcid::GLU:
        // 0=GLH (protonated carboxyl)
        return (variantIdx == 0) ? "GLH" : "GLU";
    case AminoAcid::LYS:
        // 0=LYN (neutral)
        return (variantIdx == 0) ? "LYN" : "LYS";
    case AminoAcid::CYS:
        // 0=CYX (disulfide), 1=CYM (thiolate)
        switch (variantIdx) {
        case 0:
            return "CYX";
        case 1:
            return "CYM";
        default:
            return "CYS";
        }
    case AminoAcid::TYR:
        // 0=TYM (phenolate)
        return (variantIdx == 0) ? "TYM" : "TYR";
    default:
        // Non-titratable residues have no variant table.
        return IupacResidue3LetterFor(aa);
    }
}

const char* IupacResidue3LetterFor(AminoAcid aa) {
    return GetQtAminoAcidType(aa).three_letter_code;
}

char ResidueOneLetterFor(AminoAcid aa) {
    return GetQtAminoAcidType(aa).one_letter_code;
}

}  // namespace h5reader::model

// QtAminoAcidType registry — static table mirroring the library's
// per-AA chemistry predicates. Values match nmr/src/AminoAcidType.cpp
// for the standard 20.
//
// Chi-angle counts per Markley et al. 1998 J. Biomol. NMR §2.1
// (sidechain dihedral conventions): ALA/GLY 0, SER/THR/CYS/VAL 1,
// PHE/TYR/TRP/HIS/ASN/ASP/LEU/ILE/PRO 2, MET/GLU/GLN 3, LYS/ARG 4.

#include "QtAminoAcidType.h"

#include <array>

namespace h5reader::model {

namespace {

// Indexed by AminoAcid ordinal. UNK sentinel at index 20.
// The const char* / char fields are literal addresses (program-
// lifetime stable); the bool fields are constants.
constexpr std::array<QtAminoAcidType, 21> kTable = {{
    /* ALA = 0  */ {AminoAcid::ALA, "ALA", 'A', false, false, true, 0},
    /* ARG = 1  */ {AminoAcid::ARG, "ARG", 'R', false, true, true, 4},
    /* ASN = 2  */ {AminoAcid::ASN, "ASN", 'N', false, false, true, 2},
    /* ASP = 3  */ {AminoAcid::ASP, "ASP", 'D', false, true, true, 2},
    /* CYS = 4  */ {AminoAcid::CYS, "CYS", 'C', false, true, true, 1},
    /* GLN = 5  */ {AminoAcid::GLN, "GLN", 'Q', false, false, true, 3},
    /* GLU = 6  */ {AminoAcid::GLU, "GLU", 'E', false, true, true, 3},
    /* GLY = 7  */ {AminoAcid::GLY, "GLY", 'G', false, false, true, 0},
    /* HIS = 8  */ {AminoAcid::HIS, "HIS", 'H', true, true, true, 2},
    /* ILE = 9  */ {AminoAcid::ILE, "ILE", 'I', false, false, true, 2},
    /* LEU = 10 */ {AminoAcid::LEU, "LEU", 'L', false, false, true, 2},
    /* LYS = 11 */ {AminoAcid::LYS, "LYS", 'K', false, true, true, 4},
    /* MET = 12 */ {AminoAcid::MET, "MET", 'M', false, false, true, 3},
    /* PHE = 13 */ {AminoAcid::PHE, "PHE", 'F', true, false, true, 2},
    /* PRO = 14 */ {AminoAcid::PRO, "PRO", 'P', false, false, false, 2},
    /* SER = 15 */ {AminoAcid::SER, "SER", 'S', false, false, true, 1},
    /* THR = 16 */ {AminoAcid::THR, "THR", 'T', false, false, true, 1},
    /* TRP = 17 */ {AminoAcid::TRP, "TRP", 'W', true, false, true, 2},
    /* TYR = 18 */ {AminoAcid::TYR, "TYR", 'Y', true, true, true, 2},
    /* VAL = 19 */ {AminoAcid::VAL, "VAL", 'V', false, false, true, 1},
    /* UNK = 20 */ {AminoAcid::Unknown, "UNK", 'X', false, false, false, 0},
}};

constexpr QtAminoAcidType kUnknownSentinel = {AminoAcid::Unknown, "UNK", 'X', false, false, false, 0};

}  // namespace

const QtAminoAcidType& GetQtAminoAcidType(AminoAcid aa) {
    const int idx = static_cast<int>(aa);
    if (idx < 0 || idx > 20)
        return kUnknownSentinel;
    return kTable[static_cast<size_t>(idx)];
}

}  // namespace h5reader::model

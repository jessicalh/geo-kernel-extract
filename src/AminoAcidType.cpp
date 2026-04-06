// AminoAcidType: the complete amino acid chemistry table.
// Single source of truth for all amino acid properties.

#include "AminoAcidType.h"

namespace nmr {

using A = AminoAcidAtom;
using E = Element;

// Standard backbone atoms (most amino acids)
#define BB  A{"N",E::N,true}, A{"CA",E::C,true}, A{"C",E::C,true}, A{"O",E::O,true}, A{"H",E::H,true}, A{"HA",E::H,true}
#define BB_PRO  A{"N",E::N,true}, A{"CA",E::C,true}, A{"C",E::C,true}, A{"O",E::O,true}, A{"HA",E::H,true}
#define BB_GLY  A{"N",E::N,true}, A{"CA",E::C,true}, A{"C",E::C,true}, A{"O",E::O,true}, A{"H",E::H,true}, A{"HA2",E::H,true}, A{"HA3",E::H,true}

static const std::vector<AminoAcidType> AMINO_ACID_TYPES = {

// ALA
{ AminoAcid::ALA, "ALA", 'A', false, false, true, 0, 0,
  {BB, {"CB",E::C,false},{"HB1",E::H,false},{"HB2",E::H,false},{"HB3",E::H,false}},
  {}, {}, {} },

// ARG (pKa ~12.5, almost always +1, but PROPKA predicts it)
{ AminoAcid::ARG, "ARG", 'R', false, true, true, 4, +1,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},{"HG2",E::H,false},{"HG3",E::H,false},
   {"CD",E::C,false},{"HD2",E::H,false},{"HD3",E::H,false},
   {"NE",E::N,false},{"HE",E::H,false},
   {"CZ",E::C,false},
   {"NH1",E::N,false},{"HH11",E::H,false},{"HH12",E::H,false},
   {"NH2",E::N,false},{"HH21",E::H,false},{"HH22",E::H,false}},
  {},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}, {"CB","CG","CD","NE"}, {"CG","CD","NE","CZ"}},
  {{"ARN", "deprotonated arginine", 0, "deprotonated"}} },

// ASN
{ AminoAcid::ASN, "ASN", 'N', false, false, true, 2, 0,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},{"OD1",E::O,false},
   {"ND2",E::N,false},{"HD21",E::H,false},{"HD22",E::H,false}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","OD1"}}, {} },

// ASP
{ AminoAcid::ASP, "ASP", 'D', false, true, true, 2, -1,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},{"OD1",E::O,false},{"OD2",E::O,false}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","OD1"}},
  {{"ASH", "protonated aspartate", 0, "protonated"}} },

// CYS
{ AminoAcid::CYS, "CYS", 'C', false, true, true, 1, -1,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"SG",E::S,false},{"HG",E::H,false}},
  {}, {{"N","CA","CB","SG"}},
  {{"CYX", "disulfide bonded", 0, "disulfide"}, {"CYM", "deprotonated thiolate", -1, "deprotonated"}} },

// GLN
{ AminoAcid::GLN, "GLN", 'Q', false, false, true, 3, 0,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},{"HG2",E::H,false},{"HG3",E::H,false},
   {"CD",E::C,false},{"OE1",E::O,false},
   {"NE2",E::N,false},{"HE21",E::H,false},{"HE22",E::H,false}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}, {"CB","CG","CD","OE1"}}, {} },

// GLU
{ AminoAcid::GLU, "GLU", 'E', false, true, true, 3, -1,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},{"HG2",E::H,false},{"HG3",E::H,false},
   {"CD",E::C,false},{"OE1",E::O,false},{"OE2",E::O,false}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}, {"CB","CG","CD","OE1"}},
  {{"GLH", "protonated glutamate", 0, "protonated"}} },

// GLY
{ AminoAcid::GLY, "GLY", 'G', false, false, true, 0, 0,
  {BB_GLY}, {}, {}, {} },

// HIS
{ AminoAcid::HIS, "HIS", 'H', true, true, true, 2, +1,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},
   {"ND1",E::N,false},
   {"CD2",E::C,false},{"HD2",E::H,false},
   {"CE1",E::C,false},{"HE1",E::H,false},
   {"NE2",E::N,false},{"HE2",E::H,false}},
  {{RingTypeIndex::HisImidazole, {"CG","ND1","CE1","NE2","CD2"}}},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","ND1"}},
  {{"HID", "Nd-protonated (delta)", 0, "delta"},
   {"HIE", "Ne-protonated (epsilon)", 0, "epsilon"},
   {"HIP", "doubly protonated", +1, "doubly"}} },

// ILE
{ AminoAcid::ILE, "ILE", 'I', false, false, true, 2, 0,
  {BB,
   {"CB",E::C,false},{"HB",E::H,false},
   {"CG1",E::C,false},{"HG12",E::H,false},{"HG13",E::H,false},
   {"CG2",E::C,false},{"HG21",E::H,false},{"HG22",E::H,false},{"HG23",E::H,false},
   {"CD1",E::C,false},{"HD11",E::H,false},{"HD12",E::H,false},{"HD13",E::H,false}},
  {}, {{"N","CA","CB","CG1"}, {"CA","CB","CG1","CD1"}}, {} },

// LEU
{ AminoAcid::LEU, "LEU", 'L', false, false, true, 2, 0,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},{"HG",E::H,false},
   {"CD1",E::C,false},{"HD11",E::H,false},{"HD12",E::H,false},{"HD13",E::H,false},
   {"CD2",E::C,false},{"HD21",E::H,false},{"HD22",E::H,false},{"HD23",E::H,false}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","CD1"}}, {} },

// LYS
{ AminoAcid::LYS, "LYS", 'K', false, true, true, 4, +1,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},{"HG2",E::H,false},{"HG3",E::H,false},
   {"CD",E::C,false},{"HD2",E::H,false},{"HD3",E::H,false},
   {"CE",E::C,false},{"HE2",E::H,false},{"HE3",E::H,false},
   {"NZ",E::N,false},{"HZ1",E::H,false},{"HZ2",E::H,false},{"HZ3",E::H,false}},
  {},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}, {"CB","CG","CD","CE"}, {"CG","CD","CE","NZ"}},
  {{"LYN", "deprotonated lysine", 0, "deprotonated"}} },

// MET
{ AminoAcid::MET, "MET", 'M', false, false, true, 3, 0,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},{"HG2",E::H,false},{"HG3",E::H,false},
   {"SD",E::S,false},
   {"CE",E::C,false},{"HE1",E::H,false},{"HE2",E::H,false},{"HE3",E::H,false}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","SD"}, {"CB","CG","SD","CE"}}, {} },

// PHE
{ AminoAcid::PHE, "PHE", 'F', true, false, true, 2, 0,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},
   {"CD1",E::C,false},{"HD1",E::H,false},
   {"CD2",E::C,false},{"HD2",E::H,false},
   {"CE1",E::C,false},{"HE1",E::H,false},
   {"CE2",E::C,false},{"HE2",E::H,false},
   {"CZ",E::C,false},{"HZ",E::H,false}},
  {{RingTypeIndex::PheBenzene, {"CG","CD1","CE1","CZ","CE2","CD2"}}},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD1"}}, {} },

// PRO
{ AminoAcid::PRO, "PRO", 'P', false, false, false, 2, 0,
  {BB_PRO,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},{"HG2",E::H,false},{"HG3",E::H,false},
   {"CD",E::C,false},{"HD2",E::H,false},{"HD3",E::H,false}},
  {}, {{"N","CA","CB","CG"}, {"CA","CB","CG","CD"}}, {} },

// SER
{ AminoAcid::SER, "SER", 'S', false, false, true, 1, 0,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"OG",E::O,false},{"HG",E::H,false}},
  {}, {{"N","CA","CB","OG"}}, {} },

// THR
{ AminoAcid::THR, "THR", 'T', false, false, true, 1, 0,
  {BB,
   {"CB",E::C,false},{"HB",E::H,false},
   {"OG1",E::O,false},{"HG1",E::H,false},
   {"CG2",E::C,false},{"HG21",E::H,false},{"HG22",E::H,false},{"HG23",E::H,false}},
  {}, {{"N","CA","CB","OG1"}}, {} },

// TRP
{ AminoAcid::TRP, "TRP", 'W', true, false, true, 2, 0,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},
   {"CD1",E::C,false},{"HD1",E::H,false},
   {"CD2",E::C,false},
   {"NE1",E::N,false},{"HE1",E::H,false},
   {"CE2",E::C,false},
   {"CE3",E::C,false},{"HE3",E::H,false},
   {"CZ2",E::C,false},{"HZ2",E::H,false},
   {"CZ3",E::C,false},{"HZ3",E::H,false},
   {"CH2",E::C,false},{"HH2",E::H,false}},
  {{RingTypeIndex::TrpBenzene,   {"CD2","CE2","CZ2","CH2","CZ3","CE3"}},
   {RingTypeIndex::TrpPyrrole,   {"CG","CD1","NE1","CE2","CD2"}},
   {RingTypeIndex::TrpPerimeter, {"CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2"}}},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD1"}}, {} },

// TYR
{ AminoAcid::TYR, "TYR", 'Y', true, true, true, 2, -1,
  {BB,
   {"CB",E::C,false},{"HB2",E::H,false},{"HB3",E::H,false},
   {"CG",E::C,false},
   {"CD1",E::C,false},{"HD1",E::H,false},
   {"CD2",E::C,false},{"HD2",E::H,false},
   {"CE1",E::C,false},{"HE1",E::H,false},
   {"CE2",E::C,false},{"HE2",E::H,false},
   {"CZ",E::C,false},{"OH",E::O,false},{"HH",E::H,false}},
  {{RingTypeIndex::TyrPhenol, {"CG","CD1","CE1","CZ","CE2","CD2"}}},
  {{"N","CA","CB","CG"}, {"CA","CB","CG","CD1"}},
  {{"TYM", "deprotonated tyrosinate", -1, "deprotonated"}} },

// VAL
{ AminoAcid::VAL, "VAL", 'V', false, false, true, 1, 0,
  {BB,
   {"CB",E::C,false},{"HB",E::H,false},
   {"CG1",E::C,false},{"HG11",E::H,false},{"HG12",E::H,false},{"HG13",E::H,false},
   {"CG2",E::C,false},{"HG21",E::H,false},{"HG22",E::H,false},{"HG23",E::H,false}},
  {}, {{"N","CA","CB","CG1"}}, {} },

// Unknown
{ AminoAcid::Unknown, "UNK", 'X', false, false, false, 0, 0,
  {}, {}, {}, {} },

};

#undef BB
#undef BB_PRO
#undef BB_GLY

const AminoAcidType& GetAminoAcidType(AminoAcid aa) {
    int i = static_cast<int>(aa);
    if (i < 0 || i >= static_cast<int>(AMINO_ACID_TYPES.size()))
        return AMINO_ACID_TYPES.back();
    return AMINO_ACID_TYPES[static_cast<size_t>(i)];
}

const std::vector<AminoAcidType>& AllAminoAcidTypes() {
    return AMINO_ACID_TYPES;
}

const AminoAcidType& AminoAcidTypeFromCode(const std::string& code) {
    for (const auto& type : AMINO_ACID_TYPES)
        if (code == type.three_letter_code) return type;
    return AMINO_ACID_TYPES.back();
}


void ValidateVariantIndices() {
    // Verify the variant ordering contract. A swap here silently breaks
    // every protonation assignment in the system.
    auto check = [](AminoAcid aa, int idx, const char* expected_name) {
        const auto& aat = GetAminoAcidType(aa);
        if (idx >= static_cast<int>(aat.variants.size()) ||
            std::string(aat.variants[idx].name) != expected_name) {
            std::fprintf(stderr,
                "FATAL: AminoAcidType variant index contract violated: "
                "%s variants[%d] expected \"%s\", got \"%s\"\n",
                aat.three_letter_code, idx, expected_name,
                (idx < static_cast<int>(aat.variants.size()))
                    ? aat.variants[idx].name : "(out of range)");
            std::abort();
        }
    };

    check(AminoAcid::HIS, 0, "HID");
    check(AminoAcid::HIS, 1, "HIE");
    check(AminoAcid::HIS, 2, "HIP");
    check(AminoAcid::ASP, 0, "ASH");
    check(AminoAcid::GLU, 0, "GLH");
    check(AminoAcid::CYS, 0, "CYX");
    check(AminoAcid::CYS, 1, "CYM");
    check(AminoAcid::LYS, 0, "LYN");
    check(AminoAcid::ARG, 0, "ARN");
    check(AminoAcid::TYR, 0, "TYM");
}

}  // namespace nmr

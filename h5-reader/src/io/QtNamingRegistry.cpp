#include "QtNamingRegistry.h"

namespace h5reader::io {

using model::AminoAcid;
using model::ProtonationVariant;

namespace {

// 20 standard codes in AminoAcid enum order.
constexpr const char* kStandardCodes[20] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL"
};

}  // namespace

AminoAcid ThreeLetterCodeToAminoAcid(const std::string& code) {
    for (int i = 0; i < 20; ++i) {
        if (code == kStandardCodes[i]) return static_cast<AminoAcid>(i);
    }

    // AMBER protonation variants for HIS
    if (code == "HID" || code == "HIE" || code == "HIP") return AminoAcid::HIS;
    // CHARMM protonation variants for HIS
    if (code == "HSD" || code == "HSE" || code == "HSP") return AminoAcid::HIS;

    if (code == "ASH")                 return AminoAcid::ASP;
    if (code == "GLH")                 return AminoAcid::GLU;
    if (code == "CYX" || code == "CYM") return AminoAcid::CYS;
    if (code == "LYN")                 return AminoAcid::LYS;
    if (code == "ARN")                 return AminoAcid::ARG;
    if (code == "TYM")                 return AminoAcid::TYR;

    // Selenomethionine → MET (standard substitution).
    if (code == "MSE")                 return AminoAcid::MET;

    return AminoAcid::Unknown;
}

const char* ThreeLetterCodeForAminoAcid(AminoAcid aa) {
    const int i = static_cast<int>(aa);
    if (i >= 0 && i < 20) return kStandardCodes[i];
    return "UNK";
}

const char* ProtonationVariantName(ProtonationVariant v) {
    switch (v) {
        case ProtonationVariant::Default:          return "default";
        case ProtonationVariant::HIS_delta:        return "HID";
        case ProtonationVariant::HIS_epsilon:      return "HIE";
        case ProtonationVariant::HIS_doubly:       return "HIP";
        case ProtonationVariant::ASP_protonated:   return "ASH";
        case ProtonationVariant::GLU_protonated:   return "GLH";
        case ProtonationVariant::LYS_deprotonated: return "LYN";
        case ProtonationVariant::ARG_deprotonated: return "ARN";
        case ProtonationVariant::TYR_deprotonated: return "TYM";
        case ProtonationVariant::CYS_disulfide:    return "CYX";
        case ProtonationVariant::CYS_deprotonated: return "CYM";
    }
    return "?";
}

QString CharmmAtomNameToCanonical(const QString& charmmName,
                                   const QString& /*residueName*/) {
    // Backbone amide H — the most visible CHARMM/ff14SB difference.
    if (charmmName == QLatin1String("HN")) return QStringLiteral("H");

    // CHARMM TYR hydroxyl hydrogen is "HH" in CHARMM, stays "HH" in
    // IUPAC/PDB. Terminal hydrogens (H1/H2/H3 on N-terminus, OXT O-
    // terminus) are already canonical in both.
    //
    // Extend as needed when inspector display surfaces an unexpected
    // name. Do not over-extend speculatively.
    return charmmName;
}

}  // namespace h5reader::io

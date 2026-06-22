#pragma once

#include "../model/Types.h"

#include <array>
#include <cstddef>
#include <string_view>

namespace h5reader::rediscover {

enum class LiteratureStatus { Cited, GoodEnough, Placeholder };

struct LiteratureConstant {
    const char* key;
    double value;
    const char* units;
    LiteratureStatus status;
    const char* source;
};

struct LiteratureStatusCounts {
    int cited = 0;
    int good_enough = 0;
    int placeholder = 0;
};

inline constexpr const char* LiteratureStatusName(LiteratureStatus s) {
    switch (s) {
    case LiteratureStatus::Cited: return "cited";
    case LiteratureStatus::GoodEnough: return "good_enough";
    case LiteratureStatus::Placeholder: return "placeholder";
    }
    return "placeholder";
}

inline constexpr LiteratureConstant kBuckinghamA_H{
    "buckingham.A.H", 0.20, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "Buckingham linear field coefficient; good-enough proton magnitude until per-IUPAC citation is acquired"};
inline constexpr LiteratureConstant kBuckinghamA_C{
    "buckingham.A.C", 0.06, "ppm/(V/angstrom)", LiteratureStatus::Placeholder,
    "self-placeholder heavy-atom linear field coefficient"};
inline constexpr LiteratureConstant kBuckinghamA_N{
    "buckingham.A.N", 0.08, "ppm/(V/angstrom)", LiteratureStatus::Placeholder,
    "self-placeholder heavy-atom linear field coefficient"};
inline constexpr LiteratureConstant kBuckinghamA_O{
    "buckingham.A.O", 0.05, "ppm/(V/angstrom)", LiteratureStatus::Placeholder,
    "self-placeholder heavy-atom linear field coefficient"};
inline constexpr LiteratureConstant kBuckinghamA_S{
    "buckingham.A.S", 0.10, "ppm/(V/angstrom)", LiteratureStatus::Placeholder,
    "self-placeholder heavy-atom linear field coefficient"};
inline constexpr LiteratureConstant kBuckinghamA_Generic{
    "buckingham.A.generic", 0.05, "ppm/(V/angstrom)", LiteratureStatus::Placeholder,
    "self-placeholder generic linear field coefficient"};

inline constexpr LiteratureConstant kBuckinghamB_Generic{
    "buckingham.B.generic", 0.010, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "good-enough quadratic term, about 20x smaller than the proton linear response at 1 V/angstrom"};

inline constexpr LiteratureConstant kCoulombKeVA{
    "electrostatics.coulomb_ke", 14.3996, "V*angstrom/e", LiteratureStatus::GoodEnough,
    "standard molecular electrostatic conversion used by charge-field and EFG producers"};

inline constexpr LiteratureConstant kSigma0_H{
    "sigma0.H", 30.0, "ppm", LiteratureStatus::Placeholder,
    "self-placeholder zero-field absolute shielding reference"};
inline constexpr LiteratureConstant kSigma0_C{
    "sigma0.C", 170.0, "ppm", LiteratureStatus::Placeholder,
    "self-placeholder zero-field absolute shielding reference"};
inline constexpr LiteratureConstant kSigma0_N{
    "sigma0.N", 250.0, "ppm", LiteratureStatus::Placeholder,
    "self-placeholder zero-field absolute shielding reference"};
inline constexpr LiteratureConstant kSigma0_O{
    "sigma0.O", 350.0, "ppm", LiteratureStatus::Placeholder,
    "self-placeholder zero-field absolute shielding reference"};
inline constexpr LiteratureConstant kSigma0_S{
    "sigma0.S", 700.0, "ppm", LiteratureStatus::Placeholder,
    "self-placeholder zero-field absolute shielding reference"};
inline constexpr LiteratureConstant kSigma0_Generic{
    "sigma0.generic", 0.0, "ppm", LiteratureStatus::Placeholder,
    "self-placeholder zero-field absolute shielding reference"};

inline constexpr LiteratureConstant kMcConnellPeptideCO{
    "mcconnell.delta_chi.peptide_co", 2.41, "10^-6 cm^3/mol", LiteratureStatus::Cited,
    "Hooper-Kaiser/Case lineage already carried on disk"};
inline constexpr LiteratureConstant kMcConnellPeptideCN{
    "mcconnell.delta_chi.peptide_cn", -5.42, "10^-6 cm^3/mol", LiteratureStatus::Cited,
    "Hooper-Kaiser/Case lineage already carried on disk"};
inline constexpr LiteratureConstant kMcConnellSidechainCO{
    "mcconnell.delta_chi.sidechain_co", 2.41, "10^-6 cm^3/mol", LiteratureStatus::Cited,
    "Peptide C=O value reused as good literature-side default for sidechain C=O"};
inline constexpr LiteratureConstant kMcConnellAromaticZero{
    "mcconnell.delta_chi.aromatic_zeroed", 0.0, "10^-6 cm^3/mol", LiteratureStatus::GoodEnough,
    "aromatic pi current is carried by ring-current kernels to avoid double-counting"};
inline constexpr LiteratureConstant kMcConnellMolarPrefactor{
    "mcconnell.molar_prefactor", 1.0e24 / 6.02214076e23, "angstrom^-3 per 10^-6 cm^3/mol",
    LiteratureStatus::Cited, "SI exact Avogadro constant conversion"};

inline constexpr LiteratureConstant kLarsenWaterTerm{
    "larsen.water_term", 2.07, "ppm", LiteratureStatus::Cited,
    "Larsen 2015 ProCS15 NMA-water amide-H term; DOI 10.7717/peerj.1344"};
inline constexpr LiteratureConstant kLarsenShieldingTensors{
    "larsen.hbond_shielding_tensors", 1.0, "producer ppm tensor scale", LiteratureStatus::Cited,
    "Larsen 2015 ProCS15 table/grid producer emits ppm tensors"};

inline constexpr std::array<LiteratureConstant, 9> kRingIntensityByType{{
    {"ring.intensity.PheBenzene", -12.0, "nA/T", LiteratureStatus::Cited, "Giessner-Prettre aromatic ring current"},
    {"ring.intensity.TyrPhenol", -11.28, "nA/T", LiteratureStatus::Cited, "Giessner-Prettre aromatic ring current"},
    {"ring.intensity.TrpBenzene", -12.48, "nA/T", LiteratureStatus::Cited, "Giessner-Prettre aromatic ring current"},
    {"ring.intensity.TrpPyrrole", -6.72, "nA/T", LiteratureStatus::Cited, "Giessner-Prettre aromatic ring current"},
    {"ring.intensity.TrpPerimeter", -19.2, "nA/T", LiteratureStatus::GoodEnough, "indole perimeter compound current"},
    {"ring.intensity.HisImidazole", -5.16, "nA/T", LiteratureStatus::Cited, "Giessner-Prettre aromatic ring current"},
    {"ring.intensity.HidImidazole", -5.16, "nA/T", LiteratureStatus::GoodEnough, "histidine tautomer reuse of imidazole current"},
    {"ring.intensity.HieImidazole", -5.16, "nA/T", LiteratureStatus::GoodEnough, "histidine tautomer reuse of imidazole current"},
    {"ring.intensity.ProPyrrolidine", 0.0, "nA/T", LiteratureStatus::Cited, "saturated ring, no aromatic pi current"},
}};

inline constexpr std::array<LiteratureConstant, 9> kJohnsonBoveyLobeOffsetByType{{
    {"ring.jb_lobe_offset.PheBenzene", 0.64, "angstrom", LiteratureStatus::Cited, "Johnson-Bovey two-loop offset"},
    {"ring.jb_lobe_offset.TyrPhenol", 0.64, "angstrom", LiteratureStatus::Cited, "Johnson-Bovey two-loop offset"},
    {"ring.jb_lobe_offset.TrpBenzene", 0.64, "angstrom", LiteratureStatus::Cited, "Johnson-Bovey two-loop offset"},
    {"ring.jb_lobe_offset.TrpPyrrole", 0.52, "angstrom", LiteratureStatus::Cited, "Johnson-Bovey two-loop offset"},
    {"ring.jb_lobe_offset.TrpPerimeter", 0.60, "angstrom", LiteratureStatus::GoodEnough, "compound indole perimeter offset"},
    {"ring.jb_lobe_offset.HisImidazole", 0.50, "angstrom", LiteratureStatus::Cited, "Johnson-Bovey two-loop offset"},
    {"ring.jb_lobe_offset.HidImidazole", 0.50, "angstrom", LiteratureStatus::GoodEnough, "histidine tautomer reuse of imidazole offset"},
    {"ring.jb_lobe_offset.HieImidazole", 0.50, "angstrom", LiteratureStatus::GoodEnough, "histidine tautomer reuse of imidazole offset"},
    {"ring.jb_lobe_offset.ProPyrrolidine", 0.0, "angstrom", LiteratureStatus::Cited, "saturated ring, no aromatic pi current"},
}};

inline constexpr std::array<LiteratureConstant, 39> kAllLiteratureConstants{{
    kBuckinghamA_H,
    kBuckinghamA_C,
    kBuckinghamA_N,
    kBuckinghamA_O,
    kBuckinghamA_S,
    kBuckinghamA_Generic,
    kBuckinghamB_Generic,
    kCoulombKeVA,
    kSigma0_H,
    kSigma0_C,
    kSigma0_N,
    kSigma0_O,
    kSigma0_S,
    kSigma0_Generic,
    kMcConnellPeptideCO,
    kMcConnellPeptideCN,
    kMcConnellSidechainCO,
    kMcConnellAromaticZero,
    kMcConnellMolarPrefactor,
    kLarsenWaterTerm,
    kLarsenShieldingTensors,
    kRingIntensityByType[0],
    kRingIntensityByType[1],
    kRingIntensityByType[2],
    kRingIntensityByType[3],
    kRingIntensityByType[4],
    kRingIntensityByType[5],
    kRingIntensityByType[6],
    kRingIntensityByType[7],
    kRingIntensityByType[8],
    kJohnsonBoveyLobeOffsetByType[0],
    kJohnsonBoveyLobeOffsetByType[1],
    kJohnsonBoveyLobeOffsetByType[2],
    kJohnsonBoveyLobeOffsetByType[3],
    kJohnsonBoveyLobeOffsetByType[4],
    kJohnsonBoveyLobeOffsetByType[5],
    kJohnsonBoveyLobeOffsetByType[6],
    kJohnsonBoveyLobeOffsetByType[7],
    kJohnsonBoveyLobeOffsetByType[8],
}};

inline constexpr LiteratureStatusCounts CountLiteratureConstantStatuses() {
    LiteratureStatusCounts counts;
    for (const LiteratureConstant& c : kAllLiteratureConstants) {
        switch (c.status) {
        case LiteratureStatus::Cited: ++counts.cited; break;
        case LiteratureStatus::GoodEnough: ++counts.good_enough; break;
        case LiteratureStatus::Placeholder: ++counts.placeholder; break;
        }
    }
    return counts;
}

inline LiteratureConstant BuckinghamA(model::Element element,
                                      std::string_view /*residueType*/,
                                      std::string_view /*atomName*/) {
    switch (element) {
    case model::Element::H: return kBuckinghamA_H;
    case model::Element::C: return kBuckinghamA_C;
    case model::Element::N: return kBuckinghamA_N;
    case model::Element::O: return kBuckinghamA_O;
    case model::Element::S: return kBuckinghamA_S;
    default: return kBuckinghamA_Generic;
    }
}

inline LiteratureConstant BuckinghamB(model::Element /*element*/,
                                      std::string_view /*residueType*/,
                                      std::string_view /*atomName*/) {
    return kBuckinghamB_Generic;
}

inline double CoulombKeVA() {
    return kCoulombKeVA.value;
}

inline LiteratureConstant Sigma0(model::Element element,
                                 std::string_view /*residueType*/,
                                 std::string_view /*atomName*/) {
    switch (element) {
    case model::Element::H: return kSigma0_H;
    case model::Element::C: return kSigma0_C;
    case model::Element::N: return kSigma0_N;
    case model::Element::O: return kSigma0_O;
    case model::Element::S: return kSigma0_S;
    default: return kSigma0_Generic;
    }
}

inline LiteratureConstant McConnellDeltaChi(model::BondCategory category) {
    switch (category) {
    case model::BondCategory::PeptideCO: return kMcConnellPeptideCO;
    case model::BondCategory::PeptideCN: return kMcConnellPeptideCN;
    case model::BondCategory::SidechainCO: return kMcConnellSidechainCO;
    case model::BondCategory::Aromatic: return kMcConnellAromaticZero;
    default:
        return {"mcconnell.delta_chi.inapplicable", 0.0, "10^-6 cm^3/mol",
                LiteratureStatus::Placeholder, "self-placeholder for unsupported bond category"};
    }
}

inline LiteratureConstant RingIntensity(model::RingTypeIndex type) {
    const int idx = static_cast<int>(type);
    if (idx >= 0 && idx < static_cast<int>(kRingIntensityByType.size()))
        return kRingIntensityByType[static_cast<std::size_t>(idx)];
    return {"ring.intensity.unknown", 0.0, "nA/T", LiteratureStatus::Placeholder,
            "self-placeholder for unknown ring type"};
}

inline LiteratureConstant JohnsonBoveyLobeOffset(model::RingTypeIndex type) {
    const int idx = static_cast<int>(type);
    if (idx >= 0 && idx < static_cast<int>(kJohnsonBoveyLobeOffsetByType.size()))
        return kJohnsonBoveyLobeOffsetByType[static_cast<std::size_t>(idx)];
    return {"ring.jb_lobe_offset.unknown", 0.0, "angstrom", LiteratureStatus::Placeholder,
            "self-placeholder for unknown ring type"};
}

}  // namespace h5reader::rediscover

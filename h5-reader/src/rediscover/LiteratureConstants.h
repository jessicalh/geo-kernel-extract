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
    const char* axis_convention_note = "";
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
    "buckingham.A.H.untyped_proxy", 1.50, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "Grayson & Chittenden 2000 Int. J. Mol. Sci. 1, 49-60, DOI 10.3390/ijms1030049; "
    "generic proton proxy when no amide-H class is available",
    "sigma=sigma0-A*Eparallel; generic H proxy, positive axis is the engine-selected Eparallel axis"};
inline constexpr LiteratureConstant kBuckinghamA_C{
    "buckingham.A.C.untyped_proxy", -0.96, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "Grayson & Chittenden 2000 Int. J. Mol. Sci. 1, 49-60, DOI 10.3390/ijms1030049; "
    "aliphatic sp3 13C proxy for untyped carbon",
    "sigma=sigma0-A*Eparallel; sp3 C proxy, positive axis is the engine-selected Eparallel axis"};
inline constexpr LiteratureConstant kBuckinghamA_N{
    "buckingham.A.N.untyped_proxy", 16.88, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "Jensen et al. 2015 JACS 137, 5573-5582, DOI 10.1021/ja512957g; "
    "amide 15N proxy for untyped protein nitrogen",
    "sigma=sigma0-A*Eparallel; amide-N proxy, positive axis is the engine-selected Eparallel axis"};
inline constexpr LiteratureConstant kBuckinghamA_O{
    "buckingham.A.O.untyped_proxy", -6.55, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "Rizzo et al. 1995 J. Chem. Phys. 102, 8953-8966, DOI 10.1063/1.468797; "
    "H2O/hydroxyl 17O proxy for untyped oxygen",
    "sigma=sigma0-A*Eparallel; H2O C2-axis proxy, positive axis is the engine-selected Eparallel axis"};
inline constexpr LiteratureConstant kBuckinghamA_S{
    "buckingham.A.S.structural_absent", 0.0, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "structural-absent marker: no sourced class-specific sulfur shielding-polarizability constant wired",
    "sigma=sigma0-A*Eparallel; sulfur Buckingham response intentionally zeroed until sourced"};
inline constexpr LiteratureConstant kBuckinghamA_Generic{
    "buckingham.A.generic.structural_absent", 0.0, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "structural-absent marker: no sourced class-specific generic shielding-polarizability constant wired",
    "sigma=sigma0-A*Eparallel; unknown element Buckingham response intentionally zeroed until sourced"};

inline constexpr LiteratureConstant kBuckinghamA_Amide15N{
    "buckingham.A.amide_15N", 16.88, "ppm/(V/angstrom)", LiteratureStatus::Cited,
    "Jensen et al. 2015 JACS 137, 5573-5582, DOI 10.1021/ja512957g; "
    "PMC4480918 supplement reports A_parallel=868 ppm/a.u.",
    "sigma=sigma0-A*Eparallel; amide 15N A_parallel, positive-axis agreement with engine frame must be checked"};
inline constexpr LiteratureConstant kBuckinghamA_Amide1HN{
    "buckingham.A.amide_1HN", 3.66, "ppm/(V/angstrom)", LiteratureStatus::Cited,
    "Jensen et al. 2015 JACS 137, 5573-5582, DOI 10.1021/ja512957g; "
    "PMC4480918 supplement reports A_parallel=188 ppm/a.u.",
    "sigma=sigma0-A*Eparallel; amide 1HN along N->H bond field axis"};
inline constexpr LiteratureConstant kBuckinghamA_Carbonyl13C{
    "buckingham.A.carbonyl_13C", 7.28, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "Augspurger, Dykstra & Oldfield 1991 JACS 113, 2447-2451, DOI 10.1021/ja00007a015; "
    "CO d_sigma_iso/dE sign-flipped to project convention",
    "sigma=sigma0-A*Eparallel; carbonyl sign depends on CO/source positive axis versus engine Eparallel"};
inline constexpr LiteratureConstant kBuckinghamA_AliphaticSp3_13C{
    "buckingham.A.aliphatic_sp3_13C", -0.96, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "Grayson & Chittenden 2000 Int. J. Mol. Sci. 1, 49-60, DOI 10.3390/ijms1030049; "
    "C2H6 sp3 13C A=-49.2 ppm/a.u.",
    "sigma=sigma0-A*Eparallel; tetrahedral C principal-axis proxy, engine Eparallel axis is local frame dependent"};
inline constexpr LiteratureConstant kBuckinghamA_Carbonyl17O{
    "buckingham.A.carbonyl_17O", -29.69, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "Augspurger, Dykstra & Oldfield 1991 JACS 113, 2447-2451, DOI 10.1021/ja00007a015; "
    "CO d_sigma_iso/dE sign-flipped to project convention",
    "sigma=sigma0-A*Eparallel; carbonyl sign depends on CO/source positive axis versus engine Eparallel"};
inline constexpr LiteratureConstant kBuckinghamA_Hydroxyl17O{
    "buckingham.A.hydroxyl_17O", -6.55, "ppm/(V/angstrom)", LiteratureStatus::GoodEnough,
    "Rizzo et al. 1995 J. Chem. Phys. 102, 8953-8966, DOI 10.1063/1.468797; "
    "H2O 17O A_z=-336.7 ppm/a.u.",
    "sigma=sigma0-A*Eparallel; H2O C2-axis hydroxyl/water proxy"};

inline constexpr LiteratureConstant kBuckinghamB_H{
    "buckingham.B.H.untyped_proxy", 0.030, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "Grayson & Chittenden 2000 Int. J. Mol. Sci. 1, 49-60, DOI 10.3390/ijms1030049; "
    "generic proton quadratic proxy",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; generic H proxy"};
inline constexpr LiteratureConstant kBuckinghamB_C{
    "buckingham.B.C.untyped_proxy", 0.35, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "Grayson & Chittenden 2000 Int. J. Mol. Sci. 1, 49-60, DOI 10.3390/ijms1030049; "
    "aliphatic sp3 13C proxy for untyped carbon",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; sp3 C proxy"};
inline constexpr LiteratureConstant kBuckinghamB_N{
    "buckingham.B.N.untyped_proxy", 0.0, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "Jensen et al. 2015 JACS 137, 5573-5582, DOI 10.1021/ja512957g; "
    "amide 15N B unreported and neglected in CSP literature",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; B=0 is an explicit unreported/neglected marker"};
inline constexpr LiteratureConstant kBuckinghamB_O{
    "buckingham.B.O.untyped_proxy", 0.52, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "Rizzo et al. 1995 J. Chem. Phys. 102, 8953-8966, DOI 10.1063/1.468797; "
    "H2O/hydroxyl 17O B_zz=1367.6 ppm/a.u.^2",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; H2O C2-axis proxy"};
inline constexpr LiteratureConstant kBuckinghamB_S{
    "buckingham.B.S.structural_absent", 0.0, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "structural-absent marker: no sourced class-specific sulfur quadratic Buckingham constant wired",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; sulfur Buckingham response intentionally zeroed until sourced"};

inline constexpr LiteratureConstant kBuckinghamB_Generic{
    "buckingham.B.generic.structural_absent", 0.0, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "structural-absent marker: no sourced class-specific generic quadratic Buckingham constant wired",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; unknown element Buckingham response intentionally zeroed"};

inline constexpr LiteratureConstant kBuckinghamB_Amide15N{
    "buckingham.B.amide_15N", 0.0, "ppm/(V/angstrom)^2", LiteratureStatus::Cited,
    "Jensen et al. 2015 JACS 137, 5573-5582, DOI 10.1021/ja512957g; "
    "15N B not reported and neglected in CSP model",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; B=0 is an explicit unreported/neglected marker"};
inline constexpr LiteratureConstant kBuckinghamB_Amide1HN{
    "buckingham.B.amide_1HN", 0.0, "ppm/(V/angstrom)^2", LiteratureStatus::Cited,
    "Jensen et al. 2015 JACS 137, 5573-5582, DOI 10.1021/ja512957g; "
    "1HN B not reported and neglected in CSP model",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; B=0 is an explicit unreported/neglected marker"};
inline constexpr LiteratureConstant kBuckinghamB_Carbonyl13C{
    "buckingham.B.carbonyl_13C", 0.0, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "Augspurger, Dykstra & Oldfield 1991 JACS 113, 2447-2451, DOI 10.1021/ja00007a015; "
    "CO carbonyl-like 13C B unreported",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; B=0 is an explicit unreported marker"};
inline constexpr LiteratureConstant kBuckinghamB_AliphaticSp3_13C{
    "buckingham.B.aliphatic_sp3_13C", 0.35, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "Grayson & Chittenden 2000 Int. J. Mol. Sci. 1, 49-60, DOI 10.3390/ijms1030049; "
    "working methyl/sp3 13C B proxy",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; tetrahedral C scalar B proxy"};
inline constexpr LiteratureConstant kBuckinghamB_Carbonyl17O{
    "buckingham.B.carbonyl_17O", 0.0, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "Augspurger, Dykstra & Oldfield 1991 JACS 113, 2447-2451, DOI 10.1021/ja00007a015; "
    "CO carbonyl-like 17O B unreported",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; B=0 is an explicit unreported marker"};
inline constexpr LiteratureConstant kBuckinghamB_Hydroxyl17O{
    "buckingham.B.hydroxyl_17O", 0.52, "ppm/(V/angstrom)^2", LiteratureStatus::GoodEnough,
    "Rizzo et al. 1995 J. Chem. Phys. 102, 8953-8966, DOI 10.1063/1.468797; "
    "H2O 17O B_zz=1367.6 ppm/a.u.^2",
    "sigma=sigma0-A*Eparallel-B*Eparallel^2; H2O C2-axis hydroxyl/water proxy"};

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

inline constexpr std::array<LiteratureConstant, 56> kAllLiteratureConstants{{
    kBuckinghamA_H,
    kBuckinghamA_C,
    kBuckinghamA_N,
    kBuckinghamA_O,
    kBuckinghamA_S,
    kBuckinghamA_Generic,
    kBuckinghamA_Amide15N,
    kBuckinghamA_Amide1HN,
    kBuckinghamA_Carbonyl13C,
    kBuckinghamA_AliphaticSp3_13C,
    kBuckinghamA_Carbonyl17O,
    kBuckinghamA_Hydroxyl17O,
    kBuckinghamB_H,
    kBuckinghamB_C,
    kBuckinghamB_N,
    kBuckinghamB_O,
    kBuckinghamB_S,
    kBuckinghamB_Generic,
    kBuckinghamB_Amide15N,
    kBuckinghamB_Amide1HN,
    kBuckinghamB_Carbonyl13C,
    kBuckinghamB_AliphaticSp3_13C,
    kBuckinghamB_Carbonyl17O,
    kBuckinghamB_Hydroxyl17O,
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

inline LiteratureStatusCounts CountBuckinghamConstantStatuses() {
    LiteratureStatusCounts counts;
    for (const LiteratureConstant& c : kAllLiteratureConstants) {
        if (std::string_view(c.key).rfind("buckingham.", 0) != 0) continue;
        switch (c.status) {
        case LiteratureStatus::Cited: ++counts.cited; break;
        case LiteratureStatus::GoodEnough: ++counts.good_enough; break;
        case LiteratureStatus::Placeholder: ++counts.placeholder; break;
        }
    }
    return counts;
}

inline constexpr bool CarbonylLikeFrame(std::string_view frameKind) {
    return frameKind == "backbone_carbonyl"
           || frameKind == "sidechain_carboxylate"
           || frameKind == "sidechain_carboxamide";
}

inline constexpr bool AliphaticCarbonFrame(std::string_view frameKind) {
    return frameKind == "aliphatic_carbon" || frameKind == "backbone_ca";
}

inline constexpr bool AmideNFrame(std::string_view frameKind) {
    return frameKind == "backbone_amide_n" || frameKind == "sidechain_carboxamide";
}

inline constexpr bool AmideHFrame(std::string_view frameKind) {
    return frameKind == "backbone_amide_h";
}

inline constexpr bool HydroxylOxygenFrame(std::string_view frameKind) {
    return frameKind == "hydroxyl_oxygen";
}

inline constexpr bool SidechainCarbonylCarbon(std::string_view residueType,
                                              std::string_view atomName) {
    return (residueType == "ASP" && atomName == "CG")
           || (residueType == "GLU" && atomName == "CD")
           || (residueType == "ASN" && atomName == "CG")
           || (residueType == "GLN" && atomName == "CD");
}

inline constexpr bool SidechainCarbonylOxygen(std::string_view residueType,
                                              std::string_view atomName) {
    return (residueType == "ASP" && (atomName == "OD1" || atomName == "OD2"))
           || (residueType == "GLU" && (atomName == "OE1" || atomName == "OE2"))
           || (residueType == "ASN" && atomName == "OD1")
           || (residueType == "GLN" && atomName == "OE1");
}

inline constexpr bool SidechainAmideNitrogen(std::string_view residueType,
                                             std::string_view atomName) {
    return (residueType == "ASN" && atomName == "ND2")
           || (residueType == "GLN" && atomName == "NE2");
}

inline constexpr bool HydroxylOxygenName(std::string_view atomName) {
    return atomName == "OG" || atomName == "OG1" || atomName == "OH";
}

inline LiteratureConstant BuckinghamA(model::Element element,
                                      std::string_view residueType,
                                      std::string_view atomName,
                                      std::string_view frameKind) {
    switch (element) {
    case model::Element::H:
        if (AmideHFrame(frameKind) || atomName == "H" || atomName == "HN") return kBuckinghamA_Amide1HN;
        return kBuckinghamA_H;
    case model::Element::C:
        if (CarbonylLikeFrame(frameKind) || atomName == "C"
            || SidechainCarbonylCarbon(residueType, atomName))
            return kBuckinghamA_Carbonyl13C;
        if (AliphaticCarbonFrame(frameKind) || atomName == "CA") return kBuckinghamA_AliphaticSp3_13C;
        return kBuckinghamA_C;
    case model::Element::N:
        if (AmideNFrame(frameKind) || atomName == "N"
            || SidechainAmideNitrogen(residueType, atomName))
            return kBuckinghamA_Amide15N;
        return kBuckinghamA_N;
    case model::Element::O:
        if (CarbonylLikeFrame(frameKind) || atomName == "O" || atomName == "OXT"
            || SidechainCarbonylOxygen(residueType, atomName))
            return kBuckinghamA_Carbonyl17O;
        if (HydroxylOxygenFrame(frameKind) || HydroxylOxygenName(atomName)) return kBuckinghamA_Hydroxyl17O;
        return kBuckinghamA_O;
    case model::Element::S: return kBuckinghamA_S;
    default: return kBuckinghamA_Generic;
    }
}

inline LiteratureConstant BuckinghamA(model::Element element,
                                      std::string_view residueType,
                                      std::string_view atomName) {
    return BuckinghamA(element, residueType, atomName, {});
}

inline LiteratureConstant BuckinghamB(model::Element element,
                                      std::string_view residueType,
                                      std::string_view atomName,
                                      std::string_view frameKind) {
    switch (element) {
    case model::Element::H:
        if (AmideHFrame(frameKind) || atomName == "H" || atomName == "HN") return kBuckinghamB_Amide1HN;
        return kBuckinghamB_H;
    case model::Element::C:
        if (CarbonylLikeFrame(frameKind) || atomName == "C"
            || SidechainCarbonylCarbon(residueType, atomName))
            return kBuckinghamB_Carbonyl13C;
        if (AliphaticCarbonFrame(frameKind) || atomName == "CA") return kBuckinghamB_AliphaticSp3_13C;
        return kBuckinghamB_C;
    case model::Element::N:
        if (AmideNFrame(frameKind) || atomName == "N"
            || SidechainAmideNitrogen(residueType, atomName))
            return kBuckinghamB_Amide15N;
        return kBuckinghamB_N;
    case model::Element::O:
        if (CarbonylLikeFrame(frameKind) || atomName == "O" || atomName == "OXT"
            || SidechainCarbonylOxygen(residueType, atomName))
            return kBuckinghamB_Carbonyl17O;
        if (HydroxylOxygenFrame(frameKind) || HydroxylOxygenName(atomName)) return kBuckinghamB_Hydroxyl17O;
        return kBuckinghamB_O;
    case model::Element::S: return kBuckinghamB_S;
    default: return kBuckinghamB_Generic;
    }
}

inline LiteratureConstant BuckinghamB(model::Element element,
                                      std::string_view residueType,
                                      std::string_view atomName) {
    return BuckinghamB(element, residueType, atomName, {});
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

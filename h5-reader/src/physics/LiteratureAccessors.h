#pragma once

#include "../model/Types.h"
#include "constants/LiteratureConstants.h"

#include <cstddef>
#include <string_view>

namespace h5reader::physics {

using nmr::constants::LiteratureConstant;
using nmr::constants::LiteratureStatus;
using nmr::constants::LiteratureStatusCounts;

inline LiteratureStatusCounts CountBuckinghamConstantStatuses() {
    LiteratureStatusCounts counts;
    for (const LiteratureConstant& c : nmr::constants::kAllLiteratureConstants) {
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
    using namespace nmr::constants;
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
    using namespace nmr::constants;
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
    return nmr::constants::kCoulombKeVA.value;
}

inline LiteratureConstant Sigma0(model::Element element,
                                 std::string_view /*residueType*/,
                                 std::string_view /*atomName*/) {
    using namespace nmr::constants;
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
    using namespace nmr::constants;
    switch (category) {
    case model::BondCategory::PeptideCO: return kMcConnellPeptideCO;
    case model::BondCategory::PeptideCN: return kMcConnellPeptideCN;
    case model::BondCategory::BackboneOther: return kMcConnellBackboneOther;
    case model::BondCategory::SidechainCO: return kMcConnellSidechainCO;
    case model::BondCategory::SidechainOther: return kMcConnellSidechainOther;
    case model::BondCategory::Aromatic: return kMcConnellAromaticZero;
    default: return kFallbackLiteratureConstants[0];
    }
}

inline LiteratureConstant RingIntensity(model::RingTypeIndex type) {
    const int idx = static_cast<int>(type);
    if (idx >= 0 && idx < static_cast<int>(nmr::constants::kRingIntensityByType.size()))
        return nmr::constants::kRingIntensityByType[static_cast<std::size_t>(idx)];
    return nmr::constants::kFallbackLiteratureConstants[1];
}

inline LiteratureConstant JohnsonBoveyLobeOffset(model::RingTypeIndex type) {
    const int idx = static_cast<int>(type);
    if (idx >= 0 && idx < static_cast<int>(nmr::constants::kJohnsonBoveyLobeOffsetByType.size()))
        return nmr::constants::kJohnsonBoveyLobeOffsetByType[static_cast<std::size_t>(idx)];
    return nmr::constants::kFallbackLiteratureConstants[2];
}

}  // namespace h5reader::physics

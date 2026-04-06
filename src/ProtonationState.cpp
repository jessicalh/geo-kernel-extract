#include "ProtonationState.h"
#include "AminoAcidType.h"
#include <algorithm>
#include <sstream>

namespace nmr {

void ProtonationState::AddResidue(ResidueProtonation decision) {
    // Replace if a decision for this residue already exists.
    for (auto& r : residues_) {
        if (r.residue_index == decision.residue_index) {
            r = decision;
            return;
        }
    }
    residues_.push_back(decision);
}


const ResidueProtonation* ProtonationState::ForResidue(size_t residue_index) const {
    for (const auto& r : residues_) {
        if (r.residue_index == residue_index) return &r;
    }
    return nullptr;
}


int ProtonationState::NetDecisionCharge() const {
    int charge = 0;
    for (const auto& r : residues_) {
        if (r.variant_index >= 0) {
            const AminoAcidType& aat = GetAminoAcidType(r.amino_acid);
            if (r.variant_index < static_cast<int>(aat.variants.size())) {
                charge += aat.variants[r.variant_index].formal_charge;
            }
        } else {
            // Explicit decision to keep default charged state
            charge += GetAminoAcidType(r.amino_acid).charged_formal_charge;
        }
    }
    return charge;
}


int ProtonationState::NetChargeForProtein(
        const std::vector<AminoAcid>& residue_types) const {
    int charge = 0;
    for (size_t ri = 0; ri < residue_types.size(); ++ri) {
        const AminoAcidType& aat = GetAminoAcidType(residue_types[ri]);

        // Non-titratable residues with formal charge still contribute
        // (e.g., ARG is always +1 but has no titratable variants)
        if (!aat.is_titratable) {
            charge += aat.charged_formal_charge;
            continue;
        }

        const ResidueProtonation* decision = ForResidue(ri);
        if (decision && decision->variant_index >= 0 &&
            decision->variant_index < static_cast<int>(aat.variants.size())) {
            charge += aat.variants[decision->variant_index].formal_charge;
        } else {
            // No decision or default state: use charged_formal_charge.
            // At physiological pH: LYS(+1), ASP(-1), GLU(-1),
            // HIS(+1 default but often 0 at pH 7)
            charge += aat.charged_formal_charge;
        }
    }
    return charge;
}


std::string ProtonationState::Describe() const {
    std::ostringstream ss;
    ss << name_;
    if (!std::isnan(pH_)) ss << " (pH " << pH_ << ")";
    ss << ", " << residues_.size() << " decisions";
    ss << ", net decision charge=" << NetDecisionCharge();
    return ss.str();
}


bool ProtonationState::operator==(const ProtonationState& other) const {
    if (residues_.size() != other.residues_.size()) return false;
    for (size_t i = 0; i < residues_.size(); ++i) {
        if (residues_[i].residue_index != other.residues_[i].residue_index) return false;
        if (residues_[i].variant_index != other.residues_[i].variant_index) return false;
    }
    return true;
}

}  // namespace nmr

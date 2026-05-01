#pragma once
//
// Residue: one amino acid at a specific position in a protein chain.
//
// Knows its type, sequence number, chain, and which atoms belong to it.
// Does NOT store positions. Backbone atom indices cached for fast access.
//

#include "Types.h"
#include "AminoAcidType.h"
#include <string>
#include <vector>
#include <limits>

namespace nmr {

enum class ResidueTerminalState {
    Internal,
    NTerminus,
    CTerminus,
    NAndCTerminus,
    Unknown
};

inline const char* ResidueTerminalStateName(ResidueTerminalState state) {
    switch (state) {
        case ResidueTerminalState::Internal:      return "internal";
        case ResidueTerminalState::NTerminus:     return "n_terminus";
        case ResidueTerminalState::CTerminus:     return "c_terminus";
        case ResidueTerminalState::NAndCTerminus: return "n_and_c_terminus";
        case ResidueTerminalState::Unknown:       return "unknown";
    }
    return "unknown";
}

class Residue {
public:
    AminoAcid           type = AminoAcid::Unknown;
    int                 sequence_number = 0;
    std::string         chain_id;
    std::string         insertion_code;
    std::vector<size_t> atom_indices;
    int                 protonation_variant_index = -1;
    bool                protonation_state_resolved = false;
    ResidueTerminalState terminal_state = ResidueTerminalState::Unknown;

    // Backbone atom indices (cached). NONE means not present.
    static constexpr size_t NONE = std::numeric_limits<size_t>::max();
    size_t N  = NONE;
    size_t CA = NONE;
    size_t C  = NONE;
    size_t O  = NONE;
    size_t H  = NONE;   // NONE for PRO
    size_t HA = NONE;   // HA2 for GLY
    size_t CB = NONE;   // NONE for GLY

    // Chi angle atom indices (up to 4)
    struct ChiAtoms {
        size_t a[4] = {NONE, NONE, NONE, NONE};
        bool Valid() const {
            return a[0] != NONE && a[1] != NONE && a[2] != NONE && a[3] != NONE;
        }
    };
    ChiAtoms chi[4];

    // Typed access to amino acid chemistry
    const AminoAcidType& AminoAcidInfo() const { return GetAminoAcidType(type); }

    // Query methods
    bool IsAromatic() const { return AminoAcidInfo().is_aromatic; }
    bool IsTitratable() const { return AminoAcidInfo().is_titratable; }
    bool HasAmideH() const { return AminoAcidInfo().has_amide_H; }
    int ChiAngleCount() const { return AminoAcidInfo().chi_angle_count; }

    // Sequence address (OBJECT_MODEL.md)
    struct SequenceAddress {
        std::string chain_id;
        int sequence_number;
        std::string insertion_code;
    };
    SequenceAddress GetSequenceAddress() const {
        return {chain_id, sequence_number, insertion_code};
    }
};

}  // namespace nmr

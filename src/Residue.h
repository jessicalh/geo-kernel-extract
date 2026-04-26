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

class Residue {
public:
    AminoAcid           type = AminoAcid::Unknown;
    int                 sequence_number = 0;
    std::string         chain_id;
    std::string         insertion_code;
    std::vector<size_t> atom_indices;
    int                 protonation_variant_index = -1;

    // Backbone atom indices (cached). NONE means not present.
    static constexpr size_t NONE = std::numeric_limits<size_t>::max();
    size_t N  = NONE;
    size_t CA = NONE;
    size_t C  = NONE;
    size_t O  = NONE;
    size_t H  = NONE;   // NONE for PRO
    size_t HA = NONE;   // HA2 for GLY
    size_t CB = NONE;   // NONE for GLY

    // Sequence context (set by Protein::StampResidueContext after
    // residues are loaded). Useful for stats consumers, prochirality
    // disambiguation, and N/C-terminal aware feature extraction.
    AminoAcid prev_residue_type = AminoAcid::Unknown;  // Unknown if N-terminal
    AminoAcid next_residue_type = AminoAcid::Unknown;  // Unknown if C-terminal
    bool is_n_terminal = false;
    bool is_c_terminal = false;

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

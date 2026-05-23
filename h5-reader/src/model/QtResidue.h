// QtResidue — one amino acid mirrored from residues.npy.
//
// All chemistry-typed fields (aminoAcid, protonationVariantIndex,
// terminalState, prev/next residue links, prev/next residue types,
// is_xpro_context) are typed. The five string fields the old format
// carried (amber/iupac residue 3-letter, one-letter, chain_id,
// insertion_code) are derived OR wrapped:
//
//   - 3-letter codes (amber, iupac) and 1-letter: DERIVED at display
//     time from QtAminoAcid + QtAminoAcidType registry via the
//     QtResidueNames helpers. Not stored on QtResidue.
//
//   - chain_id + insertion_code: WRAPPED in QtChainAddress, which
//     deliberately deletes operator== so a comparison is always an
//     explicit IsSameAddress() call.
//
// The convenience flags (is_proline, is_aromatic, is_titratable,
// has_amide_h) the sidecar projects are derivable from aminoAcid via
// QtAminoAcidType — store anyway since reading them once from the
// sidecar avoids 50 lookups per protein. is_xpro_context is unique
// (bond-graph derived; not in any AA table).
//
// Backbone index cache (N/CA/C/O/H/HA/CB) is rebuilt at load from
// atoms' typed BackboneRole + Locant — no string scan, no
// "first non-aromatic sidechain C is CB" heuristic. SIZE-checked at
// -1 sentinel (NONE).

#pragma once

#include "QtChainAddress.h"
#include "QtSemanticEnums.h"
#include "Types.h"

#include <cstdint>
#include <vector>

namespace h5reader::model {

struct QtResidue {
    // ----- Identity indices -----
    int32_t residueIndex = -1;

    // ----- Addressing (typed wrapper; operator== deleted) -----
    QtChainAddress address;

    // ----- Residue type + protonation variant -----
    AminoAcid aminoAcid = AminoAcid::Unknown;
    int8_t protonationVariantIndex = -1;
    TerminalState terminalState = TerminalState::Internal;

    // ----- Backbone chain links (bond-graph derived; -1 == no link) -----
    int32_t prevResidueIndex = -1;
    int32_t nextResidueIndex = -1;
    AminoAcid prevResidueType = AminoAcid::Unknown;
    AminoAcid nextResidueType = AminoAcid::Unknown;

    // ----- Sidecar-projected convenience flags (also derivable via
    // QtAminoAcidType; stored to avoid per-call lookup) -----
    int32_t atomCount = 0;
    bool isProline = false;
    bool isAromatic = false;
    bool isTitratable = false;
    bool hasAmideH = false;
    bool isXProContext = false;  // backbone successor is PRO

    // ----- Atom membership (built at load from atoms.residue_index walk) -----
    std::vector<int32_t> atomIndices;

    // ----- Backbone atom-index cache (built at load from typed
    // BackboneRole + Locant; no string scan). NONE = absent.
    static constexpr int32_t NONE = -1;
    int32_t N = NONE;
    int32_t CA = NONE;
    int32_t C = NONE;
    int32_t O = NONE;
    int32_t H = NONE;   // NONE for PRO
    int32_t HA = NONE;  // GLY HA pseudoatom covers both HA2/HA3
    int32_t CB = NONE;  // NONE for GLY

    // ----- Queries -----
    bool HasN() const { return N != NONE; }
    bool HasCA() const { return CA != NONE; }
    bool HasC() const { return C != NONE; }
    bool HasO() const { return O != NONE; }
    bool HasCB() const { return CB != NONE; }
    bool IsTerminalN() const {
        return terminalState == TerminalState::NTerminus || terminalState == TerminalState::NAndCTerminus;
    }
    bool IsTerminalC() const {
        return terminalState == TerminalState::CTerminus || terminalState == TerminalState::NAndCTerminus;
    }
};

}  // namespace h5reader::model

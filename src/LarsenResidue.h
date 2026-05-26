#pragma once
// Typed model of the five pieces perceived from a Larsen ProCS15
// tripeptide DFT record. Runtime chemistry uses AtomMechanicalIdentity;
// dft_atom_idx is preserved only for correlation with the source row.

#include "SemanticEnums.h"
#include "Types.h"

#include <array>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace nmr {

struct TripeptideDftRecord;  // forward; defined in TripeptideDftTable.h


// One residue's worth of perceived atoms in a Larsen tripeptide.
class LarsenResidue {
public:
    enum class Kind : std::uint8_t {
        AceCap   = 0,  ///< Acetyl cap atoms: -C(=O)-CH3.
        NCapAla  = 1,  ///< N-terminal flanking ALA residue.
        Central  = 2,  ///< Central residue (residue field carries which).
        CCapAla  = 3,  ///< C-terminal flanking ALA residue.
        NmeCap   = 4,  ///< N-methylamide cap: -NH-CH3.
    };

    // One perceived atom. The `identity` is the typed chemistry
    // assignment from graph isomorphism; `dft_atom_idx` is the
    // 1-based row index in the Gaussian / ORCA Standard orientation
    // output (or equivalent), preserved for log correlation only.
    struct PerAtom {
        AtomMechanicalIdentity identity = {};
        int                    dft_atom_idx = 0;
        Element                element  = Element::Unknown;
        Vec3                   position = Vec3::Zero();

        // Pre-decomposed shielding tensor (irreps split at ingest
        // and carried on the DB row; preserved here for direct use
        // by calculators).
        Mat3                  shielding_tensor = Mat3::Zero();
        double                isotropic        = 0.0;
        double                anisotropy       = 0.0;
        std::array<double, 5> t2_components    = {};

        // True iff K=3 Weisfeiler-Lehman leaves this atom in a
        // graph-automorphic class, such as PHE/TYR CD1/CD2 or ARG
        // NH1/NH2. The
        // chemical identity tuple (Element, Locant, BackboneRole) is
        // sound; only the within-class label (BranchAddress and
        // DiastereotopicIndex) is interchangeable. Matchers should
        // drop those two fields and resolve by nearest-spatial when
        // this flag is true.
        bool canonical_assignment_ambiguous = false;
    };

    // One bond in this piece's local covalent graph. External amide
    // bonds are represented only by the LarsenTripeptide-level chain.
    struct Bond {
        int          a     = -1;
        int          b     = -1;
        std::uint8_t order = 1;  ///< 1=single, 2=double, etc. (heuristic)
    };

    Kind                 kind     = Kind::Central;
    AminoAcid            residue  = AminoAcid::Unknown;
    std::vector<PerAtom> atoms;
    std::vector<Bond>    bonds;

    // Role-pinned slot cache. -1 if absent for this piece's kind.
    int N_idx  = -1;
    int H_idx  = -1;
    int CA_idx = -1;
    int HA_idx = -1;
    int CB_idx = -1;
    int C_idx  = -1;
    int O_idx  = -1;

    // Find the local index of an atom with the given typed identity.
    // Returns -1 if not found. For chemically-equivalent atom sets
    // (e.g., methyl Hs that collapse to one identity), returns the
    // first match.
    int LookupByIdentity(const AtomMechanicalIdentity& id) const;

    bool HasAllRequiredSlots() const;
};


// The full 5-piece Larsen tripeptide.
struct LarsenTripeptide {
    LarsenResidue ace;
    LarsenResidue n_cap;
    LarsenResidue central;
    LarsenResidue c_cap;
    LarsenResidue nme;

    // HIS central variant matched by perception: 0=HID, 1=HIE, 2=HIP.
    // -1 for non-HIS centrals or when no variant has been determined.
    int central_variant_index = -1;

    int TotalAtoms() const;

    // Find a (piece*, local_idx) by the global dft_atom_idx (1-based
    // in the source Gaussian / ORCA output). Returns (nullptr, -1)
    // if not found.
    std::pair<const LarsenResidue*, int> FindByDftIdx(int dft_atom_idx) const;
};


// Perceive a full LarsenTripeptide from a TripeptideDftRecord.
//
// Verifies that the central piece matches expected_central's canonical
// chemistry; mismatches return nullopt.
//
// his_variant_hint applies only when expected_central == HIS. Values
// are the canonical protonation variant indices: 0=HID, 1=HIE, 2=HIP.
// -1 means "no hint" (default). When a hint is provided, perception
// requires that variant to match. Without a hint, perception tries
// HID/HIE/HIP in order and accepts the first that fits, which can
// misassign histidine ring protons if the protein's actual variant
// differs from the perceived one. Pass
// the protein's `Residue::protonation_variant_index` to keep the
// match honest.
//
// On failure, returns nullopt. The reason is logged via OperationLog.
std::optional<LarsenTripeptide> PerceiveLarsenTripeptide(
    const TripeptideDftRecord& rec,
    AminoAcid                  expected_central,
    int                        his_variant_hint = -1);


}  // namespace nmr

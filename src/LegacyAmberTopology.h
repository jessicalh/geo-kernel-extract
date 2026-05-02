#pragma once
//
// LegacyAmberTopology: explicit topology contract for the current
// Amber-derived calculator language.
//
// This is the string barrier output. Loaders and construction code may use
// names to build this object; calculators consume its typed/indexed facts.
//
// Two scopes of invariant data live here:
//   - structural identity (atom_count, residue_count, bonds) — set at
//     construction by Protein::FinalizeConstruction; available on every
//     load path.
//   - FF-numerical (mass, partial charges read by GROMACS, atom types,
//     LJ params, exclusions, fudgeQQ, SETTLE / virtual site geometry)
//     — populated by AttachAmberFFData when the load path is a TPR
//     (FullSystemReader). Calculators that need these fields check
//     HasAmberFFData(); the legacy ff14SB-flat-table path leaves the
//     enrichment empty.
//

#include "ProteinTopology.h"
#include "CovalentTopology.h"
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace nmr {

// Pair-LJ parameters for a single atom-type pair.
// GROMACS stores c6 and c12 directly (kJ/mol·nm^6 and kJ/mol·nm^12);
// sigma/epsilon are derived if needed.
struct AmberLjPair {
    double c6 = 0.0;
    double c12 = 0.0;
};

// Per-atom invariant FF-numerical fields populated from the TPR.
// Vectors are parallel to the protein's typed atom list.
struct AmberFFData {
    // Per-atom invariant FF-numerical fields (vectors size = atom_count).
    std::vector<double> mass;                   // Da
    std::vector<int>    ff_atom_type_index;     // index into the FF atom-type table
    std::vector<int>    ptype;                  // GROMACS ParticleType enum (cast to int)
    std::vector<std::string> atomtype_string;   // FF atom-type name (e.g. "CT", "N", "H1")

    // Exclusion lists from gmx_localtop_t.excls, filtered to protein.
    // Each inner vector lists atom indices excluded from the index's
    // nonbonded interactions (1-2, 1-3, 1-4 within nrexcl).
    std::vector<std::vector<int>> exclusions;

    // Pair-LJ parameter table. Rectangular atnr × atnr; LJ is symmetric
    // in i,j. Index as lj_pair_table[i * atnr + j].
    std::vector<AmberLjPair> lj_pair_table;
    int                       atnr = 0;          // atom-type table dimension

    // Force-field-wide constants.
    double fudge_qq = 1.0;                       // 1-4 Coulomb scaling
    double rep_pow  = 12.0;                      // LJ repulsion exponent (provenance)

    // Number of non-perturbed interactions per InteractionFunction kind,
    // copied from idef.numNonperturbedInteractions. Non-zero values past
    // this index in idef.il[fn] would indicate FEP perturbations baked
    // into the TPR. Surfaced for diagnostics and logonlyloadrecord; we
    // do not currently consume FEP perturbations.
    std::array<int, 256> num_non_perturbed = {};

    // SETTLE constraint geometry (per water moltype, one entry per
    // distinct water model that appeared in the TPR). Empty if no
    // water moltype was present.
    struct SettleEntry {
        std::string moltype_name;   // e.g. "SOL"
        double doh = 0.0;           // OH bond length, nm
        double dhh = 0.0;           // HH distance, nm
    };
    std::vector<SettleEntry> settle;

    // Virtual site geometry. Empty for TIP3P fixtures (no vsites);
    // populated when a TPR has TIP4P+ water (M-site) or a vsite-
    // bearing protein moltype. For the latter case we'd flag in the
    // load record; this slot just preserves the geometric constraint
    // descriptor so a future reader has it.
    struct VirtualSiteEntry {
        size_t atom_index = 0;       // protein-local index of the vsite atom
        int    construction_type = 0; // GROMACS InteractionFunction enum value
        std::array<size_t, 5> source_atoms = {};
        int    source_atom_count = 0;
        std::array<double, 6> params = {};
    };
    std::vector<VirtualSiteEntry> virtual_sites;
};


class LegacyAmberTopology final : public ProteinTopology {
public:
    LegacyAmberTopology(size_t atom_count,
                        size_t residue_count,
                        std::unique_ptr<CovalentTopology> bonds);

    ProteinTopologyKind Kind() const override {
        return ProteinTopologyKind::LegacyAmber;
    }
    std::string_view Name() const override { return "LegacyAmberTopology"; }
    size_t AtomCount() const override { return atom_count_; }
    size_t ResidueCount() const override { return residue_count_; }

    const CovalentTopology& Bonds() const { return *bonds_; }

    size_t BondCount() const { return bonds_->BondCount(); }
    const Bond& BondAt(size_t i) const { return bonds_->BondAt(i); }
    const std::vector<Bond>& BondList() const { return bonds_->Bonds(); }
    const std::vector<size_t>& BondIndicesFor(size_t atom_index) const {
        return bonds_->BondIndicesFor(atom_index);
    }
    size_t HydrogenParentOf(size_t atom_index) const {
        return bonds_->HydrogenParentOf(atom_index);
    }

    // FF-numerical enrichment. AttachAmberFFData is a one-shot
    // post-construction populate; calling it twice is a logic error.
    // FullSystemReader calls it after Protein::FinalizeConstruction
    // creates the topology so the FF-numerical data extracted from the
    // TPR can be attached to the structural topology.
    void AttachAmberFFData(AmberFFData data);
    bool HasAmberFFData() const { return amber_ff_data_.has_value(); }
    const AmberFFData& FFData() const;

private:
    size_t atom_count_ = 0;
    size_t residue_count_ = 0;
    std::unique_ptr<CovalentTopology> bonds_;
    std::optional<AmberFFData> amber_ff_data_;
};

}  // namespace nmr

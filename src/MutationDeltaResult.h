#pragma once
//
// MutationDeltaResult: ConformationResult that stores the per-atom deltas
// between this (WT) conformation and a mutant conformation.
//
// Attaches to the WT conformation. The mutant conformation is an extra
// input to Compute, not a dependency. Stores atom correspondence and
// per-atom delta tensors (DFT shielding, APBS fields, MOPAC charges,
// DSSP, graph, charges) internally — does NOT write to ConformationAtom.
//
// Ring proximity: for each matched atom, stores cylindrical coordinates
// relative to each removed ring (the mutation site rings). Uses SVD-based
// ring normals from GeometryResult, cylindrical decomposition (z, rho,
// theta), McConnell factor, and exponential decay weighting.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>
#include <limits>

namespace nmr {

class ProteinConformation;


struct MutationSite {
    size_t residue_index = 0;
    AminoAcid wt_type = AminoAcid::Unknown;
    AminoAcid mut_type = AminoAcid::Unknown;
    std::vector<size_t> wt_ring_indices;  // rings that were removed
};


// Cylindrical coordinates of an atom relative to one removed ring
struct RingProximity {
    size_t ring_index = 0;
    RingTypeIndex ring_type = RingTypeIndex::PheBenzene;
    double distance = 0.0;          // |d| to ring center (Angstroms)
    double z = 0.0;                 // height above ring plane (signed)
    double rho = 0.0;               // radial distance in ring plane
    double theta = 0.0;             // cone angle from normal: atan2(rho, z)
    double mcconnell_factor = 0.0;  // (3cos^2 theta - 1) / r^3
    double exp_decay = 0.0;         // exp(-distance / 4.0)
};


struct MatchedAtomData {
    size_t wt_index = 0;
    size_t mut_index = 0;
    double match_distance = 0.0;

    // WT atom identity (from typed objects, not strings)
    Element element = Element::Unknown;
    AtomRole role = AtomRole::Unknown;
    bool is_backbone = false;
    size_t residue_index = 0;

    // DFT shielding delta (always present — OrcaShieldingResult required)
    Mat3 delta_shielding = Mat3::Zero();
    SphericalTensor delta_shielding_spherical;

    // APBS delta
    Vec3 delta_efield = Vec3::Zero();
    Mat3 delta_efg = Mat3::Zero();
    SphericalTensor delta_efg_spherical;
    bool has_apbs_delta = false;

    // Charge deltas
    double delta_partial_charge = 0.0;  // ff14SB charge delta
    double delta_mopac_charge = 0.0;
    bool has_mopac_delta = false;

    // DSSP deltas (per residue, stored on the atom for convenience)
    double delta_phi = 0.0;
    double delta_psi = 0.0;
    double delta_sasa = 0.0;
    bool has_dssp_delta = false;

    // Graph deltas
    int delta_graph_dist_ring = 0;  // WT - mutant (positive = was closer in WT)
    double delta_bfs_decay = 0.0;
    bool has_graph_delta = false;

    // Ring proximity to removed rings (cylindrical coords in ring frame)
    std::vector<RingProximity> removed_ring_proximity;
    size_t nearest_removed_ring = SIZE_MAX;   // index into removed_ring_proximity
    double nearest_removed_ring_dist = 99.0;
};


// Aggregate summary statistics
struct DeltaSummary {
    // By element
    struct ElementBin {
        Element element = Element::Unknown;
        int count = 0;
        double mean_delta_t0 = 0.0;
        double mean_abs_delta_t0 = 0.0;
        double max_abs_delta_t0 = 0.0;
        double mean_t2_magnitude = 0.0;
    };
    std::vector<ElementBin> by_element;

    // By distance to nearest removed ring (1A bins from 0-15A)
    struct DistanceBin {
        double bin_start = 0.0;  // Angstroms
        double bin_end = 0.0;
        int count = 0;
        double mean_abs_delta_t0 = 0.0;
        double mean_t2_magnitude = 0.0;
    };
    std::vector<DistanceBin> by_distance;

    // Backbone vs sidechain
    int backbone_count = 0;
    double backbone_mean_abs_t0 = 0.0;
    int sidechain_count = 0;
    double sidechain_mean_abs_t0 = 0.0;
};


class MutationDeltaResult : public ConformationResult {
public:
    std::string Name() const override { return "MutationDeltaResult"; }

    std::vector<std::type_index> Dependencies() const override;

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    static std::unique_ptr<MutationDeltaResult> Compute(
        ProteinConformation& wt_conf,
        const ProteinConformation& mut_conf);

    // --- Counts ---
    size_t MatchedAtomCount() const { return matched_atoms_.size(); }
    size_t UnmatchedWtAtomCount() const;

    // --- Per-atom queries (return zero for unmatched atoms) ---
    bool HasMatch(size_t wt_atom_index) const;
    const MatchedAtomData& MatchedDataAt(size_t wt_atom_index) const;
    const Mat3& DeltaShieldingAt(size_t wt_atom_index) const;
    const SphericalTensor& DeltaShieldingSphericalAt(size_t wt_atom_index) const;
    double DeltaT0At(size_t wt_atom_index) const;

    // APBS deltas
    Vec3 DeltaEFieldAt(size_t wt_atom_index) const;
    Mat3 DeltaEFGAt(size_t wt_atom_index) const;
    bool HasApbsDelta() const { return has_apbs_delta_; }

    // Charge deltas
    double DeltaPartialChargeAt(size_t wt_atom_index) const;
    double DeltaMopacChargeAt(size_t wt_atom_index) const;
    bool HasMopacDelta() const { return has_mopac_delta_; }

    // DSSP deltas
    bool HasDsspDelta() const { return has_dssp_delta_; }

    // Graph deltas
    bool HasGraphDelta() const { return has_graph_delta_; }

    // --- Mutation sites ---
    const std::vector<MutationSite>& MutationSites() const { return mutation_sites_; }

    // --- Ring proximity ---
    double NearestRemovedRingDistance(size_t wt_atom_index) const;

    // --- Summary statistics ---
    const DeltaSummary& Summary() const { return summary_; }

    // --- Atom correspondence ---
    size_t MutantAtomFor(size_t wt_atom_index) const;
    double MatchDistanceAt(size_t wt_atom_index) const;

private:
    const ProteinConformation* wt_conf_ = nullptr;

    std::vector<size_t> wt_to_matched_;  // SIZE_MAX if unmatched
    std::vector<MatchedAtomData> matched_atoms_;
    std::vector<MutationSite> mutation_sites_;
    DeltaSummary summary_;

    bool has_apbs_delta_ = false;
    bool has_mopac_delta_ = false;
    bool has_dssp_delta_ = false;
    bool has_graph_delta_ = false;

    static const Mat3 zero_mat3_;
    static const SphericalTensor zero_spherical_;
    static const MatchedAtomData empty_match_;
};

}  // namespace nmr

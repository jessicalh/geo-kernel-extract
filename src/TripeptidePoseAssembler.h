#pragma once
// Aligns a TripeptideDftRecord onto a protein residue and maps DFT
// tensors to protein atoms by typed identity. Central-residue residuals
// are emitted as features; cap residuals above the threshold are
// rejected.

#include "Types.h"
#include "TripeptideDftTable.h"

#include <cstdint>
#include <string>
#include <vector>

namespace nmr {

class Protein;
class ProteinConformation;


// NTerm reads the N-terminal ALA cap of the next-centered tripeptide;
// CTerm reads the C-terminal ALA cap of the previous-centered tripeptide.
enum class TripeptidePoseSide : std::uint8_t {
    Central = 0,
    NTerm   = 1,
    CTerm   = 2,
};


struct AlignedDftAtom {
    int         dft_atom_idx     = -1;          // into TripeptideDftRecord::atoms
    std::size_t protein_atom_idx = 0;            // into Protein::atoms
    Element     element          = Element::Unknown;

    Vec3 aligned_position = Vec3::Zero();         // R*(p_dft - cs) + ct
    Vec3 residual_vec     = Vec3::Zero();         // aligned - protein_pos
    double residual_distance = 0.0;                // residual_vec.norm()

    Mat3            shielding_tensor_aligned    = Mat3::Zero();    // R σ R^T, ppm
    SphericalTensor shielding_spherical_aligned;                    // Decompose

    bool substrate_role_agrees = false;
};


struct AssembledTripeptide {
    bool ok = false;                          // false on hard structural failure
    int  calc_id = 0;
    std::string frame_type;
    TripeptidePoseSide side = TripeptidePoseSide::Central;
    int  n_chi_axes_used = 0;                 // caller-set diagnostic

    Mat3   rotation        = Mat3::Identity();
    Vec3   source_centroid = Vec3::Zero();
    Vec3   target_centroid = Vec3::Zero();
    double backbone_kabsch_rmsd = 0.0;         // 3-point Kabsch fit, Å

    std::vector<AlignedDftAtom> aligned_atoms;

    int    n_substrate_disagreements = 0;
    int    n_above_threshold = 0;             // cap rejects; central counts only.
    double max_residual_A = 0.0;
    double mean_residual_A = 0.0;
};


// validation_threshold_A rejects cap atoms but only increments
// n_above_threshold on the central path. substrate_check_strict applies
// to cap slot role checks; central matching uses typed identity before
// emission. ok=false means the record was absent, required N/CA/C or
// cap slots were missing, perception failed, or no atom was emitted.
AssembledTripeptide AssembleTripeptide(
    const Protein&             protein,
    const ProteinConformation& conf,
    std::size_t                protein_residue_idx,
    const TripeptideDftRecord& rec,
    TripeptidePoseSide         side,
    double                     validation_threshold_A = 3.0,
    bool                       substrate_check_strict = true);


}  // namespace nmr

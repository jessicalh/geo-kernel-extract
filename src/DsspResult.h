#pragma once
//
// DsspResult: secondary structure, phi/psi, SASA from DSSP.
//
// Uses libdssp (Joosten et al. 2011, from Kabsch & Sander 1983).
// Per-residue: secondary structure char, raw libdssp/IUPAC phi (rad),
// raw libdssp/IUPAC psi (rad), SASA (A^2), H-bond acceptor and donor partners.
//
// Dependencies: none (DSSP needs only positions).
//
// CROSS-RESULT READ (writer side, 2026-05-19, PATTERNS §17):
//   AllResidues() (secondary_structure + acceptors[2] / donors[2] +
//   their energies) is read per-frame by
//   Dssp8TimeSeriesTrajectoryResult (emits SS code + H-bond partner/
//   energy timelines) and Dssp8TransitionTrajectoryResult (emits SS
//   transition stats). Both reader TRs declare no Dependencies()
//   because DsspResult attaches conditionally on skip_dssp == false;
//   per-frame source-attached gate captures absence.
//
//   Phi() / Psi() are read by MutationDeltaResult when both WT and
//   mutant conformations have DsspResult attached. They are also read
//   by DihedralTimeSeriesTrajectoryResult's cross-result numerical-
//   consistency test; libdssp returns plain IUPAC phi/psi, while the TR
//   uses the project convention (negative of IUPAC). The test therefore
//   compares the TR against -DsspResult.Phi/Psi at 1e-2 rad tolerance.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include <cstdint>
#include <vector>

namespace nmr {

struct DsspResidue {
    // observed = true iff a DSSP output row actually mapped to this
    // residue index. The result vector is resized to ResidueCount() at
    // Compute time, but only residues whose (chain_id, seq_num) lookup
    // hits in the cif++/libdssp output get populated. Without this
    // flag, unmapped residues (caps, unmapped chains, insertion-code
    // collisions, DSSP skips) would silently look like real coil "C" —
    // the secondary_structure default, biasing downstream coil
    // occupancy and helix-dwell statistics. Codex review caught this
    // 2026-05-19; the failure mode is invisible on linear-single-chain
    // test fixtures but real on multi-chain / capped / engineered
    // structures the fleet will see.
    bool   observed = false;
    char secondary_structure = 'C';  // H/G/I/E/B/T/S/C (valid only when observed)
    double phi = 0.0;                // libdssp/IUPAC radians; observed only
    double psi = 0.0;                // libdssp/IUPAC radians; observed only
    double sasa = 0.0;               // Angstroms^2 (valid only when observed)

    // H-bond partner indices (into protein residue list, SIZE_MAX if none)
    struct HBondPartner {
        size_t residue_index = SIZE_MAX;
        double energy = 0.0;         // kcal/mol
    };
    HBondPartner acceptors[2];       // residues whose C=O accepts this N-H
    HBondPartner donors[2];          // residues whose N-H donates to this C=O
};

inline bool DsspIsPpii(char ss) { return ss == 'P'; }


class DsspResult : public ConformationResult {
public:
    std::string Name() const override { return "DsspResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }
    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    // Factory: run DSSP on the conformation and return the result.
    // Returns nullptr on failure (writes diagnostic to stderr).
    static std::unique_ptr<DsspResult> Compute(ProteinConformation& conf);
    static std::unique_ptr<DsspResult> CreateForTesting(
        std::vector<DsspResidue> residues);

    // Physics query methods
    char SecondaryStructure(size_t residue_index) const;
    double Phi(size_t residue_index) const;
    double Psi(size_t residue_index) const;
    double SASA(size_t residue_index) const;

    // All per-residue data
    const std::vector<DsspResidue>& AllResidues() const { return residues_; }

private:
    void PopulateCanonicalTorsions(const ProteinConformation* conf);

    std::vector<DsspResidue> residues_;

    // Residue-major matrices, six columns in canonical order:
    // phi, psi, chi1, chi2, chi3, chi4.
    std::vector<double> torsion_angle_;
    std::vector<double> torsion_sin_;
    std::vector<double> torsion_cos_;
    std::vector<std::uint8_t> torsion_valid_;
};

}  // namespace nmr

#pragma once
//
// FullSystemReader: reads full-system GROMACS topology and trajectory
// frames, splitting atoms into protein + solvent.
//
// The TPR contains the topology for ALL atoms (protein, water, ions).
// The full-system .xtc contains positions for ALL atoms at each frame.
// This reader:
//   1. Reads the TPR once — extracts atom ranges, bonded interaction
//      parameters, and everything needed to build the Protein object
//   2. Splits full-system .xtc frames: protein → Vec3, solvent → SolventEnvironment
//
// Single parse: ReadTopology() parses the TPR once and populates
// SystemTopology + BondedParameters + stored protein topology data.
// Everything downstream uses stored results — no re-reads.
//
// Coordinates are in Angstroms (converted from nm).
//

#include "SolventEnvironment.h"
#include "BondedEnergyResult.h"
#include "LegacyAmberTopology.h"  // AmberFFData enrichment payload
#include "OperationRunner.h"      // ForceField enum

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace nmr {

struct BuildResult;                   // forward — defined in BuildResult.h
namespace detail { struct TprData; }  // opaque — holds gmx_mtop_t, defined in .cpp

// Topology layout extracted from TPR.
struct SystemTopology {
    // Atom ranges in the full-system frame.
    size_t protein_start = 0;
    size_t protein_count = 0;
    size_t water_O_start = 0;   // first water oxygen index
    size_t water_count = 0;     // number of water MOLECULES (3 atoms each)
    size_t ion_start = 0;
    size_t ion_count = 0;
    size_t total_atoms = 0;

    // Per-water-molecule charge (from first water in topology).
    double water_O_charge = 0.0;
    double water_H_charge = 0.0;

    // Per-ion charges and elements (one per ion).
    std::vector<double> ion_charges;
    std::vector<int>    ion_atomic_numbers;
};


class FullSystemReader {
public:
    FullSystemReader();
    ~FullSystemReader();

    // Not copyable (holds unique_ptr to parsed TPR data).
    FullSystemReader(const FullSystemReader&) = delete;
    FullSystemReader& operator=(const FullSystemReader&) = delete;

    // Read the TPR once. Populates SystemTopology (atom ranges),
    // BondedParameters (interaction lists + CMAP grids), and stores
    // parsed topology for BuildProtein(). Single parse — everything
    // downstream uses stored results.
    // Returns false on error (check error()).
    bool ReadTopology(const std::string& tpr_path);

    // Given a full-system XTC frame (all atoms, in nm),
    // extract protein positions (in Angstroms) and SolventEnvironment.
    // protein_positions will have protein_count entries.
    bool ExtractFrame(const std::vector<float>& full_frame_xyz,
                      std::vector<Vec3>& protein_positions,
                      SolventEnvironment& solvent) const;

    // Bonded interaction parameters extracted during ReadTopology().
    // Bond, angle, UB, proper dihedral, improper dihedral, and CMAP
    // interactions for the protein atoms. Atom indices are protein-local
    // (0-based, same as our Protein).
    const BondedParameters& BondedParams() const { return bonded_params_; }
    bool HasBondedParams() const { return !bonded_params_.interactions.empty(); }

    // FF-numerical enrichment data extracted during ReadTopology().
    // Per-atom mass, atom-type index, ptype, atomtype string; exclusion
    // lists; LJ pair table; fudgeQQ; SETTLE per water moltype; vsite
    // geometry. Consumed by the trajectory loader path: after
    // Protein::FinalizeConstruction creates the LegacyAmberTopology,
    // the loader calls protein.MutableLegacyAmber().AttachAmberFFData(
    // reader.ConsumeAmberFFData()).
    const AmberFFData& AmberFF() const { return amber_ff_data_; }
    AmberFFData ConsumeAmberFFData() { return std::move(amber_ff_data_); }

    // Build a Protein + ChargeSource from the stored TPR parse.
    // Must be called after ReadTopology(). Returns a Protein with
    // residues + atoms but no conformations (FinalizeConstruction
    // needs positions from the first XTC frame).
    BuildResult BuildProtein(const std::string& protein_id,
                             ForceField force_field = ForceField::CHARMM36m) const;

    // PBC-fix the protein slice in-place. `protein_coords` is the
    // contiguous protein-only float coordinate buffer (size =
    // 3 * Topology().protein_count, typically carved from a full XTC
    // frame). `box_in` is this frame's box (varies per frame in NPT).
    // Walks the protein-only mtop (built once at ReadTopology time)
    // so the bond-graph traversal touches only the protein atoms.
    // Returns false if ReadTopology has not run, the buffer size is
    // wrong, or the captured pbcType was unset.
    bool MakeProteinWhole(std::vector<float>& protein_coords,
                          const float box_in[3][3]) const;

    const SystemTopology& Topology() const { return topo_; }
    const std::string& error() const { return error_; }

private:
    SystemTopology topo_;
    BondedParameters bonded_params_;
    AmberFFData amber_ff_data_;
    std::unique_ptr<detail::TprData> tpr_;
    std::string error_;
};

}  // namespace nmr

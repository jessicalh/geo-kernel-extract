#pragma once
/// @file
/// Typed per-mode input specifications and the @ref ModeSpec variant.
///
/// Each supported nmr_extract mode is a distinct struct holding only
/// the fields that mode consumes. A @c std::variant binds them into
/// the @ref ModeSpec dispatch type. The runner uses @c std::visit
/// rather than a JobMode-enum switch.

#include "FrameEmission.h"
#include "OrcaRunLoader.h"

#include <filesystem>
#include <optional>
#include <variant>

namespace nmr::cli {

/// @brief Load a bare PDB, protonate with @c reduce, run all calculators.
///
/// Heavy-atom-only PDB on input. Hydrogens are added at runtime via the
/// @c reduce tool, charges from ff14SB, then the full per-conformation
/// pipeline runs and emits per-atom NPYs.
struct PdbMode {
    std::filesystem::path pdb;            ///< Path to the bare PDB.
    double                pH      = 7.0;  ///< pH passed to @c reduce.
    bool                  mopac   = true; ///< Run PM7+MOZYME + derived calcs.
    bool                  apbs    = true; ///< Run APBS solvated EFG.
    bool                  coulomb = true; ///< Run vacuum Coulomb EFG.
};

/// @brief Load a PDB that already has H atoms, run all calculators.
///
/// Skips the @c reduce step; protonation variant is detected from the
/// explicit H atoms present.
struct ProtonatedPdbMode {
    std::filesystem::path pdb;
    bool                  mopac   = true;
    bool                  apbs    = true;
    bool                  coulomb = true;
};

/// @brief Load a single tleap/AMBER-prepared pose, run all calculators.
///
/// The pose is described by a root name: @c {root}.xyz (coordinates),
/// @c {root}.prmtop (topology + charges), and optional
/// @c {root}_nmr.out (ORCA DFT shielding tensors for OrcaShieldingResult).
struct OrcaMode {
    OrcaRunFiles files;
    bool         mopac   = true;
    bool         apbs    = true;
    bool         coulomb = true;
};

/// @brief Load a WT+ALA mutant pair, run on both, compute Δ shielding.
///
/// Each side is a pose with the OrcaMode shape. MutationDeltaResult
/// attaches to the WT conformation and stores per-atom delta tensors.
struct MutantMode {
    OrcaRunFiles wt;
    OrcaRunFiles ala;
    bool         mopac   = true;
    bool         apbs    = true;
    bool         coulomb = true;
};

/// @brief Read a GROMACS relaxation run, per-frame calculators + H5.
///
/// @c dir contains @c production.tpr, @c production.trr, @c production.edr.
/// The parent of @c dir is searched for @c topol.top by convention
/// (GROMACS readback). Emits the trajectory H5 unconditionally; per-frame
/// PDBs (@c emit_pdbs) and per-frame NPYs (@c emit_npys) are independent
/// opt-ins.
///
/// MOPAC default is off here because the per-frame cost
/// (~10 min/frame for PM7+MOZYME) makes fleet-scale trajectory runs
/// infeasible.
struct TrajectoryMode {
    std::filesystem::path           dir;
    bool                            mopac = false;
    /// @brief Stride at which MOPAC runs (in TRR-frame-index units).
    ///
    /// When @c mopac is true, MOPAC is invoked only on frames where
    /// @c frame_idx % mopac_stride == 0. Other dispatched frames run
    /// the rest of the pipeline without MopacResult attached (per the
    /// conditional-attach TR discipline). Should match the NPY-emit
    /// stride so MOPAC-touched frames coincide with the disk-emitted
    /// frames. Default 1 = MOPAC on every dispatched frame (the
    /// original FullFatFrameExtraction behaviour).
    std::size_t mopac_stride = 1;
    std::optional<FramePdbEmission> emit_pdbs;
    std::optional<FrameNpyEmission> emit_npys;
};

/// @brief Discriminated union over the five supported modes.
///
/// Use @c std::visit at the dispatch site:
/// @code
/// std::visit([&](auto&& mode) { Run(mode, common, session); }, spec);
/// @endcode
using ModeSpec = std::variant<PdbMode,
                              ProtonatedPdbMode,
                              OrcaMode,
                              MutantMode,
                              TrajectoryMode>;

}  // namespace nmr::cli

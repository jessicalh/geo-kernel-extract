#pragma once
/// @file
/// Typed per-mode input specifications and the @ref ModeSpec variant.
///
/// Each supported nmr_extract mode is a distinct struct holding only
/// the fields that mode consumes. A @c std::variant binds them into
/// the @ref ModeSpec dispatch type. The runner uses @c std::visit
/// rather than a JobMode-enum switch.

#include "OrcaRunLoader.h"

#include <filesystem>
#include <variant>

namespace nmr::cli {

/// @brief Load a bare PDB, protonate with @c reduce, run the standard pipeline.
///
/// Heavy-atom-only PDB on input. Hydrogens are added at runtime via the
/// @c reduce tool, charges from ff14SB, then the full per-conformation
/// pipeline runs and emits per-atom NPYs.
struct PdbMode {
    std::filesystem::path pdb;            ///< Path to the bare PDB.
    double                pH    = 7.0;    ///< pH passed to @c reduce.
    bool                  mopac = true;   ///< Run PM7+MOZYME + derived calcs.
};

/// @brief Load a PDB that already has H atoms, run the standard pipeline.
///
/// Skips the @c reduce step; protonation variant is detected from the
/// explicit H atoms present.
struct ProtonatedPdbMode {
    std::filesystem::path pdb;
    bool                  mopac = true;
};

/// @brief Load a single tleap/AMBER-prepared pose, run the standard pipeline.
///
/// The pose is described by a root name: @c {root}.xyz (coordinates),
/// @c {root}.prmtop (topology + charges), and @c {root}_nmr.out path for
/// OrcaShieldingResult.
struct OrcaMode {
    OrcaRunFiles files;
    bool         mopac = true;
};

/// @brief Load a WT+ALA mutant pair, run on both, compute Δ shielding.
///
/// Each side is a pose with the OrcaMode shape. MutationDeltaResult
/// attaches to the WT conformation and stores per-atom delta tensors.
struct MutantMode {
    OrcaRunFiles wt;
    OrcaRunFiles ala;
    bool         mopac = true;
};

/// @brief Read a GROMACS relaxation run, per-frame calculators + H5.
///
/// @c dir contains @c production.tpr, @c production.trr, @c production.edr.
/// The parent of @c dir is searched for @c topol.top by convention
/// (GROMACS readback). Into @c --output the run always writes the
/// trajectory H5, the per-protein topology sidecars, and per-frame NPYs
/// (@c output/npys/) + per-frame PDBs (@c output/pdbs/) — one set per
/// dispatched frame.
///
/// Frame cadence is governed by exactly one knob, @c stride (CLI
/// @c --stride): MOPAC (when enabled), NPY emission and PDB emission all
/// act on the dispatched frames, never on independent sub-strides. The
/// former @c --mopac-stride / @c --npy-stride / @c --pdb-stride and the
/// emit time-windows were removed 2026-05-31 (one of them silently
/// halved production runs).
///
/// MOPAC default is off here because PM7+MOZYME is too expensive for
/// fleet-scale trajectory runs.
struct TrajectoryMode {
    std::filesystem::path dir;
    bool                  mopac  = false;  ///< @c --mopac: FullFat (MOPAC on every dispatched frame).
    std::size_t           stride = 1;      ///< @c --stride N: process every N-th TRR frame. The one cadence knob.
};

/// @brief The three GROMACS production files a @c --trajectory DIR holds.
///
/// Single source of truth for the @c production.{tpr,trr,edr} basenames,
/// which otherwise appear by hand in validation and execution. No
/// discovery — the names are fixed by convention (see CLAUDE.md mode 5
/// and PATTERNS.md "No file discovery"); @ref FromProductionDir only
/// joins them onto @c dir.
struct TrajectoryInputFiles {
    std::filesystem::path tpr;  ///< @c {dir}/production.tpr
    std::filesystem::path trr;  ///< @c {dir}/production.trr
    std::filesystem::path edr;  ///< @c {dir}/production.edr

    static TrajectoryInputFiles FromProductionDir(
            const std::filesystem::path& dir) {
        return {dir / "production.tpr",
                dir / "production.trr",
                dir / "production.edr"};
    }
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

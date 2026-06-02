/// @file
/// nmr::cli::PrintUsage implementation.

#include "Cli/PrintUsage.h"

#include <cstdio>

namespace nmr::cli {

void PrintUsage(const char* program_name) {
    (void)std::fprintf(stderr,
                       "Usage: %s MODE [common options]\n"
                       "\n"
                       "Modes (exactly one required):\n"
                       "\n"
                       "  --pdb FILE [--pH 7.0] [--no-mopac]\n"
                       "      Load a bare PDB, protonate with `reduce`, apply ff14SB charges,\n"
                       "      run every calculator, emit per-atom NPYs.\n"
                       "\n"
                       "  --protonated-pdb FILE [--no-mopac]\n"
                       "      Load a PDB that already has H atoms, apply ff14SB charges,\n"
                       "      run every calculator, emit per-atom NPYs.\n"
                       "\n"
                       "  --orca --root NAME [--no-mopac]\n"
                       "      Load a tleap/AMBER-prepared pose. NAME expands to:\n"
                       "        {NAME}.xyz      (coordinates, required)\n"
                       "        {NAME}.prmtop   (topology + charges, required)\n"
                       "        {NAME}_nmr.out  (DFT shielding tensors, optional)\n"
                       "\n"
                       "  --mutant --wt NAME --ala NAME [--no-mopac]\n"
                       "      Load a WT+ALA mutant pair (each a --root-style pose).\n"
                       "      Runs every calculator on both; computes WT-ALA delta tensors.\n"
                       "\n"
                       "  --trajectory DIR [--stride N] [--mopac]\n"
                       "      Read a GROMACS production run. DIR must contain production.tpr,\n"
                       "      production.trr, production.edr. Into --output writes the trajectory\n"
                       "      H5, the topology sidecars, per-frame PDBs (output/pdbs/) and per-frame\n"
                       "      NPYs (output/npys/).\n"
                       "      --stride N (default 1): process every N-th TRR frame. MOPAC, NPY and\n"
                       "        PDB output all act on exactly those dispatched frames -- this is the\n"
                       "        only cadence knob.\n"
                       "      --mopac: switch PerFrameExtractionSet -> FullFatFrameExtraction (MOPAC\n"
                       "        on every dispatched frame). Default off.\n"
                       "\n"
                       "Common options:\n"
                       "  --output DIR     Output directory for NPY feature arrays + trajectory H5.\n"
                       "  --config FILE    TOML with calculator parameter overrides (optional).\n"
                       "  --aimnet2 FILE   AIMNet2 .jpt model path. Required for AIMNet2-derived\n"
                       "                   calculators (charges, embedding, charge-response gradient).\n"
                       "  --help, -h       Show this message.\n",
                       program_name);
}

}  // namespace nmr::cli

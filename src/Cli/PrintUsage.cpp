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
                       "  --trajectory DIR [--mopac [--mopac-stride N]]\n"
                       "                   [--emit-frame-pdbs DIR [--pdb-stride N] [--pdb-from-ps T0]\n"
                       "                                          [--pdb-to-ps T1] [--pdb-decorator TAG]]\n"
                       "                   [--emit-frame-npys DIR [--npy-stride N] [--npy-from-ps T0]\n"
                       "                                          [--npy-to-ps T1]]\n"
                       "      Read a GROMACS production run. DIR must contain production.tpr,\n"
                       "      production.trr, production.edr. Emits trajectory H5 always;\n"
                       "      per-frame PDB and NPY sidecars are independent opt-ins.\n"
                       "      --mopac switches from PerFrameExtractionSet to FullFatFrameExtraction.\n"
                       "      --mopac-stride N (default 1): run MOPAC only on frames where\n"
                       "        frame_idx %% N == 0; other dispatched frames skip MOPAC. Should\n"
                       "        match --npy-stride and --pdb-stride so heavy work coincides\n"
                       "        with disk-emitted frames.\n"
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

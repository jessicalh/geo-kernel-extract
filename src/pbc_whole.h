/*
 * pbc_whole.h — Make protein molecules whole after PBC wrapping.
 *
 * XTC trajectories from GROMACS store coordinates wrapped into the
 * periodic box.  Atoms of a single molecule can end up on opposite
 * sides of the box boundary, producing unphysical geometries when
 * the frame is extracted as a PDB.
 *
 * Our XTC files contain only the protein (not solvent/ions), so the
 * coordinate count doesn't match the full-system TPR.  We solve this
 * by loading the full topology and trimming it in-place to keep only
 * the first molecule block (the protein).  Then we call do_pbc_mtop()
 * — the canonical GROMACS function that walks the bond graph to make
 * every molecule whole.
 *
 * This handles any box geometry (cubic, rhombic dodecahedron, triclinic)
 * and any molecular topology (branched sidechains, disulfide bridges,
 * multi-chain proteins) correctly.
 *
 * Requires linking against libgromacs.
 *
 * Usage:
 *   MoleculeWholer wholer(tpr_path);   // load once per protein
 *   wholer.make_whole(frame);           // fix each frame in-place
 */

#pragma once

#include "xtc_reader.h"

#include <filesystem>
#include <stdexcept>
#include <string>

#include <gromacs/fileio/tpxio.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/topology/topology.h>
#include <gromacs/utility/real.h>
#include <gromacs/utility/vectypes.h>

// XTC stores float coordinates.  Our reinterpret_cast to rvec* is only
// valid when GROMACS real == float (mixed-precision build, GMX_DOUBLE=0).
static_assert(std::is_same_v<real, float>,
    "pbc_whole.h requires single-precision GROMACS (GMX_DOUBLE=0)");

namespace fs = std::filesystem;

class MoleculeWholer {
public:
    MoleculeWholer(const MoleculeWholer&) = delete;
    MoleculeWholer& operator=(const MoleculeWholer&) = delete;
    explicit MoleculeWholer(const fs::path& tpr_path) {
        if (!fs::exists(tpr_path))
            throw std::runtime_error("TPR not found: " + tpr_path.string());

        // Load full system topology.
        int full_natoms = 0;
        matrix full_box;
        pbc_type_ = read_tpx(tpr_path, nullptr, full_box, &full_natoms,
                             nullptr, nullptr, &mtop_);

        if (mtop_.molblock.empty())
            throw std::runtime_error("TPR has no molecule blocks: " +
                                     tpr_path.string());

        // Identify protein molblocks: keep all leading blocks whose
        // moltype name starts with "Protein" (GROMACS pdb2gmx convention).
        // Multi-chain proteins get split into Protein_chain_A, Protein_chain_A2,
        // etc. — each is a separate molblock but the XTC contains all of them.
        // Non-protein blocks (SOL, NA, CL) are always after all protein blocks.
        size_t n_protein_blocks = 0;
        protein_natoms_ = 0;
        for (size_t b = 0; b < mtop_.molblock.size(); ++b) {
            int mt = mtop_.molblock[b].type;
            const char* name = *mtop_.moltype[mt].name;
            // GROMACS names protein moltypes "Protein_chain_X" or "Protein".
            // Non-protein moltypes are SOL, NA, CL, etc.
            if (std::string(name).rfind("Protein", 0) != 0)
                break;
            if (mtop_.molblock[b].nmol != 1)
                throw std::runtime_error(
                    "Multi-molecule protein block not supported (nmol=" +
                    std::to_string(mtop_.molblock[b].nmol) + "): " +
                    tpr_path.string());
            protein_natoms_ += mtop_.moltype[mt].atoms.nr;
            ++n_protein_blocks;
        }

        if (n_protein_blocks == 0)
            throw std::runtime_error("No protein molblocks found in: " +
                                     tpr_path.string());

        // Trim to protein blocks only so do_pbc_mtop iterates over
        // exactly the protein-only coordinate array.
        mtop_.molblock.resize(n_protein_blocks);
        mtop_.natoms = protein_natoms_;
        mtop_.finalize();
    }

    // Number of atoms in the protein molecule.
    int natoms() const { return protein_natoms_; }

    // Make the protein molecule whole in-place.
    // Uses do_pbc_mtop with the protein-only topology and the frame's
    // own box vectors (correct for NPT where box size varies).
    void make_whole(XtcFrame& frame) const {
        if (frame.natoms != protein_natoms_)
            throw std::runtime_error(
                "Atom count mismatch: protein topology=" +
                std::to_string(protein_natoms_) +
                " XTC frame=" + std::to_string(frame.natoms));

        auto* coords = reinterpret_cast<rvec*>(frame.x.data());

        // Use the frame's box (NPT simulations change box each step).
        matrix frame_box;
        for (int i = 0; i < DIM; ++i)
            for (int j = 0; j < DIM; ++j)
                frame_box[i][j] = frame.box[i][j];

        do_pbc_mtop(pbc_type_, frame_box, &mtop_, coords);
    }

    // Make whole from a raw coordinate buffer + box (dihedral pipeline).
    void make_whole(std::vector<float>& coords, const float box_in[3][3]) const {
        if (static_cast<int>(coords.size()) != protein_natoms_ * 3)
            throw std::runtime_error(
                "Coord size mismatch: expected " +
                std::to_string(protein_natoms_ * 3) +
                " got " + std::to_string(coords.size()));

        auto* rvecs = reinterpret_cast<rvec*>(coords.data());

        matrix frame_box;
        for (int i = 0; i < DIM; ++i)
            for (int j = 0; j < DIM; ++j)
                frame_box[i][j] = box_in[i][j];

        do_pbc_mtop(pbc_type_, frame_box, &mtop_, rvecs);
    }

private:
    gmx_mtop_t mtop_;           // Full topology, trimmed to protein block only
    PbcType    pbc_type_;
    int        protein_natoms_ = 0;
};

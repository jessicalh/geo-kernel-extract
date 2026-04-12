#include "FullSystemReader.h"
#include "OperationLog.h"

// GROMACS TPR reading
#include "gromacs/fileio/tpxio.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/inputrec.h"

#include <cstring>
#include <string>

namespace nmr {

// Identify water residues by name.
static bool IsWater(const char* name) {
    return strcmp(name, "SOL") == 0
        || strcmp(name, "TIP3") == 0
        || strcmp(name, "TIP3P") == 0
        || strcmp(name, "HOH") == 0
        || strcmp(name, "WAT") == 0;
}

// Identify common ions.
static bool IsIon(const char* name) {
    return strcmp(name, "NA") == 0
        || strcmp(name, "CL") == 0
        || strcmp(name, "K") == 0
        || strcmp(name, "Na+") == 0
        || strcmp(name, "Cl-") == 0
        || strcmp(name, "MG") == 0
        || strcmp(name, "CA") == 0
        || strcmp(name, "ZN") == 0
        || strcmp(name, "ION") == 0;
}

// Identify protein molblocks: name starts with "Protein" (GROMACS convention)
// or first residue is a known amino acid.
static bool IsProteinBlock(const gmx_moltype_t& mt) {
    const char* name = *mt.name;
    if (strncmp(name, "Protein", 7) == 0) return true;
    // Fallback: check first residue name
    if (mt.atoms.nres > 0) {
        const char* res = *(mt.atoms.resinfo[0].name);
        // Known amino acids (CHARMM naming)
        static const char* aa[] = {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
            "HSE","HSD","HSP","HIE","HID","HIP","ILE","LEU","LYS",
            "MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
            "NALA","NARG","CGLU","CLYS", nullptr
        };
        for (int i = 0; aa[i]; ++i)
            if (strcmp(res, aa[i]) == 0) return true;
    }
    return false;
}


bool FullSystemReader::ReadTopology(const std::string& tpr_path) {
    OperationLog::Scope scope("FullSystemReader::ReadTopology", tpr_path);

    TpxFileHeader tpx = readTpxHeader(tpr_path.c_str(), true);
    t_inputrec ir;
    t_state state;
    gmx_mtop_t mtop;
    read_tpx_state(tpr_path.c_str(), tpx.bIr ? &ir : nullptr,
                   &state, tpx.bTop ? &mtop : nullptr);

    topo_.total_atoms = mtop.natoms;

    // Walk molblocks to identify protein, water, ion atom ranges.
    size_t global_idx = 0;
    bool found_protein = false;

    for (const auto& molblock : mtop.molblock) {
        const auto& mt = mtop.moltype[molblock.type];
        size_t atoms_per_mol = mt.atoms.nr;
        size_t total_in_block = atoms_per_mol * molblock.nmol;

        if (!found_protein && IsProteinBlock(mt)) {
            // Protein block — assume one copy (nmol=1 for protein)
            topo_.protein_start = global_idx;
            topo_.protein_count = atoms_per_mol;
            found_protein = true;

        } else if (mt.atoms.nres > 0 &&
                   IsWater(*(mt.atoms.resinfo[0].name))) {
            // Water block
            if (topo_.water_count == 0)
                topo_.water_O_start = global_idx;

            topo_.water_count += molblock.nmol;

            // Extract charges from first water molecule
            if (topo_.water_O_charge == 0.0 && atoms_per_mol >= 3) {
                topo_.water_O_charge = mt.atoms.atom[0].q;
                topo_.water_H_charge = mt.atoms.atom[1].q;
            }

        } else if (mt.atoms.nres > 0 &&
                   IsIon(*(mt.atoms.resinfo[0].name))) {
            // Ion block
            if (topo_.ion_count == 0)
                topo_.ion_start = global_idx;

            for (int m = 0; m < molblock.nmol; ++m) {
                topo_.ion_charges.push_back(mt.atoms.atom[0].q);
                topo_.ion_atomic_numbers.push_back(mt.atoms.atom[0].atomnumber);
            }
            topo_.ion_count += molblock.nmol;
        }

        global_idx += total_in_block;
    }

    if (!found_protein) {
        error_ = "no protein block found in " + tpr_path;
        return false;
    }

    OperationLog::Info(LogCalcOther, "FullSystemReader",
        "total=" + std::to_string(topo_.total_atoms) +
        " protein=" + std::to_string(topo_.protein_count) +
        " water=" + std::to_string(topo_.water_count) +
        " ions=" + std::to_string(topo_.ion_count) +
        " water_O_q=" + std::to_string(topo_.water_O_charge) +
        " water_H_q=" + std::to_string(topo_.water_H_charge));

    return true;
}


bool FullSystemReader::ExtractFrame(
        const std::vector<float>& xyz,
        std::vector<Vec3>& protein_positions,
        SolventEnvironment& solvent) const {

    size_t n_total = xyz.size() / 3;
    if (n_total < topo_.total_atoms) {
        return false;  // frame doesn't have enough atoms
    }

    // nm → Angstrom conversion
    auto pos = [&](size_t idx) -> Vec3 {
        return Vec3(
            static_cast<double>(xyz[idx * 3 + 0]) * 10.0,
            static_cast<double>(xyz[idx * 3 + 1]) * 10.0,
            static_cast<double>(xyz[idx * 3 + 2]) * 10.0);
    };

    // 1. Protein positions
    protein_positions.resize(topo_.protein_count);
    for (size_t i = 0; i < topo_.protein_count; ++i)
        protein_positions[i] = pos(topo_.protein_start + i);

    // 2. Water molecules (groups of 3: O, H1, H2)
    solvent.waters.resize(topo_.water_count);
    solvent.water_O_positions.resize(topo_.water_count);
    for (size_t w = 0; w < topo_.water_count; ++w) {
        size_t base = topo_.water_O_start + w * 3;
        auto& mol = solvent.waters[w];
        mol.O_pos     = pos(base);
        mol.H1_pos    = pos(base + 1);
        mol.H2_pos    = pos(base + 2);
        mol.O_charge  = topo_.water_O_charge;
        mol.H_charge  = topo_.water_H_charge;
        solvent.water_O_positions[w] = mol.O_pos;
    }

    // 3. Ions (one atom each)
    solvent.ions.resize(topo_.ion_count);
    for (size_t i = 0; i < topo_.ion_count; ++i) {
        size_t idx = topo_.ion_start + i;
        solvent.ions[i].pos = pos(idx);
        solvent.ions[i].charge = topo_.ion_charges[i];
        solvent.ions[i].atomic_number = topo_.ion_atomic_numbers[i];
    }

    return true;
}

}  // namespace nmr

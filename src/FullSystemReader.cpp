#include "FullSystemReader.h"
#include "BuildResult.h"
#include "Protein.h"
#include "AminoAcidType.h"
#include "NamingRegistry.h"
#include "ChargeSource.h"
#include "PhysicalConstants.h"
#include "OperationLog.h"

// GROMACS TPR reading
#include "gromacs/fileio/tpxio.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/inputrec.h"

#include <cstring>
#include <set>
#include <string>

namespace nmr {

// ── Opaque stored parse ─────────────────────────────────────────
// Holds the parsed gmx_mtop_t so downstream methods (BondedParams,
// BuildProtein) use the same parse. Defined here, not in the header,
// so GROMACS includes stay in the .cpp.

namespace detail {
struct TprData {
    gmx_mtop_t mtop;
};
}  // namespace detail


FullSystemReader::FullSystemReader() = default;
FullSystemReader::~FullSystemReader() = default;


// ── PDB LOADING BOUNDARY helpers ────────────────────────────────
// Translation from GROMACS types to our typed objects. These operate
// on GROMACS t_atoms — valid only during the TPR parse, never at
// runtime. Same functions exist in GromacsEnsembleLoader.cpp for
// the fleet path.

static Element ElementFromAtomicNumber(int an) {
    switch (an) {
        case 1:  return Element::H;
        case 6:  return Element::C;
        case 7:  return Element::N;
        case 8:  return Element::O;
        case 16: return Element::S;
        default: return Element::Unknown;
    }
}

static const std::set<std::string> KNOWN_AMINO_ACIDS = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "HID","HIE","HIP","HSD","HSE","HSP","ASH","GLH","ASPP","GLUP",
    "CYX","CYM","CYS2","LYN","ARN","TYM","MSE"
};

// Detect HIS protonation variant by checking which H atoms are present.
static int DetectHisVariantFromAtoms(const t_atoms& atoms,
                                     int res_start, int res_end) {
    bool has_HD1 = false, has_HE2 = false;
    for (int ai = res_start; ai < res_end; ++ai) {
        std::string name = *(atoms.atomname[ai]);
        if (name == "HD1") has_HD1 = true;
        if (name == "HE2") has_HE2 = true;
    }
    if (has_HD1 && has_HE2) return 2;  // HIP/HSP
    if (has_HD1)            return 0;  // HID/HSD
    if (has_HE2)            return 1;  // HIE/HSE
    return -1;
}

// CHARMM residue name → protonation variant index.
static int VariantFromCharmmResidueName(const std::string& charmm_name,
                                        AminoAcid type) {
    if (type == AminoAcid::HIS) {
        if (charmm_name == "HSD" || charmm_name == "HID") return 0;
        if (charmm_name == "HSE" || charmm_name == "HIE") return 1;
        if (charmm_name == "HSP" || charmm_name == "HIP") return 2;
        return -1;
    }
    if (type == AminoAcid::ASP) {
        if (charmm_name == "ASH" || charmm_name == "ASPP") return 1;
        return 0;
    }
    if (type == AminoAcid::GLU) {
        if (charmm_name == "GLH" || charmm_name == "GLUP") return 1;
        return 0;
    }
    if (type == AminoAcid::CYS) {
        if (charmm_name == "CYX" || charmm_name == "CYS2") return 1;
        if (charmm_name == "CYM") return 2;
        return 0;
    }
    if (type == AminoAcid::LYS) {
        if (charmm_name == "LYN") return 1;
        return 0;
    }
    return 0;
}


// ── Identify molecule block types ───────────────────────────────

static bool IsWater(const char* name) {
    return strcmp(name, "SOL") == 0
        || strcmp(name, "TIP3") == 0
        || strcmp(name, "TIP3P") == 0
        || strcmp(name, "HOH") == 0
        || strcmp(name, "WAT") == 0;
}

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

static bool IsProteinBlock(const gmx_moltype_t& mt) {
    const char* name = *mt.name;
    if (strncmp(name, "Protein", 7) == 0) return true;
    if (mt.atoms.nres > 0) {
        const char* res = *(mt.atoms.resinfo[0].name);
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


// ── ReadTopology ────────────────────────────────────────────────
//
// Single TPR parse. Populates:
//   1. SystemTopology (atom ranges, water/ion charges)
//   2. BondedParameters (interaction lists, CMAP grids)
//   3. Stored gmx_mtop_t for BuildProtein()

bool FullSystemReader::ReadTopology(const std::string& tpr_path) {
    OperationLog::Scope scope("FullSystemReader::ReadTopology", tpr_path);

    tpr_ = std::make_unique<detail::TprData>();

    TpxFileHeader tpx = readTpxHeader(tpr_path.c_str(), true);
    t_inputrec ir;
    t_state state;
    read_tpx_state(tpr_path.c_str(), tpx.bIr ? &ir : nullptr,
                   &state, tpx.bTop ? &tpr_->mtop : nullptr);

    const auto& mtop = tpr_->mtop;
    topo_.total_atoms = mtop.natoms;

    // ── Pass 1: atom ranges (protein, water, ions) ──────────────
    size_t global_idx = 0;
    bool found_protein = false;

    for (const auto& molblock : mtop.molblock) {
        const auto& mt = mtop.moltype[molblock.type];
        size_t atoms_per_mol = mt.atoms.nr;
        size_t total_in_block = atoms_per_mol * molblock.nmol;

        if (!found_protein && IsProteinBlock(mt)) {
            topo_.protein_start = global_idx;
            topo_.protein_count = atoms_per_mol;
            found_protein = true;

        } else if (mt.atoms.nres > 0 &&
                   IsWater(*(mt.atoms.resinfo[0].name))) {
            if (topo_.water_count == 0)
                topo_.water_O_start = global_idx;
            topo_.water_count += molblock.nmol;
            if (topo_.water_O_charge == 0.0 && atoms_per_mol >= 3) {
                topo_.water_O_charge = mt.atoms.atom[0].q;
                topo_.water_H_charge = mt.atoms.atom[1].q;
            }

        } else if (mt.atoms.nres > 0 &&
                   IsIon(*(mt.atoms.resinfo[0].name))) {
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

    // ── Pass 2: bonded interaction parameters ───────────────────
    // Same molblock walk, extracting force field interaction lists
    // for the protein atoms. Uses the stored mtop — no re-read.

    const auto& ffp = mtop.ffparams;

    bonded_params_.cmap_grid_spacing = ffp.cmap_grid.grid_spacing;
    bonded_params_.cmap_grids.clear();
    for (const auto& cm : ffp.cmap_grid.cmapdata) {
        bonded_params_.cmap_grids.push_back(
            std::vector<double>(cm.cmap.begin(), cm.cmap.end()));
    }

    size_t global_offset = 0;
    for (const auto& molblock : mtop.molblock) {
        const auto& mt = mtop.moltype[molblock.type];
        size_t atoms_per_mol = mt.atoms.nr;

        if (!IsProteinBlock(mt)) {
            global_offset += atoms_per_mol * molblock.nmol;
            continue;
        }

        size_t protein_local_offset = global_offset - topo_.protein_start;

        auto prot_idx = [&](int local) -> size_t {
            return protein_local_offset + static_cast<size_t>(local);
        };

        const auto& ilist = mt.ilist;

        // Bonds (harmonic)
        {
            const auto& il = ilist[InteractionFunction::Bonds];
            const auto& ia = il.iatoms;
            int nra = interaction_function[InteractionFunction::Bonds].nratoms;
            for (size_t j = 0; j < ia.size(); j += nra + 1) {
                int type = ia[j];
                const auto& p = ffp.iparams[type].harmonic;
                BondedInteraction bi;
                bi.type = BondedInteraction::Bond;
                bi.n_atoms = 2;
                bi.atoms[0] = prot_idx(ia[j+1]);
                bi.atoms[1] = prot_idx(ia[j+2]);
                bi.p[0] = p.rA;
                bi.p[1] = p.krA;
                bonded_params_.interactions.push_back(bi);
            }
        }

        // Angles (harmonic) — GROMACS stores theta0 in DEGREES
        {
            const auto& il = ilist[InteractionFunction::Angles];
            const auto& ia = il.iatoms;
            int nra = interaction_function[InteractionFunction::Angles].nratoms;
            for (size_t j = 0; j < ia.size(); j += nra + 1) {
                int type = ia[j];
                const auto& p = ffp.iparams[type].harmonic;
                BondedInteraction bi;
                bi.type = BondedInteraction::Angle;
                bi.n_atoms = 3;
                bi.atoms[0] = prot_idx(ia[j+1]);
                bi.atoms[1] = prot_idx(ia[j+2]);
                bi.atoms[2] = prot_idx(ia[j+3]);
                bi.p[0] = p.rA * DEG_TO_RAD;
                bi.p[1] = p.krA;
                bonded_params_.interactions.push_back(bi);
            }
        }

        // Urey-Bradley
        {
            const auto& il = ilist[InteractionFunction::UreyBradleyPotential];
            const auto& ia = il.iatoms;
            int nra = interaction_function[InteractionFunction::UreyBradleyPotential].nratoms;
            for (size_t j = 0; j < ia.size(); j += nra + 1) {
                int type = ia[j];
                const auto& p = ffp.iparams[type].u_b;
                BondedInteraction bi_angle;
                bi_angle.type = BondedInteraction::Angle;
                bi_angle.n_atoms = 3;
                bi_angle.atoms[0] = prot_idx(ia[j+1]);
                bi_angle.atoms[1] = prot_idx(ia[j+2]);
                bi_angle.atoms[2] = prot_idx(ia[j+3]);
                bi_angle.p[0] = p.thetaA * DEG_TO_RAD;
                bi_angle.p[1] = p.kthetaA;
                bonded_params_.interactions.push_back(bi_angle);
                if (p.kUBA != 0.0) {
                    BondedInteraction bi_ub;
                    bi_ub.type = BondedInteraction::UreyBradley;
                    bi_ub.n_atoms = 3;
                    bi_ub.atoms[0] = prot_idx(ia[j+1]);
                    bi_ub.atoms[1] = prot_idx(ia[j+2]);
                    bi_ub.atoms[2] = prot_idx(ia[j+3]);
                    bi_ub.p[0] = p.r13A;
                    bi_ub.p[1] = p.kUBA;
                    bonded_params_.interactions.push_back(bi_ub);
                }
            }
        }

        // Proper dihedrals (periodic)
        {
            const auto& il = ilist[InteractionFunction::ProperDihedrals];
            const auto& ia = il.iatoms;
            int nra = interaction_function[InteractionFunction::ProperDihedrals].nratoms;
            for (size_t j = 0; j < ia.size(); j += nra + 1) {
                int type = ia[j];
                const auto& p = ffp.iparams[type].pdihs;
                BondedInteraction bi;
                bi.type = BondedInteraction::ProperDih;
                bi.n_atoms = 4;
                bi.atoms[0] = prot_idx(ia[j+1]);
                bi.atoms[1] = prot_idx(ia[j+2]);
                bi.atoms[2] = prot_idx(ia[j+3]);
                bi.atoms[3] = prot_idx(ia[j+4]);
                bi.p[0] = p.phiA * DEG_TO_RAD;
                bi.p[1] = p.cpA;
                bi.p[2] = static_cast<double>(p.mult);
                bonded_params_.interactions.push_back(bi);
            }
        }

        // Improper dihedrals (harmonic)
        {
            const auto& il = ilist[InteractionFunction::ImproperDihedrals];
            const auto& ia = il.iatoms;
            int nra = interaction_function[InteractionFunction::ImproperDihedrals].nratoms;
            for (size_t j = 0; j < ia.size(); j += nra + 1) {
                int type = ia[j];
                const auto& p = ffp.iparams[type].harmonic;
                BondedInteraction bi;
                bi.type = BondedInteraction::ImproperDih;
                bi.n_atoms = 4;
                bi.atoms[0] = prot_idx(ia[j+1]);
                bi.atoms[1] = prot_idx(ia[j+2]);
                bi.atoms[2] = prot_idx(ia[j+3]);
                bi.atoms[3] = prot_idx(ia[j+4]);
                bi.p[0] = p.rA * DEG_TO_RAD;
                bi.p[1] = p.krA;
                bonded_params_.interactions.push_back(bi);
            }
        }

        // CMAP (dihedral energy correction)
        {
            const auto& il = ilist[InteractionFunction::DihedralEnergyCorrectionMap];
            const auto& ia = il.iatoms;
            int nra = interaction_function[InteractionFunction::DihedralEnergyCorrectionMap].nratoms;
            for (size_t j = 0; j < ia.size(); j += nra + 1) {
                int type = ia[j];
                int grid_idx = ffp.iparams[type].cmap.cmapA;
                BondedInteraction bi;
                bi.type = BondedInteraction::CMAP;
                bi.n_atoms = 5;
                bi.atoms[0] = prot_idx(ia[j+1]);
                bi.atoms[1] = prot_idx(ia[j+2]);
                bi.atoms[2] = prot_idx(ia[j+3]);
                bi.atoms[3] = prot_idx(ia[j+4]);
                bi.atoms[4] = prot_idx(ia[j+5]);
                bi.p[0] = static_cast<double>(grid_idx);
                bonded_params_.interactions.push_back(bi);
            }
        }

        global_offset += atoms_per_mol * molblock.nmol;
    }

    OperationLog::Info(LogCalcOther, "FullSystemReader",
        std::to_string(bonded_params_.interactions.size()) + " bonded interactions, " +
        std::to_string(bonded_params_.cmap_grids.size()) + " CMAP grids");

    return true;
}


// ── ExtractFrame ────────────────────────────────────────────────

bool FullSystemReader::ExtractFrame(
        const std::vector<float>& xyz,
        std::vector<Vec3>& protein_positions,
        SolventEnvironment& solvent) const {

    size_t n_total = xyz.size() / 3;
    if (n_total < topo_.total_atoms) {
        return false;
    }

    auto pos = [&](size_t idx) -> Vec3 {
        return Vec3(
            static_cast<double>(xyz[idx * 3 + 0]) * 10.0,
            static_cast<double>(xyz[idx * 3 + 1]) * 10.0,
            static_cast<double>(xyz[idx * 3 + 2]) * 10.0);
    };

    protein_positions.resize(topo_.protein_count);
    for (size_t i = 0; i < topo_.protein_count; ++i)
        protein_positions[i] = pos(topo_.protein_start + i);

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

    solvent.ions.resize(topo_.ion_count);
    for (size_t i = 0; i < topo_.ion_count; ++i) {
        size_t idx = topo_.ion_start + i;
        solvent.ions[i].pos = pos(idx);
        solvent.ions[i].charge = topo_.ion_charges[i];
        solvent.ions[i].atomic_number = topo_.ion_atomic_numbers[i];
    }

    return true;
}


// ── BuildProtein ────────────────────────────────────────────────
//
// Build a Protein + ChargeSource from the stored TPR parse.
// Same logic as BuildProteinFromTpr in GromacsEnsembleLoader.cpp
// but uses the gmx_mtop_t already parsed by ReadTopology().
// PDB LOADING BOUNDARY: translates CHARMM naming to typed objects.

BuildResult FullSystemReader::BuildProtein(
        const std::string& protein_id,
        ForceField force_field) const {

    BuildResult result;

    if (!tpr_) {
        result.error = "BuildProtein called before ReadTopology";
        return result;
    }

    const auto& mtop = tpr_->mtop;

    // Find the protein moltype
    const gmx_moltype_t* protein_moltype = nullptr;
    for (const auto& molblock : mtop.molblock) {
        const auto& moltype = mtop.moltype[molblock.type];
        if (moltype.atoms.nres > 0) {
            std::string first_res = *(moltype.atoms.resinfo[0].name);
            if (KNOWN_AMINO_ACIDS.count(first_res) > 0) {
                protein_moltype = &moltype;
                break;
            }
        }
    }

    if (!protein_moltype) {
        result.error = "no protein molecule found in stored TPR for " + protein_id;
        return result;
    }

    const t_atoms& tpr_atoms = protein_moltype->atoms;
    size_t n_protein_atoms = tpr_atoms.nr;
    size_t n_residues = tpr_atoms.nres;

    auto& registry = GlobalNamingRegistry();
    auto protein = std::make_unique<Protein>();

    // Residue ranges
    struct ResRange { int start; int end; };
    std::vector<ResRange> ranges(n_residues);
    {
        std::vector<int> first_atom(n_residues, -1);
        for (int ai = 0; ai < tpr_atoms.nr; ++ai) {
            int ri = tpr_atoms.atom[ai].resind;
            if (ri >= 0 && ri < static_cast<int>(n_residues) && first_atom[ri] < 0)
                first_atom[ri] = ai;
        }
        for (size_t ri = 0; ri < n_residues; ++ri) {
            ranges[ri].start = first_atom[ri];
            ranges[ri].end = (ri + 1 < n_residues) ? first_atom[ri + 1]
                                                     : tpr_atoms.nr;
        }
    }

    // Add residues
    for (size_t ri = 0; ri < n_residues; ++ri) {
        Residue res;
        std::string charmm_name = *(tpr_atoms.resinfo[ri].name);
        std::string canonical = registry.ToCanonical(charmm_name);
        if (canonical.empty()) {
            result.error = "unknown residue '" + charmm_name +
                           "' at position " + std::to_string(ri) +
                           " in " + protein_id;
            return result;
        }
        res.type = AminoAcidFromThreeLetterCode(canonical);
        if (res.type == AminoAcid::Unknown) {
            result.error = "unrecognised amino acid '" + canonical +
                           "' from '" + charmm_name + "' at " + std::to_string(ri);
            return result;
        }
        res.sequence_number = tpr_atoms.resinfo[ri].nr;
        res.chain_id = "A";
        res.protonation_variant_index =
            VariantFromCharmmResidueName(charmm_name, res.type);
        if (res.type == AminoAcid::HIS && res.protonation_variant_index < 0) {
            res.protonation_variant_index =
                DetectHisVariantFromAtoms(tpr_atoms, ranges[ri].start, ranges[ri].end);
        }
        protein->AddResidue(std::move(res));
    }

    // Add atoms
    for (int ai = 0; ai < tpr_atoms.nr; ++ai) {
        auto atom = std::make_unique<Atom>();
        std::string charmm_atom_name = *(tpr_atoms.atomname[ai]);
        size_t ri = tpr_atoms.atom[ai].resind;
        std::string canonical_residue =
            ThreeLetterCodeForAminoAcid(protein->ResidueAt(ri).type);
        atom->pdb_atom_name = registry.TranslateAtomName(
            charmm_atom_name, canonical_residue,
            ToolContext::Charmm, ToolContext::Standard);
        atom->element = ElementFromAtomicNumber(tpr_atoms.atom[ai].atomnumber);
        atom->residue_index = ri;
        size_t idx = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(idx);
        const std::string& name = protein->AtomAt(idx).pdb_atom_name;
        Residue& res_ref = protein->MutableResidueAt(ri);
        if      (name == "N"  && res_ref.N  == Residue::NONE) res_ref.N  = idx;
        else if (name == "CA" && res_ref.CA == Residue::NONE) res_ref.CA = idx;
        else if (name == "C"  && res_ref.C  == Residue::NONE &&
                 protein->AtomAt(idx).element == Element::C &&
                 res_ref.CA != Residue::NONE) res_ref.C = idx;
        else if (name == "O"  && res_ref.O  == Residue::NONE) res_ref.O  = idx;
        else if ((name == "H" || name == "HN") &&
                 res_ref.H == Residue::NONE) res_ref.H = idx;
        else if ((name == "HA" || name == "HA2") &&
                 res_ref.HA == Residue::NONE) res_ref.HA = idx;
        else if (name == "CB" && res_ref.CB == Residue::NONE) res_ref.CB = idx;
    }

    // Build context
    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->pdb_source = protein_id;
    ctx->force_field = ForceFieldName(force_field);
    ctx->protonation_tool = "pdb2gmx";
    protein->SetBuildContext(std::move(ctx));

    // Extract charges from TPR
    std::vector<AtomChargeRadius> charge_vec(n_protein_atoms);
    double charge_sum = 0.0;
    for (size_t ai = 0; ai < n_protein_atoms; ++ai) {
        charge_vec[ai].charge = tpr_atoms.atom[ai].q;
        charge_vec[ai].radius = 1.5;
        charge_sum += tpr_atoms.atom[ai].q;
    }
    result.net_charge = static_cast<int>(
        charge_sum + (charge_sum > 0 ? 0.5 : -0.5));
    result.charges = std::make_unique<PreloadedChargeSource>(
        std::move(charge_vec), force_field);

    result.protein = std::move(protein);
    result.ok = true;

    OperationLog::Info(LogCalcOther, "FullSystemReader::BuildProtein",
        protein_id + ": " +
        std::to_string(n_protein_atoms) + " atoms, " +
        std::to_string(n_residues) + " residues");

    return result;
}

}  // namespace nmr

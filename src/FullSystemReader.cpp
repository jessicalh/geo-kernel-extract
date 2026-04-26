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
#include "gromacs/pbcutil/pbc.h"
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
// BuildProtein, MakeProteinWhole) use the same parse. Defined here,
// not in the header, so GROMACS includes stay in the .cpp.
//
// After ReadTopology completes, `mtop` is the protein-only topology:
// the original full-system mtop is trimmed in place at the end of
// the parse so only the leading protein molblocks remain and natoms
// equals topo_.protein_count. This single trimmed mtop serves
// BuildProtein, the bonded-params extraction (which already ran on
// the full one before the trim), and MakeProteinWhole (which passes
// it to do_pbc_mtop). gmx_mtop_t is not copyable, so an in-place
// trim is the only way to keep one parsed instance.
// `pbc_type` is captured from the TPR's t_inputrec at parse time.

namespace detail {
struct TprData {
    gmx_mtop_t mtop;
    PbcType    pbc_type = PbcType::Unset;
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


// ── Identify molecule block types by composition (not by name) ──
//
// The trajectory layer reads what GROMACS chose to write, including
// pdb2gmx's internal moltype names ("Protein", "Protein_chain_A1",
// "SOL", "NA", ...). Matching those names amounts to programming via
// strings — once a moltype name shifts (different water model
// renamed by a tool, force field convention drift, custom build
// pipeline) the rule is silently wrong. Composition is what we
// actually mean: water has a 3-atom O+H+H signature; an ion is a
// single-atom moltype; everything in the leading non-solvent run is
// the protein. The predicates below use only typed integer fields
// (atom counts, residue counts, atomic numbers) from gmx_moltype_t.

// Water moltype: 3 atoms, exactly one O and two H by atomic number,
// one residue. TIP3P / SPC / OPC all match (3-site models). TIP4P
// and 4-site models would need an explicit accommodation; we don't
// currently use them and would want to know if a TPR appeared with
// one rather than have it slip through.
static bool IsWaterMoltype(const gmx_moltype_t& mt) {
    if (mt.atoms.nr != 3 || mt.atoms.nres != 1) return false;
    int n_O = 0, n_H = 0, n_other = 0;
    for (int i = 0; i < 3; ++i) {
        const int z = mt.atoms.atom[i].atomnumber;
        if      (z == 8) ++n_O;
        else if (z == 1) ++n_H;
        else             ++n_other;
    }
    return n_O == 1 && n_H == 2 && n_other == 0;
}

// Ion moltype: a single atom in a single residue. Atomic number and
// charge come from the moltype's first (only) atom. We do not gate
// on which element it is — counter-ions in CHARMM36m / ff14SB
// fleets are Na+, Cl-, K+, Mg2+, Ca2+, Zn2+, but the typed-int rule
// covers any single-atom moltype that is not water-shaped.
static bool IsIonMoltype(const gmx_moltype_t& mt) {
    return mt.atoms.nr == 1 && mt.atoms.nres == 1;
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

    // Capture pbcType from the inputrec for later use by MakeProteinWhole.
    // Production TPRs always carry an inputrec (tpx.bIr). If a future
    // input lacks one, MakeProteinWhole will fall back to PbcType::Unset
    // and do_pbc_mtop will assume cubic — surfaced as an explicit error.
    if (tpx.bIr) {
        tpr_->pbc_type = ir.pbcType;
    }

    // Per-molblock diagnostic: emitted unconditionally so the GROMACS
    // moltype layout is visible in every run log. Useful both when
    // GROMACS splits one biological chain into multiple Protein_*
    // moltypes and when it doesn't — fleet-wide audit-by-grep.
    for (size_t b = 0; b < mtop.molblock.size(); ++b) {
        const auto& molblock = mtop.molblock[b];
        const auto& mt = mtop.moltype[molblock.type];
        OperationLog::Info(LogCalcOther, "FullSystemReader",
            "molblock[" + std::to_string(b) + "] name=" +
            std::string(*mt.name) +
            " nmol=" + std::to_string(molblock.nmol) +
            " atoms_per_mol=" + std::to_string(mt.atoms.nr));
    }

    // ── Pass 1: classify molblocks by composition ───────────────
    //
    // Convention: GROMACS lays the topology out in a fixed order —
    // protein molblock(s) first, then water, then ions. The protein
    // is potentially split across multiple Protein_* moltypes for
    // GROMACS-internal reasons (chain breaks, disulfide handling).
    // We accept whatever leading non-water non-ion run is present
    // as ONE biological chain (the project does not currently support
    // multi-chain biology in trajectory mode), and refuse a TPR that
    // shows a non-{water,ion} block after the solvent region — that
    // would be a cofactor / co-solute pattern we don't yet handle.
    size_t global_idx = 0;
    size_t n_protein_molblocks = 0;
    bool past_protein = false;

    for (const auto& molblock : mtop.molblock) {
        const auto& mt = mtop.moltype[molblock.type];
        const size_t atoms_per_mol = mt.atoms.nr;
        const size_t total_in_block = atoms_per_mol * molblock.nmol;

        const bool is_water = IsWaterMoltype(mt);
        const bool is_ion   = IsIonMoltype(mt);

        if (!past_protein && !is_water && !is_ion) {
            // Leading non-solvent block: protein (possibly chain-split).
            // Strict nmol=1 — pdb2gmx convention. Homo-multimers
            // (nmol>1) would need replicated residue/atom construction
            // in BuildProtein and we want to fail loudly if seen.
            if (molblock.nmol != 1) {
                error_ = "Protein-shape molblock '" +
                         std::string(*mt.name) +
                         "' has nmol=" + std::to_string(molblock.nmol) +
                         " (expected 1; pdb2gmx convention). " +
                         "Replicated protein moltypes are not supported.";
                return false;
            }
            if (n_protein_molblocks == 0) {
                topo_.protein_start = global_idx;
            }
            topo_.protein_count += atoms_per_mol;
            ++n_protein_molblocks;

        } else {
            // First non-protein block closes the protein region.
            if (n_protein_molblocks > 0) past_protein = true;

            if (is_water) {
                if (topo_.water_count == 0)
                    topo_.water_O_start = global_idx;
                topo_.water_count += molblock.nmol;
                if (topo_.water_O_charge == 0.0) {
                    topo_.water_O_charge = mt.atoms.atom[0].q;
                    topo_.water_H_charge = mt.atoms.atom[1].q;
                }
            } else if (is_ion) {
                if (topo_.ion_count == 0)
                    topo_.ion_start = global_idx;
                for (int m = 0; m < molblock.nmol; ++m) {
                    topo_.ion_charges.push_back(mt.atoms.atom[0].q);
                    topo_.ion_atomic_numbers.push_back(
                        mt.atoms.atom[0].atomnumber);
                }
                topo_.ion_count += molblock.nmol;
            } else {
                // Non-water non-ion AFTER the protein region:
                // cofactor, lipid, co-solute. Refuse explicitly so a
                // future case is surfaced rather than mis-sliced.
                error_ = "Unsupported moltype after protein/solvent region: '" +
                         std::string(*mt.name) +
                         "' (atoms=" + std::to_string(atoms_per_mol) +
                         " residues=" + std::to_string(mt.atoms.nres) +
                         " nmol=" + std::to_string(molblock.nmol) +
                         "). Cofactors and co-solutes are not currently " +
                         "supported by the trajectory loader.";
                return false;
            }
        }

        global_idx += total_in_block;
    }

    if (n_protein_molblocks == 0) {
        error_ = "no protein-shape molblock found in " + tpr_path;
        return false;
    }

    // Strong, distinct log line for the chain-split case so it is
    // unmistakable in fleet audit output.
    if (n_protein_molblocks > 1) {
        OperationLog::Warn("FullSystemReader",
            "PROTEIN SPLIT: " + std::to_string(n_protein_molblocks) +
            " leading non-solvent molblocks summed to " +
            std::to_string(topo_.protein_count) +
            " atoms (one biological chain, encoded by GROMACS as " +
            "multiple Protein_* moltypes).");
    }

    OperationLog::Info(LogCalcOther, "FullSystemReader",
        "total=" + std::to_string(topo_.total_atoms) +
        " protein=" + std::to_string(topo_.protein_count) +
        " (across " + std::to_string(n_protein_molblocks) +
        " molblock" + (n_protein_molblocks == 1 ? "" : "s") + ")" +
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

    // Iterate exactly the leading protein molblocks identified in
    // Pass 1. Each one's local atom indices are offset by the
    // running global_offset so the resulting BondedInteraction list
    // uses indices consistent with BuildProtein's concatenated atom
    // ordering.
    size_t global_offset = 0;
    for (size_t b = 0; b < n_protein_molblocks; ++b) {
        const auto& molblock = mtop.molblock[b];
        const auto& mt = mtop.moltype[molblock.type];
        size_t atoms_per_mol = mt.atoms.nr;

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

    // Trim the parsed mtop in place to the protein-only slice.
    // gmx_mtop_t is not copyable, so this is the single canonical
    // mtop from this point on. BuildProtein and MakeProteinWhole
    // both consume it; do_pbc_mtop walks only the protein bond
    // graph because solvent/ion molblocks are no longer present.
    tpr_->mtop.molblock.resize(n_protein_molblocks);
    tpr_->mtop.natoms = static_cast<int>(topo_.protein_count);
    tpr_->mtop.finalize();

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


// ── MakeProteinWhole ────────────────────────────────────────────
//
// Replaces the former MoleculeWholer dependency. do_pbc_mtop walks
// the protein-only mtop (built once at ReadTopology time) and makes
// each molecule whole within the box. Coordinates are float (XTC
// precision) and treated as rvec*; the static_assert in the GROMACS
// headers guarantees real == float in this build.

bool FullSystemReader::MakeProteinWhole(
        std::vector<float>& protein_coords,
        const float box_in[3][3]) const {
    if (!tpr_) {
        return false;
    }
    if (protein_coords.size() != topo_.protein_count * 3) {
        OperationLog::Error("FullSystemReader::MakeProteinWhole",
            "buffer size " + std::to_string(protein_coords.size()) +
            " != 3 * protein_count " +
            std::to_string(topo_.protein_count));
        return false;
    }
    if (tpr_->pbc_type == PbcType::Unset) {
        OperationLog::Error("FullSystemReader::MakeProteinWhole",
            "pbc_type unset (TPR did not carry an inputrec)");
        return false;
    }

    auto* rvecs = reinterpret_cast<rvec*>(protein_coords.data());
    matrix frame_box;
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            frame_box[i][j] = box_in[i][j];
    do_pbc_mtop(tpr_->pbc_type, frame_box, &tpr_->mtop, rvecs);
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

    // Collect leading non-water non-ion moltypes — the protein
    // chain[s] in molblock order. Composition predicates match the
    // ones in ReadTopology Pass 1; the slice they imply is identical
    // and consistent with Topology().protein_count by construction.
    std::vector<const gmx_moltype_t*> protein_moltypes;
    for (const auto& molblock : mtop.molblock) {
        const auto& moltype = mtop.moltype[molblock.type];
        if (IsWaterMoltype(moltype) || IsIonMoltype(moltype)) break;
        if (molblock.nmol != 1) {
            result.error = "Protein-shape molblock '" +
                           std::string(*moltype.name) +
                           "' has nmol=" + std::to_string(molblock.nmol) +
                           " in " + protein_id +
                           " (BuildProtein expects nmol=1).";
            return result;
        }
        protein_moltypes.push_back(&moltype);
    }

    if (protein_moltypes.empty()) {
        result.error = "no protein-shape moltype found in stored TPR for " +
                       protein_id;
        return result;
    }

    if (protein_moltypes.size() > 1) {
        OperationLog::Warn("FullSystemReader::BuildProtein",
            protein_id + ": building one biological chain from " +
            std::to_string(protein_moltypes.size()) +
            " GROMACS molblocks (chain split).");
    }

    auto& registry = GlobalNamingRegistry();
    auto protein = std::make_unique<Protein>();

    // Continuous indexing across protein moltypes. residue_offset
    // advances after each block's residues are added; atom_offset is
    // the running total of atoms so charges array indexing matches
    // the Protein's atom order.
    size_t residue_offset = 0;
    size_t atom_offset = 0;
    std::vector<AtomChargeRadius> charge_vec;
    double charge_sum = 0.0;

    for (const gmx_moltype_t* moltype : protein_moltypes) {
        const t_atoms& tpr_atoms = moltype->atoms;
        const size_t n_residues = tpr_atoms.nres;

        // Residue ranges within this moltype (atom indices local to
        // tpr_atoms — will be globalised when atoms are added below).
        struct ResRange { int start; int end; };
        std::vector<ResRange> ranges(n_residues);
        {
            std::vector<int> first_atom(n_residues, -1);
            for (int ai = 0; ai < tpr_atoms.nr; ++ai) {
                int ri = tpr_atoms.atom[ai].resind;
                if (ri >= 0 && ri < static_cast<int>(n_residues) &&
                    first_atom[ri] < 0)
                    first_atom[ri] = ai;
            }
            for (size_t ri = 0; ri < n_residues; ++ri) {
                ranges[ri].start = first_atom[ri];
                ranges[ri].end = (ri + 1 < n_residues) ? first_atom[ri + 1]
                                                          : tpr_atoms.nr;
            }
        }

        // Residues for this moltype.
        for (size_t ri = 0; ri < n_residues; ++ri) {
            Residue res;
            std::string charmm_name = *(tpr_atoms.resinfo[ri].name);
            std::string canonical = registry.ToCanonical(charmm_name);
            if (canonical.empty()) {
                result.error = "unknown residue '" + charmm_name +
                               "' at position " +
                               std::to_string(residue_offset + ri) +
                               " in " + protein_id;
                return result;
            }
            res.type = AminoAcidFromThreeLetterCode(canonical);
            if (res.type == AminoAcid::Unknown) {
                result.error = "unrecognised amino acid '" + canonical +
                               "' from '" + charmm_name + "' at " +
                               std::to_string(residue_offset + ri);
                return result;
            }
            res.sequence_number = tpr_atoms.resinfo[ri].nr;
            res.chain_id = "A";  // single biological chain (multi-molblock
                                 // is a GROMACS encoding artifact, not a
                                 // chain split at the biology layer).
            res.protonation_variant_index =
                VariantFromCharmmResidueName(charmm_name, res.type);
            if (res.type == AminoAcid::HIS &&
                res.protonation_variant_index < 0) {
                res.protonation_variant_index = DetectHisVariantFromAtoms(
                    tpr_atoms, ranges[ri].start, ranges[ri].end);
            }
            protein->AddResidue(std::move(res));
        }

        // Atoms for this moltype, with global residue index.
        for (int ai = 0; ai < tpr_atoms.nr; ++ai) {
            auto atom = std::make_unique<Atom>();
            std::string charmm_atom_name = *(tpr_atoms.atomname[ai]);
            const size_t local_ri = tpr_atoms.atom[ai].resind;
            const size_t global_ri = residue_offset + local_ri;
            std::string canonical_residue =
                ThreeLetterCodeForAminoAcid(protein->ResidueAt(global_ri).type);
            atom->iupac_name = registry.TranslateAtomName(
                charmm_atom_name, canonical_residue,
                ToolContext::Charmm, ToolContext::Standard);
            atom->element = ElementFromAtomicNumber(
                tpr_atoms.atom[ai].atomnumber);
            atom->residue_index = global_ri;
            const size_t idx = protein->AddAtom(std::move(atom));
            protein->MutableResidueAt(global_ri).atom_indices.push_back(idx);
            const IupacAtomName& name = protein->AtomAt(idx).iupac_name;
            Residue& res_ref = protein->MutableResidueAt(global_ri);
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

            AtomChargeRadius cr;
            cr.charge = tpr_atoms.atom[ai].q;
            cr.radius = 1.5;
            charge_vec.push_back(cr);
            charge_sum += tpr_atoms.atom[ai].q;
        }

        residue_offset += n_residues;
        atom_offset += static_cast<size_t>(tpr_atoms.nr);
    }

    // Build context
    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->pdb_source = protein_id;
    ctx->force_field = ForceFieldName(force_field);
    ctx->protonation_tool = "pdb2gmx";
    protein->SetBuildContext(std::move(ctx));

    result.net_charge = static_cast<int>(
        charge_sum + (charge_sum > 0 ? 0.5 : -0.5));
    result.charges = std::make_unique<PreloadedChargeSource>(
        std::move(charge_vec), force_field);

    result.protein = std::move(protein);
    result.ok = true;

    OperationLog::Info(LogCalcOther, "FullSystemReader::BuildProtein",
        protein_id + ": " +
        std::to_string(atom_offset) + " atoms, " +
        std::to_string(residue_offset) + " residues" +
        (protein_moltypes.size() > 1 ?
            " across " + std::to_string(protein_moltypes.size()) +
            " molblocks" : ""));

    return result;
}

}  // namespace nmr

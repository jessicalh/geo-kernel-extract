// tools/topology/build_semantic_tables.cpp
//
// Build-time generator for the LegacyAmberTopology semantic tables.
//
// String barrier note: this file is the STRING BOUNDARY in its
// strongest form. Strings flow IN from CCD (cifpp) and from RDKit's
// perception output. They are reconciled here against the typed-enum
// vocabulary in src/SemanticEnums.h and emitted as typed-enum
// literals into src/generated/LegacyAmberSemanticTables.cpp. After
// emission, no string survives. The runtime library libnmr_shielding
// does not link RDKit and never sees this code.
//
// Process log: this generator emits a structured log (NDJSON) of
// every decision it makes -- per (residue, atom, field), what each
// source said, what was chosen, why. The log is committed alongside
// the generated C++ source so the audit trail lives with the data.
// Inspection: the user reads the log once after each generation run;
// committed logs serve as a reproducibility record.

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// cifpp -- CCD reader at the loading boundary.
#include <cif++.hpp>

// RDKit -- chemistry perception.
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/new_canon.h>
#include <CIPLabeler/CIPLabeler.h>
#include <RDGeneral/versions.h>

// Types.h + SemanticEnums.h -- the typed-enum vocabulary the runtime sees.
#include "Types.h"
#include "SemanticEnums.h"

// Runtime substrate header -- the home of `ParseAtomName`,
// `LocantLetterToEnum`, `ElementFromAtomName`, `ParsedAtomName`. Lifted
// out of this generator into the runtime header so the generator and
// the runtime share a SINGLE copy. This file used to carry private
// duplicates; they are removed and the generator now calls the lifted
// versions inline. See spec/plan/bones/topology-and-identity-pivot-synthesis-2026-05-05.md
// §13.1 for the rationale.
#include "generated/LegacyAmberSemanticTables.h"

using namespace nmr;

// Bring the lifted parser into the generator's local scope so the
// existing call sites (SynthesisedFor*, BuildAtomSemanticEntry, the
// per-residue loop in BuildAtomSemanticEntries) compile unchanged.
// The generator's local `ParsedName` is now an alias for the runtime
// header's `ParsedAtomName`; struct shape is identical (same fields,
// same defaults, same semantics) so existing field accesses keep working.
using ParsedName = ::nmr::topology_generated::ParsedAtomName;
using ::nmr::topology_generated::ElementFromAtomName;
using ::nmr::topology_generated::LocantLetterToEnum;
using ::nmr::topology_generated::ParseAtomName;


namespace {

struct Args {
    std::string ccd_path;
    std::string output_cpp_path;
    std::string log_path;
    bool        help        = false;
};

void PrintUsage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s --ccd PATH --output PATH --log PATH\n"
        "\n"
        "Build the LegacyAmberTopology semantic tables.\n"
        "\n"
        "Required:\n"
        "  --ccd PATH      Path to wwPDB Chemical Component Dictionary\n"
        "                  (typically data/ccd/components.cif).\n"
        "  --output PATH   Path to write the generated C++ source\n"
        "                  (typically src/generated/LegacyAmberSemanticTables.cpp).\n"
        "  --log PATH      Path to write the structured process log\n"
        "                  (typically src/generated/LegacyAmberSemanticTables.log.txt).\n"
        "\n"
        "Optional:\n"
        "  --help          Show this message.\n",
        prog);
}

bool ParseArgs(int argc, char** argv, Args& args) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--help" || a == "-h") {
            args.help = true;
            return true;
        }
        if (i + 1 >= argc) {
            std::fprintf(stderr, "ERROR: %s requires a value\n", a.c_str());
            return false;
        }
        const std::string v = argv[++i];
        if      (a == "--ccd")    args.ccd_path        = v;
        else if (a == "--output") args.output_cpp_path = v;
        else if (a == "--log")    args.log_path        = v;
        else {
            std::fprintf(stderr, "ERROR: unknown argument %s\n", a.c_str());
            return false;
        }
    }
    if (args.ccd_path.empty() || args.output_cpp_path.empty() || args.log_path.empty()) {
        std::fprintf(stderr, "ERROR: --ccd, --output, --log are all required\n");
        return false;
    }
    return true;
}


// Process log: structured plain-text record. Each line is either a
// section header (==== SECTION NAME ====) or a key=value record. The
// format is grep-able and human-readable; machine consumers can
// parse with awk or a simple state machine if needed.
//
// Open once at start of run, write each event, close at end. The
// log is committed alongside the generated C++ source.
class ProcessLog {
public:
    explicit ProcessLog(const std::string& path) : out_(path) {
        if (!out_.is_open()) {
            std::fprintf(stderr, "ERROR: cannot open log file %s\n", path.c_str());
            std::exit(2);
        }
        out_ << "# LegacyAmberTopology semantic tables -- generator process log\n"
                "# Format: section headers (==== SECTION ====), key=value records,\n"
                "# blank lines for grouping. grep-able by section / key.\n"
                "#\n"
                "# This log is committed alongside src/generated/LegacyAmberSemanticTables.cpp\n"
                "# as the reproducibility audit trail. Inspect once per regeneration.\n"
                "#\n";
    }

    void Section(const std::string& name) {
        out_ << "\n==== " << name << " ====\n";
    }
    void KV(const std::string& key, const std::string& value) {
        out_ << key << " = " << value << "\n";
    }
    void KV(const std::string& key, int value) {
        out_ << key << " = " << value << "\n";
    }
    void Note(const std::string& msg) {
        out_ << "# " << msg << "\n";
    }

private:
    std::ofstream out_;
};


// ============================================================================
// CCD entry: typed-string record of one residue's atoms + bonds.
//
// These are the raw strings flowing IN from CCD via cifpp. They live
// only inside this generator binary; they are never serialised, never
// returned across an API boundary, and never reach the runtime
// library. The reconciler reads them, runs RDKit perception, and
// writes typed-enum values into the emitted C++ source.
// ============================================================================

struct CcdAtom {
    std::string atom_id;             ///< e.g. "CG", "HD11"
    std::string alt_atom_id;         ///< legacy/alternate name from CCD
    std::string type_symbol;         ///< element symbol "C", "N", ...
    int         formal_charge = 0;
    bool        aromatic_flag = false;          ///< pdbx_aromatic_flag Y/N
    bool        leaving_atom_flag = false;      ///< pdbx_leaving_atom_flag Y/N
    std::string stereo_config;       ///< pdbx_stereo_config "R"/"S"/"N"
    int         ordinal = 0;         ///< pdbx_ordinal -- 1-based atom index
};

struct CcdBond {
    std::string atom_id_1;           ///< first atom in the bond
    std::string atom_id_2;           ///< second atom in the bond
    std::string value_order;         ///< CCD value_order: SING/DOUB/TRIP/AROM/DELO/...
    bool        aromatic_flag = false;
    std::string stereo_config;       ///< pdbx_stereo_config "E"/"Z"/"N"
    int         ordinal = 0;
};

struct CcdResidue {
    std::string           three_letter;
    std::string           name;        ///< from _chem_comp.name (e.g. PHENYLALANINE)
    std::string           type;        ///< from _chem_comp.type (e.g. "L-PEPTIDE LINKING")
    int                   formal_charge_total = 0;
    std::vector<CcdAtom>  atoms;
    std::vector<CcdBond>  bonds;
};


// Load one residue's CCD entry into a typed CcdResidue. Returns
// nullopt if the entry is not in the CCD file.
//
// String boundary: this function reads strings from cifpp. The
// strings live in the returned struct; downstream consumers (RDKit
// builder, reconciler) translate them to typed-enum values and
// discard.
std::optional<CcdResidue> LoadCcdResidue(cif::file& ccd,
                                          const std::string& three_letter,
                                          ProcessLog& log) {
    log.Section(std::string("load-ccd-residue::") + three_letter);

    cif::datablock* block_ptr = nullptr;
    for (auto& db : ccd) {
        if (db.name() == three_letter) {
            block_ptr = &db;
            break;
        }
    }
    if (block_ptr == nullptr) {
        log.KV("status", "MISSING");
        log.Note("Residue not found in CCD; check spelling or update CCD.");
        return std::nullopt;
    }
    auto& block = *block_ptr;

    CcdResidue res;
    res.three_letter = three_letter;

    // Scalar header fields. cifpp categories are accessed by name;
    // null-or-unknown values ('.' / '?') yield empty item_handle, for
    // which value_or returns the default we provide.
    auto get_str = [](auto&& handle, const char* dv) -> std::string {
        return handle.template value_or<std::string>(std::string(dv));
    };
    auto get_int = [](auto&& handle, int dv) -> int {
        return handle.template value_or<int>(dv);
    };

    auto* chem_comp = block.get("chem_comp");
    if (chem_comp != nullptr && !chem_comp->empty()) {
        const auto& row = chem_comp->front();
        res.name                = get_str(row["name"], "");
        res.type                = get_str(row["type"], "");
        res.formal_charge_total = get_int(row["pdbx_formal_charge"], 0);
    }
    log.KV("name", res.name);
    log.KV("type", res.type);
    log.KV("formal_charge_total", res.formal_charge_total);

    // _chem_comp_atom loop.
    auto* atom_cat = block.get("chem_comp_atom");
    if (atom_cat != nullptr) {
        for (const auto& row : *atom_cat) {
            CcdAtom a;
            a.atom_id            = get_str(row["atom_id"], "");
            a.alt_atom_id        = get_str(row["alt_atom_id"], "");
            a.type_symbol        = get_str(row["type_symbol"], "");
            a.formal_charge      = get_int(row["charge"], 0);
            a.aromatic_flag      = (get_str(row["pdbx_aromatic_flag"], "N") == "Y");
            a.leaving_atom_flag  = (get_str(row["pdbx_leaving_atom_flag"], "N") == "Y");
            a.stereo_config      = get_str(row["pdbx_stereo_config"], "N");
            a.ordinal            = get_int(row["pdbx_ordinal"], 0);
            res.atoms.push_back(std::move(a));
        }
    }
    log.KV("atom_count", static_cast<int>(res.atoms.size()));

    // _chem_comp_bond loop.
    auto* bond_cat = block.get("chem_comp_bond");
    if (bond_cat != nullptr) {
        for (const auto& row : *bond_cat) {
            CcdBond b;
            b.atom_id_1     = get_str(row["atom_id_1"], "");
            b.atom_id_2     = get_str(row["atom_id_2"], "");
            b.value_order   = get_str(row["value_order"], "");
            b.aromatic_flag = (get_str(row["pdbx_aromatic_flag"], "N") == "Y");
            b.stereo_config = get_str(row["pdbx_stereo_config"], "N");
            b.ordinal       = get_int(row["pdbx_ordinal"], 0);
            res.bonds.push_back(std::move(b));
        }
    }
    log.KV("bond_count", static_cast<int>(res.bonds.size()));
    log.KV("status", "OK");
    return res;
}


// ============================================================================
// CCD-to-RDKit: build an RDKit RWMol from a typed CcdResidue.
//
// The RDKit mol carries the chemistry graph for perception: aromatic
// and hybridisation perception, CIP labelling, canonical ranking,
// bond order. Atom names from CCD are stored as atom properties so
// RDKit's perception output can be mapped back to the CCD atom
// identities.
//
// Returns the RWMol plus a parallel vector mapping
// rdkit_atom_idx -> CCD atom_id, used by the reconciler.
// ============================================================================

struct RdkitMolWithMap {
    std::unique_ptr<RDKit::RWMol>  mol;
    std::vector<std::string>       rdkit_idx_to_ccd_atom_id;
};

int CcdElementAtomicNumber(const std::string& type_symbol) {
    static const std::map<std::string, int> kSymbolToZ = {
        {"H",  1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B",  5},
        {"C",  6}, {"N",  7}, {"O",  8}, {"F",  9}, {"Ne",10},
        {"Na",11}, {"Mg",12}, {"Al",13}, {"Si",14}, {"P", 15},
        {"S", 16}, {"Cl",17}, {"Ar",18}, {"K", 19}, {"Ca",20},
        {"Fe",26}, {"Br",35}, {"I", 53},
    };
    auto it = kSymbolToZ.find(type_symbol);
    return it != kSymbolToZ.end() ? it->second : 0;
}

RDKit::Bond::BondType CcdBondOrderToRdkit(const std::string& value_order) {
    if (value_order == "SING") return RDKit::Bond::SINGLE;
    if (value_order == "DOUB") return RDKit::Bond::DOUBLE;
    if (value_order == "TRIP") return RDKit::Bond::TRIPLE;
    if (value_order == "AROM") return RDKit::Bond::AROMATIC;
    if (value_order == "DELO") return RDKit::Bond::ONEANDAHALF;  // delocalised
    if (value_order == "QUAD") return RDKit::Bond::QUADRUPLE;
    return RDKit::Bond::UNSPECIFIED;
}

RdkitMolWithMap BuildRdkitMol(const CcdResidue& res, ProcessLog& log) {
    log.Section(std::string("build-rdkit-mol::") + res.three_letter);

    RdkitMolWithMap out;
    out.mol = std::make_unique<RDKit::RWMol>();

    std::map<std::string, unsigned int> ccd_id_to_rdkit_idx;

    // Add atoms.
    for (const auto& ccd_atom : res.atoms) {
        const int z = CcdElementAtomicNumber(ccd_atom.type_symbol);
        if (z <= 0) {
            log.KV("warning_unknown_element", ccd_atom.type_symbol);
            continue;
        }
        auto* rd_atom = new RDKit::Atom(z);
        rd_atom->setFormalCharge(ccd_atom.formal_charge);
        // No setIsAromatic here -- aromaticity perception runs as
        // part of sanitization. Carrying CCD's flag through would
        // double-up perception. The CCD aromatic flag becomes a
        // cross-check witness during reconciliation.
        const auto idx = out.mol->addAtom(rd_atom, /*updateLabel*/false, /*takeOwnership*/true);
        ccd_id_to_rdkit_idx[ccd_atom.atom_id] = idx;
        out.rdkit_idx_to_ccd_atom_id.push_back(ccd_atom.atom_id);
    }
    log.KV("rdkit_atoms_added", static_cast<int>(out.mol->getNumAtoms()));

    // Add bonds.
    int bonds_added = 0;
    int bonds_skipped_unknown_atom = 0;
    for (const auto& b : res.bonds) {
        auto it1 = ccd_id_to_rdkit_idx.find(b.atom_id_1);
        auto it2 = ccd_id_to_rdkit_idx.find(b.atom_id_2);
        if (it1 == ccd_id_to_rdkit_idx.end() || it2 == ccd_id_to_rdkit_idx.end()) {
            ++bonds_skipped_unknown_atom;
            continue;
        }
        const auto bond_type = CcdBondOrderToRdkit(b.value_order);
        out.mol->addBond(it1->second, it2->second, bond_type);
        ++bonds_added;
    }
    log.KV("rdkit_bonds_added", bonds_added);
    if (bonds_skipped_unknown_atom > 0) {
        log.KV("bonds_skipped_unknown_atom", bonds_skipped_unknown_atom);
    }

    return out;
}


// Run RDKit's standard perception pipeline. After this call:
//  - Aromaticity is perceived (Atom::getIsAromatic, Bond::getIsAromatic).
//  - Hybridisation is perceived (Atom::getHybridization).
//  - Implicit valences are computed.
//  - Stereochemistry can be assigned by AssignStereochemistry.
//
// The CCD entries already carry consistent valences (the canonical
// free-amino-acid form has its protonated termini), so sanitization
// generally succeeds. Failures here indicate a malformed CCD entry
// or unusual element/bond combination that needs investigation.
bool SanitizeMol(RDKit::RWMol& mol, ProcessLog& log) {
    try {
        unsigned int operations_performed = 0;
        RDKit::MolOps::sanitizeMol(mol, operations_performed);
        log.KV("sanitize_status", "OK");
        log.KV("sanitize_ops_mask_hex",
               std::string("0x") + std::to_string(operations_performed));
        return true;
    } catch (const std::exception& e) {
        log.KV("sanitize_status", "FAILED");
        log.KV("sanitize_error", e.what());
        std::fprintf(stderr, "ERROR sanitizing mol: %s\n", e.what());
        return false;
    }
}


// Inspect the perceived molecule and log RDKit-derived facts: per-atom
// aromaticity, hybridisation, ring membership, and the SMILES
// (round-trip sanity check). At this stage we do not yet write to
// the typed-enum tables; we just confirm the perception layer works.
void LogRdkitPerception(const RDKit::RWMol& mol,
                         const std::vector<std::string>& idx_to_ccd_id,
                         ProcessLog& log,
                         const std::string& residue_3letter) {
    log.Section(std::string("rdkit-perception::") + residue_3letter);
    const auto* ring_info = mol.getRingInfo();
    log.KV("ring_count", static_cast<int>(ring_info->numRings()));

    // Count aromatic rings manually (RDKit 2023.09 does not expose
    // numAromaticRings on RingInfo).
    int aromatic_ring_count = 0;
    for (const auto& bond_ring : ring_info->bondRings()) {
        bool all_aromatic = !bond_ring.empty();
        for (int bidx : bond_ring) {
            if (!mol.getBondWithIdx(bidx)->getIsAromatic()) {
                all_aromatic = false;
                break;
            }
        }
        if (all_aromatic) ++aromatic_ring_count;
    }
    log.KV("aromatic_ring_count", aromatic_ring_count);

    int aromatic_atoms = 0;
    int sp_count = 0, sp2_count = 0, sp3_count = 0, other_hyb = 0;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const auto* a = mol.getAtomWithIdx(i);
        if (a->getIsAromatic()) ++aromatic_atoms;
        switch (a->getHybridization()) {
            case RDKit::Atom::SP:   ++sp_count;  break;
            case RDKit::Atom::SP2:  ++sp2_count; break;
            case RDKit::Atom::SP3:  ++sp3_count; break;
            default:                ++other_hyb; break;
        }
    }
    log.KV("aromatic_atom_count", aromatic_atoms);
    log.KV("hybridisation_sp",  sp_count);
    log.KV("hybridisation_sp2", sp2_count);
    log.KV("hybridisation_sp3", sp3_count);
    log.KV("hybridisation_other_or_unset", other_hyb);

    // Round-trip SMILES as a sanity probe. The CCD entry should
    // produce a reasonable SMILES; obviously-broken perception shows
    // up here as a malformed string or exception.
    try {
        const std::string smi = RDKit::MolToSmiles(mol);
        log.KV("smiles_roundtrip", smi);
    } catch (const std::exception& e) {
        log.KV("smiles_roundtrip_error", e.what());
    }

    // First few atoms with their CCD names: spot-check that the
    // atom-index-to-ccd-id map is sensible.
    for (unsigned int i = 0; i < mol.getNumAtoms() && i < 6u; ++i) {
        const auto* a = mol.getAtomWithIdx(i);
        const std::string ccd_id = (i < idx_to_ccd_id.size()) ? idx_to_ccd_id[i] : "?";
        std::string fact;
        fact += "Z=" + std::to_string(a->getAtomicNum());
        fact += " arom=" + std::to_string(static_cast<int>(a->getIsAromatic()));
        fact += " hyb=" + std::to_string(static_cast<int>(a->getHybridization()));
        log.KV(std::string("atom_") + std::to_string(i) + "_" + ccd_id, fact);
    }
}


// ============================================================================
// Atom-name parser
//
// Lifted to src/generated/LegacyAmberSemanticTables.h so the runtime
// composition path and this generator share a SINGLE source of truth.
// `ParseAtomName`, `ElementFromAtomName`, `LocantLetterToEnum`, and
// the `ParsedAtomName` struct live there as inline functions; the
// `using` declarations near the top of this file bring them into local
// scope under the legacy generator-side name `ParsedName` so the
// existing call sites (SynthesisedFor*, BuildAtomSemanticEntry, the
// per-residue loop) compile unchanged.
//
// String boundary: the parser consumes CCD atom-name strings and
// returns typed-enum values. The strings die at the function's
// return; callers see typed values only.
// ============================================================================


// Build a map from each atom name to its first heavy-atom neighbour
// (used to derive H -> parent for ParseAtomName).
std::map<std::string, std::string> BuildHydrogenParentMap(const CcdResidue& res) {
    std::map<std::string, std::string> parents;
    for (const auto& b : res.bonds) {
        const bool a1_is_h = !b.atom_id_1.empty() && b.atom_id_1.front() == 'H';
        const bool a2_is_h = !b.atom_id_2.empty() && b.atom_id_2.front() == 'H';
        if (a1_is_h && !a2_is_h) parents[b.atom_id_1] = b.atom_id_2;
        else if (a2_is_h && !a1_is_h) parents[b.atom_id_2] = b.atom_id_1;
    }
    return parents;
}


// ============================================================================
// RDKit perception extraction -- typed values for ProchiralStereo,
// hybridisation, aromatic flag, ring info, equivalence class.
//
// String boundary: RDKit returns its perception as typed C++ values
// already. The CIPLabeler returns "R"/"S"/etc. as a property string;
// we translate to the typed enum here.
// ============================================================================

ProchiralStereo RdkitCipCodeToEnum(const std::string& code) {
    if (code == "R") return ProchiralStereo::ProR;
    if (code == "S") return ProchiralStereo::ProS;
    return ProchiralStereo::Unassigned;
}

Hybridisation RdkitHybToEnum(RDKit::Atom::HybridizationType h) {
    switch (h) {
        case RDKit::Atom::SP:   return Hybridisation::sp;
        case RDKit::Atom::SP2:  return Hybridisation::sp2;
        case RDKit::Atom::SP3:  return Hybridisation::sp3;
        default:                return Hybridisation::Unassigned;
    }
}

BondOrderToNeighbour RdkitBondTypeToEnum(RDKit::Bond::BondType t) {
    switch (t) {
        case RDKit::Bond::SINGLE:        return BondOrderToNeighbour::Single;
        case RDKit::Bond::DOUBLE:        return BondOrderToNeighbour::Double;
        case RDKit::Bond::TRIPLE:        return BondOrderToNeighbour::Triple;
        case RDKit::Bond::AROMATIC:      return BondOrderToNeighbour::Aromatic;
        case RDKit::Bond::ONEANDAHALF:   return BondOrderToNeighbour::Delocalised;
        case RDKit::Bond::QUADRUPLE:     return BondOrderToNeighbour::Quadruple;
        default:                         return BondOrderToNeighbour::None;
    }
}


struct RdkitFacts {
    // Per-rdkit-atom-index facts.
    std::vector<bool>            aromatic;
    std::vector<Hybridisation>   hybridisation;
    std::vector<ProchiralStereo> cip_stereo;
    std::vector<uint8_t>         canonical_rank;   ///< RDKit equivalence-class index
    std::vector<BondOrderMask>   bond_orders;      ///< up to 4 bonded neighbours
    std::vector<bool>            in_ring;
};

RdkitFacts ExtractRdkitFacts(RDKit::RWMol& mol, ProcessLog& log) {
    log.Section("rdkit-facts-extraction");

    // Run CIP labelling. Annotates atoms with property "_CIPCode" = "R"/"S".
    try {
        RDKit::CIPLabeler::assignCIPLabels(mol);
        log.KV("cip_labeller_status", "OK");
    } catch (const std::exception& e) {
        log.KV("cip_labeller_status", "FAILED");
        log.KV("cip_labeller_error", e.what());
    }

    // Canonical-rank = equivalence class within molecule. Atoms with the
    // same rank are NMR-equivalent (e.g. methyl Hs collapse, ring-flip-fast
    // ortho atoms collapse).
    std::vector<unsigned int> rank;
    try {
        RDKit::Canon::rankMolAtoms(mol, rank, /*breakTies*/ false);
        log.KV("canonical_rank_status", "OK");
    } catch (const std::exception& e) {
        log.KV("canonical_rank_status", "FAILED");
        log.KV("canonical_rank_error", e.what());
        rank.assign(mol.getNumAtoms(), 0);
    }

    RdkitFacts facts;
    const auto n = mol.getNumAtoms();
    facts.aromatic.resize(n);
    facts.hybridisation.resize(n);
    facts.cip_stereo.resize(n, ProchiralStereo::NotProchiral);
    facts.canonical_rank.resize(n);
    facts.bond_orders.resize(n);
    facts.in_ring.resize(n);

    const auto* ring_info = mol.getRingInfo();
    int cip_assigned = 0;
    int aromatic_count = 0;
    for (unsigned int i = 0; i < n; ++i) {
        const auto* a = mol.getAtomWithIdx(i);
        facts.aromatic[i]      = a->getIsAromatic();
        facts.hybridisation[i] = RdkitHybToEnum(a->getHybridization());
        facts.canonical_rank[i] = static_cast<uint8_t>(
            i < rank.size() ? std::min<unsigned>(rank[i], 255u) : 0u);
        facts.in_ring[i] = ring_info->numAtomRings(i) > 0;
        if (facts.aromatic[i]) ++aromatic_count;

        std::string cip;
        if (a->getPropIfPresent(RDKit::common_properties::_CIPCode, cip) && !cip.empty()) {
            facts.cip_stereo[i] = RdkitCipCodeToEnum(cip);
            ++cip_assigned;
        }

        // Per-atom bond orders to up to 4 neighbours.
        unsigned int slot = 0;
        for (const auto& nbr : mol.atomBonds(a)) {
            if (slot >= 4) break;
            facts.bond_orders[i].orders[slot] = RdkitBondTypeToEnum(nbr->getBondType());
            ++slot;
        }
    }
    log.KV("aromatic_atom_count", aromatic_count);
    log.KV("cip_labels_assigned", cip_assigned);
    return facts;
}


// ============================================================================
// Synthesised-field tables -- per-residue chemistry encoding.
//
// PlanarGroupKind, PolarHKind, PseudoatomMembership, RingSystem and
// RingPosition labels for the standard 20 residues. Encoded as a
// per-residue lookup function.
//
// Citation pattern: each lookup carries a comment naming the chemistry
// authority. PolarHKind from Wuethrich 1986 + Englander 2008 frame;
// PlanarGroupKind from Pauling 1951 / Cantor & Schimmel 1980;
// pseudoatom set from Markley 1998 Table 1 (Wuethrich 1983 extension).
// RingPosition from Vollhardt + Joule & Mills heterocycle conventions.
// ============================================================================

struct SynthesisedFields {
    PlanarGroupKind       planar_group  = PlanarGroupKind::None;
    PlanarStereo          planar_stereo = PlanarStereo::NotApplicable;
    PolarHKind            polar_h       = PolarHKind::NotPolar;
    PseudoatomMembership  pseudoatom    = {};
    RingPosition          ring_pos      = {};

    // Markley-derived prochiral assignment (ProR / ProS for
    // prochiral methylene Hs and prochiral methyl-bearing branch
    // carbons). Hand-encoded in the residue reference because
    // RDKit CIPLabeler labels chiral centres, NOT prochiral Hs --
    // see the reference doc + dependencies §C.1. When the
    // synthesised value is not NotProchiral, it OVERRIDES the
    // RDKit-perceived value (which is None for prochiral Hs in
    // CCD's free-amino-acid forms with explicit Hs).
    ProchiralStereo       prochiral     = ProchiralStereo::NotProchiral;

    // Optional per-atom formal-charge override (variant-specific).
    // The CCD entry's per-atom charge is the default authority. When
    // a variant's protonation state differs from the CCD parent's
    // chemistry (e.g. CYM changes Sγ from neutral to -1; LYN changes
    // Nζ from +1 to 0; HID/HIE remove the imidazolium +1 from Nδ1;
    // HIP relocates the +1 from Nδ1 to Nε2), the variant patch sets
    // this override and BuildAtomSemanticEntry uses it instead of
    // the CCD value. std::optional<> empty = "use CCD"; populated =
    // "variant overrides CCD here." The standard 20 patch functions
    // never set this, so their output is unchanged.
    std::optional<int8_t> charge_override = std::nullopt;
};


// ============================================================================
// Per-residue + per-variant synthesised-field tables
// ============================================================================
//
// Each function below encodes one (residue, variant) block from the
// signed-off reference document
// `spec/plan/bones/topology-residue-reference-2026-05-05.md` Section 3
// (alphabetical per-residue blocks) + Section 2 conventions (charge
// placement, residue-specific E/Z, Markley alternation rule).
//
// Architecture: every variant has its own CCD entry (verified
// 2026-05-05: HID/HIE/HIP/CYX/CYM/LYN/ARN/TYM/ASH/GLH all present in
// data/ccd/components.cif). Therefore each variant has its own
// SynthesisedFor<Variant> function dispatched by `res.three_letter`,
// just like PHE. No variant-aware branching inside the residue-default
// function is needed.
//
// Markley alternation rule applied uniformly:
//   β-onward methylene Hs in standard L-amino acids:
//     H<locant>2 = ProS; H<locant>3 = ProR.
//   Glycine Hα is INVERTED: HA2 = ProR; HA3 = ProS.
//   Branched heavy atoms (Val, Leu, Ile) per Markley Fig 1.
//
// Citation: every block carries a single-line pointer to the
// reference doc's Section 3 entry for that residue / variant.
// ============================================================================


// Markley alternation rule for prochiral methylene Hs in standard
// L-amino acids (β-onward). Returns ProS for "...2", ProR for "...3",
// NotProchiral otherwise. Glycine Hα is NOT covered here (it inverts;
// see SynthesisedForGly). Branched-heavy-atom prochirality is
// per-residue and also encoded explicitly per function.
ProchiralStereo MarkleyMethyleneProchiral(const std::string& atom_id) {
    if (atom_id.size() < 3 || atom_id[0] != 'H') return ProchiralStereo::NotProchiral;
    // The matching atom-id forms covered by the alternation rule:
    //   HB2/HB3, HG2/HG3, HD2/HD3, HE2/HE3 (no branch-digit suffix).
    // Forms with a digit between locant and inner index (e.g. HD11
    // for Leu's δ1-methyl, HG21 for Val's γ2-methyl) are methyl-Hs,
    // not prochiral methylene Hs; do NOT label them ProR/ProS.
    if (atom_id.size() != 3) return ProchiralStereo::NotProchiral;
    const char locant = atom_id[1];
    const char digit  = atom_id[2];
    if (locant != 'B' && locant != 'G' && locant != 'D' && locant != 'E') {
        return ProchiralStereo::NotProchiral;
    }
    if (digit == '2') return ProchiralStereo::ProS;
    if (digit == '3') return ProchiralStereo::ProR;
    return ProchiralStereo::NotProchiral;
}


// ALA -- alanine. Reference Section 3 (ALA block).
SynthesisedFields SynthesisedForAla(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;

    // Backbone amide plane (N + C + O + H form the planar peptide unit).
    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    // β-methyl pseudoatom (Markley MB; 3 equivalent β-methyl Hs).
    if (atom_id == "HB1" || atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::M,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    return s;
}


// ARG -- arginine (default, charged). Reference Section 3 (ARG block).
// Charge placement (Section 2): +1 on Nε per Lewis convention; Cζ +
// Nη1 + Nη2 = 0. E/Z mapping: Nη1 cis to Cδ ⇒ Z; Nη2 trans ⇒ E.
// Hηxx labels: BMRB confirms HH11=Z (cis to Nε via Nη1), HH12=E,
// HH21=Z (cis to Nε via Nη2), HH22=E. Note this is OPPOSITE polarity
// from Asn/Gln because on Cζ, Nε IS the higher-priority substituent.
//
// CCD ARG carries +1 on NH2 (its Lewis convention places +1 on a
// nitrogen of the guanidinium); we override to put +1 on Nε per the
// project's chosen Lewis-convention placement (reference Section 2),
// and zero out NH2's CCD-default +1.
SynthesisedFields SynthesisedForArg(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β/γ/δ methylenes

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    // Lewis-convention charge placement (Section 2 of reference doc):
    // +1 on Nε; CCD default has +1 on NH2 (also a valid Lewis form, but
    // not the project's canonical placement). Override both.
    if (atom_id == "NE")  s.charge_override = static_cast<int8_t>(+1);
    if (atom_id == "NH2") s.charge_override = static_cast<int8_t>(0);

    // β methylene QB.
    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    // γ methylene QG.
    if (atom_id == "HG2" || atom_id == "HG3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Gamma), 0, false};
    }
    // δ methylene QD.
    if (atom_id == "HD2" || atom_id == "HD3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Delta), 0, false};
    }

    // Guanidinium group: NE + HE + CZ + NH1/HH11/HH12 + NH2/HH21/HH22.
    if (atom_id == "NE" || atom_id == "HE" || atom_id == "CZ" ||
        atom_id == "NH1" || atom_id == "HH11" || atom_id == "HH12" ||
        atom_id == "NH2" || atom_id == "HH21" || atom_id == "HH22") {
        s.planar_group = PlanarGroupKind::Guanidinium;
    }
    if (atom_id == "HE")    s.polar_h = PolarHKind::GuanidiniumNH;
    if (atom_id == "HH11" || atom_id == "HH12" ||
        atom_id == "HH21" || atom_id == "HH22") {
        s.polar_h = PolarHKind::GuanidiniumNH;
    }

    // E/Z labels on the Nη atoms and their Hs (Section 2 of reference).
    if (atom_id == "NH1")  s.planar_stereo = PlanarStereo::Z;   // Nη1 cis to Cδ
    if (atom_id == "NH2")  s.planar_stereo = PlanarStereo::E;   // Nη2 trans to Cδ
    if (atom_id == "HH11") s.planar_stereo = PlanarStereo::Z;   // BMRB: Z
    if (atom_id == "HH12") s.planar_stereo = PlanarStereo::E;   // BMRB: E
    if (atom_id == "HH21") s.planar_stereo = PlanarStereo::Z;   // BMRB: Z
    if (atom_id == "HH22") s.planar_stereo = PlanarStereo::E;   // BMRB: E

    // QH super-group across all four Hη; per-atom QH1 vs QH2 by outer.
    if (atom_id == "HH11" || atom_id == "HH12") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Eta), 1, true};
    }
    if (atom_id == "HH21" || atom_id == "HH22") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Eta), 2, true};
    }
    return s;
}


// ARN -- deprotonated arginine variant. Reference Section 3 (ARN
// block) + dependencies §B.1 + §F. ff14SB has no canonical ARN
// template; HE is removed by chemistry inference (lone pair on Nε).
// AmberAminoAcidVariantTable["ARN"] verified 2026-05-05: charge=0,
// no explicit atom inventory, so reference's HE-removal stands.
//
// ARG default carries +1 on Nε (HE-bearing guanidinium nitrogen) and
// 0 on NH2 per Lewis convention. ARN drops HE and the lone pair lives
// on Nε, so ARN's Nε is neutral (0). NH2 stays 0. (CCD's +1 on NH2
// was zeroed by ARG already; this remains 0 for ARN.)
SynthesisedFields SynthesisedForArn(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForArg(atom_id, parsed);
    // HE is removed at the pipeline level; defensive zeroing here.
    if (atom_id == "HE") {
        s.polar_h = PolarHKind::NotPolar;
    }
    // ARG default put +1 on Nε; ARN is neutral, so re-override to 0.
    if (atom_id == "NE")  s.charge_override = static_cast<int8_t>(0);
    if (atom_id == "NH2") s.charge_override = static_cast<int8_t>(0);
    return s;
}


// ASN -- asparagine. Reference Section 3 (ASN block).
// E/Z labels (CRITICAL, recently corrected per dependencies §E.1-2):
// On Cγ=Nδ2 planar bond, OD1 (O) is the high-priority Cγ substituent.
// HD21 = "1" = cis to Cβ = trans to OD1 = E (BMRB confirms).
// HD22 = "2" = trans to Cβ = cis to OD1 = Z (BMRB confirms).
SynthesisedFields SynthesisedForAsn(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β methylene only

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    // β methylene QB.
    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }

    // Side-chain amide group: CG + OD1 + ND2 + HD21 + HD22.
    if (atom_id == "CG"  || atom_id == "OD1" || atom_id == "ND2" ||
        atom_id == "HD21" || atom_id == "HD22") {
        s.planar_group = PlanarGroupKind::SidechainAmide;
    }
    if (atom_id == "HD21" || atom_id == "HD22") {
        s.polar_h = PolarHKind::SidechainPrimaryAmide;
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Delta), 0, false};
    }
    if (atom_id == "HD21") s.planar_stereo = PlanarStereo::E;   // BMRB
    if (atom_id == "HD22") s.planar_stereo = PlanarStereo::Z;   // BMRB

    return s;
}


// ASP -- aspartate. Reference Section 3 (ASP block).
// Carboxylate: -1 on Oδ2 per Lewis convention; Cγ + Oδ1 carry 0.
//
// CCD ASP entry is the *protonated* (neutral) form, so its OD2 has
// formal_charge=0. The deprotonated chain form (project default) needs
// an explicit -1 override on OD2 per Section 2 of the reference doc.
// (The protonated variant ASH inherits SynthesisedForAsp; ASH then
// re-overrides OD2 back to 0 -- see SynthesisedForAsh.)
SynthesisedFields SynthesisedForAsp(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β methylene only

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }

    // Side-chain carboxylate: CG + OD1 + OD2.
    if (atom_id == "CG" || atom_id == "OD1" || atom_id == "OD2") {
        s.planar_group = PlanarGroupKind::Carboxylate;
    }

    // Lewis-convention charge placement (reference Section 2):
    // -1 on Oδ2; Oδ1 = 0 (delocalised in reality, Lewis-localised here).
    if (atom_id == "OD2") s.charge_override = static_cast<int8_t>(-1);
    return s;
}


// ASH -- protonated aspartate variant. Reference Section 3 (ASH block).
// Note that the CCD's "ASP" entry IS the protonated form (it carries
// HD2 + OXT + HXT and has total formal_charge=0). The standard ASP
// emission therefore already passes through the HD2 atom; the ASH
// variant's job is to label that atom's chemistry correctly
// (CarboxylOH on HD2). No atom-add, no atom-remove.
//
// Charge: SynthesisedForAsp sets -1 on OD2 (the deprotonated chain
// form). ASH is the protonated (neutral) form, so we re-override OD2
// back to 0 here -- HD2 is bonded to OD2 in ASH, making OD2 a neutral
// hydroxyl oxygen.
SynthesisedFields SynthesisedForAsh(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForAsp(atom_id, parsed);
    // HD2 is the protonated carboxyl OH on Oδ2 in ASH.
    if (atom_id == "HD2") {
        s.planar_group = PlanarGroupKind::Carboxylate;
        s.polar_h      = PolarHKind::CarboxylOH;
    }
    // ASH is neutral: undo the ASP -1 on OD2 (now a hydroxyl-bearing O).
    if (atom_id == "OD2") s.charge_override = static_cast<int8_t>(0);
    return s;
}


// CYS -- cysteine. Reference Section 3 (CYS block).
SynthesisedFields SynthesisedForCys(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β methylene only

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    if (atom_id == "HG") s.polar_h = PolarHKind::ThiolSH;     // reduced thiol H
    return s;
}


// CYX -- disulfide-bonded cysteine variant. Reference Section 3 (CYX
// block). HG removed at the variant pipeline level (Sγ-Sγ disulfide
// bond); the inter-residue disulfide topology is captured upstream
// via CovalentTopology. Sγ formal charge stays 0 (CCD: 0).
SynthesisedFields SynthesisedForCyx(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForCys(atom_id, parsed);
    // Defensive/idempotent: HG row is removed entirely at the
    // pipeline level; this branch ensures we never label a stray HG
    // entry as ThiolSH for CYX.
    if (atom_id == "HG") s.polar_h = PolarHKind::NotPolar;
    return s;
}


// CYM -- deprotonated thiolate cysteine variant. Reference Section 3
// (CYM block). HG removed at the variant pipeline level. Sγ formal
// charge: CCD has 0 (CYS); CYM is thiolate so override SG to -1.
SynthesisedFields SynthesisedForCym(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForCys(atom_id, parsed);
    if (atom_id == "HG") s.polar_h = PolarHKind::NotPolar;
    if (atom_id == "SG") s.charge_override = static_cast<int8_t>(-1);
    return s;
}


// GLN -- glutamine. Reference Section 3 (GLN block).
// E/Z labels (CRITICAL, recently corrected per dependencies §E.3-4):
// On Cδ=Nε2 planar bond, OE1 (O) is the high-priority Cδ substituent.
// HE21 = "1" = cis to Cγ = trans to OE1 = E (BMRB confirms).
// HE22 = "2" = trans to Cγ = cis to OE1 = Z (BMRB confirms).
SynthesisedFields SynthesisedForGln(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β + γ methylenes

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    if (atom_id == "HG2" || atom_id == "HG3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Gamma), 0, false};
    }

    // Side-chain amide group: CD + OE1 + NE2 + HE21 + HE22.
    if (atom_id == "CD" || atom_id == "OE1" || atom_id == "NE2" ||
        atom_id == "HE21" || atom_id == "HE22") {
        s.planar_group = PlanarGroupKind::SidechainAmide;
    }
    if (atom_id == "HE21" || atom_id == "HE22") {
        s.polar_h = PolarHKind::SidechainPrimaryAmide;
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Epsilon), 0, false};
    }
    if (atom_id == "HE21") s.planar_stereo = PlanarStereo::E;  // BMRB
    if (atom_id == "HE22") s.planar_stereo = PlanarStereo::Z;  // BMRB

    return s;
}


// GLU -- glutamate. Reference Section 3 (GLU block).
// Carboxylate: -1 on Oε2 per Lewis convention.
//
// CCD GLU entry is the *protonated* (neutral) form (formal_charge=0;
// carries HE2 + OXT + HXT). The deprotonated chain form (project
// default) needs an explicit -1 override on OE2. (GLH inherits and
// re-overrides OE2 back to 0 -- see SynthesisedForGlh.)
SynthesisedFields SynthesisedForGlu(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β + γ methylenes

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    if (atom_id == "HG2" || atom_id == "HG3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Gamma), 0, false};
    }

    // Side-chain carboxylate: CD + OE1 + OE2.
    if (atom_id == "CD" || atom_id == "OE1" || atom_id == "OE2") {
        s.planar_group = PlanarGroupKind::Carboxylate;
    }

    // Lewis-convention charge placement (reference Section 2):
    // -1 on Oε2; Oε1 = 0.
    if (atom_id == "OE2") s.charge_override = static_cast<int8_t>(-1);
    return s;
}


// GLH -- protonated glutamate variant. Reference Section 3 (GLH
// block). Like ASP→ASH above, the CCD "GLU" entry IS the protonated
// form (it carries HE2 + OXT + HXT, total formal_charge=0). The
// GLH variant labels HE2's chemistry as CarboxylOH; the SynthesisedForGlu
// default sets -1 on OE2 (deprotonated chain form), so GLH needs to
// re-override OE2 back to 0 (HE2 bonded to OE2 makes it a neutral
// hydroxyl oxygen).
SynthesisedFields SynthesisedForGlh(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForGlu(atom_id, parsed);
    if (atom_id == "HE2") {
        s.planar_group = PlanarGroupKind::Carboxylate;
        s.polar_h      = PolarHKind::CarboxylOH;
    }
    // GLH is neutral: undo the GLU -1 on OE2.
    if (atom_id == "OE2") s.charge_override = static_cast<int8_t>(0);
    return s;
}


// GLY -- glycine. Reference Section 3 (GLY block).
// THE INVERSION: HA2 = ProR, HA3 = ProS (the only standard-20 case
// where the alternation rule reverses; Markley Fig 1 marks HA2 as R).
// Prochiral fields are populated mechanically by ParsedName for now;
// the synthesised block contributes the QA pseudoatom membership only.
// (Note: HA2/HA3 ProR/ProS assignment is out-of-scope here per the
// reference's guidance to use the prochiral column directly per atom;
// the C++ generator's `prochiral` field is sourced from RDKit
// CIPLabeler, which handles chiral centres but not glycine prochiral
// Hs. The Gly inversion is the canonical Markley result; if the
// runtime field is empty for Gly Hα atoms, the pre-CIP fallback rule
// is "ProR = HA2; ProS = HA3" -- INVERTED from typical alternation.)
SynthesisedFields SynthesisedForGly(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    // Glycine's two diastereotopic Hα form the QA pseudoatom group.
    // INVERSION: HA2 = ProR; HA3 = ProS (Markley Fig 1; the only
    // standard-20 case where the alternation rule reverses).
    if (atom_id == "HA2") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Alpha), 0, false};
        s.prochiral  = ProchiralStereo::ProR;
    }
    if (atom_id == "HA3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Alpha), 0, false};
        s.prochiral  = ProchiralStereo::ProS;
    }
    return s;
}


// HIS -- histidine. AMBER ff14SB DEFAULT = HIE (Nε2-protonated,
// neutral). The bare 3-letter "HIS" therefore maps to the HIE
// encoding: no Hδ1, has Hε2; ND1 = Heteroatom_NoH, NE2 = Heteroatom_NH.
// Reference Section 3 (HIS block). Locked decision per dependencies
// §D.1 (HIE-default-for-HIS); cross-check at src/AmberLeapInput.cpp:29.
//
// CCD HIS entry is the imidazolium form (+1 on Nδ1, with HD1 + HE2
// both present). The default chain form is the neutral HIE-equivalent:
// Nδ1 + Nε2 + ring carbons all carry 0, and HD1 is partitioned out by
// the §H.10 chain-atom whitelist filter. Override Nδ1 +1 -> 0 here.
SynthesisedFields SynthesisedForHis(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β methylene only
    // Aromatic-CH / imidazole-NH override: HD2 (aromatic CH on Cδ2)
    // and HE2 (imidazole NH on Nε2) match the helper's H[BGDE][23]
    // pattern but are NOT prochiral methylene Hs. Clear the helper's
    // false-positive ProS tag here; same cleanup pattern as PHE/TRP/TYR.
    if (atom_id == "HD2" || atom_id == "HE2") {
        s.prochiral = ProchiralStereo::NotProchiral;
    }

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }

    // Imidazole 5-ring: CG (ipso), ND1, CE1 (PyrroleAlpha), HE1, NE2,
    // HE2, CD2 (PyrroleBeta), HD2.
    auto set_ring = [&](RingPositionLabel pos, bool het_with_h) {
        s.ring_pos.primary.ring          = RingSystemKind::Imidazole_His;
        s.ring_pos.primary.position      = pos;
        s.ring_pos.primary.ring_size     = 5;
        s.ring_pos.primary.aromatic      = true;
        s.ring_pos.primary.planar        = true;
        s.ring_pos.primary.n_heteroatoms = 2;
        s.planar_group = PlanarGroupKind::Imidazole;
        (void)het_with_h;
    };

    if (atom_id == "CG")                           set_ring(RingPositionLabel::Ipso, false);
    else if (atom_id == "ND1")                     set_ring(RingPositionLabel::Heteroatom_NoH, false);
    else if (atom_id == "CE1" || atom_id == "HE1") set_ring(RingPositionLabel::PyrroleAlpha, false);
    else if (atom_id == "NE2")                     set_ring(RingPositionLabel::Heteroatom_NH, true);
    else if (atom_id == "HE2")                     set_ring(RingPositionLabel::Heteroatom_NH, true);
    else if (atom_id == "CD2" || atom_id == "HD2") set_ring(RingPositionLabel::PyrroleBeta, false);

    if (atom_id == "HE2") s.polar_h = PolarHKind::ImidazoleNH;

    // HIS = HIE default (neutral). CCD HIS has +1 on Nδ1 (imidazolium);
    // override Nδ1 to 0. NE2 remains 0 (CCD default).
    if (atom_id == "ND1") s.charge_override = static_cast<int8_t>(0);
    return s;
}


// HID -- Nδ1-protonated histidine variant. Reference Section 3 (HID
// block). HD1 retained from CCD parent; HE2 dropped at the variant
// pipeline level (see ProcessVariantResidue); ND1 becomes
// Heteroatom_NH; NE2 becomes Heteroatom_NoH. Formal charge: CCD's
// HIS entry has +1 on Nδ1 (it's the imidazolium form). HID is
// neutral, so we override Nδ1 +1 → 0.
SynthesisedFields SynthesisedForHid(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForHis(atom_id, parsed);
    // ND1 now bears Hδ1: change from Heteroatom_NoH to Heteroatom_NH.
    if (atom_id == "ND1" || atom_id == "HD1") {
        s.ring_pos.primary.position = RingPositionLabel::Heteroatom_NH;
        s.planar_group = PlanarGroupKind::Imidazole;
        s.ring_pos.primary.ring          = RingSystemKind::Imidazole_His;
        s.ring_pos.primary.ring_size     = 5;
        s.ring_pos.primary.aromatic      = true;
        s.ring_pos.primary.planar        = true;
        s.ring_pos.primary.n_heteroatoms = 2;
    }
    // NE2 now has no H: change to Heteroatom_NoH.
    if (atom_id == "NE2") {
        s.ring_pos.primary.position = RingPositionLabel::Heteroatom_NoH;
    }
    if (atom_id == "HD1") s.polar_h = PolarHKind::ImidazoleNH;
    // HE2 is removed entirely at the pipeline level (see
    // ProcessVariantResidue's atoms_to_remove), so no row will be
    // emitted for it; this branch is here for defensive/idempotent
    // patching only.
    if (atom_id == "HE2") s.polar_h = PolarHKind::NotPolar;
    // Variant formal charge: HID is neutral; CCD HIS has +1 on ND1
    // (imidazolium form). Override ND1 to 0 here.
    if (atom_id == "ND1") s.charge_override = static_cast<int8_t>(0);
    return s;
}


// HIE -- Nε2-protonated histidine variant. Reference Section 3 (HIE
// block). HD1 dropped at the variant pipeline level; HE2 retained.
// AMBER ff14SB's "HIE" matches what the existing SynthesisedForHis
// function already encodes for HE2/NE2 (HE2 = ImidazoleNH; NE2 =
// Heteroatom_NH; ND1 = Heteroatom_NoH); we therefore re-use the HIS
// patch wholesale here. CCD HIS has +1 on Nδ1 (imidazolium); HIE is
// neutral so we override Nδ1 +1 → 0.
SynthesisedFields SynthesisedForHie(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForHis(atom_id, parsed);
    // HIE is neutral; remove CCD's +1 on Nδ1.
    if (atom_id == "ND1") s.charge_override = static_cast<int8_t>(0);
    return s;
}


// HIP -- imidazolium (both NHs protonated, +1 charge). Reference
// Section 3 (HIP block). Charge placement (Section 2): +1 on Nε2 per
// chosen convention; ND1, CE1, CD2, CG carry 0. CCD's HIS entry has
// +1 on ND1; HIP relocates the +1 to NE2 per the chosen convention.
SynthesisedFields SynthesisedForHip(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForHis(atom_id, parsed);
    // Both ND1 and NE2 are NH-bearing in HIP.
    if (atom_id == "ND1" || atom_id == "HD1") {
        s.ring_pos.primary.position = RingPositionLabel::Heteroatom_NH;
        s.planar_group = PlanarGroupKind::Imidazole;
        s.ring_pos.primary.ring          = RingSystemKind::Imidazole_His;
        s.ring_pos.primary.ring_size     = 5;
        s.ring_pos.primary.aromatic      = true;
        s.ring_pos.primary.planar        = true;
        s.ring_pos.primary.n_heteroatoms = 2;
    }
    if (atom_id == "HD1") s.polar_h = PolarHKind::ImidazoleNH;
    // HE2 stays ImidazoleNH (inherited from HIS default).
    // Variant formal charge (Lewis-localised convention, Section
    // 2 / dependencies §C.4): place +1 on Nε2; clear the CCD's
    // ND1 +1.
    if (atom_id == "ND1") s.charge_override = static_cast<int8_t>(0);
    if (atom_id == "NE2") s.charge_override = static_cast<int8_t>(1);
    return s;
}


// ILE -- isoleucine. Reference Section 3 (ILE block).
// Branched aliphatic; CG2 (γ-methyl, ProS branch); CG1 (γ-methylene,
// HG12=ProS, HG13=ProR); CD1 (δ-methyl, only δ atom -> NotProchiral).
// MG and MD are NOT in higher-order Q super-groups (Ile has no QG_super
// because γ2-methyl and γ1-methylene are not equivalent in any
// pseudoatom).
SynthesisedFields SynthesisedForIle(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    // CG2: γ2-methyl carbon, ProS branch (only γ-methyl C in Ile).
    if (atom_id == "CG2") s.prochiral = ProchiralStereo::ProS;
    // γ2-methyl (CG2): MG with branch=2, super=false.
    if (atom_id == "HG21" || atom_id == "HG22" || atom_id == "HG23") {
        s.pseudoatom = {PseudoatomKind::M,
                        static_cast<uint8_t>(Locant::Gamma), 2, false};
    }
    // γ1-methylene (CG1): QG with branch=1, super=false. Hγ12=ProS,
    // Hγ13=ProR (Markley Fig 1; the asymmetric γ-methylene of Ile).
    if (atom_id == "HG12") s.prochiral = ProchiralStereo::ProS;
    if (atom_id == "HG13") s.prochiral = ProchiralStereo::ProR;
    if (atom_id == "HG12" || atom_id == "HG13") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Gamma), 1, false};
    }
    // δ1-methyl (CD1): MD with branch=1, super=false. Cδ1 is the only
    // δ atom (no diastereotopy) -> NotProchiral.
    if (atom_id == "HD11" || atom_id == "HD12" || atom_id == "HD13") {
        s.pseudoatom = {PseudoatomKind::M,
                        static_cast<uint8_t>(Locant::Delta), 1, false};
    }
    return s;
}


// LEU -- leucine. Reference Section 3 (LEU block).
// β methylene QB (NOT in QD super); γ has single Hγ (no diastereotopy);
// δ1-methyl (Cδ1=ProR) + δ2-methyl (Cδ2=ProS); MD1+MD2 both in QD
// super (all six δ-methyl Hs collapse to QD).
SynthesisedFields SynthesisedForLeu(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β methylene

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    // Branched δ-methyls: Cδ1 = ProR; Cδ2 = ProS (Markley Fig 1).
    if (atom_id == "CD1") s.prochiral = ProchiralStereo::ProR;
    if (atom_id == "CD2") s.prochiral = ProchiralStereo::ProS;
    // δ1-methyl: MD1 in QD super.
    if (atom_id == "HD11" || atom_id == "HD12" || atom_id == "HD13") {
        s.pseudoatom = {PseudoatomKind::M,
                        static_cast<uint8_t>(Locant::Delta), 1, true};
    }
    // δ2-methyl: MD2 in QD super.
    if (atom_id == "HD21" || atom_id == "HD22" || atom_id == "HD23") {
        s.pseudoatom = {PseudoatomKind::M,
                        static_cast<uint8_t>(Locant::Delta), 2, true};
    }
    return s;
}


// LYS -- lysine (default, charged). Reference Section 3 (LYS block).
// Charge placement: +1 on Nζ; HZ1/2/3 = AmmoniumNH (0). All four
// methylenes get QB/QG/QD/QE; Nζ has QZ (3 ammonium Hs).
SynthesisedFields SynthesisedForLys(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β/γ/δ/ε methylenes

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    if (atom_id == "HG2" || atom_id == "HG3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Gamma), 0, false};
    }
    if (atom_id == "HD2" || atom_id == "HD3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Delta), 0, false};
    }
    if (atom_id == "HE2" || atom_id == "HE3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Epsilon), 0, false};
    }
    // QZ over the ammonium Hs.
    if (atom_id == "HZ1" || atom_id == "HZ2" || atom_id == "HZ3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Zeta), 0, false};
        s.polar_h = PolarHKind::AmmoniumNH;
    }
    return s;
}


// LYN -- neutral lysine variant. Reference Section 3 (LYN block) +
// dependencies §B.2 + §E.5-6. CRITICAL (recently corrected): HZ1 is
// REMOVED (NOT HZ3); HZ2 + HZ3 retained as AmineNH (per the new
// PolarHKind::AmineNH enum). NZ formal charge: CCD LYS has +1 on
// NZ (charged ammonium); LYN is neutral so override NZ +1 → 0.
SynthesisedFields SynthesisedForLyn(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForLys(atom_id, parsed);
    // HZ1 removed at pipeline level; defensive/idempotent zeroing
    // here in case a row sneaks through.
    if (atom_id == "HZ1") {
        s.polar_h    = PolarHKind::NotPolar;
        s.pseudoatom = {};
    }
    if (atom_id == "HZ2" || atom_id == "HZ3") {
        s.polar_h = PolarHKind::AmineNH;
    }
    if (atom_id == "NZ") s.charge_override = static_cast<int8_t>(0);
    return s;
}


// MET -- methionine. Reference Section 3 (MET block).
// Sulphur sidechain: β methylene QB; γ methylene QG; thioether SD; ε
// methyl ME (3 equivalent Hs); no super-group.
SynthesisedFields SynthesisedForMet(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β + γ methylenes

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    if (atom_id == "HG2" || atom_id == "HG3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Gamma), 0, false};
    }
    // ε-methyl: HE1/HE2/HE3 form ME (Markley says 3 equivalent Hs).
    // Note: HE1/HE2/HE3 here are 3-H METHYL group (not methylene). The
    // helper would tag HE2/HE3 as ProS/ProR, but the methyl is not
    // prochiral -- override back to NotProchiral.
    if (atom_id == "HE1" || atom_id == "HE2" || atom_id == "HE3") {
        s.prochiral = ProchiralStereo::NotProchiral;
        s.pseudoatom = {PseudoatomKind::M,
                        static_cast<uint8_t>(Locant::Epsilon), 0, false};
    }
    return s;
}


// PHE -- phenylalanine (acceptance gate; original implementation
// from 2026-05-05 morning session). Reference Section 3 (PHE block).
// Aromatic 6-ring: CG (ipso), CD1 (ortho1), CD2 (ortho2), CE1
// (meta1), CE2 (meta2), CZ (para); Hs follow.
// PHE-only initial implementation; generalises in a follow-up commit
// once the architecture is proven on one residue.
SynthesisedFields SynthesisedForPhe(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s;
    // β methylene Hs: HB2=ProS, HB3=ProR. Ring Hs (HD1/HD2/HE1/HE2)
    // are caught by the helper but are NOT prochiral; override below.
    s.prochiral = MarkleyMethyleneProchiral(atom_id);
    if (atom_id == "HD1" || atom_id == "HD2" ||
        atom_id == "HE1" || atom_id == "HE2") {
        s.prochiral = ProchiralStereo::NotProchiral;
    }

    // Phe ring atoms: CG (ipso), CD1 (ortho1), CD2 (ortho2),
    //                 CE1 (meta1), CE2 (meta2), CZ (para).
    // Their hydrogens take the corresponding label too.
    auto set_ring = [&](RingPositionLabel label) {
        s.ring_pos.primary.ring          = RingSystemKind::Benzene_Phe;
        s.ring_pos.primary.position      = label;
        s.ring_pos.primary.ring_size     = 6;
        s.ring_pos.primary.aromatic      = true;
        s.ring_pos.primary.planar        = true;
        s.ring_pos.primary.n_heteroatoms = 0;
        s.planar_group = PlanarGroupKind::Aromatic6Ring;
    };
    // Note: PHE has no HG atom (CG is the ipso C, fully substituted).
    // The previous `|| atom_id == "HG"` clause was dead code; removed.
    if      (atom_id == "CG")                      set_ring(RingPositionLabel::Ipso);
    else if (atom_id == "CD1" || atom_id == "HD1") set_ring(RingPositionLabel::Ortho1);
    else if (atom_id == "CD2" || atom_id == "HD2") set_ring(RingPositionLabel::Ortho2);
    else if (atom_id == "CE1" || atom_id == "HE1") set_ring(RingPositionLabel::Meta1);
    else if (atom_id == "CE2" || atom_id == "HE2") set_ring(RingPositionLabel::Meta2);
    else if (atom_id == "CZ"  || atom_id == "HZ")  set_ring(RingPositionLabel::Para);

    // Polar-H classification.
    if (parsed.element == Element::H) {
        if      (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;
        else if (parsed.is_n_terminus)              s.polar_h = PolarHKind::AmmoniumNH;
        else if (parsed.is_c_terminus)              s.polar_h = PolarHKind::CarboxylOH;
        // ring Hs and CH/CHn are NotPolar.
    }

    // Backbone peptide planar group: N, CA(?, no, CA is sp3), C, O, H form
    // the planar peptide unit. Strictly: N, C, O, H of i; CA of i and i-1
    // are at the boundaries. In the canonical free-amino-acid CCD entry
    // we tag N + C + O + H as members.
    if (atom_id == "N" || atom_id == "C" || atom_id == "O" || atom_id == "H" ||
        atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }

    // C-terminal carboxylate (Phe in CCD has OXT + HXT i.e. protonated):
    // tag OXT/HXT/C as Carboxylate when both Os are present. PHE in CCD
    // is the protonated form so we use AromaticHydroxyl-style PolarH for
    // the OH proton.
    if (atom_id == "OXT" || atom_id == "C" || atom_id == "O") {
        if (s.planar_group == PlanarGroupKind::None) s.planar_group = PlanarGroupKind::Carboxylate;
    }

    // Pseudoatom membership: PHE has QB (Hβ pair), QD (Hδ ring pair),
    // QE (Hε ring pair), QR (all ring Hs).
    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q, /*locant*/ static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    if (atom_id == "HD1" || atom_id == "HD2") {
        s.pseudoatom = {PseudoatomKind::Q, /*locant*/ static_cast<uint8_t>(Locant::Delta), 0, /*super*/ true};
    }
    if (atom_id == "HE1" || atom_id == "HE2") {
        s.pseudoatom = {PseudoatomKind::Q, /*locant*/ static_cast<uint8_t>(Locant::Epsilon), 0, /*super*/ true};
    }
    if (atom_id == "HZ") {
        s.pseudoatom = {PseudoatomKind::R, /*locant*/ static_cast<uint8_t>(Locant::Zeta), 0, /*super*/ true};
    }

    return s;
}


// PRO -- proline. Reference Section 3 (PRO block).
// Saturated 5-ring (pyrrolidine): N + Cα + Cβ + Cγ + Cδ.
// Special: NO backbone H (secondary amine, in-ring N).
// Cα is in the ring AND has Locant::None (orthogonality, dependencies
// §D.2). RingMembership.aromatic=false; RingMembership.planar=false
// (saturated; ring puckering is conformation-side).
//
// Per-atom RingPositionLabel uses the chemistry-distinct Pro labels
// (ProRingNitrogen / ProRingAlphaCarbon / ProRingBeta /
// ProRingPuckerPivot / ProRingDelta) rather than the generic
// `Saturated` label. The five labels surface puckering chemistry
// (Cγ as endo/exo pivot driving the 1-2 ppm Cβ/Cδ exo/endo shift
// difference; Vega & Boyer 1979, Schubert et al. 2002), backbone
// chemistry (N as in-ring secondary amine), and the canonical
// orthogonality case (Cα in ring with Locant::None).
//
// Hydrogens attached to ring carbons inherit the parent C's typed
// label per the existing convention used for aromatic rings.
SynthesisedFields SynthesisedForPro(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β/γ/δ methylenes

    auto set_ring_pyrrolidine = [&](RingPositionLabel pos) {
        s.ring_pos.primary.ring          = RingSystemKind::Pyrrolidine_Pro;
        s.ring_pos.primary.position      = pos;
        s.ring_pos.primary.ring_size     = 5;
        s.ring_pos.primary.aromatic      = false;
        s.ring_pos.primary.planar        = false;
        s.ring_pos.primary.n_heteroatoms = 1;
    };

    // Pro N: in ring, secondary amine (no H in the ring), planar_group = None.
    if (atom_id == "N") set_ring_pyrrolidine(RingPositionLabel::ProRingNitrogen);
    // Pro Cα + HA: in ring, Locant::None per the orthogonality rule.
    if (atom_id == "CA" || atom_id == "HA") {
        set_ring_pyrrolidine(RingPositionLabel::ProRingAlphaCarbon);
    }
    // Backbone C and O still in PeptideAmide (next residue's amide partner).
    if (atom_id == "C" || atom_id == "O") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    // Cβ + HB2/HB3: sidechain methylene, distinct from the Cγ pucker pivot.
    if (atom_id == "CB" || atom_id == "HB2" || atom_id == "HB3") {
        set_ring_pyrrolidine(RingPositionLabel::ProRingBeta);
    }
    // Cγ + HG2/HG3: the puckering pivot atom (endo/exo motion).
    if (atom_id == "CG" || atom_id == "HG2" || atom_id == "HG3") {
        set_ring_pyrrolidine(RingPositionLabel::ProRingPuckerPivot);
    }
    // Cδ + HD2/HD3: methylene bonded directly to the ring N (chemistry-
    // distinct from Cβ because of N proximity).
    if (atom_id == "CD" || atom_id == "HD2" || atom_id == "HD3") {
        set_ring_pyrrolidine(RingPositionLabel::ProRingDelta);
    }
    // Methylene Hs: QB / QG / QD pseudoatoms.
    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    if (atom_id == "HG2" || atom_id == "HG3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Gamma), 0, false};
    }
    if (atom_id == "HD2" || atom_id == "HD3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Delta), 0, false};
    }
    return s;
}


// SER -- serine. Reference Section 3 (SER block).
// β methylene QB; γ-OH (HydroxylOH_Aliphatic).
SynthesisedFields SynthesisedForSer(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    s.prochiral = MarkleyMethyleneProchiral(atom_id);  // β methylene

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }
    if (atom_id == "HG") s.polar_h = PolarHKind::HydroxylOH_Aliphatic;
    return s;
}


// THR -- threonine. Reference Section 3 (THR block).
// β-CH (single Hβ); γ1-OH (HydroxylOH_Aliphatic on HG1); γ2-methyl MG
// (HG21/22/23).
SynthesisedFields SynthesisedForThr(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HG1") s.polar_h = PolarHKind::HydroxylOH_Aliphatic;
    // γ2-methyl: HG21/22/23 form MG (Markley: methyl Hs).
    if (atom_id == "HG21" || atom_id == "HG22" || atom_id == "HG23") {
        s.pseudoatom = {PseudoatomKind::M,
                        static_cast<uint8_t>(Locant::Gamma), 2, false};
    }
    return s;
}


// TRP -- tryptophan. Reference Section 3 (TRP block).
// Fused indole: pyrrole 5-ring (CG, CD1, NE1, CE2, CD2) + benzene
// 6-ring (CD2, CE2, CE3, CZ2, CZ3, CH2). Bridgeheads (CD2, CE2)
// belong to BOTH rings: primary = 5-ring (smaller-ring convention),
// secondary = 6-ring.
// Trp 6-ring perimeter labels (Section 2 + dependencies §C.3): Ortho1
// = CE3, Ortho2 = CZ2, Meta1 = CZ3, Meta2 = CH2 (synthesised; Markley
// does not publish ipso/ortho/meta/para for indole 6-ring).
// IndoleNH on HE1.
SynthesisedFields SynthesisedForTrp(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    // β methylene Hs: HB2=ProS, HB3=ProR. Ring HE3 (locant E, digit
    // 3) would be wrongly tagged ProR by the helper; HE3 is an
    // aromatic CH on the indole 6-ring. Override below.
    s.prochiral = MarkleyMethyleneProchiral(atom_id);
    if (atom_id == "HE3") s.prochiral = ProchiralStereo::NotProchiral;

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }

    // Pyrrole 5-ring atoms.
    auto set_5ring = [&](RingPositionLabel pos) {
        s.ring_pos.primary.ring          = RingSystemKind::Indole_Trp_5;
        s.ring_pos.primary.position      = pos;
        s.ring_pos.primary.ring_size     = 5;
        s.ring_pos.primary.aromatic      = true;
        s.ring_pos.primary.planar        = true;
        s.ring_pos.primary.n_heteroatoms = 1;
        s.planar_group = PlanarGroupKind::Aromatic5Ring;
    };
    // Benzene 6-ring perimeter atoms (non-bridgehead).
    auto set_6ring = [&](RingPositionLabel pos) {
        s.ring_pos.primary.ring          = RingSystemKind::Indole_Trp_6;
        s.ring_pos.primary.position      = pos;
        s.ring_pos.primary.ring_size     = 6;
        s.ring_pos.primary.aromatic      = true;
        s.ring_pos.primary.planar        = true;
        s.ring_pos.primary.n_heteroatoms = 0;
        s.planar_group = PlanarGroupKind::Aromatic6Ring;
    };
    // Bridgeheads: in 5-ring (primary, smaller) AND 6-ring (secondary).
    auto set_bridgehead = [&]() {
        s.ring_pos.primary.ring          = RingSystemKind::Indole_Trp_5;
        s.ring_pos.primary.position      = RingPositionLabel::BridgeFusion;
        s.ring_pos.primary.ring_size     = 5;
        s.ring_pos.primary.aromatic      = true;
        s.ring_pos.primary.planar        = true;
        s.ring_pos.primary.n_heteroatoms = 1;
        s.ring_pos.secondary.ring          = RingSystemKind::Indole_Trp_6;
        s.ring_pos.secondary.position      = RingPositionLabel::BridgeFusion;
        s.ring_pos.secondary.ring_size     = 6;
        s.ring_pos.secondary.aromatic      = true;
        s.ring_pos.secondary.planar        = true;
        s.ring_pos.secondary.n_heteroatoms = 0;
        s.planar_group = PlanarGroupKind::Aromatic5Ring;  // primary
    };

    // 5-ring non-bridgehead atoms: CG (ipso), CD1 (PyrroleBeta), HD1,
    // NE1 (Heteroatom_NH), HE1.
    if      (atom_id == "CG")                       set_5ring(RingPositionLabel::Ipso);
    else if (atom_id == "CD1" || atom_id == "HD1")  set_5ring(RingPositionLabel::PyrroleBeta);
    else if (atom_id == "NE1" || atom_id == "HE1")  set_5ring(RingPositionLabel::Heteroatom_NH);
    // 5-ring + 6-ring bridgeheads.
    else if (atom_id == "CE2" || atom_id == "CD2")  set_bridgehead();
    // 6-ring perimeter atoms.
    else if (atom_id == "CE3" || atom_id == "HE3")  set_6ring(RingPositionLabel::Ortho1);
    else if (atom_id == "CZ2" || atom_id == "HZ2")  set_6ring(RingPositionLabel::Ortho2);
    else if (atom_id == "CZ3" || atom_id == "HZ3")  set_6ring(RingPositionLabel::Meta1);
    else if (atom_id == "CH2" || atom_id == "HH2")  set_6ring(RingPositionLabel::Meta2);

    // Indole 9-atom perimeter (Case 1995, J. Biomol. NMR 6, 341-346):
    // a third ring system covering all 9 heavy atoms of the conjugated
    // π current circuit. Encoded in the tertiary RingMembership slot.
    // The 9-atom perimeter set is the UNION of the 5-ring and 6-ring
    // heavy atoms (the bridgeheads CE2 and CD2 ARE in the perimeter,
    // not removed); what's "missing" relative to the dual-ring walk
    // is the shared bridgehead-fusion edge between CE2 and CD2 — the
    // perimeter wraps the outside of the indole, not crossing the
    // fusion. Ring-attached H atoms inherit perimeter membership per
    // the convention used for primary/secondary slots.
    auto set_perimeter_tertiary = [&]() {
        s.ring_pos.tertiary.ring          = RingSystemKind::Indole_Trp_9;
        s.ring_pos.tertiary.position      = RingPositionLabel::PerimeterMember;
        s.ring_pos.tertiary.ring_size     = 9;
        s.ring_pos.tertiary.aromatic      = true;
        s.ring_pos.tertiary.planar        = true;
        s.ring_pos.tertiary.n_heteroatoms = 1;  // NE1 is the sole heteroatom.
    };
    if (atom_id == "CG"  || atom_id == "CD1" || atom_id == "HD1" ||
        atom_id == "NE1" || atom_id == "HE1" || atom_id == "CE2" ||
        atom_id == "CZ2" || atom_id == "HZ2" || atom_id == "CH2" ||
        atom_id == "HH2" || atom_id == "CZ3" || atom_id == "HZ3" ||
        atom_id == "CE3" || atom_id == "HE3" || atom_id == "CD2") {
        set_perimeter_tertiary();
    }

    if (atom_id == "HE1") s.polar_h = PolarHKind::IndoleNH;
    return s;
}


// TYR -- tyrosine (default, neutral phenol form). Reference Section 3
// (TYR block). Aromatic 6-ring like PHE; CZ carries the para-OH.
// HH on the OH is HydroxylOH_Aromatic. Cζ-O-H rotation is conformation-
// side. QR super-group covers all four ring Hs (HD1/HD2/HE1/HE2), but
// not HH (HH is in AromaticHydroxyl planar group).
SynthesisedFields SynthesisedForTyr(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;
    // β methylene Hs: HB2=ProS, HB3=ProR. Ring Hs (HD1/HD2/HE1/HE2)
    // are caught by the helper as digit '1'/'2' but only HD2/HE2
    // (digit 2) match -> would be wrongly tagged ProS. Override below.
    s.prochiral = MarkleyMethyleneProchiral(atom_id);
    if (atom_id == "HD1" || atom_id == "HD2" ||
        atom_id == "HE1" || atom_id == "HE2") {
        s.prochiral = ProchiralStereo::NotProchiral;
    }

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    if (atom_id == "HB2" || atom_id == "HB3") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Beta), 0, false};
    }

    auto set_ring = [&](RingPositionLabel pos) {
        s.ring_pos.primary.ring          = RingSystemKind::Benzene_Tyr;
        s.ring_pos.primary.position      = pos;
        s.ring_pos.primary.ring_size     = 6;
        s.ring_pos.primary.aromatic      = true;
        s.ring_pos.primary.planar        = true;
        s.ring_pos.primary.n_heteroatoms = 0;
        s.planar_group = PlanarGroupKind::Aromatic6Ring;
    };
    if      (atom_id == "CG")                       set_ring(RingPositionLabel::Ipso);
    else if (atom_id == "CD1" || atom_id == "HD1")  set_ring(RingPositionLabel::Ortho1);
    else if (atom_id == "CD2" || atom_id == "HD2")  set_ring(RingPositionLabel::Ortho2);
    else if (atom_id == "CE1" || atom_id == "HE1")  set_ring(RingPositionLabel::Meta1);
    else if (atom_id == "CE2" || atom_id == "HE2")  set_ring(RingPositionLabel::Meta2);
    else if (atom_id == "CZ")                       set_ring(RingPositionLabel::Para);

    // Pseudoatom membership: QD / QE / QR.
    if (atom_id == "HD1" || atom_id == "HD2") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Delta), 0, true};
    }
    if (atom_id == "HE1" || atom_id == "HE2") {
        s.pseudoatom = {PseudoatomKind::Q,
                        static_cast<uint8_t>(Locant::Epsilon), 0, true};
    }

    // Para-OH: AromaticHydroxyl planar group spans Cζ-O-H. Note that
    // CZ retains its Aromatic6Ring assignment (set above); only OH and
    // HH belong to the AromaticHydroxyl group exclusively.
    if (atom_id == "OH") s.planar_group = PlanarGroupKind::AromaticHydroxyl;
    if (atom_id == "HH") {
        s.planar_group = PlanarGroupKind::AromaticHydroxyl;
        s.polar_h      = PolarHKind::HydroxylOH_Aromatic;
    }
    return s;
}


// TYM -- deprotonated tyrosine (aryloxide) variant. Reference Section
// 3 (TYM block) + dependencies §E.8. HH removed at the pipeline
// level; OH planar_group AromaticHydroxyl → AromaticOxide (per the
// new PlanarGroupKind::AromaticOxide enum, 2026-05-05). Oη formal
// charge: CCD TYR has 0 on OH (neutral phenol); TYM is phenolate so
// override OH 0 → -1.
SynthesisedFields SynthesisedForTym(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s = SynthesisedForTyr(atom_id, parsed);
    // Oη: AromaticHydroxyl -> AromaticOxide (phenolate; π-conjugation
    // qualitatively different from neutral phenol-OH).
    if (atom_id == "OH") s.planar_group = PlanarGroupKind::AromaticOxide;
    // HH removed at pipeline level; defensive zeroing here.
    if (atom_id == "HH") {
        s.planar_group = PlanarGroupKind::None;
        s.polar_h      = PolarHKind::NotPolar;
    }
    if (atom_id == "OH") s.charge_override = static_cast<int8_t>(-1);
    return s;
}


// VAL -- valine. Reference Section 3 (VAL block).
// Branched aliphatic: β-CH (single Hβ); γ1-methyl (Cγ1 = ProR);
// γ2-methyl (Cγ2 = ProS); MG1 + MG2 both in QG super (all six γ-methyl
// Hs collapse to QG).
SynthesisedFields SynthesisedForVal(const std::string& atom_id, const ParsedName& /*parsed*/) {
    SynthesisedFields s;

    if (atom_id == "N" || atom_id == "C" || atom_id == "O" ||
        atom_id == "H" || atom_id == "HN") {
        s.planar_group = PlanarGroupKind::PeptideAmide;
    }
    if (atom_id == "H" || atom_id == "HN") s.polar_h = PolarHKind::BackboneAmide;

    // Branched γ-methyls: Cγ1 = ProR; Cγ2 = ProS (Markley Fig 1).
    // Methyl Hs themselves are NotProchiral; only the parent C atoms
    // carry the prochiral label.
    if (atom_id == "CG1") s.prochiral = ProchiralStereo::ProR;
    if (atom_id == "CG2") s.prochiral = ProchiralStereo::ProS;

    // γ1-methyl: MG1 in QG super.
    if (atom_id == "HG11" || atom_id == "HG12" || atom_id == "HG13") {
        s.pseudoatom = {PseudoatomKind::M,
                        static_cast<uint8_t>(Locant::Gamma), 1, true};
    }
    // γ2-methyl: MG2 in QG super.
    if (atom_id == "HG21" || atom_id == "HG22" || atom_id == "HG23") {
        s.pseudoatom = {PseudoatomKind::M,
                        static_cast<uint8_t>(Locant::Gamma), 2, true};
    }
    return s;
}


// ============================================================================
// AtomSemanticEntry -- the typed record for one atom in the residue.
// Combines mechanical (parser), algorithmic (RDKit), and synthesised
// (per-residue table) sources into a single typed record with
// provenance witnesses.
// ============================================================================

struct AtomSemanticEntry {
    // Identification (in the generator scope only -- the runtime
    // table is keyed by typed structural matching against
    // AtomMechanicalIdentity (Element + Locant + Branch + DiIndex +
    // BackboneRole), per §H of
    // spec/plan/bones/topology-encoding-dependencies-2026-05-05.md).
    std::string               atom_id_for_log;

    // Mechanical-identity fields (the typed lookup key).
    Element                   element           = Element::Unknown;
    Locant                    locant            = Locant::None;
    BranchAddress             branch            = {};
    DiastereotopicIndex       di_index          = DiastereotopicIndex::None;
    BackboneRole              backbone_role     = BackboneRole::None;

    // Chemistry-substrate fields.
    ProchiralStereo           prochiral         = ProchiralStereo::NotProchiral;
    PlanarGroupKind           planar_group      = PlanarGroupKind::None;
    PlanarStereo              planar_stereo     = PlanarStereo::NotApplicable;
    PseudoatomMembership      pseudoatom        = {};
    PolarHKind                polar_h           = PolarHKind::NotPolar;
    RingPosition              ring_position     = {};
    bool                      aromatic          = false;
    Hybridisation             hybridisation     = Hybridisation::Unassigned;
    BondOrderMask             bond_orders       = {};
    int8_t                    formal_charge     = 0;
    bool                      is_exchangeable   = false;
    uint8_t                   equivalence_class = 0;

    // Provenance per field. Keyed by a small enum; populated as
    // each source contributes.
    SemanticProvenance        prov_locant;
    SemanticProvenance        prov_prochiral;
    SemanticProvenance        prov_planar_group;
    SemanticProvenance        prov_polar_h;
    SemanticProvenance        prov_ring_position;
    SemanticProvenance        prov_aromatic;
    SemanticProvenance        prov_hybridisation;
    // (Other provenance fields default-constructed; expanded as
    // additional sources are added in subsequent commits.)
};

// Combine sources with the precedence-table policy.
AtomSemanticEntry BuildAtomSemanticEntry(const std::string& atom_id,
                                          const ParsedName& parsed,
                                          unsigned int rdkit_idx,
                                          const RdkitFacts& facts,
                                          const SynthesisedFields& syn,
                                          int formal_charge_from_ccd,
                                          const std::string& ccd_type_symbol) {
    AtomSemanticEntry e;
    e.atom_id_for_log = atom_id;

    // Mechanical-identity fields. Element comes from the CCD record
    // (`_chem_comp_atom.type_symbol`) -- the parser's first-letter rule
    // is unsuitable for this slot because it conflates 'H' as an element
    // with 'H' as the eta-locant suffix. The CCD authority is canonical.
    // Locant + Branch + DiastereotopicIndex + BackboneRole come from the
    // parser. Together these five form the AtomMechanicalIdentity tuple
    // used by the generated `LookupBy` function. See §H of
    // `spec/plan/bones/topology-encoding-dependencies-2026-05-05.md`.
    e.element       = ElementFromSymbol(ccd_type_symbol);
    e.locant        = parsed.locant;
    e.branch        = parsed.branch;
    e.di_index      = parsed.di_index;
    e.backbone_role = parsed.backbone_role;
    e.prov_locant.witnesses[0] = {SemanticSource::IUPAC_1969,
                                  static_cast<uint8_t>(parsed.locant)};
    e.prov_locant.confidence   = SemanticConfidence::AlgorithmAuthoritative;
    e.prov_locant.sources_agree = true;

    // Algorithmic fields from RDKit.
    e.aromatic      = facts.aromatic[rdkit_idx];
    e.hybridisation = facts.hybridisation[rdkit_idx];
    e.bond_orders   = facts.bond_orders[rdkit_idx];
    e.equivalence_class = facts.canonical_rank[rdkit_idx];
    // Formal charge: CCD per-atom value is the default authority. A
    // SynthesisedFor<Residue|Variant> patch may set syn.charge_override
    // to redirect a per-atom charge for chemistry that differs from
    // CCD. Variants use this for protonation deltas (CYM Sγ -1, LYN Nζ
    // 0, HIP Nε2 +1, TYM Oη -1, ARN HE removal). Standard 20 patches
    // also use this to enforce Lewis-localised conventions per
    // dependencies §C.1: ASP OD2=-1, GLU OE2=-1, ARG NE=+1 + NH2=0,
    // HIS ND1=0 (HIE-default per AmberLeapInput.cpp:29). Atoms not
    // overridden track CCD untouched.
    e.formal_charge = syn.charge_override.has_value()
        ? *syn.charge_override
        : static_cast<int8_t>(formal_charge_from_ccd);

    // Prochiral assignment combines RDKit (chiral centres) with the
    // synthesised Markley-derived per-residue table (prochiral Hs
    // + prochiral branch carbons). When the synthesised value is
    // non-None, it WINS -- RDKit CIPLabeler does not label prochiral
    // Hs in the CCD free-amino-acid form, so the synthesised table
    // is the authoritative source for those atoms.
    if (syn.prochiral != ProchiralStereo::NotProchiral) {
        e.prochiral = syn.prochiral;
        e.prov_prochiral.witnesses[0] = {SemanticSource::Markley1998_Figure1,
                                         static_cast<uint8_t>(syn.prochiral)};
        e.prov_prochiral.confidence    = SemanticConfidence::AlgorithmAuthoritative;
        e.prov_prochiral.cip_verified  = false;
        e.prov_prochiral.sources_agree = true;
    } else {
        e.prochiral = facts.cip_stereo[rdkit_idx];
        e.prov_prochiral.witnesses[0] = {SemanticSource::RDKit_CIPLabeler,
                                         static_cast<uint8_t>(facts.cip_stereo[rdkit_idx])};
        e.prov_prochiral.confidence   = (facts.cip_stereo[rdkit_idx] == ProchiralStereo::NotProchiral)
            ? SemanticConfidence::AlgorithmAuthoritative
            : SemanticConfidence::CIPVerified;
        e.prov_prochiral.cip_verified = (facts.cip_stereo[rdkit_idx] != ProchiralStereo::NotProchiral &&
                                          facts.cip_stereo[rdkit_idx] != ProchiralStereo::Unassigned);
        e.prov_prochiral.sources_agree = true;
    }

    e.prov_aromatic.witnesses[0] = {SemanticSource::RDKit_AromaticPerception,
                                    static_cast<uint8_t>(facts.aromatic[rdkit_idx] ? 1 : 0)};
    e.prov_aromatic.confidence   = SemanticConfidence::AlgorithmAuthoritative;
    e.prov_aromatic.sources_agree = true;

    e.prov_hybridisation.witnesses[0] = {SemanticSource::RDKit_Hybridisation,
                                          static_cast<uint8_t>(facts.hybridisation[rdkit_idx])};
    e.prov_hybridisation.confidence   = SemanticConfidence::AlgorithmAuthoritative;
    e.prov_hybridisation.sources_agree = true;

    // Synthesised fields from per-residue table.
    e.planar_group  = syn.planar_group;
    e.planar_stereo = syn.planar_stereo;
    e.polar_h       = syn.polar_h;
    e.pseudoatom    = syn.pseudoatom;
    e.ring_position = syn.ring_pos;
    e.is_exchangeable = (e.polar_h != PolarHKind::NotPolar);

    // DiastereotopicIndex docstring (`SemanticEnums.h:123`) defines None
    // as "atom is not part of a prochiral methylene pair." The mechanical
    // parser keys on the digit in the atom name (HB2 -> Position2 etc.)
    // without knowing whether the atom belongs to a methyl or methylene.
    // For methyl-pseudoatom members (Ala HB2/HB3, Val MG[12]'s methyl Hs,
    // etc.), the digit reflects the H ordering within the methyl, NOT a
    // prochiral pair. Clear di_index here so the field's semantics match
    // its docstring.
    if (e.pseudoatom.kind == PseudoatomKind::M) {
        e.di_index = DiastereotopicIndex::None;
    }

    e.prov_planar_group.witnesses[0] = {SemanticSource::SynthesizedFromChemistry,
                                         static_cast<uint8_t>(syn.planar_group)};
    e.prov_planar_group.confidence   = SemanticConfidence::AlgorithmAuthoritative;

    e.prov_polar_h.witnesses[0] = {SemanticSource::SynthesizedFromChemistry,
                                    static_cast<uint8_t>(syn.polar_h)};
    e.prov_polar_h.confidence   = SemanticConfidence::AlgorithmAuthoritative;

    e.prov_ring_position.witnesses[0] = {SemanticSource::SynthesizedFromChemistry,
                                          static_cast<uint8_t>(syn.ring_pos.primary.position)};
    e.prov_ring_position.confidence   = SemanticConfidence::AlgorithmAuthoritative;

    return e;
}


// Dispatch the chemistry-source plug-in for one (residue-or-variant
// code, atom_id) pair. Single source-of-truth for which
// SynthesisedFor<...> function applies. Standard 20 codes route to
// their per-residue function; AMBER variant codes
// (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN/ARN/TYM) route to their
// per-variant patch (which itself layers on top of the parent
// SynthesisedFor<Parent>).
//
// Returns default-constructed SynthesisedFields if the code is
// unknown -- the dispatch is permissive so unknown codes do not fail
// hard at this layer; the main loop owns the residue-set decision.
SynthesisedFields DispatchSynthesised(const std::string& code,
                                       const std::string& atom_id,
                                       const ParsedName& parsed) {
    // Standard 20.
    if      (code == "ALA") return SynthesisedForAla(atom_id, parsed);
    else if (code == "ARG") return SynthesisedForArg(atom_id, parsed);
    else if (code == "ASN") return SynthesisedForAsn(atom_id, parsed);
    else if (code == "ASP") return SynthesisedForAsp(atom_id, parsed);
    else if (code == "CYS") return SynthesisedForCys(atom_id, parsed);
    else if (code == "GLN") return SynthesisedForGln(atom_id, parsed);
    else if (code == "GLU") return SynthesisedForGlu(atom_id, parsed);
    else if (code == "GLY") return SynthesisedForGly(atom_id, parsed);
    else if (code == "HIS") return SynthesisedForHis(atom_id, parsed);
    else if (code == "ILE") return SynthesisedForIle(atom_id, parsed);
    else if (code == "LEU") return SynthesisedForLeu(atom_id, parsed);
    else if (code == "LYS") return SynthesisedForLys(atom_id, parsed);
    else if (code == "MET") return SynthesisedForMet(atom_id, parsed);
    else if (code == "PHE") return SynthesisedForPhe(atom_id, parsed);
    else if (code == "PRO") return SynthesisedForPro(atom_id, parsed);
    else if (code == "SER") return SynthesisedForSer(atom_id, parsed);
    else if (code == "THR") return SynthesisedForThr(atom_id, parsed);
    else if (code == "TRP") return SynthesisedForTrp(atom_id, parsed);
    else if (code == "TYR") return SynthesisedForTyr(atom_id, parsed);
    else if (code == "VAL") return SynthesisedForVal(atom_id, parsed);
    // AMBER protonation variants.
    else if (code == "HID") return SynthesisedForHid(atom_id, parsed);
    else if (code == "HIE") return SynthesisedForHie(atom_id, parsed);
    else if (code == "HIP") return SynthesisedForHip(atom_id, parsed);
    else if (code == "ASH") return SynthesisedForAsh(atom_id, parsed);
    else if (code == "GLH") return SynthesisedForGlh(atom_id, parsed);
    else if (code == "CYX") return SynthesisedForCyx(atom_id, parsed);
    else if (code == "CYM") return SynthesisedForCym(atom_id, parsed);
    else if (code == "LYN") return SynthesisedForLyn(atom_id, parsed);
    else if (code == "ARN") return SynthesisedForArn(atom_id, parsed);
    else if (code == "TYM") return SynthesisedForTym(atom_id, parsed);
    return SynthesisedFields{};
}


// Cap-atom names that always belong to a terminal-state table, never
// to a per-residue chain table. The CCD entries for the canonical
// free-amino-acid forms include these atoms (the CCD records the
// fully-protonated zwitterion or its protonated-carboxyl variant);
// the runtime composition layer reads them from the dedicated
// terminal-state tables (`kCapNtermCharged`, `kCapNtermNeutral`,
// `kCapCtermDeprotonated`, `kCapCtermProtonated`) only when a residue
// is actually at a chain end. Per §H.4 of
// `spec/plan/bones/topology-encoding-dependencies-2026-05-05.md`.
//
// The five names below are residue-independent: the chemistry of
// NTERM_CHARGED Hs (H1/H2/H3) is the same on every residue's N
// terminus, etc. Centralising them eliminates per-residue duplication
// and makes the chain table semantically clean.
const std::set<std::string>& CapAtomNames() {
    static const std::set<std::string> kCaps = {
        "OXT",   // C-terminal carboxylate / carboxylic acid second oxygen.
        "HXT",   // C-terminal carboxylic-acid OH proton.
        "H1",    // N-terminal first H (NTERM_CHARGED + NTERM_NEUTRAL).
        "H2",    // N-terminal second H (NTERM_CHARGED + NTERM_NEUTRAL).
        "H3",    // N-terminal third H (NTERM_CHARGED only).
    };
    return kCaps;
}


// Chain-atom whitelist per residue (Fix 4 / §H.10 of dependencies file).
//
// The substrate's per-residue tables must match AMBER's chain-residue
// inventory exactly: CCD's free-amino-acid form is a superset (e.g.
// CCD ASP carries the protonated HD2; CCD GLU carries HE2; CCD HIS
// carries HD1; CCD PRO carries the backbone amide H, but Pro is a
// secondary amine in chain context). Filter the CCD atom list against
// this whitelist for the standard 20.
//
// Source: parallel encoding of `AmberAminoAcidVariantTable` in
// `src/AminoAcidType.cpp::AMINO_ACID_TYPES`. Duplicated here (acceptable
// for build-time utility) because the generator cannot link
// AminoAcidType.cpp without pulling in significant runtime deps. Verify
// against AminoAcidType.cpp when a residue's atom inventory changes.
//
// Variant tables (kHisAtoms_HID etc.) are NOT filtered here; they retain
// their full atom set. The variant-specific atoms (HD1 in HID/HIP, HD2
// in ASH, HE2 in GLH) are intentionally excluded from the chain table
// and intentionally included in the variant table.
const std::set<std::string>& ChainAtomWhitelist(const std::string& residue_3letter) {
    // ALA: 10 atoms.
    static const std::set<std::string> kALA = {
        "N","CA","C","O","H","HA","CB","HB1","HB2","HB3"
    };
    // ARG: 24 atoms (charged form, includes HE).
    static const std::set<std::string> kARG = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","HG2","HG3","CD","HD2","HD3",
        "NE","HE","CZ","NH1","HH11","HH12","NH2","HH21","HH22"
    };
    // ASN: 14 atoms.
    static const std::set<std::string> kASN = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","OD1","ND2","HD21","HD22"
    };
    // ASP: 12 atoms (deprotonated; drops HD2 from CCD's protonated form).
    static const std::set<std::string> kASP = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","OD1","OD2"
    };
    // CYS: 11 atoms (reduced thiol).
    static const std::set<std::string> kCYS = {
        "N","CA","C","O","H","HA","CB","HB2","HB3","SG","HG"
    };
    // GLN: 17 atoms.
    static const std::set<std::string> kGLN = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","HG2","HG3","CD","OE1","NE2","HE21","HE22"
    };
    // GLU: 15 atoms (deprotonated; drops HE2 from CCD's protonated form).
    static const std::set<std::string> kGLU = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","HG2","HG3","CD","OE1","OE2"
    };
    // GLY: 7 atoms.
    static const std::set<std::string> kGLY = {
        "N","CA","C","O","H","HA2","HA3"
    };
    // HIS: 17 atoms (HIE-default; drops HD1 from CCD's imidazolium form).
    static const std::set<std::string> kHIS = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","ND1","CD2","HD2","CE1","HE1","NE2","HE2"
    };
    // ILE: 19 atoms.
    static const std::set<std::string> kILE = {
        "N","CA","C","O","H","HA",
        "CB","HB","CG1","HG12","HG13","CG2","HG21","HG22","HG23",
        "CD1","HD11","HD12","HD13"
    };
    // LEU: 19 atoms.
    static const std::set<std::string> kLEU = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","HG",
        "CD1","HD11","HD12","HD13","CD2","HD21","HD22","HD23"
    };
    // LYS: 22 atoms.
    static const std::set<std::string> kLYS = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","HG2","HG3","CD","HD2","HD3",
        "CE","HE2","HE3","NZ","HZ1","HZ2","HZ3"
    };
    // MET: 17 atoms.
    static const std::set<std::string> kMET = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","HG2","HG3","SD","CE","HE1","HE2","HE3"
    };
    // PHE: 20 atoms.
    static const std::set<std::string> kPHE = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","CD1","HD1","CD2","HD2",
        "CE1","HE1","CE2","HE2","CZ","HZ"
    };
    // PRO: 14 atoms (drops backbone H; Pro is a secondary amine).
    static const std::set<std::string> kPRO = {
        "N","CA","C","O","HA",
        "CB","HB2","HB3","CG","HG2","HG3","CD","HD2","HD3"
    };
    // SER: 11 atoms.
    static const std::set<std::string> kSER = {
        "N","CA","C","O","H","HA","CB","HB2","HB3","OG","HG"
    };
    // THR: 14 atoms.
    static const std::set<std::string> kTHR = {
        "N","CA","C","O","H","HA",
        "CB","HB","OG1","HG1","CG2","HG21","HG22","HG23"
    };
    // TRP: 24 atoms.
    static const std::set<std::string> kTRP = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","CD1","HD1","CD2",
        "NE1","HE1","CE2","CE3","HE3",
        "CZ2","HZ2","CZ3","HZ3","CH2","HH2"
    };
    // TYR: 21 atoms.
    static const std::set<std::string> kTYR = {
        "N","CA","C","O","H","HA",
        "CB","HB2","HB3","CG","CD1","HD1","CD2","HD2",
        "CE1","HE1","CE2","HE2","CZ","OH","HH"
    };
    // VAL: 16 atoms.
    static const std::set<std::string> kVAL = {
        "N","CA","C","O","H","HA",
        "CB","HB","CG1","HG11","HG12","HG13","CG2","HG21","HG22","HG23"
    };
    static const std::set<std::string> kEMPTY;

    if (residue_3letter == "ALA") return kALA;
    if (residue_3letter == "ARG") return kARG;
    if (residue_3letter == "ASN") return kASN;
    if (residue_3letter == "ASP") return kASP;
    if (residue_3letter == "CYS") return kCYS;
    if (residue_3letter == "GLN") return kGLN;
    if (residue_3letter == "GLU") return kGLU;
    if (residue_3letter == "GLY") return kGLY;
    if (residue_3letter == "HIS") return kHIS;
    if (residue_3letter == "ILE") return kILE;
    if (residue_3letter == "LEU") return kLEU;
    if (residue_3letter == "LYS") return kLYS;
    if (residue_3letter == "MET") return kMET;
    if (residue_3letter == "PHE") return kPHE;
    if (residue_3letter == "PRO") return kPRO;
    if (residue_3letter == "SER") return kSER;
    if (residue_3letter == "THR") return kTHR;
    if (residue_3letter == "TRP") return kTRP;
    if (residue_3letter == "TYR") return kTYR;
    if (residue_3letter == "VAL") return kVAL;
    return kEMPTY;
}


// Build the per-atom semantic entries for one residue. The
// `synthesis_code` is the dispatch key for `DispatchSynthesised`; it
// is normally `res.three_letter` for the standard 20 but for variants
// it is the variant code (HID/HIE/HIP/...). The `atoms_to_remove`
// set names CCD-parent atoms that should be dropped from the variant
// output (e.g. HE2 for HID; HZ1 for LYN; HG for CYX/CYM; HH for TYM;
// HE for ARN). Empty set means emit every CCD atom (standard 20
// path).
//
// Per §H.4 of the dependencies file, cap atoms (OXT, HXT, H1, H2, H3)
// are partitioned out of the per-residue table irrespective of variant
// -- they live only in the terminal-state tables emitted separately.
// This means the standard 20 tables are SMALLER than they were before
// the partitioning landed (drop 2-4 atoms each).
//
// Per §H.10 (chain-atom whitelist filter, landed by Fix 4 of the
// 2026-05-05 OpenAI evaluation): the standard 20 tables are filtered
// against the chain-residue inventory in `ChainAtomWhitelist`. CCD's
// free-amino-acid form may carry extra atoms (HD2 in ASP, HE2 in GLU,
// HD1 in HIS, backbone H in PRO) that are NOT in AMBER's chain-residue
// inventory; the whitelist drops them. Empty whitelist (variant path)
// disables filtering.
std::vector<AtomSemanticEntry> BuildResidueEntries(
        const CcdResidue& res,
        const RdkitMolWithMap& built,
        const RdkitFacts& facts,
        const std::string& synthesis_code,
        const std::set<std::string>& atoms_to_remove,
        const std::set<std::string>& chain_whitelist,
        ProcessLog& log) {
    log.Section(std::string("build-entries::") + res.three_letter +
                (synthesis_code != res.three_letter
                    ? "::variant=" + synthesis_code
                    : std::string()));
    std::vector<AtomSemanticEntry> entries;
    entries.reserve(res.atoms.size());

    auto parents = BuildHydrogenParentMap(res);
    std::map<std::string, unsigned int> name_to_rdkit_idx;
    for (unsigned int i = 0; i < built.rdkit_idx_to_ccd_atom_id.size(); ++i) {
        name_to_rdkit_idx[built.rdkit_idx_to_ccd_atom_id[i]] = i;
    }

    int dropped = 0;
    int dropped_caps = 0;
    int dropped_non_chain = 0;
    const auto& cap_names = CapAtomNames();
    const bool apply_whitelist = !chain_whitelist.empty();
    for (const auto& ccd_atom : res.atoms) {
        if (atoms_to_remove.count(ccd_atom.atom_id) > 0) {
            ++dropped;
            log.KV("variant_dropped_atom", ccd_atom.atom_id);
            continue;
        }
        if (cap_names.count(ccd_atom.atom_id) > 0) {
            // Per §H.4 partitioning: cap atoms never live in a
            // per-residue table; they live in the four terminal-state
            // cap tables emitted separately at the end of the run.
            ++dropped_caps;
            log.KV("partitioned_cap_atom", ccd_atom.atom_id);
            continue;
        }
        if (apply_whitelist && chain_whitelist.count(ccd_atom.atom_id) == 0) {
            // Per §H.10 chain-atom whitelist: CCD-only atoms not in
            // AMBER's chain-residue inventory (HD2 in ASP, HE2 in GLU,
            // HD1 in HIS, backbone H in PRO) are dropped from the per-
            // residue table. They remain in the variant tables where
            // appropriate.
            ++dropped_non_chain;
            log.KV("non_chain_atom_filtered", ccd_atom.atom_id);
            continue;
        }
        std::string parent_name;
        if (ccd_atom.type_symbol == "H") {
            auto it = parents.find(ccd_atom.atom_id);
            if (it != parents.end()) parent_name = it->second;
        }
        const ParsedName parsed = ParseAtomName(ccd_atom.atom_id, parent_name);

        unsigned int rdkit_idx = 0;
        auto it = name_to_rdkit_idx.find(ccd_atom.atom_id);
        if (it != name_to_rdkit_idx.end()) rdkit_idx = it->second;

        SynthesisedFields syn = DispatchSynthesised(synthesis_code,
                                                     ccd_atom.atom_id,
                                                     parsed);

        AtomSemanticEntry e = BuildAtomSemanticEntry(ccd_atom.atom_id, parsed,
                                                      rdkit_idx, facts, syn,
                                                      ccd_atom.formal_charge,
                                                      ccd_atom.type_symbol);
        entries.push_back(std::move(e));
    }
    if (dropped > 0) log.KV("variant_atoms_dropped_total", dropped);
    if (dropped_caps > 0) log.KV("partitioned_cap_atoms_total", dropped_caps);
    if (dropped_non_chain > 0) log.KV("non_chain_atoms_filtered_total", dropped_non_chain);
    log.KV("entries_built", static_cast<int>(entries.size()));
    return entries;
}


// Spot-log a few entries for inspection.
void LogEntriesSpotCheck(const std::vector<AtomSemanticEntry>& entries,
                         ProcessLog& log,
                         const std::string& residue_3letter) {
    log.Section(std::string("entries-spot-check::") + residue_3letter);
    for (size_t i = 0; i < entries.size(); ++i) {
        const auto& e = entries[i];
        std::string fact;
        fact += "elem=" + std::to_string(static_cast<int>(e.element));
        fact += " loc=" + std::to_string(static_cast<int>(e.locant));
        fact += " branch=" + std::to_string(e.branch.outer) + "/" + std::to_string(e.branch.inner);
        fact += " di=" + std::to_string(static_cast<int>(e.di_index));
        fact += " bbrole=" + std::to_string(static_cast<int>(e.backbone_role));
        fact += " prochiral=" + std::to_string(static_cast<int>(e.prochiral));
        fact += " planar=" + std::to_string(static_cast<int>(e.planar_group));
        fact += " polarH=" + std::to_string(static_cast<int>(e.polar_h));
        fact += " ring_pos=" + std::to_string(static_cast<int>(e.ring_position.primary.position));
        fact += " arom=" + std::to_string(static_cast<int>(e.aromatic));
        fact += " hyb=" + std::to_string(static_cast<int>(e.hybridisation));
        fact += " eqcls=" + std::to_string(static_cast<int>(e.equivalence_class));
        log.KV(std::to_string(i) + "_" + e.atom_id_for_log, fact);
    }
}


// ============================================================================
// C++ emitter -- writes the typed-enum table file consumed at runtime.
//
// Format: a generated .cpp file with constexpr std::array of
// AtomSemanticTable per residue and per AMBER protonation variant,
// plus four cap tables (one per TerminalState) carrying the residue-
// independent terminal-state chemistry, plus two Lookup functions
// (`LookupBy` for per-(residue, variant) tables; `LookupCap` for cap
// tables) keyed on AtomMechanicalIdentity (Element + Locant +
// BranchAddress + DiastereotopicIndex + BackboneRole) -- typed
// structural matching, NOT atom_local_idx and NOT atom_id strings.
// Per §H of spec/plan/bones/topology-encoding-dependencies-2026-05-05.md.
// The output includes only typed-enum references; no std::string
// literals carrying chemistry data, no gemmi/RDKit/cifpp symbols, no
// Eigen.
//
// The string barrier from the generator's perspective: this Emit
// function is the only place where chemistry-derived enum names are
// converted back to text -- but only as C++ source code identifiers
// (e.g. "nmr::Locant::Beta"), not as data. The compiler then turns
// those identifier strings back into the typed-enum integer values
// at compile time. No runtime string survives.
// ============================================================================

const char* ElementLiteral(Element e) {
    switch (e) {
        case Element::H:       return "nmr::Element::H";
        case Element::C:       return "nmr::Element::C";
        case Element::N:       return "nmr::Element::N";
        case Element::O:       return "nmr::Element::O";
        case Element::S:       return "nmr::Element::S";
        case Element::Unknown: return "nmr::Element::Unknown";
    }
    return "nmr::Element::Unknown";
}

const char* BackboneRoleLiteral(BackboneRole r) {
    switch (r) {
        case BackboneRole::None:           return "nmr::BackboneRole::None";
        case BackboneRole::Nitrogen:       return "nmr::BackboneRole::Nitrogen";
        case BackboneRole::AlphaCarbon:    return "nmr::BackboneRole::AlphaCarbon";
        case BackboneRole::CarbonylCarbon: return "nmr::BackboneRole::CarbonylCarbon";
        case BackboneRole::CarbonylOxygen: return "nmr::BackboneRole::CarbonylOxygen";
        case BackboneRole::AmideHydrogen:  return "nmr::BackboneRole::AmideHydrogen";
        case BackboneRole::AlphaHydrogen:  return "nmr::BackboneRole::AlphaHydrogen";
    }
    return "nmr::BackboneRole::None";
}

const char* TerminalStateLiteral(TerminalState t) {
    switch (t) {
        case TerminalState::Internal:          return "nmr::TerminalState::Internal";
        case TerminalState::NtermCharged:      return "nmr::TerminalState::NtermCharged";
        case TerminalState::NtermNeutral:      return "nmr::TerminalState::NtermNeutral";
        case TerminalState::CtermDeprotonated: return "nmr::TerminalState::CtermDeprotonated";
        case TerminalState::CtermProtonated:   return "nmr::TerminalState::CtermProtonated";
    }
    return "nmr::TerminalState::Internal";
}

const char* LocantLiteral(Locant l) {
    switch (l) {
        case Locant::None:    return "nmr::Locant::None";
        case Locant::Alpha:   return "nmr::Locant::Alpha";
        case Locant::Beta:    return "nmr::Locant::Beta";
        case Locant::Gamma:   return "nmr::Locant::Gamma";
        case Locant::Delta:   return "nmr::Locant::Delta";
        case Locant::Epsilon: return "nmr::Locant::Epsilon";
        case Locant::Zeta:    return "nmr::Locant::Zeta";
        case Locant::Eta:     return "nmr::Locant::Eta";
    }
    return "nmr::Locant::None";
}

const char* DiastereotopicLiteral(DiastereotopicIndex d) {
    switch (d) {
        case DiastereotopicIndex::None:      return "nmr::DiastereotopicIndex::None";
        case DiastereotopicIndex::Position2: return "nmr::DiastereotopicIndex::Position2";
        case DiastereotopicIndex::Position3: return "nmr::DiastereotopicIndex::Position3";
    }
    return "nmr::DiastereotopicIndex::None";
}

const char* ProchiralLiteral(ProchiralStereo p) {
    switch (p) {
        case ProchiralStereo::NotProchiral: return "nmr::ProchiralStereo::NotProchiral";
        case ProchiralStereo::ProR:         return "nmr::ProchiralStereo::ProR";
        case ProchiralStereo::ProS:         return "nmr::ProchiralStereo::ProS";
        case ProchiralStereo::Unassigned:   return "nmr::ProchiralStereo::Unassigned";
    }
    return "nmr::ProchiralStereo::NotProchiral";
}

const char* PlanarGroupLiteral(PlanarGroupKind p) {
    switch (p) {
        case PlanarGroupKind::None:             return "nmr::PlanarGroupKind::None";
        case PlanarGroupKind::PeptideAmide:     return "nmr::PlanarGroupKind::PeptideAmide";
        case PlanarGroupKind::SidechainAmide:   return "nmr::PlanarGroupKind::SidechainAmide";
        case PlanarGroupKind::Guanidinium:      return "nmr::PlanarGroupKind::Guanidinium";
        case PlanarGroupKind::Imidazole:        return "nmr::PlanarGroupKind::Imidazole";
        case PlanarGroupKind::Aromatic6Ring:    return "nmr::PlanarGroupKind::Aromatic6Ring";
        case PlanarGroupKind::Aromatic5Ring:    return "nmr::PlanarGroupKind::Aromatic5Ring";
        case PlanarGroupKind::Carboxylate:      return "nmr::PlanarGroupKind::Carboxylate";
        case PlanarGroupKind::AromaticHydroxyl: return "nmr::PlanarGroupKind::AromaticHydroxyl";
        case PlanarGroupKind::AromaticOxide:    return "nmr::PlanarGroupKind::AromaticOxide";
    }
    return "nmr::PlanarGroupKind::None";
}

const char* PlanarStereoLiteral(PlanarStereo p) {
    switch (p) {
        case PlanarStereo::NotApplicable: return "nmr::PlanarStereo::NotApplicable";
        case PlanarStereo::E:             return "nmr::PlanarStereo::E";
        case PlanarStereo::Z:             return "nmr::PlanarStereo::Z";
        case PlanarStereo::Unspecified:   return "nmr::PlanarStereo::Unspecified";
    }
    return "nmr::PlanarStereo::NotApplicable";
}

const char* PseudoatomKindLiteral(PseudoatomKind k) {
    switch (k) {
        case PseudoatomKind::None: return "nmr::PseudoatomKind::None";
        case PseudoatomKind::M:    return "nmr::PseudoatomKind::M";
        case PseudoatomKind::Q:    return "nmr::PseudoatomKind::Q";
        case PseudoatomKind::R:    return "nmr::PseudoatomKind::R";
    }
    return "nmr::PseudoatomKind::None";
}

const char* PolarHLiteral(PolarHKind p) {
    switch (p) {
        case PolarHKind::NotPolar:              return "nmr::PolarHKind::NotPolar";
        case PolarHKind::BackboneAmide:         return "nmr::PolarHKind::BackboneAmide";
        case PolarHKind::SidechainPrimaryAmide: return "nmr::PolarHKind::SidechainPrimaryAmide";
        case PolarHKind::IndoleNH:              return "nmr::PolarHKind::IndoleNH";
        case PolarHKind::AmmoniumNH:            return "nmr::PolarHKind::AmmoniumNH";
        case PolarHKind::GuanidiniumNH:         return "nmr::PolarHKind::GuanidiniumNH";
        case PolarHKind::ImidazoleNH:           return "nmr::PolarHKind::ImidazoleNH";
        case PolarHKind::CarboxylOH:            return "nmr::PolarHKind::CarboxylOH";
        case PolarHKind::HydroxylOH_Aliphatic:  return "nmr::PolarHKind::HydroxylOH_Aliphatic";
        case PolarHKind::HydroxylOH_Aromatic:   return "nmr::PolarHKind::HydroxylOH_Aromatic";
        case PolarHKind::ThiolSH:               return "nmr::PolarHKind::ThiolSH";
        case PolarHKind::AmineNH:               return "nmr::PolarHKind::AmineNH";
        case PolarHKind::OtherPolarH:           return "nmr::PolarHKind::OtherPolarH";
    }
    return "nmr::PolarHKind::NotPolar";
}

const char* RingSystemLiteral(RingSystemKind r) {
    switch (r) {
        case RingSystemKind::NotInRing:       return "nmr::RingSystemKind::NotInRing";
        case RingSystemKind::Benzene_Phe:     return "nmr::RingSystemKind::Benzene_Phe";
        case RingSystemKind::Benzene_Tyr:     return "nmr::RingSystemKind::Benzene_Tyr";
        case RingSystemKind::Imidazole_His:   return "nmr::RingSystemKind::Imidazole_His";
        case RingSystemKind::Indole_Trp_5:    return "nmr::RingSystemKind::Indole_Trp_5";
        case RingSystemKind::Indole_Trp_6:    return "nmr::RingSystemKind::Indole_Trp_6";
        case RingSystemKind::Pyrrolidine_Pro: return "nmr::RingSystemKind::Pyrrolidine_Pro";
        case RingSystemKind::Indole_Trp_9:    return "nmr::RingSystemKind::Indole_Trp_9";
    }
    return "nmr::RingSystemKind::NotInRing";
}

const char* RingPositionLabelLiteral(RingPositionLabel r) {
    switch (r) {
        case RingPositionLabel::NotInRing:           return "nmr::RingPositionLabel::NotInRing";
        case RingPositionLabel::Ipso:                return "nmr::RingPositionLabel::Ipso";
        case RingPositionLabel::Ortho1:              return "nmr::RingPositionLabel::Ortho1";
        case RingPositionLabel::Ortho2:              return "nmr::RingPositionLabel::Ortho2";
        case RingPositionLabel::Meta1:               return "nmr::RingPositionLabel::Meta1";
        case RingPositionLabel::Meta2:               return "nmr::RingPositionLabel::Meta2";
        case RingPositionLabel::Para:                return "nmr::RingPositionLabel::Para";
        case RingPositionLabel::PyrroleAlpha:        return "nmr::RingPositionLabel::PyrroleAlpha";
        case RingPositionLabel::PyrroleBeta:         return "nmr::RingPositionLabel::PyrroleBeta";
        case RingPositionLabel::BridgeFusion:        return "nmr::RingPositionLabel::BridgeFusion";
        case RingPositionLabel::Heteroatom_NH:       return "nmr::RingPositionLabel::Heteroatom_NH";
        case RingPositionLabel::Heteroatom_NoH:      return "nmr::RingPositionLabel::Heteroatom_NoH";
        case RingPositionLabel::Heteroatom_OH:       return "nmr::RingPositionLabel::Heteroatom_OH";
        case RingPositionLabel::Saturated:           return "nmr::RingPositionLabel::Saturated";
        case RingPositionLabel::ProRingNitrogen:     return "nmr::RingPositionLabel::ProRingNitrogen";
        case RingPositionLabel::ProRingAlphaCarbon:  return "nmr::RingPositionLabel::ProRingAlphaCarbon";
        case RingPositionLabel::ProRingBeta:         return "nmr::RingPositionLabel::ProRingBeta";
        case RingPositionLabel::ProRingPuckerPivot:  return "nmr::RingPositionLabel::ProRingPuckerPivot";
        case RingPositionLabel::ProRingDelta:        return "nmr::RingPositionLabel::ProRingDelta";
        case RingPositionLabel::PerimeterMember:     return "nmr::RingPositionLabel::PerimeterMember";
    }
    return "nmr::RingPositionLabel::NotInRing";
}

std::string RingMembershipLiteral(const RingMembership& rm) {
    std::ostringstream o;
    o << "{" << RingSystemLiteral(rm.ring)
      << ", " << RingPositionLabelLiteral(rm.position)
      << ", " << static_cast<int>(rm.ring_size)
      << ", " << (rm.aromatic ? "true" : "false")
      << ", " << (rm.planar ? "true" : "false")
      << ", " << static_cast<int>(rm.n_heteroatoms)
      << "}";
    return o.str();
}

std::string PseudoatomLiteral(const PseudoatomMembership& p) {
    std::ostringstream o;
    o << "{" << PseudoatomKindLiteral(p.kind)
      << ", " << static_cast<int>(p.locant)
      << ", " << static_cast<int>(p.branch)
      << ", " << (p.in_super_group ? "true" : "false")
      << "}";
    return o.str();
}

// Emit one AtomSemanticTable record as a brace-initialiser. Field
// order MUST match SemanticEnums.h's AtomSemanticTable definition.
// The first five fields (Element, Locant, BranchAddress,
// DiastereotopicIndex, BackboneRole) form the AtomMechanicalIdentity
// tuple used by the generated `LookupBy` function -- see §H of
// `spec/plan/bones/topology-encoding-dependencies-2026-05-05.md`.
std::string EmitEntryLiteral(const AtomSemanticEntry& e) {
    std::ostringstream o;
    o << "    { " << ElementLiteral(e.element)
      << ", " << LocantLiteral(e.locant)
      << ", {" << static_cast<int>(e.branch.outer)
      <<  "," << static_cast<int>(e.branch.inner) << "}"
      << ", " << DiastereotopicLiteral(e.di_index)
      << ", " << BackboneRoleLiteral(e.backbone_role)
      << ", " << ProchiralLiteral(e.prochiral)
      << ", " << PlanarGroupLiteral(e.planar_group)
      << ", " << PlanarStereoLiteral(e.planar_stereo)
      << ", " << PseudoatomLiteral(e.pseudoatom)
      << ", " << PolarHLiteral(e.polar_h)
      << ", {" << RingMembershipLiteral(e.ring_position.primary)
      <<  ", " << RingMembershipLiteral(e.ring_position.secondary)
      <<  ", " << RingMembershipLiteral(e.ring_position.tertiary) << "}"
      << ", " << (e.aromatic ? "true" : "false")
      << ", " << static_cast<int>(e.formal_charge)
      << ", " << (e.is_exchangeable ? "true" : "false")
      << ", " << static_cast<int>(e.equivalence_class)
      << " }";
    return o.str();
}


// Holds entries emitted for one (residue, variant) combination.
struct ResidueEntries {
    std::string residue_3letter;
    std::string variant_3letter;   ///< "" for canonical; "HID" / "ASH" / etc for variants
    std::vector<std::string>      atom_id_for_log;   ///< For comments alongside each row.
    std::vector<AtomSemanticEntry> entries;
};


// Holds the entries for one cap-table (terminal-state) emission. Cap
// atom rows are residue-independent: the chemistry of NTERM_CHARGED
// H1/H2/H3 is the same on every residue's N terminus.
struct CapTableEntries {
    TerminalState                  state;
    std::string                    array_name;       ///< "kCapNtermCharged" etc.
    std::vector<std::string>       atom_id_for_log;
    std::vector<AtomSemanticEntry> entries;
};


// Build one cap-atom AtomSemanticEntry. Cap chemistry is hand-encoded
// from §4 of `spec/plan/bones/topology-residue-reference-2026-05-05.md`.
// The mechanical-identity tuple (Element, Locant, Branch, DiIndex,
// BackboneRole) for cap atoms uses Locant::None and BackboneRole::None
// for terminus-added Hs / Os (they are neither sidechain Greek-letter
// atoms nor backbone peptide-amide-unit atoms). For backbone OVERRIDE
// rows (the N row in NTERM caps, the C and O rows in CTERM caps) the
// BackboneRole IS populated -- these rows compose with the per-residue
// table by matching the same (Element, Locant=None, Branch=0, DiIndex,
// BackboneRole) identity as the chain entry. See `MakeCapBackboneOverride`.
AtomSemanticEntry MakeCapAtomEntry(const std::string& atom_id,
                                    Element element,
                                    PlanarGroupKind planar_group,
                                    PolarHKind polar_h,
                                    int8_t formal_charge,
                                    bool is_exchangeable,
                                    PseudoatomMembership pseudoatom = {}) {
    AtomSemanticEntry e;
    e.atom_id_for_log = atom_id;
    e.element         = element;
    e.locant          = Locant::None;
    e.branch          = {};
    e.di_index        = DiastereotopicIndex::None;
    e.backbone_role   = BackboneRole::None;
    e.prochiral       = ProchiralStereo::NotProchiral;
    e.planar_group    = planar_group;
    e.planar_stereo   = PlanarStereo::NotApplicable;
    e.pseudoatom      = pseudoatom;
    e.polar_h         = polar_h;
    e.ring_position   = {};
    e.aromatic        = false;
    e.formal_charge   = formal_charge;
    e.is_exchangeable = is_exchangeable;
    e.equivalence_class = 0;
    return e;
}


// Build a cap-table OVERRIDE entry for a backbone atom whose chemistry
// differs at the terminus from the canonical internal-chain chemistry.
// These rows have a non-None BackboneRole so they share their mechanical
// identity (Element, Locant::None, Branch={0,0}, DiIndex::None,
// backbone_role) with the chain table's same-role entry; the runtime
// composition rule is "if LookupCap returns non-null, override the
// chain entry's fields with the cap entry's fields for that atom."
// See the comment block at the top of LookupCap in the generated .cpp
// for the full composition rule.
AtomSemanticEntry MakeCapBackboneOverride(const std::string& atom_id,
                                           Element element,
                                           BackboneRole backbone_role,
                                           PlanarGroupKind planar_group,
                                           PolarHKind polar_h,
                                           int8_t formal_charge,
                                           bool is_exchangeable) {
    AtomSemanticEntry e;
    e.atom_id_for_log = atom_id;
    e.element         = element;
    e.locant          = Locant::None;
    e.branch          = {};
    e.di_index        = DiastereotopicIndex::None;
    e.backbone_role   = backbone_role;
    e.prochiral       = ProchiralStereo::NotProchiral;
    e.planar_group    = planar_group;
    e.planar_stereo   = PlanarStereo::NotApplicable;
    e.pseudoatom      = {};
    e.polar_h         = polar_h;
    e.ring_position   = {};
    e.aromatic        = false;
    e.formal_charge   = formal_charge;
    e.is_exchangeable = is_exchangeable;
    e.equivalence_class = 0;
    return e;
}


// Synthesise the four cap tables per §4 of the reference doc.
//
// NTERM_CHARGED (NH3+):
//   N (override) -- planar_group=None (terminal N is not a peptide-amide
//     nitrogen); polarH=NotPolar (label is on the H, not the N); formal_chg=+1.
//   H1, H2, H3 -- ammonium NH; PolarHKind::AmmoniumNH; formal_chg=0;
//   exchangeable=true; Pseudoatom = Q with Locant::None (Markley Table 1
//   does not list a pseudoatom for terminus-specific Hs; Q-with-locant=None
//   is the Q-like equivalent-group encoding per the reference doc).
// NTERM_NEUTRAL (NH2):
//   N (override) -- planar_group=None; polarH=NotPolar; formal_chg=0.
//   H1, H2 -- amine NH; PolarHKind::AmineNH; formal_chg=0;
//   exchangeable=true.
// CTERM_DEPROTONATED (COO-):
//   C (override) -- planar_group=Carboxylate (was PeptideAmide on chain);
//     polarH=NotPolar; formal_chg=0.
//   O (override = O' in Markley convention) -- planar_group=Carboxylate
//     (was PeptideAmide on chain); polarH=NotPolar; formal_chg=0
//     (the -1 is on OXT/O'').
//   OXT (= O'') -- planar_group=Carboxylate; polarH=NotPolar;
//   formal_chg=-1; exchangeable=false.
// CTERM_PROTONATED (COOH):
//   C (override) -- planar_group=Carboxylate; polarH=NotPolar; formal_chg=0.
//   O (override) -- planar_group=Carboxylate; polarH=NotPolar; formal_chg=0.
//   OXT (= O'') -- planar_group=Carboxylate; polarH=NotPolar;
//   formal_chg=0; exchangeable=false.
//   HXT (= H'') -- planar_group=Carboxylate; polarH=CarboxylOH;
//   formal_chg=0; exchangeable=true.
//
// Backbone H removal in NTERM is handled by AminoAcidType (terminal
// residues use H1/H2/H3, not H) -- no sentinel needed in the cap table.
//
// Composition rule (documented at LookupCap site below): if
// LookupCap(state, identity) returns non-null, the cap entry overrides
// the per-residue entry for that atom. Backbone overrides have the
// same (Element, Locant::None, Branch={0,0}, DiIndex::None,
// backbone_role) identity as the chain row for the corresponding
// backbone slot, so they shadow it under that rule.
std::vector<CapTableEntries> BuildCapTables(ProcessLog& log) {
    log.Section("build-cap-tables");

    // Per-cap-H pseudoatom encoding (Markley Table 1 has no entry for
    // terminus-specific Hs; the Q/None/0/false encoding is the
    // closest-fit "equivalent group" per spec Section 4 + §A.3 of the
    // dependencies file).
    const PseudoatomMembership q_terminus = {PseudoatomKind::Q, 0, 0, false};

    std::vector<CapTableEntries> caps;

    {
        CapTableEntries c;
        c.state      = TerminalState::NtermCharged;
        c.array_name = "kCapNtermCharged";
        // Backbone N override: chain-N planar_group is PeptideAmide; the
        // terminal NH3+ is a primary amine, not a peptide-amide unit.
        // formal_charge moves to +1.
        c.entries.push_back(MakeCapBackboneOverride(
            "N", Element::N, BackboneRole::Nitrogen,
            PlanarGroupKind::None, PolarHKind::NotPolar,
            static_cast<int8_t>(+1), false));
        c.entries.push_back(MakeCapAtomEntry(
            "H1", Element::H, PlanarGroupKind::None,
            PolarHKind::AmmoniumNH, 0, true, q_terminus));
        c.entries.push_back(MakeCapAtomEntry(
            "H2", Element::H, PlanarGroupKind::None,
            PolarHKind::AmmoniumNH, 0, true, q_terminus));
        c.entries.push_back(MakeCapAtomEntry(
            "H3", Element::H, PlanarGroupKind::None,
            PolarHKind::AmmoniumNH, 0, true, q_terminus));
        for (const auto& e : c.entries) c.atom_id_for_log.push_back(e.atom_id_for_log);
        log.KV("kCapNtermCharged_size", static_cast<int>(c.entries.size()));
        caps.push_back(std::move(c));
    }

    {
        CapTableEntries c;
        c.state      = TerminalState::NtermNeutral;
        c.array_name = "kCapNtermNeutral";
        // Backbone N override: chain-N planar_group is PeptideAmide; the
        // terminal NH2 is a primary amine. formal_charge stays 0.
        c.entries.push_back(MakeCapBackboneOverride(
            "N", Element::N, BackboneRole::Nitrogen,
            PlanarGroupKind::None, PolarHKind::NotPolar,
            static_cast<int8_t>(0), false));
        c.entries.push_back(MakeCapAtomEntry(
            "H1", Element::H, PlanarGroupKind::None,
            PolarHKind::AmineNH, 0, true, q_terminus));
        c.entries.push_back(MakeCapAtomEntry(
            "H2", Element::H, PlanarGroupKind::None,
            PolarHKind::AmineNH, 0, true, q_terminus));
        for (const auto& e : c.entries) c.atom_id_for_log.push_back(e.atom_id_for_log);
        log.KV("kCapNtermNeutral_size", static_cast<int>(c.entries.size()));
        caps.push_back(std::move(c));
    }

    {
        CapTableEntries c;
        c.state      = TerminalState::CtermDeprotonated;
        c.array_name = "kCapCtermDeprotonated";
        // Backbone C and O overrides: chain planar_group is PeptideAmide;
        // the terminal COO- is a carboxylate. -1 is on OXT/O'' per the
        // CTERM_DEPROTONATED entry below, so this O (= O') stays formal_charge 0.
        c.entries.push_back(MakeCapBackboneOverride(
            "C", Element::C, BackboneRole::CarbonylCarbon,
            PlanarGroupKind::Carboxylate, PolarHKind::NotPolar,
            static_cast<int8_t>(0), false));
        c.entries.push_back(MakeCapBackboneOverride(
            "O", Element::O, BackboneRole::CarbonylOxygen,
            PlanarGroupKind::Carboxylate, PolarHKind::NotPolar,
            static_cast<int8_t>(0), false));
        c.entries.push_back(MakeCapAtomEntry(
            "OXT", Element::O, PlanarGroupKind::Carboxylate,
            PolarHKind::NotPolar, -1, false));
        for (const auto& e : c.entries) c.atom_id_for_log.push_back(e.atom_id_for_log);
        log.KV("kCapCtermDeprotonated_size", static_cast<int>(c.entries.size()));
        caps.push_back(std::move(c));
    }

    {
        CapTableEntries c;
        c.state      = TerminalState::CtermProtonated;
        c.array_name = "kCapCtermProtonated";
        // Backbone C and O overrides: chain planar_group PeptideAmide ->
        // Carboxylate at the protonated terminus. formal_charges stay 0
        // (no charged form here).
        c.entries.push_back(MakeCapBackboneOverride(
            "C", Element::C, BackboneRole::CarbonylCarbon,
            PlanarGroupKind::Carboxylate, PolarHKind::NotPolar,
            static_cast<int8_t>(0), false));
        c.entries.push_back(MakeCapBackboneOverride(
            "O", Element::O, BackboneRole::CarbonylOxygen,
            PlanarGroupKind::Carboxylate, PolarHKind::NotPolar,
            static_cast<int8_t>(0), false));
        c.entries.push_back(MakeCapAtomEntry(
            "OXT", Element::O, PlanarGroupKind::Carboxylate,
            PolarHKind::NotPolar, 0, false));
        c.entries.push_back(MakeCapAtomEntry(
            "HXT", Element::H, PlanarGroupKind::Carboxylate,
            PolarHKind::CarboxylOH, 0, true));
        for (const auto& e : c.entries) c.atom_id_for_log.push_back(e.atom_id_for_log);
        log.KV("kCapCtermProtonated_size", static_cast<int>(c.entries.size()));
        caps.push_back(std::move(c));
    }

    return caps;
}


// Format the residue's identifier for use as a C++ array name.
// Canonical -> kPheAtoms, variant -> kPheAtoms_HID.
std::string CppArrayName(const std::string& residue, const std::string& variant) {
    auto title_case = [](const std::string& s) {
        if (s.empty()) return s;
        std::string out = s;
        for (size_t i = 1; i < out.size(); ++i) {
            out[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(out[i])));
        }
        return out;
    };
    std::string n = "k" + title_case(residue) + "Atoms";
    if (!variant.empty()) n += "_" + variant;
    return n;
}


// Map an AMBER 3-letter code to its AminoAcid enum identifier (as a
// C++ source token). Used by the LookupBy emitter to dispatch on
// (residue, variant_idx).
const char* AminoAcidEnumLiteral(const std::string& three_letter) {
    if (three_letter == "ALA") return "nmr::AminoAcid::ALA";
    if (three_letter == "ARG") return "nmr::AminoAcid::ARG";
    if (three_letter == "ASN") return "nmr::AminoAcid::ASN";
    if (three_letter == "ASP") return "nmr::AminoAcid::ASP";
    if (three_letter == "CYS") return "nmr::AminoAcid::CYS";
    if (three_letter == "GLN") return "nmr::AminoAcid::GLN";
    if (three_letter == "GLU") return "nmr::AminoAcid::GLU";
    if (three_letter == "GLY") return "nmr::AminoAcid::GLY";
    if (three_letter == "HIS") return "nmr::AminoAcid::HIS";
    if (three_letter == "ILE") return "nmr::AminoAcid::ILE";
    if (three_letter == "LEU") return "nmr::AminoAcid::LEU";
    if (three_letter == "LYS") return "nmr::AminoAcid::LYS";
    if (three_letter == "MET") return "nmr::AminoAcid::MET";
    if (three_letter == "PHE") return "nmr::AminoAcid::PHE";
    if (three_letter == "PRO") return "nmr::AminoAcid::PRO";
    if (three_letter == "SER") return "nmr::AminoAcid::SER";
    if (three_letter == "THR") return "nmr::AminoAcid::THR";
    if (three_letter == "TRP") return "nmr::AminoAcid::TRP";
    if (three_letter == "TYR") return "nmr::AminoAcid::TYR";
    if (three_letter == "VAL") return "nmr::AminoAcid::VAL";
    return "nmr::AminoAcid::Unknown";
}


// Map a variant 3-letter code to its variant index per
// `src/AminoAcidType.h` canonical contract (asserted by
// ValidateVariantIndices()):
//   HIS:  HID=0, HIE=1, HIP=2
//   ASP:  ASH=0
//   GLU:  GLH=0
//   CYS:  CYX=0, CYM=1
//   LYS:  LYN=0
//   ARG:  ARN=0
//   TYR:  TYM=0
// Default chain (no variant) has variant_idx = 255 in this generator's
// emitter; the runtime conventionally uses 0 for "no variant" but the
// AminoAcidType variant table is 1-indexed in spirit. We use 255 as a
// distinct "default chain" marker in the emitted lookup; runtime
// callers pass the same convention. Standard 20 default chain emits
// variant_idx=255; HIS chain (which IS HIE in AMBER) also emits
// variant_idx=255 in the default block.
int VariantIndexForEmitter(const std::string& variant_3letter) {
    // 255 = "default chain (no variant)"; matches the usage in
    // SynthesisedFor<Residue> default-block emission.
    if (variant_3letter.empty()) return 255;
    if (variant_3letter == "HID") return 0;
    if (variant_3letter == "HIE") return 1;
    if (variant_3letter == "HIP") return 2;
    if (variant_3letter == "ASH") return 0;
    if (variant_3letter == "GLH") return 0;
    if (variant_3letter == "CYX") return 0;
    if (variant_3letter == "CYM") return 1;
    if (variant_3letter == "LYN") return 0;
    if (variant_3letter == "ARN") return 0;
    if (variant_3letter == "TYM") return 0;
    return 255;
}


// Write the generated .cpp file. Strings that appear in the output:
// (a) the C++ identifier for each enum value (as the typed-enum name);
// (b) atom-id comments next to each entry, for human readability.
// (a) is compile-time identifier, not runtime data; (b) is decoration
// in a // comment, also not runtime data.
void EmitCppFile(const std::string& output_path,
                 const std::vector<ResidueEntries>& all,
                 const std::vector<CapTableEntries>& caps,
                 ProcessLog& log) {
    log.Section("emit-cpp");
    log.KV("output_path", output_path);

    std::ofstream out(output_path);
    if (!out.is_open()) {
        log.KV("status", "FAILED");
        log.KV("error", "cannot open output file");
        std::fprintf(stderr, "ERROR: cannot open %s for writing\n", output_path.c_str());
        std::exit(3);
    }

    out << "// src/generated/LegacyAmberSemanticTables.cpp\n";
    out << "//\n";
    out << "// AUTOGENERATED by tools/topology/build_semantic_tables.\n";
    out << "// Do NOT edit by hand. Regenerate with the steps in\n";
    out << "// tools/topology/README.md.\n";
    out << "//\n";
    out << "// Lookup architecture: per-residue and cap tables are queried at\n";
    out << "// runtime composition time via typed structural matching against\n";
    out << "// AtomMechanicalIdentity (Element + Locant + BranchAddress +\n";
    out << "// DiastereotopicIndex + BackboneRole). Atom positions in the\n";
    out << "// arrays are NOT used as runtime keys; the LookupBy / LookupCap\n";
    out << "// functions emitted at the bottom of this file own the lookup.\n";
    out << "// See section H of\n";
    out << "// `spec/plan/bones/topology-encoding-dependencies-2026-05-05.md` for\n";
    out << "// the architectural rationale (typed identity instead of\n";
    out << "// atom_local_idx; cap atoms partitioned out of per-residue tables).\n";
    out << "//\n";
    out << "// String barrier: this file contains only typed-enum\n";
    out << "// identifiers (compile-time names) and atom-id comments.\n";
    out << "// No std::string literals, no gemmi / RDKit / cifpp\n";
    out << "// symbols, no chemistry-string round-tripping at runtime.\n";
    out << "//\n";
    out << "// Variant-index sentinel: the base (chain, no protonation\n";
    out << "// variant) table is selected by variant_idx == kBaseVariantIdx.\n";
    out << "// Any other invalid index returns nullptr (fail-fast).\n";
    out << "//\n";
    out << "// Audit trail: see src/generated/LegacyAmberSemanticTables.log.txt\n";
    out << "// for the structured generation log committed alongside.\n";
    out << "\n";
    out << "#include \"../SemanticEnums.h\"\n";
    out << "#include \"../Types.h\"\n";
    out << "\n";
    out << "namespace nmr::topology_generated {\n";
    out << "\n";
    // Emit the typed sentinel constant for the base (chain) table.
    out << "// Sentinel for the base (chain, non-variant) table. Runtime\n";
    out << "// callers pass this value when there is no protonation variant\n";
    out << "// in play; LookupBy dispatches the base table on this case and\n";
    out << "// returns nullptr for any other invalid index (fail-fast).\n";
    out << "constexpr std::uint8_t kBaseVariantIdx = 255;\n";
    out << "\n";

    int total_atoms = 0;
    for (const auto& re : all) {
        const auto array_name = CppArrayName(re.residue_3letter, re.variant_3letter);
        out << "// === " << re.residue_3letter;
        if (!re.variant_3letter.empty()) out << " (variant " << re.variant_3letter << ")";
        out << " ===\n";
        out << "constexpr std::array<AtomSemanticTable, " << re.entries.size() << "> "
            << array_name << " = {{\n";
        for (size_t i = 0; i < re.entries.size(); ++i) {
            out << EmitEntryLiteral(re.entries[i]);
            const std::string atom_id = (i < re.atom_id_for_log.size())
                                            ? re.atom_id_for_log[i] : "";
            out << ",  // " << i << ": " << atom_id << "\n";
        }
        out << "}};\n\n";
        total_atoms += static_cast<int>(re.entries.size());
        log.KV(std::string("emitted_") + array_name + "_size",
               static_cast<int>(re.entries.size()));
    }

    // Cap tables (per §H.4 of the dependencies file). One per
    // TerminalState; chain atoms live in the per-residue tables above.
    out << "// === Terminal-state cap tables ===\n";
    out << "// Per spec/plan/bones/topology-residue-reference-2026-05-05.md Section 4.\n";
    out << "// Cap atoms (OXT, HXT, H1, H2, H3) are residue-independent; one\n";
    out << "// table per terminal state covers all 20 standard residues.\n";
    out << "\n";
    int cap_total = 0;
    for (const auto& c : caps) {
        out << "constexpr std::array<AtomSemanticTable, " << c.entries.size() << "> "
            << c.array_name << " = {{\n";
        for (size_t i = 0; i < c.entries.size(); ++i) {
            out << EmitEntryLiteral(c.entries[i]);
            const std::string atom_id = (i < c.atom_id_for_log.size())
                                            ? c.atom_id_for_log[i] : "";
            out << ",  // " << i << ": " << atom_id << "\n";
        }
        out << "}};\n\n";
        cap_total += static_cast<int>(c.entries.size());
        log.KV(std::string("emitted_") + c.array_name + "_size",
               static_cast<int>(c.entries.size()));
    }

    // Lookup function emission. The runtime calls these at protein-
    // construction time once per atom; the result drives ApplySemantic
    // onto the typed semantic-record slots. Linear scan over a
    // residue's table; the typed mechanical identity is unique within
    // a residue's table by construction.
    //
    // No strings cross this barrier: enum types only.
    out << "// ============================================================================\n";
    out << "// LookupBy -- canonical runtime lookup over a residue + variant table\n";
    out << "// ============================================================================\n";
    out << "//\n";
    out << "// Returns the AtomSemanticTable entry whose mechanical identity\n";
    out << "// (Element + Locant + Branch + DiIndex + BackboneRole) matches\n";
    out << "// the queried `identity`. Linear scan; tables are small (each\n";
    out << "// residue is ~10-25 entries). Returns nullptr if no entry matches.\n";
    out << "//\n";
    out << "// `variant_idx` follows `src/AminoAcidType.h` ValidateVariantIndices()\n";
    out << "// contract: HIS HID=0, HIE=1, HIP=2; ASP ASH=0; GLU GLH=0; CYS CYX=0,\n";
    out << "// CYM=1; LYS LYN=0; ARG ARN=0; TYR TYM=0. The runtime convention\n";
    out << "// `Residue::protonation_variant_index = -1` (no variant) is mapped\n";
    out << "// to variant_idx=kBaseVariantIdx (255) by the runtime caller (cast\n";
    out << "// from -1) or passed as 255 directly; the function handles\n";
    out << "// kBaseVariantIdx as \"chain (no protonation variant)\" via each\n";
    out << "// residue's explicit case -- for HIS this returns the CCD-derived\n";
    out << "// HIS table (which corresponds to the AMBER ff14SB \"HIS = HIE\"\n";
    out << "// convention).\n";
    out << "//\n";
    out << "// FAIL-FAST: any variant_idx that is neither one of the known\n";
    out << "// variant indices for this residue nor kBaseVariantIdx returns\n";
    out << "// nullptr. Callers must pass kBaseVariantIdx for chain residues;\n";
    out << "// silently returning the base table on garbage input would mask\n";
    out << "// upstream bugs.\n";
    out << "//\n";
    out << "// THIS IS THE CANONICAL RUNTIME LOOKUP. atom_local_idx is NOT used\n";
    out << "// (index spaces don't align across CCD, AmberAminoAcidVariantTable,\n";
    out << "// and other producers; typed structural matching gives unambiguous\n";
    out << "// lookup that survives any reordering of either side). Per §H.2 of\n";
    out << "// spec/plan/bones/topology-encoding-dependencies-2026-05-05.md.\n";
    out << "namespace detail {\n";
    out << "    constexpr const AtomSemanticTable*\n";
    out << "    LookupInArray(const AtomSemanticTable* base, std::size_t n,\n";
    out << "                  const AtomMechanicalIdentity& q) {\n";
    out << "        for (std::size_t i = 0; i < n; ++i) {\n";
    out << "            const auto& e = base[i];\n";
    out << "            if (e.element == q.element\n";
    out << "                && e.locant == q.locant\n";
    out << "                && e.branch == q.branch\n";
    out << "                && e.di_index == q.di_index\n";
    out << "                && e.backbone_role == q.backbone_role) {\n";
    out << "                return &e;\n";
    out << "            }\n";
    out << "        }\n";
    out << "        return nullptr;\n";
    out << "    }\n";
    out << "}  // namespace detail\n";
    out << "\n";

    out << "const AtomSemanticTable*\n";
    out << "LookupBy(nmr::AminoAcid residue, std::uint8_t variant_idx,\n";
    out << "         const AtomMechanicalIdentity& identity) {\n";

    // Group entries by AminoAcid first; emit a switch on residue.
    // Within each residue branch, dispatch on variant_idx.
    std::map<std::string, std::vector<const ResidueEntries*>> by_residue;
    for (const auto& re : all) {
        by_residue[re.residue_3letter].push_back(&re);
    }

    out << "    switch (residue) {\n";
    for (const auto& kv : by_residue) {
        const std::string& res = kv.first;
        out << "        case " << AminoAcidEnumLiteral(res) << ": {\n";
        // Find default-chain entry (variant_3letter empty).
        const ResidueEntries* default_entry = nullptr;
        std::vector<const ResidueEntries*> variants_only;
        for (const ResidueEntries* re : kv.second) {
            if (re->variant_3letter.empty()) default_entry = re;
            else variants_only.push_back(re);
        }
        out << "            switch (variant_idx) {\n";
        // Explicit base-chain case via the typed sentinel. Replaces the
        // older "default: -> base table" pattern; default now returns
        // nullptr (fail-fast on invalid indices).
        if (default_entry != nullptr) {
            const std::string array_name =
                CppArrayName(default_entry->residue_3letter, "");
            out << "                case kBaseVariantIdx:\n";
            out << "                    return detail::LookupInArray("
                << array_name << ".data(), "
                << array_name << ".size(), identity);\n";
        }
        for (const ResidueEntries* re : variants_only) {
            const int vidx = VariantIndexForEmitter(re->variant_3letter);
            const std::string array_name =
                CppArrayName(re->residue_3letter, re->variant_3letter);
            out << "                case " << vidx << ":\n";
            out << "                    return detail::LookupInArray("
                << array_name << ".data(), "
                << array_name << ".size(), identity);\n";
        }
        out << "                default: return nullptr;\n";
        out << "            }\n";
        out << "        }\n";
    }
    out << "        default: return nullptr;\n";
    out << "    }\n";
    out << "}\n";
    out << "\n";

    // LookupCap emission.
    out << "// ============================================================================\n";
    out << "// LookupCap -- runtime lookup over a terminal-state cap table\n";
    out << "// ============================================================================\n";
    out << "//\n";
    out << "// Returns the cap-atom AtomSemanticTable entry whose mechanical\n";
    out << "// identity matches `identity`, for the requested terminal state.\n";
    out << "// Returns nullptr if no entry matches (typical when the caller\n";
    out << "// queries an atom that is not actually a cap atom on this terminus,\n";
    out << "// or when the terminal state is `Internal`).\n";
    out << "//\n";
    out << "// Override-composition rule (per §H.5 of dependencies file):\n";
    out << "// When a residue is at a terminus, `Protein::FinalizeConstruction`\n";
    out << "// MUST query both `LookupBy` (chain) and `LookupCap` (terminal-\n";
    out << "// state). If `LookupCap` returns non-null, the cap entry's fields\n";
    out << "// override the chain entry's fields for that atom. The cap table\n";
    out << "// MAY contain entries for backbone atoms (N for NTERM, C/O for\n";
    out << "// CTERM) whose chemistry differs at the terminus; those entries\n";
    out << "// have priority. Backbone overrides have non-None BackboneRole\n";
    out << "// matching the chain's same-role entry; terminus-added Hs / Os\n";
    out << "// (H1, H2, H3, OXT, HXT) carry BackboneRole::None and only match\n";
    out << "// via the cap-table lookup.\n";
    out << "//\n";
    out << "// Per §H.4 + §H.5 of spec/plan/bones/topology-encoding-dependencies-2026-05-05.md.\n";
    out << "const AtomSemanticTable*\n";
    out << "LookupCap(nmr::TerminalState state,\n";
    out << "          const AtomMechanicalIdentity& identity) {\n";
    out << "    switch (state) {\n";
    for (const auto& c : caps) {
        out << "        case " << TerminalStateLiteral(c.state) << ":\n";
        out << "            return detail::LookupInArray("
            << c.array_name << ".data(), "
            << c.array_name << ".size(), identity);\n";
    }
    out << "        case nmr::TerminalState::Internal:\n";
    out << "        default:\n";
    out << "            return nullptr;\n";
    out << "    }\n";
    out << "}\n";
    out << "\n";

    out << "}  // namespace nmr::topology_generated\n";
    out.flush();
    log.KV("total_atoms_emitted", total_atoms);
    log.KV("cap_atoms_emitted", cap_total);
    log.KV("status", "OK");
}


// Validation pass: load a residue, build mol, sanitize, extract
// facts, build typed entries, log them. Returns the entries on
// success for downstream emission.
//
// Standard-20 path: synthesis_code = residue_3letter,
// atoms_to_remove empty.
std::optional<ResidueEntries> ProcessOneResidue(cif::file& ccd,
                                                  const std::string& residue_3letter,
                                                  ProcessLog& log) {
    log.Section(std::string("process-residue::") + residue_3letter);
    auto residue_opt = LoadCcdResidue(ccd, residue_3letter, log);
    if (!residue_opt) {
        log.KV("status", "FAILED_TO_LOAD");
        return std::nullopt;
    }
    auto built = BuildRdkitMol(*residue_opt, log);
    if (built.mol == nullptr || built.mol->getNumAtoms() == 0) {
        log.KV("status", "FAILED_TO_BUILD_RDKIT_MOL");
        return std::nullopt;
    }
    if (!SanitizeMol(*built.mol, log)) {
        log.KV("status", "FAILED_TO_SANITIZE");
        return std::nullopt;
    }
    LogRdkitPerception(*built.mol, built.rdkit_idx_to_ccd_atom_id, log, residue_3letter);

    auto facts = ExtractRdkitFacts(*built.mol, log);
    // Standard-20 path: apply the chain-atom whitelist (§H.10) so the
    // per-residue table reflects AMBER's chain-residue inventory.
    auto entries = BuildResidueEntries(*residue_opt, built, facts,
                                        residue_3letter, /*atoms_to_remove*/ {},
                                        ChainAtomWhitelist(residue_3letter),
                                        log);
    LogEntriesSpotCheck(entries, log, residue_3letter);

    ResidueEntries re;
    re.residue_3letter = residue_3letter;
    re.variant_3letter = "";  // canonical form for the standard 20
    re.atom_id_for_log.reserve(entries.size());
    for (const auto& e : entries) re.atom_id_for_log.push_back(e.atom_id_for_log);
    re.entries = std::move(entries);
    log.KV("status", "OK");
    return re;
}


// Variant path: parent_3letter = the CCD entry name that owns the
// chemistry graph (HIS for HID/HIE/HIP; ASP for ASH; GLU for GLH;
// CYS for CYX/CYM; LYS for LYN; ARG for ARN; TYR for TYM). The
// variant_3letter is the AMBER protonation-variant code emitted in
// the table name (kHisAtoms_HID etc.). atoms_to_remove names the
// CCD-parent atoms that the variant lacks (e.g. HE2 for HID). The
// SynthesisedFor<Variant> patch sets remaining field deltas
// (charge_override + planar_group + polar_h + ring_position
// changes); see the per-variant SynthesisedFor<Variant> comments.
//
// Architectural note: AMBER protonation-variant codes are NOT
// canonical CCD entries -- they are AMBER-internal names. The CCD
// happens to have data_HID, data_HIE, data_HIP, data_ASH, data_GLH,
// data_CYX, data_CYM, data_LYN, data_ARN, data_TYM blocks, but
// every one of them is an unrelated small molecule that shares the
// 3-letter code by historical accident. We MUST NOT iterate variant
// codes against the CCD; we synthesise variants from the parent's
// chemistry graph.
std::optional<ResidueEntries> ProcessVariantResidue(
        cif::file& ccd,
        const std::string& variant_3letter,
        const std::string& parent_3letter,
        const std::set<std::string>& atoms_to_remove,
        ProcessLog& log) {
    log.Section(std::string("process-variant::") + variant_3letter +
                "::parent=" + parent_3letter);
    auto residue_opt = LoadCcdResidue(ccd, parent_3letter, log);
    if (!residue_opt) {
        log.KV("status", "FAILED_TO_LOAD_PARENT");
        return std::nullopt;
    }
    auto built = BuildRdkitMol(*residue_opt, log);
    if (built.mol == nullptr || built.mol->getNumAtoms() == 0) {
        log.KV("status", "FAILED_TO_BUILD_RDKIT_MOL");
        return std::nullopt;
    }
    if (!SanitizeMol(*built.mol, log)) {
        log.KV("status", "FAILED_TO_SANITIZE");
        return std::nullopt;
    }
    LogRdkitPerception(*built.mol, built.rdkit_idx_to_ccd_atom_id, log,
                       parent_3letter + "::variant=" + variant_3letter);

    auto facts = ExtractRdkitFacts(*built.mol, log);
    // Variant path: empty whitelist disables §H.10 filtering. Variant
    // tables retain their full atom set (including atoms partitioned
    // out of the chain table for that residue, e.g. HD2 in ASH, HE2 in
    // GLH, HD1 in HID/HIP).
    auto entries = BuildResidueEntries(*residue_opt, built, facts,
                                        variant_3letter, atoms_to_remove,
                                        /*chain_whitelist*/ {},
                                        log);
    LogEntriesSpotCheck(entries, log, parent_3letter + "_" + variant_3letter);

    ResidueEntries re;
    re.residue_3letter = parent_3letter;
    re.variant_3letter = variant_3letter;
    re.atom_id_for_log.reserve(entries.size());
    for (const auto& e : entries) re.atom_id_for_log.push_back(e.atom_id_for_log);
    re.entries = std::move(entries);
    log.KV("status", "OK");
    log.KV("variant_atoms_emitted", static_cast<int>(re.entries.size()));
    log.KV("variant_atoms_removed",
           static_cast<int>(residue_opt->atoms.size() - re.entries.size()));
    return re;
}


}  // namespace


int main(int argc, char** argv) {
    Args args;
    if (!ParseArgs(argc, argv, args) || args.help) {
        PrintUsage(argv[0]);
        return args.help ? 0 : 2;
    }

    ProcessLog log(args.log_path);
    log.Section("run-info");
    log.KV("generator", "build_semantic_tables");
    log.KV("ccd_path", args.ccd_path);
    log.KV("output_cpp", args.output_cpp_path);
    log.KV("log_path", args.log_path);
    log.KV("rdkit_version", RDKit::rdkitVersion);
    log.KV("step", "2 of N -- cifpp + RDKit pipeline validated on PHE");
    log.Note("Subsequent steps generalise to all 20 residues + variants,");
    log.Note("add per-field population, reconciliation per the precedence");
    log.Note("table, and emit the generated C++ source.");

    log.Section("ccd-load");
    log.KV("ccd_path", args.ccd_path);
    cif::file ccd;
    try {
        ccd.load(args.ccd_path);
    } catch (const std::exception& e) {
        log.KV("status", "FAILED");
        log.KV("error", e.what());
        std::fprintf(stderr, "ERROR loading CCD: %s\n", e.what());
        return 1;
    }
    log.KV("status", "OK");

    // Process the standard 20 amino acids first; emit one
    // ResidueEntries per residue, indexed by `kFooAtoms` in the
    // generated C++.
    std::vector<ResidueEntries> all_entries;
    for (const std::string& res : {
            std::string("ALA"),
            std::string("ARG"),
            std::string("ASN"),
            std::string("ASP"),
            std::string("CYS"),
            std::string("GLN"),
            std::string("GLU"),
            std::string("GLY"),
            std::string("HIS"),  // AMBER default; HIS in CCD has both HD1 and HE2
            std::string("ILE"),
            std::string("LEU"),
            std::string("LYS"),
            std::string("MET"),
            std::string("PHE"),
            std::string("PRO"),
            std::string("SER"),
            std::string("THR"),
            std::string("TRP"),
            std::string("TYR"),
            std::string("VAL"),
         }) {
        auto entries = ProcessOneResidue(ccd, res, log);
        if (!entries) {
            log.Section("done");
            log.KV("status", "FAILED");
            std::fprintf(stderr, "ERROR: failed to process %s; see log %s\n",
                         res.c_str(), args.log_path.c_str());
            return 1;
        }
        all_entries.push_back(std::move(*entries));
    }

    // AMBER protonation variants. Each variant is synthesised from
    // its parent CCD entry plus a per-atom delta:
    //   - atoms_to_remove: CCD-parent atoms absent from the variant
    //   - SynthesisedFor<Variant> patch: field deltas (charge,
    //     planar_group, polar_h, ring_position) for retained atoms
    //
    // Architectural note (recorded for archaeology): AMBER's variant
    // 3-letter codes (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN/ARN/TYM) are
    // NOT canonical CCD entries. The CCD does have data_X blocks for
    // those codes but they refer to unrelated small molecules that
    // share the code by historical accident (e.g. data_HID =
    // 5-hydroxy-1H-indol-3-yl acetic acid; data_LYN = lysine amide).
    // Iterating variant codes against the CCD would emit the wrong
    // chemistry. Variants are synthesised from the parent only.
    //
    // CCD parent inventory note: the CCD's HIS entry includes BOTH
    // HD1 AND HE2 (it's the imidazolium form, +1 on Nδ1); the CCD's
    // ASP includes HD2 (it's the protonated/ASH form); the CCD's GLU
    // includes HE2 (it's the protonated/GLH form). This means
    // ASH/GLH/HIP need NO atom removal -- the parent CCD already has
    // the atoms they require. HID and HIE each remove one atom (HE2
    // and HD1 respectively). Other variants remove their respective
    // atoms (HG for CYX/CYM; HZ1 for LYN; HE for ARN; HH for TYM).
    struct VariantSpec {
        std::string variant_code;
        std::string parent_code;
        std::set<std::string> atoms_to_remove;
    };
    const std::vector<VariantSpec> variants = {
        // Histidine family. CCD HIS = imidazolium (HD1 + HE2 both
        // present, +1 on Nδ1). HID drops HE2; HIE drops HD1; HIP
        // keeps both atoms but relocates +1 from Nδ1 to Nε2 (handled
        // in SynthesisedForHip).
        {"HID", "HIS", {"HE2"}},
        {"HIE", "HIS", {"HD1"}},
        {"HIP", "HIS", {}},
        // Aspartate. CCD ASP = neutral protonated (has HD2 + OXT +
        // HXT, total formal_charge 0). ASH variant just labels HD2
        // as CarboxylOH; no atom delta.
        {"ASH", "ASP", {}},
        // Glutamate. CCD GLU = neutral protonated (has HE2 + OXT +
        // HXT). GLH labels HE2 as CarboxylOH; no atom delta.
        {"GLH", "GLU", {}},
        // Cysteine family. CCD CYS = thiol (has HG). CYX (disulfide-
        // bonded) and CYM (thiolate) both drop HG; CYM also sets
        // SG charge to -1 (handled in SynthesisedForCym).
        {"CYX", "CYS", {"HG"}},
        {"CYM", "CYS", {"HG"}},
        // Lysine. CCD LYS = ammonium (HZ1, HZ2, HZ3, +1 on Nζ).
        // LYN drops HZ1 and clears the +1 on Nζ; HZ2/HZ3 become
        // AmineNH (handled in SynthesisedForLyn).
        {"LYN", "LYS", {"HZ1"}},
        // Arginine. CCD ARG = guanidinium (HE + HH11 + HH12 + HH21 +
        // HH22, +1 on NH2). ARN drops HE and clears the +1 on NH2
        // (handled in SynthesisedForArn). Note: AMBER ff14SB has no
        // canonical ARN template; HE is removed by chemistry
        // inference (lone pair on Nε after deprotonation). See
        // dependencies §B.1, §F.
        {"ARN", "ARG", {"HE"}},
        // Tyrosine. CCD TYR = neutral phenol (has HH on OH). TYM
        // drops HH; OH becomes AromaticOxide and gets formal_chg=-1
        // (handled in SynthesisedForTym).
        {"TYM", "TYR", {"HH"}},
    };
    for (const auto& v : variants) {
        auto entries = ProcessVariantResidue(ccd, v.variant_code,
                                              v.parent_code,
                                              v.atoms_to_remove, log);
        if (!entries) {
            log.Section("done");
            log.KV("status", "FAILED_VARIANT");
            log.KV("failed_variant", v.variant_code);
            std::fprintf(stderr, "ERROR: failed to process variant %s "
                                  "(parent %s); see log %s\n",
                         v.variant_code.c_str(), v.parent_code.c_str(),
                         args.log_path.c_str());
            return 1;
        }
        all_entries.push_back(std::move(*entries));
    }

    // Build the four terminal-state cap tables. Per §H.4 of
    // spec/plan/bones/topology-encoding-dependencies-2026-05-05.md, cap
    // chemistry is residue-independent; one table per terminal state
    // covers all 20 standard residues.
    auto cap_tables = BuildCapTables(log);

    EmitCppFile(args.output_cpp_path, all_entries, cap_tables, log);

    log.Section("done");
    log.KV("status", "OK");
    log.KV("residues_emitted", static_cast<int>(all_entries.size()));
    log.Note("Standard 20 + 10 AMBER protonation variants emitted.");
    log.Note("Standard 20 byte-identical regeneration is the contract:");
    log.Note("variant emission appends after the 20 and never modifies");
    log.Note("the 20's tables.");

    std::printf("OK -- generated %s with %d residue table(s).\n"
                "Log: %s\n",
                args.output_cpp_path.c_str(),
                static_cast<int>(all_entries.size()),
                args.log_path.c_str());
    return 0;
}

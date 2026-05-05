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

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <optional>
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
#include <RDGeneral/versions.h>


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


// Scaffolded validation pass: load PHE from CCD, build an RDKit mol,
// sanitize, log perception. This is Step 2-3 of the generator
// pipeline; subsequent steps add per-field population, reconciliation,
// and emission.
int ValidatePheRoundTrip(cif::file& ccd, ProcessLog& log) {
    log.Section("phe-roundtrip-validation");
    auto residue_opt = LoadCcdResidue(ccd, "PHE", log);
    if (!residue_opt) {
        log.KV("status", "FAILED_TO_LOAD_PHE");
        return 1;
    }
    auto built = BuildRdkitMol(*residue_opt, log);
    if (built.mol == nullptr || built.mol->getNumAtoms() == 0) {
        log.KV("status", "FAILED_TO_BUILD_RDKIT_MOL");
        return 1;
    }
    if (!SanitizeMol(*built.mol, log)) {
        log.KV("status", "FAILED_TO_SANITIZE");
        return 1;
    }
    LogRdkitPerception(*built.mol, built.rdkit_idx_to_ccd_atom_id, log, "PHE");
    log.KV("status", "OK");
    return 0;
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

    if (int rc = ValidatePheRoundTrip(ccd, log); rc != 0) {
        log.Section("done");
        log.KV("status", "FAILED");
        std::fprintf(stderr, "ERROR: PHE round-trip validation failed; see log %s\n",
                     args.log_path.c_str());
        return rc;
    }

    log.Section("done");
    log.KV("status", "OK");
    log.Note("Step 2 of N done: cifpp -> RDKit pipeline validated on PHE.");
    log.Note("No tables emitted yet.");

    std::printf("OK -- step 2 (PHE round-trip via cifpp + RDKit) complete.\n"
                "See log at %s\n",
                args.log_path.c_str());
    return 0;
}

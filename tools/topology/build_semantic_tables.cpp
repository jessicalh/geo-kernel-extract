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

using namespace nmr;


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
// Atom-name parser -- mechanical fields (Locant, BranchAddress,
// DiastereotopicIndex) from the CCD atom_id string plus the residue's
// bond graph (to map H atoms to their heavy-atom parent).
//
// String boundary: this parser consumes CCD atom-name strings and
// returns typed-enum values. The strings die at this function's
// return; callers see typed values only.
// ============================================================================

struct ParsedName {
    Element             element;
    bool                is_backbone     = false;  ///< N/CA/C/O/H/HA (and HN alias)
    bool                is_n_terminus   = false;  ///< H1/H2/H3 (extra ammonium Hs)
    bool                is_c_terminus   = false;  ///< OXT/HXT
    Locant              locant          = Locant::None;
    BranchAddress       branch          = {};
    DiastereotopicIndex di_index        = DiastereotopicIndex::None;
};

Element ElementFromAtomName(const std::string& name) {
    if (name.empty()) return Element::Unknown;
    switch (name[0]) {
        case 'H': return Element::H;
        case 'C': return Element::C;
        case 'N': return Element::N;
        case 'O': return Element::O;
        case 'S': return Element::S;
        default:  return Element::Unknown;
    }
}

Locant LocantLetterToEnum(char c) {
    switch (c) {
        case 'A': return Locant::Alpha;
        case 'B': return Locant::Beta;
        case 'G': return Locant::Gamma;
        case 'D': return Locant::Delta;
        case 'E': return Locant::Epsilon;
        case 'Z': return Locant::Zeta;
        case 'H': return Locant::Eta;     // Note: 'H' as locant suffix (Arg/Tyr eta), NOT element.
        default:  return Locant::None;
    }
}

// Parse an atom name plus its parent map (H -> heavy parent name).
// For non-H atoms parent_name is empty.
ParsedName ParseAtomName(const std::string& name, const std::string& parent_name) {
    ParsedName p;
    p.element = ElementFromAtomName(name);

    // Pure backbone names.
    if (name == "N" || name == "CA" || name == "C" || name == "O" ||
        name == "H" || name == "HN" || name == "HA") {
        p.is_backbone = true;
        if (name == "HA") p.locant = Locant::Alpha;  // HA is a Cα-H
        return p;
    }
    // N-terminal extra ammonium Hs.
    if (name == "H1" || name == "H2" || name == "H3" || name == "H2N") {
        p.is_n_terminus = true;
        return p;
    }
    // C-terminal carboxyl atoms.
    if (name == "OXT" || name == "HXT" || name == "OT1" || name == "OT2") {
        p.is_c_terminus = true;
        return p;
    }
    // Glycine alpha-Hs (only diastereotopic CH2 in the standard 20).
    if (name == "HA2" || name == "HA3") {
        p.is_backbone = true;
        p.locant = Locant::Alpha;
        p.di_index = (name.back() == '2') ? DiastereotopicIndex::Position2
                                          : DiastereotopicIndex::Position3;
        return p;
    }

    // Sidechain: <element-letter><locant-letter><digits...>.
    // The locant letter is the second character of the atom name.
    if (name.size() < 2) return p;  // unparseable, give up
    p.locant = LocantLetterToEnum(name[1]);
    if (p.locant == Locant::None) return p;

    // Trailing digits form a numeric suffix: "" / "1" / "2" / "12" / "21" etc.
    const std::string suffix = name.substr(2);
    if (suffix.empty()) return p;  // single sidechain atom (e.g. CB, SG, OH)

    // For non-H atoms with a digit suffix: the digit is the BranchAddress::outer
    // (Val CG1/CG2, Leu CD1/CD2, Phe CD1/CD2, Asn OD1/ND2, Asn OD1, etc.).
    if (p.element != Element::H) {
        // Single digit: simple branch.
        if (suffix.size() == 1 && std::isdigit(static_cast<unsigned char>(suffix[0]))) {
            p.branch.outer = static_cast<uint8_t>(suffix[0] - '0');
            return p;
        }
        // Multi-digit on heavy atom is unusual; record as outer for now.
        p.branch.outer = static_cast<uint8_t>(suffix[0] - '0');
        return p;
    }

    // For H atoms: use the heavy-atom parent map to disambiguate.
    //   - parent has no branch digit (e.g. CB, CG): suffix encodes
    //     diastereotopic pair (HB2/HB3).
    //   - parent has branch digit (e.g. CD1, NH1): suffix's first
    //     digit is the parent's branch (outer); remaining digits
    //     are inner (e.g. HD11 = on CD1, inner=1).
    bool parent_has_digit = false;
    char parent_digit = 0;
    if (!parent_name.empty()) {
        for (char c : parent_name) {
            if (std::isdigit(static_cast<unsigned char>(c))) {
                parent_has_digit = true;
                parent_digit = c;
                break;
            }
        }
    }

    if (!parent_has_digit) {
        // Diastereotopic methylene-style: HB2 vs HB3, HG2 vs HG3, etc.
        // Methyl-on-unbranched (Ala HB1/HB2/HB3) also lands here; we
        // record the digit as DiastereotopicIndex but the methyl-pair
        // distinction is downstream (PseudoatomMembership).
        if (suffix.size() == 1 && std::isdigit(static_cast<unsigned char>(suffix[0]))) {
            const char d = suffix[0];
            if (d == '2') p.di_index = DiastereotopicIndex::Position2;
            else if (d == '3') p.di_index = DiastereotopicIndex::Position3;
            // d == '1' (Ala HB1) leaves di_index = None; the membership
            // of a methyl pseudoatom (MB) is what carries the meaning.
        }
        return p;
    }

    // Parent has branch digit: outer = parent's branch, inner = remaining
    // digit on the H name (if any).
    p.branch.outer = static_cast<uint8_t>(parent_digit - '0');
    // Extract the H atom's inner digit: search the suffix for a digit
    // that isn't the parent's outer match.
    for (char c : suffix) {
        if (!std::isdigit(static_cast<unsigned char>(c))) continue;
        if (c == parent_digit && p.branch.inner == 0 && suffix.size() == 1) {
            // Single matching digit (e.g. HG1 on CG1 means just "the H on CG1").
            return p;
        }
        if (c != parent_digit) {
            p.branch.inner = static_cast<uint8_t>(c - '0');
            break;
        }
    }
    // Special case: HG21 on Thr CG2 -> outer=2, inner=1.
    // The above loop catches that because digits after the first match are inner.
    if (p.branch.inner == 0) {
        // Suffix has multiple matching digits or only one: handle by
        // taking the last digit as inner.
        p.branch.inner = static_cast<uint8_t>(suffix.back() - '0');
        if (p.branch.inner == p.branch.outer && suffix.size() == 1) {
            // Suffix was just the parent digit; H is the only one on
            // that branch (e.g. Thr HG1 on OG1).
            p.branch.inner = 0;
        }
    }
    return p;
}


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
};


// PHE-only initial implementation; generalises in a follow-up commit
// once the architecture is proven on one residue.
SynthesisedFields SynthesisedForPhe(const std::string& atom_id, const ParsedName& parsed) {
    SynthesisedFields s;

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
    if      (atom_id == "CG"  || atom_id == "HG")  set_ring(RingPositionLabel::Ipso);
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


// ============================================================================
// AtomSemanticEntry -- the typed record for one atom in the residue.
// Combines mechanical (parser), algorithmic (RDKit), and synthesised
// (per-residue table) sources into a single typed record with
// provenance witnesses.
// ============================================================================

struct AtomSemanticEntry {
    // Identification (in the generator scope only -- the runtime
    // table indexes by (residue, atom_local_idx), not by name).
    std::string               atom_id_for_log;

    // Typed fields, all 14 of them.
    Locant                    locant            = Locant::None;
    BranchAddress             branch            = {};
    DiastereotopicIndex       di_index          = DiastereotopicIndex::None;
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
                                          int formal_charge_from_ccd) {
    AtomSemanticEntry e;
    e.atom_id_for_log = atom_id;

    // Mechanical fields from atom-name parser.
    e.locant   = parsed.locant;
    e.branch   = parsed.branch;
    e.di_index = parsed.di_index;
    e.prov_locant.witnesses[0] = {SemanticSource::IUPAC_1969,
                                  static_cast<uint8_t>(parsed.locant)};
    e.prov_locant.confidence   = SemanticConfidence::AlgorithmAuthoritative;
    e.prov_locant.sources_agree = true;

    // Algorithmic fields from RDKit.
    e.prochiral     = facts.cip_stereo[rdkit_idx];
    e.aromatic      = facts.aromatic[rdkit_idx];
    e.hybridisation = facts.hybridisation[rdkit_idx];
    e.bond_orders   = facts.bond_orders[rdkit_idx];
    e.equivalence_class = facts.canonical_rank[rdkit_idx];
    e.formal_charge = static_cast<int8_t>(formal_charge_from_ccd);

    e.prov_prochiral.witnesses[0] = {SemanticSource::RDKit_CIPLabeler,
                                     static_cast<uint8_t>(facts.cip_stereo[rdkit_idx])};
    e.prov_prochiral.confidence   = (facts.cip_stereo[rdkit_idx] == ProchiralStereo::NotProchiral)
        ? SemanticConfidence::AlgorithmAuthoritative
        : SemanticConfidence::CIPVerified;
    e.prov_prochiral.cip_verified = (facts.cip_stereo[rdkit_idx] != ProchiralStereo::NotProchiral &&
                                      facts.cip_stereo[rdkit_idx] != ProchiralStereo::Unassigned);
    e.prov_prochiral.sources_agree = true;

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


// Build the per-atom semantic entries for one residue.
std::vector<AtomSemanticEntry> BuildResidueEntries(const CcdResidue& res,
                                                    const RdkitMolWithMap& built,
                                                    const RdkitFacts& facts,
                                                    ProcessLog& log) {
    log.Section(std::string("build-entries::") + res.three_letter);
    std::vector<AtomSemanticEntry> entries;
    entries.reserve(res.atoms.size());

    auto parents = BuildHydrogenParentMap(res);
    std::map<std::string, unsigned int> name_to_rdkit_idx;
    for (unsigned int i = 0; i < built.rdkit_idx_to_ccd_atom_id.size(); ++i) {
        name_to_rdkit_idx[built.rdkit_idx_to_ccd_atom_id[i]] = i;
    }

    for (const auto& ccd_atom : res.atoms) {
        std::string parent_name;
        if (ccd_atom.type_symbol == "H") {
            auto it = parents.find(ccd_atom.atom_id);
            if (it != parents.end()) parent_name = it->second;
        }
        const ParsedName parsed = ParseAtomName(ccd_atom.atom_id, parent_name);

        unsigned int rdkit_idx = 0;
        auto it = name_to_rdkit_idx.find(ccd_atom.atom_id);
        if (it != name_to_rdkit_idx.end()) rdkit_idx = it->second;

        SynthesisedFields syn;
        if (res.three_letter == "PHE") {
            syn = SynthesisedForPhe(ccd_atom.atom_id, parsed);
        }

        AtomSemanticEntry e = BuildAtomSemanticEntry(ccd_atom.atom_id, parsed,
                                                      rdkit_idx, facts, syn,
                                                      ccd_atom.formal_charge);
        entries.push_back(std::move(e));
    }
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
        fact += "loc=" + std::to_string(static_cast<int>(e.locant));
        fact += " branch=" + std::to_string(e.branch.outer) + "/" + std::to_string(e.branch.inner);
        fact += " di=" + std::to_string(static_cast<int>(e.di_index));
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
// AtomSemanticTable per residue, plus a Lookup function indexed by
// (AminoAcid, atom_local_idx). The output includes only typed-enum
// references; no std::string literals carrying chemistry data, no
// gemmi/RDKit/cifpp symbols, no Eigen.
//
// The string barrier from the generator's perspective: this Emit
// function is the only place where chemistry-derived enum names are
// converted back to text -- but only as C++ source code identifiers
// (e.g. "nmr::Locant::Beta"), not as data. The compiler then turns
// those identifier strings back into the typed-enum integer values
// at compile time. No runtime string survives.
// ============================================================================

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
    }
    return "nmr::RingSystemKind::NotInRing";
}

const char* RingPositionLabelLiteral(RingPositionLabel r) {
    switch (r) {
        case RingPositionLabel::NotInRing:      return "nmr::RingPositionLabel::NotInRing";
        case RingPositionLabel::Ipso:           return "nmr::RingPositionLabel::Ipso";
        case RingPositionLabel::Ortho1:         return "nmr::RingPositionLabel::Ortho1";
        case RingPositionLabel::Ortho2:         return "nmr::RingPositionLabel::Ortho2";
        case RingPositionLabel::Meta1:          return "nmr::RingPositionLabel::Meta1";
        case RingPositionLabel::Meta2:          return "nmr::RingPositionLabel::Meta2";
        case RingPositionLabel::Para:           return "nmr::RingPositionLabel::Para";
        case RingPositionLabel::PyrroleAlpha:   return "nmr::RingPositionLabel::PyrroleAlpha";
        case RingPositionLabel::PyrroleBeta:    return "nmr::RingPositionLabel::PyrroleBeta";
        case RingPositionLabel::BridgeFusion:   return "nmr::RingPositionLabel::BridgeFusion";
        case RingPositionLabel::Heteroatom_NH:  return "nmr::RingPositionLabel::Heteroatom_NH";
        case RingPositionLabel::Heteroatom_NoH: return "nmr::RingPositionLabel::Heteroatom_NoH";
        case RingPositionLabel::Heteroatom_OH:  return "nmr::RingPositionLabel::Heteroatom_OH";
        case RingPositionLabel::Saturated:      return "nmr::RingPositionLabel::Saturated";
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
std::string EmitEntryLiteral(const AtomSemanticEntry& e) {
    std::ostringstream o;
    o << "    { " << LocantLiteral(e.locant)
      << ", {" << static_cast<int>(e.branch.outer)
      <<  "," << static_cast<int>(e.branch.inner) << "}"
      << ", " << DiastereotopicLiteral(e.di_index)
      << ", " << ProchiralLiteral(e.prochiral)
      << ", " << PlanarGroupLiteral(e.planar_group)
      << ", " << PlanarStereoLiteral(e.planar_stereo)
      << ", " << PseudoatomLiteral(e.pseudoatom)
      << ", " << PolarHLiteral(e.polar_h)
      << ", {" << RingMembershipLiteral(e.ring_position.primary)
      <<  ", " << RingMembershipLiteral(e.ring_position.secondary) << "}"
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


// Write the generated .cpp file. Strings that appear in the output:
// (a) the C++ identifier for each enum value (as the typed-enum name);
// (b) atom-id comments next to each entry, for human readability.
// (a) is compile-time identifier, not runtime data; (b) is decoration
// in a // comment, also not runtime data.
void EmitCppFile(const std::string& output_path,
                 const std::vector<ResidueEntries>& all,
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
    out << "// String barrier: this file contains only typed-enum\n";
    out << "// identifiers (compile-time names) and atom-id comments.\n";
    out << "// No std::string literals, no gemmi / RDKit / cifpp\n";
    out << "// symbols, no chemistry-string round-tripping at runtime.\n";
    out << "//\n";
    out << "// Audit trail: see src/generated/LegacyAmberSemanticTables.log.txt\n";
    out << "// for the structured generation log committed alongside.\n";
    out << "\n";
    out << "#include \"../SemanticEnums.h\"\n";
    out << "\n";
    out << "namespace nmr::topology_generated {\n";
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

    out << "}  // namespace nmr::topology_generated\n";
    out.flush();
    log.KV("total_atoms_emitted", total_atoms);
    log.KV("status", "OK");
}


// Validation pass: load a residue, build mol, sanitize, extract
// facts, build typed entries, log them. Returns the entries on
// success for downstream emission.
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
    auto entries = BuildResidueEntries(*residue_opt, built, facts, log);
    LogEntriesSpotCheck(entries, log, residue_3letter);

    ResidueEntries re;
    re.residue_3letter = residue_3letter;
    re.variant_3letter = "";  // canonical form for now
    re.atom_id_for_log.reserve(entries.size());
    for (const auto& e : entries) re.atom_id_for_log.push_back(e.atom_id_for_log);
    re.entries = std::move(entries);
    log.KV("status", "OK");
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

    // Process residues one by one. Currently only PHE is in the
    // synthesised-fields lookup (SynthesisedForPhe); other residues
    // get NotInRing / NotPolar defaults until per-residue tables
    // are added in the generalisation step. The architecture is
    // proven on PHE end-to-end; subsequent residues land mechanically
    // by extending the synthesised-fields lookup.
    std::vector<ResidueEntries> all_entries;
    for (const std::string& res : {std::string("PHE")}) {
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

    EmitCppFile(args.output_cpp_path, all_entries, log);

    log.Section("done");
    log.KV("status", "OK");
    log.KV("residues_emitted", static_cast<int>(all_entries.size()));
    log.Note("Step 5 of N done: emitter writes generated C++ source.");
    log.Note("Architecture proven end-to-end on PHE.");
    log.Note("Generalisation to remaining 19 residues + 10 variants is the");
    log.Note("next session's work; the synthesised-field tables in");
    log.Note("SynthesisedForPhe extend mechanically per residue.");

    std::printf("OK -- generated %s with %d residue table(s).\n"
                "Log: %s\n",
                args.output_cpp_path.c_str(),
                static_cast<int>(all_entries.size()),
                args.log_path.c_str());
    return 0;
}

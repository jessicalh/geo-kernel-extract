#include "PdbFileReader.h"
#include "ReduceProtonation.h"
#include "ChargeSource.h"
#include "AmberChargeResolver.h"
#include "ForceFieldChargeTable.h"
#include "RuntimeEnvironment.h"
#include "AminoAcidType.h"
#include "NamingRegistry.h"
#include "OperationLog.h"
// RuntimeEnvironment::Ff14sbParams() resolves the param file location.
// Discovery: NMR_FF14SB_PARAMS env → data/ff14sb_params.dat → TOML config.
#include <cif++.hpp>
#include <cif++/pdb/pdb2cif.hpp>
#include <dssp.hpp>
#include <fstream>
#include <sstream>
#include <map>
#include <filesystem>
#include <cstdio>
#include <tuple>
#include <random>

namespace fs = std::filesystem;

namespace nmr {


// ============================================================================
// Internal: parse a PDB string into a Protein.
// This is the cif++ parsing logic. Used by BuildFromPdb on the
// protonated PDB output from reduce.
// ============================================================================

static std::unique_ptr<Protein> ParsePdb(const std::string& pdb_text,
                                          std::string& error) {
    // Fix PDB text for cif++ compatibility
    std::istringstream rawStream(pdb_text);
    std::stringstream fixedPdb;
    std::string rawLine;
    bool hasEnd = false;
    while (std::getline(rawStream, rawLine)) {
        if (rawLine.empty()) continue;
        std::string rec = rawLine.substr(0, std::min(rawLine.size(), size_t(6)));

        // Skip records that confuse cif++ strict parser
        if (rec == "TITLE " || rec == "REMARK" || rec == "MODEL " ||
            rec == "ENDMDL") continue;
        if (rec == "USER  ") continue;  // reduce adds USER MOD lines
        if (rec == "END   " || rawLine == "END") hasEnd = true;

        // Fix blank chain ID (column 21, 0-indexed)
        if ((rec == "ATOM  " || rec == "HETATM") &&
            rawLine.size() > 21 && rawLine[21] == ' ') {
            rawLine[21] = 'A';
        }
        fixedPdb << rawLine << '\n';
    }
    if (!hasEnd) fixedPdb << "END\n";

    // Parse PDB via cif++
    std::istringstream pdbStream(fixedPdb.str());
    cif::file cifFile;
    try {
        cif::pdb::ReadPDBFile(pdbStream, cifFile);
    } catch (const std::exception& e) {
        error = std::string("cif++ parse error: ") + e.what();
        return nullptr;
    }

    if (cifFile.empty()) {
        error = "Empty CIF file from PDB input";
        return nullptr;
    }

    // Build Protein atoms and residues
    cif::datablock& db = cifFile.front();
    cif::mm::structure structure(db);

    auto protein = std::make_unique<Protein>();
    std::vector<Vec3> positions;

    struct ResKey {
        std::string chain;
        std::string seq_id;
        bool operator<(const ResKey& o) const {
            return std::tie(chain, seq_id) < std::tie(o.chain, o.seq_id);
        }
    };
    std::map<ResKey, size_t> res_map;

    for (auto& poly : structure.polymers()) {
        std::string chain_id = poly.get_auth_asym_id();

        for (auto& mono : poly) {
            std::string comp_id = mono.get_compound_id();
            std::string seq_id = mono.get_auth_seq_id();

            AminoAcid aa_type = AminoAcidFromThreeLetterCode(comp_id);
            if (aa_type == AminoAcid::Unknown) continue;

            ResKey key{chain_id, seq_id};
            size_t res_idx;

            auto it = res_map.find(key);
            if (it == res_map.end()) {
                Residue res;
                res.type = aa_type;
                res.chain_id = chain_id;
                try {
                    res.sequence_number = std::stoi(seq_id);
                } catch (...) {
                    res.sequence_number = 0;
                }
                res_idx = protein->AddResidue(res);
                res_map[key] = res_idx;
            } else {
                res_idx = it->second;
            }

            // PDB LOADING BOUNDARY (the cifpp surface).
            //
            // Two-pass per-residue canonicalisation via the
            // NamingApplicator. The first pass collects every atom's
            // raw input name to form the sibling-set snapshot; the
            // second pass calls applicator.Apply(ctx) per atom with
            // the snapshot in NamingContext::sibling_input_names.
            // Sibling-aware predicates examine the snapshot to decide
            // whether to fire (e.g. PRO HD1->HD2 fires only when
            // siblings show {HD1,HD2,no HD3}; on canonical {HD2,HD3}
            // input the predicate evaluates false and the rule does
            // not fire — the input passes through unchanged).
            //
            // source = NamingSource::CifppPdbInput. cifpp-loaded PDB
            // files carry IUPAC / canonical-AMBER-ff14SB names; the
            // applicator's pass-through branch in Resolve() returns
            // them unchanged. Pre-Markley fixtures (LYN HZ1/HZ2 only,
            // GLY collapsed-HA) have rules tagged AmberFf14SBCanonical
            // that fire independent of source and produce canonical
            // outputs.
            //
            // KNOWN GAP (deferred to PROPKA wiring at the
            // PdbFileReader TODO, separate commit): residues labelled
            // "LYS" but structurally LYN (only HZ1/HZ2 NZ-H, no HZ3)
            // canonicalise correctly during this load-time pass via
            // the LYN HZ1/HZ2 sibling-aware shift rules (siblings
            // {HZ1,HZ2,no HZ3} matches the LYN-pre-Markley pattern).
            // The post-protonation second applicator pass at
            // FinalizeConstruction is idempotent on the now-canonical
            // names.
            const auto& applicator = GlobalNamingApplicator();

            // First pass: harvest raw atom names for the sibling
            // snapshot, plus the per-atom (element, position, alt-A
            // gating) data we'll need for the second pass.
            struct PendingAtom {
                std::string raw_name;
                Element     elem;
                Vec3        pos;
            };
            std::vector<PendingAtom> pending;
            pending.reserve(mono.atoms().size());
            for (auto& atom : mono.atoms()) {
                std::string alt_id = atom.get_label_alt_id();
                if (!alt_id.empty() && alt_id != "A") continue;

                Element elem = ElementFromSymbol(
                    cif::atom_type_traits(atom.get_type()).symbol());
                if (elem == Element::Unknown) continue;

                auto [x, y, z] = atom.get_location();
                pending.push_back({atom.get_label_atom_id(), elem, Vec3(x, y, z)});
            }

            // Compute the sibling snapshot once for the residue.
            std::vector<std::string> input_names;
            input_names.reserve(pending.size());
            for (const PendingAtom& p : pending) {
                input_names.push_back(p.raw_name);
            }
            // Parent-input-name vector is empty (PdbFileReader doesn't
            // resolve heavy-atom parents at load time; this happens in
            // CovalentTopology::Resolve later).
            const std::vector<std::string> parent_names;

            const Residue& res_now = protein->ResidueAt(res_idx);
            const auto canonical_names = applicator.ApplyResidue(
                input_names,
                parent_names,
                res_now.type,
                /*variant_index=*/-1,
                TerminalState::Internal,
                NamingSource::CifppPdbInput,
                res_now.sequence_number,
                res_now.chain_id);

            // Second pass: write atoms with canonical names.
            for (size_t k = 0; k < pending.size(); ++k) {
                auto new_atom = Atom::Create(pending[k].elem);
                new_atom->pdb_atom_name = canonical_names[k];
                new_atom->residue_index = res_idx;
                size_t atom_idx = protein->AddAtom(std::move(new_atom));
                protein->MutableResidueAt(res_idx).atom_indices.push_back(atom_idx);
                positions.push_back(pending[k].pos);
            }
        }
    }

    if (positions.empty()) {
        error = "No protein atoms found in PDB input";
        return nullptr;
    }

    protein->FinalizeConstruction(positions);
    protein->AddConformation(std::move(positions), "protonated");

    return protein;
}


// ============================================================================
// Generate a unique temp file path (GUID-like to avoid collisions)
// ============================================================================

static std::string UniqueTempPath(const std::string& stem) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dist;
    char buf[256];
    std::snprintf(buf, sizeof(buf), "/tmp/nmr_%s_%016lx.pdb",
                  stem.c_str(), dist(gen));
    return buf;
}


// ============================================================================
// BuildFromPdb: the public entry point.
//
// Read PDB → protonate with reduce → parse → assign charges → net charge.
// ============================================================================

BuildResult BuildFromPdb(const std::string& path, double pH) {
    BuildResult result;

    OperationLog::Scope scope("BuildFromPdb",
        path + " pH=" + std::to_string(pH));

    // 1. Read original PDB content
    std::ifstream f(path);
    if (!f.is_open()) {
        result.error = "Cannot open '" + path + "'";
        return result;
    }
    std::string pdb_content((std::istreambuf_iterator<char>(f)),
                             std::istreambuf_iterator<char>());
    f.close();

    // 2. Protonate with reduce
    std::string protonated = ProtonateWithReduce(pdb_content);
    if (protonated.empty()) {
        result.error = "reduce protonation failed for " + path;
        return result;
    }

    // 3. Parse protonated structure
    std::string parse_error;
    result.protein = ParsePdb(protonated, parse_error);
    if (!result.protein) {
        result.error = "failed to parse protonated PDB: " + parse_error;
        return result;
    }

    // 4. Build context
    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->pdb_source = fs::path(path).filename().string();
    ctx->protonation_tool = "reduce";
    ctx->protonation_pH = pH;
    ctx->force_field = "ff14SB";
    result.protein->SetBuildContext(std::move(ctx));

    // 5. Resolve charge source through the typed dispatch.
    AmberSourceConfig source_config;
    source_config.flat_table_path = RuntimeEnvironment::Ff14sbParams();
    source_config.preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;
    if (source_config.flat_table_path.empty() ||
            !fs::exists(source_config.flat_table_path)) {
        result.error = "ff14SB params not found (RuntimeEnvironment::Ff14sbParams()="
                       + (source_config.flat_table_path.empty()
                              ? std::string("<not set>")
                              : source_config.flat_table_path) + ")";
        return result;
    }

    std::string charge_err;
    result.charges = ResolveAmberChargeSource(
        *result.protein, result.protein->BuildContext(),
        source_config, charge_err);
    if (!result.charges) {
        result.error = charge_err;
        return result;
    }

    // Compute net charge from charge sum
    auto& conf = result.protein->Conformation();
    if (!result.protein->PrepareForceFieldCharges(*result.charges, conf, charge_err)) {
        result.error = "charge preparation failed: " + charge_err;
        return result;
    }
    double charge_sum = result.protein->ForceFieldCharges().TotalCharge();
    result.net_charge = static_cast<int>(
        charge_sum + (charge_sum > 0 ? 0.5 : -0.5));

    result.ok = true;

    OperationLog::Info(LogCharges, "BuildFromPdb",
        "protonated " + fs::path(path).filename().string() + ": " +
        std::to_string(result.protein->AtomCount()) + " atoms, " +
        "net_charge=" + std::to_string(result.net_charge));

    return result;
}


// ============================================================================
// BuildFromProtonatedPdb: load an already-protonated PDB, assign charges.
// No reduce call. For test data and pre-processed inputs.
// ============================================================================

BuildResult BuildFromProtonatedPdb(const std::string& path) {
    BuildResult result;

    OperationLog::Scope scope("BuildFromProtonatedPdb", path);

    // Read and parse directly — no reduce step
    std::ifstream f(path);
    if (!f.is_open()) {
        result.error = "Cannot open '" + path + "'";
        return result;
    }
    std::string pdb_content((std::istreambuf_iterator<char>(f)),
                             std::istreambuf_iterator<char>());
    f.close();

    std::string parse_error;
    result.protein = ParsePdb(pdb_content, parse_error);
    if (!result.protein) {
        result.error = "failed to parse PDB: " + parse_error;
        return result;
    }

    // Build context
    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->pdb_source = fs::path(path).filename().string();
    ctx->force_field = "ff14SB";
    result.protein->SetBuildContext(std::move(ctx));

    // Resolve charge source through the typed dispatch.
    AmberSourceConfig source_config;
    source_config.flat_table_path = RuntimeEnvironment::Ff14sbParams();
    source_config.preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;
    if (source_config.flat_table_path.empty() ||
            !fs::exists(source_config.flat_table_path)) {
        result.error = "ff14SB params not found (RuntimeEnvironment::Ff14sbParams()="
                       + (source_config.flat_table_path.empty()
                              ? std::string("<not set>")
                              : source_config.flat_table_path) + ")";
        return result;
    }

    std::string charge_err;
    result.charges = ResolveAmberChargeSource(
        *result.protein, result.protein->BuildContext(),
        source_config, charge_err);
    if (!result.charges) {
        result.error = charge_err;
        return result;
    }

    auto& conf = result.protein->Conformation();
    if (!result.protein->PrepareForceFieldCharges(*result.charges, conf, charge_err)) {
        result.error = "charge preparation failed: " + charge_err;
        return result;
    }
    double charge_sum = result.protein->ForceFieldCharges().TotalCharge();
    result.net_charge = static_cast<int>(
        charge_sum + (charge_sum > 0 ? 0.5 : -0.5));

    result.ok = true;
    return result;
}


}  // namespace nmr

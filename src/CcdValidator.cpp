#include "CcdValidator.h"
#include "Protein.h"
#include "Atom.h"
#include "Residue.h"
#include "AminoAcidType.h"
#include "RuntimeEnvironment.h"
#include "OperationLog.h"
#include "NpyWriter.h"

#include <gemmi/cif.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/chemcomp.hpp>
#include <gemmi/elem.hpp>

#include <filesystem>
#include <fstream>
#include <mutex>
#include <set>
#include <unordered_map>
#include <utility>

namespace fs = std::filesystem;

namespace nmr {

// ============================================================================
// CcdStore: process-global, lazy-loaded chemical component dictionary.
//
// Reading components.cif.gz takes a few seconds and ~half a gigabyte of RAM.
// The Document is parsed once and held for the program's lifetime. Lookups
// are by three-letter code into a name-keyed block index; ChemComp objects
// are constructed on first lookup and cached.
// ============================================================================

namespace {

class CcdStore {
public:
    static CcdStore& Instance() {
        static CcdStore inst;
        return inst;
    }

    // Returns nullptr if the residue code is not present in CCD.
    // CCD load itself is fatal-on-failure (see Load), so reaching this
    // point implies the dictionary is loaded.
    const gemmi::ChemComp* Find(const std::string& residue_three_letter) {
        std::call_once(load_once_, [this]() { Load(); });

        {
            std::lock_guard<std::mutex> guard(cache_mutex_);
            auto it = cache_.find(residue_three_letter);
            if (it != cache_.end()) return &it->second;
        }

        auto block_it = block_index_.find(residue_three_letter);
        if (block_it == block_index_.end()) return nullptr;

        gemmi::ChemComp cc = gemmi::make_chemcomp_from_block(*block_it->second);

        std::lock_guard<std::mutex> guard(cache_mutex_);
        auto [it, inserted] = cache_.emplace(residue_three_letter, std::move(cc));
        return &it->second;
    }

    const std::string& SourcePath() const { return source_path_; }

private:
    CcdStore() = default;

    // CCD load is an environment requirement, not a soft dependency.
    // Failure to load means the batch run cannot proceed safely (silent
    // degradation across hundreds of proteins is worse than a hard stop).
    // Every load-time failure routes to a controlled fatal abort with a
    // diagnostic; gemmi parser exceptions are caught so they never unwind
    // out of FinalizeConstruction.
    void Load() {
        // Path resolution: RuntimeEnvironment::CcdPath() handles
        // TOML > env var > NMR_DATA_DIR fallback. No file discovery.
        source_path_ = RuntimeEnvironment::CcdPath();

        auto fatal = [this](const std::string& detail) {
            OperationLog::Error("CcdValidator::Load", detail);
            fprintf(stderr,
                "FATAL: CcdValidator::Load — %s\n"
                "  Configure ccd_path in ~/.nmr_tools.toml,\n"
                "  set NMR_CCD_PATH, or run scripts/fetch_ccd.sh to populate\n"
                "  data/ccd/components.cif.gz.\n",
                detail.c_str());
            std::abort();
        };

        if (source_path_.empty()) {
            fatal("no CCD path configured (RuntimeEnvironment::CcdPath empty)");
        }
        if (!fs::exists(source_path_)) {
            fatal("CCD file not found at " + source_path_);
        }

        OperationLog::Info("CcdValidator::Load",
            "reading CCD from " + source_path_);

        // gemmi::cif::read_file is header-only and reads plain (uncompressed)
        // CIF. The CCD ships gzipped from RCSB; scripts/fetch_ccd.sh
        // gunzips after download. Catch parser exceptions and convert to
        // a controlled abort so they never escape into the ctor path.
        try {
            doc_ = gemmi::cif::read_file(source_path_);
        } catch (const std::exception& e) {
            fatal("gemmi failed to parse " + source_path_ + ": " + e.what());
        } catch (...) {
            fatal("gemmi failed to parse " + source_path_ + ": unknown exception");
        }

        block_index_.reserve(doc_.blocks.size());
        for (const gemmi::cif::Block& block : doc_.blocks) {
            // Some monomer libraries prefix names with "comp_"; CCD doesn't.
            const std::string& name = block.name;
            if (name.size() >= 5 && name.compare(0, 5, "comp_") == 0)
                block_index_[name.substr(5)] = &block;
            else
                block_index_[name] = &block;
        }
        loaded_ = true;
        OperationLog::Info("CcdValidator::Load",
            "loaded " + std::to_string(block_index_.size()) + " CCD components");
    }

    std::once_flag load_once_;
    bool loaded_ = false;
    std::string source_path_;
    gemmi::cif::Document doc_;
    std::unordered_map<std::string, const gemmi::cif::Block*> block_index_;
    std::unordered_map<std::string, gemmi::ChemComp> cache_;
    std::mutex cache_mutex_;
};

}  // namespace


// ============================================================================
// DescribeAtom
// ============================================================================

CcdAtomDescription CcdValidator::DescribeAtom(
        const std::string& residue_three_letter,
        const std::string& atom_name) {
    CcdAtomDescription desc;
    desc.atom_name = atom_name;

    const gemmi::ChemComp* cc = CcdStore::Instance().Find(residue_three_letter);
    if (!cc) return desc;

    auto it = cc->find_atom(atom_name);
    if (it == cc->atoms.end()) return desc;

    desc.element_symbol = it->el.name();
    desc.is_aromatic = false;
    for (const auto& bond : cc->rt.bonds) {
        if (bond.id1.atom == atom_name || bond.id2.atom == atom_name) {
            if (bond.aromatic || bond.type == gemmi::BondType::Aromatic) {
                desc.is_aromatic = true;
                break;
            }
        }
    }
    desc.leaving_atom = false;  // CCD's pdbx_leaving_atom_flag not in ChemComp::Atom;
                                // the column is rarely used and we don't need it for diagnostics.
    desc.found = true;
    return desc;
}


// ============================================================================
// Check
// ============================================================================

CcdValidationReport CcdValidator::Check(const Protein& protein) {
    CcdValidationReport report;
    // Trigger load on first call. If load fails, CcdStore aborts with
    // a diagnostic — control does not return here. Reaching the next
    // line means the dictionary is loaded.
    (void)CcdStore::Instance().Find("ALA");
    report.ccd_path = CcdStore::Instance().SourcePath();

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        const AminoAcidType& aatype = res.AminoAcidInfo();
        const std::string three = aatype.three_letter_code;
        report.n_residues_checked++;

        const gemmi::ChemComp* cc = CcdStore::Instance().Find(three);
        if (!cc) {
            CcdValidationFinding f;
            f.kind = CcdValidationFinding::ResidueTypeUnknownToCcd;
            f.residue_three_letter = three;
            f.residue_position = res.sequence_number;
            f.chain_id = res.chain_id;
            f.detail = "CCD has no entry for residue type " + three;
            OperationLog::Warn("CcdValidator::Check", f.detail);
            report.findings.push_back(std::move(f));
            // Still count the atoms for n_atoms_checked even if no CCD match.
            report.n_atoms_checked += static_cast<int>(res.atom_indices.size());
            continue;
        }

        // Build name-keyed maps from each side for fast cross-checks.
        std::unordered_map<std::string, size_t> protein_atoms;
        protein_atoms.reserve(res.atom_indices.size());
        for (size_t ai : res.atom_indices) {
            protein_atoms.emplace(protein.AtomAt(ai).iupac_name.AsString(), ai);
        }

        // --- Protein-side check: every atom present in Protein. ---
        for (size_t ai : res.atom_indices) {
            const Atom& atom = protein.AtomAt(ai);
            const std::string aname = atom.iupac_name.AsString();
            report.n_atoms_checked++;

            auto cc_it = cc->find_atom(aname);
            if (cc_it == cc->atoms.end()) {
                CcdValidationFinding f;
                f.kind = CcdValidationFinding::AtomNotInCcd;
                f.residue_three_letter = three;
                f.residue_position = res.sequence_number;
                f.chain_id = res.chain_id;
                f.atom_name = aname;
                f.detail = "atom " + aname + " not listed in CCD entry " + three;
                OperationLog::Warn("CcdValidator::Check",
                    "residue " + three + std::to_string(res.sequence_number) +
                    ": " + f.detail);
                report.findings.push_back(std::move(f));
                continue;
            }

            // Element cross-check (CCD element symbol vs typed atom element)
            const std::string ccd_elem = cc_it->el.name();
            const std::string protein_elem = SymbolForElement(atom.element);
            if (!ccd_elem.empty() && !protein_elem.empty() && ccd_elem != protein_elem) {
                CcdValidationFinding f;
                f.kind = CcdValidationFinding::ElementMismatch;
                f.residue_three_letter = three;
                f.residue_position = res.sequence_number;
                f.chain_id = res.chain_id;
                f.atom_name = aname;
                f.detail = "Protein element " + protein_elem +
                           ", CCD element " + ccd_elem;
                OperationLog::Warn("CcdValidator::Check",
                    "residue " + three + std::to_string(res.sequence_number) +
                    " atom " + aname + ": " + f.detail);
                report.findings.push_back(std::move(f));
            }
        }

        // --- CCD-side check: every CCD heavy atom present in Protein. ---
        //
        // Skip hydrogens (commonly absent from crystal sources, constant-pH
        // variants, etc.) AND skip atoms that the CCD lists as conditionally
        // present — terminal carboxylate atoms (OXT, OT1, OT2) appear only
        // on the C-terminal residue, terminal H atoms (HXT) only on the
        // protonated form. Without filtering, CCD's complete-residue
        // listing produces ~75 false flags per protein for a 76-residue
        // chain (every residue except the last gets "missing OXT").
        //
        // gemmi's ChemComp::Atom does not expose pdbx_leaving_atom_flag
        // in its data model, so we hardcode the small set of conditional
        // heavy/H names we know about. Keep this list tight — it is
        // safer to ship a noisy diagnostic than to silently mask a real
        // structure error.
        static const std::set<std::string> kConditionalCcdAtoms = {
            "OXT",  // C-terminal carboxylate oxygen
            "OT1",  // alternate name for C-terminal O
            "OT2",  // alternate name for OXT
            "HXT",  // protonated C-terminal carboxyl H
        };

        for (const auto& cc_atom : cc->atoms) {
            if (cc_atom.is_hydrogen()) continue;
            if (kConditionalCcdAtoms.count(cc_atom.id)) continue;
            if (protein_atoms.count(cc_atom.id) == 0) {
                CcdValidationFinding f;
                f.kind = CcdValidationFinding::AtomMissingFromProtein;
                f.residue_three_letter = three;
                f.residue_position = res.sequence_number;
                f.chain_id = res.chain_id;
                f.atom_name = cc_atom.id;
                f.detail = "CCD lists heavy atom " + cc_atom.id +
                           " (" + cc_atom.el.name() + ")" +
                           " not present in Protein residue " + three;
                OperationLog::Warn("CcdValidator::Check",
                    "residue " + three + std::to_string(res.sequence_number) +
                    ": " + f.detail);
                report.findings.push_back(std::move(f));
            }
        }
    }

    OperationLog::Info("CcdValidator::Check",
        "checked " + std::to_string(report.n_atoms_checked) + " atoms in " +
        std::to_string(report.n_residues_checked) + " residues; " +
        std::to_string(report.findings.size()) + " findings");

    return report;
}


// ============================================================================
// WriteCcdValidationDiagnostics
//
// Two files in <output_dir>/ccd_validation/:
//   summary.npy   (1, 5) int32: [n_atoms_checked, n_residues_checked,
//                                 n_findings, n_atom_not_in_ccd,
//                                 n_atom_missing_from_protein]
//   findings.csv  one row per finding (kind, residue, position, chain, atom, detail)
//
// CSV is the natural type for findings — gemmi-side strings (atom names CCD
// lists that are not in Protein) have no atom_index in the typed object model
// to reference, and string packing into NPY is awkward for diagnostic data.
// summary.npy keeps the typed counts queryable from the SDK.
// ============================================================================

static const char* KindString(CcdValidationFinding::Kind k) {
    switch (k) {
        case CcdValidationFinding::AtomNotInCcd:           return "AtomNotInCcd";
        case CcdValidationFinding::AtomMissingFromProtein: return "AtomMissingFromProtein";
        case CcdValidationFinding::ElementMismatch:        return "ElementMismatch";
        case CcdValidationFinding::ResidueTypeUnknownToCcd:return "ResidueTypeUnknownToCcd";
    }
    return "Unknown";
}

static std::string CsvEscape(const std::string& s) {
    bool needs_quotes = false;
    for (char c : s) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            needs_quotes = true;
            break;
        }
    }
    if (!needs_quotes) return s;
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) {
        if (c == '"') out.push_back('"');
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}

int WriteCcdValidationDiagnostics(const CcdValidationReport& report,
                                   const std::string& output_dir) {
    const std::string subdir = output_dir + "/ccd_validation";
    std::error_code ec;
    fs::create_directories(subdir, ec);
    if (ec) {
        OperationLog::Error("WriteCcdValidationDiagnostics",
            "failed to create " + subdir + ": " + ec.message());
        return 0;
    }

    int written = 0;

    // --- summary.npy ---
    int n_atom_not_in_ccd = 0;
    int n_atom_missing_from_protein = 0;
    int n_element_mismatch = 0;
    int n_residue_unknown = 0;
    for (const auto& f : report.findings) {
        switch (f.kind) {
            case CcdValidationFinding::AtomNotInCcd:           n_atom_not_in_ccd++; break;
            case CcdValidationFinding::AtomMissingFromProtein: n_atom_missing_from_protein++; break;
            case CcdValidationFinding::ElementMismatch:        n_element_mismatch++; break;
            case CcdValidationFinding::ResidueTypeUnknownToCcd:n_residue_unknown++; break;
        }
    }
    int32_t summary[5] = {
        report.n_atoms_checked,
        report.n_residues_checked,
        static_cast<int32_t>(report.findings.size()),
        n_atom_not_in_ccd,
        n_atom_missing_from_protein
    };
    if (NpyWriter::WriteInt32(subdir + "/summary.npy", summary, 5))
        written++;

    // --- findings.csv ---
    std::ofstream csv(subdir + "/findings.csv");
    if (csv.is_open()) {
        csv << "kind,residue_three_letter,residue_position,chain_id,atom_name,detail\n";
        for (const auto& f : report.findings) {
            csv << KindString(f.kind) << ","
                << CsvEscape(f.residue_three_letter) << ","
                << f.residue_position << ","
                << CsvEscape(f.chain_id) << ","
                << CsvEscape(f.atom_name) << ","
                << CsvEscape(f.detail) << "\n";
        }
        if (csv.good()) written++;
    } else {
        OperationLog::Error("WriteCcdValidationDiagnostics",
            "could not open " + subdir + "/findings.csv for write");
    }

    // Counts of finding kinds we don't surface in summary.npy still go
    // to the log so they're not lost.
    OperationLog::Info("WriteCcdValidationDiagnostics",
        "summary written to " + subdir +
        " (atom_not_in_ccd=" + std::to_string(n_atom_not_in_ccd) +
        " atom_missing_from_protein=" + std::to_string(n_atom_missing_from_protein) +
        " element_mismatch=" + std::to_string(n_element_mismatch) +
        " residue_unknown=" + std::to_string(n_residue_unknown) + ")");

    return written;
}

}  // namespace nmr

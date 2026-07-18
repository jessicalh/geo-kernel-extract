#include "OrcaRunLoader.h"
#include "Of3Loader.h"
#include "AminoAcidType.h"
#include "NamingRegistry.h"
#include "ChargeSource.h"
#include "AmberChargeResolver.h"
#include "ForceFieldChargeTable.h"
#include "OperationLog.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdio>

namespace fs = std::filesystem;

namespace nmr {

// Internal result type (pre-BuildResult). LoadWithPrmtop produces a Protein;
// the public BuildFromOrca wraps it with charges and net_charge.
struct OrcaLoadInternal {
    std::unique_ptr<Protein> protein;
    bool ok = false;
    std::string error;
};

// File-local neutral carrier: the paths + provenance the heavy structural
// parse (LoadWithPrmtop) actually consumes. It has NO nmr/DFT field —
// structural loading never needs one. Both public entries populate this and
// share LoadWithPrmtop verbatim:
//   - BuildFromOrca  -> NamingSource::OrcaEcho,        "AlphaFold+tleap"
//   - BuildFromOf3   -> NamingSource::Of3PrmtopInput,  "OpenFold+tleap"
// This is the provenance seam that lets the two modes reuse one heavy parse
// without either claiming the other's origin. Introducing it does not change
// BuildFromOrca's observable behavior or public signature.
struct PreparedAmberLoad {
    std::string pdb_path;
    std::string xyz_path;
    std::string prmtop_path;
    std::string tleap_script_path;
    NamingSource naming_source;
    std::string prediction_method;
};

// ============================================================================
// Prmtop section readers (PDB LOADING BOUNDARY)
// ============================================================================

static std::vector<std::string> ReadPrmtopStrings(
        const std::string& path, const std::string& flag_name, int field_width) {
    std::vector<std::string> values;
    std::ifstream in(path);
    if (!in.is_open()) return values;

    std::string line;
    bool in_section = false, found_format = false;

    while (std::getline(in, line)) {
        if (line.find("%FLAG " + flag_name) != std::string::npos) {
            in_section = true; found_format = false; continue;
        }
        if (in_section && !found_format && line.find("%FORMAT") != std::string::npos) {
            found_format = true; continue;
        }
        if (in_section && found_format) {
            if (line.find("%FLAG") != std::string::npos) break;
            for (size_t pos = 0; pos + static_cast<size_t>(field_width) <= line.size();
                 pos += field_width) {
                std::string s = line.substr(pos, field_width);
                size_t start = s.find_first_not_of(' ');
                size_t end = s.find_last_not_of(' ');
                if (start != std::string::npos)
                    values.push_back(s.substr(start, end - start + 1));
            }
        }
    }
    return values;
}

static std::vector<int> ReadPrmtopInts(
        const std::string& path, const std::string& flag_name) {
    std::vector<int> values;
    std::ifstream in(path);
    if (!in.is_open()) return values;

    std::string line;
    bool in_section = false, found_format = false;

    while (std::getline(in, line)) {
        if (line.find("%FLAG " + flag_name) != std::string::npos) {
            in_section = true; found_format = false; continue;
        }
        if (in_section && !found_format && line.find("%FORMAT") != std::string::npos) {
            found_format = true; continue;
        }
        if (in_section && found_format) {
            if (line.find("%FLAG") != std::string::npos) break;
            std::istringstream iss(line);
            int v;
            while (iss >> v) values.push_back(v);
        }
    }
    return values;
}


// ============================================================================
// XYZ reader (PDB LOADING BOUNDARY)
// ============================================================================

struct XyzAtom { std::string element; double x, y, z; };

static std::vector<XyzAtom> ReadXyz(const std::string& path) {
    std::vector<XyzAtom> atoms;
    std::ifstream in(path);
    if (!in.is_open()) return atoms;

    std::string line;
    if (!std::getline(in, line)) return atoms;
    int n = 0;
    std::sscanf(line.c_str(), "%d", &n);
    if (n <= 0) return atoms;

    std::getline(in, line);  // title

    atoms.reserve(n);
    while (std::getline(in, line) && static_cast<int>(atoms.size()) < n) {
        XyzAtom a;
        char elem[4] = {};
        if (std::sscanf(line.c_str(), " %3s %lf %lf %lf",
                         elem, &a.x, &a.y, &a.z) == 4) {
            a.element = elem;
            atoms.push_back(a);
        }
    }
    return atoms;
}


// ============================================================================
// AMBER rst7 / inpcrd reader (PDB LOADING BOUNDARY)
//
// --of3 only. tleap emits coordinates as an AMBER restart (rst7/inpcrd) beside
// the prmtop, so OF3 reads geometry straight from it (no xyz conversion). This
// mirrors ReadXyz as a file-local static and returns the same XyzAtom vector so
// the shared LoadWithPrmtop consumes it unchanged. The ORCA/mutant xyz path
// (ReadXyz above) is untouched.
//
// Format: line 1 = title; line 2 = atom count (+ optional time); then the
// coordinate block as 6F12.7 — six fixed-width 12-char floats per line (2
// atoms) — over exactly ceil(N/2) lines, optionally followed by a velocity
// block and/or a box line (both ignored). rst7 carries NO element symbols
// (they come from the prmtop ATOMIC_NUMBER), so XyzAtom::element is left empty;
// the caller resolves elements from the prmtop exactly as on the xyz path.
//
// FAIL-LOUD hardening (silent-wrong-geometry is unacceptable). On success the
// returned vector has N atoms and `error` is empty; on any violation the vector
// is empty and `error` explains why:
//   - Record shape: read EXACTLY ceil(N/2) coordinate lines, each holding
//     EXACTLY its expected count of 12-char fields (6, or 3 on the final line
//     when N is odd). A coordinate is NEVER pulled from a following velocity/box
//     block to complete a short record.
//   - Over-width columns: a coordinate wider than 12 chars (|coord| overflowing
//     F12.7, e.g. a negative value <= -1000) shifts fixed-column slicing, so the
//     line width no longer matches 12*fields — rejected (never re-sliced).
//   - Full-field numeric consumption: std::stod must consume the entire field
//     (only surrounding whitespace may remain), so a numeric prefix like
//     "1.0garbage" or a Fortran "D"-exponent is rejected, not silently truncated.
//   - Finiteness: every coordinate must be std::isfinite (rejects nan/inf).
// Fields are parsed by fixed 12-char column, NOT whitespace: adjacent
// large-magnitude coordinates (e.g. two values <= -100) touch with no separator.
// ============================================================================
static std::vector<XyzAtom> ReadInpcrd(const std::string& path,
                                       std::string& error) {
    std::vector<XyzAtom> atoms;
    std::ifstream in(path);
    if (!in.is_open()) { error = "cannot open inpcrd: " + path; return atoms; }

    std::string line;
    if (!std::getline(in, line)) { error = "empty inpcrd (no title line)"; return atoms; }
    if (!std::getline(in, line)) { error = "inpcrd truncated (no atom-count line)"; return atoms; }
    int n = 0;
    if (std::sscanf(line.c_str(), "%d", &n) != 1 || n <= 0) {
        error = "inpcrd: unreadable atom count on line 2: '" + line + "'";
        return atoms;
    }

    // rst7 coordinates: 6F12.7 (2 atoms per line) over exactly ceil(N/2) lines.
    // Read that many lines, validating each RECORD's shape — never accumulate
    // "3*N values from anywhere", which is what let a short record borrow from a
    // trailing velocity/box block.
    const long total_values = 3L * n;
    const long coord_lines  = (static_cast<long>(n) + 1) / 2;  // ceil(N/2)
    constexpr size_t kFieldW = 12;

    std::vector<double> coords;
    coords.reserve(total_values);

    for (long li = 0; li < coord_lines; ++li) {
        if (!std::getline(in, line)) {
            error = "inpcrd: coordinate block truncated at line " +
                    std::to_string(li) + " of " + std::to_string(coord_lines) +
                    " (need " + std::to_string(total_values) + " coords for " +
                    std::to_string(n) + " atoms)";
            return {};
        }
        const long remaining = total_values - static_cast<long>(coords.size());
        const long fields_here = std::min<long>(6, remaining);

        // The record must be EXACTLY fields_here fixed-width 12-char columns.
        // Trim only trailing whitespace/CR (leading spaces belong to field 1).
        // A wider record (over-width coordinate) or a shorter/blank record is
        // rejected, never blindly re-sliced or completed from the next block.
        const size_t last = line.find_last_not_of(" \t\r\n");
        const size_t width = (last == std::string::npos) ? 0 : last + 1;
        const size_t expect_w = static_cast<size_t>(fields_here) * kFieldW;
        if (width != expect_w) {
            error = "inpcrd coordinate line " + std::to_string(li) +
                    " has width " + std::to_string(width) + ", expected " +
                    std::to_string(expect_w) + " (" + std::to_string(fields_here) +
                    " x F12.7). An over-width coordinate (|coord| overflowing "
                    "F12.7) or a short/blank record shifts the fixed-column parse.";
            return {};
        }

        for (long f = 0; f < fields_here; ++f) {
            const std::string field = line.substr(static_cast<size_t>(f) * kFieldW, kFieldW);
            size_t consumed = 0;
            double v = 0.0;
            try {
                v = std::stod(field, &consumed);
            } catch (...) {
                error = "inpcrd: non-numeric coordinate field '" + field +
                        "' on line " + std::to_string(li);
                return {};
            }
            // Require FULL-field numeric consumption: only whitespace may follow
            // the number. Catches "1.0garbage", Fortran "D"-exponents, and any
            // partial parse that std::stod would otherwise accept.
            while (consumed < field.size() &&
                   std::isspace(static_cast<unsigned char>(field[consumed]))) ++consumed;
            if (consumed != field.size()) {
                error = "inpcrd: coordinate field '" + field + "' is not fully "
                        "numeric (partial parse) on line " + std::to_string(li);
                return {};
            }
            if (!std::isfinite(v)) {
                error = "inpcrd: non-finite coordinate '" + field + "' on line " +
                        std::to_string(li);
                return {};
            }
            coords.push_back(v);
        }
    }

    if (static_cast<long>(coords.size()) != total_values) {
        error = "inpcrd: read " + std::to_string(coords.size()) +
                " coordinates, expected " + std::to_string(total_values);
        return {};
    }

    atoms.resize(n);
    for (int i = 0; i < n; ++i) {
        atoms[i].x = coords[3 * i + 0];
        atoms[i].y = coords[3 * i + 1];
        atoms[i].z = coords[3 * i + 2];
        // element left empty on purpose: resolved from prmtop ATOMIC_NUMBER.
    }
    error.clear();
    return atoms;
}

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


// ============================================================================
// LoadOrcaRun: WITH prmtop path
// ============================================================================

static OrcaLoadInternal LoadWithPrmtop(const PreparedAmberLoad& load,
                                      const std::vector<XyzAtom>& xyz) {
    OrcaLoadInternal result;

    auto atom_names = ReadPrmtopStrings(load.prmtop_path, "ATOM_NAME", 4);
    auto res_labels = ReadPrmtopStrings(load.prmtop_path, "RESIDUE_LABEL", 4);
    auto res_pointers = ReadPrmtopInts(load.prmtop_path, "RESIDUE_POINTER");
    auto atomic_numbers = ReadPrmtopInts(load.prmtop_path, "ATOMIC_NUMBER");

    if (atom_names.empty() || res_labels.empty() || res_pointers.empty()) {
        result.error = "incomplete prmtop: " + load.prmtop_path;
        return result;
    }

    size_t n_atoms = atom_names.size();
    size_t n_residues = res_labels.size();

    if (xyz.size() != n_atoms) {
        result.error = "XYZ has " + std::to_string(xyz.size()) +
                       " atoms, prmtop has " + std::to_string(n_atoms);
        return result;
    }

    auto& registry = GlobalNamingRegistry();

    // Build Protein using Protein's public API
    auto protein = std::make_unique<Protein>();

    // Build residue ranges (1-based → 0-based)
    struct ResRange { size_t start; size_t end; };
    std::vector<ResRange> ranges(n_residues);
    for (size_t ri = 0; ri < n_residues; ++ri) {
        ranges[ri].start = res_pointers[ri] - 1;
        ranges[ri].end = (ri + 1 < n_residues)
                         ? static_cast<size_t>(res_pointers[ri + 1] - 1)
                         : n_atoms;
    }

    // Add residues
    for (size_t ri = 0; ri < n_residues; ++ri) {
        Residue res;

        // PDB LOADING BOUNDARY: AMBER residue name → canonical
        std::string canonical = registry.ToCanonical(res_labels[ri]);
        if (canonical.empty()) {
            result.error = "unknown residue: " + res_labels[ri] +
                           " at position " + std::to_string(ri);
            return result;
        }
        res.type = AminoAcidFromThreeLetterCode(canonical);
        res.sequence_number = static_cast<int>(ri + 1);
        res.chain_id = "A";

        // Detect protonation variant from AMBER label
        if (res_labels[ri] != canonical) {
            const AminoAcidType& aat = GetAminoAcidType(res.type);
            for (int vi = 0; vi < static_cast<int>(aat.variants.size()); ++vi) {
                if (res_labels[ri] == aat.variants[vi].name) {
                    res.protonation_variant_index = vi;
                    break;
                }
            }
        }

        protein->AddResidue(std::move(res));
    }

    // PDB LOADING BOUNDARY (the prmtop ATOM_NAME surface).
    //
    // prmtop already carries AMBER ff14SB canonical names, so most
    // atom-name strings pass through unchanged. The applicator's
    // pass-through branch in Resolve() returns canonical inputs
    // unchanged. Any non-canonical residual (rare on prmtop loads) is
    // covered by sibling-aware shift rules.
    //
    // The N-terminal cap atom H1 is preserved unchanged: H1 is a
    // distinct cap-only atom in AMBER (kCapNtermCharged: H1, H2, H3
    // on the +1 ammonium nitrogen), NOT a backbone amide H. The
    // canonicality oracle whitelists H1/H2/H3 as cap-only canonical
    // names; no rule shifts them.
    //
    // source = load.naming_source (OrcaEcho for --orca/--mutant,
    // Of3PrmtopInput for --of3 — both the prmtop ATOM_NAME surface, differing
    // only in provenance). No canonicalisation RULE gates on either of those
    // two source values, so the two modes produce identical canonical output;
    // rules tagged AmberFf14SBCanonical (LYN HZ shifts, GLY HA collapse) fire
    // independent of source when the sibling pattern matches.
    const auto& applicator = GlobalNamingApplicator();

    // Per-residue sibling snapshot: collect raw atom names by residue.
    std::vector<std::vector<std::string>> per_residue_raw_names(n_residues);
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        size_t resolved_ri = 0;
        for (size_t ri = 0; ri < n_residues; ++ri) {
            if (ai >= ranges[ri].start && ai < ranges[ri].end) {
                resolved_ri = ri;
                break;
            }
        }
        per_residue_raw_names[resolved_ri].push_back(atom_names[ai]);
    }

    // Canonicalise per residue.
    std::vector<std::vector<std::string>> per_residue_canonical(n_residues);
    for (size_t ri = 0; ri < n_residues; ++ri) {
        const Residue& res_now = protein->ResidueAt(ri);
        const std::vector<std::string> parent_names;
        per_residue_canonical[ri] = applicator.ApplyResidue(
            per_residue_raw_names[ri],
            parent_names,
            res_now.type,
            res_now.protonation_variant_index,
            TerminalState::Internal,
            load.naming_source,
            res_now.sequence_number,
            res_now.chain_id);
    }

    // Per-residue atom counter: tracks position within each residue's
    // canonical-names vector as we walk the global atom list in order.
    std::vector<size_t> per_residue_cursor(n_residues, 0);

    // Add atoms
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto atom = std::make_unique<Atom>();

        // Find residue for this atom (range membership).
        size_t resolved_ri = 0;
        for (size_t ri = 0; ri < n_residues; ++ri) {
            if (ai >= ranges[ri].start && ai < ranges[ri].end) {
                resolved_ri = ri;
                break;
            }
        }
        atom->pdb_atom_name = per_residue_canonical[resolved_ri]
                                  [per_residue_cursor[resolved_ri]++];

        if (ai < atomic_numbers.size())
            atom->element = ElementFromAtomicNumber(atomic_numbers[ai]);
        else
            atom->element = ElementFromSymbol(xyz[ai].element);

        atom->residue_index = resolved_ri;

        // H parent: in tleap ordering, H follows its parent heavy atom
        if (atom->element == Element::H && ai > 0) {
            for (size_t prev = ai; prev > 0; --prev) {
                size_t pi = prev - 1;
                if (pi < atomic_numbers.size() && atomic_numbers[pi] != 1) {
                    atom->parent_atom_index = pi;
                    break;
                }
            }
        }

        size_t idx = protein->AddAtom(std::move(atom));

        // Update residue atom_indices and backbone cache
        Residue& res = protein->MutableResidueAt(
            protein->AtomAt(idx).residue_index);
        res.atom_indices.push_back(idx);

        // Backbone index cache from atom name
        // PDB LOADING BOUNDARY: string name → typed backbone index
        //
        // After canonicalisation, the name is canonical AMBER ff14SB.
        // Backbone amide H (canonical "H"); HN never appears here
        // because it canonicalises to H upstream. H1 is intentionally
        // NOT aliased to H -- it is a distinct N-terminal cap atom.
        // HA2 stays because Gly's backbone canonical IS HA2/HA3.
        const std::string& name = protein->AtomAt(idx).pdb_atom_name;
        if (name == "N" && res.N == Residue::NONE) res.N = idx;
        else if (name == "CA" && res.CA == Residue::NONE) res.CA = idx;
        else if (name == "C" && res.C == Residue::NONE &&
                 protein->AtomAt(idx).element == Element::C &&
                 res.CA != Residue::NONE) res.C = idx;
        else if (name == "O" && res.O == Residue::NONE) res.O = idx;
        else if (name == "H" && res.H == Residue::NONE) res.H = idx;
        else if ((name == "HA" || name == "HA2") &&
                 res.HA == Residue::NONE) res.HA = idx;
        else if (name == "CB" && res.CB == Residue::NONE) res.CB = idx;
    }

    // Build positions and create conformation
    std::vector<Vec3> positions;
    positions.reserve(n_atoms);
    for (const auto& a : xyz)
        positions.push_back(Vec3(a.x, a.y, a.z));

    // Set build context
    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->pdb_source = load.pdb_path;
    ctx->force_field = "ff14SB";
    ctx->protonation_tool = "tleap";
    ctx->prmtop_path = load.prmtop_path;
    ctx->tleap_script_path = load.tleap_script_path;
    protein->SetBuildContext(std::move(ctx));

    // Finalize: backbone indices, bonds, rings — same as PdbFileReader.
    // Must run BEFORE creating the conformation.
    protein->FinalizeConstruction(positions);

    // Create the single conformation from XYZ positions. The prediction
    // label records the honest per-mode provenance supplied at the mode
    // boundary ("AlphaFold+tleap" for --orca/--mutant, "OpenFold+tleap" for
    // --of3), never a hardcoded convention.
    protein->AddPrediction(std::move(positions), load.prediction_method);

    result.protein = std::move(protein);
    result.ok = true;
    return result;
}


// ============================================================================
// LoadOrcaRun: dispatcher
// ============================================================================

BuildResult BuildFromOrca(const OrcaRunFiles& files) {
    BuildResult result;

    OperationLog::Scope scope("BuildFromOrca",
        "pdb=" + files.pdb_path + " xyz=" + files.xyz_path);

    if (!fs::exists(files.xyz_path)) {
        result.error = "XYZ not found: " + files.xyz_path;
        return result;
    }

    auto xyz = ReadXyz(files.xyz_path);
    if (xyz.empty()) {
        result.error = "failed to read XYZ: " + files.xyz_path;
        return result;
    }

    // ORCA / mutant paths require an upstream PRMTOP. Missing or
    // unreadable PRMTOP is a hard load error — there is no fall-through
    // to the ff14SB flat table on the ORCA paths.
    if (files.prmtop_path.empty() || !fs::exists(files.prmtop_path)) {
        result.error = "no prmtop available for " + files.pdb_path +
                       ". --orca/--mutant require an upstream PRMTOP "
                       "(provide files.prmtop_path).";
        return result;
    }

    // Map the ORCA input onto the shared structural carrier. The ORCA
    // provenance is preserved exactly: OrcaEcho atom-name surface and the
    // "AlphaFold+tleap" prediction label. LoadWithPrmtop is the heavy parse
    // shared verbatim with BuildFromOf3; this mapping does not change
    // BuildFromOrca's observable behavior or public signature.
    PreparedAmberLoad load;
    load.pdb_path          = files.pdb_path;
    load.xyz_path          = files.xyz_path;
    load.prmtop_path       = files.prmtop_path;
    load.tleap_script_path = files.tleap_script_path;
    load.naming_source     = NamingSource::OrcaEcho;
    load.prediction_method = "AlphaFold+tleap";

    auto internal = LoadWithPrmtop(load, xyz);
    if (!internal.ok) {
        result.error = internal.error;
        return result;
    }
    result.protein = std::move(internal.protein);

    // Charges from prmtop, via the resolver (branch 1 short-circuit).
    AmberSourceConfig source_config;
    source_config.preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;
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

    OperationLog::Info(LogCharges, "BuildFromOrca",
        "loaded " + std::to_string(result.protein->AtomCount()) + " atoms, " +
        std::to_string(result.protein->ResidueCount()) + " residues, " +
        "net_charge=" + std::to_string(result.net_charge));

    return result;
}


// ============================================================================
// BuildFromOf3: OpenFold-prepared pose (peer of BuildFromOrca)
// ============================================================================
//
// Shares the file-local heavy parse (ReadXyz + LoadWithPrmtop + the prmtop
// section readers) and the Branch-1 charge resolution with BuildFromOrca. The
// only differences are honest OpenFold provenance (Of3PrmtopInput atom-name
// surface, input.prediction_method conformation label) and the deliberate
// ABSENCE of any DFT/NMR input. Defined here so OF3 reuses the same
// file-local heavy parse verbatim; only the thin per-mode wrapper below is
// duplicated so the ORCA entry point above stays observably identical.
BuildResult BuildFromOf3(const Of3Input& input) {
    BuildResult result;

    OperationLog::Scope scope("BuildFromOf3",
        "pdb=" + input.pdb_path + " inpcrd=" + input.inpcrd_path);

    // Honest provenance guard: the conformation records input.prediction_method
    // verbatim. An empty method would silently record false/blank provenance.
    if (input.prediction_method.empty()) {
        result.error = "BuildFromOf3: empty prediction_method; the OF3 mode "
                       "boundary must supply conformation provenance "
                       "(\"OpenFold+tleap\").";
        return result;
    }

    if (!fs::exists(input.inpcrd_path)) {
        result.error = "inpcrd not found: " + input.inpcrd_path;
        return result;
    }

    // Geometry comes straight from the tleap-emitted AMBER restart (rst7),
    // in prmtop atom order. ReadInpcrd fails loud (never silent-wrong geometry)
    // on over-width columns, partial/non-numeric fields, non-finite values, or a
    // short record borrowing from a velocity/box block. The prmtop-vs-coords
    // atom-count equality check is enforced by the shared LoadWithPrmtop below,
    // exactly as on the xyz path.
    std::string inpcrd_err;
    auto coords = ReadInpcrd(input.inpcrd_path, inpcrd_err);
    if (coords.empty()) {
        result.error = "failed to read inpcrd: " + input.inpcrd_path +
                       (inpcrd_err.empty() ? "" : " (" + inpcrd_err + ")");
        return result;
    }

    // OF3 requires an upstream PRMTOP. Missing or unreadable PRMTOP is a hard
    // load error — there is no fall-through to the ff14SB flat table and no
    // runtime tleap preparation. The supplied prmtop is the only charge/radii
    // authority.
    if (input.prmtop_path.empty() || !fs::exists(input.prmtop_path)) {
        result.error = "no prmtop available for OF3 input " + input.inpcrd_path +
                       ". --of3 requires an upstream PRMTOP "
                       "(input.prmtop_path).";
        return result;
    }

    // Map onto the shared structural carrier with honest OpenFold provenance,
    // then run the exact same heavy parse ORCA uses. The carrier's coords-path
    // slot (xyz_path) is not consumed by LoadWithPrmtop — positions flow via
    // the coords vector — so it is left unset for the inpcrd path.
    PreparedAmberLoad load;
    load.pdb_path          = input.pdb_path;
    load.prmtop_path       = input.prmtop_path;
    load.tleap_script_path = input.tleap_script_path;
    load.naming_source     = NamingSource::Of3PrmtopInput;
    load.prediction_method = input.prediction_method;

    auto internal = LoadWithPrmtop(load, coords);
    if (!internal.ok) {
        result.error = internal.error;
        return result;
    }
    result.protein = std::move(internal.protein);

    // Charges from prmtop via the resolver's Branch-1 short-circuit — the same
    // authoritative path ORCA uses. No new charge code, no flat-table fallback.
    AmberSourceConfig source_config;
    source_config.preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;
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

    OperationLog::Info(LogCharges, "BuildFromOf3",
        "loaded " + std::to_string(result.protein->AtomCount()) + " atoms, " +
        std::to_string(result.protein->ResidueCount()) + " residues, " +
        "net_charge=" + std::to_string(result.net_charge));

    return result;
}

}  // namespace nmr

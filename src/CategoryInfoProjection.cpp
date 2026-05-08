#include "CategoryInfoProjection.h"

#include "AminoAcidType.h"
#include "Atom.h"
#include "LegacyAmberTopology.h"
#include "OperationLog.h"
#include "Protein.h"
#include "Residue.h"
#include "SemanticEnums.h"
#include "Types.h"
#include "generated/LegacyAmberSemanticTables.h"

#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace nmr {

// ============================================================================
// Internal state and helpers (anonymous namespace)
// ============================================================================

namespace {

// ── Naming-provenance enum ──────────────────────────────────────
//
// Mirrored as int8 in the NPY. No "DriftLogged" — substrate fields are
// invariants; substrate-vs-atom_nom.tbl disagreement is a chemistry
// error, not a routine fallback.
enum class NamingProvenance : int8_t {
    Match      = 0,
    MissLogged = 1,
};

// ── Parsed atom_nom.tbl row ────────────────────────────────────
//
// We keep only the columns we use: BMRB (col 2) and SC stereo (col 3).
// IUPAC == BMRB for the standard 20 (the BMRB column IS the IUPAC/IUB
// 1970 / Markley 1998 standard).
struct AtomNomRow {
    std::string bmrb_name;   // canonical IUPAC name
    std::string sc_stereo;   // "pro-R" / "pro-S" / "Z" / "E" / ""
};

// ── File-local state ───────────────────────────────────────────
struct State {
    bool configured = false;
    // Per-residue table keyed by (1-letter residue code, AMBER atom name).
    std::map<std::pair<char, std::string>, AtomNomRow> per_residue;
    // Cap-atom table (atom_nom.tbl rows with 1-letter == 'X'); keyed by
    // atom name only. Looked up as fallback after per-residue miss.
    std::map<std::string, AtomNomRow> cap_atoms;
    // "RES:NAME" -> miss count. Includes provenance prefix when needed.
    std::map<std::string, int> miss_log;
};
State& Get() {
    static State s;
    return s;
}

// ── atom_nom.tbl parser ────────────────────────────────────────
//
// Format: tab-separated columns; lines starting with '#' are comments.
// The file's actual columnar layout (verified by `cat -A`) has an
// empty column between the 1-letter residue code and the BMRB name —
// so split-on-tab gives:
//
//   cols[0]: 1-letter residue code (A=Ala, R=Arg, ..., X=universal-cap)
//   cols[1]: empty (legacy padding column)
//   cols[2]: BMRB / IUPAC atom name
//   cols[3]: SC stereo designation (pro-R / pro-S / Z / E / blank)
//   cols[4..]: other naming systems (UCSF, XPLOR, MSI, PDB, SYBYL, MIDAS)
//
// We consume cols 0, 2, 3 only.
//
void ParseAtomNomTable(const fs::path& path, State& s) {
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr,
            "FATAL: CategoryInfoProjection::Configure -- could not open "
            "atom_nom.tbl at \"%s\".\n",
            path.string().c_str());
        std::abort();
    }

    std::string line;
    size_t row_count = 0;
    while (std::getline(in, line)) {
        // Strip comments.
        if (!line.empty() && line[0] == '#') continue;
        // Skip blank lines.
        bool blank = true;
        for (char c : line) {
            if (!std::isspace(static_cast<unsigned char>(c))) { blank = false; break; }
        }
        if (blank) continue;

        // Tokenise on tabs/whitespace runs. Tabs are the documented
        // separator; some rows have variable spacing.
        std::vector<std::string> cols;
        std::string cur;
        bool in_token = false;
        for (char c : line) {
            if (c == '\t') {
                cols.push_back(cur);
                cur.clear();
                in_token = false;
                continue;
            }
            if (std::isspace(static_cast<unsigned char>(c))) {
                if (in_token) {
                    // Inside a token, treat space as separator only if
                    // tabs aren't present elsewhere — but the file
                    // mixes them. Conservative: end token at first space.
                    cols.push_back(cur);
                    cur.clear();
                    in_token = false;
                }
                continue;
            }
            cur.push_back(c);
            in_token = true;
        }
        if (!cur.empty()) cols.push_back(cur);

        if (cols.size() < 3) continue;
        const std::string& res_code = cols[0];
        const std::string& bmrb_name = cols[2];
        std::string sc = cols.size() >= 4 ? cols[3] : "";
        if (bmrb_name.empty()) continue;

        // SC may be blank in the source file — strip.
        if (sc.size() == 1 && std::isspace(static_cast<unsigned char>(sc[0]))) sc.clear();

        AtomNomRow row{bmrb_name, sc};

        if (res_code.size() == 1 && res_code[0] == 'X') {
            s.cap_atoms[bmrb_name] = row;
        } else if (res_code.size() == 1 && std::isupper(static_cast<unsigned char>(res_code[0]))) {
            s.per_residue[{res_code[0], bmrb_name}] = row;
        }
        // Multi-character first-column rows (rare; e.g. "C-S-S") are skipped.

        ++row_count;
    }

    OperationLog::Info("CategoryInfoProjection::Configure",
        "parsed atom_nom.tbl rows=" + std::to_string(row_count) +
        " per_residue=" + std::to_string(s.per_residue.size()) +
        " cap_atoms=" + std::to_string(s.cap_atoms.size()));

    // ── Configure-time deep-consistency sanity check ──────────────
    //
    // Fail loud HERE rather than mid-run if the file is shorter or
    // malformed than expected. The standard 20 produce ~329 per-residue
    // rows (each residue has ~10-25 named atoms). Anything substantially
    // smaller indicates a parser or file-format regression.
    //
    constexpr size_t kMinPerResidue = 280;   // generous floor; real value ~329
    constexpr size_t kMinCapAtoms   = 4;     // H1, H2, H3, OXT minimally
    if (s.per_residue.size() < kMinPerResidue) {
        std::fprintf(stderr,
            "FATAL: CategoryInfoProjection::Configure -- atom_nom.tbl "
            "produced only %zu per-residue rows (expected >= %zu). "
            "File may be truncated or parser regression.\n",
            s.per_residue.size(), kMinPerResidue);
        std::abort();
    }
    if (s.cap_atoms.size() < kMinCapAtoms) {
        std::fprintf(stderr,
            "FATAL: CategoryInfoProjection::Configure -- atom_nom.tbl "
            "produced only %zu cap-atom rows (expected >= %zu).\n",
            s.cap_atoms.size(), kMinCapAtoms);
        std::abort();
    }

    // Spot-check: every standard 1-letter residue code must have at least
    // a backbone N atom keyed for it.
    static constexpr const char* kStandard20 =
        "ARNDCQEGHILKMFPSTWYV";
    for (const char* p = kStandard20; *p; ++p) {
        if (s.per_residue.find({*p, "N"}) == s.per_residue.end()) {
            std::fprintf(stderr,
                "FATAL: CategoryInfoProjection::Configure -- atom_nom.tbl "
                "missing backbone N row for residue '%c'. File integrity "
                "check failed.\n", *p);
            std::abort();
        }
    }
}

// ── 1-letter / 3-letter / variant projections ─────────────────
//
// AminoAcidType already carries `three_letter_code` (canonical IUPAC)
// and `one_letter_code`. Variant residues (CYX/HID/HIE/HIP/ASH/GLH/
// LYN/ARN/TYM) use `variants[i].name` for the AMBER variant string.
//
const char* IupacThreeLetter(AminoAcid type) {
    if (type == AminoAcid::Unknown) return "UNK";
    return GetAminoAcidType(type).three_letter_code;
}
char OneLetter(AminoAcid type) {
    if (type == AminoAcid::Unknown) return 'X';
    return GetAminoAcidType(type).one_letter_code;
}
std::string AmberThreeLetter(AminoAcid type, int variant_index) {
    if (type == AminoAcid::Unknown) return "UNK";
    const AminoAcidType& aat = GetAminoAcidType(type);
    if (variant_index < 0) return aat.three_letter_code;
    if (static_cast<size_t>(variant_index) < aat.variants.size()) {
        return aat.variants[variant_index].name;
    }
    return aat.three_letter_code;
}

// ── AMBER → BMRB cap-name map ────────────────────────────────
//
// atom_nom.tbl indexes cap atoms under their BMRB names, which differ
// from AMBER ff14SB for the C-terminal extras:
//
//   AMBER  →  BMRB
//   ─────     ─────
//   OXT       O''
//   HXT       H''
//   H1        H1   (matches; no projection needed)
//   H2        H2
//   H3        H3
//
// Without this map, AMBER OXT would miss the cap-atom table and emit
// MissLogged for every C-terminus on the fleet — silently breaking
// any BMRB-keyed join. Verified by python/tests/test_atom_nom_consistency.py
// (test_amber_cap_names_NOT_in_table); if atom_nom.tbl ever adds OXT/HXT
// directly the python test will fail and remind us to revisit this.
//
const std::string& AmberToBmrbCapName(const std::string& amber_name) {
    static const std::string oxt_bmrb = "O''";
    static const std::string hxt_bmrb = "H''";
    if (amber_name == "OXT") return oxt_bmrb;
    if (amber_name == "HXT") return hxt_bmrb;
    return amber_name;
}

// ── atom_nom.tbl lookup ──────────────────────────────────────
//
// Returns nullptr on miss; caller logs to miss_log.
//
// Tries per-residue first (e.g. ('A', "HB1") for ALA HB1), then falls
// back to the cap-atoms table for residue-agnostic atoms (H1, OXT, etc.).
// Caps go through AmberToBmrbCapName for AMBER → BMRB name projection.
//
const AtomNomRow* LookupAtomNom(const State& s,
                                  AminoAcid residue_type,
                                  const std::string& atom_name) {
    if (!s.configured) return nullptr;
    const char one_letter = OneLetter(residue_type);
    auto it = s.per_residue.find({one_letter, atom_name});
    if (it != s.per_residue.end()) return &it->second;
    const std::string& bmrb_lookup = AmberToBmrbCapName(atom_name);
    auto it_cap = s.cap_atoms.find(bmrb_lookup);
    if (it_cap != s.cap_atoms.end()) return &it_cap->second;
    return nullptr;
}

// ── S-column packing with strict-length fail-loud guard ───────
//
// Per locked decision: any string column that exceeds its width aborts
// with a diagnostic naming the offending row. PATTERNS §9: fprintf +
// std::abort, no exceptions.
//
void PackString(char* dst, size_t dst_size,
                const std::string& src,
                const char* field_name,
                size_t row_index) {
    if (src.size() >= dst_size) {
        // dst_size includes no null terminator (numpy S<N> = N raw bytes,
        // possibly null-padded). We require src.size() <= dst_size; if
        // src.size() == dst_size we still pack (no terminator needed).
        if (src.size() > dst_size) {
            std::fprintf(stderr,
                "FATAL: CategoryInfoProjection -- field \"%s\" length %zu "
                "exceeds NPY column width %zu at row %zu (value=\"%s\").\n",
                field_name, src.size(), dst_size, row_index, src.c_str());
            std::abort();
        }
    }
    std::memset(dst, 0, dst_size);
    std::memcpy(dst, src.data(), src.size());
}

// ── Packed structured-record schema ─────────────────────────
//
// Layout exactly matches the NPY structured dtype declared at write
// time. C++ packing must match numpy's packed (align=False) layout.
// We use byte-by-byte serialisation through a local buffer to avoid
// any compiler alignment or struct-padding subtleties.
//
constexpr size_t kS8 = 8;
constexpr size_t kS4 = 4;
constexpr size_t kS1 = 1;

constexpr size_t kRecordSize =
    /* atom_index             */ 4 +
    /* residue_index          */ 4 +
    /* element                */ 1 +
    /* amber_atom_name        */ kS8 +
    /* iupac_atom_name        */ kS8 +
    /* bmrb_atom_name         */ kS8 +
    /* amber_residue_3letter  */ kS4 +
    /* iupac_residue_3letter  */ kS4 +
    /* bmrb_residue_3letter   */ kS4 +
    /* residue_1letter        */ kS1 +
    /* residue_type           */ 1 +
    /* residue_variant_index  */ 1 +
    /* terminal_state         */ 1 +
    /* locant                 */ 1 +
    /* branch_outer           */ 1 +
    /* branch_inner           */ 1 +
    /* di_index               */ 1 +
    /* backbone_role          */ 1 +
    /* prochiral              */ 1 +
    /* planar_group           */ 1 +
    /* planar_stereo          */ 1 +
    /* polar_h_kind           */ 1 +
    /* ring_position_primary  */ 1 +
    /* ring_position_secondary*/ 1 +
    /* pseudoatom_kind        */ 1 +
    /* in_super_group         */ 1 +
    /* aromatic               */ 1 +
    /* formal_charge          */ 1 +
    /* is_exchangeable        */ 1 +
    /* iupac_naming_provenance*/ 1 +
    /* bmrb_naming_provenance */ 1;

constexpr const char* kStructuredDtypeDescr =
    "[('atom_index', '<i4'),"
    " ('residue_index', '<i4'),"
    " ('element', 'i1'),"
    " ('amber_atom_name', '|S8'),"
    " ('iupac_atom_name', '|S8'),"
    " ('bmrb_atom_name', '|S8'),"
    " ('amber_residue_3letter', '|S4'),"
    " ('iupac_residue_3letter', '|S4'),"
    " ('bmrb_residue_3letter', '|S4'),"
    " ('residue_1letter', '|S1'),"
    " ('residue_type', 'i1'),"
    " ('residue_variant_index', 'i1'),"
    " ('terminal_state', 'i1'),"
    " ('locant', 'i1'),"
    " ('branch_outer', 'i1'),"
    " ('branch_inner', 'i1'),"
    " ('di_index', 'i1'),"
    " ('backbone_role', 'i1'),"
    " ('prochiral', 'i1'),"
    " ('planar_group', 'i1'),"
    " ('planar_stereo', 'i1'),"
    " ('polar_h_kind', 'i1'),"
    " ('ring_position_primary', 'i1'),"
    " ('ring_position_secondary', 'i1'),"
    " ('pseudoatom_kind', 'i1'),"
    " ('in_super_group', 'i1'),"
    " ('aromatic', 'i1'),"
    " ('formal_charge', 'i1'),"
    " ('is_exchangeable', 'i1'),"
    " ('iupac_naming_provenance', 'i1'),"
    " ('bmrb_naming_provenance', 'i1')]";

// ── Record builder ──────────────────────────────────────────
//
// Build the per-atom structured-record byte buffer in the exact order
// declared in kStructuredDtypeDescr. Returns the buffer.
//
std::vector<unsigned char> BuildRecords(const Protein& protein, State& s) {
    const size_t N = protein.AtomCount();
    std::vector<unsigned char> buf(N * kRecordSize, 0);

    const bool has_substrate = protein.LegacyAmber().HasAtomSemantic();

    for (size_t ai = 0; ai < N; ++ai) {
        unsigned char* row = buf.data() + ai * kRecordSize;
        size_t off = 0;

        const Atom& atom = protein.AtomAt(ai);
        const Residue& res = protein.ResidueAt(atom.residue_index);

        // ── atom_index, residue_index (little-endian int32) ──
        const int32_t atom_idx_val = static_cast<int32_t>(ai);
        std::memcpy(row + off, &atom_idx_val, 4); off += 4;
        const int32_t res_idx_val = static_cast<int32_t>(atom.residue_index);
        std::memcpy(row + off, &res_idx_val, 4); off += 4;

        // ── element (int8) ──
        row[off++] = static_cast<int8_t>(atom.element);

        // ── Atom names ──
        const std::string& amber_name = atom.pdb_atom_name;
        PackString(reinterpret_cast<char*>(row + off), kS8, amber_name,
                   "amber_atom_name", ai);
        off += kS8;

        std::string iupac_name;
        std::string bmrb_name;
        std::string sc_stereo;
        NamingProvenance iupac_prov = NamingProvenance::MissLogged;
        NamingProvenance bmrb_prov  = NamingProvenance::MissLogged;

        if (s.configured) {
            const AtomNomRow* nom = LookupAtomNom(s, res.type, amber_name);
            if (nom) {
                iupac_name = nom->bmrb_name;   // BMRB column == IUPAC for std-20
                bmrb_name  = nom->bmrb_name;
                sc_stereo  = nom->sc_stereo;
                iupac_prov = NamingProvenance::Match;
                bmrb_prov  = NamingProvenance::Match;
            } else {
                // Log once per (residue 3-letter, atom name).
                std::string key = std::string(IupacThreeLetter(res.type)) + ":" + amber_name;
                ++s.miss_log[key];
                // Fallback: emit AMBER name as IUPAC/BMRB so consumers
                // still get a usable string (logged-fallback discipline).
                iupac_name = amber_name;
                bmrb_name  = amber_name;
                sc_stereo  = "";
            }
        } else {
            // No projection configured. Emit AMBER name as fallback.
            iupac_name = amber_name;
            bmrb_name  = amber_name;
        }

        PackString(reinterpret_cast<char*>(row + off), kS8, iupac_name, "iupac_atom_name", ai); off += kS8;
        PackString(reinterpret_cast<char*>(row + off), kS8, bmrb_name,  "bmrb_atom_name",  ai); off += kS8;

        // ── Residue 3-letters ──
        const std::string amber_3 = AmberThreeLetter(res.type, res.protonation_variant_index);
        const std::string iupac_3 = IupacThreeLetter(res.type);
        const std::string bmrb_3  = iupac_3;  // == IUPAC for std-20

        PackString(reinterpret_cast<char*>(row + off), kS4, amber_3, "amber_residue_3letter", ai); off += kS4;
        PackString(reinterpret_cast<char*>(row + off), kS4, iupac_3, "iupac_residue_3letter", ai); off += kS4;
        PackString(reinterpret_cast<char*>(row + off), kS4, bmrb_3,  "bmrb_residue_3letter",  ai); off += kS4;

        // ── Residue 1-letter ──
        char one = OneLetter(res.type);
        PackString(reinterpret_cast<char*>(row + off), kS1, std::string(1, one),
                   "residue_1letter", ai);
        off += kS1;

        // ── Residue type / variant / terminal ──
        row[off++] = static_cast<int8_t>(res.type);
        row[off++] = static_cast<int8_t>(res.protonation_variant_index);
        row[off++] = static_cast<int8_t>(res.terminal_state);

        // ── Mechanical identity / chemistry from substrate ──
        if (has_substrate) {
            const AtomSemanticTable& sem = protein.LegacyAmber().SemanticAt(ai);
            row[off++] = static_cast<int8_t>(sem.locant);
            row[off++] = static_cast<int8_t>(sem.branch.outer);
            row[off++] = static_cast<int8_t>(sem.branch.inner);
            row[off++] = static_cast<int8_t>(sem.di_index);
            row[off++] = static_cast<int8_t>(sem.backbone_role);
            row[off++] = static_cast<int8_t>(sem.prochiral);
            row[off++] = static_cast<int8_t>(sem.planar_group);
            row[off++] = static_cast<int8_t>(sem.planar_stereo);
            row[off++] = static_cast<int8_t>(sem.polar_h);
            row[off++] = static_cast<int8_t>(sem.ring_position.primary.position);
            row[off++] = static_cast<int8_t>(sem.ring_position.secondary.position);
            row[off++] = static_cast<int8_t>(sem.pseudoatom.kind);
            row[off++] = sem.pseudoatom.in_super_group ? 1 : 0;
            row[off++] = sem.aromatic ? 1 : 0;
            row[off++] = static_cast<int8_t>(sem.formal_charge);
            row[off++] = sem.is_exchangeable ? 1 : 0;
        } else {
            // Substrate empty: zero out the 16 substrate-derived bytes.
            std::memset(row + off, 0, 16);
            off += 16;
        }

        // ── Provenance ──
        row[off++] = static_cast<int8_t>(iupac_prov);
        row[off++] = static_cast<int8_t>(bmrb_prov);

        // Sanity: we should have written exactly kRecordSize bytes.
        if (off != kRecordSize) {
            std::fprintf(stderr,
                "FATAL: CategoryInfoProjection record size mismatch at "
                "atom %zu: wrote %zu bytes, expected %zu.\n",
                ai, off, kRecordSize);
            std::abort();
        }
    }

    return buf;
}

// ── NPY 1.0 structured-array writer ─────────────────────────
//
// File layout (numpy NPY format spec):
//
//   bytes 0..5   magic "\x93NUMPY"
//   byte  6      major version (1)
//   byte  7      minor version (0)
//   bytes 8..9   header_len (little-endian uint16)
//   bytes 10..   header (Python dict literal as ASCII), padded with
//                spaces to 64-byte total alignment, ending with '\n'.
//   bytes after  raw record bytes (n_records × kRecordSize)
//
bool WriteStructuredNpy(const fs::path& path,
                          size_t n_records,
                          const std::vector<unsigned char>& bytes) {
    if (bytes.size() != n_records * kRecordSize) {
        std::fprintf(stderr,
            "FATAL: CategoryInfoProjection -- buffer size mismatch: "
            "%zu bytes for %zu records (expected %zu).\n",
            bytes.size(), n_records, n_records * kRecordSize);
        std::abort();
    }

    std::ostringstream hdr;
    hdr << "{'descr': " << kStructuredDtypeDescr
        << ", 'fortran_order': False, 'shape': (" << n_records << ",), }";
    std::string header = hdr.str();

    // Pre-header bytes: 6 magic + 2 version + 2 header_len = 10.
    constexpr size_t kPreHeader = 10;
    // Pad header so total file offset to data is a multiple of 64.
    size_t total = kPreHeader + header.size() + 1;  // +1 for trailing '\n'
    size_t pad = (64 - (total % 64)) % 64;
    header.append(pad, ' ');
    header.push_back('\n');

    if (header.size() > 0xFFFF) {
        std::fprintf(stderr,
            "FATAL: CategoryInfoProjection -- NPY header too large for "
            "v1.0 format (%zu > 65535).\n", header.size());
        std::abort();
    }
    const uint16_t header_len = static_cast<uint16_t>(header.size());

    std::ofstream out(path, std::ios::binary);
    if (!out) {
        OperationLog::Error("CategoryInfoProjection::WriteFeatures",
            "could not open " + path.string() + " for write");
        return false;
    }
    const char magic[] = {'\x93','N','U','M','P','Y'};
    out.write(magic, 6);
    const char ver[] = {'\x01','\x00'};
    out.write(ver, 2);
    out.write(reinterpret_cast<const char*>(&header_len), 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    out.write(reinterpret_cast<const char*>(bytes.data()),
              static_cast<std::streamsize>(bytes.size()));
    return out.good();
}

}  // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

void CategoryInfoProjection::Configure(Config config) {
    State& s = Get();
    s.configured = false;
    s.per_residue.clear();
    s.cap_atoms.clear();
    s.miss_log.clear();

    if (config.atom_nom_tbl.empty()) {
        OperationLog::Info("CategoryInfoProjection::Configure",
            "atom_nom_tbl path empty; projection runs inert");
        return;
    }
    if (!fs::exists(config.atom_nom_tbl)) {
        std::fprintf(stderr,
            "FATAL: CategoryInfoProjection::Configure -- atom_nom_tbl "
            "path \"%s\" does not exist.\n",
            config.atom_nom_tbl.string().c_str());
        std::abort();
    }

    ParseAtomNomTable(config.atom_nom_tbl, s);
    s.configured = true;
}

void CategoryInfoProjection::Reset() {
    State& s = Get();
    s.configured = false;
    s.per_residue.clear();
    s.cap_atoms.clear();
    s.miss_log.clear();
}

bool CategoryInfoProjection::IsActive() {
    return Get().configured;
}

int CategoryInfoProjection::WriteFeatures(const Protein& protein,
                                            const std::string& output_dir) {
    OperationLog::Scope scope("CategoryInfoProjection::WriteFeatures",
        "atoms=" + std::to_string(protein.AtomCount()) +
        " dir=" + output_dir);

    fs::create_directories(output_dir);
    State& s = Get();
    std::vector<unsigned char> records = BuildRecords(protein, s);

    fs::path out_path = fs::path(output_dir) / "atoms_category_info.npy";
    if (!WriteStructuredNpy(out_path, protein.AtomCount(), records)) {
        OperationLog::Error("CategoryInfoProjection::WriteFeatures",
            "write failed: " + out_path.string());
        return 0;
    }
    OperationLog::Info("CategoryInfoProjection::WriteFeatures",
        "wrote " + out_path.string() +
        " (records=" + std::to_string(protein.AtomCount()) +
        " bytes/record=" + std::to_string(kRecordSize) +
        " misses=" + std::to_string(s.miss_log.size()) + ")");
    return 1;
}

// ── Per-atom queries ─────────────────────────────────────────

std::string CategoryInfoProjection::IupacAtomName(const Protein& protein,
                                                    std::size_t atom_index) {
    const State& s = Get();
    if (!s.configured) return {};
    const Atom& atom = protein.AtomAt(atom_index);
    const Residue& res = protein.ResidueAt(atom.residue_index);
    const AtomNomRow* row = LookupAtomNom(s, res.type, atom.pdb_atom_name);
    return row ? row->bmrb_name : std::string{};
}

std::string CategoryInfoProjection::BmrbAtomName(const Protein& protein,
                                                   std::size_t atom_index) {
    return IupacAtomName(protein, atom_index);  // == IUPAC for std-20
}

std::string CategoryInfoProjection::BmrbStereoLabel(const Protein& protein,
                                                      std::size_t atom_index) {
    const State& s = Get();
    if (!s.configured) return {};
    const Atom& atom = protein.AtomAt(atom_index);
    const Residue& res = protein.ResidueAt(atom.residue_index);
    const AtomNomRow* row = LookupAtomNom(s, res.type, atom.pdb_atom_name);
    return row ? row->sc_stereo : std::string{};
}

// ── Per-residue queries ──────────────────────────────────────

std::string CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid type,
                                                              int variant_index) {
    return AmberThreeLetter(type, variant_index);
}

std::string CategoryInfoProjection::IupacResidueThreeLetter(AminoAcid type) {
    return IupacThreeLetter(type);
}

char CategoryInfoProjection::ResidueOneLetter(AminoAcid type) {
    return OneLetter(type);
}

const std::map<std::string, int>& CategoryInfoProjection::MissLog() {
    return Get().miss_log;
}

}  // namespace nmr

#include "TripeptideDftTable.h"

#include "AminoAcidType.h"
#include "OperationLog.h"

#include <libpq-fe.h>

#include <set>
#include <unordered_map>
#include <unordered_set>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace nmr {


// ============================================================================
// Minimal JSONB lexer
//
// Tuned to the tensorcs15 schema. The two payloads we parse are
//   geometry : [{element, x, y, z, atom_idx}, …]
//   tensors  : [{element, atom_idx, isotropic, anisotropy, tensor_3x3,
//                eigenvalues, t2_components}, …]
// Numbers (signed, possibly with decimal/exponent), strings ("..."), and
// the structural tokens [, ], {, }, :, ,. No general JSON support.
// Ported verbatim from the gotham TripeptideDatabase.cpp parser; the
// schema isn't going to grow.
// ============================================================================

namespace {

struct JsonToken {
    enum Type {
        Number, String, ArrayStart, ArrayEnd,
        ObjStart, ObjEnd, Colon, Comma, Eof
    };
    Type type = Eof;
    std::string text;
    double num = 0.0;
};


class JsonLexer {
public:
    explicit JsonLexer(const char* s) : p_(s) {}

    JsonToken Next() {
        while (*p_ && (*p_ == ' ' || *p_ == '\n' ||
                       *p_ == '\r' || *p_ == '\t')) ++p_;
        if (!*p_) return {JsonToken::Eof, "", 0.0};

        switch (*p_) {
            case '[': ++p_; return {JsonToken::ArrayStart, "[", 0.0};
            case ']': ++p_; return {JsonToken::ArrayEnd,   "]", 0.0};
            case '{': ++p_; return {JsonToken::ObjStart,   "{", 0.0};
            case '}': ++p_; return {JsonToken::ObjEnd,     "}", 0.0};
            case ':': ++p_; return {JsonToken::Colon,      ":", 0.0};
            case ',': ++p_; return {JsonToken::Comma,      ",", 0.0};
            case '"': {
                ++p_;
                const char* start = p_;
                while (*p_ && *p_ != '"') ++p_;
                std::string s(start, p_);
                if (*p_ == '"') ++p_;
                return {JsonToken::String, std::move(s), 0.0};
            }
            default: {
                const char* start = p_;
                if (*p_ == '-') ++p_;
                while (*p_ && ((*p_ >= '0' && *p_ <= '9') ||
                               *p_ == '.' || *p_ == 'e' || *p_ == 'E' ||
                               *p_ == '+' || *p_ == '-')) ++p_;
                std::string s(start, p_);
                double v = 0.0;
                if (!s.empty()) {
                    try { v = std::stod(s); } catch (...) { v = 0.0; }
                }
                return {JsonToken::Number, std::move(s), v};
            }
        }
    }

    // Consume tokens until the matching ']' that closes an array we are
    // currently sitting at (ArrayStart already consumed by the caller).
    // Used to skip eigenvalues / t2_components arrays we don't care
    // about (when called from the ParseTensors per-key dispatch).
    void SkipArrayBody() {
        int depth = 1;
        while (depth > 0) {
            JsonToken t = Next();
            if (t.type == JsonToken::ArrayStart) ++depth;
            else if (t.type == JsonToken::ArrayEnd) --depth;
            else if (t.type == JsonToken::Eof) return;
        }
    }

private:
    const char* p_;
};


// Read t2_components: [a, b, c, d, e]. The outer ArrayStart is already
// consumed. Reads exactly 5 numbers; returns them in sphericart order.
std::array<double, 5> ParseT2Array(JsonLexer& lex) {
    std::array<double, 5> out = {};
    for (int i = 0; i < 5; ++i) {
        JsonToken v = lex.Next();
        if (v.type == JsonToken::Number) out[i] = v.num;
        if (i < 4) lex.Next();   // consume comma
    }
    JsonToken end = lex.Next();
    (void)end;  // expect ArrayEnd
    return out;
}


// geometry JSONB: [{element, x, y, z, atom_idx}, …]
struct GeometryEntry {
    int     atom_idx = 0;
    Element element  = Element::Unknown;
    Vec3    pos      = Vec3::Zero();
};


std::vector<GeometryEntry> ParseGeometryJson(const std::string& json_text) {
    std::vector<GeometryEntry> out;
    JsonLexer lex(json_text.c_str());
    lex.Next();  // [

    while (true) {
        JsonToken t = lex.Next();
        if (t.type == JsonToken::ArrayEnd || t.type == JsonToken::Eof) break;
        if (t.type != JsonToken::ObjStart) continue;

        GeometryEntry g;
        while (true) {
            JsonToken key = lex.Next();
            if (key.type == JsonToken::ObjEnd) break;
            lex.Next();  // :
            JsonToken val = lex.Next();
            if      (key.text == "atom_idx") g.atom_idx = static_cast<int>(val.num);
            else if (key.text == "element")  g.element  = ElementFromSymbol(val.text);
            else if (key.text == "x")        g.pos.x()  = val.num;
            else if (key.text == "y")        g.pos.y()  = val.num;
            else if (key.text == "z")        g.pos.z()  = val.num;
            JsonToken sep = lex.Next();
            if (sep.type == JsonToken::ObjEnd) break;
        }
        out.push_back(std::move(g));
        JsonToken sep = lex.Next();
        if (sep.type == JsonToken::ArrayEnd || sep.type == JsonToken::Eof) break;
    }
    return out;
}


// tensors JSONB: [{element, atom_idx, isotropic, anisotropy, tensor_3x3,
//                  eigenvalues, t2_components}, …]
struct TensorEntry {
    int                   atom_idx       = 0;
    Element               element        = Element::Unknown;
    double                isotropic      = 0.0;
    double                anisotropy     = 0.0;
    Mat3                  tensor         = Mat3::Zero();
    std::array<double, 5> t2_components  = {};
};


std::vector<TensorEntry> ParseTensorJson(const std::string& json_text) {
    std::vector<TensorEntry> out;
    JsonLexer lex(json_text.c_str());
    lex.Next();  // [

    while (true) {
        JsonToken t = lex.Next();
        if (t.type == JsonToken::ArrayEnd || t.type == JsonToken::Eof) break;
        if (t.type != JsonToken::ObjStart) continue;

        TensorEntry te;
        while (true) {
            JsonToken key = lex.Next();
            if (key.type == JsonToken::ObjEnd) break;
            lex.Next();  // :

            if (key.text == "atom_idx") {
                JsonToken v = lex.Next();
                te.atom_idx = static_cast<int>(v.num);
            }
            else if (key.text == "element") {
                JsonToken v = lex.Next();
                te.element = ElementFromSymbol(v.text);
            }
            else if (key.text == "isotropic") {
                JsonToken v = lex.Next();
                te.isotropic = v.num;
            }
            else if (key.text == "anisotropy") {
                JsonToken v = lex.Next();
                te.anisotropy = v.num;
            }
            else if (key.text == "tensor_3x3") {
                lex.Next();  // outer [
                for (int row = 0; row < 3; ++row) {
                    lex.Next();  // inner [
                    for (int col = 0; col < 3; ++col) {
                        JsonToken v = lex.Next();
                        te.tensor(row, col) = v.num;
                        if (col < 2) lex.Next();  // comma
                    }
                    lex.Next();  // inner ]
                    if (row < 2) lex.Next();  // comma
                }
                lex.Next();  // outer ]
            }
            else if (key.text == "t2_components") {
                lex.Next();  // [
                te.t2_components = ParseT2Array(lex);
            }
            else {
                // Skip eigenvalues (length 3 array) or any unknown key.
                JsonToken v = lex.Next();
                if (v.type == JsonToken::ArrayStart) lex.SkipArrayBody();
            }

            JsonToken sep = lex.Next();
            if (sep.type == JsonToken::ObjEnd) break;
        }
        out.push_back(std::move(te));
        JsonToken sep = lex.Next();
        if (sep.type == JsonToken::ArrayEnd || sep.type == JsonToken::Eof) break;
    }
    return out;
}


// Round to nearest grid multiple in [-180, 180). 180 normalises to -180
// for the 20° grid (the grid points are -180, -160, …, 160). The 2°
// ALA-baseline grid keeps 180 as-is: it's a real grid point.
int RoundToGrid(double angle, int grid_spacing) {
    int rounded = static_cast<int>(std::round(angle / grid_spacing)) *
                  grid_spacing;
    while (rounded >  180) rounded -= 360;
    while (rounded < -180) rounded += 360;
    if (rounded == 180 && grid_spacing != 2) rounded = -180;
    return rounded;
}


// Number of chi angles the central residue uses. Sets the depth of
// chi columns we include in the WHERE clause.
//
// Per the AminoAcidType per-residue tables (cross-checked against
// per-residue tensorcs15 row counts):
//   0 chi: A, G, P
//   1 chi: S, T, C, V
//   2 chi: L, I, N, D, H, F, Y, W
//   3 chi: E, M, Q
//   4 chi: K, R
int ChiAngleCountForResidue(char letter) {
    switch (letter) {
        case 'A': case 'G': case 'P':                   return 0;
        case 'S': case 'T': case 'C': case 'V':         return 1;
        case 'L': case 'I': case 'N': case 'D':
        case 'H': case 'F': case 'Y': case 'W':         return 2;
        case 'E': case 'M': case 'Q':                   return 3;
        case 'K': case 'R':                             return 4;
        default:                                        return 0;
    }
}


// Merge geometry + tensor entries by atom_idx (NOT by parallel-vector
// position). Both JSONB payloads carry an atom_idx per record; the
// ingest pipeline is expected to emit them in the same order, but
// trusting that without validation is exactly the silent-corruption
// vector that perception is meant to guard against on the downstream
// end — a future ingest script that reordered one array would leave
// every tensor on the wrong atom and same-element swaps would be
// invisible. We key by atom_idx structurally and fail loud on
// mismatched indices, missing tensors for a geometry atom, or duplicate
// atom_idx values.
//
// Same-element swap detection: even when atom_idx values match across
// arrays, we cross-check that the tensor's element matches the
// geometry's element. A divergence means the JSON layer is producing
// inconsistent records — fail loud rather than silently mis-assigning.
std::vector<TripeptideDftAtom> MergeGeometryAndTensors(
        const std::vector<GeometryEntry>& geom,
        const std::vector<TensorEntry>&   tens) {
    std::vector<TripeptideDftAtom> out;
    if (geom.size() != tens.size()) {
        OperationLog::Warn(
            "TripeptideDftTable::MergeGeometryAndTensors",
            "geometry/tensor count mismatch: geom=" +
                std::to_string(geom.size()) +
                " tensor=" + std::to_string(tens.size()) +
                " — proceeding with set intersection");
    }

    // Build atom_idx → TensorEntry map. Duplicate atom_idx in tensors
    // payload is a data-corruption signal; fail loud.
    std::unordered_map<int, const TensorEntry*> tens_by_idx;
    tens_by_idx.reserve(tens.size());
    for (const TensorEntry& t : tens) {
        auto [it, inserted] = tens_by_idx.emplace(t.atom_idx, &t);
        if (!inserted) {
            OperationLog::Error(
                "TripeptideDftTable::MergeGeometryAndTensors",
                "duplicate atom_idx in tensors JSONB: " +
                std::to_string(t.atom_idx) +
                " — silent corruption risk; aborting merge");
            return {};  // empty result; downstream sees rec.atoms.empty()
        }
    }

    out.reserve(geom.size());
    std::unordered_set<int> seen_geom_idx;
    for (const GeometryEntry& g : geom) {
        if (!seen_geom_idx.insert(g.atom_idx).second) {
            OperationLog::Error(
                "TripeptideDftTable::MergeGeometryAndTensors",
                "duplicate atom_idx in geometry JSONB: " +
                std::to_string(g.atom_idx) +
                " — silent corruption risk; aborting merge");
            return {};
        }
        auto it = tens_by_idx.find(g.atom_idx);
        if (it == tens_by_idx.end()) {
            OperationLog::Warn(
                "TripeptideDftTable::MergeGeometryAndTensors",
                "geometry atom_idx=" + std::to_string(g.atom_idx) +
                " has no matching tensor row; skipping");
            continue;
        }
        const TensorEntry& t = *it->second;

        // Same-element cross-check. atom_idx alignment without element
        // agreement means the JSON payloads disagree about chemistry
        // for this atom — a serious upstream bug we want to surface.
        if (g.element != t.element) {
            OperationLog::Error(
                "TripeptideDftTable::MergeGeometryAndTensors",
                "element mismatch at atom_idx=" + std::to_string(g.atom_idx) +
                ": geometry=" + std::to_string(static_cast<int>(g.element)) +
                " tensor=" + std::to_string(static_cast<int>(t.element)) +
                " — silent corruption risk; aborting merge");
            return {};
        }

        TripeptideDftAtom a;
        a.atom_idx         = g.atom_idx;
        a.element          = g.element;
        a.position         = g.pos;
        a.shielding_tensor = t.tensor;
        a.isotropic        = t.isotropic;
        a.anisotropy       = t.anisotropy;
        a.t2_components    = t.t2_components;
        out.push_back(std::move(a));
    }
    return out;
}

}  // anonymous namespace


// ============================================================================
// TripeptideDftTable lifecycle
// ============================================================================

// Redact sensitive kv pairs from a libpq DSN before it reaches the
// operation log. Uses libpq's own PQconninfoParse for structural
// tokenization so the redaction is robust to:
//   - case-insensitive keys (libpq DSN keys are case-insensitive;
//     "Password=..." is valid)
//   - quoted values with embedded whitespace ('secret with spaces')
//   - URI form (postgresql://user:secret@host/db)
//   - non-canonical whitespace (tabs, multiple spaces between kv pairs)
//
// Keys redacted: any key in `kSensitive`. Returns a canonical
// "key='value'" reconstruction of the parsed DSN. On parse failure,
// returns "<dsn redaction failed: ...>" — never falls back to logging
// the raw DSN, since a parser failure is the worst-case scenario for
// silent password leakage.
static std::string RedactDsnForLog(const std::string& dsn) {
    static const std::set<std::string> kSensitive = {
        "password", "passfile",
    };

    char* err = nullptr;
    PQconninfoOption* opts = PQconninfoParse(dsn.c_str(), &err);
    if (!opts) {
        std::string msg = err ? err : "unknown libpq parse error";
        if (err) PQfreemem(err);
        return "<dsn redaction failed: " + msg + ">";
    }

    std::string out;
    for (PQconninfoOption* o = opts; o->keyword; ++o) {
        if (!o->val) continue;
        // libpq sets `val` to NULL for unset options. We only include
        // those the user provided. Output key in lowercase canonical
        // form; libpq emits all known options lowercased already.
        std::string key = o->keyword;
        const std::string val = kSensitive.count(key)
            ? std::string("<redacted>")
            : std::string(o->val);
        if (!out.empty()) out += ' ';
        out += key + "='" + val + "'";
    }
    PQconninfoFree(opts);
    return out;
}


TripeptideDftTable::TripeptideDftTable(const std::string& conn_str) {
    conn_ = PQconnectdb(conn_str.c_str());
    if (PQstatus(conn_) != CONNECTION_OK) {
        std::string err = PQerrorMessage(conn_);
        PQfinish(conn_);
        conn_ = nullptr;
        throw std::runtime_error(
            "TripeptideDftTable: connection failed: " + err);
    }
    OperationLog::Info(LogCalcOther, "TripeptideDftTable",
        "connected to tensorcs15 (conn_str=\"" + RedactDsnForLog(conn_str) + "\")");
}


TripeptideDftTable::~TripeptideDftTable() {
    if (conn_) PQfinish(conn_);
}


bool TripeptideDftTable::IsConnected() const {
    return conn_ != nullptr && PQstatus(conn_) == CONNECTION_OK;
}


// ============================================================================
// QueryNearest
// ============================================================================

TripeptideDftRecord TripeptideDftTable::QueryNearest(
        char residue_letter,
        double phi, double psi,
        double chi1, double chi2,
        double chi3, double chi4,
        int n_chi_axes,
        int his_variant_hint) const {

    TripeptideDftRecord entry;
    if (!conn_) {
        throw std::runtime_error(
            "TripeptideDftTable::QueryNearest: not connected");
    }

    // Grid: 2° for AAA baseline, 20° for everything else.
    const int grid = (residue_letter == 'A') ? 2 : 20;
    const std::string tripeptide =
        std::string("A") + residue_letter + "A";
    const int g_phi = RoundToGrid(phi, grid);
    const int g_psi = RoundToGrid(psi, grid);

    // Caller override of chi depth: when n_chi_axes ≥ 0, use that
    // exact depth (capped to the residue's natural depth so we never
    // include more chi columns than the residue table holds).
    const int natural_n_chi = ChiAngleCountForResidue(residue_letter);
    int n_chi = (n_chi_axes < 0)
                    ? natural_n_chi
                    : std::min(n_chi_axes, natural_n_chi);
    if (n_chi < 0) n_chi = 0;
    const int g_chi1 = (n_chi >= 1) ? RoundToGrid(chi1, 20) : 0;
    const int g_chi2 = (n_chi >= 2) ? RoundToGrid(chi2, 20) : 0;
    const int g_chi3 = (n_chi >= 3) ? RoundToGrid(chi3, 20) : 0;
    const int g_chi4 = (n_chi >= 4) ? RoundToGrid(chi4, 20) : 0;

    // Build the WHERE clause one chi axis at a time so we can fall
    // back from the deepest match. Caller-side fallback (re-run with
    // fewer axes) is handled per-residue by the calculator; here we
    // do the single attempt at the requested depth.
    std::string where = "tripeptide='" + tripeptide + "' AND phi=" +
                        std::to_string(g_phi) + " AND psi=" +
                        std::to_string(g_psi);
    if (n_chi >= 1) where += " AND chi1=" + std::to_string(g_chi1);
    if (n_chi >= 2) where += " AND chi2=" + std::to_string(g_chi2);
    if (n_chi >= 3) where += " AND chi3=" + std::to_string(g_chi3);
    if (n_chi >= 4) where += " AND chi4=" + std::to_string(g_chi4);

    const std::string query =
        "SELECT calc_id, tripeptide, phi, psi, chi1, chi2, chi3, chi4, "
        "n_atoms, frame_type, geometry::text, tensors::text "
        "FROM raw_dft_calculations WHERE " + where + " LIMIT 1";

    PGresult* res = PQexec(conn_, query.c_str());
    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        std::string err = PQerrorMessage(conn_);
        PQclear(res);
        throw std::runtime_error(
            "TripeptideDftTable::QueryNearest: " + err);
    }

    if (PQntuples(res) == 0) {
        PQclear(res);
        OperationLog::Warn(
            "TripeptideDftTable::QueryNearest",
            "no match for " + tripeptide +
                " phi=" + std::to_string(g_phi) +
                " psi=" + std::to_string(g_psi) +
                " n_chi=" + std::to_string(n_chi));
        return entry;
    }

    entry.calc_id   = std::stoi(PQgetvalue(res, 0, 0));
    entry.tripeptide = PQgetvalue(res, 0, 1);
    entry.phi  = std::stoi(PQgetvalue(res, 0, 2));
    entry.psi  = std::stoi(PQgetvalue(res, 0, 3));
    entry.chi1 = PQgetisnull(res, 0, 4) ? 0 : std::stoi(PQgetvalue(res, 0, 4));
    entry.chi2 = PQgetisnull(res, 0, 5) ? 0 : std::stoi(PQgetvalue(res, 0, 5));
    entry.chi3 = PQgetisnull(res, 0, 6) ? 0 : std::stoi(PQgetvalue(res, 0, 6));
    entry.chi4 = PQgetisnull(res, 0, 7) ? 0 : std::stoi(PQgetvalue(res, 0, 7));
    entry.n_atoms = std::stoi(PQgetvalue(res, 0, 8));
    entry.frame_type =
        PQgetisnull(res, 0, 9) ? "" : PQgetvalue(res, 0, 9);

    const std::string geom_json = PQgetvalue(res, 0, 10);
    const std::string tens_json = PQgetvalue(res, 0, 11);
    PQclear(res);

    auto geom_entries = ParseGeometryJson(geom_json);
    auto tens_entries = ParseTensorJson(tens_json);
    entry.atoms = MergeGeometryAndTensors(geom_entries, tens_entries);

    // Perceive typed LarsenTripeptide. Resolve expected central residue
    // from the 1-letter code.
    AminoAcid expected_central = AminoAcid::Unknown;
    for (const auto& t : AllAminoAcidTypes()) {
        if (t.one_letter_code == residue_letter) {
            expected_central = t.index; break;
        }
    }
    entry.larsen = PerceiveLarsenTripeptide(entry, expected_central,
                                              his_variant_hint);
    const std::string perception_tag = entry.larsen.has_value() ? "ok" : "FAILED";

    OperationLog::Info(LogCalcOther,
        "TripeptideDftTable::QueryNearest",
        entry.tripeptide +
            " phi=" + std::to_string(entry.phi) +
            " psi=" + std::to_string(entry.psi) +
            " chi1=" + std::to_string(entry.chi1) +
            " chi2=" + std::to_string(entry.chi2) +
            " chi3=" + std::to_string(entry.chi3) +
            " chi4=" + std::to_string(entry.chi4) +
            " calc_id=" + std::to_string(entry.calc_id) +
            " frame_type=" + entry.frame_type +
            " atoms=" + std::to_string(entry.atoms.size()) +
            " perception=" + perception_tag);

    // Perception self-logs its specific failure reason via
    // OperationLog::Warn("PerceiveLarsenTripeptide", ...) at the
    // exact failure site (amide count, segmentation, graph-iso,
    // slot-cache completeness). We do not double-log here.

    return entry;
}


}  // namespace nmr

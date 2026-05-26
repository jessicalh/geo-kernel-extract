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

#include <nlohmann/json.hpp>

namespace nmr {


namespace {


struct GeometryEntry {
    int     atom_idx = 0;
    Element element  = Element::Unknown;
    Vec3    pos      = Vec3::Zero();
};


std::vector<GeometryEntry> ParseGeometryJson(const std::string& json_text) {
    // Missing coordinates or atom_idx would otherwise default to a
    // plausible but wrong atom; .at() makes the row malformed instead.
    std::vector<GeometryEntry> out;
    try {
        for (const auto& o : nlohmann::json::parse(json_text)) {
            GeometryEntry g;
            g.atom_idx = o.at("atom_idx").get<int>();
            g.element  = ElementFromSymbol(o.at("element").get<std::string>());
            g.pos.x()  = o.at("x").get<double>();
            g.pos.y()  = o.at("y").get<double>();
            g.pos.z()  = o.at("z").get<double>();
            out.push_back(std::move(g));
        }
    } catch (const std::exception& e) {
        OperationLog::Error("TripeptideDftTable::ParseGeometryJson",
            std::string("malformed or incomplete geometry JSONB: ") + e.what() +
            " — aborting; downstream sees no atoms");
        return {};
    }
    return out;
}


struct TensorEntry {
    int                   atom_idx       = 0;
    Element               element        = Element::Unknown;
    double                isotropic      = 0.0;
    double                anisotropy     = 0.0;
    Mat3                  tensor         = Mat3::Zero();
    std::array<double, 5> t2_components  = {};
};


std::vector<TensorEntry> ParseTensorJson(const std::string& json_text) {
    // Missing tensor fields would otherwise become silent zeroes; exact
    // shapes are checked because nlohmann::json::at() does not reject
    // extra rows or components.
    std::vector<TensorEntry> out;
    try {
        for (const auto& o : nlohmann::json::parse(json_text)) {
            TensorEntry te;
            te.atom_idx   = o.at("atom_idx").get<int>();
            te.element    = ElementFromSymbol(o.at("element").get<std::string>());
            te.isotropic  = o.at("isotropic").get<double>();
            te.anisotropy = o.at("anisotropy").get<double>();

            const auto& tensor = o.at("tensor_3x3");
            if (tensor.size() != 3)
                throw std::runtime_error("tensor_3x3 must have exactly 3 rows");
            for (int row = 0; row < 3; ++row) {
                if (tensor.at(row).size() != 3)
                    throw std::runtime_error("tensor_3x3 row must have exactly 3 cols");
                for (int col = 0; col < 3; ++col)
                    te.tensor(row, col) = tensor.at(row).at(col).get<double>();
            }

            const auto& t2 = o.at("t2_components");
            if (t2.size() != 5)
                throw std::runtime_error("t2_components must have exactly 5 entries");
            for (int i = 0; i < 5; ++i)
                te.t2_components[i] = t2.at(i).get<double>();

            out.push_back(std::move(te));
        }
    } catch (const std::exception& e) {
        OperationLog::Error("TripeptideDftTable::ParseTensorJson",
            std::string("malformed or incomplete tensor JSONB: ") + e.what() +
            " — aborting; downstream sees no atoms");
        return {};
    }
    return out;
}


// For the 20-degree grid, 180 normalizes to -180; the 2-degree ALA
// baseline keeps 180 as a real grid point.
int RoundToGrid(double angle, int grid_spacing) {
    int rounded = static_cast<int>(std::round(angle / grid_spacing)) *
                  grid_spacing;
    while (rounded >  180) rounded -= 360;
    while (rounded < -180) rounded += 360;
    if (rounded == 180 && grid_spacing != 2) rounded = -180;
    return rounded;
}


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


// Key by atom_idx rather than JSON array position. A reordered geometry
// or tensor array would otherwise attach shielding tensors to the
// wrong atoms while still passing same-length checks.
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


// Use libpq's own parser so password redaction handles URI DSNs,
// quoted values, case-insensitive keys, and odd whitespace.
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

    const int grid = (residue_letter == 'A') ? 2 : 20;
    const std::string tripeptide =
        std::string("A") + residue_letter + "A";
    const int g_phi = RoundToGrid(phi, grid);
    const int g_psi = RoundToGrid(psi, grid);

    // The override is capped to the residue's natural depth so extra
    // chi columns are never queried for shallow residue tables.
    const int natural_n_chi = ChiAngleCountForResidue(residue_letter);
    int n_chi = (n_chi_axes < 0)
                    ? natural_n_chi
                    : std::min(n_chi_axes, natural_n_chi);
    if (n_chi < 0) n_chi = 0;
    const int g_chi1 = (n_chi >= 1) ? RoundToGrid(chi1, 20) : 0;
    const int g_chi2 = (n_chi >= 2) ? RoundToGrid(chi2, 20) : 0;
    const int g_chi3 = (n_chi >= 3) ? RoundToGrid(chi3, 20) : 0;
    const int g_chi4 = (n_chi >= 4) ? RoundToGrid(chi4, 20) : 0;

    std::string where = "tripeptide='" + tripeptide + "' AND phi=" +
                        std::to_string(g_phi) + " AND psi=" +
                        std::to_string(g_psi);
    if (n_chi >= 1) where += " AND chi1=" + std::to_string(g_chi1);
    if (n_chi >= 2) where += " AND chi2=" + std::to_string(g_chi2);
    if (n_chi >= 3) where += " AND chi3=" + std::to_string(g_chi3);
    if (n_chi >= 4) where += " AND chi4=" + std::to_string(g_chi4);

    // Caller-side chi fallback can make several rows match. calc_id is
    // the stable tie-breaker; otherwise PostgreSQL row order is not
    // replay-stable.
    const std::string query =
        "SELECT calc_id, tripeptide, phi, psi, chi1, chi2, chi3, chi4, "
        "n_atoms, frame_type, geometry::text, tensors::text "
        "FROM raw_dft_calculations WHERE " + where +
        " ORDER BY calc_id ASC LIMIT 1";

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

    return entry;
}


}  // namespace nmr

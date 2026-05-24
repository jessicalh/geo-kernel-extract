#include "BlessCompare.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace nmr::test {

const char* VerdictName(BlessVerdict v) {
    switch (v) {
        case BlessVerdict::Identical:            return "Identical";
        case BlessVerdict::WithinTolerance:      return "WithinTolerance";
        case BlessVerdict::Drifted:              return "Drifted";
        case BlessVerdict::ZeroOutput:           return "ZeroOutput";
        case BlessVerdict::DtypeOrShapeMismatch: return "DtypeOrShapeMismatch";
        case BlessVerdict::ReadFailed:           return "ReadFailed";
    }
    return "Unknown";
}


// ============================================================================
// Minimal NPY reader for the dtypes NpyWriter emits.
//
// NPY v1 layout: 6-byte magic '\x93NUMPY', 1-byte major, 1-byte minor,
// 2-byte little-endian header_len, ASCII Python dict literal of length
// header_len, raw column-major-or-row-major data.
//
// We only consume what NpyWriter writes: fortran_order=False, shapes as
// (n,) or (rows, cols) or (rows, d1, d2, ...). dtype one of <f8 <f4 <i4 <i8.
// Header parsing is deliberately tiny — this is a file format we author.
// ============================================================================

namespace {

struct NpyArray {
    std::string descr;             // e.g. "<f8"
    std::vector<size_t> shape;     // dimensions
    std::vector<double> as_double; // all values cast to double for compare
    size_t element_count = 0;
    bool ok = false;
    std::string error;
};

// Extract the substring between the first occurrence of needle..close in s.
// Returns "" on miss. close defaults to "'".
static std::string ExtractBetween(const std::string& s,
                                  const std::string& needle,
                                  char close) {
    auto pos = s.find(needle);
    if (pos == std::string::npos) return {};
    pos += needle.size();
    auto end = s.find(close, pos);
    if (end == std::string::npos) return {};
    return s.substr(pos, end - pos);
}

// Parse "shape': (1231, 9, )" into [1231, 9]. Returns empty vector on parse
// failure — caller treats that as a structural error.
static std::vector<size_t> ParseShape(const std::string& header) {
    std::vector<size_t> result;
    auto pos = header.find("'shape':");
    if (pos == std::string::npos) return result;
    auto open = header.find('(', pos);
    auto close = header.find(')', pos);
    if (open == std::string::npos || close == std::string::npos
        || close <= open) return result;

    std::string inside = header.substr(open + 1, close - open - 1);
    std::stringstream ss(inside);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        // trim
        while (!tok.empty() && (tok.front() == ' ' || tok.front() == '\t'))
            tok.erase(tok.begin());
        while (!tok.empty() && (tok.back() == ' ' || tok.back() == '\t'))
            tok.pop_back();
        if (tok.empty()) continue;
        char* end = nullptr;
        unsigned long long v = std::strtoull(tok.c_str(), &end, 10);
        if (end == tok.c_str()) { result.clear(); return result; }
        result.push_back(static_cast<size_t>(v));
    }
    return result;
}

static NpyArray ReadNpy(const std::string& path) {
    NpyArray out;
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) {
        out.error = "cannot open " + path;
        return out;
    }

    char magic[6];
    if (!in.read(magic, 6)
        || magic[0] != '\x93' || magic[1] != 'N' || magic[2] != 'U'
        || magic[3] != 'M' || magic[4] != 'P' || magic[5] != 'Y') {
        out.error = "bad NPY magic";
        return out;
    }

    uint8_t major = 0, minor = 0;
    in.read(reinterpret_cast<char*>(&major), 1);
    in.read(reinterpret_cast<char*>(&minor), 1);

    size_t hlen = 0;
    if (major == 1) {
        uint16_t h;
        in.read(reinterpret_cast<char*>(&h), 2);
        hlen = h;
    } else if (major == 2 || major == 3) {
        uint32_t h;
        in.read(reinterpret_cast<char*>(&h), 4);
        hlen = h;
    } else {
        out.error = "unsupported NPY version " + std::to_string(major);
        return out;
    }

    std::string header(hlen, '\0');
    in.read(header.data(), static_cast<std::streamsize>(hlen));

    out.descr = ExtractBetween(header, "'descr': '", '\'');
    out.shape = ParseShape(header);

    if (out.descr.empty() || out.shape.empty()) {
        out.error = "header parse failed: " + header.substr(0, 64);
        return out;
    }

    if (header.find("'fortran_order': False") == std::string::npos) {
        out.error = "fortran_order=True not supported";
        return out;
    }

    size_t n = 1;
    for (size_t d : out.shape) n *= d;
    out.element_count = n;

    // Decode element type and read all into a double buffer for
    // uniform comparison. This buffers the entire array; smoke-test
    // sizes (~10k doubles per file) make that a non-concern.
    out.as_double.resize(n);
    if (out.descr == "<f8") {
        std::vector<double> buf(n);
        in.read(reinterpret_cast<char*>(buf.data()),
                static_cast<std::streamsize>(n * sizeof(double)));
        if (in.gcount() != static_cast<std::streamsize>(n * sizeof(double))) {
            out.error = "short read on <f8 data";
            return out;
        }
        out.as_double = std::move(buf);
    } else if (out.descr == "<f4") {
        std::vector<float> buf(n);
        in.read(reinterpret_cast<char*>(buf.data()),
                static_cast<std::streamsize>(n * sizeof(float)));
        if (in.gcount() != static_cast<std::streamsize>(n * sizeof(float))) {
            out.error = "short read on <f4 data";
            return out;
        }
        for (size_t i = 0; i < n; ++i) out.as_double[i] = buf[i];
    } else if (out.descr == "<i4") {
        std::vector<int32_t> buf(n);
        in.read(reinterpret_cast<char*>(buf.data()),
                static_cast<std::streamsize>(n * sizeof(int32_t)));
        if (in.gcount() != static_cast<std::streamsize>(n * sizeof(int32_t))) {
            out.error = "short read on <i4 data";
            return out;
        }
        for (size_t i = 0; i < n; ++i) out.as_double[i] = buf[i];
    } else if (out.descr == "<i8") {
        std::vector<int64_t> buf(n);
        in.read(reinterpret_cast<char*>(buf.data()),
                static_cast<std::streamsize>(n * sizeof(int64_t)));
        if (in.gcount() != static_cast<std::streamsize>(n * sizeof(int64_t))) {
            out.error = "short read on <i8 data";
            return out;
        }
        for (size_t i = 0; i < n; ++i) out.as_double[i] = static_cast<double>(buf[i]);
    } else {
        out.error = "unsupported dtype " + out.descr;
        return out;
    }

    out.ok = true;
    return out;
}

// Check whether two files are bit-identical without loading both fully.
// Mirrors test_smoke.cpp's FilesIdentical so callers can short-circuit
// the expensive numeric path when nothing has changed at all.
static bool FilesByteIdentical(const std::string& a, const std::string& b) {
    std::ifstream fa(a, std::ios::binary | std::ios::ate);
    std::ifstream fb(b, std::ios::binary | std::ios::ate);
    if (!fa.is_open() || !fb.is_open()) return false;
    if (fa.tellg() != fb.tellg()) return false;
    fa.seekg(0); fb.seekg(0);
    constexpr size_t BUF = 8192;
    char ba[BUF], bb[BUF];
    while (fa && fb) {
        fa.read(ba, BUF);
        fb.read(bb, BUF);
        auto n = fa.gcount();
        if (n != fb.gcount()) return false;
        if (std::memcmp(ba, bb, static_cast<size_t>(n)) != 0) return false;
    }
    return true;
}

}  // namespace


// ============================================================================
// Comparison
// ============================================================================

BlessResult CompareNpy(const std::string& run_path,
                       const std::string& blessed_path,
                       const BlessPolicy& policy) {
    BlessResult out;

    if (FilesByteIdentical(run_path, blessed_path)) {
        out.verdict = BlessVerdict::Identical;
        return out;
    }

    NpyArray run     = ReadNpy(run_path);
    NpyArray blessed = ReadNpy(blessed_path);

    if (!run.ok || !blessed.ok) {
        out.verdict = BlessVerdict::ReadFailed;
        out.diagnostic = "read failure — run: "
                       + (run.ok ? std::string("ok") : run.error)
                       + "; blessed: "
                       + (blessed.ok ? std::string("ok") : blessed.error);
        return out;
    }

    if (run.descr != blessed.descr || run.shape != blessed.shape) {
        out.verdict = BlessVerdict::DtypeOrShapeMismatch;
        std::ostringstream msg;
        msg << "run dtype=" << run.descr << " shape=(";
        for (size_t i = 0; i < run.shape.size(); ++i) {
            if (i) msg << ",";
            msg << run.shape[i];
        }
        msg << "); blessed dtype=" << blessed.descr << " shape=(";
        for (size_t i = 0; i < blessed.shape.size(); ++i) {
            if (i) msg << ",";
            msg << blessed.shape[i];
        }
        msg << ")";
        out.diagnostic = msg.str();
        return out;
    }

    const size_t n = run.element_count;

    // Nonzero-fraction sanity check (per pattern §11 numerical-noise: a
    // calculator that quietly stopped writing should be a loud failure,
    // not a bit-identical-zero pass).
    auto count_nonzero = [](const std::vector<double>& v) {
        size_t c = 0;
        for (double x : v) if (std::fabs(x) > 0.0) ++c;
        return c;
    };
    size_t run_nz     = count_nonzero(run.as_double);
    size_t blessed_nz = count_nonzero(blessed.as_double);
    double run_nz_frac     = n > 0 ? static_cast<double>(run_nz)     / n : 0.0;
    double blessed_nz_frac = n > 0 ? static_cast<double>(blessed_nz) / n : 0.0;

    // Only flag ZeroOutput if blessed had data (so the array IS supposed
    // to carry signal) but run lost most of it. A blessed-empty array
    // staying empty is fine — that's the per-array policy's job to allow.
    if (blessed_nz_frac >= policy.min_nonzero_fraction
        && run_nz_frac < policy.min_nonzero_fraction) {
        out.verdict = BlessVerdict::ZeroOutput;
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "run nonzero %.3f below min %.3f (blessed had %.3f); "
            "run_nz=%zu blessed_nz=%zu of %zu elements",
            run_nz_frac, policy.min_nonzero_fraction, blessed_nz_frac,
            run_nz, blessed_nz, n);
        out.diagnostic = buf;
        return out;
    }

    // Element-wise tolerance check.
    // |Δ| ≤ atol + rtol·|blessed|  → within tolerance.
    size_t exceeded = 0;
    double max_abs_diff = 0.0;
    double max_rel_diff = 0.0;
    size_t max_idx = 0;
    for (size_t i = 0; i < n; ++i) {
        double d = std::fabs(run.as_double[i] - blessed.as_double[i]);
        double thresh = policy.atol + policy.rtol * std::fabs(blessed.as_double[i]);
        if (d > thresh) ++exceeded;
        if (d > max_abs_diff) {
            max_abs_diff = d;
            max_idx = i;
        }
        double denom = std::fabs(blessed.as_double[i]);
        double rel = denom > 0.0 ? d / denom : 0.0;
        if (rel > max_rel_diff) max_rel_diff = rel;
    }

    if (exceeded == 0) {
        out.verdict = BlessVerdict::WithinTolerance;
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "max |Δ|=%.3e max |Δ|/|val|=%.3e across %zu elements "
            "(rtol=%.0e atol=%.0e)",
            max_abs_diff, max_rel_diff, n, policy.rtol, policy.atol);
        out.diagnostic = buf;
        return out;
    }

    out.verdict = BlessVerdict::Drifted;
    char buf[384];
    std::snprintf(buf, sizeof(buf),
        "%zu of %zu elements exceed tolerance "
        "(rtol=%.0e atol=%.0e); max |Δ|=%.3e at index [%zu]; "
        "max |Δ|/|val|=%.3e; run_nz_frac=%.3f blessed_nz_frac=%.3f",
        exceeded, n, policy.rtol, policy.atol,
        max_abs_diff, max_idx, max_rel_diff,
        run_nz_frac, blessed_nz_frac);
    out.diagnostic = buf;
    return out;
}


// ============================================================================
// Policy loader
//
// One-time read of bless_policy.toml. Section parser handles [default]
// and [arrays.<stem>] only; values are key = number with optional
// scientific notation. This is a closed format we author — keeping the
// parser purpose-built avoids pulling a TOML library into the test tree.
// ============================================================================

namespace {

struct PolicyTable {
    BlessPolicy default_policy;
    std::unordered_map<std::string, BlessPolicy> per_array;
    bool loaded = false;
    std::string source_path;
};

static std::mutex g_policy_mutex;
static std::unordered_map<std::string, PolicyTable> g_policy_tables;

static void TrimInPlace(std::string& s) {
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t'))
        s.erase(s.begin());
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t'
                          || s.back() == '\r'))
        s.pop_back();
}

static bool ParseDouble(const std::string& v, double& out) {
    char* end = nullptr;
    double d = std::strtod(v.c_str(), &end);
    if (end == v.c_str()) return false;
    out = d;
    return true;
}

static void ApplyKv(BlessPolicy& p, const std::string& key,
                    const std::string& val) {
    double d;
    if (!ParseDouble(val, d)) return;
    if      (key == "rtol")                 p.rtol = d;
    else if (key == "atol")                 p.atol = d;
    else if (key == "min_nonzero_fraction") p.min_nonzero_fraction = d;
    // unknown keys silently ignored — forward-compatible
}

static const PolicyTable& LoadTable(const std::string& path) {
    std::lock_guard<std::mutex> lk(g_policy_mutex);
    auto it = g_policy_tables.find(path);
    if (it != g_policy_tables.end()) return it->second;

    PolicyTable& t = g_policy_tables[path];
    t.source_path = path;

    if (path.empty() || !fs::exists(path)) {
        // Default table; per-array map empty.
        t.loaded = true;
        return t;
    }

    std::ifstream in(path);
    std::string line;
    std::string section;       // "default" or "arrays.<stem>" or ""
    std::string array_stem;    // populated when section == "arrays.<stem>"
    BlessPolicy current;

    auto flush = [&]() {
        if (section == "default") {
            t.default_policy = current;
        } else if (section.rfind("arrays.", 0) == 0) {
            t.per_array[array_stem] = current;
        }
    };

    while (std::getline(in, line)) {
        // strip comments
        auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        TrimInPlace(line);
        if (line.empty()) continue;

        if (line.front() == '[' && line.back() == ']') {
            // close prior section
            flush();
            section = line.substr(1, line.size() - 2);
            TrimInPlace(section);
            array_stem.clear();
            if (section.rfind("arrays.", 0) == 0) {
                array_stem = section.substr(7);
                // Per-array sections inherit from [default] and override
                // only keys they specify. Convention: place [default]
                // before any [arrays.X] in the file.
                current = t.default_policy;
            } else {
                current = BlessPolicy{};  // will be filled by following kv lines
            }
            continue;
        }

        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);
        TrimInPlace(key);
        TrimInPlace(val);
        ApplyKv(current, key, val);
    }
    flush();

    t.loaded = true;
    return t;
}

}  // namespace


BlessPolicy LoadPolicy(const std::string& policy_toml_path,
                       const std::string& array_stem) {
    const PolicyTable& t = LoadTable(policy_toml_path);
    auto it = t.per_array.find(array_stem);
    if (it != t.per_array.end()) return it->second;
    return t.default_policy;
}

}  // namespace nmr::test

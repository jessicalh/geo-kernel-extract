#include "MopacResult.h"
#include "Protein.h"
#include "RuntimeEnvironment.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <regex>
#include <algorithm>
#include <thread>
#include <cstdlib>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <set>

namespace fs = std::filesystem;

namespace nmr {

// Log channel for MOPAC computation

namespace {

double QuietNaN() {
    return std::numeric_limits<double>::quiet_NaN();
}

std::string Trim(std::string s) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

std::string Upper(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::string> SplitWhitespace(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream ss(line);
    std::string tok;
    while (ss >> tok) out.push_back(tok);
    return out;
}

bool ParseMopacDouble(std::string token, double& value) {
    token = Trim(token);
    if (token.empty()) return false;
    while (!token.empty() && (token.back() == ',' || token.back() == ';'))
        token.pop_back();
    for (char& c : token) {
        if (c == 'D' || c == 'd') c = 'E';
    }
    char* end = nullptr;
    value = std::strtod(token.c_str(), &end);
    return end != token.c_str() && *end == '\0';
}

double ParseMopacDoubleOrNaN(const std::string& token) {
    double v = QuietNaN();
    return ParseMopacDouble(token, v) ? v : QuietNaN();
}

bool IsIntegerToken(const std::string& token) {
    if (token.empty()) return false;
    size_t i = (token[0] == '+' || token[0] == '-') ? 1 : 0;
    if (i == token.size()) return false;
    for (; i < token.size(); ++i)
        if (!std::isdigit(static_cast<unsigned char>(token[i]))) return false;
    return true;
}

std::vector<std::string> SplitLines(const std::string& text) {
    std::vector<std::string> lines;
    std::istringstream ss(text);
    std::string line;
    while (std::getline(ss, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    return lines;
}

std::vector<size_t> LineOffsets(const std::string& text) {
    std::vector<size_t> offsets;
    offsets.push_back(0);
    for (size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '\n' && i + 1 < text.size())
            offsets.push_back(i + 1);
    }
    return offsets;
}

size_t LineNumberForOffset(const std::vector<size_t>& offsets, size_t byte_offset) {
    auto it = std::upper_bound(offsets.begin(), offsets.end(), byte_offset);
    if (it == offsets.begin()) return 0;
    return static_cast<size_t>((it - offsets.begin()) - 1);
}

uint64_t LocalPairKey(size_t a, size_t b) {
    size_t lo = (a < b) ? a : b;
    size_t hi = (a < b) ? b : a;
    return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
}

std::string Hex32(uint32_t v) {
    std::ostringstream os;
    os << std::hex << std::setfill('0') << std::setw(8) << v;
    return os.str();
}

uint32_t Rotr(uint32_t x, uint32_t n) {
    return (x >> n) | (x << (32 - n));
}

std::string Sha256Hex(const std::vector<unsigned char>& bytes) {
    static constexpr uint32_t k[64] = {
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
        0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
        0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
        0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
        0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
        0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
        0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
        0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
        0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
        0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
        0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
        0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
        0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
        0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
        0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
        0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
    };

    uint32_t h[8] = {
        0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
        0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
    };

    std::vector<unsigned char> msg = bytes;
    const uint64_t bit_len = static_cast<uint64_t>(bytes.size()) * 8u;
    msg.push_back(0x80u);
    while ((msg.size() % 64) != 56) msg.push_back(0u);
    for (int i = 7; i >= 0; --i)
        msg.push_back(static_cast<unsigned char>((bit_len >> (i * 8)) & 0xffu));

    for (size_t chunk = 0; chunk < msg.size(); chunk += 64) {
        uint32_t w[64];
        for (int i = 0; i < 16; ++i) {
            const size_t off = chunk + static_cast<size_t>(i) * 4;
            w[i] = (static_cast<uint32_t>(msg[off]) << 24) |
                   (static_cast<uint32_t>(msg[off + 1]) << 16) |
                   (static_cast<uint32_t>(msg[off + 2]) << 8) |
                   static_cast<uint32_t>(msg[off + 3]);
        }
        for (int i = 16; i < 64; ++i) {
            const uint32_t s0 = Rotr(w[i - 15], 7) ^ Rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const uint32_t s1 = Rotr(w[i - 2], 17) ^ Rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }

        uint32_t a = h[0], b = h[1], c = h[2], d = h[3];
        uint32_t e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t S1 = Rotr(e, 6) ^ Rotr(e, 11) ^ Rotr(e, 25);
            const uint32_t ch = (e & f) ^ ((~e) & g);
            const uint32_t temp1 = hh + S1 + ch + k[i] + w[i];
            const uint32_t S0 = Rotr(a, 2) ^ Rotr(a, 13) ^ Rotr(a, 22);
            const uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t temp2 = S0 + maj;
            hh = g;
            g = f;
            f = e;
            e = d + temp1;
            d = c;
            c = b;
            b = a;
            a = temp1 + temp2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d;
        h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }

    std::string out;
    for (uint32_t v : h) out += Hex32(v);
    return out;
}

bool LooksTextArtifact(const fs::path& p) {
    std::string ext = Upper(p.extension().string());
    return ext == ".MOP" || ext == ".OUT" || ext == ".ARC" ||
           ext == ".AUX" || ext == ".HTML" || ext == ".HTM" ||
           ext == ".PDB" || ext == ".TXT" || ext == ".LOG";
}

std::string BytesToText(const std::vector<unsigned char>& bytes) {
    return std::string(reinterpret_cast<const char*>(bytes.data()), bytes.size());
}

bool ReadWholeFile(const fs::path& path, std::vector<unsigned char>& bytes) {
    std::error_code fec;
    const size_t file_size = fs::file_size(path, fec);
    if (fec) return false;
    bytes.resize(file_size);
    FILE* fp = fopen(path.string().c_str(), "rb");
    if (!fp) return false;
    const size_t nread = fread(bytes.data(), 1, file_size, fp);
    fclose(fp);
    bytes.resize(nread);
    return nread == file_size;
}

std::vector<MopacRawArtifact> CaptureSameStemArtifacts(const std::string& stem) {
    std::vector<MopacRawArtifact> artifacts;
    const fs::path stem_path(stem);
    const fs::path dir = stem_path.parent_path();
    const std::string base = stem_path.filename().string();
    std::set<fs::path> paths;

    for (const char* ext : {".mop", ".out", ".arc", ".den", ".aux", ".html", ".htm", ".pdb"}) {
        fs::path p = stem + ext;
        if (fs::exists(p)) paths.insert(fs::absolute(p));
    }

    std::error_code ec;
    if (fs::exists(dir, ec)) {
        for (const auto& de : fs::directory_iterator(dir, ec)) {
            if (ec) break;
            if (!de.is_regular_file(ec)) continue;
            const fs::path p = de.path();
            if (p.stem().string() == base)
                paths.insert(fs::absolute(p));
        }
    }

    for (const auto& p : paths) {
        MopacRawArtifact artifact;
        artifact.filename = p.filename().string();
        artifact.extension = p.extension().string();
        artifact.is_text = LooksTextArtifact(p);
        if (ReadWholeFile(p, artifact.bytes)) {
            artifact.size_bytes = static_cast<std::uint64_t>(artifact.bytes.size());
            artifact.sha256 = Sha256Hex(artifact.bytes);
            artifacts.push_back(std::move(artifact));
        }
    }
    return artifacts;
}

const MopacRawArtifact* FindArtifact(const MopacRunRecord& record, const std::string& extension) {
    const std::string wanted = Upper(extension);
    for (const auto& artifact : record.artifacts) {
        if (Upper(artifact.extension) == wanted) return &artifact;
    }
    return nullptr;
}

bool IsLikelySectionHeading(const std::string& line) {
    const std::string t = Trim(line);
    if (t.size() < 6) return false;
    if (t.find("----") != std::string::npos || t.find("****") != std::string::npos)
        return true;
    int letters = 0, upper = 0, lower = 0, digits = 0;
    for (unsigned char c : t) {
        if (std::isalpha(c)) {
            ++letters;
            if (std::isupper(c)) ++upper;
            if (std::islower(c)) ++lower;
        } else if (std::isdigit(c)) {
            ++digits;
        }
    }
    if (letters < 4 || lower > 0) return false;
    if (digits > letters * 2) return false;
    if (std::isdigit(static_cast<unsigned char>(t.front()))) return false;
    return upper >= letters * 8 / 10;
}

void CaptureGenericOutSections(const std::string& text, MopacRunRecord& record) {
    const auto lines = SplitLines(text);
    const auto offsets = LineOffsets(text);
    std::vector<size_t> starts;
    starts.push_back(0);
    for (size_t i = 0; i < lines.size(); ++i) {
        if (i != 0 && IsLikelySectionHeading(lines[i]))
            starts.push_back(i);
    }
    starts.push_back(lines.size());

    record.sections.clear();
    for (size_t si = 0; si + 1 < starts.size(); ++si) {
        const size_t begin = starts[si];
        const size_t end = starts[si + 1];
        if (begin >= end) continue;
        MopacTextSection section;
        section.id = "out_section_" + std::to_string(record.sections.size());
        section.label = Trim(lines[begin]).empty() ? "OUT preamble" : Trim(lines[begin]);
        section.start_line = begin;
        section.end_line = end;
        section.start_byte = begin < offsets.size() ? offsets[begin] : text.size();
        section.end_byte = end < offsets.size() ? offsets[end] : text.size();
        section.parser_status = "raw_lines_preserved";
        section.heading_lines.push_back(lines[begin]);
        for (size_t i = begin; i < end; ++i)
            section.raw_lines.push_back(lines[i]);
        record.sections.push_back(std::move(section));
    }
}

void CaptureGenericNumericTables(MopacRunRecord& record) {
    for (const auto& section : record.sections) {
        size_t run_start = SIZE_MAX;
        std::vector<std::vector<std::string>> original;
        std::vector<std::vector<double>> numeric;

        auto flush = [&]() {
            if (run_start == SIZE_MAX || original.size() < 2) {
                run_start = SIZE_MAX;
                original.clear();
                numeric.clear();
                return;
            }
            MopacGenericTable table;
            table.id = "out_numeric_table_" + std::to_string(record.tables.size());
            table.section_id = section.id;
            table.axis = "unknown";
            table.label = section.label + " numeric block";
            table.original_strings = std::move(original);
            table.numeric_values = std::move(numeric);
            table.parse_warnings.push_back(
                "generic numeric block: column labels unavailable; raw section is authoritative");
            record.tables.push_back(std::move(table));
            run_start = SIZE_MAX;
            original.clear();
            numeric.clear();
        };

        for (size_t i = 0; i < section.raw_lines.size(); ++i) {
            const auto toks = SplitWhitespace(section.raw_lines[i]);
            std::vector<std::string> raw_row;
            std::vector<double> num_row;
            for (const auto& tok : toks) {
                double v = 0.0;
                if (ParseMopacDouble(tok, v)) {
                    raw_row.push_back(tok);
                    num_row.push_back(v);
                }
            }
            if (num_row.size() >= 2) {
                if (run_start == SIZE_MAX) run_start = i;
                original.push_back(std::move(raw_row));
                numeric.push_back(std::move(num_row));
            } else {
                flush();
            }
        }
        flush();
    }
}

bool IsAuxAssignmentLine(const std::string& line) {
    const std::string t = Trim(line);
    if (t.empty() || t[0] == '#') return false;
    const auto eq = t.find('=');
    if (eq == std::string::npos || eq == 0) return false;
    const std::string lhs = t.substr(0, eq);
    bool has_alpha = false;
    for (unsigned char c : lhs) {
        if (std::islower(c)) return false;
        if (std::isalpha(c)) has_alpha = true;
    }
    return has_alpha;
}

MopacAuxRecord ParseAuxHeader(const std::string& line, size_t line_index) {
    MopacAuxRecord rec;
    rec.source_line = line_index;
    rec.raw_lines.push_back(line);

    std::string t = Trim(line);
    const auto eq = t.find('=');
    std::string lhs = Trim(t.substr(0, eq));
    std::string rhs = Trim(t.substr(eq + 1));

    std::smatch m;
    std::regex count_re(R"(\[(\d+)\]\s*$)");
    if (std::regex_search(lhs, m, count_re)) {
        rec.declared_count = static_cast<size_t>(std::stoull(m[1].str()));
        lhs = Trim(lhs.substr(0, static_cast<size_t>(m.position())));
    }

    const auto colon = lhs.find(':');
    if (colon != std::string::npos) {
        rec.key = Trim(lhs.substr(0, colon));
        rec.unit = Trim(lhs.substr(colon + 1));
    } else {
        rec.key = Trim(lhs);
    }

    if (!rhs.empty()) {
        if (rhs.size() >= 2 && rhs.front() == '"' && rhs.back() == '"') {
            rec.values.push_back(rhs.substr(1, rhs.size() - 2));
        } else {
            auto toks = SplitWhitespace(rhs);
            rec.values.insert(rec.values.end(), toks.begin(), toks.end());
        }
    }
    return rec;
}

void FinalizeAuxRecord(MopacRunRecord& record, MopacAuxRecord& rec) {
    rec.numeric_values.clear();
    for (const auto& value : rec.values) {
        double v = 0.0;
        if (ParseMopacDouble(value, v))
            rec.numeric_values.push_back(v);
    }
    record.aux_record_index[rec.key] = record.aux_records.size();
    record.aux_records.push_back(std::move(rec));
}

void ParseAuxText(const std::string& aux_text, MopacRunRecord& record) {
    const auto lines = SplitLines(aux_text);
    bool active = false;
    MopacAuxRecord current;

    for (size_t i = 0; i < lines.size(); ++i) {
        if (IsAuxAssignmentLine(lines[i])) {
            if (active) FinalizeAuxRecord(record, current);
            current = ParseAuxHeader(lines[i], i);
            active = true;
            continue;
        }
        if (!active) continue;
        current.raw_lines.push_back(lines[i]);
        const std::string t = Trim(lines[i]);
        if (t.empty() || t[0] == '#') continue;
        auto toks = SplitWhitespace(t);
        current.values.insert(current.values.end(), toks.begin(), toks.end());
    }
    if (active) FinalizeAuxRecord(record, current);
}

const MopacAuxRecord* GetAux(const MopacRunRecord& record, const std::string& key) {
    auto it = record.aux_record_index.find(key);
    if (it == record.aux_record_index.end()) return nullptr;
    return &record.aux_records[it->second];
}

std::vector<double> AuxNumeric(const MopacRunRecord& record, const std::string& key) {
    const MopacAuxRecord* rec = GetAux(record, key);
    return rec ? rec->numeric_values : std::vector<double>{};
}

std::vector<std::string> AuxValues(const MopacRunRecord& record, const std::string& key) {
    const MopacAuxRecord* rec = GetAux(record, key);
    return rec ? rec->values : std::vector<std::string>{};
}

void CaptureProvenanceFromAux(MopacRunRecord& record) {
    auto scalar_string = [&](const std::string& key) -> std::string {
        auto vals = AuxValues(record, key);
        return vals.empty() ? "" : vals.front();
    };
    if (record.version.empty()) record.version = scalar_string("MOPAC_VERSION");
    if (record.run_date.empty()) record.run_date = scalar_string("DATE");
    if (record.method.empty()) record.method = scalar_string("METHOD");
    if (record.empirical_formula.empty()) record.empirical_formula = scalar_string("EMPIRICAL_FORMULA");
    if (record.point_group.empty()) record.point_group = scalar_string("POINT_GROUP");
}

void CaptureScalarTermsFromAux(MopacRunRecord& record) {
    for (const auto& rec : record.aux_records) {
        if (rec.numeric_values.size() != 1) continue;
        if (rec.declared_count > 1) continue;
        MopacScalarTerm term;
        term.label = rec.key;
        term.value = rec.numeric_values.front();
        term.unit = rec.unit;
        term.source = "aux";
        term.original_string = rec.values.empty() ? "" : rec.values.front();
        record.scalar_terms.push_back(std::move(term));
    }
}

void CaptureAoBasisFromAux(MopacRunRecord& record) {
    auto atom_index = AuxNumeric(record, "AO_ATOMINDEX");
    auto sym_type = AuxValues(record, "ATOM_SYMTYPE");
    auto zeta = AuxNumeric(record, "AO_ZETA");
    auto pqn = AuxNumeric(record, "ATOM_PQN");
    auto populations = AuxNumeric(record, "AO_POPULATIONS");
    if (populations.empty()) populations = AuxNumeric(record, "AO_CHARGES");
    if (atom_index.empty() && sym_type.empty() && zeta.empty() && pqn.empty()) return;

    const size_t n = std::max({atom_index.size(), sym_type.size(), zeta.size(), pqn.size(), populations.size()});
    record.ao_basis.clear();
    record.ao_basis.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        MopacAOBasisFunction ao;
        ao.ao_index = i;
        if (i < atom_index.size() && atom_index[i] > 0)
            ao.atom_index = static_cast<size_t>(atom_index[i] - 1);
        if (i < sym_type.size()) ao.symmetry_type = sym_type[i];
        if (i < zeta.size()) ao.zeta = zeta[i];
        if (i < pqn.size()) ao.principal_quantum_number = static_cast<int>(pqn[i]);
        if (i < populations.size()) ao.population = populations[i];
        record.ao_basis.push_back(std::move(ao));
    }
}

void AddMatrixBlockFromAux(MopacRunRecord& record,
                           const std::string& key,
                           const std::string& storage,
                           size_t rows,
                           size_t cols) {
    const MopacAuxRecord* rec = GetAux(record, key);
    if (!rec || rec->numeric_values.empty()) return;
    MopacMatrixBlock block;
    block.name = key;
    block.unit = rec->unit;
    block.storage = storage;
    block.rows = rows;
    block.cols = cols;
    block.values = rec->numeric_values;
    block.original_strings = rec->values;
    record.matrix_blocks.push_back(std::move(block));
}

void CaptureMatricesFromAux(MopacRunRecord& record) {
    const size_t nao = record.ao_basis.size();
    const size_t packed = nao == 0 ? 0 : (nao * (nao + 1)) / 2;
    AddMatrixBlockFromAux(record, "OVERLAP_MATRIX", "lower_triangle_packed", packed, 1);
    AddMatrixBlockFromAux(record, "DENSITY_MATRIX", "lower_triangle_packed", packed, 1);
    AddMatrixBlockFromAux(record, "TOTAL_DENSITY_MATRIX", "lower_triangle_packed", packed, 1);
    AddMatrixBlockFromAux(record, "LMO_VECTORS", "row_major_lmo_by_ao", nao, nao);
    AddMatrixBlockFromAux(record, "EIGENVECTORS", "row_major_mo_by_ao", nao, nao);
    AddMatrixBlockFromAux(record, "OVERLAP_COEFFICIENTS", "compressed_coefficients_with_OVERLAP_INDICES", 0, 1);
    AddMatrixBlockFromAux(record, "DENSITY_MATRIX_COEFFICIENTS", "compressed_coefficients_with_DENSITY_MATRIX_INDICES", 0, 1);
    AddMatrixBlockFromAux(record, "MO_COEFFICIENTS", "compressed_coefficients_with_MO_INDICES", 0, 1);

    const MopacAuxRecord* density = GetAux(record, "DENSITY_MATRIX");
    if (!density) density = GetAux(record, "TOTAL_DENSITY_MATRIX");
    const MopacAuxRecord* overlap = GetAux(record, "OVERLAP_MATRIX");
    if (density && overlap &&
        !density->numeric_values.empty() &&
        density->numeric_values.size() == overlap->numeric_values.size()) {
        MopacMatrixBlock block;
        block.name = "MULLIKEN_AO_OVERLAP_POPULATION";
        block.storage = "lower_triangle_packed_density_times_overlap";
        block.rows = density->numeric_values.size();
        block.cols = 1;
        block.values.resize(density->numeric_values.size());
        for (size_t i = 0; i < block.values.size(); ++i)
            block.values[i] = density->numeric_values[i] * overlap->numeric_values[i];
        record.matrix_blocks.push_back(std::move(block));
    }
}

void CaptureMolecularOrbitalsFromAux(MopacRunRecord& record) {
    auto eigenvalues = AuxNumeric(record, "EIGENVALUES");
    if (eigenvalues.empty()) eigenvalues = AuxNumeric(record, "LMO_ENERGY_LEVELS");
    auto occupations = AuxNumeric(record, "MOLECULAR_ORBITAL_OCCUPANCIES");
    auto labels = AuxValues(record, "M.O.SYMMETRY_LABELS");
    const size_t n = std::max({eigenvalues.size(), occupations.size(), labels.size()});
    if (n == 0) return;
    record.molecular_orbitals.resize(n);
    for (size_t i = 0; i < n; ++i) {
        auto& mo = record.molecular_orbitals[i];
        mo.orbital_index = i;
        if (i < eigenvalues.size()) mo.energy = eigenvalues[i];
        if (i < occupations.size()) mo.occupation = occupations[i];
        if (i < labels.size()) mo.label = labels[i];
        if (mo.label.empty() && !eigenvalues.empty()) mo.label = "LMO";
    }
}

void CaptureInputFromMopArtifact(MopacRunRecord& record) {
    const MopacRawArtifact* mop = FindArtifact(record, ".mop");
    if (!mop) return;
    const auto lines = SplitLines(BytesToText(mop->bytes));
    if (!lines.empty()) record.keyword_line = lines[0];
    if (lines.size() > 1) record.title_line = lines[1];
    if (lines.size() > 2) record.comment_line = lines[2];
    for (size_t li = 3; li < lines.size(); ++li) {
        const auto toks = SplitWhitespace(lines[li]);
        if (toks.size() < 7) continue;
        MopacInputAtomLine atom_line;
        atom_line.atom_index = record.input_atoms.size();
        atom_line.element = toks[0];
        atom_line.position = Vec3(ParseMopacDoubleOrNaN(toks[1]),
                                  ParseMopacDoubleOrNaN(toks[3]),
                                  ParseMopacDoubleOrNaN(toks[5]));
        atom_line.x_flag = IsIntegerToken(toks[2]) ? std::stoi(toks[2]) : 0;
        atom_line.y_flag = IsIntegerToken(toks[4]) ? std::stoi(toks[4]) : 0;
        atom_line.z_flag = IsIntegerToken(toks[6]) ? std::stoi(toks[6]) : 0;
        atom_line.raw_line = lines[li];
        record.input_atoms.push_back(std::move(atom_line));
    }
}

void CaptureOutProvenance(const std::string& text, MopacRunRecord& record) {
    const auto lines = SplitLines(text);
    for (const auto& line : lines) {
        const std::string u = Upper(line);
        if (record.version.empty() && u.find("MOPAC") != std::string::npos)
            record.version = Trim(line);
        if (u.find("JOB ENDED NORMALLY") != std::string::npos ||
            u.find("ENDED NORMALLY") != std::string::npos)
            record.termination_status = Trim(line);
        if (u.find("WARNING") != std::string::npos)
            record.warnings.push_back(Trim(line));
        if (u.find("ERROR") != std::string::npos ||
            u.find("ABNORMAL") != std::string::npos)
            record.errors.push_back(Trim(line));
        if (u.find("SCF") != std::string::npos ||
            u.find("MOZYME") != std::string::npos ||
            u.find("ITER") != std::string::npos ||
            u.find("CONVERG") != std::string::npos)
            record.convergence_records.push_back(Trim(line));
        if (u.find("CPU") != std::string::npos ||
            u.find("TIME") != std::string::npos ||
            u.find("SECONDS") != std::string::npos ||
            u.find("MINUTES") != std::string::npos)
            record.timing_records.push_back(Trim(line));
    }
}

void CaptureScalarTermsFromOut(const std::string& text, MopacRunRecord& record) {
    std::regex scalar_re(R"(^\s*([A-Z][A-Z0-9 ._/\-\(\)]+?)\s*=\s*([+-]?\d+(?:\.\d*)?(?:[DEde][+-]?\d+)?)\s*([^=\n\r]*))");
    const auto lines = SplitLines(text);
    for (const auto& line : lines) {
        std::smatch m;
        if (!std::regex_search(line, m, scalar_re)) continue;
        MopacScalarTerm term;
        term.label = Trim(m[1].str());
        term.value = ParseMopacDoubleOrNaN(m[2].str());
        term.unit = Trim(m[3].str());
        term.source = "out";
        term.original_string = Trim(m[2].str());
        record.scalar_terms.push_back(std::move(term));
    }
}

void CaptureDipoleRowsFromOut(const std::string& text, MopacRunRecord& record) {
    auto pos = text.find("DIPOLE           X");
    if (pos == std::string::npos) return;
    std::istringstream section(text.substr(pos));
    std::string line;
    std::getline(section, line);
    while (std::getline(section, line)) {
        const auto toks = SplitWhitespace(line);
        if (toks.size() < 5) {
            if (!Trim(line).empty()) continue;
            if (!record.dipole_rows.empty()) break;
            continue;
        }
        double x = 0.0, y = 0.0, z = 0.0, total = 0.0;
        if (!ParseMopacDouble(toks[1], x) ||
            !ParseMopacDouble(toks[2], y) ||
            !ParseMopacDouble(toks[3], z) ||
            !ParseMopacDouble(toks[4], total)) {
            continue;
        }
        MopacDipoleRow row;
        row.label = toks[0];
        row.x = x;
        row.y = y;
        row.z = z;
        row.total = total;
        row.raw_line = line;
        record.dipole_rows.push_back(std::move(row));
        if (Upper(toks[0]) == "SUM") break;
    }
}

void CaptureAtomPopulationsFromOut(const std::string& text, size_t natoms, MopacRunRecord& record) {
    record.atom_populations.resize(natoms);
    for (size_t i = 0; i < natoms; ++i)
        record.atom_populations[i].atom_index = i;

    auto pos = text.find("NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS");
    if (pos == std::string::npos) return;

    std::istringstream section(text.substr(pos));
    std::string line;
    std::vector<std::string> headers;
    std::getline(section, line);
    for (int i = 0; i < 2 && std::getline(section, line); ++i) {
        auto toks = SplitWhitespace(line);
        if (!toks.empty()) headers.insert(headers.end(), toks.begin(), toks.end());
    }

    while (std::getline(section, line)) {
        if (line.empty() || line.find("DIPOLE") != std::string::npos) break;
        const auto toks = SplitWhitespace(line);
        if (toks.size() < 5 || !IsIntegerToken(toks[0])) continue;
        const size_t atom = static_cast<size_t>(std::stoul(toks[0]) - 1);
        if (atom >= natoms) continue;
        auto& pop = record.atom_populations[atom];
        pop.element = toks[1];
        pop.column_labels = headers;
        pop.original_cells = toks;
        pop.net_charge = ParseMopacDoubleOrNaN(toks[2]);
        pop.electron_density = ParseMopacDoubleOrNaN(toks[3]);
        pop.s_population = ParseMopacDoubleOrNaN(toks[4]);
        if (toks.size() > 5) pop.p_population = ParseMopacDoubleOrNaN(toks[5]);
        if (toks.size() > 6) pop.d_population = ParseMopacDoubleOrNaN(toks[6]);
        if (toks.size() > 7) pop.f_population = ParseMopacDoubleOrNaN(toks[7]);
    }
}

void CaptureAtomicOrbitalPopulationsFromOut(const std::string& text, size_t natoms, MopacRunRecord& record) {
    auto pos = text.find("ATOMIC ORBITAL ELECTRON POPULATIONS");
    if (pos == std::string::npos) return;

    constexpr size_t kCols = 9;
    std::istringstream section(text.substr(pos));
    std::string line;
    bool saw_header = false;
    bool saw_row = false;
    std::vector<std::string> labels = {"s", "px", "py", "pz", "x^2-y^2", "xz", "z^2", "yz", "xy"};

    while (std::getline(section, line)) {
        const auto toks = SplitWhitespace(line);
        if (toks.empty()) {
            if (saw_row) break;
            continue;
        }
        const std::string u = Upper(line);
        if (u.find("BOND ORDERS") != std::string::npos ||
            u.find("VALENCIES") != std::string::npos ||
            u.find("MULLIKEN POPULATION") != std::string::npos) {
            break;
        }
        if (!saw_header) {
            if (toks.size() >= 2 && Upper(toks[0]) == "ATOM") {
                labels.clear();
                for (size_t i = 1; i < toks.size(); ++i)
                    labels.push_back(toks[i]);
                labels.resize(kCols);
                saw_header = true;
            }
            continue;
        }
        if (toks.size() < 3 || !IsIntegerToken(toks[0])) continue;
        const size_t atom = static_cast<size_t>(std::stoul(toks[0]) - 1);
        if (atom >= natoms) continue;

        MopacAtomicOrbitalPopulation row;
        row.atom_index = atom;
        row.element = toks[1];
        row.column_labels = labels;
        row.original_cells = toks;
        for (size_t c = 0; c < kCols && c + 2 < toks.size(); ++c)
            row.populations[c] = ParseMopacDoubleOrNaN(toks[c + 2]);
        record.atomic_orbital_populations.push_back(std::move(row));
        saw_row = true;
    }
}

void CapturePrintedBondOrdersFromOut(const std::string& text, size_t natoms, MopacRunRecord& record) {
    auto pos = text.find("BOND ORDERS");
    if (pos == std::string::npos) return;
    const auto offsets = LineOffsets(text);
    const size_t start_line = LineNumberForOffset(offsets, pos);
    std::istringstream section(text.substr(pos));
    std::string line;
    std::getline(section, line);
    std::getline(section, line);

    std::regex line_re(R"(\s*(\d+)\s+([A-Za-z]+)\s+\(([+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?)\)(.*))");
    std::regex pair_re(R"((\d+)\s+([A-Za-z]+)\s+([+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?))");
    size_t source_line = start_line + 2;
    size_t row_order = 0;
    bool have_row = false;
    size_t row_atom = 0;
    size_t current_row_order = 0;
    std::string row_element;
    double row_valency = QuietNaN();

    auto consume_pairs = [&](const std::string& rest, const std::string& raw, size_t line_no) {
        if (!have_row) return;
        auto begin = std::sregex_iterator(rest.begin(), rest.end(), pair_re);
        auto end = std::sregex_iterator();
        for (auto it = begin; it != end; ++it) {
            const int nbr_atom_0 = std::stoi((*it)[1].str()) - 1;
            if (nbr_atom_0 < 0 || static_cast<size_t>(nbr_atom_0) >= natoms) continue;
            MopacPrintedBondEntry entry;
            entry.printed_entry_index = record.printed_bond_entries.size();
            entry.row_order = current_row_order;
            entry.row_atom = row_atom;
            entry.row_element = row_element;
            entry.row_valency = row_valency;
            entry.neighbour_atom = static_cast<size_t>(nbr_atom_0);
            entry.neighbour_element = (*it)[2].str();
            entry.order = ParseMopacDoubleOrNaN((*it)[3].str());
            entry.source_line = line_no;
            entry.raw_line = raw;
            record.printed_bond_entries.push_back(std::move(entry));
        }
    };

    while (std::getline(section, line)) {
        ++source_line;
        if (Trim(line).empty()) {
            have_row = false;
            continue;
        }
        if (line.find("BOND ORDERS") != std::string::npos) continue;
        if (line.find("CARTESIAN") != std::string::npos ||
            line.find("ATOMIC") != std::string::npos ||
            line.find("MULLIKEN") != std::string::npos ||
            line.find("JOB ENDED") != std::string::npos ||
            line.find("COMPUTATION") != std::string::npos) break;

        std::smatch lm;
        if (!std::regex_match(line, lm, line_re)) {
            consume_pairs(line, line, source_line);
            continue;
        }
        const int row_atom_1 = std::stoi(lm[1].str());
        const int row_atom_0 = row_atom_1 - 1;
        if (row_atom_0 < 0 || static_cast<size_t>(row_atom_0) >= natoms) continue;

        have_row = true;
        current_row_order = row_order++;
        row_atom = static_cast<size_t>(row_atom_0);
        row_element = lm[2].str();
        row_valency = ParseMopacDoubleOrNaN(lm[3].str());
        const std::string rest = lm[4].str();

        if (row_atom < record.atom_populations.size()) {
            auto& pop = record.atom_populations[row_atom];
            pop.mopac_valency = row_valency;
            if (pop.element.empty()) pop.element = row_element;
        }

        consume_pairs(rest, line, source_line);
    }
}

void CaptureUniqueBondOrders(MopacRunRecord& record) {
    struct Accum {
        double max_order = -std::numeric_limits<double>::infinity();
        double sum_order = 0.0;
        size_t count = 0;
        std::vector<size_t> refs;
    };
    std::map<std::pair<size_t, size_t>, Accum> accum;
    for (const auto& entry : record.printed_bond_entries) {
        if (!std::isfinite(entry.order)) continue;
        const size_t a = std::min(entry.row_atom, entry.neighbour_atom);
        const size_t b = std::max(entry.row_atom, entry.neighbour_atom);
        auto& acc = accum[{a, b}];
        acc.max_order = std::max(acc.max_order, entry.order);
        acc.sum_order += entry.order;
        acc.count += 1;
        acc.refs.push_back(entry.printed_entry_index);
    }

    record.unique_bond_orders.clear();
    for (auto& [key, acc] : accum) {
        MopacUniqueBondOrder u;
        u.atom_a = key.first;
        u.atom_b = key.second;
        u.max_order = acc.max_order;
        u.mean_order = acc.count ? acc.sum_order / static_cast<double>(acc.count) : QuietNaN();
        u.printed_entry_indices = std::move(acc.refs);
        record.unique_bond_orders.push_back(std::move(u));
    }
}

void CaptureMoBondingContributionsFromOut(const std::string& text, MopacRunRecord& record) {
    const auto lines = SplitLines(text);
    bool in_section = false;
    std::string section_id;
    std::vector<std::string> section_lines;
    std::regex row_re(R"(^\s*(\d+)\s+(.*?)([+-]?\d*\.?\d+(?:[DEde][+-]?\d+)?)\s*$)");

    auto is_mo_bonding_heading = [](const std::string& line) {
        const std::string u = Upper(line);
        const bool orbital = u.find("ORBITAL") != std::string::npos ||
                             u.find("M.O") != std::string::npos ||
                             u.find("MO ") != std::string::npos ||
                             u.find("LMO") != std::string::npos ||
                             u.find("LOCALIZED") != std::string::npos;
        const bool bonding = u.find("BOND") != std::string::npos ||
                             u.find("VALENC") != std::string::npos ||
                             u.find("CONTRIBUT") != std::string::npos;
        return orbital && bonding;
    };

    auto flush_section = [&]() {
        if (!in_section) return;
        MopacGenericTable table;
        table.id = "mo_bonding_contribution_" + std::to_string(record.tables.size());
        table.section_id = section_id;
        table.axis = "mo";
        table.label = "MO/LMO bonding contribution table";
        table.column_labels = {"orbital_index", "bonding_contribution"};
        table.parse_warnings.push_back(
            "heading-bounded parser: original section retained in output sections; "
            "numeric rows use the final numeric token as the contribution");

        for (const auto& raw : section_lines) {
            std::smatch m;
            if (!std::regex_match(raw, m, row_re)) continue;
            const size_t idx = static_cast<size_t>(std::stoul(m[1].str()) - 1);
            const double contribution = ParseMopacDoubleOrNaN(m[3].str());
            if (idx >= record.molecular_orbitals.size())
                record.molecular_orbitals.resize(idx + 1);
            auto& mo = record.molecular_orbitals[idx];
            mo.orbital_index = idx;
            mo.bonding_contribution = contribution;
            if (mo.label.empty()) mo.label = Trim(m[2].str());
            table.row_labels.push_back(std::to_string(idx));
            table.original_strings.push_back({m[1].str(), Trim(m[2].str()), m[3].str()});
            table.numeric_values.push_back({static_cast<double>(idx), contribution});
        }
        if (!table.numeric_values.empty())
            record.tables.push_back(std::move(table));
        section_lines.clear();
        in_section = false;
    };

    for (size_t i = 0; i < lines.size(); ++i) {
        const std::string u = Upper(lines[i]);
        if (u.find("BOND ORDERS") != std::string::npos) {
            flush_section();
            return;
        }
        if (!in_section && is_mo_bonding_heading(lines[i])) {
            in_section = true;
            section_id = "out_line_" + std::to_string(i);
            section_lines.clear();
            continue;
        }
        if (!in_section) continue;
        if (Trim(lines[i]).empty()) {
            if (!section_lines.empty()) flush_section();
            continue;
        }
        if (IsLikelySectionHeading(lines[i]) && !is_mo_bonding_heading(lines[i])) {
            flush_section();
            continue;
        }
        section_lines.push_back(lines[i]);
    }
    flush_section();
}

void PopulateFullRecordFromParsedLegacy(MopacRunRecord& record,
                                        const std::vector<double>& charges,
                                        const std::vector<double>& s_pop,
                                        const std::vector<double>& p_pop,
                                        const std::vector<double>& project_valencies) {
    if (record.atom_populations.size() < charges.size())
        record.atom_populations.resize(charges.size());
    for (size_t i = 0; i < charges.size(); ++i) {
        auto& pop = record.atom_populations[i];
        const bool had_typed_population_row =
            !pop.original_cells.empty() ||
            std::isfinite(pop.net_charge) ||
            std::isfinite(pop.electron_density);
        pop.atom_index = i;
        if (!std::isfinite(pop.net_charge)) pop.net_charge = charges[i];
        if (!std::isfinite(pop.s_population)) pop.s_population = s_pop[i];
        if (!std::isfinite(pop.p_population) &&
            (!had_typed_population_row || p_pop[i] != 0.0)) {
            pop.p_population = p_pop[i];
        }
        if (i < project_valencies.size()) pop.project_valency = project_valencies[i];
    }
}

void BuildTopologyBondRecords(const Protein& protein, MopacRunRecord& record) {
    std::map<uint64_t, size_t> unique_by_key;
    for (size_t ui = 0; ui < record.unique_bond_orders.size(); ++ui) {
        const auto& u = record.unique_bond_orders[ui];
        unique_by_key[LocalPairKey(u.atom_a, u.atom_b)] = ui;
    }

    record.topology_bond_records.clear();
    record.topology_bond_records.reserve(protein.BondCount());
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const auto& bond = protein.BondAt(bi);
        const uint64_t key = LocalPairKey(bond.atom_index_a, bond.atom_index_b);
        MopacTopologyBondOrderRecord row;
        row.bond_index = bi;
        row.atom_a = bond.atom_index_a;
        row.atom_b = bond.atom_index_b;
        auto it = unique_by_key.find(key);
        if (it != unique_by_key.end()) {
            const auto& u = record.unique_bond_orders[it->second];
            row.present = true;
            row.order = u.max_order;
            row.unique_pair_index = it->second;
            row.printed_entry_count = u.printed_entry_indices.size();
            record.unique_bond_orders[it->second].topology_bond_index = bi;
        } else {
            row.present = false;
            row.order = QuietNaN();
            row.absence_reason =
                "not present in MOPAC printed bond-order rows; below print threshold, not printed by MOZYME, or absent from electronic bond graph";
        }
        record.topology_bond_records.push_back(std::move(row));
    }
}

void FillElementLabelsFromProtein(const Protein& protein, MopacRunRecord& record) {
    if (record.atom_populations.size() < protein.AtomCount())
        record.atom_populations.resize(protein.AtomCount());
    for (size_t i = 0; i < protein.AtomCount(); ++i) {
        auto& pop = record.atom_populations[i];
        pop.atom_index = i;
        if (pop.element.empty()) pop.element = SymbolForElement(protein.AtomAt(i).element);
    }
}

}  // anonymous namespace


// ============================================================================
// Write .mop input file for MOPAC PM7+MOZYME 1SCF.
//
// Keywords: PM7 MOZYME 1SCF CHARGE=N ALLBONDS MULLIK AUX(...) LARGE LET GEO-OK THREADS=<threads>
//   PM7:     semiempirical Hamiltonian
//   MOZYME:  linear-scaling SCF (~45s for 889 atoms)
//   1SCF:    single-point, no geometry optimisation
//   CHARGE:  net formal charge
//   ALLBONDS: print Wiberg bond order matrix including hydrogen/low-order rows
//   MULLIK:  print Mulliken population analysis
//   AUX:     write the self-describing machine-readable result stream
//   LARGE:   request full population/MO detail where MOPAC supports it
//   LET:     allow unusual geometries without aborting
//   GEO-OK:  suppress geometry warnings
//   THREADS: OpenMP parallelism
// ============================================================================

static std::string WriteMopFile(const Protein& protein,
                                 const ProteinConformation& conf,
                                 int net_charge,
                                 int threads,
                                 const std::string& path) {
    std::ofstream out(path);
    if (!out.is_open()) return "Cannot open " + path + " for writing";

    out << "PM7 MOZYME 1SCF CHARGE=" << net_charge
        << " ALLBONDS MULLIK AUX(PRECISION=9,MOS=99999) LARGE LET GEO-OK THREADS="
        << threads << "\n";

    std::string name = protein.BuildContext().pdb_source;
    if (!name.empty()) name = fs::path(name).stem().string();
    if (name.empty()) name = "protein";
    out << name << " " << conf.AtomCount() << " atoms\n";
    out << "\n";

    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        Vec3 pos = conf.PositionAt(i);
        std::string elem = SymbolForElement(protein.AtomAt(i).element);
        char line[128];
        snprintf(line, sizeof(line), "  %-2s  %14.8f 0  %14.8f 0  %14.8f 0\n",
                 elem.c_str(), pos.x(), pos.y(), pos.z());
        out << line;
    }
    out.close();
    return "";
}


// ============================================================================
// Parse MOPAC .out file for Mulliken charges, orbital populations,
// bond orders, heat of formation.
//
// This is the text-parsing path needed because the subprocess writes
// the .out file. The parsed fields feed the current mopac_* NPY schema.
// ============================================================================

struct MopacParsed {
    std::vector<double> charges;
    std::vector<double> s_pop;
    std::vector<double> p_pop;
    std::vector<MopacBondOrder> bond_orders;
    double heat_of_formation = 0.0;
    Vec3 dipole = Vec3::Zero();
    MopacRunRecord run_record;
    bool success = false;
    std::string error;
};

static MopacParsed ParseMopacOutput(const std::string& out_path, size_t natoms) {
    MopacParsed result;
    result.charges.resize(natoms, 0.0);
    result.s_pop.resize(natoms, 0.0);
    result.p_pop.resize(natoms, 0.0);

    // Read the .out file with C stdio after the subprocess exits, so
    // MOPAC's Fortran I/O has closed its file handles.
    std::string text;
    {
        std::error_code fec;
        size_t file_size = fs::file_size(out_path, fec);
        if (fec || file_size == 0) {
            result.error = "Cannot stat " + out_path;
            return result;
        }
        text.resize(file_size);
        FILE* fp = fopen(out_path.c_str(), "rb");
        if (!fp) {
            result.error = "Cannot open " + out_path;
            return result;
        }
        size_t nread = fread(text.data(), 1, file_size, fp);
        fclose(fp);
        text.resize(nread);
    }

    // Check for normal termination
    if (text.find("JOB ENDED NORMALLY") == std::string::npos &&
        text.find("ended normally") == std::string::npos) {
        result.error = "MOPAC abnormal termination";
        return result;
    }

    // --- Heat of formation ---
    {
        std::regex hof_re(R"(FINAL HEAT OF FORMATION\s*=\s*([-\d.]+))");
        std::smatch m;
        if (std::regex_search(text, m, hof_re)) {
            result.heat_of_formation = std::stod(m[1].str());
        }
    }

    // --- NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS ---
    // This section has: index, element, charge, ?, s_pop, p_pop
    {
        auto pos = text.find("NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS");
        if (pos != std::string::npos) {
            std::istringstream section(text.substr(pos));
            std::string line;
            // Skip header lines (title + blank + column headers)
            std::getline(section, line);  // title
            std::getline(section, line);  // blank or underline
            std::getline(section, line);  // column headers

            while (std::getline(section, line)) {
                if (line.empty() || line.find("DIPOLE") != std::string::npos) break;
                std::istringstream ss(line);
                int idx;
                std::string elem;
                double charge, dummy, sp;
                if (ss >> idx >> elem >> charge >> dummy >> sp) {
                    size_t i = static_cast<size_t>(idx - 1);
                    if (i < natoms) {
                        result.charges[i] = charge;
                        result.s_pop[i] = sp;
                        // p-Pop is absent for hydrogen (only s orbital)
                        double pp = 0.0;
                        if (ss >> pp) {
                            result.p_pop[i] = pp;
                        }
                    }
                }
            }
        }
    }

    // --- BOND ORDERS ---
    // Format: "  4  C     (3.119)    17  C 1.014     6  C 0.968 ..."
    {
        auto pos = text.find("BOND ORDERS");
        if (pos != std::string::npos) {
            std::istringstream section(text.substr(pos));
            std::string line;
            std::getline(section, line);  // "BOND ORDERS" title
            std::getline(section, line);  // blank or header

            std::regex line_re(R"(\s+(\d+)\s+\w+\s+\([\d.]+\)(.*))");
            std::regex pair_re(R"((\d+)\s+\w+\s+([\d.]+))");

            while (std::getline(section, line)) {
                if (line.empty()) continue;
                if (line.find("BOND ORDERS") != std::string::npos) continue;
                // Stop at next section
                if (line.find("CARTESIAN") != std::string::npos ||
                    line.find("ATOMIC") != std::string::npos ||
                    line.find("COMPUTATION") != std::string::npos) break;

                std::smatch lm;
                if (std::regex_match(line, lm, line_re)) {
                    int atom_i = std::stoi(lm[1].str()) - 1;
                    std::string rest = lm[2].str();

                    auto begin = std::sregex_iterator(rest.begin(), rest.end(), pair_re);
                    auto end = std::sregex_iterator();
                    for (auto it = begin; it != end; ++it) {
                        int atom_j = std::stoi((*it)[1].str()) - 1;
                        double bo = std::stod((*it)[2].str());
                        if (bo > 0.01 &&
                            atom_i >= 0 && static_cast<size_t>(atom_i) < natoms &&
                            atom_j >= 0 && static_cast<size_t>(atom_j) < natoms) {
                            MopacBondOrder mbo;
                            mbo.atom_a = static_cast<size_t>(atom_i);
                            mbo.atom_b = static_cast<size_t>(atom_j);
                            mbo.wiberg_order = bo;
                            result.bond_orders.push_back(mbo);
                        }
                    }
                }
            }
        }
    }

    // --- DIPOLE (Debye) ---
    // Table format:
    //           DIPOLE           X         Y         Z       TOTAL
    //  POINT-CHG.        ...
    //  HYBRID            ...
    //  SUM               dx        dy        dz        total
    {
        auto pos = text.find("DIPOLE           X");
        if (pos != std::string::npos) {
            std::istringstream section(text.substr(pos));
            std::string line;
            std::getline(section, line);  // header
            while (std::getline(section, line)) {
                if (line.find("SUM") != std::string::npos) {
                    std::istringstream ss(line);
                    std::string label;
                    double dx, dy, dz;
                    if (ss >> label >> dx >> dy >> dz) {
                        result.dipole = Vec3(dx, dy, dz);
                    }
                    break;
                }
            }
        }
    }

    // --- Full-run additive capture ---
    // The legacy fields above remain the compatibility projection. The
    // following pass preserves the broader output and parses typed views from
    // raw .out/.aux records without feeding back into the legacy regexes.
    {
        std::string stem = out_path;
        const auto dot = stem.rfind('.');
        if (dot != std::string::npos) stem = stem.substr(0, dot);
        result.run_record.atom_count = natoms;
        result.run_record.temp_stem = stem;
        result.run_record.artifacts = CaptureSameStemArtifacts(stem);
        CaptureInputFromMopArtifact(result.run_record);

        CaptureOutProvenance(text, result.run_record);
        CaptureGenericOutSections(text, result.run_record);
        CaptureGenericNumericTables(result.run_record);
        CaptureScalarTermsFromOut(text, result.run_record);
        CaptureDipoleRowsFromOut(text, result.run_record);
        CaptureAtomPopulationsFromOut(text, natoms, result.run_record);
        CaptureAtomicOrbitalPopulationsFromOut(text, natoms, result.run_record);
        CapturePrintedBondOrdersFromOut(text, natoms, result.run_record);
        CaptureUniqueBondOrders(result.run_record);

        if (const MopacRawArtifact* aux = FindArtifact(result.run_record, ".aux")) {
            ParseAuxText(BytesToText(aux->bytes), result.run_record);
            CaptureProvenanceFromAux(result.run_record);
            CaptureScalarTermsFromAux(result.run_record);
            CaptureAoBasisFromAux(result.run_record);
            CaptureMolecularOrbitalsFromAux(result.run_record);
            CaptureMatricesFromAux(result.run_record);
        }
        CaptureMoBondingContributionsFromOut(text, result.run_record);
    }

    result.success = true;
    return result;
}


// ============================================================================
// Compute: write .mop, call MOPAC as a subprocess, parse .out
// ============================================================================

std::unique_ptr<MopacResult> MopacResult::Compute(
        ProteinConformation& conf,
        int net_charge,
        int threads) {

    const Protein& protein = conf.ProteinRef();
    const size_t natoms = conf.AtomCount();

    // Resolve thread count: 0 = auto = 3/4 of hardware concurrency.
    // This machine is dedicated to this work — give MOPAC most of the cores.
    if (threads <= 0) {
        int hw = static_cast<int>(std::thread::hardware_concurrency());
        threads = std::max(4, (hw * 3) / 4);
    }

    // Set OMP environment for MOPAC's Fortran/OpenMP internals.
    // OMP_STACKSIZE: Fortran stack overflow at ~500 atoms with the default
    // 8 MB. 2 GB is generous and safe for proteins up to ~4000 atoms.
    // OMP_NUM_THREADS: belt-and-suspenders with the THREADS keyword.
    // These are process-global but MOPAC is the only OpenMP consumer here
    // and we run one protein at a time (no concurrent mozyme_scf calls).
    setenv("OMP_STACKSIZE", "2G", 1);
    setenv("OMP_NUM_THREADS", std::to_string(threads).c_str(), 1);

    OperationLog::Scope scope("MopacResult::Compute",
        "atoms=" + std::to_string(natoms) +
        " charge=" + std::to_string(net_charge) +
        " threads=" + std::to_string(threads));

    // Generate guid-unique temp file paths
    std::string protein_name = protein.BuildContext().pdb_source;
    if (!protein_name.empty())
        protein_name = fs::path(protein_name).stem().string();
    if (protein_name.empty()) protein_name = "protein";

    std::string mop_path = RuntimeEnvironment::TempFilePath(protein_name, "mopac.mop");

    // Write .mop input
    std::string err = WriteMopFile(protein, conf, net_charge, threads, mop_path);
    if (!err.empty()) {
        OperationLog::Error("MopacResult::Compute", err);
        return nullptr;
    }

    // Call MOPAC as a subprocess. We link libmopac but use the binary
    // because run_mopac_from_input() leaves Fortran I/O buffers unflushed —
    // the .out file is incomplete when we try to read it. The subprocess
    // exits cleanly, all file handles close, the output is complete.
    // The binary stays resident in the page cache after the first call.
    const std::string& mopac_bin = RuntimeEnvironment::Mopac();
    if (mopac_bin.empty() || !fs::exists(mopac_bin)) {
        OperationLog::Error("MopacResult::Compute",
            "MOPAC binary not found (configured: " +
            (mopac_bin.empty() ? "<not set>" : mopac_bin) + ")");
        return nullptr;
    }

    OperationLog::Log(OperationLog::Level::Info, LogMopac,
        "MopacResult::Compute",
        "running PM7+MOZYME 1SCF atoms=" + std::to_string(natoms));

    std::string cmd = "OMP_STACKSIZE=2G OMP_NUM_THREADS=" +
        std::to_string(threads) + " " + mopac_bin + " " + mop_path +
        " > /dev/null 2>&1";
    int rc = std::system(cmd.c_str());

    // MOPAC writes .out alongside .mop (same stem, different extension)
    std::string out_path = mop_path.substr(0, mop_path.size() - 4) + ".out";

    std::error_code ec;

    if (rc != 0) {
        OperationLog::Error("MopacResult::Compute",
            "MOPAC returned error code " + std::to_string(rc));
        fs::remove(mop_path, ec);
        fs::remove(out_path, ec);
        return nullptr;
    }

    // Parse the .out file
    MopacParsed parsed = ParseMopacOutput(out_path, natoms);

    if (!parsed.success) {
        OperationLog::Error("MopacResult::Compute", parsed.error);
        // Leave temp files for debugging
        OperationLog::Error("MopacResult::Compute",
            "MOPAC output preserved at: " + out_path);
        return nullptr;
    }

    // Clean up temp files on success
    std::string stem = mop_path.substr(0, mop_path.size() - 4);
    fs::remove(mop_path, ec);
    fs::remove(out_path, ec);
    fs::remove(stem + ".arc", ec);
    fs::remove(stem + ".den", ec);
    fs::remove(stem + ".aux", ec);

    // Build result
    auto result = std::make_unique<MopacResult>();
    result->charges_ = std::move(parsed.charges);
    result->s_pop_ = std::move(parsed.s_pop);
    result->p_pop_ = std::move(parsed.p_pop);
    result->heat_of_formation_ = parsed.heat_of_formation;
    result->dipole_ = parsed.dipole;
    result->bond_orders_ = std::move(parsed.bond_orders);
    result->run_record_ = std::move(parsed.run_record);
    result->run_record_.net_charge = net_charge;
    result->run_record_.threads = threads;
    result->run_record_.mopac_binary = mopac_bin;

    // Build valencies from bond orders (sum of orders per atom = valency)
    result->valencies_.resize(natoms, 0.0);
    for (const auto& bo : result->bond_orders_) {
        result->valencies_[bo.atom_a] += bo.wiberg_order;
        // Bond orders are listed once per pair from the MOPAC output
        // (atom_i lists atom_j, but atom_j's line also lists atom_i).
        // We only add for atom_a here; atom_b's contribution comes from
        // the reversed entry elsewhere in the bond order list.
    }

    // Build O(1) bond order lookup map
    for (const auto& bo : result->bond_orders_) {
        uint64_t key = PairKey(bo.atom_a, bo.atom_b);
        // MOPAC lists each pair from both sides; keep the max reported order.
        auto it = result->bond_order_map_.find(key);
        if (it == result->bond_order_map_.end() || bo.wiberg_order > it->second) {
            result->bond_order_map_[key] = bo.wiberg_order;
        }
    }

    PopulateFullRecordFromParsedLegacy(result->run_record_,
                                       result->charges_,
                                       result->s_pop_,
                                       result->p_pop_,
                                       result->valencies_);
    FillElementLabelsFromProtein(protein, result->run_record_);

    // Build topology bridge: parallel to protein.Bonds()
    // Map each covalent bond's atom pair to its index for reverse lookup
    std::unordered_map<uint64_t, size_t> bond_pair_to_index;
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const auto& bond = protein.BondAt(bi);
        bond_pair_to_index[PairKey(bond.atom_index_a, bond.atom_index_b)] = bi;
    }

    result->topology_bond_orders_.resize(protein.BondCount(), 0.0);
    for (const auto& [key, order] : result->bond_order_map_) {
        auto it = bond_pair_to_index.find(key);
        if (it != bond_pair_to_index.end()) {
            result->topology_bond_orders_[it->second] = order;
        }
    }
    BuildTopologyBondRecords(protein, result->run_record_);

    // Store per-atom data on ConformationAtom
    for (size_t i = 0; i < natoms; ++i) {
        auto& ca = conf.MutableAtomAt(i);
        ca.mopac_charge = result->charges_[i];
        ca.mopac_s_pop = result->s_pop_[i];
        ca.mopac_p_pop = result->p_pop_[i];
        ca.mopac_valency = result->valencies_[i];
    }

    // Build per-atom bond neighbour lists on ConformationAtom
    // First, collect per atom
    std::vector<std::vector<MopacBondNeighbour>> per_atom(natoms);
    for (const auto& [key, order] : result->bond_order_map_) {
        size_t a = static_cast<size_t>(key >> 32);
        size_t b = static_cast<size_t>(key & 0xFFFFFFFF);

        // Look up topology bond index
        size_t topo_idx = SIZE_MAX;
        auto tit = bond_pair_to_index.find(key);
        if (tit != bond_pair_to_index.end()) topo_idx = tit->second;

        MopacBondNeighbour nb_a;
        nb_a.other_atom = b;
        nb_a.wiberg_order = order;
        nb_a.topology_bond_index = topo_idx;
        per_atom[a].push_back(nb_a);

        MopacBondNeighbour nb_b;
        nb_b.other_atom = a;
        nb_b.wiberg_order = order;
        nb_b.topology_bond_index = topo_idx;
        per_atom[b].push_back(nb_b);
    }

    // Sort descending by bond order and store
    for (size_t i = 0; i < natoms; ++i) {
        std::sort(per_atom[i].begin(), per_atom[i].end(),
            [](const MopacBondNeighbour& a, const MopacBondNeighbour& b) {
                return a.wiberg_order > b.wiberg_order;
            });
        conf.MutableAtomAt(i).mopac_bond_neighbours = std::move(per_atom[i]);
    }

    OperationLog::Log(OperationLog::Level::Info, LogMopac,
        "MopacResult::Compute",
        "heat=" + std::to_string(result->heat_of_formation_) +
        " kcal/mol, " + std::to_string(result->bond_order_map_.size()) +
        " bond orders");

    return result;
}


// ============================================================================
// Query methods
// ============================================================================

double MopacResult::ChargeAt(size_t i) const {
    return (i < charges_.size()) ? charges_[i] : 0.0;
}

double MopacResult::SPopAt(size_t i) const {
    return (i < s_pop_.size()) ? s_pop_[i] : 0.0;
}

double MopacResult::PPopAt(size_t i) const {
    return (i < p_pop_.size()) ? p_pop_[i] : 0.0;
}

double MopacResult::ValencyAt(size_t i) const {
    return (i < valencies_.size()) ? valencies_[i] : 0.0;
}

double MopacResult::BondOrder(size_t atom_a, size_t atom_b) const {
    auto it = bond_order_map_.find(PairKey(atom_a, atom_b));
    return (it != bond_order_map_.end()) ? it->second : 0.0;
}

double MopacResult::TopologyBondOrder(size_t bond_index) const {
    return (bond_index < topology_bond_orders_.size())
        ? topology_bond_orders_[bond_index] : 0.0;
}

const MopacAtomPopulation* MopacResult::AtomPopulationAt(size_t atom_index) const {
    return atom_index < run_record_.atom_populations.size()
        ? &run_record_.atom_populations[atom_index] : nullptr;
}

const MopacAOBasisFunction* MopacResult::AOBasisAt(size_t ao_index) const {
    return ao_index < run_record_.ao_basis.size()
        ? &run_record_.ao_basis[ao_index] : nullptr;
}

const MopacAuxRecord* MopacResult::AuxRecordByKey(const std::string& key) const {
    auto it = run_record_.aux_record_index.find(key);
    if (it == run_record_.aux_record_index.end()) return nullptr;
    return &run_record_.aux_records[it->second];
}


// ============================================================================
// WriteFeatures: current mopac_* NPY output schema.
//
//   mopac_charges.npy     — (N,)   Mulliken charges
//   mopac_scalars.npy     — (N, 4) [charge, s_pop, p_pop, valency]
//   mopac_bond_orders.npy — (B, 3) [atom_i, atom_j, bond_order]
//   mopac_global.npy      — (4,)   [heat_of_formation, dipole_x, dipole_y, dipole_z]
// ============================================================================

int MopacResult::WriteFeatures(const ProteinConformation& conf,
                                const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    int written = 0;

    // mopac_charges: (N,)
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).mopac_charge;
        NpyWriter::WriteFloat64(output_dir + "/mopac_charges.npy", data.data(), N);
        ++written;
    }

    // mopac_scalars: (N, 4) — [charge, s_pop, p_pop, valency]
    {
        std::vector<double> data(N * 4);
        for (size_t i = 0; i < N; ++i) {
            data[i*4 + 0] = conf.AtomAt(i).mopac_charge;
            data[i*4 + 1] = conf.AtomAt(i).mopac_s_pop;
            data[i*4 + 2] = conf.AtomAt(i).mopac_p_pop;
            data[i*4 + 3] = conf.AtomAt(i).mopac_valency;
        }
        NpyWriter::WriteFloat64(output_dir + "/mopac_scalars.npy", data.data(), N, 4);
        ++written;
    }

    // mopac_bond_orders: (B, 3) — [atom_i, atom_j, bond_order]
    // Only unique pairs (from bond_order_map_, not the raw list which has duplicates)
    {
        std::vector<double> data;
        data.reserve(bond_order_map_.size() * 3);
        for (const auto& [key, order] : bond_order_map_) {
            size_t a = static_cast<size_t>(key >> 32);
            size_t b = static_cast<size_t>(key & 0xFFFFFFFF);
            data.push_back(static_cast<double>(a));
            data.push_back(static_cast<double>(b));
            data.push_back(order);
        }
        size_t nbonds = bond_order_map_.size();
        if (nbonds > 0) {
            NpyWriter::WriteFloat64(output_dir + "/mopac_bond_orders.npy",
                                    data.data(), nbonds, 3);
        } else {
            // Write empty (0, 3) array
            NpyWriter::WriteFloat64(output_dir + "/mopac_bond_orders.npy",
                                    nullptr, 0, 3);
        }
        ++written;
    }

    // mopac_bond_neighbors: (M, 4) — [atom_i, atom_j, order, topology_bond_index]
    // Directed per-atom rows, sorted descending by order within each atom. A
    // missing covalent topology bond is represented as -1.
    {
        constexpr size_t C = 4;
        size_t rows = 0;
        for (size_t i = 0; i < N; ++i) {
            rows += conf.AtomAt(i).mopac_bond_neighbours.size();
        }

        std::vector<double> data(rows * C, 0.0);
        size_t row = 0;
        for (size_t i = 0; i < N; ++i) {
            for (const auto& nb : conf.AtomAt(i).mopac_bond_neighbours) {
                data[row*C + 0] = static_cast<double>(i);
                data[row*C + 1] = static_cast<double>(nb.other_atom);
                data[row*C + 2] = nb.wiberg_order;
                data[row*C + 3] = nb.topology_bond_index == SIZE_MAX
                    ? -1.0
                    : static_cast<double>(nb.topology_bond_index);
                ++row;
            }
        }

        NpyWriter::WriteFloat64(output_dir + "/mopac_bond_neighbors.npy",
                                data.empty() ? nullptr : data.data(), rows, C);
        ++written;
    }

    // mopac_global: (4,) — [heat_of_formation, dipole_x, dipole_y, dipole_z]
    {
        double data[4] = { heat_of_formation_, dipole_.x(), dipole_.y(), dipole_.z() };
        NpyWriter::WriteFloat64(output_dir + "/mopac_global.npy", data, 4);
        ++written;
    }

    // mopac_atom_populations: (N, 12)
    {
        constexpr size_t C = 12;
        std::vector<double> data(N * C, QuietNaN());
        for (size_t i = 0; i < N; ++i) {
            MopacAtomPopulation fallback;
            fallback.atom_index = i;
            fallback.net_charge = ChargeAt(i);
            fallback.s_population = SPopAt(i);
            fallback.p_population = PPopAt(i);
            fallback.project_valency = ValencyAt(i);
            const MopacAtomPopulation* p = AtomPopulationAt(i);
            const MopacAtomPopulation& row = p ? *p : fallback;
            data[i*C + 0] = row.net_charge;
            data[i*C + 1] = row.electron_density;
            data[i*C + 2] = row.s_population;
            data[i*C + 3] = row.p_population;
            data[i*C + 4] = row.d_population;
            data[i*C + 5] = row.f_population;
            data[i*C + 6] = row.dipole_x;
            data[i*C + 7] = row.dipole_y;
            data[i*C + 8] = row.dipole_z;
            data[i*C + 9] = row.dipole_total;
            data[i*C + 10] = row.mopac_valency;
            data[i*C + 11] = row.project_valency;
        }
        NpyWriter::WriteFloat64(output_dir + "/mopac_atom_populations.npy",
                                data.data(), N, C);
        ++written;
    }

    // mopac_atomic_orbital_populations: (K, 9)
    {
        constexpr size_t C = 9;
        const size_t K = run_record_.atomic_orbital_populations.size();
        std::vector<double> data(K * C, QuietNaN());
        for (size_t i = 0; i < K; ++i) {
            const auto& row = run_record_.atomic_orbital_populations[i];
            for (size_t c = 0; c < C; ++c)
                data[i*C + c] = row.populations[c];
        }
        NpyWriter::WriteFloat64(output_dir + "/mopac_atomic_orbital_populations.npy",
                                data.empty() ? nullptr : data.data(), K, C);
        ++written;
    }

    // mopac_bond_valencies: (N,) — MOPAC diagonal valencies, not recomputed.
    {
        std::vector<double> data(N, QuietNaN());
        for (size_t i = 0; i < N; ++i) {
            if (const MopacAtomPopulation* p = AtomPopulationAt(i))
                data[i] = p->mopac_valency;
        }
        NpyWriter::WriteFloat64(output_dir + "/mopac_bond_valencies.npy",
                                data.data(), N);
        ++written;
    }

    // mopac_bond_orders_unique: (U, 8)
    {
        constexpr size_t C = 8;
        const size_t U = run_record_.unique_bond_orders.size();
        std::vector<double> data(U * C, QuietNaN());
        for (size_t i = 0; i < U; ++i) {
            const auto& u = run_record_.unique_bond_orders[i];
            data[i*C + 0] = static_cast<double>(u.atom_a);
            data[i*C + 1] = static_cast<double>(u.atom_b);
            data[i*C + 2] = u.max_order;
            data[i*C + 3] = u.mean_order;
            data[i*C + 4] = static_cast<double>(u.printed_entry_indices.size());
            data[i*C + 5] = u.printed_entry_indices.empty()
                ? QuietNaN() : static_cast<double>(u.printed_entry_indices[0]);
            data[i*C + 6] = u.printed_entry_indices.size() < 2
                ? QuietNaN() : static_cast<double>(u.printed_entry_indices[1]);
            data[i*C + 7] = u.topology_bond_index == SIZE_MAX
                ? QuietNaN() : static_cast<double>(u.topology_bond_index);
        }
        NpyWriter::WriteFloat64(output_dir + "/mopac_bond_orders_unique.npy",
                                data.empty() ? nullptr : data.data(), U, C);
        ++written;
    }

    // mopac_topology_bond_orders_full: (B_topo, 8)
    {
        constexpr size_t C = 8;
        const size_t B = run_record_.topology_bond_records.size();
        std::vector<double> data(B * C, QuietNaN());
        for (size_t i = 0; i < B; ++i) {
            const auto& row = run_record_.topology_bond_records[i];
            data[i*C + 0] = static_cast<double>(row.bond_index);
            data[i*C + 1] = static_cast<double>(row.atom_a);
            data[i*C + 2] = static_cast<double>(row.atom_b);
            data[i*C + 3] = row.order;
            data[i*C + 4] = row.present ? 1.0 : 0.0;
            data[i*C + 5] = row.unique_pair_index == SIZE_MAX
                ? QuietNaN() : static_cast<double>(row.unique_pair_index);
            data[i*C + 6] = row.present ? 0.0 : 1.0;
            data[i*C + 7] = static_cast<double>(row.printed_entry_count);
        }
        NpyWriter::WriteFloat64(output_dir + "/mopac_topology_bond_orders_full.npy",
                                data.empty() ? nullptr : data.data(), B, C);
        ++written;
    }

    return written;
}

}  // namespace nmr

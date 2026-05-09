#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

// Compile the generated implementation into this test translation unit so the
// internal constexpr arrays are visible to the audit harness. This target does
// not link libnmr_shielding, avoiding duplicate LookupBy/LookupCap symbols.
#include "generated/LegacyAmberSemanticTables.cpp"

namespace {

using nmr::AminoAcid;
using nmr::AtomMechanicalIdentity;
using nmr::AtomSemanticTable;
using nmr::BackboneRole;
using nmr::BranchAddress;
using nmr::DiastereotopicIndex;
using nmr::Element;
using nmr::Locant;
using nmr::PlanarGroupKind;
using nmr::PlanarStereo;
using nmr::PolarHKind;
using nmr::ProchiralStereo;
using nmr::PseudoatomKind;
using nmr::RingPositionLabel;
using nmr::RingSystemKind;
using nmr::TerminalState;

namespace gen = nmr::topology_generated;

struct TableRef {
    const char* name;
    const char* label;
    const AtomSemanticTable* rows;
    std::size_t size;
    AminoAcid residue;
    std::uint8_t variant_idx;
    TerminalState cap_state;
    bool is_cap;
};

template <std::size_t N>
TableRef ResidueTable(const char* name,
                      const char* label,
                      AminoAcid residue,
                      std::uint8_t variant_idx,
                      const std::array<AtomSemanticTable, N>& rows) {
    return TableRef{name, label, rows.data(), rows.size(), residue,
                    variant_idx, TerminalState::Internal, false};
}

template <std::size_t N>
TableRef CapTable(const char* name,
                  const char* label,
                  TerminalState state,
                  const std::array<AtomSemanticTable, N>& rows) {
    return TableRef{name, label, rows.data(), rows.size(), AminoAcid::Unknown,
                    0, state, true};
}

std::vector<TableRef> AllTables() {
    return {
        ResidueTable("kAlaAtoms", "ALA", AminoAcid::ALA, gen::kBaseVariantIdx, gen::kAlaAtoms),
        ResidueTable("kArgAtoms", "ARG", AminoAcid::ARG, gen::kBaseVariantIdx, gen::kArgAtoms),
        ResidueTable("kAsnAtoms", "ASN", AminoAcid::ASN, gen::kBaseVariantIdx, gen::kAsnAtoms),
        ResidueTable("kAspAtoms", "ASP", AminoAcid::ASP, gen::kBaseVariantIdx, gen::kAspAtoms),
        ResidueTable("kCysAtoms", "CYS", AminoAcid::CYS, gen::kBaseVariantIdx, gen::kCysAtoms),
        ResidueTable("kGlnAtoms", "GLN", AminoAcid::GLN, gen::kBaseVariantIdx, gen::kGlnAtoms),
        ResidueTable("kGluAtoms", "GLU", AminoAcid::GLU, gen::kBaseVariantIdx, gen::kGluAtoms),
        ResidueTable("kGlyAtoms", "GLY", AminoAcid::GLY, gen::kBaseVariantIdx, gen::kGlyAtoms),
        ResidueTable("kHisAtoms", "HIS", AminoAcid::HIS, gen::kBaseVariantIdx, gen::kHisAtoms),
        ResidueTable("kIleAtoms", "ILE", AminoAcid::ILE, gen::kBaseVariantIdx, gen::kIleAtoms),
        ResidueTable("kLeuAtoms", "LEU", AminoAcid::LEU, gen::kBaseVariantIdx, gen::kLeuAtoms),
        ResidueTable("kLysAtoms", "LYS", AminoAcid::LYS, gen::kBaseVariantIdx, gen::kLysAtoms),
        ResidueTable("kMetAtoms", "MET", AminoAcid::MET, gen::kBaseVariantIdx, gen::kMetAtoms),
        ResidueTable("kPheAtoms", "PHE", AminoAcid::PHE, gen::kBaseVariantIdx, gen::kPheAtoms),
        ResidueTable("kProAtoms", "PRO", AminoAcid::PRO, gen::kBaseVariantIdx, gen::kProAtoms),
        ResidueTable("kSerAtoms", "SER", AminoAcid::SER, gen::kBaseVariantIdx, gen::kSerAtoms),
        ResidueTable("kThrAtoms", "THR", AminoAcid::THR, gen::kBaseVariantIdx, gen::kThrAtoms),
        ResidueTable("kTrpAtoms", "TRP", AminoAcid::TRP, gen::kBaseVariantIdx, gen::kTrpAtoms),
        ResidueTable("kTyrAtoms", "TYR", AminoAcid::TYR, gen::kBaseVariantIdx, gen::kTyrAtoms),
        ResidueTable("kValAtoms", "VAL", AminoAcid::VAL, gen::kBaseVariantIdx, gen::kValAtoms),
        ResidueTable("kHisAtoms_HID", "HIS/HID", AminoAcid::HIS, 0, gen::kHisAtoms_HID),
        ResidueTable("kHisAtoms_HIE", "HIS/HIE", AminoAcid::HIS, 1, gen::kHisAtoms_HIE),
        ResidueTable("kHisAtoms_HIP", "HIS/HIP", AminoAcid::HIS, 2, gen::kHisAtoms_HIP),
        ResidueTable("kAspAtoms_ASH", "ASP/ASH", AminoAcid::ASP, 0, gen::kAspAtoms_ASH),
        ResidueTable("kGluAtoms_GLH", "GLU/GLH", AminoAcid::GLU, 0, gen::kGluAtoms_GLH),
        ResidueTable("kCysAtoms_CYX", "CYS/CYX", AminoAcid::CYS, 0, gen::kCysAtoms_CYX),
        ResidueTable("kCysAtoms_CYM", "CYS/CYM", AminoAcid::CYS, 1, gen::kCysAtoms_CYM),
        ResidueTable("kLysAtoms_LYN", "LYS/LYN", AminoAcid::LYS, 0, gen::kLysAtoms_LYN),
        ResidueTable("kArgAtoms_ARN", "ARG/ARN", AminoAcid::ARG, 0, gen::kArgAtoms_ARN),
        ResidueTable("kTyrAtoms_TYM", "TYR/TYM", AminoAcid::TYR, 0, gen::kTyrAtoms_TYM),
        CapTable("kCapNtermCharged", "NTERM_CHARGED", TerminalState::NtermCharged, gen::kCapNtermCharged),
        CapTable("kCapNtermNeutral", "NTERM_NEUTRAL", TerminalState::NtermNeutral, gen::kCapNtermNeutral),
        CapTable("kCapCtermDeprotonated", "CTERM_DEPROTONATED", TerminalState::CtermDeprotonated, gen::kCapCtermDeprotonated),
        CapTable("kCapCtermProtonated", "CTERM_PROTONATED", TerminalState::CtermProtonated, gen::kCapCtermProtonated),
    };
}

std::string SourceRoot() {
#ifdef NMR_SOURCE_DIR
    return NMR_SOURCE_DIR;
#else
    return ".";
#endif
}

std::string Trim(std::string s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) {
        s.erase(s.begin());
    }
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) {
        s.pop_back();
    }
    return s;
}

std::vector<std::string> SplitWhitespace(const std::string& s) {
    std::istringstream iss(s);
    std::vector<std::string> fields;
    for (std::string field; iss >> field;) {
        fields.push_back(field);
    }
    return fields;
}

std::string Slurp(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        ADD_FAILURE() << "Cannot open " << path;
        return {};
    }
    std::ostringstream oss;
    oss << in.rdbuf();
    return oss.str();
}

std::map<std::string, std::vector<std::string>> ParseGeneratedAtomNames() {
    const std::string path = SourceRoot() + "/src/generated/LegacyAmberSemanticTables.cpp";
    std::ifstream in(path);
    if (!in.is_open()) {
        ADD_FAILURE() << "Cannot open generated table source: " << path;
        return {};
    }

    const std::regex table_re(
        R"(^constexpr std::array<AtomSemanticTable, [0-9]+> (k[A-Za-z0-9_]+) = \{\{)");
    const std::regex atom_re(R"(//[[:space:]]*[0-9]+:[[:space:]]*([A-Z0-9]+)[[:space:]]*$)");

    std::map<std::string, std::vector<std::string>> out;
    std::string current;
    for (std::string line; std::getline(in, line);) {
        std::smatch m;
        if (std::regex_search(line, m, table_re)) {
            current = m[1].str();
            out[current] = {};
            continue;
        }
        if (current.empty()) continue;
        if (std::regex_search(line, m, atom_re)) {
            out[current].push_back(m[1].str());
        }
        if (line.find("}};") != std::string::npos) {
            current.clear();
        }
    }
    return out;
}

std::map<std::string, std::set<std::string>> ParseStandardCcdElementSymbols() {
    const std::set<std::string> standard_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL",
    };
    const std::string path = SourceRoot() + "/data/ccd/components.cif";
    std::ifstream in(path);
    if (!in.is_open()) {
        ADD_FAILURE() << "Cannot open CCD component file: " << path;
        return {};
    }

    std::map<std::string, std::set<std::string>> symbols;
    std::string current;
    bool reading_atom_headers = false;
    bool reading_atom_rows = false;
    std::vector<std::string> headers;
    std::size_t type_symbol_idx = static_cast<std::size_t>(-1);

    for (std::string raw; std::getline(in, raw);) {
        const std::string line = Trim(raw);
        if (line.rfind("data_", 0) == 0) {
            const std::string comp = line.substr(5);
            current = standard_residues.count(comp) ? comp : "";
            reading_atom_headers = false;
            reading_atom_rows = false;
            headers.clear();
            type_symbol_idx = static_cast<std::size_t>(-1);
            continue;
        }
        if (current.empty()) continue;

        if (line == "loop_") {
            reading_atom_headers = true;
            reading_atom_rows = false;
            headers.clear();
            type_symbol_idx = static_cast<std::size_t>(-1);
            continue;
        }

        bool process_as_row = false;
        if (reading_atom_headers) {
            if (line.rfind("_chem_comp_atom.", 0) == 0) {
                headers.push_back(line);
                if (line == "_chem_comp_atom.type_symbol") {
                    type_symbol_idx = headers.size() - 1;
                }
                continue;
            }

            if (!headers.empty() &&
                type_symbol_idx != static_cast<std::size_t>(-1)) {
                reading_atom_headers = false;
                reading_atom_rows = true;
                process_as_row = true;
            } else {
                reading_atom_headers = false;
            }
        }

        if (reading_atom_rows || process_as_row) {
            if (line.empty() || line[0] == '#' || line == "loop_" ||
                line[0] == '_' || line.rfind("data_", 0) == 0) {
                reading_atom_rows = false;
                continue;
            }
            const auto fields = SplitWhitespace(line);
            if (fields.size() > type_symbol_idx && fields.front() == current) {
                symbols[current].insert(fields[type_symbol_idx]);
            }
        }
    }
    return symbols;
}

const TableRef* FindTable(const std::vector<TableRef>& tables,
                          const std::string& name) {
    for (const auto& table : tables) {
        if (name == table.name) return &table;
    }
    return nullptr;
}

const std::vector<std::string>& NamesFor(
    const std::map<std::string, std::vector<std::string>>& parsed,
    const TableRef& table) {
    static const std::vector<std::string> empty;
    const auto it = parsed.find(table.name);
    if (it == parsed.end()) {
        ADD_FAILURE() << "Generated source did not expose atom comments for "
                      << table.name;
        return empty;
    }
    return it->second;
}

std::size_t IndexOfAtom(const std::vector<std::string>& names,
                        const std::string& atom) {
    const auto it = std::find(names.begin(), names.end(), atom);
    if (it == names.end()) return names.size();
    return static_cast<std::size_t>(std::distance(names.begin(), it));
}

const AtomSemanticTable* FindAtom(
    const std::map<std::string, std::vector<std::string>>& parsed,
    const std::vector<TableRef>& tables,
    const std::string& table_name,
    const std::string& atom) {
    const TableRef* table = FindTable(tables, table_name);
    if (table == nullptr) {
        ADD_FAILURE() << "Unknown test table " << table_name;
        return nullptr;
    }
    const auto& names = NamesFor(parsed, *table);
    const std::size_t idx = IndexOfAtom(names, atom);
    if (idx >= table->size) {
        ADD_FAILURE() << "Atom " << atom << " not found in " << table_name;
        return nullptr;
    }
    return &table->rows[idx];
}

void ExpectHasAtom(const std::map<std::string, std::vector<std::string>>& parsed,
                   const std::string& table_name,
                   const std::string& atom) {
    const auto it = parsed.find(table_name);
    ASSERT_NE(it, parsed.end()) << table_name;
    EXPECT_NE(IndexOfAtom(it->second, atom), it->second.size())
        << table_name << " should contain " << atom;
}

void ExpectNoAtom(const std::map<std::string, std::vector<std::string>>& parsed,
                  const std::string& table_name,
                  const std::string& atom) {
    const auto it = parsed.find(table_name);
    ASSERT_NE(it, parsed.end()) << table_name;
    EXPECT_EQ(IndexOfAtom(it->second, atom), it->second.size())
        << table_name << " should not contain " << atom;
}

AtomMechanicalIdentity IdentityOf(const AtomSemanticTable& row) {
    return AtomMechanicalIdentity{row.element, row.locant, row.branch,
                                  row.di_index, row.backbone_role};
}

struct IdentityKey {
    Element element;
    Locant locant;
    BranchAddress branch;
    DiastereotopicIndex di_index;
    BackboneRole backbone_role;

    explicit IdentityKey(const AtomSemanticTable& row)
        : element(row.element),
          locant(row.locant),
          branch(row.branch),
          di_index(row.di_index),
          backbone_role(row.backbone_role) {}

    bool operator<(const IdentityKey& other) const {
        return std::tie(element, locant, branch.outer, branch.inner, di_index,
                        backbone_role) <
               std::tie(other.element, other.locant, other.branch.outer,
                        other.branch.inner, other.di_index,
                        other.backbone_role);
    }
};

bool SamePseudoatom(const nmr::PseudoatomMembership& a,
                    const nmr::PseudoatomMembership& b) {
    return a.kind == b.kind && a.locant == b.locant &&
           a.branch == b.branch && a.in_super_group == b.in_super_group;
}

bool SameRingMembership(const nmr::RingMembership& a,
                        const nmr::RingMembership& b) {
    return a.ring == b.ring && a.position == b.position &&
           a.ring_size == b.ring_size && a.aromatic == b.aromatic &&
           a.planar == b.planar && a.n_heteroatoms == b.n_heteroatoms;
}

bool SameRingPosition(const nmr::RingPosition& a,
                      const nmr::RingPosition& b) {
    return SameRingMembership(a.primary, b.primary) &&
           SameRingMembership(a.secondary, b.secondary);
}

bool SameSemanticRow(const AtomSemanticTable& a, const AtomSemanticTable& b) {
    return a.element == b.element && a.locant == b.locant &&
           a.branch == b.branch && a.di_index == b.di_index &&
           a.backbone_role == b.backbone_role &&
           a.prochiral == b.prochiral &&
           a.planar_group == b.planar_group &&
           a.planar_stereo == b.planar_stereo &&
           SamePseudoatom(a.pseudoatom, b.pseudoatom) &&
           a.polar_h == b.polar_h &&
           SameRingPosition(a.ring_position, b.ring_position) &&
           a.aromatic == b.aromatic &&
           a.formal_charge == b.formal_charge &&
           a.is_exchangeable == b.is_exchangeable &&
           a.equivalence_class == b.equivalence_class;
}

std::string JoinNames(const std::vector<std::string>& names,
                      const std::vector<std::size_t>& indices) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < indices.size(); ++i) {
        if (i != 0) oss << ",";
        oss << names[indices[i]];
    }
    return oss.str();
}

bool IsAllowedEquivalentHydrogenCollision(const TableRef& table,
                                          const std::vector<std::string>& names,
                                          const std::vector<std::size_t>& indices) {
    bool all_h = true;
    bool all_methyl = true;
    std::set<std::string> atom_names;
    for (const std::size_t idx : indices) {
        const auto& row = table.rows[idx];
        all_h = all_h && row.element == Element::H;
        all_methyl = all_methyl && row.pseudoatom.kind == PseudoatomKind::M;
        atom_names.insert(names[idx]);
    }
    if (all_h && all_methyl) return true;
    if (std::string(table.name) == "kCapNtermCharged") {
        return atom_names == std::set<std::string>{"H1", "H2", "H3"};
    }
    if (std::string(table.name) == "kCapNtermNeutral") {
        return atom_names == std::set<std::string>{"H1", "H2"};
    }
    return false;
}

TEST(TopologySemanticTables, GeneratedSourceNamesStayAlignedWithCompiledArrays) {
    const auto parsed = ParseGeneratedAtomNames();
    for (const auto& table : AllTables()) {
        SCOPED_TRACE(table.label);
        const auto it = parsed.find(table.name);
        ASSERT_NE(it, parsed.end()) << table.name;
        EXPECT_EQ(table.size, it->second.size())
            << table.name << " compiled row count and generated atom-name "
            << "comments diverged";
    }
}

TEST(TopologySemanticTables, AllRowsHavePopulatedElementsAndPolarConsistency) {
    for (const auto& table : AllTables()) {
        SCOPED_TRACE(table.label);
        for (std::size_t i = 0; i < table.size; ++i) {
            const auto& row = table.rows[i];
            EXPECT_NE(Element::Unknown, row.element) << table.name << "[" << i << "]";
            if (row.polar_h == PolarHKind::NotPolar) {
                EXPECT_FALSE(row.is_exchangeable) << table.name << "[" << i << "]";
            } else {
                EXPECT_EQ(Element::H, row.element) << table.name << "[" << i << "]";
                EXPECT_TRUE(row.is_exchangeable) << table.name << "[" << i << "]";
            }
        }
    }
}

TEST(TopologySemanticTables, CcdStandardResidueElementsStayInsideRuntimeScope) {
    const auto symbols = ParseStandardCcdElementSymbols();
    const std::set<std::string> expected_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL",
    };
    const std::set<std::string> runtime_scope = {"H", "C", "N", "O", "S"};

    EXPECT_EQ(expected_residues.size(), symbols.size())
        << "The CCD scope parser should find all 20 standard residues";
    std::set<std::string> observed_symbols;
    for (const auto& residue : expected_residues) {
        const auto it = symbols.find(residue);
        ASSERT_NE(it, symbols.end()) << residue;
        EXPECT_FALSE(it->second.empty()) << residue;
        for (const auto& symbol : it->second) {
            observed_symbols.insert(symbol);
            EXPECT_TRUE(runtime_scope.count(symbol))
                << residue << " contains CCD type_symbol=" << symbol
                << "; runtime Element currently encodes only H/C/N/O/S";
        }
    }
    EXPECT_EQ(runtime_scope, observed_symbols)
        << "The standard-20 CCD grounding should exercise the whole "
        << "protein Element runtime scope";
}

TEST(TopologySemanticTables, GeneratorLoggedNoUnknownCcdElementsForCurrentTables) {
    const std::string log =
        Slurp(SourceRoot() + "/src/generated/LegacyAmberSemanticTables.log.txt");
    ASSERT_FALSE(log.empty());
    EXPECT_EQ(std::string::npos, log.find("warning_unknown_element"))
        << "The generator skipped at least one CCD atom because its "
        << "type_symbol was not mapped at the RDKit build boundary";
}

TEST(TopologySemanticTables, MechanicalIdentitiesAreUniqueExceptEquivalentHydrogenSets) {
    const auto parsed = ParseGeneratedAtomNames();
    for (const auto& table : AllTables()) {
        SCOPED_TRACE(table.label);
        const auto& names = NamesFor(parsed, table);
        ASSERT_EQ(table.size, names.size());

        std::map<IdentityKey, std::vector<std::size_t>> groups;
        for (std::size_t i = 0; i < table.size; ++i) {
            groups[IdentityKey(table.rows[i])].push_back(i);
        }

        for (const auto& item : groups) {
            const auto& indices = item.second;
            if (indices.size() <= 1) continue;

            EXPECT_TRUE(IsAllowedEquivalentHydrogenCollision(table, names, indices))
                << table.name << " has an unwhitelisted mechanical-identity "
                << "collision among " << JoinNames(names, indices);

            for (std::size_t i = 1; i < indices.size(); ++i) {
                EXPECT_TRUE(SameSemanticRow(table.rows[indices[0]], table.rows[indices[i]]))
                    << table.name << " equivalent-H collision rows differ for "
                    << JoinNames(names, indices);
            }
        }
    }
}

TEST(TopologySemanticTables, LookupByDispatchesEveryResidueAndVariantRow) {
    for (const auto& table : AllTables()) {
        if (table.is_cap) continue;
        SCOPED_TRACE(table.label);
        for (std::size_t i = 0; i < table.size; ++i) {
            const auto hit = gen::LookupBy(table.residue, table.variant_idx,
                                           IdentityOf(table.rows[i]));
            ASSERT_NE(nullptr, hit) << table.name << "[" << i << "]";
            EXPECT_EQ(IdentityOf(table.rows[i]), IdentityOf(*hit))
                << table.name << "[" << i << "]";
        }

        AtomMechanicalIdentity missing = IdentityOf(table.rows[0]);
        missing.element = Element::Unknown;
        EXPECT_EQ(nullptr, gen::LookupBy(table.residue, table.variant_idx, missing))
            << table.name << " should fail lookup on an impossible identity";
    }
}

TEST(TopologySemanticTables, LookupByFailsFastOnInvalidResidueOrVariantIndex) {
    const auto id = IdentityOf(gen::kAlaAtoms[0]);
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::Unknown, gen::kBaseVariantIdx, id));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::ALA, 0, id))
        << "ALA has no variant 0; callers must pass kBaseVariantIdx for chain";
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::ALA, 254, id));

    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 0, IdentityOf(gen::kHisAtoms_HID[0])));
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 1, IdentityOf(gen::kHisAtoms_HIE[0])));
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 2, IdentityOf(gen::kHisAtoms_HIP[0])));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::HIS, 3, IdentityOf(gen::kHisAtoms[0])));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::HIS, 254, IdentityOf(gen::kHisAtoms[0])));
}

TEST(TopologySemanticTables, LookupCapCoversAllTerminalStatesAndMissesCleanly) {
    for (const auto& table : AllTables()) {
        if (!table.is_cap) continue;
        SCOPED_TRACE(table.label);
        for (std::size_t i = 0; i < table.size; ++i) {
            const auto hit = gen::LookupCap(table.cap_state, IdentityOf(table.rows[i]));
            ASSERT_NE(nullptr, hit) << table.name << "[" << i << "]";
            EXPECT_EQ(IdentityOf(table.rows[i]), IdentityOf(*hit))
                << table.name << "[" << i << "]";
            EXPECT_EQ(nullptr, gen::LookupCap(TerminalState::Internal,
                                              IdentityOf(table.rows[i])));
        }

        AtomMechanicalIdentity missing = IdentityOf(table.rows[0]);
        missing.element = Element::Unknown;
        EXPECT_EQ(nullptr, gen::LookupCap(table.cap_state, missing));
    }
}

TEST(TopologySemanticTables, CapTablesMatchSectionFourChemistry) {
    const auto parsed = ParseGeneratedAtomNames();
    const auto tables = AllTables();

    const auto* nterm_charged_n = FindAtom(parsed, tables, "kCapNtermCharged", "N");
    ASSERT_NE(nullptr, nterm_charged_n);
    EXPECT_EQ(1, nterm_charged_n->formal_charge);
    for (const char* atom : {"H1", "H2", "H3"}) {
        const auto* h = FindAtom(parsed, tables, "kCapNtermCharged", atom);
        ASSERT_NE(nullptr, h) << atom;
        EXPECT_EQ(Element::H, h->element) << atom;
        EXPECT_EQ(PolarHKind::AmmoniumNH, h->polar_h) << atom;
        EXPECT_EQ(0, h->formal_charge) << atom;
        EXPECT_TRUE(h->is_exchangeable) << atom;
    }

    const auto* nterm_neutral_n = FindAtom(parsed, tables, "kCapNtermNeutral", "N");
    ASSERT_NE(nullptr, nterm_neutral_n);
    EXPECT_EQ(0, nterm_neutral_n->formal_charge);
    for (const char* atom : {"H1", "H2"}) {
        const auto* h = FindAtom(parsed, tables, "kCapNtermNeutral", atom);
        ASSERT_NE(nullptr, h) << atom;
        EXPECT_EQ(Element::H, h->element) << atom;
        EXPECT_EQ(PolarHKind::AmineNH, h->polar_h) << atom;
        EXPECT_EQ(0, h->formal_charge) << atom;
        EXPECT_TRUE(h->is_exchangeable) << atom;
    }
    ExpectNoAtom(parsed, "kCapNtermNeutral", "H3");

    const auto* cterm_deprot_oxt =
        FindAtom(parsed, tables, "kCapCtermDeprotonated", "OXT");
    ASSERT_NE(nullptr, cterm_deprot_oxt);
    EXPECT_EQ(Element::O, cterm_deprot_oxt->element);
    EXPECT_EQ(PlanarGroupKind::Carboxylate, cterm_deprot_oxt->planar_group);
    EXPECT_EQ(-1, cterm_deprot_oxt->formal_charge);
    EXPECT_FALSE(cterm_deprot_oxt->is_exchangeable);

    const auto* cterm_prot_oxt =
        FindAtom(parsed, tables, "kCapCtermProtonated", "OXT");
    ASSERT_NE(nullptr, cterm_prot_oxt);
    EXPECT_EQ(Element::O, cterm_prot_oxt->element);
    EXPECT_EQ(PlanarGroupKind::Carboxylate, cterm_prot_oxt->planar_group);
    EXPECT_EQ(0, cterm_prot_oxt->formal_charge);

    const auto* hxt = FindAtom(parsed, tables, "kCapCtermProtonated", "HXT");
    ASSERT_NE(nullptr, hxt);
    EXPECT_EQ(Element::H, hxt->element);
    EXPECT_EQ(PolarHKind::CarboxylOH, hxt->polar_h);
    EXPECT_EQ(0, hxt->formal_charge);
    EXPECT_TRUE(hxt->is_exchangeable);
}

TEST(TopologySemanticTables, BackboneRoleConventionsStayPinned) {
    const auto parsed = ParseGeneratedAtomNames();
    const auto tables = AllTables();

    const auto* pro_n = FindAtom(parsed, tables, "kProAtoms", "N");
    ASSERT_NE(nullptr, pro_n);
    EXPECT_EQ(BackboneRole::Nitrogen, pro_n->backbone_role);
    EXPECT_EQ(Locant::None, pro_n->locant);
    EXPECT_EQ(RingSystemKind::Pyrrolidine_Pro, pro_n->ring_position.primary.ring)
        << "PRO N remains a backbone N even though it is also in the "
        << "pyrrolidine ring; ring membership is orthogonal to BackboneRole";

    const auto* gly_ha2 = FindAtom(parsed, tables, "kGlyAtoms", "HA2");
    const auto* gly_ha3 = FindAtom(parsed, tables, "kGlyAtoms", "HA3");
    ASSERT_NE(nullptr, gly_ha2);
    ASSERT_NE(nullptr, gly_ha3);
    EXPECT_EQ(Locant::Alpha, gly_ha2->locant);
    EXPECT_EQ(Locant::Alpha, gly_ha3->locant);
    EXPECT_EQ(DiastereotopicIndex::Position2, gly_ha2->di_index);
    EXPECT_EQ(DiastereotopicIndex::Position3, gly_ha3->di_index);
    EXPECT_EQ(BackboneRole::None, gly_ha2->backbone_role);
    EXPECT_EQ(BackboneRole::None, gly_ha3->backbone_role);

    const auto* ala_ha = FindAtom(parsed, tables, "kAlaAtoms", "HA");
    ASSERT_NE(nullptr, ala_ha);
    EXPECT_EQ(Locant::Alpha, ala_ha->locant);
    EXPECT_EQ(BackboneRole::AlphaHydrogen, ala_ha->backbone_role);
}

TEST(TopologySemanticTables, StandardResidueTablesExcludeTerminalAndVariantOnlyAtoms) {
    const auto parsed = ParseGeneratedAtomNames();
    for (const auto& table : AllTables()) {
        if (table.is_cap || table.variant_idx != gen::kBaseVariantIdx) continue;
        SCOPED_TRACE(table.label);
        for (const char* cap_atom : {"H1", "H2", "H3", "OXT", "HXT"}) {
            ExpectNoAtom(parsed, table.name, cap_atom);
        }
    }

    ExpectNoAtom(parsed, "kAspAtoms", "HD2");
    ExpectNoAtom(parsed, "kGluAtoms", "HE2");
    ExpectNoAtom(parsed, "kHisAtoms", "HD1");
    ExpectNoAtom(parsed, "kProAtoms", "H");
}

TEST(TopologySemanticTables, VariantInventoriesAndChemistryStayPinned) {
    const auto parsed = ParseGeneratedAtomNames();
    const auto tables = AllTables();

    ExpectHasAtom(parsed, "kAspAtoms_ASH", "HD2");
    ExpectNoAtom(parsed, "kAspAtoms", "HD2");
    const auto* asp_od2 = FindAtom(parsed, tables, "kAspAtoms", "OD2");
    const auto* ash_od2 = FindAtom(parsed, tables, "kAspAtoms_ASH", "OD2");
    const auto* ash_hd2 = FindAtom(parsed, tables, "kAspAtoms_ASH", "HD2");
    ASSERT_NE(nullptr, asp_od2);
    ASSERT_NE(nullptr, ash_od2);
    ASSERT_NE(nullptr, ash_hd2);
    EXPECT_EQ(-1, asp_od2->formal_charge);
    EXPECT_EQ(0, ash_od2->formal_charge);
    EXPECT_EQ(PolarHKind::CarboxylOH, ash_hd2->polar_h);

    ExpectHasAtom(parsed, "kGluAtoms_GLH", "HE2");
    ExpectNoAtom(parsed, "kGluAtoms", "HE2");
    const auto* glu_oe2 = FindAtom(parsed, tables, "kGluAtoms", "OE2");
    const auto* glh_oe2 = FindAtom(parsed, tables, "kGluAtoms_GLH", "OE2");
    const auto* glh_he2 = FindAtom(parsed, tables, "kGluAtoms_GLH", "HE2");
    ASSERT_NE(nullptr, glu_oe2);
    ASSERT_NE(nullptr, glh_oe2);
    ASSERT_NE(nullptr, glh_he2);
    EXPECT_EQ(-1, glu_oe2->formal_charge);
    EXPECT_EQ(0, glh_oe2->formal_charge);
    EXPECT_EQ(PolarHKind::CarboxylOH, glh_he2->polar_h);

    ExpectHasAtom(parsed, "kHisAtoms_HID", "HD1");
    ExpectNoAtom(parsed, "kHisAtoms_HID", "HE2");
    ExpectNoAtom(parsed, "kHisAtoms_HIE", "HD1");
    ExpectHasAtom(parsed, "kHisAtoms_HIE", "HE2");
    ExpectHasAtom(parsed, "kHisAtoms_HIP", "HD1");
    ExpectHasAtom(parsed, "kHisAtoms_HIP", "HE2");
    const auto* hid_hd1 = FindAtom(parsed, tables, "kHisAtoms_HID", "HD1");
    const auto* hie_he2 = FindAtom(parsed, tables, "kHisAtoms_HIE", "HE2");
    const auto* hip_hd1 = FindAtom(parsed, tables, "kHisAtoms_HIP", "HD1");
    const auto* hip_he2 = FindAtom(parsed, tables, "kHisAtoms_HIP", "HE2");
    ASSERT_NE(nullptr, hid_hd1);
    ASSERT_NE(nullptr, hie_he2);
    ASSERT_NE(nullptr, hip_hd1);
    ASSERT_NE(nullptr, hip_he2);
    EXPECT_EQ(PolarHKind::ImidazoleNH, hid_hd1->polar_h);
    EXPECT_EQ(PolarHKind::ImidazoleNH, hie_he2->polar_h);
    EXPECT_EQ(PolarHKind::ImidazoleNH, hip_hd1->polar_h);
    EXPECT_EQ(PolarHKind::ImidazoleNH, hip_he2->polar_h);

    ExpectHasAtom(parsed, "kCysAtoms", "HG");
    ExpectNoAtom(parsed, "kCysAtoms_CYX", "HG");
    ExpectNoAtom(parsed, "kCysAtoms_CYM", "HG");

    ExpectHasAtom(parsed, "kLysAtoms", "HZ1");
    ExpectHasAtom(parsed, "kLysAtoms_LYN", "HZ2");
    ExpectHasAtom(parsed, "kLysAtoms_LYN", "HZ3");
    ExpectNoAtom(parsed, "kLysAtoms_LYN", "HZ1");
    const auto* lyn_hz2 = FindAtom(parsed, tables, "kLysAtoms_LYN", "HZ2");
    const auto* lyn_hz3 = FindAtom(parsed, tables, "kLysAtoms_LYN", "HZ3");
    ASSERT_NE(nullptr, lyn_hz2);
    ASSERT_NE(nullptr, lyn_hz3);
    EXPECT_EQ(PolarHKind::AmineNH, lyn_hz2->polar_h);
    EXPECT_EQ(PolarHKind::AmineNH, lyn_hz3->polar_h);

    ExpectHasAtom(parsed, "kArgAtoms", "HE");
    ExpectNoAtom(parsed, "kArgAtoms_ARN", "HE");

    ExpectHasAtom(parsed, "kTyrAtoms", "HH");
    ExpectNoAtom(parsed, "kTyrAtoms_TYM", "HH");
    const auto* tym_oh = FindAtom(parsed, tables, "kTyrAtoms_TYM", "OH");
    ASSERT_NE(nullptr, tym_oh);
    EXPECT_EQ(-1, tym_oh->formal_charge);
    EXPECT_EQ(PolarHKind::NotPolar, tym_oh->polar_h);
}

TEST(TopologySemanticTables, VariantChargeOverridesRunAfterParentSynthesis) {
    const auto parsed = ParseGeneratedAtomNames();
    const auto tables = AllTables();

    const auto* asp_od2 = FindAtom(parsed, tables, "kAspAtoms", "OD2");
    const auto* ash_od2 = FindAtom(parsed, tables, "kAspAtoms_ASH", "OD2");
    ASSERT_NE(nullptr, asp_od2);
    ASSERT_NE(nullptr, ash_od2);
    EXPECT_EQ(-1, asp_od2->formal_charge);
    EXPECT_EQ(0, ash_od2->formal_charge)
        << "ASH must patch ASP's OD2=-1 override back to neutral after "
        << "delegating to the ASP synthesis path";

    const auto* glu_oe2 = FindAtom(parsed, tables, "kGluAtoms", "OE2");
    const auto* glh_oe2 = FindAtom(parsed, tables, "kGluAtoms_GLH", "OE2");
    ASSERT_NE(nullptr, glu_oe2);
    ASSERT_NE(nullptr, glh_oe2);
    EXPECT_EQ(-1, glu_oe2->formal_charge);
    EXPECT_EQ(0, glh_oe2->formal_charge)
        << "GLH must patch GLU's OE2=-1 override back to neutral after "
        << "delegating to the GLU synthesis path";

    const auto* his_nd1 = FindAtom(parsed, tables, "kHisAtoms", "ND1");
    const auto* his_ne2 = FindAtom(parsed, tables, "kHisAtoms", "NE2");
    const auto* hip_nd1 = FindAtom(parsed, tables, "kHisAtoms_HIP", "ND1");
    const auto* hip_ne2 = FindAtom(parsed, tables, "kHisAtoms_HIP", "NE2");
    ASSERT_NE(nullptr, his_nd1);
    ASSERT_NE(nullptr, his_ne2);
    ASSERT_NE(nullptr, hip_nd1);
    ASSERT_NE(nullptr, hip_ne2);
    EXPECT_EQ(0, his_nd1->formal_charge);
    EXPECT_EQ(0, his_ne2->formal_charge);
    EXPECT_EQ(0, hip_nd1->formal_charge);
    EXPECT_EQ(1, hip_ne2->formal_charge)
        << "HIP must retain its variant-specific positive imidazolium "
        << "charge after inheriting the HIS ring synthesis";
}

TEST(TopologySemanticTables, VariantDispatchUsesTheVariantSpecificTables) {
    const auto parsed = ParseGeneratedAtomNames();
    const auto tables = AllTables();

    const auto* ash_hd2 = FindAtom(parsed, tables, "kAspAtoms_ASH", "HD2");
    ASSERT_NE(nullptr, ash_hd2);
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::ASP, 0, IdentityOf(*ash_hd2)));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::ASP, gen::kBaseVariantIdx,
                                     IdentityOf(*ash_hd2)));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::ASP, 1, IdentityOf(*ash_hd2)));

    const auto* glh_he2 = FindAtom(parsed, tables, "kGluAtoms_GLH", "HE2");
    ASSERT_NE(nullptr, glh_he2);
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::GLU, 0, IdentityOf(*glh_he2)));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::GLU, gen::kBaseVariantIdx,
                                     IdentityOf(*glh_he2)));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::GLU, 1, IdentityOf(*glh_he2)));

    const auto* hid_hd1 = FindAtom(parsed, tables, "kHisAtoms_HID", "HD1");
    ASSERT_NE(nullptr, hid_hd1);
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 0, IdentityOf(*hid_hd1)));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::HIS, 1, IdentityOf(*hid_hd1)));
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 2, IdentityOf(*hid_hd1)));

    const auto* hie_he2 = FindAtom(parsed, tables, "kHisAtoms_HIE", "HE2");
    ASSERT_NE(nullptr, hie_he2);
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::HIS, 0, IdentityOf(*hie_he2)));
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 1, IdentityOf(*hie_he2)));
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 2, IdentityOf(*hie_he2)));
}

}  // namespace

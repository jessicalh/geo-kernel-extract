// tests/test_string_barrier.cpp
//
// String barrier discipline test for LegacyAmberTopology semantic
// substrate. Asserts the structural property that chemistry-string
// libraries (gemmi, RDKit) cannot be invoked from the runtime
// library, and that cifpp's use is bounded to the existing
// PDB-loading boundary.
//
// The barrier is enforced at multiple layers; this test backstops
// the discipline:
//
//   1. Linker-level: libnmr_shielding.a contains no RDKit symbols.
//      This is the strongest guarantee -- runtime code physically
//      cannot call RDKit because the symbols are not present.
//
//   2. Header-level: zero src/ headers include cifpp / RDKit / gemmi.
//
//   3. Module-level: cifpp #includes appear in exactly two files
//      (PdbFileReader.cpp, DsspResult.cpp) -- the documented
//      PDB-loading boundary; RDKit and gemmi appear in zero src/
//      files.
//
//   4. Generated-output-level: src/generated/*.cpp contains no
//      std::string literals carrying chemistry data (only typed-enum
//      identifier text, which is compile-time).
//
// See spec/plan/topology-substrate-implementation-plan-2026-05-05.md
// for the architectural rationale.

#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

namespace fs = std::filesystem;

namespace {

// Run a shell command, return stdout as a string. Used to grep
// the source tree.
std::string ShellCapture(const std::string& cmd) {
    std::string out;
    FILE* p = ::popen(cmd.c_str(), "r");
    if (p == nullptr) return out;
    char buf[4096];
    while (std::fgets(buf, sizeof(buf), p) != nullptr) {
        out += buf;
    }
    ::pclose(p);
    return out;
}

// Source root: provided via the build system's NMR_SOURCE_DIR define
// (set in CMakeLists for the test target).
std::string SourceRoot() {
#ifdef NMR_SOURCE_DIR
    return NMR_SOURCE_DIR;
#else
    return ".";
#endif
}

// Whitelisted files where cifpp use is sanctioned (the PDB-loading
// boundary established before this work). New entries require an
// architectural decision, not a code-review approval.
const std::set<std::string>& CifppAllowList() {
    static const std::set<std::string> allow = {
        "src/PdbFileReader.cpp",
        "src/DsspResult.cpp",
    };
    return allow;
}

// Files matching a relative pattern under the source root.
std::vector<std::string> GrepSrcRecursive(const std::string& pattern) {
    const std::string root = SourceRoot();
    // -l: just filenames; -E: extended regex; --include filters by
    // extension so we don't grep build artifacts.
    const std::string cmd =
        "cd " + root + " && grep -lE '" + pattern + "' "
        "--include='*.h' --include='*.cpp' --include='*.hpp' "
        "-r src/ 2>/dev/null || true";
    const std::string out = ShellCapture(cmd);
    std::vector<std::string> lines;
    std::istringstream iss(out);
    for (std::string line; std::getline(iss, line); ) {
        if (!line.empty()) lines.push_back(line);
    }
    return lines;
}


// ============================================================================
// Layer 1 -- Linker-level: RDKit symbols absent from libnmr_shielding.a
// ============================================================================

TEST(StringBarrier, RuntimeLibraryHasZeroRdkitSymbols) {
    const std::string archive = SourceRoot() + "/build/libnmr_shielding.a";
    if (!fs::exists(archive)) {
        GTEST_SKIP() << "build/libnmr_shielding.a not present (run cmake build); "
                     << "test asserts the structural barrier when artifact exists.";
    }
    const std::string cmd = "nm -A " + archive +
                            " 2>/dev/null | grep -c '_ZN5RDKit' || true";
    const std::string out = ShellCapture(cmd);
    int count = std::atoi(out.c_str());
    EXPECT_EQ(0, count) << "RDKit symbols present in libnmr_shielding.a; the "
                          << "string barrier was breached. Inspect with:\n"
                          << "  nm -A " << archive << " | grep _ZN5RDKit";
}


// ============================================================================
// Layer 2 -- Header-level: no chemistry-string libraries in headers
// ============================================================================

TEST(StringBarrier, NoCifppRdkitGemmiInHeaders) {
    // Walk src/ headers (.h / .hpp) only; .cpp use of cifpp is
    // governed by Layer 3. Pattern matches the #include directives
    // (the only way these libraries reach the compilation unit),
    // not identifiers that happen to contain "cifpp" / "rdkit" /
    // "gemmi" as substrings (e.g. SemanticSource enum values).
    const std::string root = SourceRoot();
    const std::string cmd =
        "cd " + root + " && grep -lE "
        "'^[[:space:]]*#[[:space:]]*include[[:space:]]+[<\"]("
        "cif\\+\\+|rdkit/|GraphMol/|RDGeneral/|CIPLabeler/|gemmi/)' "
        "--include='*.h' --include='*.hpp' "
        "-r src/ 2>/dev/null || true";
    const std::string out = ShellCapture(cmd);
    EXPECT_TRUE(out.empty())
        << "Chemistry-string libraries #include'd in src/ headers. The "
           "string barrier requires these libraries appear only at the "
           "loading boundary (specific .cpp files). Offenders:\n" << out;
}


// ============================================================================
// Layer 3 -- Module-level: cifpp confined to allow-list; RDKit/gemmi nowhere
// ============================================================================

TEST(StringBarrier, CifppOnlyInAllowListedTranslationUnits) {
    // Match cifpp #include directives or the cif:: namespace usage.
    // Does NOT match identifiers like "cifpp_CCD_pdbx_aromatic" (an
    // enum value name in SemanticEnums.h that describes the source
    // of a witness, not a chemistry-string boundary breach).
    const auto cifpp_users = GrepSrcRecursive(
        "^[[:space:]]*#[[:space:]]*include[[:space:]]+[<\"]cif\\+\\+|"
        "cif::");
    const auto& allow = CifppAllowList();

    std::vector<std::string> violations;
    for (const auto& f : cifpp_users) {
        if (allow.find(f) == allow.end()) {
            violations.push_back(f);
        }
    }
    EXPECT_TRUE(violations.empty())
        << "cifpp used outside the documented PDB-loading boundary. "
           "Allow-listed files: PdbFileReader.cpp, DsspResult.cpp. "
           "Adding a new cifpp consumer requires an architectural "
           "decision -- update CifppAllowList() in this test plus "
           "the architecture rationale in spec/plan/. Violations:\n"
        << [&violations]{ std::string s; for (const auto& v : violations) s += v + "\n"; return s; }();
}

TEST(StringBarrier, RdkitAndGemmiNeverInRuntimeSrc) {
    // Same regex strategy: match #include directives + namespace
    // usage; do not match bare substrings in identifier names.
    const auto rdkit_users = GrepSrcRecursive(
        "^[[:space:]]*#[[:space:]]*include[[:space:]]+[<\"]("
        "rdkit/|GraphMol/|RDGeneral/|CIPLabeler/)|RDKit::");
    const auto gemmi_users = GrepSrcRecursive(
        "^[[:space:]]*#[[:space:]]*include[[:space:]]+[<\"]gemmi/|gemmi::");
    EXPECT_TRUE(rdkit_users.empty())
        << "RDKit appears in src/ runtime files. RDKit is restricted to "
           "tools/topology/ (the build-time generator). The runtime library "
           "must never link RDKit. Violations:\n"
        << [&rdkit_users]{ std::string s; for (const auto& v : rdkit_users) s += v + "\n"; return s; }();
    EXPECT_TRUE(gemmi_users.empty())
        << "gemmi appears in src/ runtime files. gemmi is not currently used "
           "by the project; if a new chemistry source is being added, the "
           "decision should be reviewed against the string barrier in "
           "spec/plan/topology-substrate-implementation-plan-2026-05-05.md. "
           "Violations:\n"
        << [&gemmi_users]{ std::string s; for (const auto& v : gemmi_users) s += v + "\n"; return s; }();
}


// ============================================================================
// Layer 4 -- Generated output: zero std::string literals carrying chemistry
// ============================================================================

TEST(StringBarrier, GeneratedTablesContainOnlyTypedEnumIdentifiers) {
    const std::string gen = SourceRoot() + "/src/generated/LegacyAmberSemanticTables.cpp";
    if (!fs::exists(gen)) {
        GTEST_SKIP() << "Generated table not present (run "
                        "tools/topology/build_semantic_tables); test "
                        "asserts the typed-only property when artifact exists.";
    }
    std::ifstream f(gen);
    ASSERT_TRUE(f.is_open()) << "cannot open " << gen;

    int line_no = 0;
    int violations = 0;
    std::string violation_msg;
    for (std::string line; std::getline(f, line); ) {
        ++line_no;
        // Skip leading-comment lines (// ...), blank lines, and
        // #include directives (which legitimately carry quoted
        // header paths). The runtime data table itself must contain
        // no string literals; only typed-enum identifiers.
        std::string stripped = line;
        const auto first = stripped.find_first_not_of(" \t");
        if (first == std::string::npos) continue;
        stripped = stripped.substr(first);
        if (stripped.substr(0, 2) == "//") continue;
        if (stripped[0] == '#') continue;
        // The closing-brace lines with trailing comments (`}};  // ...`)
        // and entry rows (`{...},  // 0: N`) carry an end-of-line
        // comment after the literal. Strip the //... tail before
        // scanning for unexpected string literals.
        const auto cmt = stripped.find("//");
        if (cmt != std::string::npos) stripped = stripped.substr(0, cmt);

        // Now scan for any "..." string-literal pattern. The generated
        // file should never contain string literals (only enum
        // identifiers, integer literals, true/false, and braces).
        for (size_t i = 0; i < stripped.size(); ++i) {
            if (stripped[i] == '"') {
                ++violations;
                if (violation_msg.size() < 1000) {
                    violation_msg += "  line " + std::to_string(line_no) + ": " + line + "\n";
                }
                break;  // one report per line
            }
        }
    }
    EXPECT_EQ(0, violations)
        << "Generated table contains string literals (chemistry data leaked "
           "into runtime code). The emitter must produce only typed-enum "
           "identifiers, never std::string literals. Offending lines:\n"
        << violation_msg;
}

}  // namespace

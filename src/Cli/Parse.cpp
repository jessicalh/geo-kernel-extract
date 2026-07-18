/// @file
/// nmr::cli::Parse implementation.

#include "Cli/Parse.h"

#include <cstdlib>
#include <cstring>

namespace nmr::cli {

namespace {

constexpr const char* kFlagPdb           = "--pdb";
constexpr const char* kFlagProtonatedPdb = "--protonated-pdb";
constexpr const char* kFlagOf3           = "--of3";
constexpr const char* kFlagOrca          = "--orca";
constexpr const char* kFlagMutant        = "--mutant";
constexpr const char* kFlagTrajectory    = "--trajectory";

constexpr const char* kFlagHelp1 = "--help";
constexpr const char* kFlagHelp2 = "-h";

bool HasFlag(int argc, char* argv[], const char* name) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], name) == 0) return true;
    }
    return false;
}

/// Return the next argv token after @c name, or empty if @c name is
/// absent or is the final token.
std::string GetArg(int argc, char* argv[], const char* name) {
    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], name) == 0) return argv[i + 1];
    }
    return "";
}

/// Return the first positional argument (no leading dash) after @c name,
/// or empty if not present. Distinct from GetArg because @c --trajectory
/// is followed by a directory positional rather than a flagged value:
/// @c --trajectory DIR --aimnet2 MODEL must read DIR, not @c --aimnet2.
std::string GetPositionalAfter(int argc, char* argv[], const char* name) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], name) != 0) continue;
        for (int j = i + 1; j < argc; ++j) {
            if (argv[j][0] != '-') return argv[j];
        }
    }
    return "";
}

CommonOptions ParseCommon(int argc, char* argv[]) {
    CommonOptions opts;
    opts.output_dir         = GetArg(argc, argv, "--output");
    opts.config_path        = GetArg(argc, argv, "--config");
    opts.aimnet2_model_path = GetArg(argc, argv, "--aimnet2");
    return opts;
}

/// Apply the single-conformation skip flags. MOPAC is the one orthogonal
/// toggle (@c --no-mopac); APBS is always on and home-rolled Coulomb is
/// retired in favour of it, so neither is switchable.
template <typename Mode>
void ApplySingleConfFlags(int argc, char* argv[], Mode& m) {
    if (HasFlag(argc, argv, "--no-mopac")) m.mopac = false;
}

OrcaRunFiles ExpandOrcaRoot(const std::string& root) {
    OrcaRunFiles f;
    f.xyz_path     = root + ".xyz";
    f.prmtop_path  = root + ".prmtop";
    f.nmr_out_path = root + "_nmr.out";
    return f;
}

// --of3 root expands to EXACTLY {root}.prmtop + {root}.inpcrd. Geometry comes
// from the tleap-emitted AMBER restart (rst7/inpcrd) that already sits beside
// the prmtop, so no separate .xyz is required. It never derives, probes,
// validates, or remembers {root}_nmr.out: DFT is not part of OF3's object
// shape. Conformation provenance is honest OpenFold+tleap, supplied here at
// the mode boundary (never hardcoded in the shared loader).
Of3Input ExpandOf3Root(const std::string& root) {
    Of3Input in;
    in.prmtop_path       = root + ".prmtop";
    in.inpcrd_path       = root + ".inpcrd";
    in.prediction_method = "OpenFold+tleap";
    return in;
}

/// Identify the mode flag in argv. Returns the matched flag string or
/// empty if none. Sets @c error if more than one is present (the old
/// parser silently picked the first; this is a deliberate behaviour
/// change).
const char* IdentifyModeFlag(int argc, char* argv[], std::string& error) {
    static constexpr const char* kModes[] = {
        kFlagPdb, kFlagProtonatedPdb, kFlagOf3, kFlagOrca, kFlagMutant,
        kFlagTrajectory,
    };
    const char* matched = nullptr;
    for (const char* m : kModes) {
        if (!HasFlag(argc, argv, m)) continue;
        if (matched != nullptr) {
            error = std::string("multiple mode flags present: ") + matched + " and " + m
                  + " (pick exactly one)";
            return nullptr;
        }
        matched = m;
    }
    return matched;
}

ParseResult ParsePdb(int argc, char* argv[]) {
    ParseResult r;
    PdbMode m;
    m.pdb = GetArg(argc, argv, kFlagPdb);
    if (m.pdb.empty()) {
        r.error = "--pdb requires a file path: --pdb FILE";
        return r;
    }
    const std::string pH_str = GetArg(argc, argv, "--pH");
    if (!pH_str.empty()) m.pH = std::strtod(pH_str.c_str(), nullptr);
    ApplySingleConfFlags(argc, argv, m);
    r.spec   = ModeSpec{std::move(m)};
    r.common = ParseCommon(argc, argv);
    return r;
}

ParseResult ParseProtonatedPdb(int argc, char* argv[]) {
    ParseResult r;
    ProtonatedPdbMode m;
    m.pdb = GetArg(argc, argv, kFlagProtonatedPdb);
    if (m.pdb.empty()) {
        r.error = "--protonated-pdb requires a file path: --protonated-pdb FILE";
        return r;
    }
    ApplySingleConfFlags(argc, argv, m);
    r.spec   = ModeSpec{std::move(m)};
    r.common = ParseCommon(argc, argv);
    return r;
}

ParseResult ParseOf3(int argc, char* argv[]) {
    ParseResult r;
    Of3Mode m;
    const std::string root = GetArg(argc, argv, "--root");
    if (root.empty()) {
        r.error =
            "--of3 requires --root NAME\n"
            "  root expands to only: {root}.prmtop, {root}.inpcrd "
            "(no {root}_nmr.out; --of3 carries no DFT input)";
        return r;
    }
    m.input = ExpandOf3Root(root);
    ApplySingleConfFlags(argc, argv, m);
    r.spec   = ModeSpec{std::move(m)};
    r.common = ParseCommon(argc, argv);
    return r;
}

ParseResult ParseOrca(int argc, char* argv[]) {
    ParseResult r;
    OrcaMode m;
    const std::string root = GetArg(argc, argv, "--root");
    if (root.empty()) {
        r.error =
            "--orca requires --root NAME\n"
            "  root expands to: {root}.xyz, {root}.prmtop, {root}_nmr.out";
        return r;
    }
    m.files = ExpandOrcaRoot(root);
    ApplySingleConfFlags(argc, argv, m);
    r.spec   = ModeSpec{std::move(m)};
    r.common = ParseCommon(argc, argv);
    return r;
}

ParseResult ParseMutant(int argc, char* argv[]) {
    ParseResult r;
    MutantMode m;
    const std::string wt_root  = GetArg(argc, argv, "--wt");
    const std::string ala_root = GetArg(argc, argv, "--ala");
    if (wt_root.empty() || ala_root.empty()) {
        r.error =
            "--mutant requires --wt NAME --ala NAME\n"
            "  each root expands to: {root}.xyz, {root}.prmtop, {root}_nmr.out";
        return r;
    }
    m.wt  = ExpandOrcaRoot(wt_root);
    m.ala = ExpandOrcaRoot(ala_root);
    ApplySingleConfFlags(argc, argv, m);
    r.spec   = ModeSpec{std::move(m)};
    r.common = ParseCommon(argc, argv);
    return r;
}

ParseResult ParseTrajectory(int argc, char* argv[]) {
    ParseResult r;
    TrajectoryMode m;
    m.dir = GetPositionalAfter(argc, argv, kFlagTrajectory);
    if (m.dir.empty()) {
        r.error =
            "--trajectory requires a directory path:\n"
            "  --trajectory DIR\n"
            "  DIR must contain production.tpr, production.trr, production.edr";
        return r;
    }
    m.mopac = HasFlag(argc, argv, "--mopac");
    // The one cadence knob. --stride N processes every N-th TRR frame and
    // emits NPY+PDB on each; default 1 = every frame. No separate
    // mopac/npy/pdb stride exists any more.
    const std::string s_stride = GetArg(argc, argv, "--stride");
    if (!s_stride.empty()) {
        m.stride = std::strtoull(s_stride.c_str(), nullptr, 10);
        if (m.stride == 0) m.stride = 1;
    }
    // Optional window over raw TRR frame indices: [window_start,
    // window_start + window_len). Bounds the SINGLE dispatch loop; stride
    // still applies within. Both-or-neither: presence is tested with
    // HasFlag (a flag may be present yet be the final token, which GetArg
    // reports as empty) so "exactly one supplied" is caught as an error
    // rather than silently treated as no-window.
    const bool has_win_start = HasFlag(argc, argv, "--window-start");
    const bool has_win_len   = HasFlag(argc, argv, "--window-len");
    if (has_win_start != has_win_len) {
        r.error =
            "--window-start N and --window-len M must be given together "
            "(both or neither); window = TRR frames [N, N+M), window-local";
        return r;
    }
    if (has_win_start && has_win_len) {
        // HasFlag confirmed presence, but a flag can be the FINAL token with no
        // value, which GetArg reports as empty -> strtoull("") == 0 ->
        // window_len 0 -> HasWindow() false -> a silently full run. Require a
        // digit string for each value and a window of >= 1 frame. Loud, not
        // silent: this feature was removed once for silent degradation.
        const std::string s_start = GetArg(argc, argv, "--window-start");
        const std::string s_len   = GetArg(argc, argv, "--window-len");
        const auto is_uint = [](const std::string& s) {
            return !s.empty() &&
                   s.find_first_not_of("0123456789") == std::string::npos;
        };
        if (!is_uint(s_start) || !is_uint(s_len)) {
            r.error =
                "--window-start and --window-len each require a non-negative "
                "integer: --window-start N --window-len M "
                "(window = TRR frames [N, N+M))";
            return r;
        }
        m.window_start = std::strtoull(s_start.c_str(), nullptr, 10);
        m.window_len   = std::strtoull(s_len.c_str(), nullptr, 10);
        if (m.window_len == 0) {
            r.error =
                "--window-len must be >= 1 (window = TRR frames "
                "[start, start+len); 0 is an empty window)";
            return r;
        }
    }
    r.spec   = ModeSpec{std::move(m)};
    r.common = ParseCommon(argc, argv);
    return r;
}

}  // namespace


ParseResult Parse(int argc, char* argv[]) {
    ParseResult r;

    if (argc < 2 || HasFlag(argc, argv, kFlagHelp1) || HasFlag(argc, argv, kFlagHelp2)) {
        r.help_requested = true;
        return r;
    }

    const char* mode = IdentifyModeFlag(argc, argv, r.error);
    if (!r.error.empty()) return r;
    if (mode == nullptr) {
        r.error =
            "no mode flag found\n"
            "  valid modes: --pdb, --protonated-pdb, --of3, --orca, --mutant, "
            "--trajectory";
        return r;
    }

    if (std::strcmp(mode, kFlagPdb)           == 0) return ParsePdb(argc, argv);
    if (std::strcmp(mode, kFlagProtonatedPdb) == 0) return ParseProtonatedPdb(argc, argv);
    if (std::strcmp(mode, kFlagOf3)           == 0) return ParseOf3(argc, argv);
    if (std::strcmp(mode, kFlagOrca)          == 0) return ParseOrca(argc, argv);
    if (std::strcmp(mode, kFlagMutant)        == 0) return ParseMutant(argc, argv);
    if (std::strcmp(mode, kFlagTrajectory)    == 0) return ParseTrajectory(argc, argv);

    r.error = std::string("internal: unhandled mode flag ") + mode;
    return r;
}

}  // namespace nmr::cli

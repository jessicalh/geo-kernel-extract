// tools/topology/build_semantic_tables.cpp
//
// Build-time generator for the LegacyAmberTopology semantic tables.
//
// String barrier note: this file is the STRING BOUNDARY in its
// strongest form. Strings flow IN from CCD (cifpp) and from RDKit's
// perception output. They are reconciled here against the typed-enum
// vocabulary in src/SemanticEnums.h and emitted as typed-enum
// literals into src/generated/LegacyAmberSemanticTables.cpp. After
// emission, no string survives. The runtime library libnmr_shielding
// does not link RDKit and never sees this code.
//
// Process log: this generator emits a structured log (NDJSON) of
// every decision it makes -- per (residue, atom, field), what each
// source said, what was chosen, why. The log is committed alongside
// the generated C++ source so the audit trail lives with the data.
// Inspection: the user reads the log once after each generation run;
// committed logs serve as a reproducibility record.

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

// cifpp -- CCD reader at the loading boundary.
#include <cif++.hpp>

// RDKit -- chemistry perception.
#include <GraphMol/GraphMol.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/versions.h>


namespace {

struct Args {
    std::string ccd_path;
    std::string output_cpp_path;
    std::string log_path;
    bool        help        = false;
};

void PrintUsage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s --ccd PATH --output PATH --log PATH\n"
        "\n"
        "Build the LegacyAmberTopology semantic tables.\n"
        "\n"
        "Required:\n"
        "  --ccd PATH      Path to wwPDB Chemical Component Dictionary\n"
        "                  (typically data/ccd/components.cif).\n"
        "  --output PATH   Path to write the generated C++ source\n"
        "                  (typically src/generated/LegacyAmberSemanticTables.cpp).\n"
        "  --log PATH      Path to write the structured process log\n"
        "                  (typically src/generated/LegacyAmberSemanticTables.log.txt).\n"
        "\n"
        "Optional:\n"
        "  --help          Show this message.\n",
        prog);
}

bool ParseArgs(int argc, char** argv, Args& args) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--help" || a == "-h") {
            args.help = true;
            return true;
        }
        if (i + 1 >= argc) {
            std::fprintf(stderr, "ERROR: %s requires a value\n", a.c_str());
            return false;
        }
        const std::string v = argv[++i];
        if      (a == "--ccd")    args.ccd_path        = v;
        else if (a == "--output") args.output_cpp_path = v;
        else if (a == "--log")    args.log_path        = v;
        else {
            std::fprintf(stderr, "ERROR: unknown argument %s\n", a.c_str());
            return false;
        }
    }
    if (args.ccd_path.empty() || args.output_cpp_path.empty() || args.log_path.empty()) {
        std::fprintf(stderr, "ERROR: --ccd, --output, --log are all required\n");
        return false;
    }
    return true;
}


// Process log: structured plain-text record. Each line is either a
// section header (==== SECTION NAME ====) or a key=value record. The
// format is grep-able and human-readable; machine consumers can
// parse with awk or a simple state machine if needed.
//
// Open once at start of run, write each event, close at end. The
// log is committed alongside the generated C++ source.
class ProcessLog {
public:
    explicit ProcessLog(const std::string& path) : out_(path) {
        if (!out_.is_open()) {
            std::fprintf(stderr, "ERROR: cannot open log file %s\n", path.c_str());
            std::exit(2);
        }
        out_ << "# LegacyAmberTopology semantic tables -- generator process log\n"
                "# Format: section headers (==== SECTION ====), key=value records,\n"
                "# blank lines for grouping. grep-able by section / key.\n"
                "#\n"
                "# This log is committed alongside src/generated/LegacyAmberSemanticTables.cpp\n"
                "# as the reproducibility audit trail. Inspect once per regeneration.\n"
                "#\n";
    }

    void Section(const std::string& name) {
        out_ << "\n==== " << name << " ====\n";
    }
    void KV(const std::string& key, const std::string& value) {
        out_ << key << " = " << value << "\n";
    }
    void KV(const std::string& key, int value) {
        out_ << key << " = " << value << "\n";
    }
    void Note(const std::string& msg) {
        out_ << "# " << msg << "\n";
    }

private:
    std::ofstream out_;
};


// Quick smoke test: open CCD via cifpp, count datablocks. Validates
// that the cifpp link is live.
int CcdSmokeTest(const std::string& ccd_path, ProcessLog& log) {
    log.Section("ccd-smoke-test");
    log.KV("ccd_path", ccd_path);

    cif::file ccd;
    try {
        ccd.load(ccd_path);
    } catch (const std::exception& e) {
        log.KV("status", "FAILED");
        log.KV("error", e.what());
        std::fprintf(stderr, "ERROR loading CCD: %s\n", e.what());
        return 1;
    }

    int n_blocks = 0;
    for (const auto& db : ccd) {
        (void)db;
        ++n_blocks;
        if (n_blocks > 10000) break;  // sanity guard
    }
    log.KV("datablocks_seen_first_10000", n_blocks);
    log.KV("status", "OK");
    return 0;
}


// Quick smoke test: instantiate an RDKit molecule. Validates that
// the RDKit link is live and runtime is callable.
int RdkitSmokeTest(ProcessLog& log) {
    log.Section("rdkit-smoke-test");
    log.KV("rdkit_version", RDKit::rdkitVersion);
    log.KV("rdkit_build", RDKit::rdkitBuild);

    RDKit::RWMol mol;
    auto a = new RDKit::Atom(6);
    mol.addAtom(a, true, true);
    log.KV("trivial_mol_atom_count", static_cast<int>(mol.getNumAtoms()));
    log.KV("status", "OK");
    return 0;
}


}  // namespace


int main(int argc, char** argv) {
    Args args;
    if (!ParseArgs(argc, argv, args) || args.help) {
        PrintUsage(argv[0]);
        return args.help ? 0 : 2;
    }

    ProcessLog log(args.log_path);
    log.Section("run-info");
    log.KV("generator", "build_semantic_tables");
    log.KV("ccd_path", args.ccd_path);
    log.KV("output_cpp", args.output_cpp_path);
    log.KV("log_path", args.log_path);
    log.KV("scaffolding_step", "1 of N");
    log.Note("Step 1 lays scaffolding only: cifpp + RDKit smoke tests, log shape.");
    log.Note("Subsequent steps add CCD entry parsing, RDKit perception, reconciliation,");
    log.Note("and generated-C++ emission. Each step lands as its own commit.");

    int rc = CcdSmokeTest(args.ccd_path, log);
    if (rc != 0) return rc;

    rc = RdkitSmokeTest(log);
    if (rc != 0) return rc;

    log.Section("done");
    log.KV("status", "OK");
    log.Note("Scaffolding only. No tables generated yet.");

    std::printf("OK -- scaffolding step complete. See log at %s\n",
                args.log_path.c_str());
    return 0;
}

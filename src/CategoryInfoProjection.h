#pragma once
//
// CategoryInfoProjection -- output-side per-atom comprehensive categorical
// record. Emits one structured NPY per protein carrying every invariant
// categorical fact about each atom: names across naming systems
// (AMBER / IUPAC / BMRB at atom AND residue level), mechanical identity
// (locant, branch, di_index, backbone_role), chemistry classification
// (prochiral, planar group/stereo, polar-H kind, ring position,
// pseudoatom membership, aromatic, formal_charge, is_exchangeable), and
// provenance for any external-table lookup.
//
// Singleton, fixed-shape (no inheritance, no virtuals). All state is
// file-local in the .cpp; the public surface is static methods. Same
// shape as FramePdbEmitter.
//
//   Configure(config) -- assemble: parse atom_nom.tbl into the lookup
//                         table when provided. Called once at startup from
//                         Session::LoadFromToml.
//   WriteFeatures(protein, output_dir)
//                         -- emit one structured NPY per protein. Without
//                         atom_nom.tbl, emits AMBER names as lookup
//                         fallbacks; without AtomSemanticTable substrate,
//                         emits zero-valued substrate fields.
//   Reset()             -- clear configuration; for tests.
//
// Primary reads (no model mutation):
//   - protein.LegacyAmber().AtomSemantic()  (typed substrate per atom)
//   - protein.LegacyAmber().AtomtypeString()
//   - protein.AtomAt(i)                     (element, residue_index, pdb_atom_name)
//   - protein.ResidueAt(j)                  (identity, sequence address, terminal_state)
//   - GetAminoAcidType(...)                 (three_letter_code, variants)
//   - parsed atom_nom.tbl                   (BMRB column, SC stereo column)
//
// Writes:
//   - one structured NPY at output_dir/atoms_category_info.npy. Shape
//     (N,) with a packed structured dtype. Mixed S8 / S4 / S1 / i1 / i4
//     columns.
//
// Architectural rule (memory feedback_naming_input_output_asymmetry):
// input-side and output-side naming systems are NEVER glued together.
// Input feeds physics (no tolerance, fail-loud). Output feeds ML
// matching (logged-fallback acceptable).
//
// The S8 / S4 / S1 string columns have a fail-loud length guard at
// emission: a value exceeding the column's char width aborts with a
// diagnostic naming the offending row. Per PATTERNS §9 fprintf+abort,
// no exceptions.
//
// Deliberately not a ConformationResult / TrajectoryResult. Holds no
// per-frame state; participates in no dependency graph; emits no H5.
// Per-output-directory projection for a Protein, called from entry
// points and trajectory output boundaries.
//

#include <cstddef>
#include <filesystem>
#include <map>
#include <string>

namespace nmr {

class Protein;
enum class AminoAcid : int;

class CategoryInfoProjection {
public:
    struct Config {
        std::filesystem::path atom_nom_tbl;  // empty = no external lookup table
    };

    // Called once at startup. Reads + parses atom_nom.tbl when configured.
    // Subsequent calls clear the existing state before parsing; same idiom
    // as FramePdbEmitter.
    static void Configure(Config config);

    // Emit the structured NPY for one protein. Returns the number of
    // arrays written (1 when emission succeeded, 0 on write failure).
    // Idempotent; safe to call multiple times.
    static int WriteFeatures(const Protein& protein,
                              const std::string& output_dir);

    // Clear configuration. For tests.
    static void Reset();

    // Test introspection.
    static bool IsActive();

    // Per-atom queries — for tests and downstream callers (viewer,
    // h5-reader, debugging). All const, all read-only on the protein.
    // Empty string when no projection match (and Configure has been
    // called). WriteFeatures records misses; these query methods are
    // read-only and do not update MissLog.
    static std::string IupacAtomName(const Protein&, std::size_t atom_index);
    static std::string BmrbAtomName(const Protein&, std::size_t atom_index);
    static std::string BmrbStereoLabel(const Protein&, std::size_t atom_index);

    // Per-residue 3-letter / 1-letter projections. Pure functions;
    // no atom_nom.tbl involvement. Static-map lookups.
    static std::string AmberResidueThreeLetter(AminoAcid type, int protonation_variant_index);
    static std::string IupacResidueThreeLetter(AminoAcid type);   // == BMRB for standard 20
    static char        ResidueOneLetter(AminoAcid type);

    // Diagnostic — total miss tally for this projection's lifetime.
    // Keyed by "RES:NAME" for atom-name misses. Cleared on Reset().
    static const std::map<std::string, int>& MissLog();

    CategoryInfoProjection() = delete;
    CategoryInfoProjection(const CategoryInfoProjection&) = delete;
    CategoryInfoProjection& operator=(const CategoryInfoProjection&) = delete;
};

}  // namespace nmr

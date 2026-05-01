#pragma once
//
// AmberChargeResolver: typed dispatch for AMBER charge sources.
//
// Single entry point for AMBER-source loaders to obtain a ChargeSource.
// Loaders populate Protein + ProteinBuildContext, then call
// ResolveAmberChargeSource(...) to get the right concrete source for
// the typed state. Three branches:
//
//   1. build_context.prmtop_path non-empty
//      → PrmtopChargeSource over that file (--orca / --mutant paths).
//   2. AnalyzeFlatTableCoverage returns Satisfiable
//      → ParamFileChargeSource over the flat ff14SB table
//        (standard --pdb / --protonated-pdb proteins).
//   3. Verdict non-satisfiable AND policy permits preparation
//      → AmberPreparedChargeSource (runtime tleap pipeline).
//
// Branch 3 under FailOnUnsupportedTerminalVariants returns nullptr
// with a typed error.
//
// AnalyzeFlatTableCoverage is the single source of truth for "can
// ParamFileChargeSource satisfy this protein." Loaders should not
// duplicate this logic.
//
// AmberPreparationPolicy values:
//   FailOnUnsupportedTerminalVariants                — default; loud fail
//   UseStockTermini                                  — reserved
//   UseCappedFragmentsForUnsupportedTerminalVariants — ACE/NME caps
//                                                       around ASH/CYM/GLH/LYN
//                                                       at unsupported termini
//

#include "ChargeSource.h"
#include "Protein.h"
#include <string>
#include <vector>

namespace nmr {

// ============================================================================
// AmberPreparationPolicy
//
// Selects what AmberPreparedChargeSource does when the flat ff14SB table
// cannot represent the protein. The default for --pdb / --protonated-pdb
// is FailOnUnsupportedTerminalVariants (current behaviour).
// ============================================================================

enum class AmberPreparationPolicy {
    UseStockTermini,
    UseCappedFragmentsForUnsupportedTerminalVariants,
    FailOnUnsupportedTerminalVariants
};

const char* AmberPreparationPolicyName(AmberPreparationPolicy policy);


// ============================================================================
// AmberSourceConfig
//
// What the resolver needs to make a decision. flat_table_path is required;
// tleap_path / work_dir are required only when the resolver chooses
// AmberPreparedChargeSource (step 4+).
// ============================================================================

struct AmberSourceConfig {
    std::string flat_table_path;
    AmberPreparationPolicy preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;
    std::string tleap_path;
    std::string work_dir;
};


// ============================================================================
// AmberFlatTableCoverageVerdict
//
// Pure-function output of AnalyzeFlatTableCoverage. Carries the kind of
// the first failure (or Satisfiable) and the full list of failures, so
// step-6 capping can decide per-end without re-walking the protein.
// ============================================================================

enum class AmberFlatTableCoverageKind {
    Satisfiable,
    UnsupportedTerminalVariant,
    UnsupportedResidue,
    MissingAtomName
};

const char* AmberFlatTableCoverageKindName(AmberFlatTableCoverageKind kind);

struct AmberFlatTableCoverageFailure {
    AmberFlatTableCoverageKind kind = AmberFlatTableCoverageKind::Satisfiable;
    std::string terminal_token;          // INTERNAL / NTERM / CTERM / NCTERM
    std::string ff_residue_name;         // ALA / HID / CYX / ASH / etc.
    std::string atom_name;               // populated for MissingAtomName
    int residue_sequence_number = 0;
    std::string chain_id;
};

struct AmberFlatTableCoverageVerdict {
    bool ok = true;
    std::vector<AmberFlatTableCoverageFailure> failures;

    bool Ok() const { return ok; }

    // First failure, when Ok() is false. Behaviour is undefined when Ok().
    const AmberFlatTableCoverageFailure& FirstFailure() const;

    // Human-readable summary suitable for error_out / methods text.
    // Always contains the phrase "no canonical fallback" when !Ok().
    std::string Detail() const;
};


// ============================================================================
// AnalyzeFlatTableCoverage
//
// Pure function: given a Protein with typed terminal_state and
// protonation_variant_index already resolved, and a path to the ff14SB
// flat parameter file, decide whether ParamFileChargeSource can satisfy
// every atom. Reports the first failure plus any later failures it can
// see in the same single walk.
//
// Construction-boundary note: this function performs string lookups on
// (terminal_token, ff_resname, atom_name) keys. That is allowed because
// it operates against a force-field table that is itself string-keyed
// — the same boundary discipline as ParamFileChargeSource. Calculators
// must never call this function.
// ============================================================================

AmberFlatTableCoverageVerdict AnalyzeFlatTableCoverage(
    const Protein& protein,
    const std::string& flat_table_path);


// ============================================================================
// ResolveAmberChargeSource
//
// Single dispatch point for AMBER charge sources. Loaders that handle
// AMBER chemistry (BuildFromPdb, BuildFromProtonatedPdb, BuildFromOrca,
// BuildFromMutant) call this instead of constructing ChargeSource
// subclasses directly.
//
// Branch 1: build_context.prmtop_path is non-empty
//     → PrmtopChargeSource over that file. Never re-tleap'd.
//       --orca / --mutant always take this branch by precondition;
//       a missing prmtop_path on those paths must already have been
//       rejected by the loader.
//
// Branch 2: AnalyzeFlatTableCoverage returns Satisfiable
//     → ParamFileChargeSource over config.flat_table_path.
//
// Branch 3: verdict is non-satisfiable AND policy permits preparation
//     → AmberPreparedChargeSource (lands in step 4; returns nullptr
//       with a step-numbered error until then).
//
// Returns nullptr on dispatch error, populating error_out with a
// typed description.
// ============================================================================

class ChargeSource;

std::unique_ptr<ChargeSource> ResolveAmberChargeSource(
    const Protein& protein,
    const ProteinBuildContext& build_context,
    const AmberSourceConfig& config,
    std::string& error_out);

}  // namespace nmr

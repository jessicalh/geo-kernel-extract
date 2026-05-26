#pragma once
//
// ResolveAmberChargeSource is the dispatch point for AMBER charge sources:
// upstream PRMTOP first, flat ff14SB table when coverage is complete, and
// runtime tleap preparation only when policy permits it.
//

#include "ChargeSource.h"
#include "Protein.h"
#include <string>
#include <vector>

namespace nmr {

enum class AmberPreparationPolicy {
    UseStockTermini,
    UseCappedFragmentsForUnsupportedTerminalVariants,
    FailOnUnsupportedTerminalVariants
};

const char* AmberPreparationPolicyName(AmberPreparationPolicy policy);


// flat_table_path is always required; tleap_path and work_dir are used only
// when runtime PRMTOP preparation is selected.
struct AmberSourceConfig {
    std::string flat_table_path;
    AmberPreparationPolicy preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;
    std::string tleap_path;
    std::string work_dir;
};


// Coverage verdict for the flat ff14SB table. Failures retain enough typed
// residue/end information for capped-fragment preparation.
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


// Pure function: given a Protein with typed terminal_state and
// protonation_variant_index already resolved, and a path to the ff14SB
// flat parameter file, decide whether ParamFileChargeSource can satisfy
// every atom. Reports the first failure plus any later failures it can
// see in the same single walk.
//
// Construction-boundary note: this function performs string lookups on
// (terminal_token, ff_resname, atom_name) keys. That is allowed because
// it operates against a force-field table that is itself string-keyed
// — the same boundary discipline as ParamFileChargeSource.

AmberFlatTableCoverageVerdict AnalyzeFlatTableCoverage(
    const Protein& protein,
    const std::string& flat_table_path);


// Dispatch rules:
//   - build_context.prmtop_path selects PrmtopChargeSource and is not re-tleap'd.
//   - satisfiable flat-table coverage selects ParamFileChargeSource.
//   - non-satisfiable coverage selects AmberPreparedChargeSource only when
//     the policy permits preparation.
// Returns nullptr on dispatch error, populating error_out with a
// typed description.

class ChargeSource;

std::unique_ptr<ChargeSource> ResolveAmberChargeSource(
    const Protein& protein,
    const ProteinBuildContext& build_context,
    const AmberSourceConfig& config,
    std::string& error_out);

}  // namespace nmr

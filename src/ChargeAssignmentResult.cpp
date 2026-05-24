#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "Protein.h"
#include "OperationLog.h"
#include <cstdio>
#include <cstdlib>

namespace nmr {


// ============================================================================
// Typed factory: assign charges from any ChargeSource.
// ============================================================================

std::unique_ptr<ChargeAssignmentResult> ChargeAssignmentResult::Compute(
        ProteinConformation& conf,
        const ChargeSource& source) {

    OperationLog::Scope scope("ChargeAssignmentResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " source=" + source.Describe());

    const Protein& protein = conf.ProteinRef();

    if (!protein.HasForceFieldCharges()) {
        std::string error;
        Protein& mutable_protein = const_cast<Protein&>(protein);
        if (!mutable_protein.PrepareForceFieldCharges(source, conf, error)) {
            OperationLog::Error("ChargeAssignmentResult::Compute",
                "charge table preparation failed: " + error);
            return nullptr;
        }
    }

    const ForceFieldChargeTable& table = protein.ForceFieldCharges();
    if (table.AtomCount() != conf.AtomCount()) {
        OperationLog::Error("ChargeAssignmentResult::Compute",
            "charge table has " + std::to_string(table.AtomCount()) +
            " entries, expected " + std::to_string(conf.AtomCount()));
        return nullptr;
    }

    auto result = std::make_unique<ChargeAssignmentResult>();
    result->conf_ = &conf;
    result->charge_table_ = &table;
    result->source_ = table.SourceDescription();

    // Projection only: prepared charges/PB radii live on Protein.
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        ca.partial_charge = table.PartialChargeAt(ai);
        ca.pb_radius = table.PbRadiusAt(ai);
    }

    result->total_charge_ = table.TotalCharge();
    result->assigned_count_ = table.AssignedCount();
    result->unassigned_count_ = table.UnassignedCount();

    OperationLog::Info(LogCharges, "ChargeAssignmentResult::Compute",
        "source=" + result->source_ +
        " assigned=" + std::to_string(result->assigned_count_) +
        " unassigned=" + std::to_string(result->unassigned_count_) +
        " total_charge=" + std::to_string(result->total_charge_));

    return result;
}

// ============================================================================
// Compute: ff14SB from param file (convenience, delegates to typed path)
// ============================================================================

std::unique_ptr<ChargeAssignmentResult> ChargeAssignmentResult::Compute(
        ProteinConformation& conf,
        const std::string& param_file_path) {

    ParamFileChargeSource source(param_file_path);
    return Compute(conf, source);
}


double ChargeAssignmentResult::ChargeAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).partial_charge;
}

double ChargeAssignmentResult::PbRadiusAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).pb_radius;
}

const ForceFieldChargeTable& ChargeAssignmentResult::ChargeTable() const {
    if (!charge_table_) {
        fprintf(stderr, "FATAL: ChargeAssignmentResult has no charge table.\n");
        std::abort();
    }
    return *charge_table_;
}

}  // namespace nmr

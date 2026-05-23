#include "ProteinConformation.h"
#include "OperationLog.h"
#include <cstdio>

namespace nmr {

ProteinConformation::ProteinConformation(
    const Protein* protein,
    std::vector<Vec3> positions,
    std::string description)
    : protein_(protein)
    , description_(std::move(description))
    , positions_(std::move(positions))
{
    // Build ConformationAtom vector from positions
    atoms_.reserve(positions_.size());
    for (const auto& pos : positions_) {
        atoms_.emplace_back(ConformationAtom(pos));
}
}


bool ProteinConformation::AttachResult(std::unique_ptr<ConformationResult> result) {
    if (!result) return false;

    std::string const name = result->Name();
    // Extract the referent for typeid: typeid(*result) is potentially-
    // evaluated; the named reference makes the polymorphic dispatch explicit.
    const ConformationResult& result_ref = *result;
    std::type_index const tid(typeid(result_ref));

    // Singleton check: already attached?
    if (results_.find(tid) != results_.end()) {
        OperationLog::Log(OperationLog::Level::Warning, LogResultAttach,
                          "AttachResult",
                          "rejected " + name + ": already attached");
        return false;
    }

    // Check dependencies
    for (const auto& dep : result->Dependencies()) {
        if (results_.find(dep) == results_.end()) {
            // Build diagnostic message listing what IS attached
            std::string attached;
            for (const auto& kv : results_) {
                if (!attached.empty()) attached += ", ";
                attached += kv.second->Name();
            }
            std::string msg;
            msg.reserve(name.size() + attached.size() + 64);
            msg.append("rejected ");
            msg.append(name);
            msg.append(": missing dependency. Attached: [");
            msg.append(attached);
            msg.append("]");
            OperationLog::Log(OperationLog::Level::Warning, LogResultAttach, "AttachResult", msg);
            return false;
        }
    }

    OperationLog::Log(OperationLog::Level::Info, LogResultAttach,
                      "AttachResult",
                      "attached " + name + " to conformation");

    results_[tid] = std::move(result);
    return true;
}


void ProteinConformation::ForceAttachResultForTesting(
        std::unique_ptr<ConformationResult> result) {
    if (!result) return;
    // Extract the referent for typeid: see AttachResult comment.
    const ConformationResult& result_ref = *result;
    const std::type_index tid(typeid(result_ref));
    results_[tid] = std::move(result);
}


// ============================================================================
// CrystalConformation
// ============================================================================

CrystalConformation::CrystalConformation(
    const Protein* protein,
    std::vector<Vec3> positions,
    double resolution,
    double r_factor,
    double temperature,
    std::string pdb_id)
    : ProteinConformation(protein, std::move(positions), "crystal:" + pdb_id)
    , resolution_angstroms_(resolution)
    , r_factor_(r_factor)
    , temperature_kelvin_(temperature)
    , pdb_id_(std::move(pdb_id))
{}


// ============================================================================
// PredictionConformation
// ============================================================================

PredictionConformation::PredictionConformation(
    const Protein* protein,
    std::vector<Vec3> positions,
    std::string method,
    double confidence)
    : ProteinConformation(protein, std::move(positions),
                          "prediction:" + method)
    , method_(std::move(method))
    , confidence_(confidence)
{}


// ============================================================================
// MDFrameConformation
// ============================================================================

MDFrameConformation::MDFrameConformation(
    const Protein* protein,
    std::vector<Vec3> positions,
    int walker,
    double time_ps,
    double weight,
    double rmsd_nm,
    double rg_nm)
    : ProteinConformation(protein, std::move(positions),
                          "md_frame:w" + std::to_string(walker) +
                          ":t" + std::to_string(time_ps))
    , walker_(walker)
    , time_ps_(time_ps)
    , weight_(weight)
    , rmsd_nm_(rmsd_nm)
    , rg_nm_(rg_nm)
{}


// ============================================================================
// DerivedConformation
// ============================================================================

DerivedConformation::DerivedConformation(
    const Protein* protein,
    std::vector<Vec3> positions,
    std::string derivation_description)
    : ProteinConformation(protein, std::move(positions),
                          "derived:" + derivation_description)
    , derivation_description_(std::move(derivation_description))
{}

}  // namespace nmr

#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "Protein.h"
#include "Atom.h"
#include "Bond.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include <cstdio>
#include <cstdlib>
#include <vector>

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

    // Projection: prepared charges/PB radii live on Protein. For atoms carrying
    // the flat placeholder PB radius (TPR path, which ships no real radii),
    // substitute mbondi2 — the same set our leap/prmtop path uses — so APBS gets
    // a physical dielectric boundary instead of the 1.5 A placeholder. mbondi2's
    // hydrogen rule needs connectivity, which is live here on the built Protein.
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        ca.partial_charge = table.PartialChargeAt(ai);
        if (table.Values()[ai].status == ChargeAssignmentStatus::PlaceholderPbRadius) {
            const Atom& atom = protein.AtomAt(ai);
            bool h_on_n = false;
            if (atom.element == Element::H) {
                for (size_t bi : atom.bond_indices) {
                    const Bond& b = protein.BondAt(bi);
                    size_t other = (b.atom_index_a == ai) ? b.atom_index_b : b.atom_index_a;
                    if (protein.AtomAt(other).element == Element::N) { h_on_n = true; break; }
                }
            }
            ca.pb_radius = Mbondi2Radius(atom.element, h_on_n);
        } else {
            ca.pb_radius = table.PbRadiusAt(ai);
        }
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


int ChargeAssignmentResult::WriteFeatures(const ProteinConformation& conf,
                                          const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    std::vector<double> charges(N, 0.0);
    std::vector<double> radii(N, 0.0);
    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        charges[i] = ca.partial_charge;
        radii[i] = ca.pb_radius;
    }

    int written = 0;
    if (NpyWriter::WriteFloat64(output_dir + "/ff_partial_charge.npy",
                                charges.data(), N)) ++written;
    if (NpyWriter::WriteFloat64(output_dir + "/ff_pb_radius.npy",
                                radii.data(), N)) ++written;
    return written;
}

}  // namespace nmr

#pragma once
//
// Raw audit of typed sidechain polar-H donors that Larsen 2015 Table 2
// does not model.  This result deliberately does not add shielding: it
// records donor/acceptor geometry so the omitted chemistry remains visible
// after the extraction freeze.

#include "ConformationResult.h"
#include "LarsenHBondGrid.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;

class LarsenSidechainDonorAuditResult : public ConformationResult {
public:
    struct DonorAtomRecord {
        std::int32_t is_sidechain_polar_h = 0;
        std::int32_t polar_h_kind = 0;
        std::int32_t parent_atom = -1;
        std::int32_t residue_index = -1;
        std::int32_t n_acceptor_candidates_3p5A = 0;
        std::int32_t n_geometry_pass = 0;
    };

    struct CandidateRecord {
        std::size_t donor_atom = SIZE_MAX;
        std::size_t donor_residue = SIZE_MAX;
        std::uint8_t polar_h_kind = 0;
        std::size_t parent_atom = SIZE_MAX;
        std::size_t acceptor_atom = SIZE_MAX;
        std::size_t acceptor_residue = SIZE_MAX;
        HBondAcceptorClass acceptor_class =
            HBondAcceptorClass::BackboneCarbonyl;
        double h_acceptor_distance_A = 0.0;
        double parent_acceptor_distance_A = 0.0;
        double angle_parent_h_acceptor_deg = 0.0;
        double rho_deg = 0.0;
        bool passes_geometry = false;
        bool modeled_by_larsen_table2 = false;
    };

    std::string Name() const override {
        return "LarsenSidechainDonorAuditResult";
    }
    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<LarsenSidechainDonorAuditResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    const std::vector<DonorAtomRecord>& DonorAtoms() const {
        return donor_atoms_;
    }
    const std::vector<CandidateRecord>& Candidates() const {
        return candidates_;
    }

private:
    std::vector<DonorAtomRecord> donor_atoms_;
    std::vector<CandidateRecord> candidates_;
};

}  // namespace nmr


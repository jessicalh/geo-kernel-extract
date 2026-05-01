#pragma once
//
// LegacyAmberTopology: explicit topology contract for the current
// Amber-derived calculator language.
//
// This is the string barrier output. Loaders and construction code may use
// names to build this object; calculators consume its typed/indexed facts.
//

#include "ProteinTopology.h"
#include "CovalentTopology.h"
#include <memory>

namespace nmr {

class LegacyAmberTopology final : public ProteinTopology {
public:
    LegacyAmberTopology(size_t atom_count,
                        size_t residue_count,
                        std::unique_ptr<CovalentTopology> bonds);

    ProteinTopologyKind Kind() const override {
        return ProteinTopologyKind::LegacyAmber;
    }
    std::string_view Name() const override { return "LegacyAmberTopology"; }
    size_t AtomCount() const override { return atom_count_; }
    size_t ResidueCount() const override { return residue_count_; }

    const CovalentTopology& Bonds() const { return *bonds_; }

    size_t BondCount() const { return bonds_->BondCount(); }
    const Bond& BondAt(size_t i) const { return bonds_->BondAt(i); }
    const std::vector<Bond>& BondList() const { return bonds_->Bonds(); }
    const std::vector<size_t>& BondIndicesFor(size_t atom_index) const {
        return bonds_->BondIndicesFor(atom_index);
    }
    size_t HydrogenParentOf(size_t atom_index) const {
        return bonds_->HydrogenParentOf(atom_index);
    }

private:
    size_t atom_count_ = 0;
    size_t residue_count_ = 0;
    std::unique_ptr<CovalentTopology> bonds_;
};

}  // namespace nmr

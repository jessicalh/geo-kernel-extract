#include "LegacyAmberTopology.h"

#include <cstdio>
#include <cstdlib>

namespace nmr {

LegacyAmberTopology::LegacyAmberTopology(
        size_t atom_count,
        size_t residue_count,
        std::unique_ptr<CovalentTopology> bonds)
    : atom_count_(atom_count)
    , residue_count_(residue_count)
    , bonds_(std::move(bonds)) {
    if (!bonds_) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology requires a CovalentTopology.\n");
        std::abort();
    }
}

void LegacyAmberTopology::AttachAmberFFData(AmberFFData data) {
    if (amber_ff_data_.has_value()) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology::AttachAmberFFData called twice. "
            "Enrichment is one-shot.\n");
        std::abort();
    }
    amber_ff_data_ = std::move(data);
}

const AmberFFData& LegacyAmberTopology::FFData() const {
    if (!amber_ff_data_.has_value()) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology::FFData called without prior "
            "AttachAmberFFData. Check HasAmberFFData() first.\n");
        std::abort();
    }
    return *amber_ff_data_;
}

}  // namespace nmr

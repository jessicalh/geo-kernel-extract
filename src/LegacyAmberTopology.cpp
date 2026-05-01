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

}  // namespace nmr

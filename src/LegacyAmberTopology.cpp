#include "LegacyAmberTopology.h"

#include <cstdio>
#include <cstdlib>

namespace nmr {

LegacyAmberTopology::LegacyAmberTopology(
        size_t atom_count,
        size_t residue_count,
        std::unique_ptr<CovalentTopology> bonds,
        LegacyAmberInvariants invariants)
    : atom_count_(atom_count)
    , residue_count_(residue_count)
    , bonds_(std::move(bonds))
    , mass_(std::move(invariants.mass))
    , ff_atom_type_index_(std::move(invariants.ff_atom_type_index))
    , ptype_(std::move(invariants.ptype))
    , atomtype_string_(std::move(invariants.atomtype_string))
    , exclusions_(std::move(invariants.exclusions))
    , fudge_qq_(invariants.fudge_qq)
    , rep_pow_(invariants.rep_pow)
    , atnr_(invariants.atnr)
    , num_non_perturbed_(invariants.num_non_perturbed) {
    if (!bonds_) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology requires a CovalentTopology.\n");
        std::abort();
    }
}

}  // namespace nmr

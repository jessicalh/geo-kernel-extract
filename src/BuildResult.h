#pragma once
//
// BuildResult: what a protein-build loader returns — the Protein, a
// ChargeSource matched to its atoms, the integer net_charge, and
// ok/error status. Check Ok() before use.
//

#include "ChargeSource.h"
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class Protein;

struct BuildResult {
    std::unique_ptr<Protein> protein;
    std::unique_ptr<ChargeSource> charges;
    int net_charge = 0;
    bool ok = false;
    std::string error;

    // Per-conformation names from source filenames. Currently unpopulated:
    // no active builder sets it.
    std::vector<std::string> pose_names;

    bool Ok() const { return ok && protein != nullptr; }
    explicit operator bool() const { return Ok(); }
};

}  // namespace nmr

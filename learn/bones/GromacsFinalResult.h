#pragma once
//
// GromacsFinalResult: finalization after trajectory processing.
//
// Runs at the end of a trajectory. Writes the {protein_id}.h5
// master file from conformation 0 + Protein topology. Moves
// winner frames from temp to output. Cleans up temp.
//
// For now: logs what was processed. .h5 writing comes when
// HighFive is integrated.
//

#include "GromacsProtein.h"
#include <string>

namespace nmr {

class GromacsFinalResult {
public:
    explicit GromacsFinalResult(GromacsProtein& gp);

    // Finalize: log summary, move winners, write .h5, clean temp.
    // For now: just logs what happened.
    bool Finalize(const std::string& output_dir,
                  const std::string& temp_dir);

    const std::string& error() const { return error_; }

private:
    GromacsProtein& gp_;
    std::string error_;
};

}  // namespace nmr

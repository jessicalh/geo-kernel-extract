#include "GromacsFinalResult.h"
#include "OperationLog.h"

#include <filesystem>

namespace fs = std::filesystem;

namespace nmr {

GromacsFinalResult::GromacsFinalResult(GromacsProtein& gp)
    : gp_(gp)
{}


bool GromacsFinalResult::Finalize(
        const std::string& output_dir,
        const std::string& temp_dir) {

    OperationLog::Scope scope("GromacsFinalResult::Finalize", gp_.protein_id());

    const auto& paths = gp_.frame_paths();

    OperationLog::Info(LogCalcOther, "GromacsFinalResult::Finalize",
        gp_.protein_id() + ": " +
        std::to_string(paths.size()) + " frames processed, " +
        std::to_string(gp_.protein().AtomCount()) + " atoms, " +
        "conformation 0 in Protein vector, " +
        std::to_string(gp_.protein().ConformationCount()) +
        " conformations in Protein");

    // TODO: write {protein_id}.h5 from conformation 0 + topology
    // TODO: select winners, move winner dirs from temp to output
    // TODO: clean up temp

    return true;
}

}  // namespace nmr

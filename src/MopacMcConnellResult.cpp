#include "MopacMcConnellResult.h"

#include "McConnellResult.h"
#include "MopacResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"

namespace nmr {

std::vector<std::type_index> MopacMcConnellResult::Dependencies() const {
    return {
        std::type_index(typeid(MopacResult)),
        std::type_index(typeid(McConnellResult))
    };
}


std::unique_ptr<MopacMcConnellResult> MopacMcConnellResult::Compute(
        ProteinConformation& conf) {
    OperationLog::Scope scope("MopacMcConnellResult::Compute",
        "compatibility wrapper; atoms=" + std::to_string(conf.AtomCount()) +
        " bonds=" + std::to_string(conf.ProteinRef().BondCount()));

    if (!conf.HasResult<McConnellResult>()) {
        OperationLog::Warn(
            "MopacMcConnellResult::Compute",
            "McConnellResult is not attached; BO-channel legacy fields may "
            "be unset. AttachResult dependency checks will reject this result "
            "in production.");
    }

    auto result_ptr = std::make_unique<MopacMcConnellResult>();
    result_ptr->conf_ = &conf;
    OperationLog::Info(LogCalcMcConnell, "MopacMcConnellResult::Compute",
        "no separate tensor computation; uses McConnellResult BO channel "
        "(Wiberg-weighted D(r)Qhat)");
    return result_ptr;
}


int MopacMcConnellResult::WriteFeatures(const ProteinConformation&,
                                         const std::string&) const {
    return 0;
}

}  // namespace nmr

#include "Session.h"

#include "AIMNet2Result.h"        // AIMNet2Model::Load
#include "OperationLog.h"
#include "RuntimeEnvironment.h"

#include <filesystem>

namespace fs = std::filesystem;

namespace nmr {


Session::Session() = default;
Session::~Session() = default;


Status Session::LoadFromToml() {
    // Subsystems still populate static state internally — the object
    // model's instance-ising of them is a separate refactor. Session
    // is the named startup step that orchestrates them in the correct
    // order.
    //
    // RuntimeEnvironment: machine paths from ~/.nmr_tools.toml.
    // OperationLog: channel mask + session-start log line.
    // CalculatorConfig: NOT loaded here — the caller owns the decision
    // of which params file to use (default in data/ vs JobSpec's
    // --config override), and CalculatorConfig::Load is idempotent.
    RuntimeEnvironment::Load();
    OperationLog::LoadChannelConfig();
    OperationLog::LogSessionStart();
    return kOk;
}


Status Session::LoadAimnet2Model(const std::string& path) {
    if (path.empty()) {
        last_error_ = "Session::LoadAimnet2Model: path is empty";
        return kSessionAimnet2LoadFailed;
    }
    if (!fs::exists(path)) {
        last_error_ =
            "Session::LoadAimnet2Model: model file not found: " + path;
        OperationLog::Error("Session::LoadAimnet2Model", last_error_);
        return kSessionAimnet2LoadFailed;
    }
    aimnet2_model_ = AIMNet2Model::Load(path);
    if (!aimnet2_model_) {
        last_error_ =
            "Session::LoadAimnet2Model: AIMNet2Model::Load returned null "
            "for " + path;
        OperationLog::Error("Session::LoadAimnet2Model", last_error_);
        return kSessionAimnet2LoadFailed;
    }
    OperationLog::Info(LogCalcOther, "Session::LoadAimnet2Model",
                       "loaded AIMNet2 model from " + path);
    return kOk;
}

}  // namespace nmr

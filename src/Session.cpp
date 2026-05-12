#include "Session.h"

#include "AIMNet2Result.h"        // AIMNet2Model::Load
#include "CategoryInfoProjection.h"
#include "LarsenHBondGrid.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"
#include "TripeptideDftTable.h"

#include <exception>
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

    // CategoryInfoProjection: one-shot setup of the output-side per-atom
    // categorical record (atom_nom.tbl). Inert when bmrb_atom_nom path is
    // not set in TOML; emission then runs but emits AMBER names as
    // fallback (provenance=MissLogged for every atom).
    CategoryInfoProjection::Config cic_cfg;
    cic_cfg.atom_nom_tbl = RuntimeEnvironment::BmrbAtomNom();
    CategoryInfoProjection::Configure(std::move(cic_cfg));
    return kOk;
}


Status Session::LoadTripeptideDftTable() {
    const std::string& dsn = RuntimeEnvironment::TensorCs15Dsn();
    if (dsn.empty()) {
        OperationLog::Info(LogCalcOther, "Session::LoadTripeptideDftTable",
            "no [databases].tensorcs15 DSN — table not loaded "
            "(TripeptideBackboneShieldingResult will skip).");
        return kOk;
    }
    try {
        tripeptide_dft_table_ =
            std::make_unique<TripeptideDftTable>(dsn);
    } catch (const std::exception& e) {
        last_error_ =
            "Session::LoadTripeptideDftTable: " + std::string(e.what());
        OperationLog::Error("Session::LoadTripeptideDftTable",
                            last_error_);
        tripeptide_dft_table_.reset();
        return kSessionTripeptideDbLoadFailed;
    }
    OperationLog::Info(LogCalcOther, "Session::LoadTripeptideDftTable",
                       "tensorcs15 connection open");
    return kOk;
}


Status Session::LoadLarsenHBondGrid() {
    const std::string& dir = RuntimeEnvironment::LarsenHBondGridDir();
    if (dir.empty()) {
        OperationLog::Info(LogCalcOther, "Session::LoadLarsenHBondGrid",
            "no larsen_hbond_grids path — grid not loaded "
            "(LarsenHBondShieldingResult will skip).");
        return kOk;
    }
    try {
        larsen_hbond_grid_ =
            std::make_unique<LarsenHBondGrid>(dir);
    } catch (const std::exception& e) {
        last_error_ =
            "Session::LoadLarsenHBondGrid: " + std::string(e.what());
        OperationLog::Error("Session::LoadLarsenHBondGrid",
                            last_error_);
        larsen_hbond_grid_.reset();
        return kSessionLarsenHBondGridLoadFailed;
    }
    OperationLog::Info(LogCalcOther, "Session::LoadLarsenHBondGrid",
                       "6 dense H-bond grids loaded from " + dir);
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

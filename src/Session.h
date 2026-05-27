#pragma once
//
// Session: the named entity that holds what the process brings into
// every run. One object, built in main(), passed by const reference
// to whichever RunXxx dispatches. Replaces the scattering of global
// singletons (RuntimeEnvironment, OperationLog channel mask,
// g_aimnet2_model) with one place the caller can point at.
//
// Why an object when the subsystems are currently static-populated
// classes: because the subsystems have verbs the library uses —
// config lookups, log emission, path validation, model loading —
// and naming the holder makes those verbs findable. The static
// populations continue to live inside each subsystem for now (the
// refactor to instance members is a separate, larger motion);
// Session's LoadFromToml orchestrates the startup statics it owns, and
// Session's accessors carry the resources that aren't statics (the
// loaded AIMNet2 model and optional lookup tables).
//

#include "errors.h"

#include <filesystem>
#include <memory>
#include <string>

namespace nmr {

struct AIMNet2Model;  // defined in AIMNet2Result.h; forward-declared
                      // here to avoid pulling torch into every header
                      // that sees Session.

class TripeptideDftTable;  // defined in TripeptideDftTable.h;
                           // forward-declared so consumers don't need
                           // libpq-fe.h transitively.

class LarsenHBondGrid;     // defined in LarsenHBondGrid.h;
                           // forward-declared so consumers don't need
                           // HighFive headers transitively.

class Session {
public:
    Session();
    ~Session();

    Session(const Session&) = delete;
    Session& operator=(const Session&) = delete;

    // Load ~/.nmr_tools.toml and initialise the startup subsystems:
    // RuntimeEnvironment (mopac path, ff14sb path, tmpdir, process GUID),
    // OperationLog channel mask, and output category projection. The
    // caller loads CalculatorConfig because it owns the default-vs-CLI
    // params-file choice.
    //
    // Returns kOk on success, or an ErrorCode with a matching
    // OperationLog::Error diagnostic emitted at the failure site.
    Status LoadFromToml();

    // Optional. If the caller wants AIMNet2 in the pipeline, load the
    // .jpt model once here; Session holds it for the whole process.
    // Pass the resolved path (CLI + fallback from TOML).
    Status LoadAimnet2Model(const std::string& path);

    // Optional. Open a libpq connection to the local tensorcs15 replica
    // and own the TripeptideDftTable for the Session lifetime. The DSN
    // (libpq kv-pair string) comes from RuntimeEnvironment::TensorCs15Dsn.
    // Empty DSN → kOk with table left null; OperationRunner skips the
    // tripeptide calculators when no table pointer is supplied.
    Status LoadTripeptideDftTable();

    // Optional. Load the 6 dense Larsen H-bond DFT grids (HDF5) and
    // own the LarsenHBondGrid for the Session lifetime. The directory
    // comes from RuntimeEnvironment::LarsenHBondGridDir. Empty dir →
    // kOk with grid left null; OperationRunner skips the Larsen H-bond
    // calculator when no grid pointer is supplied.
    Status LoadLarsenHBondGrid();

    // Accessors. AIMNet2Model pointer is the persistent resource;
    // RuntimeEnvironment and OperationLog continue to be read via their
    // existing static accessors at consumer sites — Session's role for
    // those is initialisation and lifetime discipline, not per-call routing.
    AIMNet2Model* Aimnet2Model() const { return aimnet2_model_.get(); }
    bool HasAimnet2Model() const { return aimnet2_model_ != nullptr; }

    // Tripeptide DFT table accessor. Const because calculators only
    // read; the table internally serialises libpq access.
    const TripeptideDftTable* TripeptideDftTablePtr() const {
        return tripeptide_dft_table_.get();
    }
    bool HasTripeptideDftTable() const {
        return tripeptide_dft_table_ != nullptr;
    }

    // Larsen H-bond grid accessor. Read-only after load.
    const LarsenHBondGrid* LarsenHBondGridPtr() const {
        return larsen_hbond_grid_.get();
    }
    bool HasLarsenHBondGrid() const {
        return larsen_hbond_grid_ != nullptr;
    }

    // Error string corresponding to the last non-ok status from one
    // of this Session's Load calls. Empty when status was ok.
    const std::string& LastError() const { return last_error_; }

private:
    std::unique_ptr<AIMNet2Model> aimnet2_model_;
    std::unique_ptr<TripeptideDftTable> tripeptide_dft_table_;
    std::unique_ptr<LarsenHBondGrid> larsen_hbond_grid_;
    std::string last_error_;
};

}  // namespace nmr

#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "TrajectoryResult.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "AIMNet2Result.h"          // AIMNet2Model (forward use)
#include "GromacsFrameHandler.h"
#include "ProteinConformation.h"
#include "OperationRunner.h"
#include "OperationLog.h"
#include "errors.h"

// GROMACS EDR reader (copied layout from former GromacsRunContext).
#include "gromacs/fileio/enxio.h"
#include "gromacs/trajectory/energyframe.h"

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <typeindex>
#include <typeinfo>

namespace nmr {


Trajectory::Trajectory(std::filesystem::path xtc_path,
                       std::filesystem::path tpr_path,
                       std::filesystem::path edr_path)
    : xtc_path_(std::move(xtc_path)),
      tpr_path_(std::move(tpr_path)),
      edr_path_(std::move(edr_path)) {
    // Preload EDR at construction; consumed per-frame in Phases 6 + 7
    // via EnergyAtTime.
    if (!edr_path_.empty()) {
        if (!LoadEdr(edr_path_)) {
            OperationLog::Warn("Trajectory",
                "EDR load failed for " + edr_path_.string() +
                " (continuing without energy)");
        }
    }
}


Trajectory::~Trajectory() = default;


// ── Run (eight phases) ───────────────────────────────────────────

Status Trajectory::Run(TrajectoryProtein& tp,
                       const RunConfiguration& config,
                       const Session& session,
                       std::vector<std::unique_ptr<TrajectoryResult>> extras,
                       std::filesystem::path output_dir) {
    if (state_ == State::Complete) {
        OperationLog::Error("Trajectory::Run",
            "called on already-completed run");
        return kCalculatorPipelineFailed;
    }
    if (state_ == State::Running) {
        OperationLog::Error("Trajectory::Run",
            "called on active run");
        return kCalculatorPipelineFailed;
    }
    state_ = State::Running;
    output_dir_ = std::move(output_dir);

    OperationLog::Info(LogCalcOther, "Trajectory::Run",
        "starting " + config.Name() + " on " + xtc_path_.string());

    // =========================================================
    // Phase 1: open handler
    // =========================================================
    // Mount XTC + build PBC fixer from TPR. No frame is read here;
    // the handler is a pure reader that advances frames on demand.

    handler_ = std::make_unique<GromacsFrameHandler>(tp);
    if (!handler_->Open(xtc_path_.string(), tpr_path_.string())) {
        OperationLog::Error("Trajectory::Run",
            "Phase 1 handler open failed: " + handler_->error());
        return kXtcOpenFailed;
    }

    // =========================================================
    // Phase 2: read frame 0 + seed
    // =========================================================
    // Seed the canonical conformation (conf0). After Seed, Protein is
    // finalized (bonds + rings detected) and TrajectoryAtoms are
    // allocated — factories in Phase 3 can rely on all of it.

    if (!handler_->ReadNextFrame()) {
        OperationLog::Error("Trajectory::Run",
            "Phase 2 failed to read first frame: " + handler_->error());
        return kFrameReadFailed;
    }
    tp.Seed(handler_->ProteinPositions(), handler_->Time());

    // =========================================================
    // Phase 3: attach TrajectoryResults
    // =========================================================
    // Configuration-declared factories first (their order is the
    // dispatch order), then caller-supplied extras.

    for (const auto& factory : config.TrajectoryResultFactories()) {
        auto result = factory(tp);
        const std::string name = result ? result->Name() : std::string("(null)");
        if (!tp.AttachResult(std::move(result))) {
            OperationLog::Error("Trajectory::Run",
                "Phase 3 attach rejected " + name);
            return kAttachRejectedSingleton;
        }
    }
    for (auto& extra : extras) {
        const std::string name = extra ? extra->Name() : std::string("(null)");
        if (!tp.AttachResult(std::move(extra))) {
            OperationLog::Error("Trajectory::Run",
                "Phase 3 attach rejected extra " + name);
            return kAttachRejectedSingleton;
        }
    }

    // =========================================================
    // Phase 4: validate deps + caller-supplied resources
    // =========================================================
    // Each attached TrajectoryResult declares type_index deps. A dep
    // is satisfied if another TrajectoryResult of that type is already
    // attached OR the type is in the config's required-
    // ConformationResult set (i.e. the per-frame pipeline will run it).
    // Mandatory caller-supplied resources (AIMNet2 model if
    // config.RequiresAimnet2()) come from Session.

    for (TrajectoryResult* r : tp.ResultsInAttachOrder()) {
        for (const std::type_index& dep : r->Dependencies()) {
            if (tp.AllResults().count(dep) > 0) continue;
            if (config.RequiresConformationResult(dep)) continue;
            OperationLog::Error("Trajectory::Run",
                std::string("Phase 4: ") + r->Name() +
                " declares unmet dependency");
            return kAttachDependencyUnmet;
        }
    }

    if (config.RequiresAimnet2() && !session.HasAimnet2Model()) {
        OperationLog::Error("Trajectory::Run",
            config.Name() +
            " requires AIMNet2 but session has no model loaded");
        return kConfigRequiresAimnet2;
    }

    // =========================================================
    // Phase 5: build per-frame RunOptions template
    // =========================================================
    // Base set once (config flags + tp's run-constant resources + the
    // session's AIMNet2 model). Per-frame fields (solvent, frame
    // energy) are pointed at env below; env is updated each frame.

    RunOptions base_opts = config.PerFrameRunOptions();
    base_opts.charge_source = tp.Charges();
    if (tp.HasBondedParams()) base_opts.bonded_params = &tp.BondedParams();
    base_opts.aimnet2_model = session.Aimnet2Model();

    // =========================================================
    // Phase 6: frame 0 — env + calc + dispatch + record
    // =========================================================

    env_.solvent            = handler_->Solvent();
    env_.current_energy     = EnergyAtTime(handler_->Time());
    env_.current_frame_idx  = 0;
    env_.current_frame_time = handler_->Time();

    RunOptions frame_opts = base_opts;
    frame_opts.solvent      = &env_.solvent;
    frame_opts.frame_energy = env_.current_energy;

    {
        auto& conf0 = tp.CanonicalConformation();
        RunResult rr = OperationRunner::Run(conf0, frame_opts);
        if (!rr.Ok()) {
            OperationLog::Error("Trajectory::Run",
                "frame 0 calculator pipeline failed: " + rr.error);
            return kCalculatorPipelineFailed;
        }
        tp.DispatchCompute(conf0, *this, /*frame_idx=*/0, handler_->Time());
        frame_times_.push_back(handler_->Time());
        frame_indices_.push_back(0);
        frame_count_ = 1;
    }

    // =========================================================
    // Phase 7: per-frame loop
    // =========================================================
    // Stride N: dispatch every N-th frame, skip the intervening N-1.
    // Each iteration is uniform: read → tick → update env → run
    // calculators → dispatch TrajectoryResults → record.

    const std::size_t stride = config.Stride();
    bool eof = false;
    while (!eof) {
        for (std::size_t s = 0; s + 1 < stride; ++s) {
            if (!handler_->Skip()) { eof = true; break; }
        }
        if (eof) break;

        if (!handler_->ReadNextFrame()) break;

        // Update env for this frame.
        env_.solvent            = handler_->Solvent();
        env_.current_energy     = EnergyAtTime(handler_->Time());
        env_.current_frame_idx  = handler_->Index();
        env_.current_frame_time = handler_->Time();

        // frame_opts pointers follow env; energy pointer resolves fresh.
        frame_opts.solvent      = &env_.solvent;
        frame_opts.frame_energy = env_.current_energy;

        auto conf = tp.TickConformation(handler_->ProteinPositions());
        RunResult rr = OperationRunner::Run(*conf, frame_opts);
        if (!rr.Ok()) {
            OperationLog::Error("Trajectory::Run",
                std::string("frame ") + std::to_string(handler_->Index()) +
                " calculator pipeline failed: " + rr.error);
            return kCalculatorPipelineFailed;
        }
        tp.DispatchCompute(*conf, *this, handler_->Index(), handler_->Time());

        frame_times_.push_back(handler_->Time());
        frame_indices_.push_back(handler_->Index());
        ++frame_count_;
    }

    if (!handler_->error().empty()) {
        OperationLog::Warn("Trajectory::Run",
            "frame loop ended with handler error: " + handler_->error());
    }

    // =========================================================
    // Phase 8: Finalize
    // =========================================================
    // Selections are emitted directly into traj.Selections() by
    // whichever TrajectoryResults push during Compute or Finalize.
    // No dynamic_cast sweep here — the bag is the source of truth.

    tp.FinalizeAllResults(*this);

    state_ = State::Complete;
    handler_.reset();

    OperationLog::Info(LogCalcOther, "Trajectory::Run",
        "complete: " + std::to_string(frame_count_) + " frames, " +
        std::to_string(selections_.Count()) + " selections");
    return kOk;
}


// ── Record accessors ─────────────────────────────────────────────

double Trajectory::TotalTimePs() const {
    if (frame_times_.empty()) return 0.0;
    return frame_times_.back() - frame_times_.front();
}


// ── EnergyAtTime ─────────────────────────────────────────────────
//
// Binary search into preloaded EDR by time. Returns the nearest
// frame; clamps at ends. Null if no EDR loaded.

const GromacsEnergy* Trajectory::EnergyAtTime(double time_ps) const {
    if (edr_frames_.empty()) return nullptr;

    auto it = std::lower_bound(
        edr_frames_.begin(), edr_frames_.end(), time_ps,
        [](const GromacsEnergy& ge, double t) { return ge.time_ps < t; });

    if (it == edr_frames_.end()) return &edr_frames_.back();
    if (it == edr_frames_.begin()) return &edr_frames_.front();

    auto prev = std::prev(it);
    if (std::fabs(it->time_ps - time_ps) < std::fabs(prev->time_ps - time_ps))
        return &(*it);
    return &(*prev);
}


// ── WriteH5: /trajectory/ group ──────────────────────────────────

void Trajectory::WriteH5(HighFive::File& file) const {
    auto grp = file.createGroup("/trajectory");

    // Source paths as attributes (provenance).
    grp.createAttribute("xtc_path", xtc_path_.string());
    grp.createAttribute("tpr_path", tpr_path_.string());
    grp.createAttribute("edr_path", edr_path_.string());
    grp.createAttribute("configuration", std::string(
        state_ == State::Complete ? "complete" : "incomplete"));

    // Frame metadata.
    auto frames = file.createGroup("/trajectory/frames");
    frames.createDataSet("time_ps", frame_times_);
    frames.createDataSet("original_index", frame_indices_);
    frames.createAttribute("n_frames", frame_count_);

    // Selections — one group per kind present.
    //
    // The bag holds all pushed records regardless of emitter; here we
    // walk the distinct kinds and emit a group per kind with its
    // records' frame_idx, time_ps, and reason arrays. Metadata stays
    // as a free-form adjunct — emitters that care about metadata
    // schema emit it themselves in their Result's own WriteH5Group.
    for (const std::type_index& kind : selections_.Kinds()) {
        auto records = selections_.ByKind(kind);
        if (records.empty()) continue;

        // Group path uses the mangled type name (compiler-dependent
        // but stable within a build). Human-readable names can be
        // added via an attribute if consumers need them.
        std::string group_path = std::string("/trajectory/selections/") +
                                 kind.name();
        auto grp = file.createGroup(group_path);

        std::vector<std::size_t> idx(records.size());
        std::vector<double> tps(records.size());
        std::vector<std::string> reason(records.size());
        for (std::size_t i = 0; i < records.size(); ++i) {
            idx[i]    = records[i]->frame_idx;
            tps[i]    = records[i]->time_ps;
            reason[i] = records[i]->reason;
        }
        grp.createDataSet("frame_idx", idx);
        grp.createDataSet("time_ps", tps);
        grp.createDataSet("reason", reason);
        grp.createAttribute("n_records", records.size());
    }
}


// ── LoadEdr (copied from the dissolved GromacsRunContext::LoadEdr) ─

bool Trajectory::LoadEdr(const std::filesystem::path& edr_path) {
    ener_file_t ef = open_enx(edr_path.string(), "r");
    if (!ef) return false;

    int nre = 0;
    gmx_enxnm_t* enms = nullptr;
    do_enxnms(ef, &nre, &enms);

    std::map<std::string, int> name_to_idx;
    for (int i = 0; i < nre; ++i) name_to_idx[enms[i].name] = i;

    auto idx = [&](const char* name) -> int {
        auto it = name_to_idx.find(name);
        return (it != name_to_idx.end()) ? it->second : -1;
    };

    int i_coul_sr = idx("Coulomb (SR)"), i_coul_recip = idx("Coul. recip.");
    int i_coul_14 = idx("Coulomb-14");
    int i_bond = idx("Bond"), i_angle = idx("Angle"), i_ub = idx("U-B");
    int i_proper = idx("Proper Dih."), i_improper = idx("Improper Dih.");
    int i_cmap = idx("CMAP Dih.");
    int i_lj_sr = idx("LJ (SR)"), i_lj_14 = idx("LJ-14");
    int i_dispcorr = idx("Disper. corr.");
    int i_potential = idx("Potential"), i_kinetic = idx("Kinetic En.");
    int i_total = idx("Total Energy"), i_enthalpy = idx("Enthalpy");
    int i_temperature = idx("Temperature"), i_pressure = idx("Pressure");
    int i_volume = idx("Volume"), i_density = idx("Density");
    int i_box_x = idx("Box-X"), i_box_y = idx("Box-Y"), i_box_z = idx("Box-Z");
    int i_vir[9] = {
        idx("Vir-XX"), idx("Vir-XY"), idx("Vir-XZ"),
        idx("Vir-YX"), idx("Vir-YY"), idx("Vir-YZ"),
        idx("Vir-ZX"), idx("Vir-ZY"), idx("Vir-ZZ") };
    int i_pres[9] = {
        idx("Pres-XX"), idx("Pres-XY"), idx("Pres-XZ"),
        idx("Pres-YX"), idx("Pres-YY"), idx("Pres-YZ"),
        idx("Pres-ZX"), idx("Pres-ZY"), idx("Pres-ZZ") };
    int i_T_prot = idx("T-Protein"), i_T_nonprot = idx("T-non-Protein");

    auto e = [](const t_enxframe& fr, int i) -> double {
        return (i >= 0) ? fr.ener[i].e : 0.0;
    };

    t_enxframe fr;
    init_enxframe(&fr);

    while (do_enx(ef, &fr)) {
        GromacsEnergy ge;
        ge.time_ps       = fr.t;
        ge.coulomb_sr    = e(fr, i_coul_sr);
        ge.coulomb_recip = e(fr, i_coul_recip);
        ge.coulomb_14    = e(fr, i_coul_14);
        ge.bond          = e(fr, i_bond);
        ge.angle         = e(fr, i_angle);
        ge.urey_bradley  = e(fr, i_ub);
        ge.proper_dih    = e(fr, i_proper);
        ge.improper_dih  = e(fr, i_improper);
        ge.cmap_dih      = e(fr, i_cmap);
        ge.lj_sr         = e(fr, i_lj_sr);
        ge.lj_14         = e(fr, i_lj_14);
        ge.disper_corr   = e(fr, i_dispcorr);
        ge.potential     = e(fr, i_potential);
        ge.kinetic       = e(fr, i_kinetic);
        ge.total_energy  = e(fr, i_total);
        ge.enthalpy      = e(fr, i_enthalpy);
        ge.temperature   = e(fr, i_temperature);
        ge.pressure      = e(fr, i_pressure);
        ge.volume        = e(fr, i_volume);
        ge.density       = e(fr, i_density);
        ge.box_x         = e(fr, i_box_x);
        ge.box_y         = e(fr, i_box_y);
        ge.box_z         = e(fr, i_box_z);
        for (int k = 0; k < 9; ++k) ge.vir[k]  = e(fr, i_vir[k]);
        for (int k = 0; k < 9; ++k) ge.pres[k] = e(fr, i_pres[k]);
        ge.T_protein     = e(fr, i_T_prot);
        ge.T_non_protein = e(fr, i_T_nonprot);
        edr_frames_.push_back(ge);
    }

    free_enxframe(&fr);
    free_enxnms(nre, enms);
    close_enx(ef);

    return !edr_frames_.empty();
}

}  // namespace nmr

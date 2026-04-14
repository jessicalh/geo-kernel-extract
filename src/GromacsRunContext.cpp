#include "GromacsRunContext.h"
#include "FullSystemReader.h"
#include "OperationLog.h"

// GROMACS EDR reading
#include "gromacs/fileio/enxio.h"
#include "gromacs/trajectory/energyframe.h"

#include <cmath>
#include <cstring>
#include <map>
#include <string>

namespace nmr {

// ── Build ───────────────────────────────────────────────────────

bool GromacsRunContext::Build(const FullSystemReader& reader,
                              const std::string& edr_path) {
    // Bonded interaction parameters from the reader's single TPR parse.
    // ReadTopology() already extracted these — no re-read.
    bonded_params = reader.BondedParams();

    // EDR energy frames (optional)
    if (!edr_path.empty()) {
        if (!LoadEdr(edr_path)) {
            OperationLog::Warn("GromacsRunContext",
                "EDR load failed: " + error + " (continuing without energy)");
            error.clear();
        }
    }

    OperationLog::Info(LogCalcOther, "GromacsRunContext::Build",
        std::to_string(bonded_params.interactions.size()) + " bonded interactions" +
        (HasEdr() ? ", " + std::to_string(edr_frames.size()) + " EDR frames" : ""));

    return true;
}


// ── AdvanceFrame ────────────────────────────────────────────────

void GromacsRunContext::AdvanceFrame(size_t idx, double time_ps) {
    frame_index = idx;
    frame_time_ps = time_ps;
    current_energy = EnergyAtTime(time_ps);
}


// ── EnergyAtTime ────────────────────────────────────────────────

const GromacsEnergy* GromacsRunContext::EnergyAtTime(double time_ps) const {
    if (edr_frames.empty()) return nullptr;

    auto it = std::lower_bound(
        edr_frames.begin(), edr_frames.end(), time_ps,
        [](const GromacsEnergy& ge, double t) { return ge.time_ps < t; });

    if (it == edr_frames.end()) return &edr_frames.back();
    if (it == edr_frames.begin()) return &edr_frames.front();

    auto prev = std::prev(it);
    if (std::fabs(it->time_ps - time_ps) < std::fabs(prev->time_ps - time_ps))
        return &(*it);
    return &(*prev);
}


// ── EDR bulk read ───────────────────────────────────────────────

bool GromacsRunContext::LoadEdr(const std::string& edr_path) {
    ener_file_t ef = open_enx(edr_path, "r");
    if (!ef) {
        error = "cannot open EDR: " + edr_path;
        return false;
    }

    int nre = 0;
    gmx_enxnm_t* enms = nullptr;
    do_enxnms(ef, &nre, &enms);

    std::map<std::string, int> name_to_idx;
    for (int i = 0; i < nre; ++i)
        name_to_idx[enms[i].name] = i;

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
        edr_frames.push_back(ge);
    }

    free_enxframe(&fr);
    free_enxnms(nre, enms);
    close_enx(ef);

    return !edr_frames.empty();
}

}  // namespace nmr

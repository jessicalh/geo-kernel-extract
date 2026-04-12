#include "GromacsEnergyResult.h"
#include "OperationLog.h"
#include "NpyWriter.h"

#include "gromacs/fileio/enxio.h"
#include "gromacs/trajectory/energyframe.h"

#include <cmath>
#include <cstring>
#include <map>
#include <string>

namespace nmr {

// ── .edr reading ─────────────────────────────────────────────────

// Read a .edr file and return the frame nearest to target_time_ps.
// Returns false if file cannot be read or is empty.
static bool ReadEdrFrame(const std::string& edr_path,
                         double target_time_ps,
                         GromacsEnergy& out) {

    ener_file_t ef = open_enx(edr_path, "r");
    if (!ef) return false;

    // Read energy term names
    int nre = 0;
    gmx_enxnm_t* enms = nullptr;
    do_enxnms(ef, &nre, &enms);

    // Build name → index map
    std::map<std::string, int> name_to_idx;
    for (int i = 0; i < nre; ++i) {
        name_to_idx[enms[i].name] = i;
    }

    // Find indices for terms we want
    // GROMACS .edr stores names with spaces (e.g. "Coulomb (SR)")
    // but gmx energy displays with dashes. Try the actual stored names.
    auto idx = [&](const std::string& name) -> int {
        auto it = name_to_idx.find(name);
        return (it != name_to_idx.end()) ? it->second : -1;
    };

    int i_coul_sr    = idx("Coulomb (SR)");
    int i_coul_recip = idx("Coul. recip.");
    int i_coul_14    = idx("Coulomb-14");
    int i_lj_sr      = idx("LJ (SR)");
    int i_potential   = idx("Potential");
    int i_temperature = idx("Temperature");
    int i_pressure    = idx("Pressure");
    int i_volume      = idx("Volume");

    // Log which terms were found
    OperationLog::Info(LogCalcOther, "ReadEdrFrame",
        "terms: Coul_SR=" + std::to_string(i_coul_sr) +
        " Coul_recip=" + std::to_string(i_coul_recip) +
        " LJ_SR=" + std::to_string(i_lj_sr) +
        " Potential=" + std::to_string(i_potential));

    // Scan frames for closest match to target time
    t_enxframe fr;
    init_enxframe(&fr);

    bool found = false;
    double best_dt = 1e30;

    while (do_enx(ef, &fr)) {
        double dt = std::fabs(fr.t - target_time_ps);
        if (dt < best_dt) {
            best_dt = dt;
            out.time_ps       = fr.t;
            out.coulomb_sr    = (i_coul_sr    >= 0) ? fr.ener[i_coul_sr].e    : 0.0;
            out.coulomb_recip = (i_coul_recip >= 0) ? fr.ener[i_coul_recip].e : 0.0;
            out.coulomb_14    = (i_coul_14    >= 0) ? fr.ener[i_coul_14].e    : 0.0;
            out.lj_sr         = (i_lj_sr      >= 0) ? fr.ener[i_lj_sr].e      : 0.0;
            out.potential     = (i_potential   >= 0) ? fr.ener[i_potential].e   : 0.0;
            out.temperature   = (i_temperature >= 0) ? fr.ener[i_temperature].e : 0.0;
            out.pressure      = (i_pressure    >= 0) ? fr.ener[i_pressure].e    : 0.0;
            out.volume        = (i_volume      >= 0) ? fr.ener[i_volume].e      : 0.0;
            found = true;
        }
        // Once we've passed the target time by more than the best match,
        // stop scanning (frames are chronological)
        if (fr.t > target_time_ps && dt > best_dt) break;
    }

    free_enxframe(&fr);
    free_enxnms(nre, enms);
    close_enx(ef);

    return found;
}


// ── ConformationResult interface ─────────────────────────────────

std::unique_ptr<GromacsEnergyResult> GromacsEnergyResult::Compute(
        ProteinConformation& conf,
        const std::string& edr_path,
        double target_time_ps) {

    OperationLog::Scope scope("GromacsEnergyResult::Compute",
        "edr=" + edr_path + " t=" + std::to_string(target_time_ps) + "ps");

    if (edr_path.empty()) {
        OperationLog::Error("GromacsEnergyResult", "no .edr path provided");
        return nullptr;
    }

    auto result = std::make_unique<GromacsEnergyResult>();

    if (!ReadEdrFrame(edr_path, target_time_ps, result->energy_)) {
        OperationLog::Error("GromacsEnergyResult",
            "failed to read .edr or no matching frame for t=" +
            std::to_string(target_time_ps) + "ps");
        return nullptr;
    }

    double dt = std::fabs(result->energy_.time_ps - target_time_ps);
    OperationLog::Info(LogCalcOther, "GromacsEnergyResult",
        "matched t=" + std::to_string(result->energy_.time_ps) +
        "ps (dt=" + std::to_string(dt) + "ps)" +
        " Coul_total=" + std::to_string(result->energy_.CoulombTotal()) +
        " potential=" + std::to_string(result->energy_.potential) +
        " T=" + std::to_string(result->energy_.temperature) + "K");

    return result;
}


int GromacsEnergyResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    // Write as a (1, 9) array: one row of per-frame scalars.
    // Columns: time_ps, coulomb_sr, coulomb_recip, coulomb_14,
    //          lj_sr, potential, temperature, pressure, volume
    const int COLS = 9;
    std::vector<double> row = {
        energy_.time_ps,
        energy_.coulomb_sr,
        energy_.coulomb_recip,
        energy_.coulomb_14,
        energy_.lj_sr,
        energy_.potential,
        energy_.temperature,
        energy_.pressure,
        energy_.volume
    };

    NpyWriter::WriteFloat64(output_dir + "/gromacs_energy.npy", row.data(), 1, COLS);
    return 1;
}

}  // namespace nmr

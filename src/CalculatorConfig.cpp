#include "CalculatorConfig.h"
#include "OperationLog.h"

#include <fstream>
#include <filesystem>
#include <cstdio>
#include <cstdlib>

namespace fs = std::filesystem;

namespace nmr {

// Static members
std::unordered_map<std::string, CalculatorConfig::ParamEntry> CalculatorConfig::defaults_;
std::unordered_map<std::string, double> CalculatorConfig::overrides_;
std::unordered_map<std::string, std::string> CalculatorConfig::string_overrides_;
bool CalculatorConfig::loaded_ = false;
bool CalculatorConfig::defaults_initialised_ = false;


void CalculatorConfig::InitDefaults() {
    if (defaults_initialised_) return;

    auto add = [](const char* key, double val, const char* unit, const char* desc) {
        defaults_[key] = {val, unit, desc};
    };

    // Ring current intensities (nA/T)
    // Giessner-Prettre & Pullman 1969; TRP indole perimeter estimated
    add("phe_benzene_ring_current_intensity",          -12.0,  "nA/T", "PHE benzene ring current");
    add("tyr_phenol_ring_current_intensity",           -11.28, "nA/T", "TYR phenol ring current");
    add("trp_benzene_ring_current_intensity",          -12.48, "nA/T", "TRP benzene ring current");
    add("trp_pyrrole_ring_current_intensity",          -6.72,  "nA/T", "TRP pyrrole ring current");
    add("his_imidazole_ring_current_intensity",        -5.16,  "nA/T", "HIS imidazole ring current");
    add("hid_imidazole_ring_current_intensity",        -5.16,  "nA/T", "HID imidazole ring current");
    add("hie_imidazole_ring_current_intensity",        -5.16,  "nA/T", "HIE imidazole ring current");
    add("trp_indole_perimeter_ring_current_intensity", -19.2,  "nA/T", "TRP indole perimeter ring current");

    // Johnson-Bovey pi-electron lobe offsets above/below ring plane (Angstroms)
    // Johnson & Bovey 1958; TRP indole perimeter estimated
    add("phe_benzene_jb_lobe_offset",          0.64, "A", "PHE benzene JB lobe offset");
    add("tyr_phenol_jb_lobe_offset",           0.64, "A", "TYR phenol JB lobe offset");
    add("trp_benzene_jb_lobe_offset",          0.64, "A", "TRP benzene JB lobe offset");
    add("trp_pyrrole_jb_lobe_offset",          0.52, "A", "TRP pyrrole JB lobe offset");
    add("his_imidazole_jb_lobe_offset",        0.50, "A", "HIS imidazole JB lobe offset");
    add("hid_imidazole_jb_lobe_offset",        0.50, "A", "HID imidazole JB lobe offset");
    add("hie_imidazole_jb_lobe_offset",        0.50, "A", "HIE imidazole JB lobe offset");
    add("trp_indole_perimeter_jb_lobe_offset", 0.60, "A", "TRP indole perimeter JB lobe offset");

    // Spatial cutoffs — define the reach of each geometric kernel (Angstroms)
    add("ring_current_spatial_cutoff",              15.0, "A", "ring current spatial cutoff");
    add("mcconnell_bond_anisotropy_cutoff",         10.0, "A", "McConnell bond anisotropy cutoff");
    add("mopac_mcconnell_bond_anisotropy_cutoff",   10.0, "A", "MOPAC McConnell bond anisotropy cutoff");
    add("coulomb_efield_cutoff",                    20.0, "A", "Coulomb E-field spatial cutoff");
    add("hbond_dipolar_max_distance",               50.0, "A", "H-bond dipolar maximum distance");
    add("hbond_counting_radius",                     3.5, "A", "H-bond counting radius");
    add("dispersion_vertex_distance_cutoff",         5.0, "A", "dispersion vertex distance cutoff");
    add("dispersion_switching_onset_distance",       4.3, "A", "dispersion switching function onset");
    add("singularity_guard_distance",                0.1, "A", "singularity guard distance");

    // Ring proximity counting shells (Angstroms)
    add("ring_proximity_shell_1",  3.0, "A", "ring proximity shell 1");
    add("ring_proximity_shell_2",  5.0, "A", "ring proximity shell 2");
    add("ring_proximity_shell_3",  8.0, "A", "ring proximity shell 3");
    add("ring_proximity_shell_4", 12.0, "A", "ring proximity shell 4");

    // Haigh-Mallion adaptive quadrature subdivision thresholds (Angstroms)
    add("haigh_mallion_subdivision_threshold_l1", 2.0, "A", "HM level-1 subdivision threshold");
    add("haigh_mallion_subdivision_threshold_l2", 1.0, "A", "HM level-2 subdivision threshold");

    // Guard thresholds — algorithmic boundaries
    add("efield_magnitude_sanity_clamp",        100.0, "V/A",     "E-field magnitude clamp");
    add("hbond_sequential_exclusion_residues",    2.0, "residues", "H-bond sequence exclusion");
    add("near_field_exclusion_ratio",             0.5, "",         "near-field exclusion ratio");

    // AIMNet2 neural network charge calculator
    add("aimnet2_cutoff_lr",                       15.0, "A",    "AIMNet2 long-range DSF Coulomb cutoff");
    add("aimnet2_max_nb",                         128.0, "",     "AIMNet2 max short-range neighbours");
    add("aimnet2_max_nb_lr",                     4096.0, "",     "AIMNet2 max long-range neighbours");
    add("aimnet2_coulomb_efg_cutoff",              20.0, "A",    "AIMNet2 Coulomb EFG cutoff");
    // charge sensitivity: no calculator params — computed by
    // GromacsFrameHandler as per-atom charge variance across ensemble

    // SASA (Shrake-Rupley per-atom solvent-accessible surface area)
    add("sasa_probe_radius",                        1.4, "A",    "SASA water probe radius (Bondi)");
    add("sasa_n_points",                           92.0, "",     "SASA Fibonacci sphere point count");

    // WaterFieldResult — explicit solvent E-field and EFG
    add("water_efield_cutoff",                     15.0, "A",    "water E-field summation cutoff (oxygen distance)");
    add("water_first_shell_cutoff",                 3.5, "A",    "first hydration shell boundary (WaterField + HydrationShell)");
    add("water_second_shell_cutoff",                5.5, "A",    "second hydration shell boundary");

    // HydrationShellResult — hydration shell geometry
    add("hydration_ion_cutoff",                    20.0, "A",    "nearest-ion search distance");

    // EEQ — Extended Electronegativity Equilibration (Caldeweyher et al. 2019)
    add("eeq_total_charge",                         0.0, "e",    "EEQ net system charge (0 = neutral protein)");
    add("eeq_cn_steepness",                         7.5, "",     "EEQ coordination number error function steepness");
    add("eeq_cn_cutoff",                           25.0, "A",    "EEQ coordination number pair cutoff distance");
    add("eeq_charge_clamp",                         2.0, "e",    "EEQ per-atom charge magnitude clamp");

    // Numerical noise floors — below these the value is treated as zero
    add("near_zero_vector_norm_threshold",     1e-10, "",    "near-zero vector norm threshold");
    add("coulomb_charge_noise_floor",          1e-15, "",    "Coulomb charge noise floor");
    add("mopac_bond_order_noise_floor",        1e-6,  "",    "MOPAC bond order noise floor");
    add("biot_savart_wire_endpoint_guard",     1e-25, "m",   "Biot-Savart wire endpoint guard (SI)");
    add("biot_savart_wire_axis_guard",         1e-70, "m^2", "Biot-Savart wire on-axis guard (SI)");
    add("haigh_mallion_triangle_area_guard",   1e-20, "",    "Haigh-Mallion triangle area guard");
    add("dispersion_switching_noise_floor",    1e-15, "",    "dispersion switching noise floor");

    defaults_initialised_ = true;
}


void CalculatorConfig::Load(const std::string& path) {
    InitDefaults();

    if (path.empty()) {
        OperationLog::Info("CalculatorConfig::Load",
            "no config path — using compiled defaults for all "
            + std::to_string(defaults_.size()) + " parameters");
        loaded_ = true;
        return;
    }

    if (!fs::exists(path)) {
        OperationLog::Warn("CalculatorConfig::Load",
            "config file not found: " + path + " — using compiled defaults");
        loaded_ = true;
        return;
    }

    std::ifstream in(path);
    std::string line;
    int line_num = 0;

    while (std::getline(in, line)) {
        ++line_num;

        // Strip comments
        auto pos = line.find('#');
        if (pos != std::string::npos) line = line.substr(0, pos);

        if (line.find('=') == std::string::npos) continue;

        auto eq = line.find('=');
        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);

        // Trim whitespace
        auto trim = [](std::string& s) {
            while (!s.empty() && (s.front() == ' ' || s.front() == '\t'))
                s.erase(s.begin());
            while (!s.empty() && (s.back() == ' ' || s.back() == '\t'))
                s.pop_back();
        };
        trim(key);
        trim(val);

        if (key.empty() || val.empty()) continue;

        char* end = nullptr;
        double d = std::strtod(val.c_str(), &end);
        if (end == val.c_str()) {
            // Not a number — store as string (strip surrounding quotes if present)
            std::string s = val;
            if (s.size() >= 2 && s.front() == '"' && s.back() == '"')
                s = s.substr(1, s.size() - 2);
            string_overrides_[key] = s;
            continue;
        }

        overrides_[key] = d;
    }

    OperationLog::Info("CalculatorConfig::Load",
        "read " + path + ": " + std::to_string(overrides_.size())
        + " overrides, " + std::to_string(defaults_.size()) + " defaults");

    loaded_ = true;
}


double CalculatorConfig::Get(const std::string& key) {
    if (!defaults_initialised_) InitDefaults();

    // Check override first
    auto ov = overrides_.find(key);
    if (ov != overrides_.end()) return ov->second;

    // Then default
    auto df = defaults_.find(key);
    if (df != defaults_.end()) return df->second.value;

    // Unknown key = programming error
    fprintf(stderr,
        "FATAL: CalculatorConfig::Get(\"%s\") — unknown parameter key.\n",
        key.c_str());
    std::abort();
    return 0.0;  // unreachable
}


std::string CalculatorConfig::GetString(const std::string& key,
                                        const std::string& defaultVal) {
    auto it = string_overrides_.find(key);
    if (it != string_overrides_.end()) return it->second;
    return defaultVal;
}


std::vector<std::string> CalculatorConfig::Validate() {
    if (!defaults_initialised_) InitDefaults();

    std::vector<std::string> unknown_keys;

    // Log every parameter
    for (const auto& [key, entry] : defaults_) {
        auto ov = overrides_.find(key);
        const char* source = (ov != overrides_.end()) ? "toml" : "default";
        double val = (ov != overrides_.end()) ? ov->second : entry.value;

        char buf[512];
        snprintf(buf, sizeof(buf), "%-50s = %12g  %-8s [%s] %s",
                 key.c_str(), val, entry.unit, source, entry.description);
        OperationLog::Info("CalculatorConfig", buf);
    }

    // Check for unknown keys in overrides
    for (const auto& [key, val] : overrides_) {
        if (defaults_.find(key) == defaults_.end()) {
            unknown_keys.push_back(key);
            OperationLog::Warn("CalculatorConfig",
                "unknown key in TOML: " + key + " (typo?)");
        }
    }

    return unknown_keys;
}

}  // namespace nmr

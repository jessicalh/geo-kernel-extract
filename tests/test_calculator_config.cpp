#include <gtest/gtest.h>
#include "CalculatorConfig.h"
#include "PhysicalConstants.h"

using namespace nmr;


// --- Defaults match current hardcoded values ---

TEST(CalculatorConfig, DefaultsMatchRingIntensities) {
    // Ring.h: PheBenzeneRing::Intensity() etc.
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("phe_benzene_ring_current_intensity"),          -12.0);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("tyr_phenol_ring_current_intensity"),           -11.28);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("trp_benzene_ring_current_intensity"),          -12.48);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("trp_pyrrole_ring_current_intensity"),          -6.72);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("his_imidazole_ring_current_intensity"),        -5.16);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("hid_imidazole_ring_current_intensity"),        -5.16);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("hie_imidazole_ring_current_intensity"),        -5.16);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("trp_indole_perimeter_ring_current_intensity"), -19.2);
}

TEST(CalculatorConfig, DefaultsMatchRingLobeOffsets) {
    // Ring.h: PheBenzeneRing::JBLobeOffset() etc.
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("phe_benzene_jb_lobe_offset"),          0.64);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("tyr_phenol_jb_lobe_offset"),           0.64);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("trp_benzene_jb_lobe_offset"),          0.64);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("trp_pyrrole_jb_lobe_offset"),          0.52);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("his_imidazole_jb_lobe_offset"),        0.50);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("hid_imidazole_jb_lobe_offset"),        0.50);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("hie_imidazole_jb_lobe_offset"),        0.50);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("trp_indole_perimeter_jb_lobe_offset"), 0.60);
}

TEST(CalculatorConfig, DefaultsMatchSpatialCutoffs) {
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("ring_current_spatial_cutoff"),            RING_CALC_CUTOFF);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("singularity_guard_distance"),             MIN_DISTANCE);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("hbond_counting_radius"),                  HBOND_COUNT_RADIUS);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("hbond_dipolar_max_distance"),             HBOND_MAX_DIST);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("efield_magnitude_sanity_clamp"),          APBS_SANITY_LIMIT);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("hbond_sequential_exclusion_residues"),
                     static_cast<double>(SEQUENTIAL_EXCLUSION_THRESHOLD));
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("near_zero_vector_norm_threshold"),        NEAR_ZERO_NORM);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("ring_proximity_shell_1"),                 RING_COUNT_SHELL_1);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("ring_proximity_shell_2"),                 RING_COUNT_SHELL_2);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("ring_proximity_shell_3"),                 RING_COUNT_SHELL_3);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("ring_proximity_shell_4"),                 RING_COUNT_SHELL_4);
}

TEST(CalculatorConfig, DefaultsMatchCalculatorConstants) {
    // From calculator .h / .cpp files
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("mcconnell_bond_anisotropy_cutoff"),        10.0);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("mopac_mcconnell_bond_anisotropy_cutoff"),  10.0);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("dispersion_vertex_distance_cutoff"),        5.0);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("dispersion_switching_onset_distance"),      4.3);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("haigh_mallion_subdivision_threshold_l1"),   2.0);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("haigh_mallion_subdivision_threshold_l2"),   1.0);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("near_field_exclusion_ratio"),               0.5);
}

TEST(CalculatorConfig, DefaultsMatchNoiseFloors) {
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("coulomb_charge_noise_floor"),          1e-15);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("mopac_bond_order_noise_floor"),        1e-6);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("biot_savart_wire_endpoint_guard"),     1e-25);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("biot_savart_wire_axis_guard"),         1e-70);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("haigh_mallion_triangle_area_guard"),   1e-20);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("dispersion_switching_noise_floor"),    1e-15);
}


// --- TOML loading round-trip ---

TEST(CalculatorConfig, TomlRoundTrip) {
    std::string toml_path = std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml";
    CalculatorConfig::Load(toml_path);

    // After loading TOML with default values, Get() returns the same numbers
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("phe_benzene_ring_current_intensity"), -12.0);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("trp_indole_perimeter_jb_lobe_offset"), 0.60);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("ring_current_spatial_cutoff"), 15.0);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("near_field_exclusion_ratio"), 0.5);
    EXPECT_DOUBLE_EQ(CalculatorConfig::Get("biot_savart_wire_axis_guard"), 1e-70);

    // Validate should find no unknown keys
    auto unknown = CalculatorConfig::Validate();
    EXPECT_TRUE(unknown.empty());
}

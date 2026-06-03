// PerAtomSubstrate -- #58 lean all-atom aggregate substrate.
//
// The emitter is a carrier on RelationshipEngine::RunTraversal's typed
// overload. It writes one row per (DFT frame slot, atom) and never emits a
// default source CSV.

#pragma once

#include "AnalysisBody.h"
#include "ExtractionSupport.h"

#include <QMap>
#include <QString>
#include <QStringList>

#include <array>
#include <cstddef>
#include <cstdint>

namespace h5reader::rediscover {

struct PerAtomSubstrateConfig {
    double ring_cutoff_A = 8.0;
    double bond_cutoff_A = 10.0;
    double charge_cutoff_A = 6.0;
    double mc_near_field_ratio = 0.5;
    bool emit_embedding = true;
    std::size_t embedding_dims = 256;

    // Bounded named query outputs. These are deliberately small defaults, not a
    // resident pair dump.
    std::size_t query_frame_slots = 1;
    std::size_t reader_pair_atom = 0;
    std::size_t reader_pair_frame_slots = 3;
    int top_k = 3;
};

struct PerAtomChannelAudit {
    std::size_t present = 0;
    double min = 0.0;
    double max = 0.0;
    bool has_range = false;
};

struct PerAtomSubstrateStats {
    std::size_t atom_count = 0;
    std::size_t dft_rows = 0;
    std::size_t rows = 0;
    std::size_t dft_present = 0;
    std::size_t ring_present = 0;
    std::size_t charge_present = 0;
    std::size_t mc_lit_valid_present = 0;
    std::size_t ff14sb_field_present = 0;
    std::size_t apbs_efield_present = 0;
    std::size_t apbs_efg_present = 0;
    std::size_t aimnet2_charge_present = 0;
    std::size_t aimnet2_crg_present = 0;
    std::size_t aimnet2_embedding_present = 0;
    std::size_t ff14sb_charge_present = 0;
    std::size_t mopac_welford_mean_charge_present = 0;
    std::size_t charge_complete = 0;
    std::size_t mopac_coulomb_shielding_present = 0;
    std::size_t mopac_mc_shielding_present = 0;
    std::size_t hbond_shielding_present = 0;
    std::size_t hbond_count_present = 0;
    std::size_t hbond_geometry_present = 0;
    std::size_t pi_quadrupole_present = 0;
    std::size_t dispersion_present = 0;
    std::size_t hm_shielding_present = 0;
    std::size_t ringchi_shielding_present = 0;
    std::size_t water_field_present = 0;
    std::size_t hydration_shell_present = 0;
    std::size_t sasa_present = 0;
    std::size_t sasa_normal_present = 0;
    std::size_t eeq_charge_present = 0;
    std::size_t eeq_coordination_number_present = 0;
    std::size_t pair_query_rows = 0;
    std::size_t top_source_query_rows = 0;
    std::size_t dominance_query_rows = 0;
    std::size_t reader_pair_query_rows = 0;
    QMap<QString, PerAtomChannelAudit> new_channel_audit;
    QStringList absent_new_channel_slabs;
};

constexpr std::size_t kPerAtomClassicalCols = 89;
constexpr std::size_t kPerAtomConditioningCols = 26;
constexpr std::size_t kPerAtomDriverMagnitudeCols = 9;
constexpr std::size_t kPerAtomBackboneAuditCols = 14;
constexpr std::size_t kPerAtomTargetDecompositionCols = 21;
constexpr std::size_t kPerAtomRingPathCols = 226;

PerAtomSubstrateStats RunPerAtomSubstrateEmit(const Body& body,
                                              const QString& outDir,
                                              const PerAtomSubstrateConfig& config,
                                              const DftFrameAlignment& alignment);

QStringList PerAtomSubstrateSidecars(const PerAtomSubstrateConfig& config);

QMap<QString, std::size_t> PerAtomSubstrateFeatureSupport(const PerAtomSubstrateStats& stats);

}  // namespace h5reader::rediscover

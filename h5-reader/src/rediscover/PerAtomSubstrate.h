// PerAtomSubstrate -- #58 lean all-atom aggregate substrate.
//
// The emitter is a carrier on RelationshipEngine::RunTraversal's typed
// overload. It writes one row per (DFT frame slot, atom) and never emits a
// default source CSV.

#pragma once

#include "AnalysisBody.h"
#include "ExtractionSupport.h"

#include <QMap>
#include <QJsonObject>
#include <QString>
#include <QStringList>

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

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

struct PerAtomIsolationScalars {
    double gap_to_2nd_ring_r = std::numeric_limits<double>::quiet_NaN();
    double gap_to_2nd_charge_r = std::numeric_limits<double>::quiet_NaN();
    double gap_to_2nd_bond_r = std::numeric_limits<double>::quiet_NaN();
    double gap_to_2nd_field_r = std::numeric_limits<double>::quiet_NaN();
    double gap_to_2nd_hbond_r = std::numeric_limits<double>::quiet_NaN();
    double dominant_fraction_ring = std::numeric_limits<double>::quiet_NaN();
    double dominant_fraction_charge = std::numeric_limits<double>::quiet_NaN();
    double dominant_fraction_mc = std::numeric_limits<double>::quiet_NaN();
    double dominant_fraction_field = std::numeric_limits<double>::quiet_NaN();
    double dominant_fraction_hbond = std::numeric_limits<double>::quiet_NaN();
};

struct PairContribution {
    QString mechanism;
    QString source_kind;
    int32_t source_id = -1;
    int32_t source_atom_index = -1;
    int32_t source_cloud_index = -1;
    int source_category_ord = -1;
    int pointer_flags = 0;
    Vec3 disp = Vec3::Zero();
    double r = std::numeric_limits<double>::quiet_NaN();
    double inv_r3 = std::numeric_limits<double>::quiet_NaN();
    double cos_theta = std::numeric_limits<double>::quiet_NaN();
    double dipolar = std::numeric_limits<double>::quiet_NaN();
    double kernel_T0 = std::numeric_limits<double>::quiet_NaN();
    std::array<double, 5> kernel_T2 = {};
    double contribution = std::numeric_limits<double>::quiet_NaN();
};

enum PointerFlags : int {
    PresentFlag = 1 << 0,
    SelfOrBondedFlag = 1 << 1,
    ProducerValidFlag = 1 << 2,
    NearFieldFlag = 1 << 3,
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
    std::size_t dominance_build4_query_rows = 0;
    std::size_t reader_pair_query_rows = 0;
    QMap<QString, std::size_t> hunter_candidate_counts;
    bool hunter_anti_circular_assertion = false;
    QJsonObject partition_bin_manifest;
    QJsonObject dominance_bin_manifest;
    QMap<QString, PerAtomChannelAudit> new_channel_audit;
    QStringList absent_new_channel_slabs;
};

constexpr std::size_t kPerAtomClassicalCols = 89;
constexpr std::size_t kPerAtomConditioningCols = 32;
constexpr std::size_t kPerAtomDominanceCols = 10;
constexpr std::size_t kPerAtomDriverMagnitudeCols = 9;
constexpr std::size_t kPerAtomPartitionBinCols = 25;
constexpr std::size_t kPerAtomDominanceBinCols = 5;
constexpr std::size_t kPerAtomBackboneAuditCols = 14;
constexpr std::size_t kPerAtomTargetDecompositionCols = 21;
constexpr std::size_t kPerAtomRingPathCols = 226;
constexpr std::size_t kPerAtomMethodPathCols = 111;
constexpr std::size_t kPerAtomHbondConditioningCols = 73;

PerAtomSubstrateStats RunPerAtomSubstrateEmit(const Body& body,
                                              const QString& outDir,
                                              const PerAtomSubstrateConfig& config,
                                              const DftFrameAlignment& alignment);

QStringList PerAtomSubstrateSidecars(const PerAtomSubstrateConfig& config);

QMap<QString, std::size_t> PerAtomSubstrateFeatureSupport(const PerAtomSubstrateStats& stats);

std::vector<PairContribution> PerAtomRowPairContributions(const Body& body,
                                                          std::size_t atom,
                                                          std::size_t row,
                                                          const PerAtomSubstrateConfig& cfg,
                                                          const LocalFrame& frame);

PerAtomIsolationScalars PerAtomIsolationScalarsForRow(const Body& body,
                                                      std::size_t atom,
                                                      std::size_t row,
                                                      const PerAtomSubstrateConfig& cfg);

}  // namespace h5reader::rediscover

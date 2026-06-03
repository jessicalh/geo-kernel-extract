// CaseHunter -- compact input-side case selection for partition/equation strips.

#pragma once

#include "AnalysisBody.h"
#include "PerAtomSubstrate.h"

#include <QMap>
#include <QString>

#include <cstddef>

namespace h5reader::rediscover {

struct CaseHunterConfig {
    std::size_t window_before = 3;
    std::size_t window_after = 3;
    std::size_t frame_stride = 20;
    std::size_t top_n = 24;
    double theta_dom = 0.35;
    double theta_gap_A = 0.02;
    double min_isolated_fraction = 0.5;
    double min_driver_variance = 1.0e-12;
    double min_score = 1.0e-3;
};

struct CaseHunterStats {
    QMap<QString, std::size_t> candidate_counts;
    bool anti_circular_assertion = false;
};

CaseHunterStats RunCaseHunter(const Body& body,
                              const PerAtomSubstrateConfig& substrate_config,
                              const QString& out_dir,
                              const CaseHunterConfig& hunter_config = {});

}  // namespace h5reader::rediscover

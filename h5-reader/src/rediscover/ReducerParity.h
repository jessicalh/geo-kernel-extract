// ReducerParity -- gates between the two current reduction truth sources.
//
// This does not choose an authority. It runs the composed Relationship reducer
// and the existing per_atom_substrate pair reducer on the same (atom, frame)
// cases and reports drift for the migration gate.

#pragma once

#include "AnalysisBody.h"
#include "PerAtomSubstrate.h"
#include "Relationship.h"

#include <QString>
#include <QStringList>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace h5reader::rediscover {

struct ReducerParityOptions {
    double abs_tolerance = 1.0e-9;
    double rel_tolerance = 1.0e-9;
    std::size_t max_cases = 0;  // 0 means all requested cases.
};

struct ReducerAggregateSnapshot {
    bool present = false;
    double sum_all = 0.0;
    double sum_valid = 0.0;
    int n_valid = 0;
    std::vector<double> per_type;
};

struct ReducerParityStats {
    std::size_t cases_checked = 0;
    std::size_t ring_cases_checked = 0;
    std::size_t mc_cases_checked = 0;
    QStringList mismatches;

    bool ok() const { return mismatches.isEmpty(); }
};

inline bool ReducerNear(double a, double b, const ReducerParityOptions& options) {
    if (std::isnan(a) && std::isnan(b)) return true;
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    const double diff = std::abs(a - b);
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return diff <= options.abs_tolerance || diff <= options.rel_tolerance * scale;
}

inline void CompareReducerAggregateSnapshots(const ReducerAggregateSnapshot& lhs,
                                             const ReducerAggregateSnapshot& rhs,
                                             const QString& label,
                                             const ReducerParityOptions& options,
                                             QStringList* mismatches) {
    if (lhs.present != rhs.present) {
        *mismatches << QStringLiteral("%1 present mismatch: %2 vs %3")
                           .arg(label)
                           .arg(lhs.present ? 1 : 0)
                           .arg(rhs.present ? 1 : 0);
    }
    if (!ReducerNear(lhs.sum_all, rhs.sum_all, options)) {
        *mismatches << QStringLiteral("%1 sum_all mismatch: %2 vs %3")
                           .arg(label)
                           .arg(lhs.sum_all, 0, 'g', 17)
                           .arg(rhs.sum_all, 0, 'g', 17);
    }
    if (!ReducerNear(lhs.sum_valid, rhs.sum_valid, options)) {
        *mismatches << QStringLiteral("%1 sum_valid mismatch: %2 vs %3")
                           .arg(label)
                           .arg(lhs.sum_valid, 0, 'g', 17)
                           .arg(rhs.sum_valid, 0, 'g', 17);
    }
    if (lhs.n_valid != rhs.n_valid) {
        *mismatches << QStringLiteral("%1 n_valid mismatch: %2 vs %3")
                           .arg(label)
                           .arg(lhs.n_valid)
                           .arg(rhs.n_valid);
    }
    if (lhs.per_type.size() != rhs.per_type.size()) {
        *mismatches << QStringLiteral("%1 per_type size mismatch: %2 vs %3")
                           .arg(label)
                           .arg(static_cast<qulonglong>(lhs.per_type.size()))
                           .arg(static_cast<qulonglong>(rhs.per_type.size()));
        return;
    }
    for (std::size_t i = 0; i < lhs.per_type.size(); ++i) {
        if (!ReducerNear(lhs.per_type[i], rhs.per_type[i], options)) {
            *mismatches << QStringLiteral("%1 per_type[%2] mismatch: %3 vs %4")
                               .arg(label)
                               .arg(static_cast<qulonglong>(i))
                               .arg(lhs.per_type[i], 0, 'g', 17)
                               .arg(rhs.per_type[i], 0, 'g', 17);
        }
    }
}

ReducerParityStats AuditReducerParity(const Body& body,
                                      const std::vector<std::size_t>& atoms,
                                      const std::vector<std::size_t>& frames,
                                      const PerAtomSubstrateConfig& config = {},
                                      const ReducerParityOptions& options = {});

}  // namespace h5reader::rediscover

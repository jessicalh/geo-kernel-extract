// SubstrateParity -- old-output parity gates for the query migration.
//
// These helpers intentionally audit the consumed reader outputs as outputs,
// never as producer FieldKind entries. Producer truth is read through the
// Catalog FieldKind door, matching SpineReachabilityProbe's model.

#pragma once

#include "AnalysisBody.h"

#include "../model/Types.h"

#include <QString>
#include <QStringList>

#include <cstddef>
#include <functional>
#include <optional>

namespace h5reader::rediscover {

struct SubstrateParityOptions {
    double abs_tolerance = 1.0e-9;
    double rel_tolerance = 1.0e-9;
    std::size_t max_rows = 0;  // 0 means all rows.
};

struct SubstrateParityStats {
    std::size_t rows_seen = 0;
    std::size_t rows_checked = 0;
    std::size_t target_raw_components_checked = 0;
    std::size_t target_T0_checked = 0;
    std::size_t target_T1_checked = 0;
    std::size_t target_T2_checked = 0;
    QStringList errors;

    bool ok() const { return errors.isEmpty(); }
};

using ProducerTargetLookup =
    std::function<std::optional<model::Mat3>(std::size_t atom, std::size_t frame,
                                             QString* reason_out)>;

SubstrateParityStats AuditPerAtomSubstrateSidecars(
    const QString& outDir,
    std::size_t atomCount,
    const ProducerTargetLookup& producerTarget,
    const SubstrateParityOptions& options = {});

SubstrateParityStats AuditAllAtomEquivariantSidecars(
    const QString& outDir,
    const ProducerTargetLookup& producerTarget,
    const SubstrateParityOptions& options = {});

SubstrateParityStats AuditPerAtomSubstrateAgainstProducer(
    const Body& body,
    const QString& outDir,
    const SubstrateParityOptions& options = {});

SubstrateParityStats AuditAllAtomEquivariantAgainstProducer(
    const Body& body,
    const QString& outDir,
    const SubstrateParityOptions& options = {});

}  // namespace h5reader::rediscover

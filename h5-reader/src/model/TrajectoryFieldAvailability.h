#pragma once

#include "Conformation.h"
#include "DashboardSignal.h"
#include "QtConformationSnapshot.h"
#include "TrajectoryConformation.h"

#include "../io/QtFieldCatalog.gen.h"
#include "../io/QtTrajectoryH5.h"

#include <QHash>
#include <QString>
#include <QVector>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <type_traits>
#include <vector>

namespace h5reader::model {

enum class TrajectoryFieldAvailabilityState : std::uint8_t {
    Absent = 0,
    NoFramePayload,
    AllMissing,
    AllZeroStructural,
    AllZeroObserved,
    Available,
};

inline const char* ToString(TrajectoryFieldAvailabilityState state) {
    switch (state) {
        case TrajectoryFieldAvailabilityState::Absent:            return "Absent";
        case TrajectoryFieldAvailabilityState::NoFramePayload:    return "NoFramePayload";
        case TrajectoryFieldAvailabilityState::AllMissing:        return "AllMissing";
        case TrajectoryFieldAvailabilityState::AllZeroStructural: return "AllZeroStructural";
        case TrajectoryFieldAvailabilityState::AllZeroObserved:   return "AllZeroObserved";
        case TrajectoryFieldAvailabilityState::Available:         return "Available";
    }
    return "?";
}

struct TrajectoryFieldAvailabilityRecord {
    TrajectoryFieldAvailabilityState state = TrajectoryFieldAvailabilityState::Absent;
    qsizetype finiteSamples = 0;
    qsizetype nonZeroSamples = 0;
};

class TrajectoryFieldAvailability {
public:
    // Sizes of the startup-loaded topology spine's tables, snapshotted by the
    // caller from the loaded QtProtein. The gate uses these to classify
    // Topology-source descriptors honestly: the topology is loaded ONCE at
    // startup (not as a per-frame NPY column), so the per-frame probe can never
    // see it. A zero count is an honest empty table (e.g. rings on a ring-less
    // peptide), not a missing source.
    struct TopologyExtent {
        qsizetype atoms = 0;
        qsizetype bonds = 0;
        qsizetype residues = 0;
        qsizetype rings = 0;
        qsizetype ringMemberships = 0;
    };

    // conformation: dense-H5 + per-frame-NPY probing, as before.
    // topology:     table sizes of the startup-loaded topology spine.
    // dftJobCount:  number of DFT jobs registered for this run (== the .LGS
    //               dft.frames count == DftShieldingStore::jobCount()); 0 when
    //               there is no DFT campaign. The live ORCA shielding lives in
    //               .out files behind DftShieldingStore, NOT in the per-frame
    //               NPYs, so it must be classified from the job count — the old
    //               code routed it through the NPY probe, which can only ever
    //               report it Absent.
    static TrajectoryFieldAvailability Build(Conformation* conformation,
                                             const TopologyExtent& topology,
                                             std::size_t dftJobCount,
                                             bool experimentalShieldingMlAvailable,
                                             const QVector<SignalDescriptor>& descriptors) {
        TrajectoryFieldAvailability out;
        const auto* traj = conformation ? conformation->asTrajectory() : nullptr;
        const auto* h5 = traj ? traj->h5() : nullptr;

        QHash<QString, TrajectoryFieldAvailabilityRecord> denseByStorage;
        QHash<QString, TrajectoryFieldAvailabilityRecord> frameByStorage;
        QHash<QString, io::FieldKind> frameFields;

        for (const SignalDescriptor& descriptor : descriptors) {
            if (descriptor.sourceKind == SignalSourceKind::DenseH5Trajectory) {
                const QString key = descriptor.storagePath;
                if (!denseByStorage.contains(key))
                    denseByStorage.insert(key, classifyDenseH5(h5, descriptor));
            } else if (descriptor.sourceKind == SignalSourceKind::FrameNpySnapshot) {
                const QString key = descriptor.storagePath;
                if (key.isEmpty()) {
                    frameByStorage.insert(key, {TrajectoryFieldAvailabilityState::Absent, 0, 0});
                    continue;
                }
                const std::optional<io::FieldKind> kind =
                    io::FindFieldByStem(key.toStdString());
                if (!kind) {
                    frameByStorage.insert(key, {TrajectoryFieldAvailabilityState::Absent, 0, 0});
                    continue;
                }
                frameFields.insert(key, *kind);
            }
        }

        if (!frameFields.isEmpty()) {
            const auto classified = classifyFramePayloads(conformation, traj, frameFields);
            for (auto it = classified.constBegin(); it != classified.constEnd(); ++it)
                frameByStorage.insert(it.key(), it.value());
        }

        // The live ORCA shielding is classified once from the registered job
        // count (it feeds the CSA glyph via DftShieldingStore; it is never a
        // per-frame NPY column). DerivedGeometry/SelectionEvents are offerable
        // affordances computed/created on demand — always offerable once a run
        // is loaded. Marking these EXPLICITLY (rather than via the old silent
        // no-record default) is what lets the missing-record default be Absent.
        const TrajectoryFieldAvailabilityRecord dftRecord = classifyDft(dftJobCount);
        const TrajectoryFieldAvailabilityRecord experimentalMlRecord =
            classifyExperimentalShieldingMl(experimentalShieldingMlAvailable);
        const TrajectoryFieldAvailabilityRecord affordanceRecord{
            TrajectoryFieldAvailabilityState::Available, 1, 1};

        for (const SignalDescriptor& descriptor : descriptors) {
            TrajectoryFieldAvailabilityRecord record;
            bool hasRecord = false;
            switch (descriptor.sourceKind) {
                case SignalSourceKind::DenseH5Trajectory: {
                    const auto it = denseByStorage.constFind(descriptor.storagePath);
                    if (it != denseByStorage.constEnd()) { record = it.value(); hasRecord = true; }
                    break;
                }
                case SignalSourceKind::FrameNpySnapshot: {
                    const auto it = frameByStorage.constFind(descriptor.storagePath);
                    if (it != frameByStorage.constEnd()) { record = it.value(); hasRecord = true; }
                    break;
                }
                case SignalSourceKind::OrcaDftFrame:
                    record = dftRecord;
                    hasRecord = true;
                    break;
                case SignalSourceKind::Topology:
                    record = classifyTopology(descriptor.storagePath, topology);
                    hasRecord = true;
                    break;
                case SignalSourceKind::DerivedGeometry:
                case SignalSourceKind::SelectionEvents:
                    record = affordanceRecord;
                    hasRecord = true;
                    break;
                case SignalSourceKind::ExperimentalShieldingMl:
                    record = experimentalMlRecord;
                    hasRecord = true;
                    break;
            }
            if (hasRecord) {
                out.byDescriptor_.insert(descriptor.id, record);
                upgradeStorage(out.byStorage_, descriptor.storagePath, record);
            }
        }

        return out;
    }

    // Classify the ORCA DFT source from the registered job count: a campaign of
    // N>0 jobs is live (Available); no campaign (N==0) is honestly Absent. Pure;
    // unit-tested directly (Build itself is integration-tested via REST because
    // it is coupled to the H5/NPY I/O it probes).
    static TrajectoryFieldAvailabilityRecord classifyDft(std::size_t jobCount) {
        if (jobCount == 0)
            return {TrajectoryFieldAvailabilityState::Absent, 0, 0};
        const auto n = static_cast<qsizetype>(jobCount);
        return {TrajectoryFieldAvailabilityState::Available, n, n};
    }

    static TrajectoryFieldAvailabilityRecord classifyExperimentalShieldingMl(bool available) {
        if (!available)
            return {TrajectoryFieldAvailabilityState::Absent, 0, 0};
        // Runtime files alone do not give the reader a per-frame payload. Keep
        // the descriptors quarantined until a graph builder + sample cache exist.
        return {TrajectoryFieldAvailabilityState::NoFramePayload, 0, 0};
    }

    // Classify a Topology-source descriptor from the loaded spine's table size.
    // A non-empty table is live (Available); an empty one is an honest
    // AllMissing (the spine loaded, but this protein has no rows for it). Pure.
    static TrajectoryFieldAvailabilityRecord classifyTopology(const QString& storagePath,
                                                              const TopologyExtent& t) {
        qsizetype n = -1;
        if (storagePath == QLatin1String("atoms")) n = t.atoms;
        else if (storagePath == QLatin1String("residues")) n = t.residues;
        else if (storagePath == QLatin1String("bonds")
                 || storagePath == QLatin1String("bond_length")) n = t.bonds;
        else if (storagePath == QLatin1String("rings")) n = t.rings;
        else if (storagePath == QLatin1String("ring_membership")) n = t.ringMemberships;
        if (n < 0)  // unknown topology stem: live iff the spine itself loaded.
            n = t.atoms;
        if (n <= 0)
            return {TrajectoryFieldAvailabilityState::AllMissing, 0, 0};
        return {TrajectoryFieldAvailabilityState::Available, n, n};
    }

    const TrajectoryFieldAvailabilityRecord* recordForDescriptor(const QString& descriptorId) const {
        const auto it = byDescriptor_.constFind(descriptorId);
        return it == byDescriptor_.constEnd() ? nullptr : &it.value();
    }

    const TrajectoryFieldAvailabilityRecord* recordForStoragePath(const QString& path) const {
        const auto it = byStorage_.constFind(path);
        return it == byStorage_.constEnd() ? nullptr : &it.value();
    }

    TrajectoryFieldAvailabilityState stateForDescriptor(const QString& descriptorId) const {
        const auto* record = recordForDescriptor(descriptorId);
        // Build() now records every source kind, so a missing record means an id
        // not in the catalog (or a future, unhandled kind) -- report it
        // conservatively as Absent. The old default of Available was the root of
        // the topology/ORCA "lies-Available" leak: an unhandled source kind got
        // no record and was silently reported live.
        return record ? record->state : TrajectoryFieldAvailabilityState::Absent;
    }

    bool allowsDescriptor(const SignalDescriptor& descriptor) const {
        const auto* record = recordForDescriptor(descriptor.id);
        if (!record)
            return true;
        return isVisibleState(record->state);
    }

    bool canSampleDescriptor(const SignalDescriptor& descriptor) const {
        const auto* record = recordForDescriptor(descriptor.id);
        if (!record)
            return true;
        return isSampleableState(record->state);
    }

    static bool isVisibleState(TrajectoryFieldAvailabilityState state) {
        switch (state) {
            case TrajectoryFieldAvailabilityState::Available:
            case TrajectoryFieldAvailabilityState::AllZeroObserved:
                return true;
            case TrajectoryFieldAvailabilityState::Absent:
            case TrajectoryFieldAvailabilityState::NoFramePayload:
            case TrajectoryFieldAvailabilityState::AllMissing:
            case TrajectoryFieldAvailabilityState::AllZeroStructural:
                return false;
        }
        return false;
    }

    static bool isSampleableState(TrajectoryFieldAvailabilityState state) {
        return isVisibleState(state);
    }

private:
    struct ScanAccum {
        bool sourcePresent = false;
        bool payloadPresent = false;
        qsizetype finiteSamples = 0;
        qsizetype nonZeroSamples = 0;
    };

    static QString trajectoryGroupName(QString storagePath) {
        const QString prefix = QStringLiteral("/trajectory/");
        if (storagePath.startsWith(prefix))
            storagePath.remove(0, prefix.size());
        return storagePath;
    }

    static bool structuralZeroDescriptor(const SignalDescriptor& /*descriptor*/) {
        // Placeholder for fields whose producer contract states that an
        // all-zero payload is structural rather than observed. Current loaded
        // buffers either have masks/absence, or all-zero is still a value the
        // user may need to see, so the conservative default is observed zero.
        return false;
    }

    static TrajectoryFieldAvailabilityRecord finish(const ScanAccum& acc,
                                                    bool structuralZero) {
        if (!acc.sourcePresent)
            return {TrajectoryFieldAvailabilityState::Absent, 0, 0};
        if (!acc.payloadPresent || acc.finiteSamples == 0)
            return {TrajectoryFieldAvailabilityState::AllMissing, 0, 0};
        if (acc.nonZeroSamples == 0) {
            return {structuralZero
                        ? TrajectoryFieldAvailabilityState::AllZeroStructural
                        : TrajectoryFieldAvailabilityState::AllZeroObserved,
                    acc.finiteSamples,
                    0};
        }
        return {TrajectoryFieldAvailabilityState::Available,
                acc.finiteSamples,
                acc.nonZeroSamples};
    }

    template <typename T>
    static void scanValues(ScanAccum& acc, const std::vector<T>& values) {
        if (values.empty())
            return;
        acc.payloadPresent = true;
        for (const T& raw : values) {
            const double v = static_cast<double>(raw);
            if constexpr (std::is_floating_point_v<T>) {
                if (!std::isfinite(v))
                    continue;
            }
            ++acc.finiteSamples;
            if (std::abs(v) > 1e-12)
                ++acc.nonZeroSamples;
        }
    }

    static bool anyAttached(const QtTimeSeriesFrameMeta& meta) {
        return meta.source_attached.empty()
               || std::any_of(meta.source_attached.begin(), meta.source_attached.end(),
                              [](std::uint8_t v) { return v != 0; });
    }

    static bool anyAttached(const QtPerResidueFrameMeta& meta) {
        return meta.source_attached.empty()
               || std::any_of(meta.source_attached.begin(), meta.source_attached.end(),
                              [](std::uint8_t v) { return v != 0; });
    }

    static bool anyAttached(const std::vector<std::uint8_t>& mask) {
        return mask.empty()
               || std::any_of(mask.begin(), mask.end(),
                              [](std::uint8_t v) { return v != 0; });
    }

    template <typename Buffer, typename Values>
    static void scanTimeSeries(ScanAccum& acc, const Buffer* buffer,
                               const Values& values,
                               const QtTimeSeriesFrameMeta& meta) {
        if (!buffer)
            return;
        acc.sourcePresent = true;
        if (!anyAttached(meta))
            return;
        scanValues(acc, values);
    }

    static TrajectoryFieldAvailabilityRecord classifyDenseH5(
        const io::QtTrajectoryH5* h5,
        const SignalDescriptor& descriptor) {
        if (!h5)
            return {TrajectoryFieldAvailabilityState::Absent, 0, 0};

        const QString path = descriptor.storagePath;
        const QString group = trajectoryGroupName(path);
        ScanAccum acc;
        acc.sourcePresent = (group == QStringLiteral("positions")) || h5->hasGroup(group);
        if (!acc.sourcePresent)
            return finish(acc, structuralZeroDescriptor(descriptor));

        if (path == QStringLiteral("/trajectory/positions")) {
            if (const auto* b = h5->positions()) {
                acc.sourcePresent = true;
                scanValues(acc, b->xyz);
            }
        } else if (path == QStringLiteral("/trajectory/bs_shielding_time_series")) {
            scanTimeSeries(acc, h5->bsShielding(), h5->bsShielding() ? h5->bsShielding()->xyz : emptyDouble(), h5->bsShielding() ? h5->bsShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/hm_shielding_time_series")) {
            scanTimeSeries(acc, h5->hmShielding(), h5->hmShielding() ? h5->hmShielding()->xyz : emptyDouble(), h5->hmShielding() ? h5->hmShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/mc_shielding_time_series")) {
            scanTimeSeries(acc, h5->mcShielding(), h5->mcShielding() ? h5->mcShielding()->xyz : emptyDouble(), h5->mcShielding() ? h5->mcShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/mopac_coulomb_shielding_time_series")) {
            scanTimeSeries(acc, h5->mopacCoulombShielding(), h5->mopacCoulombShielding() ? h5->mopacCoulombShielding()->t2 : emptyDouble(), h5->mopacCoulombShielding() ? h5->mopacCoulombShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/mopac_mc_shielding_time_series")) {
            scanTimeSeries(acc, h5->mopacMcShielding(), h5->mopacMcShielding() ? h5->mopacMcShielding()->xyz : emptyDouble(), h5->mopacMcShielding() ? h5->mopacMcShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/mopac_vs_ff14sb_reconciliation")) {
            scanTimeSeries(acc, h5->mopacVsFf14sbReconciliation(), h5->mopacVsFf14sbReconciliation() ? h5->mopacVsFf14sbReconciliation()->data : emptyDouble(), h5->mopacVsFf14sbReconciliation() ? h5->mopacVsFf14sbReconciliation()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/tripeptide_bb_shielding_time_series")) {
            scanTimeSeries(acc, h5->tripeptideBbShielding(), h5->tripeptideBbShielding() ? h5->tripeptideBbShielding()->xyz : emptyDouble(), h5->tripeptideBbShielding() ? h5->tripeptideBbShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/tripeptide_neighbor_shielding_time_series")) {
            scanTimeSeries(acc, h5->tripeptideNeighborShielding(), h5->tripeptideNeighborShielding() ? h5->tripeptideNeighborShielding()->xyz : emptyDouble(), h5->tripeptideNeighborShielding() ? h5->tripeptideNeighborShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/larsen_hbond_1pHB_shielding_time_series")) {
            scanTimeSeries(acc, h5->larsenHBond1pHBShielding(), h5->larsenHBond1pHBShielding() ? h5->larsenHBond1pHBShielding()->xyz : emptyDouble(), h5->larsenHBond1pHBShielding() ? h5->larsenHBond1pHBShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/larsen_hbond_1pHaB_shielding_time_series")) {
            scanTimeSeries(acc, h5->larsenHBond1pHaBShielding(), h5->larsenHBond1pHaBShielding() ? h5->larsenHBond1pHaBShielding()->xyz : emptyDouble(), h5->larsenHBond1pHaBShielding() ? h5->larsenHBond1pHaBShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/larsen_hbond_2pHB_shielding_time_series")) {
            scanTimeSeries(acc, h5->larsenHBond2pHBShielding(), h5->larsenHBond2pHBShielding() ? h5->larsenHBond2pHBShielding()->xyz : emptyDouble(), h5->larsenHBond2pHBShielding() ? h5->larsenHBond2pHBShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/larsen_hbond_2pHaB_shielding_time_series")) {
            scanTimeSeries(acc, h5->larsenHBond2pHaBShielding(), h5->larsenHBond2pHaBShielding() ? h5->larsenHBond2pHaBShielding()->xyz : emptyDouble(), h5->larsenHBond2pHaBShielding() ? h5->larsenHBond2pHaBShielding()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/sasa_time_series")) {
            scanTimeSeries(acc, h5->sasa(), h5->sasa() ? h5->sasa()->data : emptyDouble(), h5->sasa() ? h5->sasa()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/aimnet2_charge_time_series")) {
            scanTimeSeries(acc, h5->aimnet2Charge(), h5->aimnet2Charge() ? h5->aimnet2Charge()->data : emptyDouble(), h5->aimnet2Charge() ? h5->aimnet2Charge()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/larsen_hbond_count_time_series")) {
            scanTimeSeries(acc, h5->larsenHBondCount(), h5->larsenHBondCount() ? h5->larsenHBondCount()->data : emptyDouble(), h5->larsenHBondCount() ? h5->larsenHBondCount()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/larsen_hbond_water_term_time_series")) {
            scanTimeSeries(acc, h5->larsenHBondWaterTerm(), h5->larsenHBondWaterTerm() ? h5->larsenHBondWaterTerm()->data : emptyDouble(), h5->larsenHBondWaterTerm() ? h5->larsenHBondWaterTerm()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/bonded_energy_time_series")) {
            scanTimeSeries(acc, h5->bondedEnergyTotal(), h5->bondedEnergyTotal() ? h5->bondedEnergyTotal()->data : emptyDouble(), h5->bondedEnergyTotal() ? h5->bondedEnergyTotal()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/apbs_efield_time_series")) {
            scanTimeSeries(acc, h5->apbsEfield(), h5->apbsEfield() ? h5->apbsEfield()->xyz : emptyDouble(), h5->apbsEfield() ? h5->apbsEfield()->meta : emptyMeta());
        } else if (path == QStringLiteral("/trajectory/apbs_efg_time_series")) {
            scanTimeSeries(acc, h5->apbsEfg(), h5->apbsEfg() ? h5->apbsEfg()->t2 : emptyDouble(), h5->apbsEfg() ? h5->apbsEfg()->meta : emptyMeta());
        } else {
            // For composite/special groups the typed object exists only when
            // the group parsed successfully. Group presence is already enough
            // to remove the old "implemented in code but absent in run" leak;
            // detailed per-channel scans can be added here as a narrow follow-up.
            acc.payloadPresent = true;
            acc.finiteSamples = 1;
            acc.nonZeroSamples = 1;
        }
        return finish(acc, structuralZeroDescriptor(descriptor));
    }

    static QHash<QString, TrajectoryFieldAvailabilityRecord> classifyFramePayloads(
        Conformation* conformation,
        const TrajectoryConformation* traj,
        const QHash<QString, io::FieldKind>& fields) {
        QHash<QString, TrajectoryFieldAvailabilityRecord> out;
        if (!conformation) {
            for (auto it = fields.constBegin(); it != fields.constEnd(); ++it)
                out.insert(it.key(), {TrajectoryFieldAvailabilityState::NoFramePayload, 0, 0});
            return out;
        }

        // Availability only needs to know whether a field EVER carries data;
        // per-frame NPY presence is uniform across frames (the producer emits
        // the same field set every frame), so a few representative sampled
        // snapshots answer it. Scanning every sampled row would do
        // K x (~50-100 NPY parses) synchronously on the GUI thread at open —
        // a trajectory-size-scaling open stall (adversarial review H1). Cap at
        // a small probe; the dense-H5 half stays exhaustive (it is cheap).
        constexpr qsizetype kMaxAvailabilityProbeRows = 4;
        QVector<std::size_t> rows;
        if (traj) {
            for (std::size_t row : traj->sampledFrameRows()) {
                rows.push_back(row);
                if (rows.size() >= kMaxAvailabilityProbeRows)
                    break;
            }
        } else if (conformation->frameCount() > 0) {
            rows.push_back(0);
        }

        if (rows.isEmpty()) {
            for (auto it = fields.constBegin(); it != fields.constEnd(); ++it)
                out.insert(it.key(), {TrajectoryFieldAvailabilityState::NoFramePayload, 0, 0});
            return out;
        }

        QHash<QString, ScanAccum> accByStorage;
        for (auto it = fields.constBegin(); it != fields.constEnd(); ++it)
            accByStorage.insert(it.key(), {});

        bool anySnapshot = false;
        for (std::size_t row : rows) {
            conformation->requestSnapshot(row);
            const auto snap = conformation->snapshot(row);
            if (!snap)
                continue;
            anySnapshot = true;
            for (auto it = fields.constBegin(); it != fields.constEnd(); ++it) {
                ScanAccum& acc = accByStorage[it.key()];
                if (!snap->has(it.value()))
                    continue;
                const NpyColumn& col = snap->column(it.value());
                acc.sourcePresent = true;
                if (col.present && col.rows > 0 && col.cols > 0)
                    scanValues(acc, col.data);
            }
        }

        for (auto it = fields.constBegin(); it != fields.constEnd(); ++it) {
            if (!anySnapshot) {
                out.insert(it.key(), {TrajectoryFieldAvailabilityState::NoFramePayload, 0, 0});
                continue;
            }
            out.insert(it.key(), finish(accByStorage.value(it.key()), false));
        }
        return out;
    }

    static const std::vector<double>& emptyDouble() {
        static const std::vector<double> empty;
        return empty;
    }

    static const QtTimeSeriesFrameMeta& emptyMeta() {
        static const QtTimeSeriesFrameMeta empty;
        return empty;
    }

    // Merge a record into the by-storage map, preferring the MORE-available one
    // when two descriptors alias the same storage stem (the orca_total /
    // residues duplication: a live `orca_dft:*` and a dead `npy:orca_*` share
    // the stem). recordForStoragePath is the secondary gate in
    // VisualizationDescriptorDataAvailable, so a dead duplicate's Absent must
    // not mask the live aliased descriptor here.
    static void upgradeStorage(QHash<QString, TrajectoryFieldAvailabilityRecord>& byStorage,
                               const QString& path,
                               const TrajectoryFieldAvailabilityRecord& record) {
        if (path.isEmpty())
            return;
        const auto it = byStorage.constFind(path);
        if (it == byStorage.constEnd()
            || (isVisibleState(record.state) && !isVisibleState(it.value().state)))
            byStorage.insert(path, record);
    }

    QHash<QString, TrajectoryFieldAvailabilityRecord> byDescriptor_;
    QHash<QString, TrajectoryFieldAvailabilityRecord> byStorage_;
};

}  // namespace h5reader::model

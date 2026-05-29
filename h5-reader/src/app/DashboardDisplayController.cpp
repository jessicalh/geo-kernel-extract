#include "DashboardDisplayController.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/QtFieldCatalog.gen.h"
#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/ConformationGeometry.h"
#include "../model/QtConformationSnapshot.h"
#include "../model/DashboardPanelModel.h"
#include "../model/DashboardSignalModel.h"
#include "../model/DftShieldingStore.h"
#include "../model/QtProtein.h"
#include "../model/TrajectoryConformation.h"
#include "../model/TrajectorySignalCatalog.h"

#include <QAbstractItemModel>
#include <QSet>
#include <QStringList>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>

namespace h5reader::app {

namespace {
bool isStripMode(const QString& mode) {
    return mode.startsWith(QStringLiteral("strip."));
}

bool hasStripMode(const QStringList& modes) {
    return std::any_of(modes.begin(), modes.end(), isStripMode);
}

QString canonicalModeChannel(const QString& mode) {
    if (mode.startsWith(QStringLiteral("strip.tensor.")))
        return mode.mid(QStringLiteral("strip.tensor.").size());
    if (mode == QStringLiteral("strip.vector.magnitude"))
        return QStringLiteral("magnitude");
    return {};
}

bool modeWantsChannel(const QString& mode, const model::ChannelDescriptor& channel) {
    const QString channelId = canonicalModeChannel(mode);
    if (!channelId.isEmpty())
        return channel.id.compare(channelId, Qt::CaseInsensitive) == 0;
    if (mode == QStringLiteral("strip.vector.component"))
        return channel.id != QStringLiteral("magnitude");
    return true;
}

model::SignalBinding bindingFromAnchor(const model::SignalDescriptor& descriptor,
                                       const model::SignalAnchor& anchor,
                                       bool followsFocus) {
    model::SignalBinding binding;
    binding.descriptorId = descriptor.id;
    binding.conceptKey = descriptor.conceptKey;
    binding.anchor = anchor;
    binding.followsFocus = followsFocus;
    return binding;
}

bool bindingHasRevealTarget(const model::SignalBinding& binding) {
    switch (model::AxisForAnchor(binding.anchor)) {
    case model::SignalAxis::Atom:
    case model::SignalAxis::Residue:
    case model::SignalAxis::Bond:
    case model::SignalAxis::Ring:
    case model::SignalAxis::AromaticRing:
    case model::SignalAxis::SaturatedRing:
    case model::SignalAxis::RingMembership:
        return true;
    case model::SignalAxis::AtomTuple:
        if (const auto* tuple = std::get_if<model::AtomTupleAnchor>(&binding.anchor))
            return !tuple->atoms.empty();
        return false;
    case model::SignalAxis::None:
    case model::SignalAxis::RingContributionPair:
    case model::SignalAxis::MutationMatchPair:
    case model::SignalAxis::Protein:
    case model::SignalAxis::System:
    case model::SignalAxis::Event:
        return false;
    }
    return false;
}

bool unitSpecPresent(const model::UnitSpec& units) {
    return !units.sourceSymbol.isEmpty() || !units.displaySymbol.isEmpty()
           || units.scaleToDisplay != 1.0 || units.offsetToDisplay != 0.0;
}

double applyDisplayUnits(double value,
                         const model::SignalDescriptor& descriptor,
                         const model::ChannelDescriptor& channel) {
    const model::UnitSpec units = unitSpecPresent(channel.defaultDisplayUnits)
                                      ? channel.defaultDisplayUnits
                                      : descriptor.defaultDisplayUnits;
    return value * units.scaleToDisplay + units.offsetToDisplay;
}

std::optional<std::size_t> firstAtomInResidue(const model::QtProtein* protein, std::size_t residue) {
    if (!protein)
        return std::nullopt;
    for (std::size_t atom = 0; atom < protein->atomCount(); ++atom) {
        if (protein->atom(atom).residueIndex >= 0
            && static_cast<std::size_t>(protein->atom(atom).residueIndex) == residue) {
            return atom;
        }
    }
    return std::nullopt;
}

std::optional<std::size_t> anchorRow(const model::SignalAnchor& anchor,
                                     model::SignalAxis axis,
                                     const model::QtProtein* protein,
                                     int rows) {
    if (rows <= 0)
        return std::nullopt;
    auto inRange = [rows](std::size_t row) -> std::optional<std::size_t> {
        return row < static_cast<std::size_t>(rows) ? std::optional<std::size_t>(row) : std::nullopt;
    };

    if (const auto* atom = std::get_if<model::AtomAnchor>(&anchor))
        return inRange(atom->atom);
    if (const auto* residue = std::get_if<model::ResidueAnchor>(&anchor)) {
        if (protein && rows == static_cast<int>(protein->atomCount())) {
            if (const auto atom = firstAtomInResidue(protein, residue->residue))
                return atom;
        }
        return inRange(residue->residue);
    }
    if (const auto* tuple = std::get_if<model::AtomTupleAnchor>(&anchor)) {
        if (!tuple->atoms.empty())
            return inRange(tuple->atoms.front());
        return std::nullopt;
    }
    if (const auto* bond = std::get_if<model::BondAnchor>(&anchor))
        return inRange(bond->bond);
    if (const auto* ring = std::get_if<model::RingAnchor>(&anchor))
        return inRange(ring->ring);
    if (const auto* ring = std::get_if<model::AromaticRingAnchor>(&anchor)) {
        if (axis == model::SignalAxis::Ring && protein) {
            const auto absolute = protein->topology().absoluteRingIndex(model::QtRingAxis::AromaticRing, ring->ring);
            return absolute ? inRange(*absolute) : std::nullopt;
        }
        return inRange(ring->ring);
    }
    if (const auto* ring = std::get_if<model::SaturatedRingAnchor>(&anchor)) {
        if (axis == model::SignalAxis::Ring && protein) {
            const auto absolute = protein->topology().absoluteRingIndex(model::QtRingAxis::SaturatedRing, ring->ring);
            return absolute ? inRange(*absolute) : std::nullopt;
        }
        return inRange(ring->ring);
    }
    if (const auto* pair = std::get_if<model::RingContributionPairAnchor>(&anchor))
        return inRange(pair->pair);
    if (const auto* membership = std::get_if<model::RingMembershipAnchor>(&anchor))
        return inRange(membership->membership);
    if (const auto* pair = std::get_if<model::MutationMatchPairAnchor>(&anchor))
        return inRange(pair->pair);

    switch (axis) {
    case model::SignalAxis::None:
    case model::SignalAxis::Protein:
    case model::SignalAxis::System:
    case model::SignalAxis::Event:
        return 0;
    case model::SignalAxis::Atom:
    case model::SignalAxis::Residue:
    case model::SignalAxis::AtomTuple:
    case model::SignalAxis::Bond:
    case model::SignalAxis::Ring:
    case model::SignalAxis::AromaticRing:
    case model::SignalAxis::SaturatedRing:
    case model::SignalAxis::RingContributionPair:
    case model::SignalAxis::RingMembership:
    case model::SignalAxis::MutationMatchPair:
        break;
    }
    return rows == 1 ? std::optional<std::size_t>(0) : std::nullopt;
}

double magnitude(const double* row, int first, int count, int cols) {
    double sum = 0.0;
    int used = 0;
    for (int i = 0; i < count; ++i) {
        const int col = first + i;
        if (col < 0 || col >= cols)
            continue;
        const double value = row[col];
        if (!std::isfinite(value))
            return std::numeric_limits<double>::quiet_NaN();
        sum += value * value;
        ++used;
    }
    return used > 0 ? std::sqrt(sum) : std::numeric_limits<double>::quiet_NaN();
}

std::optional<double> sampleNpyValue(const model::NpyColumn& column,
                                     std::size_t rowIndex,
                                     const model::ChannelDescriptor& channel,
                                     const QString& displayModeId) {
    if (rowIndex >= static_cast<std::size_t>(column.rows) || column.cols <= 0)
        return std::nullopt;

    const double* row = column.row(rowIndex);
    const QString id = channel.id.toLower();

    if (id == QStringLiteral("x"))
        return column.cols > 0 ? std::optional<double>(row[0]) : std::nullopt;
    if (id == QStringLiteral("y"))
        return column.cols > 1 ? std::optional<double>(row[1]) : std::nullopt;
    if (id == QStringLiteral("z"))
        return column.cols > 2 ? std::optional<double>(row[2]) : std::nullopt;
    if (id == QStringLiteral("magnitude") || displayModeId == QStringLiteral("strip.vector.magnitude"))
        return magnitude(row, 0, std::min(3, column.cols), column.cols);
    if (id == QStringLiteral("t0") || displayModeId == QStringLiteral("strip.tensor.T0"))
        return column.cols > 0 ? std::optional<double>(row[0]) : std::nullopt;
    if (id == QStringLiteral("t1") || displayModeId == QStringLiteral("strip.tensor.T1"))
        return column.cols >= 4 ? std::optional<double>(magnitude(row, 1, 3, column.cols)) : std::nullopt;
    if (id == QStringLiteral("t2") || displayModeId == QStringLiteral("strip.tensor.T2")) {
        if (column.cols >= 9)
            return magnitude(row, 4, 5, column.cols);
        if (column.cols >= 5)
            return magnitude(row, 0, 5, column.cols);
        return std::nullopt;
    }

    return std::optional<double>(row[0]);
}

model::FrameSignalSample finiteSample(double value) {
    if (!std::isfinite(value))
        return model::FrameSignalSample::Gap(model::GapReason::NaNSentinel);
    return model::FrameSignalSample::Valid(value);
}

template <std::size_t N>
double arrayMagnitude(const std::array<double, N>& values) {
    double sum = 0.0;
    for (double value : values) {
        if (!std::isfinite(value))
            return std::numeric_limits<double>::quiet_NaN();
        sum += value * value;
    }
    return std::sqrt(sum);
}

std::optional<double> sampleTensorValue(const model::SphericalTensor& tensor,
                                        const model::ChannelDescriptor& channel,
                                        const QString& displayModeId) {
    const QString id = channel.id.toLower();
    if (id == QStringLiteral("t0") || displayModeId == QStringLiteral("strip.tensor.T0"))
        return tensor.T0;
    if (id == QStringLiteral("t1") || displayModeId == QStringLiteral("strip.tensor.T1"))
        return arrayMagnitude(tensor.T1);
    if (id == QStringLiteral("t2")
        || id == QStringLiteral("magnitude")
        || displayModeId == QStringLiteral("strip.tensor.T2")) {
        return tensor.T2Magnitude();
    }
    if (id == QStringLiteral("component") || displayModeId == QStringLiteral("strip.tensor.component"))
        return tensor.T2[0];
    return tensor.T0;
}

std::optional<double> sampleT2Value(const std::array<double, 5>& t2,
                                    const model::ChannelDescriptor& channel,
                                    const QString& displayModeId) {
    const QString id = channel.id.toLower();
    if (id == QStringLiteral("component") || displayModeId == QStringLiteral("strip.tensor.component"))
        return t2[0];
    return arrayMagnitude(t2);
}

std::optional<double> sampleVecValue(const model::Vec3& value,
                                     const model::ChannelDescriptor& channel,
                                     const QString& displayModeId) {
    const QString id = channel.id.toLower();
    if (id == QStringLiteral("x"))
        return value.x();
    if (id == QStringLiteral("y"))
        return value.y();
    if (id == QStringLiteral("z"))
        return value.z();
    if (id == QStringLiteral("magnitude") || displayModeId == QStringLiteral("strip.vector.magnitude"))
        return value.norm();
    return value.x();
}

std::optional<std::size_t> atomFromAnchor(const model::SignalAnchor& anchor,
                                          const model::QtProtein* protein) {
    if (const auto* atom = std::get_if<model::AtomAnchor>(&anchor))
        return (!protein || atom->atom < protein->atomCount()) ? std::optional<std::size_t>(atom->atom)
                                                               : std::nullopt;
    if (const auto* residue = std::get_if<model::ResidueAnchor>(&anchor))
        return firstAtomInResidue(protein, residue->residue);
    if (const auto* tuple = std::get_if<model::AtomTupleAnchor>(&anchor)) {
        if (!tuple->atoms.empty() && (!protein || tuple->atoms.front() < protein->atomCount()))
            return tuple->atoms.front();
    }
    return std::nullopt;
}

std::optional<std::vector<std::size_t>> atomTupleFromAnchor(const model::SignalAnchor& anchor,
                                                            const model::QtProtein* protein) {
    if (const auto* tuple = std::get_if<model::AtomTupleAnchor>(&anchor)) {
        if (tuple->atoms.size() < 2)
            return std::nullopt;
        if (protein) {
            for (std::size_t atom : tuple->atoms) {
                if (atom >= protein->atomCount())
                    return std::nullopt;
            }
        }
        return tuple->atoms;
    }
    return std::nullopt;
}

struct SamplePlan {
    std::function<model::FrameSignalSample(std::size_t frame)> sample;
    bool needsFrameSnapshot = false;
    bool needsDftFrame = false;
};

SamplePlan pendingPlan() {
    SamplePlan plan;
    plan.sample = [](std::size_t) {
        return model::FrameSignalSample::Gap(model::GapReason::Pending);
    };
    return plan;
}

SamplePlan frameNpyPlan(const model::SignalDescriptor& descriptor,
                        const model::ChannelDescriptor& channel,
                        const QString& displayModeId,
                        const model::SignalAnchor& anchor,
                        const model::QtProtein* protein,
                        const QPointer<model::Conformation>& conformation) {
    SamplePlan plan;
    plan.needsFrameSnapshot = true;

    const std::string stem = descriptor.storagePath.toStdString();
    const std::optional<h5reader::io::FieldKind> fieldKind = h5reader::io::FindFieldByStem(stem);
    if (!fieldKind) {
        plan.sample = [](std::size_t) {
            return model::FrameSignalSample::Gap(model::GapReason::SourceAbsent);
        };
        return plan;
    }

    plan.sample = [descriptor, channel, displayModeId, anchor, protein, conformation, fieldKind](std::size_t frame) {
        if (!conformation)
            return model::FrameSignalSample::Gap(model::GapReason::SourceAbsent);
        const std::shared_ptr<const model::QtConformationSnapshot> snapshot = conformation->snapshot(frame);
        if (!snapshot)
            return model::FrameSignalSample::Gap(model::GapReason::FrameSourceAbsent);
        const model::NpyColumn& column = snapshot->column(*fieldKind);
        if (!column.present)
            return model::FrameSignalSample::Gap(model::GapReason::SourceAbsent);
        const std::optional<std::size_t> row = anchorRow(anchor, descriptor.requiredAnchor, protein, column.rows);
        if (!row)
            return model::FrameSignalSample::Gap(model::GapReason::AnchorUnavailable);
        const std::optional<double> raw = sampleNpyValue(column, *row, channel, displayModeId);
        if (!raw)
            return model::FrameSignalSample::Gap(model::GapReason::MalformedSource);
        const double value = applyDisplayUnits(*raw, descriptor, channel);
        if (!std::isfinite(value))
            return model::FrameSignalSample::Gap(model::GapReason::NaNSentinel);
        return model::FrameSignalSample::Valid(value);
    };
    return plan;
}

// The dense H5 router is deliberately explicit: it is the visible handoff from
// catalog descriptors to concrete typed H5 buffers. Keep these storagePath cases
// in step with TrajectorySignalCatalog and the dashboard smoke coverage; if this
// grows further, the next move is a table of descriptor-path samplers, not hidden
// reflection or ad hoc dynamic lookup.
SamplePlan denseH5Plan(const model::SignalDescriptor& descriptor,
                       const model::ChannelDescriptor& channel,
                       const QString& displayModeId,
                       const model::SignalAnchor& anchor,
                       const model::QtProtein* protein,
                       const QPointer<model::Conformation>& conformation) {
    SamplePlan plan;
    plan.sample = [descriptor, channel, displayModeId, anchor, protein, conformation](std::size_t frame) {
        auto gap = [](model::GapReason reason) {
            return model::FrameSignalSample::Gap(reason);
        };
        auto finish = [&](std::optional<double> raw) {
            if (!raw)
                return gap(model::GapReason::MalformedSource);
            return finiteSample(applyDisplayUnits(*raw, descriptor, channel));
        };

        const model::TrajectoryConformation* trajectory = conformation ? conformation->asTrajectory() : nullptr;
        const auto* h5 = trajectory ? trajectory->h5() : nullptr;
        if (!h5)
            return gap(model::GapReason::SourceAbsent);

        auto rowFor = [&](model::SignalAxis axis, std::size_t rows) -> std::optional<std::size_t> {
            if (rows > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                return std::nullopt;
            return anchorRow(anchor, axis, protein, static_cast<int>(rows));
        };
        auto rowForDescriptor = [&](std::size_t rows) -> std::optional<std::size_t> {
            return rowFor(descriptor.requiredAnchor, rows);
        };
        auto sourceMaskOff = [](const std::vector<uint8_t>& mask, std::size_t t) {
            return !mask.empty() && (t >= mask.size() || mask[t] == 0);
        };

        auto samplePositions = [&]() {
            const model::QtPositionsTimeSeries* ts = h5->positions();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(sampleVecValue(ts->at(*row, frame), channel, displayModeId));
        };

        auto sampleTensorSeries = [&](const model::QtShieldingTimeSeries* ts) {
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (!ts->sourceAttachedAt(frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(sampleTensorValue(ts->at(*row, frame), channel, displayModeId));
        };

        auto sampleT2Series = [&](const model::QtT2TimeSeries* ts) {
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (!ts->sourceAttachedAt(frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(sampleT2Value(ts->at(*row, frame), channel, displayModeId));
        };

        auto sampleScalarSeries = [&](const model::QtScalarTimeSeries* ts) {
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (!ts->sourceAttachedAt(frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(ts->at(*row, frame));
        };

        auto sampleVecSeries = [&](const model::QtVec3TimeSeries* ts) {
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (!ts->sourceAttachedAt(frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(sampleVecValue(ts->at(*row, frame), channel, displayModeId));
        };

        const QString path = descriptor.storagePath;

        if (path == QStringLiteral("/trajectory/positions"))
            return samplePositions();

        if (path == QStringLiteral("/trajectory/bs_shielding_time_series"))
            return sampleTensorSeries(h5->bsShielding());
        if (path == QStringLiteral("/trajectory/hm_shielding_time_series"))
            return sampleTensorSeries(h5->hmShielding());
        if (path == QStringLiteral("/trajectory/mc_shielding_time_series"))
            return sampleTensorSeries(h5->mcShielding());
        if (path == QStringLiteral("/trajectory/piquad_shielding_time_series"))
            return sampleTensorSeries(h5->piQuadShielding());
        if (path == QStringLiteral("/trajectory/ringchi_shielding_time_series"))
            return sampleTensorSeries(h5->ringChiShielding());
        if (path == QStringLiteral("/trajectory/disp_shielding_time_series"))
            return sampleTensorSeries(h5->dispShielding());
        if (path == QStringLiteral("/trajectory/hbond_shielding_time_series"))
            return sampleTensorSeries(h5->hbondShielding());
        if (path == QStringLiteral("/trajectory/mopac_coulomb_shielding_time_series"))
            return sampleT2Series(h5->mopacCoulombShielding());
        if (path == QStringLiteral("/trajectory/mopac_mc_shielding_time_series"))
            return sampleTensorSeries(h5->mopacMcShielding());
        if (path == QStringLiteral("/trajectory/tripeptide_bb_shielding_time_series"))
            return sampleTensorSeries(h5->tripeptideBbShielding());
        if (path == QStringLiteral("/trajectory/tripeptide_neighbor_shielding_time_series"))
            return sampleTensorSeries(h5->tripeptideNeighborShielding());
        if (path == QStringLiteral("/trajectory/larsen_hbond_1pHB_shielding_time_series"))
            return sampleTensorSeries(h5->larsenHBond1pHBShielding());
        if (path == QStringLiteral("/trajectory/larsen_hbond_1pHaB_shielding_time_series"))
            return sampleTensorSeries(h5->larsenHBond1pHaBShielding());
        if (path == QStringLiteral("/trajectory/larsen_hbond_2pHB_shielding_time_series"))
            return sampleTensorSeries(h5->larsenHBond2pHBShielding());
        if (path == QStringLiteral("/trajectory/larsen_hbond_2pHaB_shielding_time_series"))
            return sampleTensorSeries(h5->larsenHBond2pHaBShielding());

        if (path == QStringLiteral("/trajectory/sasa_time_series"))
            return sampleScalarSeries(h5->sasa());
        if (path == QStringLiteral("/trajectory/aimnet2_charge_time_series"))
            return sampleScalarSeries(h5->aimnet2Charge());
        if (path == QStringLiteral("/trajectory/larsen_hbond_count_time_series"))
            return sampleScalarSeries(h5->larsenHBondCount());
        if (path == QStringLiteral("/trajectory/larsen_hbond_water_term_time_series"))
            return sampleScalarSeries(h5->larsenHBondWaterTerm());
        if (path == QStringLiteral("/trajectory/bonded_energy_time_series"))
            return sampleScalarSeries(h5->bondedEnergyTotal());
        if (path == QStringLiteral("/trajectory/mopac_vs_ff14sb_reconciliation"))
            return sampleScalarSeries(h5->mopacVsFf14sbReconciliation());

        if (path == QStringLiteral("/trajectory/apbs_efield_time_series"))
            return sampleVecSeries(h5->apbsEfield());
        if (path == QStringLiteral("/trajectory/apbs_efg_time_series"))
            return sampleT2Series(h5->apbsEfg());
        if (path == QStringLiteral("/trajectory/tripeptide_bb_residual_vec_time_series"))
            return sampleVecSeries(h5->tripeptideBbResidualVec());
        if (path == QStringLiteral("/trajectory/tripeptide_neighbor_residual_vec_prev_time_series"))
            return sampleVecSeries(h5->tripeptideNeighborResidualVecPrev());
        if (path == QStringLiteral("/trajectory/tripeptide_neighbor_residual_vec_next_time_series"))
            return sampleVecSeries(h5->tripeptideNeighborResidualVecNext());

        if (path == QStringLiteral("/trajectory/water_field_time_series")) {
            const model::QtWaterFieldTimeSeries* ts = h5->waterFieldTimeSeries();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (sourceMaskOff(ts->source_attached, frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            if (descriptor.conceptKey == QStringLiteral("water_efield"))
                return finish(sampleVecValue(ts->efieldAt(*row, frame), channel, displayModeId));
            if (descriptor.conceptKey == QStringLiteral("water_efg"))
                return finish(sampleT2Value(ts->efgAt(*row, frame), channel, displayModeId));
            return finish(static_cast<double>(ts->nFirstAt(*row, frame)));
        }

        if (path == QStringLiteral("/trajectory/hydration_shell_time_series")) {
            const model::QtHydrationShellTimeSeries* ts = h5->hydrationShellTimeSeries();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (sourceMaskOff(ts->source_attached, frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(ts->halfShellAsymmetryAt(*row, frame));
        }

        if (path == QStringLiteral("/trajectory/hydration_geometry_time_series")) {
            const model::QtHydrationGeometryTimeSeries* ts = h5->hydrationGeometryTimeSeries();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (sourceMaskOff(ts->source_attached, frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(ts->dipole_alignment[*row * ts->n_frames + frame]);
        }

        if (path == QStringLiteral("/trajectory/aimnet2_embedding_time_series")) {
            const model::QtEmbeddingTimeSeries* ts = h5->aimnet2Embedding();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (sourceMaskOff(ts->meta.source_attached, frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            const float* values = ts->dataAt(*row, frame);
            if (!values || ts->n_dims == 0)
                return gap(model::GapReason::MalformedSource);
            return finish(static_cast<double>(values[0]));
        }

        if (path == QStringLiteral("/trajectory/aimnet2_charge_response_gradient_time_series")) {
            const model::QtAimnet2ChargeResponseGradientTimeSeries* ts = h5->aimnet2ChargeResponseGradient();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (sourceMaskOff(ts->meta.source_attached, frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(sampleVecValue(ts->vecAt(*row, frame), channel, displayModeId));
        }

        if (path == QStringLiteral("/trajectory/tripeptide_bb_method_tag_time_series")) {
            const model::QtTagTimeSeries* ts = h5->tripeptideBbMethodTag();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (sourceMaskOff(ts->meta.source_attached, frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(static_cast<double>(ts->at(*row, frame)));
        }

        if (path == QStringLiteral("/trajectory/dihedral_time_series")) {
            const model::QtDihedralTimeSeries* ts = h5->dihedrals();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (!ts->sourceAttachedAt(frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_residues);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            const QString id = channel.id.toLower();
            if (id == QStringLiteral("psi"))
                return finish(ts->psiAt(*row, frame));
            if (id == QStringLiteral("omega"))
                return finish(ts->omegaAt(*row, frame));
            if (id == QStringLiteral("chi"))
                return finish(ts->chiAt(*row, frame, 0));
            return finish(ts->phiAt(*row, frame));
        }

        if (path == QStringLiteral("/trajectory/dssp8_time_series")) {
            const model::QtDssp8TimeSeries* ts = h5->dssp8();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (!ts->sourceAttachedAt(frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_residues);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(static_cast<double>(static_cast<int>(ts->codeAt(*row, frame))));
        }

        if (path == QStringLiteral("/trajectory/j_coupling_time_series")) {
            const model::QtJCouplingTimeSeries* ts = h5->jCouplingTimeSeries();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (!ts->meta.sourceAttachedAt(frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_residues);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            if (*row < ts->J_Cprime_Cgamma_exists.size() && ts->J_Cprime_Cgamma_exists[*row] == 0)
                return gap(model::GapReason::NotApplicable);
            return finish(ts->at(ts->J_Cprime_Cgamma, *row, frame));
        }

        if (path == QStringLiteral("/trajectory/ring_pucker_time_series")) {
            const model::QtRingPuckerTimeSeries* ts = h5->ringPucker();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (sourceMaskOff(ts->source_attached, frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_aromatic_rings + ts->n_saturated_rings);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            if (*row < ts->n_aromatic_rings)
                return finish(ts->aromaticChi2At(*row, frame));
            const std::size_t saturated = *row - ts->n_aromatic_rings;
            if (saturated < ts->n_saturated_rings)
                return finish(ts->puckerQAt(saturated, frame));
            return gap(model::GapReason::AnchorUnavailable);
        }

        if (path == QStringLiteral("/trajectory/ring_neighbourhood_trajectory_stats")) {
            const model::QtRingNeighbourhoodTimeSeries* ts = h5->ringNeighbourhood();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (sourceMaskOff(ts->source_attached, frame))
                return gap(model::GapReason::SourceMaskOff);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            if (ts->n_slots == 0 || ts->ringIndexAt(*row, 0) < 0)
                return gap(model::GapReason::NotApplicable);
            const std::array<double, 4> values = ts->at(*row, frame, 0);
            const QString id = channel.id.toLower();
            if (id == QStringLiteral("rho"))
                return finish(values[1]);
            if (id == QStringLiteral("z"))
                return finish(values[2]);
            if (id == QStringLiteral("in_plane_angle"))
                return finish(values[3]);
            return finish(values[0]);
        }

        if (path == QStringLiteral("/trajectory/gromacs_energy_time_series")) {
            const model::QtSystemEnergyTimeSeries* ts = h5->gromacsEnergy();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (!ts->sourceAttachedAt(frame))
                return gap(model::GapReason::SourceMaskOff);
            const QString id = channel.id.toLower();
            if (id == QStringLiteral("temperature"))
                return finish(frame < ts->temperature.size() ? std::optional<double>(ts->temperature[frame]) : std::nullopt);
            if (id == QStringLiteral("pressure"))
                return finish(frame < ts->pressure.size() ? std::optional<double>(ts->pressure[frame]) : std::nullopt);
            if (id == QStringLiteral("volume"))
                return finish(frame < ts->volume.size() ? std::optional<double>(ts->volume[frame]) : std::nullopt);
            return finish(frame < ts->total_energy.size() ? std::optional<double>(ts->total_energy[frame]) : std::nullopt);
        }

        if (path == QStringLiteral("/trajectory/rmsd_tracking")) {
            const model::QtRmsdTracking* ts = h5->rmsdTracking();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            if (frame >= ts->n_frames)
                return gap(model::GapReason::FrameSourceAbsent);
            if (sourceMaskOff(ts->source_attached, frame))
                return gap(model::GapReason::SourceMaskOff);
            return finish(frame < ts->rmsd.size() ? std::optional<double>(ts->rmsd[frame]) : std::nullopt);
        }

        auto sampleShieldingRollup = [&](const model::QtShieldingWelford* ts) {
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            if (*row >= ts->t0.size())
                return gap(model::GapReason::MalformedSource);
            return finish(ts->t0[*row].mean);
        };
        auto sampleScalarRollup = [&](const model::QtScalarWelford* ts) {
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            if (*row >= ts->value.size())
                return gap(model::GapReason::MalformedSource);
            return finish(ts->value[*row].mean);
        };
        auto sampleVecRollup = [&](const model::QtVec3Welford* ts) {
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            if (*row >= ts->magnitude.size())
                return gap(model::GapReason::MalformedSource);
            return finish(ts->magnitude[*row].mean);
        };

        if (path == QStringLiteral("/trajectory/bond_length_stats")) {
            const model::QtBondLengthStats* ts = h5->bondLengthStats();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_bonds);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(*row < ts->length_mean.size() ? std::optional<double>(ts->length_mean[*row]) : std::nullopt);
        }
        if (path == QStringLiteral("/trajectory/bs_welford"))
            return sampleShieldingRollup(h5->bsWelford());
        if (path == QStringLiteral("/trajectory/hm_welford"))
            return sampleShieldingRollup(h5->hmWelford());
        if (path == QStringLiteral("/trajectory/mc_welford"))
            return sampleShieldingRollup(h5->mcWelford());
        if (path == QStringLiteral("/trajectory/sasa_welford"))
            return sampleScalarRollup(h5->sasaWelford());
        if (path == QStringLiteral("/trajectory/eeq_welford"))
            return sampleScalarRollup(h5->eeqWelford());
        if (path == QStringLiteral("/trajectory/hbond_count_welford"))
            return sampleScalarRollup(h5->hbondCountWelford());
        if (path == QStringLiteral("/trajectory/mopac_charge_welford"))
            return sampleScalarRollup(h5->mopacChargeWelford());
        if (path == QStringLiteral("/trajectory/mopac_bond_order_welford")) {
            const model::QtBondOrderWelford* ts = h5->mopacBondOrderWelford();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_bonds);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            if (*row >= ts->bond_order.size())
                return gap(model::GapReason::MalformedSource);
            return finish(ts->bond_order[*row].mean);
        }
        if (path == QStringLiteral("/trajectory/water_field_welford"))
            return sampleVecRollup(h5->waterFieldWelford());
        if (path == QStringLiteral("/trajectory/aimnet2_charge_response_gradient_welford"))
            return sampleVecRollup(h5->aimnet2ChargeResponseGradientWelford());
        if (path == QStringLiteral("/trajectory/hydration_shell_welford")
            || path == QStringLiteral("/trajectory/hydration_geometry_welford")) {
            const model::QtHydrationWelford* ts = path == QStringLiteral("/trajectory/hydration_shell_welford")
                                                     ? h5->hydrationShellWelford()
                                                     : h5->hydrationGeometryWelford();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            if (ts->channels.empty() || *row >= ts->channels.front().moments.size())
                return gap(model::GapReason::MalformedSource);
            return finish(ts->channels.front().moments[*row].mean);
        }
        if (path == QStringLiteral("/trajectory/bs_t0_autocorrelation")) {
            const model::QtAutocorrelation* ts = h5->bsT0Autocorrelation();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(ts->at(*row, 0));
        }

        if (path == QStringLiteral("/trajectory/dssp8_transition")) {
            const model::QtDssp8Transitions* ts = h5->dssp8Transitions();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_residues);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(*row < ts->ss8_transition_count.size()
                              ? std::optional<double>(static_cast<double>(ts->ss8_transition_count[*row]))
                              : std::nullopt);
        }
        if (path == QStringLiteral("/trajectory/dihedral_bin_transition")) {
            const model::QtDihedralBinTransitions* ts = h5->dihedralBinTransitions();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowForDescriptor(ts->n_residues);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(*row < ts->backbone_transition_count.size()
                              ? std::optional<double>(static_cast<double>(ts->backbone_transition_count[*row]))
                              : std::nullopt);
        }

        return gap(model::GapReason::Pending);
    };
    return plan;
}

SamplePlan dftPlan(const model::SignalDescriptor& descriptor,
                   const QString& displayModeId,
                   const model::SignalAnchor& anchor,
                   const model::QtProtein* protein,
                   const QPointer<model::Conformation>& conformation,
                   const QPointer<model::DftShieldingStore>& dftStore) {
    const std::optional<std::size_t> atom = atomFromAnchor(anchor, protein);
    if (!atom)
        return pendingPlan();

    SamplePlan plan;
    plan.needsDftFrame = true;

    model::DftPart part = model::DftPart::Total;
    if (descriptor.conceptKey.contains(QStringLiteral("diamagnetic"), Qt::CaseInsensitive))
        part = model::DftPart::Dia;
    else if (descriptor.conceptKey.contains(QStringLiteral("paramagnetic"), Qt::CaseInsensitive))
        part = model::DftPart::Para;

    const model::DftScalar scalar = displayModeId == QStringLiteral("strip.tensor.T2")
                                        ? model::DftScalar::AnisotropyT2
                                        : model::DftScalar::IsotropicT0;

    plan.sample = [conformation, dftStore, atom = *atom, part, scalar](std::size_t frame) {
        if (!conformation || !dftStore)
            return model::FrameSignalSample::Gap(model::GapReason::SourceAbsent);
        const std::size_t original = conformation->originalFrameIndex(frame);
        if (!dftStore->hasJob(original))
            return model::FrameSignalSample::Gap(model::GapReason::FrameSourceAbsent);
        const std::optional<double> value = dftStore->sample(original, atom, part, scalar);
        if (!value)
            return model::FrameSignalSample::Gap(model::GapReason::FrameSourceAbsent);
        if (!std::isfinite(*value))
            return model::FrameSignalSample::Gap(model::GapReason::NaNSentinel);
        return model::FrameSignalSample::Valid(*value);
    };
    return plan;
}

SamplePlan geometryPlan(const model::SignalDescriptor& descriptor,
                        const model::SignalAnchor& anchor,
                        const model::QtProtein* protein,
                        const QPointer<model::Conformation>& conformation) {
    if (descriptor.storagePath == QStringLiteral("atom_displacement")) {
        const std::optional<std::size_t> atom = atomFromAnchor(anchor, protein);
        if (!atom)
            return pendingPlan();

        SamplePlan plan;
        plan.sample = [conformation, atom = *atom](std::size_t frame) {
            if (!conformation)
                return model::FrameSignalSample::Gap(model::GapReason::SourceAbsent);
            if (frame >= conformation->frameCount())
                return model::FrameSignalSample::Gap(model::GapReason::FrameSourceAbsent);
            const model::Vec3 delta = conformation->atomPosition(frame, atom) - conformation->atomPosition(0, atom);
            return finiteSample(delta.norm());
        };
        return plan;
    }

    const std::optional<std::vector<std::size_t>> atoms = atomTupleFromAnchor(anchor, protein);
    if (!atoms)
        return pendingPlan();

    SamplePlan plan;
    plan.sample = [descriptor, conformation, atoms = *atoms](std::size_t frame) {
        if (!conformation)
            return model::FrameSignalSample::Gap(model::GapReason::SourceAbsent);
        const model::GeometryMeasurement measurement = model::Measure(*conformation, frame, atoms);
        if (!measurement.valid)
            return model::FrameSignalSample::Gap(model::GapReason::AnchorUnavailable);
        if (descriptor.storagePath == QStringLiteral("distance")
            && measurement.kind != model::GeometryKind::Distance) {
            return model::FrameSignalSample::Gap(model::GapReason::NotApplicable);
        }
        if (descriptor.storagePath == QStringLiteral("angle")
            && measurement.kind != model::GeometryKind::Angle) {
            return model::FrameSignalSample::Gap(model::GapReason::NotApplicable);
        }
        if (descriptor.storagePath == QStringLiteral("dihedral")
            && measurement.kind != model::GeometryKind::Dihedral) {
            return model::FrameSignalSample::Gap(model::GapReason::NotApplicable);
        }
        return model::FrameSignalSample::Valid(measurement.value);
    };
    return plan;
}

SamplePlan topologyPlan(const model::SignalDescriptor& descriptor,
                        const model::SignalAnchor& anchor,
                        const model::QtProtein* protein,
                        const QPointer<model::Conformation>& conformation) {
    if (descriptor.storagePath != QStringLiteral("bond_length"))
        return pendingPlan();
    const auto* bond = std::get_if<model::BondAnchor>(&anchor);
    if (!bond)
        return pendingPlan();

    SamplePlan plan;
    plan.sample = [protein, conformation, bondIndex = bond->bond](std::size_t frame) {
        if (!protein || !conformation)
            return model::FrameSignalSample::Gap(model::GapReason::SourceAbsent);
        if (bondIndex >= protein->bondCount())
            return model::FrameSignalSample::Gap(model::GapReason::AnchorUnavailable);
        if (frame >= conformation->frameCount())
            return model::FrameSignalSample::Gap(model::GapReason::FrameSourceAbsent);
        const model::QtBond& bond = protein->bond(bondIndex);
        if (bond.atomIndexA < 0 || bond.atomIndexB < 0)
            return model::FrameSignalSample::Gap(model::GapReason::MalformedSource);
        const auto atomA = static_cast<std::size_t>(bond.atomIndexA);
        const auto atomB = static_cast<std::size_t>(bond.atomIndexB);
        if (atomA >= protein->atomCount() || atomB >= protein->atomCount())
            return model::FrameSignalSample::Gap(model::GapReason::MalformedSource);
        const model::Vec3 delta = conformation->atomPosition(frame, atomA) - conformation->atomPosition(frame, atomB);
        return finiteSample(delta.norm());
    };
    return plan;
}

SamplePlan selectionEventsPlan(const model::SignalDescriptor& descriptor,
                               const QPointer<model::Conformation>& conformation,
                               const QPointer<model::AtomSelection>& selection) {
    SamplePlan plan;
    plan.sample = [descriptor, conformation, selection](std::size_t frame) {
        if (descriptor.storagePath == QStringLiteral("/trajectory/selections")) {
            const model::TrajectoryConformation* trajectory = conformation ? conformation->asTrajectory() : nullptr;
            const auto* h5 = trajectory ? trajectory->h5() : nullptr;
            if (!h5)
                return model::FrameSignalSample::Gap(model::GapReason::SourceAbsent);
            const std::size_t originalFrame = conformation ? conformation->originalFrameIndex(frame) : frame;
            int count = 0;
            for (const model::QtSelectionEvent& event : h5->selections().events()) {
                if (static_cast<std::size_t>(event.frame_idx) == originalFrame)
                    ++count;
            }
            return model::FrameSignalSample::Valid(static_cast<double>(count));
        }

        if (descriptor.storagePath == QStringLiteral("selection_timeline")
            || descriptor.storagePath == QStringLiteral("selection_counts")) {
            return model::FrameSignalSample::Valid(selection ? static_cast<double>(selection->count()) : 0.0);
        }

        return model::FrameSignalSample::Gap(model::GapReason::Pending);
    };
    return plan;
}

SamplePlan samplePlanFor(const model::DashboardSignal& signal,
                         const model::SignalDescriptor& descriptor,
                         const model::ChannelDescriptor& channel,
                         const QString& displayModeId,
                         const model::SignalAnchor& anchor,
                         const model::QtProtein* protein,
                         const QPointer<model::Conformation>& conformation,
                         const QPointer<model::DftShieldingStore>& dftStore,
                         const QPointer<model::AtomSelection>& selection) {
    (void)signal;
    switch (descriptor.sourceKind) {
    case model::SignalSourceKind::DenseH5Trajectory:
        return denseH5Plan(descriptor, channel, displayModeId, anchor, protein, conformation);
    case model::SignalSourceKind::FrameNpySnapshot:
        return frameNpyPlan(descriptor, channel, displayModeId, anchor, protein, conformation);
    case model::SignalSourceKind::OrcaDftFrame:
        return dftPlan(descriptor, displayModeId, anchor, protein, conformation, dftStore);
    case model::SignalSourceKind::DerivedGeometry:
        return geometryPlan(descriptor, anchor, protein, conformation);
    case model::SignalSourceKind::Topology:
        return topologyPlan(descriptor, anchor, protein, conformation);
    case model::SignalSourceKind::SelectionEvents:
        return selectionEventsPlan(descriptor, conformation, selection);
    }
    return pendingPlan();
}
}  // namespace

DashboardDisplayController::DashboardDisplayController(QObject* parent)
    : QObject(parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("DashboardDisplayController"));
}

void DashboardDisplayController::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    ASSERT_THREAD(this);
    protein_ = protein;
    conformation_ = conformation;
    rebuild();
}

void DashboardDisplayController::setSignalModels(model::TrajectorySignalCatalog* catalog,
                                                 model::DashboardSignalModel* activeModel) {
    ASSERT_THREAD(this);
    if (activeModel_)
        disconnect(activeModel_, nullptr, this, nullptr);
    catalog_ = catalog;
    activeModel_ = activeModel;
    if (activeModel_) {
        ACONNECT(activeModel_.data(), &model::DashboardSignalModel::signalAdded,
                 this, [this](const QUuid&) { rebuild(); });
        ACONNECT(activeModel_.data(), &model::DashboardSignalModel::signalRemoved,
                 this, [this](const QUuid&) { rebuild(); });
        ACONNECT(activeModel_.data(), &model::DashboardSignalModel::signalChanged,
                 this, [this](const QUuid&) { rebuild(); });
        ACONNECT(activeModel_.data(), &QAbstractItemModel::modelReset,
                 this, &DashboardDisplayController::rebuild);
    }
    rebuild();
}

void DashboardDisplayController::setPanelModel(model::DashboardPanelModel* panelModel) {
    ASSERT_THREAD(this);
    if (panelModel_)
        disconnect(panelModel_, nullptr, this, nullptr);
    panelModel_ = panelModel;
    if (panelModel_) {
        ACONNECT(panelModel_.data(), &model::DashboardPanelModel::activePanelChanged,
                 this, [this](const QUuid&) { refreshPanelVisibility(); });
        ACONNECT(panelModel_.data(), &model::DashboardPanelModel::displayRefsChanged,
                 this, [this](const QUuid&) { refreshPanelVisibility(); });
        ACONNECT(panelModel_.data(), &model::DashboardPanelModel::panelAdded,
                 this, [this](const QUuid&) { refreshPanelVisibility(); });
        ACONNECT(panelModel_.data(), &model::DashboardPanelModel::panelRemoved,
                 this, [this](const QUuid&, const QVector<model::DashboardDisplayRef>&) {
                     refreshPanelVisibility();
                 });
        ACONNECT(panelModel_.data(), &QAbstractItemModel::modelReset,
                 this, &DashboardDisplayController::refreshPanelVisibility);
    }
    refreshPanelVisibility();
}

void DashboardDisplayController::setSelection(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    if (selection_)
        disconnect(selection_, nullptr, this, nullptr);
    selection_ = selection;
    if (selection_) {
        ACONNECT(selection_.data(), &model::AtomSelection::changed,
                 this, &DashboardDisplayController::rebuild);
        ACONNECT(selection_.data(), &model::AtomSelection::focusChanged,
                 this, [this](std::size_t) { rebuild(); });
        ACONNECT(selection_.data(), &model::AtomSelection::cleared,
                 this, &DashboardDisplayController::rebuild);
    }
    rebuild();
}

void DashboardDisplayController::setDftStore(model::DftShieldingStore* store) {
    ASSERT_THREAD(this);
    dftStore_ = store;
    rebuild();
}

void DashboardDisplayController::setFrame(int frame) {
    ASSERT_THREAD(this);
    frame_ = std::max(0, frame);
    extendToFrame(frame_);
    emit stripTracksChanged();
}

QVector<DashboardDisplayController::StripTrack> DashboardDisplayController::stripTracks() const {
    QVector<StripTrack> out;
    out.reserve(series_.size());
    for (const ActiveSeries& series : series_) {
        if (!seriesIsVisibleInActivePanel(series))
            continue;
        StripTrack item;
        item.buffer = &series.buffer.channel;
        item.color = series.color;
        item.hasBinding = series.hasBinding;
        item.binding = series.binding;
        out.push_back(item);
    }
    return out;
}

DashboardSmokeSummary DashboardDisplayController::smokeSummary() const {
    int lastFrame = -1;
    for (const ActiveSeries& series : series_) {
        if (!series.buffer.statuses.empty()) {
            lastFrame = std::max(lastFrame,
                                 static_cast<int>(series.buffer.statuses.size() - 1));
        }
    }
    return smokeSummary(0, lastFrame);
}

DashboardSmokeSummary DashboardDisplayController::smokeSummary(int firstFrame, int lastFrame) const {
    DashboardSmokeSummary summary;
    summary.seriesCount = static_cast<int>(series_.size());
    summary.seriesSparseness.reserve(series_.size());

    const int first = std::max(0, firstFrame);
    const bool hasWindow = lastFrame >= first;
    const std::size_t begin = hasWindow ? static_cast<std::size_t>(first) : 0;
    const std::size_t end = hasWindow ? static_cast<std::size_t>(lastFrame) + 1 : 0;
    auto rangeCount = [begin, end](std::size_t size) -> long long {
        if (end <= begin || size <= begin)
            return 0;
        return static_cast<long long>(std::min(size, end) - begin);
    };

    for (const ActiveSeries& series : series_) {
        const std::vector<model::SampleStatus>& statuses = series.buffer.statuses;
        const std::vector<model::GapReason>& gapReasons = series.buffer.gapReasons;
        DashboardSmokeSummary::SeriesSparseness sparseness;
        sparseness.signalLabel = series.signal.label;
        sparseness.descriptorId = series.descriptor.id;
        sparseness.conceptKey = series.descriptor.conceptKey;
        sparseness.sourceKind = model::ToString(series.descriptor.sourceKind);
        sparseness.storagePath = series.descriptor.storagePath;
        sparseness.displayModeId = series.displayModeId;
        sparseness.channelId = series.channel.id;
        sparseness.channelLabel = series.channel.label;
        sparseness.samples = rangeCount(statuses.size());

        summary.channelValues += rangeCount(series.buffer.channel.values.size());
        summary.channelValidity += rangeCount(series.buffer.channel.valid.size());
        const bool bufferSizeMismatch =
            series.buffer.channel.values.size() != statuses.size()
            || series.buffer.channel.valid.size() != statuses.size()
            || gapReasons.size() != statuses.size();
        const bool requestedWindowMissing =
            hasWindow
            && (statuses.size() < end
                || series.buffer.channel.values.size() < end
                || series.buffer.channel.valid.size() < end
                || gapReasons.size() < end);
        if (bufferSizeMismatch || requestedWindowMissing) {
            ++summary.seriesWithMismatchedBuffers;
        }
        if (sparseness.samples > 0)
            ++summary.seriesWithSamples;

        bool hasValid = false;
        bool pendingOnly = sparseness.samples > 0;
        summary.samples += sparseness.samples;
        int currentValidRun = 0;
        int currentGapRun = 0;

        const std::size_t loopEnd = hasWindow ? std::min(statuses.size(), end) : 0;
        for (std::size_t i = begin; i < loopEnd; ++i) {
            const model::SampleStatus status = statuses[i];
            const model::GapReason reason = i < gapReasons.size() ? gapReasons[i] : model::GapReason::None;

            switch (status) {
            case model::SampleStatus::Valid:
                ++summary.validSamples;
                ++sparseness.validSamples;
                hasValid = true;
                pendingOnly = false;
                if (sparseness.firstValidFrame < 0)
                    sparseness.firstValidFrame = static_cast<int>(i);
                sparseness.lastValidFrame = static_cast<int>(i);
                ++currentValidRun;
                currentGapRun = 0;
                sparseness.longestValidRun = std::max(sparseness.longestValidRun, currentValidRun);
                break;
            case model::SampleStatus::Gap:
                ++summary.gapSamples;
                ++sparseness.gapSamples;
                ++currentGapRun;
                currentValidRun = 0;
                sparseness.longestGapRun = std::max(sparseness.longestGapRun, currentGapRun);
                break;
            case model::SampleStatus::NotAvailable:
                ++summary.gapSamples;
                ++sparseness.gapSamples;
                pendingOnly = false;
                ++currentGapRun;
                currentValidRun = 0;
                sparseness.longestGapRun = std::max(sparseness.longestGapRun, currentGapRun);
                break;
            case model::SampleStatus::Invalid:
                ++summary.invalidSamples;
                ++sparseness.invalidSamples;
                pendingOnly = false;
                ++currentGapRun;
                currentValidRun = 0;
                sparseness.longestGapRun = std::max(sparseness.longestGapRun, currentGapRun);
                break;
            }

            switch (reason) {
            case model::GapReason::Pending:
                ++summary.pendingGapSamples;
                ++sparseness.pendingGapSamples;
                break;
            case model::GapReason::SourceAbsent:
                ++summary.sourceAbsentGapSamples;
                ++sparseness.sourceAbsentGapSamples;
                pendingOnly = false;
                break;
            case model::GapReason::FrameSourceAbsent:
                ++summary.frameSourceAbsentGapSamples;
                ++sparseness.frameSourceAbsentGapSamples;
                pendingOnly = false;
                break;
            case model::GapReason::AnchorUnavailable:
                ++summary.anchorUnavailableGapSamples;
                ++sparseness.anchorUnavailableGapSamples;
                pendingOnly = false;
                break;
            case model::GapReason::SourceMaskOff:
                ++sparseness.sourceMaskOffGapSamples;
                pendingOnly = false;
                break;
            case model::GapReason::NotApplicable:
                ++sparseness.notApplicableGapSamples;
                pendingOnly = false;
                break;
            case model::GapReason::NaNSentinel:
                ++sparseness.nanSentinelGapSamples;
                pendingOnly = false;
                break;
            case model::GapReason::MalformedSource:
                ++sparseness.malformedSourceGapSamples;
                pendingOnly = false;
                break;
            case model::GapReason::None:
                pendingOnly = false;
                break;
            }
        }

        if (hasValid)
            ++summary.seriesWithValidSamples;
        if (pendingOnly)
            ++summary.seriesPendingOnly;
        if (sparseness.samples > 0 && sparseness.validSamples == sparseness.samples)
            ++summary.denseSeries;
        if (sparseness.samples > 0 && sparseness.validSamples == 0)
            ++summary.allGapSeries;
        if (sparseness.validSamples > 0 && sparseness.gapSamples > 0)
            ++summary.sparseSeries;
        if (sparseness.frameSourceAbsentGapSamples > 0) {
            ++summary.seriesWithFrameSourceAbsentGaps;
            if (series.descriptor.sourceKind == model::SignalSourceKind::FrameNpySnapshot) {
                ++summary.frameNpySeriesWithFrameSourceAbsentGaps;
                summary.frameNpyFrameSourceAbsentGapSamples += sparseness.frameSourceAbsentGapSamples;
            } else if (series.descriptor.sourceKind == model::SignalSourceKind::OrcaDftFrame) {
                ++summary.orcaDftSeriesWithFrameSourceAbsentGaps;
                summary.orcaDftFrameSourceAbsentGapSamples += sparseness.frameSourceAbsentGapSamples;
            }
        }
        if (sparseness.sourceAbsentGapSamples > 0)
            ++summary.seriesWithSourceAbsentGaps;
        if (sparseness.anchorUnavailableGapSamples > 0)
            ++summary.seriesWithAnchorUnavailableGaps;
        summary.maxLongestGapRun = std::max(summary.maxLongestGapRun, sparseness.longestGapRun);
        summary.seriesSparseness.push_back(std::move(sparseness));
    }

    return summary;
}

void DashboardDisplayController::rebuild() {
    ASSERT_THREAD(this);

    QVector<ActiveSeries> next;
    activeStripSignalCount_ = 0;

    if (catalog_ && activeModel_) {
        for (const model::DashboardSignal& signal : activeModel_->activeSignals()) {
            if (!signal.enabled || !hasStripMode(signal.displayModeIds))
                continue;
            ++activeStripSignalCount_;
            const model::SignalDescriptor* descriptor = catalog_->findDescriptor(signal.binding.descriptorId);
            if (!descriptor)
                continue;
            buildGenericTracks(signal, *descriptor, next);
        }
    }

    for (int i = 0; i < next.size(); ++i)
        next[i].color = colorForIndex(i);

    series_ = std::move(next);
    extendToFrame(frame_);

    updateStatusText();
    emit stripTracksChanged();
}

void DashboardDisplayController::refreshPanelVisibility() {
    ASSERT_THREAD(this);
    updateStatusText();
    emit stripTracksChanged();
}

void DashboardDisplayController::updateStatusText() {
    if (!catalog_ || !activeModel_) {
        statusText_ = QStringLiteral("Dashboard signal model is not connected.");
    } else if (activeStripSignalCount_ == 0) {
        statusText_ = QStringLiteral("No active strip display modes.");
    } else {
        const int displayedSeries = activePanelSeriesCount();
        statusText_ = QStringLiteral("%1 displayed signal time series from %2 active strip signal%3.")
                          .arg(displayedSeries)
                          .arg(activeStripSignalCount_)
                          .arg(activeStripSignalCount_ == 1 ? QString() : QStringLiteral("s"));
        statusText_ += QStringLiteral(" Unimplemented sources append explicit pending gaps.");
    }
}

void DashboardDisplayController::buildGenericTracks(const model::DashboardSignal& signal,
                                                    const model::SignalDescriptor& descriptor,
                                                    QVector<ActiveSeries>& series) const {
    QSet<QString> emitted;
    const model::SignalAnchor anchor = resolvedAnchorForSignal(signal, descriptor);
    const model::SignalBinding reveal = bindingFromAnchor(descriptor, anchor, signal.binding.followsFocus);

    for (const QString& mode : signal.displayModeIds) {
        if (!isStripMode(mode))
            continue;
        const QVector<model::ChannelDescriptor> channels = channelsForMode(descriptor, mode);
        for (const model::ChannelDescriptor& channel : channels) {
            const QString key = QStringLiteral("%1|%2|%3").arg(descriptor.id, mode, channel.id);
            if (emitted.contains(key))
                continue;
            emitted.insert(key);

            ActiveSeries active;
            active.signal = signal;
            active.descriptor = descriptor;
            active.channel = channel;
            active.displayModeId = mode;
            active.buffer.key.signalId = signal.id;
            active.buffer.key.descriptorId = descriptor.id;
            active.buffer.key.displayModeId = mode;
            active.buffer.key.channelId = channel.id;
            active.buffer.channel.id = active.buffer.key.stableId();
            active.buffer.channel.label = channelLabel(signal, descriptor, channel, mode);
            active.buffer.channel.unit = unitsLabel(descriptor, channel);
            SamplePlan plan = samplePlanFor(signal,
                                            descriptor,
                                            channel,
                                            mode,
                                            anchor,
                                            protein_,
                                            conformation_,
                                            dftStore_,
                                            selection_);
            active.sample = std::move(plan.sample);
            active.needsFrameSnapshot = plan.needsFrameSnapshot;
            active.needsDftFrame = plan.needsDftFrame;
            active.hasBinding = bindingHasRevealTarget(reveal);
            active.binding = reveal;
            series.push_back(std::move(active));
        }
    }
}

QVector<model::ChannelDescriptor>
DashboardDisplayController::channelsForMode(const model::SignalDescriptor& descriptor,
                                            const QString& displayModeId) const {
    QVector<model::ChannelDescriptor> channels;

    if (descriptor.channels.isEmpty()) {
        model::ChannelDescriptor channel;
        channel.id = QStringLiteral("value");
        channel.label = descriptor.valueShape == model::SignalValueShape::EventRecord
                            ? QStringLiteral("Event")
                            : QStringLiteral("Value");
        channel.valueShape = descriptor.valueShape;
        channel.sourceUnits = descriptor.sourceUnits;
        channel.defaultDisplayUnits = descriptor.defaultDisplayUnits;
        channels.push_back(channel);
        return channels;
    }

    for (const model::ChannelDescriptor& channel : descriptor.channels) {
        if (modeWantsChannel(displayModeId, channel))
            channels.push_back(channel);
    }
    if (channels.isEmpty())
        channels = descriptor.channels;
    return channels;
}

model::SignalAnchor
DashboardDisplayController::resolvedAnchorForSignal(const model::DashboardSignal& signal,
                                                    const model::SignalDescriptor& descriptor) const {
    if (!signal.binding.followsFocus)
        return signal.binding.anchor;

    switch (descriptor.requiredAnchor) {
    case model::SignalAxis::Atom:
        if (selection_ && selection_->hasFocus())
            return model::AtomAnchor{selection_->focus()};
        return model::NoneAnchor{};
    case model::SignalAxis::Residue:
        if (selection_ && selection_->hasFocus() && protein_) {
            const std::size_t atom = selection_->focus();
            if (atom < protein_->atomCount()) {
                const int residue = protein_->atom(atom).residueIndex;
                if (residue >= 0)
                    return model::ResidueAnchor{static_cast<std::size_t>(residue)};
            }
        }
        return model::NoneAnchor{};
    case model::SignalAxis::AtomTuple:
        if (selection_ && !selection_->atoms().empty())
            return model::AtomTupleAnchor{selection_->atoms()};
        return model::NoneAnchor{};
    case model::SignalAxis::None:
    case model::SignalAxis::Bond:
    case model::SignalAxis::Ring:
    case model::SignalAxis::AromaticRing:
    case model::SignalAxis::SaturatedRing:
    case model::SignalAxis::RingContributionPair:
    case model::SignalAxis::RingMembership:
    case model::SignalAxis::MutationMatchPair:
    case model::SignalAxis::Protein:
    case model::SignalAxis::System:
    case model::SignalAxis::Event:
        break;
    }

    return signal.binding.anchor;
}

QString DashboardDisplayController::channelLabel(const model::DashboardSignal& signal,
                                                 const model::SignalDescriptor& descriptor,
                                                 const model::ChannelDescriptor& channel,
                                                 const QString& displayModeId) const {
    QString label = signal.label.isEmpty() ? descriptor.label : signal.label;
    if (label.isEmpty())
        label = descriptor.id;
    if (!channel.label.isEmpty() && channel.id != QStringLiteral("value"))
        label += QStringLiteral(" / %1").arg(channel.label);
    label += QStringLiteral(" [%1]").arg(displayModeId.mid(QStringLiteral("strip.").size()));

    const model::SignalAnchor anchor = resolvedAnchorForSignal(signal, descriptor);
    const QString anchorText = model::AnchorLabel(anchor);
    if (!anchorText.isEmpty() && anchorText != QStringLiteral("No anchor"))
        label += QStringLiteral(" — %1").arg(anchorText);
    return label;
}

QString DashboardDisplayController::unitsLabel(const model::SignalDescriptor& descriptor,
                                               const model::ChannelDescriptor& channel) const {
    if (!channel.defaultDisplayUnits.displaySymbol.isEmpty())
        return channel.defaultDisplayUnits.displaySymbol;
    if (!channel.sourceUnits.displaySymbol.isEmpty())
        return channel.sourceUnits.displaySymbol;
    if (!descriptor.defaultDisplayUnits.displaySymbol.isEmpty())
        return descriptor.defaultDisplayUnits.displaySymbol;
    if (!descriptor.sourceUnits.displaySymbol.isEmpty())
        return descriptor.sourceUnits.displaySymbol;
    if (!channel.sourceUnits.sourceSymbol.isEmpty())
        return channel.sourceUnits.sourceSymbol;
    return descriptor.sourceUnits.sourceSymbol;
}

bool DashboardDisplayController::seriesIsVisibleInActivePanel(const ActiveSeries& series) const {
    if (!panelModel_)
        return true;
    const model::DashboardPanel* panel = panelModel_->activePanel();
    if (!panel)
        return false;
    const model::DashboardDisplayRef ref{
        series.signal.id,
        series.displayModeId,
        series.channel.id,
    };
    return panel->displays.contains(ref);
}

int DashboardDisplayController::activePanelSeriesCount() const {
    return static_cast<int>(std::count_if(series_.begin(), series_.end(), [this](const ActiveSeries& series) {
        return seriesIsVisibleInActivePanel(series);
    }));
}

void DashboardDisplayController::extendToFrame(int frame) {
    if (frame < 0)
        return;
    const bool needsSnapshot = std::any_of(series_.begin(), series_.end(), [](const ActiveSeries& series) {
        return series.needsFrameSnapshot;
    });
    const bool needsDft = std::any_of(series_.begin(), series_.end(), [](const ActiveSeries& series) {
        return series.needsDftFrame;
    });

    const long long startFrame = [&]() {
        long long start = frame;
        for (const ActiveSeries& series : series_)
            start = std::min(start, series.buffer.lastFrame() + 1);
        return std::max<long long>(0, start);
    }();

    for (long long f = startFrame; f <= frame; ++f) {
        const std::size_t sampleFrame = static_cast<std::size_t>(f);
        if (needsSnapshot && conformation_)
            conformation_->requestSnapshot(sampleFrame);
        if (needsDft && conformation_ && dftStore_)
            dftStore_->requestFrame(conformation_->originalFrameIndex(sampleFrame));

        for (ActiveSeries& series : series_) {
            if (series.buffer.lastFrame() >= f)
                continue;
            if (series.sample)
                series.buffer.append(series.sample(sampleFrame));
            else
                series.buffer.append(model::FrameSignalSample::Gap(model::GapReason::Pending));
        }
    }
}

QColor DashboardDisplayController::colorForIndex(int index) const {
    static const QColor colors[] = {
        QColor(86, 166, 244),
        QColor(255, 179, 87),
        QColor(120, 184, 92),
        QColor(180, 131, 230),
        QColor(229, 99, 99),
        QColor(94, 170, 220),
        QColor(215, 190, 85),
    };
    return colors[static_cast<std::size_t>(index) % (sizeof(colors) / sizeof(colors[0]))];
}

}  // namespace h5reader::app

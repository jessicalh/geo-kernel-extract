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
#include "../model/QtBondVectorBuffers.h"
#include "../model/QtPerAtomChannelBuffers.h"
#include "../model/TrajectoryConformation.h"
#include "../model/TrajectorySignalCatalog.h"
#include "ChordCouplingPanel.h"
#include "FixedFreqPanel.h"
#include "LagDecayPanel.h"
#include "PowerSpectrumPanel.h"
#include "SceneRevealOverlay.h"
#include "SequenceBarPanel.h"

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

struct ScopedScrubReleaseFlag {
    bool& flag;

    explicit ScopedScrubReleaseFlag(bool& f) : flag(f) { flag = true; }
    ~ScopedScrubReleaseFlag() { flag = false; }
};

// Static-display modes that map to AbstractStripPanel subclasses
// rendered via setOwnedPanels (NOT via the temporal-strip ChannelBuffer
// path). Each new panel kind landing in Phases C-G appends its mode
// here.
bool isPanelMode(const QString& mode) {
    return mode == QStringLiteral("static.bar.sequence")
        || mode == QStringLiteral("static.spectrum.power")
        || mode == QStringLiteral("static.curve.lag.animated")
        || mode == QStringLiteral("static.chord.coupling")
        || mode == QStringLiteral("static.fixed_freq");
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
    case model::SignalAxis::BondVector:
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
    case model::SignalAxis::BondVector:
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

        if (path == QStringLiteral("/trajectory/kernel_dynamics")) {
            // Five descriptors share this path. Curve-shaped ones
            // (ACF/PSD) bypass denseH5Plan and render via the owned-panels
            // path; only the 3 scalar reductions land here. Dispatch on
            // descriptor.conceptKey to pick the right (N, C) buffer, then
            // on channel.id to pick the column.
            const model::QtKernelDynamics* kd = h5->kernelDynamics();
            if (!kd)
                return gap(model::GapReason::SourceAbsent);
            const std::optional<std::size_t> row = rowFor(model::SignalAxis::Atom, kd->n_atoms);
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            const QString chan = channel.id;
            const int colIdx = kd->channel_names.indexOf(chan);
            if (colIdx < 0)
                return gap(model::GapReason::AnchorUnavailable);
            const std::size_t col = static_cast<std::size_t>(colIdx);
            const QString conceptKey = descriptor.conceptKey;
            if (conceptKey == QStringLiteral("kernel_dynamics.decay_time"))
                return finish(kd->decay_time_ps.at(*row, col));
            if (conceptKey == QStringLiteral("kernel_dynamics.peak_freq"))
                return finish(kd->peak_freq_per_ps.at(*row, col));
            if (conceptKey == QStringLiteral("kernel_dynamics.spectral_centroid"))
                return finish(kd->spectral_centroid_per_ps.at(*row, col));
            return gap(model::GapReason::Pending);
        }

        if (path == QStringLiteral("/trajectory/dihedral_autocorrelation")) {
            // Per-residue phi/psi corr_time scalars. ACF curves go via
            // the owned-panels path. Anchor is ResidueAnchor.
            const model::QtDihedralAutocorrelation* da = h5->dihedralAutocorrelation();
            if (!da)
                return gap(model::GapReason::SourceAbsent);
            const auto* res = std::get_if<model::ResidueAnchor>(&anchor);
            if (!res || res->residue >= da->n_residues)
                return gap(model::GapReason::AnchorUnavailable);
            const QString conceptKey = descriptor.conceptKey;
            if (conceptKey == QStringLiteral("dihedral.phi_corr_time"))
                return finish(da->phi_corr_time.at(res->residue));
            if (conceptKey == QStringLiteral("dihedral.psi_corr_time"))
                return finish(da->psi_corr_time.at(res->residue));
            return gap(model::GapReason::Pending);
        }

        if (path == QStringLiteral("/trajectory/reorientational_dynamics")) {
            // Scalar reductions only (the 2 TCF curves render via the
            // owned-panels path). Dispatch on conceptKey to pick which
            // QtPerBondVectorScalar; same BondVectorAnchor / Residue
            // widening as iRED.
            const model::QtReorientationalDynamics* rd = h5->reorientationalDynamics();
            if (!rd)
                return gap(model::GapReason::SourceAbsent);
            std::optional<std::size_t> row;
            if (const auto* vec = std::get_if<model::BondVectorAnchor>(&anchor)) {
                row = rd->identity.rowFor(vec->residue, vec->kind);
            } else if (const auto* res = std::get_if<model::ResidueAnchor>(&anchor)) {
                // Wildcard kind = first matching row (NH preferred by
                // identity-table order, since producer adds NH before
                // CaHa before CO per residue).
                row = rd->identity.rowFor(res->residue, /*kind=*/0);
            }
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            const QString conceptKey = descriptor.conceptKey;
            if (conceptKey == QStringLiteral("reorient.s2"))
                return finish(rd->s2.at(*row));
            if (conceptKey == QStringLiteral("reorient.tau_e"))
                return finish(rd->tau_e.at(*row));
            if (conceptKey == QStringLiteral("reorient.r1"))
                return finish(rd->r1.at(*row));
            if (conceptKey == QStringLiteral("reorient.r2"))
                return finish(rd->r2.at(*row));
            if (conceptKey == QStringLiteral("reorient.noe"))
                return finish(rd->noe.at(*row));
            return gap(model::GapReason::Pending);
        }

        if (path == QStringLiteral("/trajectory/ired_order_parameters")) {
            // Static per-bond-vector S² — same value every frame
            // (matches the bs_t0_autocorrelation pattern above). Anchor
            // is either a BondVectorAnchor(residue, kind) or — via the
            // BondVector ↔ Residue widening rule — a ResidueAnchor,
            // in which case we resolve to the first matching row
            // (kind=0 wildcard).
            const model::QtIRedOrderParameters* ts = h5->iredOrderParameters();
            if (!ts)
                return gap(model::GapReason::SourceAbsent);
            std::optional<std::size_t> row;
            if (const auto* vec = std::get_if<model::BondVectorAnchor>(&anchor)) {
                row = ts->identity.rowFor(vec->residue, vec->kind);
            } else if (const auto* res = std::get_if<model::ResidueAnchor>(&anchor)) {
                row = ts->identity.rowFor(res->residue, /*kind wildcard=*/0);
            }
            if (!row)
                return gap(model::GapReason::AnchorUnavailable);
            return finish(ts->s2_ired[*row]);
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
        // Codex NOW-1 (2026-05-29): displayRefsChanged must trigger a
        // full rebuild() rather than the lightweight
        // refreshPanelVisibility(). The dialog adds a static panel
        // signal in two steps — addSignal() (sync, fires the
        // controller's signalAdded handler before any refs exist) +
        // addDisplayRef() (fires displayRefsChanged after). The
        // earlier rebuild on signalAdded sees no refs and filters the
        // panel out; without a rebuild here the panel never appears
        // until some unrelated event triggers another rebuild.
        ACONNECT(panelModel_.data(), &model::DashboardPanelModel::displayRefsChanged,
                 this, [this](const QUuid&) { rebuild(); });
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

void DashboardDisplayController::setSceneOverlay(SceneRevealOverlay* overlay) {
    ASSERT_THREAD(this);
    // If the overlay is being detached, clear any tensor reveal we
    // pushed to the previous one so it doesn't linger after rebuild.
    if (sceneOverlay_ && sceneOverlay_ != overlay)
        sceneOverlay_->clearTensor();
    sceneOverlay_ = overlay;
    rebuild();
}

void DashboardDisplayController::setFrame(int frame) {
    ASSERT_THREAD(this);
    frame_ = std::max(0, frame);
    extendToFrame(frame_);
    emit stripTracksChanged();
}

// Slider drags defer per-frame snapshot fetches so one valueChanged
// cascade does not stack hundreds of synchronous NPY directory reads
// on the GUI thread. extendToFrame() bails while scrubActive_ is true;
// release runs one catch-up extendToFrame(frame_). See
// h5-reader/notes/ROBUSTNESS_BACKLOG_2026-05-30.md item #3.
void DashboardDisplayController::setScrubActive(bool active) {
    ASSERT_THREAD(this);
    if (scrubActive_ == active)
        return;
    scrubActive_ = active;
    if (!scrubActive_) {
        const ScopedScrubReleaseFlag guard(scrubReleasePending_);
        extendToFrame(frame_);
        emit stripTracksChanged();
    }
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
    std::vector<std::unique_ptr<AbstractStripPanel>> nextPanels;
    activeStripSignalCount_ = 0;

    // L-4 (2026-05-29): auto-compose. Pre-scan for Reorient scalar
    // signals with static.bar.sequence mode in the active panel.
    // If 2+ are present, they collapse into ONE composite panel
    // (built post-loop) instead of one panel per signal. Absorbed
    // signals get added to `absorbedSignals` so the per-signal
    // dispatch below skips their static.bar.sequence mode (other
    // modes on the same signal still emit normally).
    QSet<QUuid> absorbedSignals;
    QVector<model::DashboardSignal> reorientCompositeGroup;
    if (catalog_ && activeModel_) {
        QVector<model::DashboardSignal> candidates;
        const auto reorientScalarConcepts = QSet<QString>{
            QStringLiteral("reorient.s2"),
            QStringLiteral("reorient.tau_e"),
            QStringLiteral("reorient.r1"),
            QStringLiteral("reorient.r2"),
            QStringLiteral("reorient.noe"),
        };
        for (const model::DashboardSignal& signal : activeModel_->activeSignals()) {
            if (!signal.enabled) continue;
            const model::SignalDescriptor* d =
                catalog_->findDescriptor(signal.binding.descriptorId);
            if (!d) continue;
            if (d->storagePath != QStringLiteral("/trajectory/reorientational_dynamics"))
                continue;
            if (!reorientScalarConcepts.contains(d->conceptKey))
                continue;
            if (!signal.displayModeIds.contains(QStringLiteral("static.bar.sequence")))
                continue;
            if (panelModel_) {
                const model::DashboardPanel* activePanel = panelModel_->activePanel();
                if (!activePanel) continue;
                const model::DashboardDisplayRef ref{
                    signal.id,
                    QStringLiteral("static.bar.sequence"),
                    QStringLiteral("panel")};
                if (!activePanel->displays.contains(ref)) continue;
            }
            candidates.push_back(signal);
        }
        if (candidates.size() >= 2) {
            reorientCompositeGroup = candidates;
            for (const auto& s : candidates)
                absorbedSignals.insert(s.id);
        }
    }

    if (catalog_ && activeModel_) {
        for (const model::DashboardSignal& signal : activeModel_->activeSignals()) {
            if (!signal.enabled)
                continue;
            const model::SignalDescriptor* descriptor =
                catalog_->findDescriptor(signal.binding.descriptorId);
            if (!descriptor)
                continue;

            // Static-display path: build an AbstractStripPanel directly.
            // Mode + descriptor.storagePath dispatch — one branch per
            // (mode, source) pair landing in Phases C-G. We loop over
            // ALL panel-mode entries in displayModeIds (matches the
            // temporal-strip loop just below) so a signal carrying
            // both static.bar.sequence and static.table emits the
            // SequenceBarPanel regardless of which mode happens to be
            // the binding's primary.
            for (const QString& mode : signal.displayModeIds) {
                if (!isPanelMode(mode))
                    continue;
                // Active-panel filter: same scope rule the temporal
                // strip path uses (see seriesIsVisibleInActivePanel).
                // If the user has multiple dashboard tabs, the panel
                // emits only for the tab whose displays include this
                // (signal, mode, channel="panel") ref.
                if (panelModel_) {
                    const model::DashboardPanel* activePanel = panelModel_->activePanel();
                    if (!activePanel) continue;
                    const model::DashboardDisplayRef ref{
                        signal.id, mode, QStringLiteral("panel")};
                    if (!activePanel->displays.contains(ref))
                        continue;
                }
                const QString& path = descriptor->storagePath;
                if (path == QStringLiteral("/trajectory/ired_order_parameters")
                    && mode == QStringLiteral("static.bar.sequence")) {
                    if (auto panel = buildIRedSequenceBarPanel(signal, *descriptor))
                        nextPanels.push_back(std::move(panel));
                } else if (path == QStringLiteral("/trajectory/kernel_dynamics")
                           && mode == QStringLiteral("static.spectrum.power")) {
                    if (auto panel = buildKernelDynamicsPowerSpectrumPanel(signal, *descriptor))
                        nextPanels.push_back(std::move(panel));
                } else if (path == QStringLiteral("/trajectory/kernel_dynamics")
                           && mode == QStringLiteral("static.curve.lag.animated")) {
                    if (auto panel = buildKernelDynamicsLagDecayPanel(signal, *descriptor))
                        nextPanels.push_back(std::move(panel));
                } else if (path == QStringLiteral("/trajectory/reorientational_dynamics")
                           && mode == QStringLiteral("static.bar.sequence")) {
                    // L-4: skip if this signal is part of an
                    // auto-composed Reorient group (the composite
                    // panel is built post-loop below). Other modes
                    // on the same signal still emit normally.
                    if (!absorbedSignals.contains(signal.id)) {
                        if (auto panel = buildReorientSequenceBarPanel(signal, *descriptor))
                            nextPanels.push_back(std::move(panel));
                    }
                } else if (path == QStringLiteral("/trajectory/reorientational_dynamics")
                           && mode == QStringLiteral("static.curve.lag.animated")) {
                    if (auto panel = buildReorientLagDecayPanel(signal, *descriptor))
                        nextPanels.push_back(std::move(panel));
                } else if (path == QStringLiteral("/trajectory/dihedral_autocorrelation")
                           && mode == QStringLiteral("static.bar.sequence")) {
                    // Covers dihedral.phi/psi_corr_time AND
                    // dihedral.chi_corr_time (L-2a) via conceptKey
                    // dispatch inside the builder.
                    if (auto panel = buildDihedralSequenceBarPanel(signal, *descriptor))
                        nextPanels.push_back(std::move(panel));
                } else if (path == QStringLiteral("/trajectory/dihedral_autocorrelation")
                           && mode == QStringLiteral("static.curve.lag.animated")) {
                    // Covers dihedral.phi/psi_acf AND dihedral.chi_acf
                    // (L-2a) — the 4-channel chi variant fans into a
                    // single multi-curve LagDecayPanel.
                    if (auto panel = buildDihedralLagDecayPanel(signal, *descriptor))
                        nextPanels.push_back(std::move(panel));
                } else if (path == QStringLiteral("/trajectory/kernel_coherence")
                           && mode == QStringLiteral("static.chord.coupling")) {
                    if (auto panel = buildKernelCoherenceChordPanel(signal, *descriptor))
                        nextPanels.push_back(std::move(panel));
                } else if (path == QStringLiteral("/trajectory/reorientational_dynamics")
                           && mode == QStringLiteral("static.fixed_freq")) {
                    // L-3b (2026-05-29): J(ω) at 5 KTB Larmor combinations.
                    if (auto panel = buildReorientFixedFreqPanel(signal, *descriptor))
                        nextPanels.push_back(std::move(panel));
                }
                // Future per-phase branches land here.
            }

            // L-3a tensor-glyph trigger INTENTIONALLY OMITTED here
            // (decision 2026-05-29 mid-session, per user request).
            // The earlier draft auto-fired sceneOverlay_->revealTensor
            // for the first active Reorient orientation_tensor signal
            // in the active panel. That was the wrong scope: the user
            // explicitly chose "no UI yet, defer trigger to follow-up"
            // because auto-fire on signal addition is too eager — the
            // ellipsoid should appear only on an explicit user
            // gesture (a Reveal button, context menu, or panel
            // interaction) that does not yet exist. The follow-up
            // session designs and lands that gesture.
            //
            // The supporting infrastructure stays in place:
            //   - revealTensor / clearTensor on SceneRevealOverlay
            //   - setSceneOverlay on the controller + dock wiring
            //   - TensorGlyphMath + h5:reorient_orientation_tensor
            //     catalog descriptor
            // When the follow-up adds the gesture, the trigger code
            // lives here and calls sceneOverlay_->revealTensor(...).

            // Temporal-strip path: existing ChannelBuffer pipeline.
            if (hasStripMode(signal.displayModeIds)) {
                ++activeStripSignalCount_;
                buildGenericTracks(signal, *descriptor, next);
            }
        }
    }

    for (int i = 0; i < next.size(); ++i)
        next[i].color = colorForIndex(i);

    // L-4 (2026-05-29): build the auto-compose Reorient composite
    // AFTER the per-signal loop has built all the other panels.
    // The composite folds the 2+ absorbed scalar signals into one
    // SequenceBarPanel with overlays.
    if (!reorientCompositeGroup.isEmpty()) {
        if (auto panel = buildReorientCompositeBarPanel(reorientCompositeGroup))
            nextPanels.push_back(std::move(panel));
    }

    activeOwnedPanelCount_ = static_cast<int>(nextPanels.size());
    series_ = std::move(next);
    ownedPanels_ = std::move(nextPanels);
    extendToFrame(frame_);

    // L-3a: tensor-glyph trigger intentionally not wired here yet —
    // see the comment block above in the per-signal loop. The
    // sceneOverlay_ pointer + revealTensor / clearTensor API exist
    // and are reachable; the call site is the deferred follow-up.

    updateStatusText();
    emit stripTracksChanged();
    emit ownedPanelsChanged();
}

std::vector<std::unique_ptr<AbstractStripPanel>>
DashboardDisplayController::takeOwnedPanels() {
    return std::move(ownedPanels_);
}

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildIRedSequenceBarPanel(
        const model::DashboardSignal& signal,
        const model::SignalDescriptor& descriptor) const {
    if (!conformation_)
        return nullptr;
    const model::TrajectoryConformation* trajectory = conformation_->asTrajectory();
    if (!trajectory)
        return nullptr;
    const auto* h5 = trajectory->h5();
    if (!h5)
        return nullptr;
    const model::QtIRedOrderParameters* ired = h5->iredOrderParameters();
    if (!ired || ired->identity.n_vectors == 0)
        return nullptr;

    std::vector<SequenceBarRow> rows;
    rows.reserve(ired->identity.n_vectors);
    for (std::size_t i = 0; i < ired->identity.n_vectors; ++i) {
        SequenceBarRow row;
        row.residue_index = ired->identity.residue_index[i];
        row.value = (i < ired->s2_ired.size()) ? ired->s2_ired[i] : 0.0;
        row.kind = ired->identity.kind[i];
        rows.push_back(row);
    }

    // Bar-click → reveal a BondVectorAnchor for the bar's (residue, kind).
    // Captured by value so the panel is self-contained.
    SequenceBarPanel::BindingForRow bindingForRow =
        [rows, descriptorId = descriptor.id](std::size_t row) -> model::SignalBinding {
        model::SignalBinding b;
        b.descriptorId = descriptorId;
        if (row < rows.size()) {
            model::BondVectorAnchor anchor;
            anchor.residue = static_cast<std::size_t>(rows[row].residue_index);
            anchor.kind = rows[row].kind;
            b.anchor = anchor;
        }
        return b;
    };

    return std::make_unique<SequenceBarPanel>(
        signal.label.isEmpty() ? descriptor.label : signal.label,
        QStringLiteral("S²"),
        std::move(rows),
        std::move(bindingForRow),
        QColor(115, 229, 214),
        0.0,
        1.0);
}

// Helper: resolve the signal's anchor → atom row (atom axis only,
// no need for the BondVector machinery here).
static std::optional<std::size_t> atomRowForKernelSignal(
        const model::DashboardSignal& signal,
        std::size_t n_atoms) {
    if (const auto* a = std::get_if<model::AtomAnchor>(&signal.binding.anchor)) {
        if (a->atom < n_atoms) return a->atom;
    }
    return std::nullopt;
}

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildKernelDynamicsPowerSpectrumPanel(
        const model::DashboardSignal& signal,
        const model::SignalDescriptor& descriptor) const {
    if (!conformation_) return nullptr;
    const auto* trajectory = conformation_->asTrajectory();
    const auto* h5 = trajectory ? trajectory->h5() : nullptr;
    if (!h5) return nullptr;
    const auto* kd = h5->kernelDynamics();
    if (!kd || kd->n_atoms == 0) return nullptr;
    const auto row = atomRowForKernelSignal(signal, kd->n_atoms);
    if (!row) return nullptr;

    model::SignalBinding reveal;
    reveal.descriptorId = descriptor.id;
    reveal.anchor = model::AtomAnchor{*row};
    return std::make_unique<PowerSpectrumPanel>(
        signal.label.isEmpty() ? descriptor.label : signal.label,
        &kd->power_spectrum,
        *row,
        std::move(reveal));
}

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildKernelDynamicsLagDecayPanel(
        const model::DashboardSignal& signal,
        const model::SignalDescriptor& descriptor) const {
    if (!conformation_) return nullptr;
    const auto* trajectory = conformation_->asTrajectory();
    const auto* h5 = trajectory ? trajectory->h5() : nullptr;
    if (!h5) return nullptr;
    const auto* kd = h5->kernelDynamics();
    if (!kd || kd->n_atoms == 0) return nullptr;
    const auto row = atomRowForKernelSignal(signal, kd->n_atoms);
    if (!row) return nullptr;

    model::SignalBinding reveal;
    reveal.descriptorId = descriptor.id;
    reveal.anchor = model::AtomAnchor{*row};
    return std::make_unique<LagDecayPanel>(
        signal.label.isEmpty() ? descriptor.label : signal.label,
        &kd->acf,
        *row,
        std::move(reveal));
}

// Resolve a Reorient signal's anchor (BondVectorAnchor or widened
// ResidueAnchor) → identity-table row.
static std::optional<std::size_t> reorientRowFor(
        const model::DashboardSignal& signal,
        const model::QtReorientationalDynamics& rd) {
    if (const auto* vec = std::get_if<model::BondVectorAnchor>(&signal.binding.anchor))
        return rd.identity.rowFor(vec->residue, vec->kind);
    if (const auto* res = std::get_if<model::ResidueAnchor>(&signal.binding.anchor))
        return rd.identity.rowFor(res->residue, /*kind=*/0);
    return std::nullopt;
}

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildReorientSequenceBarPanel(
        const model::DashboardSignal& signal,
        const model::SignalDescriptor& descriptor) const {
    if (!conformation_) return nullptr;
    const auto* trajectory = conformation_->asTrajectory();
    const auto* h5 = trajectory ? trajectory->h5() : nullptr;
    if (!h5) return nullptr;
    const auto* rd = h5->reorientationalDynamics();
    if (!rd || rd->identity.n_vectors == 0) return nullptr;

    // Which scalar to plot? conceptKey selects.
    const QString conceptKey = descriptor.conceptKey;
    const model::QtPerBondVectorScalar* scalar = nullptr;
    QString unit = QStringLiteral("dimensionless");
    std::optional<double> yMin, yMax;
    if (conceptKey == QStringLiteral("reorient.s2"))   { scalar = &rd->s2;    yMin = 0.0; yMax = 1.0; }
    else if (conceptKey == QStringLiteral("reorient.tau_e")) { scalar = &rd->tau_e; unit = QStringLiteral("ps"); }
    else if (conceptKey == QStringLiteral("reorient.r1"))    { scalar = &rd->r1;    unit = QStringLiteral("s⁻¹"); }
    else if (conceptKey == QStringLiteral("reorient.r2"))    { scalar = &rd->r2;    unit = QStringLiteral("s⁻¹"); }
    else if (conceptKey == QStringLiteral("reorient.noe"))   { scalar = &rd->noe;   yMin = -1.0; yMax = 1.0; }
    if (!scalar) return nullptr;

    std::vector<SequenceBarRow> rows;
    rows.reserve(rd->identity.n_vectors);
    for (std::size_t i = 0; i < rd->identity.n_vectors; ++i) {
        SequenceBarRow row;
        row.residue_index = rd->identity.residue_index[i];
        row.value = scalar->at(i);
        row.kind = rd->identity.kind[i];
        if (!std::isfinite(row.value))
            continue;  // NaN rows (e.g. R1/R2/NOE on non-NH) drop out cleanly
        rows.push_back(row);
    }
    if (rows.empty()) return nullptr;

    SequenceBarPanel::BindingForRow bindingForRow =
        [rows, descriptorId = descriptor.id](std::size_t row) -> model::SignalBinding {
        model::SignalBinding b;
        b.descriptorId = descriptorId;
        if (row < rows.size()) {
            model::BondVectorAnchor anchor;
            anchor.residue = static_cast<std::size_t>(rows[row].residue_index);
            anchor.kind = rows[row].kind;
            b.anchor = anchor;
        }
        return b;
    };

    return std::make_unique<SequenceBarPanel>(
        signal.label.isEmpty() ? descriptor.label : signal.label,
        unit,
        std::move(rows),
        std::move(bindingForRow),
        QColor(255, 175, 76),  // amber — distinguishes Reorient from iRED's teal
        yMin,
        yMax);
}

namespace {

// Same conceptKey → (scalar buffer pointer, unit, y-bounds) lookup
// used by buildReorientSequenceBarPanel, factored out so the
// composite builder can call it per signal without duplicating the
// switch. Returns std::nullopt for unknown conceptKeys.
struct ReorientScalarPick {
    const model::QtPerBondVectorScalar* scalar = nullptr;
    QString unit;
    std::optional<double> yMin;
    std::optional<double> yMax;
};

std::optional<ReorientScalarPick> reorientScalarFor(
        const QString& conceptKey,
        const model::QtReorientationalDynamics& rd) {
    ReorientScalarPick p;
    if (conceptKey == QStringLiteral("reorient.s2"))   { p.scalar = &rd.s2;    p.unit = QStringLiteral("dimensionless"); p.yMin = 0.0; p.yMax = 1.0; }
    else if (conceptKey == QStringLiteral("reorient.tau_e")) { p.scalar = &rd.tau_e; p.unit = QStringLiteral("ps"); }
    else if (conceptKey == QStringLiteral("reorient.r1"))    { p.scalar = &rd.r1;    p.unit = QStringLiteral("s⁻¹"); }
    else if (conceptKey == QStringLiteral("reorient.r2"))    { p.scalar = &rd.r2;    p.unit = QStringLiteral("s⁻¹"); }
    else if (conceptKey == QStringLiteral("reorient.noe"))   { p.scalar = &rd.noe;   p.unit = QStringLiteral("dimensionless"); p.yMin = -1.0; p.yMax = 1.0; }
    if (!p.scalar) return std::nullopt;
    return p;
}

std::vector<SequenceBarRow> reorientScalarRows(
        const model::QtReorientationalDynamics& rd,
        const model::QtPerBondVectorScalar& scalar) {
    std::vector<SequenceBarRow> rows;
    rows.reserve(rd.identity.n_vectors);
    for (std::size_t i = 0; i < rd.identity.n_vectors; ++i) {
        SequenceBarRow row;
        row.residue_index = rd.identity.residue_index[i];
        row.value = scalar.at(i);
        row.kind = rd.identity.kind[i];
        if (!std::isfinite(row.value)) continue;
        rows.push_back(row);
    }
    return rows;
}

QColor compositeOverlayColor(std::size_t index) {
    // Distinct overlay colours, chosen for legibility against the
    // amber primary used by buildReorientSequenceBarPanel.
    static const QColor palette[] = {
        QColor( 74, 184, 220),  // cyan
        QColor(143, 100, 220),  // violet
        QColor(245,  90, 130),  // pink
        QColor(110, 200, 110),  // green
    };
    return palette[index % (sizeof(palette) / sizeof(palette[0]))];
}

}  // namespace

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildReorientCompositeBarPanel(
        const QVector<model::DashboardSignal>& group) const {
    if (group.size() < 2 || !conformation_ || !catalog_) return nullptr;
    const auto* trajectory = conformation_->asTrajectory();
    const auto* h5 = trajectory ? trajectory->h5() : nullptr;
    if (!h5) return nullptr;
    const auto* rd = h5->reorientationalDynamics();
    if (!rd || rd->identity.n_vectors == 0) return nullptr;

    // Primary: first signal in the group. Resolve its
    // descriptor → scalar + units + y-bounds.
    const model::DashboardSignal& primarySignal = group[0];
    const model::SignalDescriptor* primaryDescriptor =
        catalog_->findDescriptor(primarySignal.binding.descriptorId);
    if (!primaryDescriptor) return nullptr;
    const auto primaryPick = reorientScalarFor(primaryDescriptor->conceptKey, *rd);
    if (!primaryPick) return nullptr;
    auto primaryRows = reorientScalarRows(*rd, *primaryPick->scalar);
    if (primaryRows.empty()) return nullptr;

    SequenceBarPanel::BindingForRow bindingForRow =
        [primaryRows, descriptorId = primaryDescriptor->id](std::size_t row) -> model::SignalBinding {
        model::SignalBinding b;
        b.descriptorId = descriptorId;
        if (row < primaryRows.size()) {
            model::BondVectorAnchor anchor;
            anchor.residue = static_cast<std::size_t>(primaryRows[row].residue_index);
            anchor.kind = primaryRows[row].kind;
            b.anchor = anchor;
        }
        return b;
    };

    // Composite label: list all conceptKeys in the group, primary
    // first. Helps the user understand which signals collapsed
    // into this panel.
    QStringList allLabels;
    allLabels.push_back(primaryDescriptor->label);
    for (int i = 1; i < group.size(); ++i) {
        if (auto* d = catalog_->findDescriptor(group[i].binding.descriptorId))
            allLabels.push_back(d->label);
    }
    const QString panelLabel = QStringLiteral("Reorient composite: %1")
                                   .arg(allLabels.join(QStringLiteral(" + ")));

    auto panel = std::make_unique<SequenceBarPanel>(
        panelLabel,
        primaryPick->unit,
        std::move(primaryRows),
        std::move(bindingForRow),
        QColor(255, 175, 76),  // amber — same as single-signal Reorient
        primaryPick->yMin,
        primaryPick->yMax);

    // Add overlays for the remaining group members.
    for (int i = 1; i < group.size(); ++i) {
        const model::DashboardSignal& s = group[i];
        const model::SignalDescriptor* d =
            catalog_->findDescriptor(s.binding.descriptorId);
        if (!d) continue;
        const auto pick = reorientScalarFor(d->conceptKey, *rd);
        if (!pick) continue;
        auto rows = reorientScalarRows(*rd, *pick->scalar);
        if (rows.empty()) continue;
        panel->addOverlay(std::move(rows),
                          compositeOverlayColor(static_cast<std::size_t>(i - 1)),
                          d->label,
                          pick->unit);
    }
    return panel;
}

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildReorientLagDecayPanel(
        const model::DashboardSignal& signal,
        const model::SignalDescriptor& descriptor) const {
    if (!conformation_) return nullptr;
    const auto* trajectory = conformation_->asTrajectory();
    const auto* h5 = trajectory ? trajectory->h5() : nullptr;
    if (!h5) return nullptr;
    const auto* rd = h5->reorientationalDynamics();
    if (!rd || rd->identity.n_vectors == 0) return nullptr;
    const auto row = reorientRowFor(signal, *rd);
    if (!row) return nullptr;

    // Pick body or lab TCF per descriptor.
    const model::QtPerBondVectorCurve* tcf = nullptr;
    if (descriptor.conceptKey == QStringLiteral("reorient.acf_internal"))
        tcf = &rd->acf_internal;
    else if (descriptor.conceptKey == QStringLiteral("reorient.acf_lab"))
        tcf = &rd->acf_lab;
    if (!tcf || tcf->n_samples < 2) return nullptr;

    // Synthesise a single-row, single-channel view from the per-bond-vector
    // curve so LagDecayPanel can render it (one polyline). The panel owns
    // this view; no lifetime coupling to the H5 buffer.
    auto view = std::make_unique<model::QtPerAtomChannelCurve>();
    view->n_atoms = 1;
    view->n_channels = 1;
    view->n_samples = tcf->n_samples;
    view->data.assign(tcf->n_samples, 0.0);
    for (std::size_t s = 0; s < tcf->n_samples; ++s)
        view->data[s] = tcf->at(*row, s);
    view->axis_values = tcf->axis_values;
    view->axis_unit = tcf->axis_unit;
    view->axis_label = QStringLiteral("lag");
    view->units = tcf->units;
    view->channel_names = QStringList{descriptor.conceptKey};

    model::SignalBinding reveal;
    reveal.descriptorId = descriptor.id;
    model::BondVectorAnchor anchor;
    anchor.residue = static_cast<std::size_t>(rd->identity.residue_index[*row]);
    anchor.kind = rd->identity.kind[*row];
    reveal.anchor = anchor;

    return std::make_unique<LagDecayPanel>(
        signal.label.isEmpty() ? descriptor.label : signal.label,
        std::move(view),
        /*atomRow=*/0,
        std::move(reveal));
}

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildKernelCoherenceChordPanel(
        const model::DashboardSignal& signal,
        const model::SignalDescriptor& descriptor) const {
    if (!conformation_) return nullptr;
    const auto* trajectory = conformation_->asTrajectory();
    const auto* h5 = trajectory ? trajectory->h5() : nullptr;
    if (!h5) return nullptr;
    const auto* kc = h5->kernelCoherence();
    if (!kc || kc->matrix.n_atoms == 0) return nullptr;
    const auto row = atomRowForKernelSignal(signal, kc->matrix.n_atoms);
    if (!row) return nullptr;

    model::SignalBinding reveal;
    reveal.descriptorId = descriptor.id;
    reveal.anchor = model::AtomAnchor{*row};
    return std::make_unique<ChordCouplingPanel>(
        signal.label.isEmpty() ? descriptor.label : signal.label,
        &kc->matrix,
        *row,
        /*threshold=*/0.3,
        std::move(reveal));
}

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildDihedralSequenceBarPanel(
        const model::DashboardSignal& signal,
        const model::SignalDescriptor& descriptor) const {
    if (!conformation_) return nullptr;
    const auto* trajectory = conformation_->asTrajectory();
    const auto* h5 = trajectory ? trajectory->h5() : nullptr;
    if (!h5) return nullptr;
    const auto* da = h5->dihedralAutocorrelation();
    if (!da || da->n_residues == 0) return nullptr;

    const QString conceptKey = descriptor.conceptKey;
    const model::QtPerResidueScalar* scalar = nullptr;
    if (conceptKey == QStringLiteral("dihedral.phi_corr_time"))      scalar = &da->phi_corr_time;
    else if (conceptKey == QStringLiteral("dihedral.psi_corr_time")) scalar = &da->psi_corr_time;

    std::vector<SequenceBarRow> rows;
    if (scalar) {
        rows.reserve(scalar->n_residues);
        for (std::size_t i = 0; i < scalar->n_residues; ++i) {
            if (!scalar->isDefined(i)) continue;
            SequenceBarRow row;
            row.residue_index = static_cast<std::int32_t>(i);
            row.value = scalar->at(i);
            if (!std::isfinite(row.value)) continue;
            rows.push_back(row);
        }
    } else if (conceptKey == QStringLiteral("dihedral.chi_corr_time")) {
        // L-2a (2026-05-29): chi composite fanned across 4 sub-slots
        // (kinds 4..7) per residue. Residues with fewer than 4 chi
        // torsions just emit fewer rows (chi_defined gates emission).
        rows.reserve(da->n_residues * 4);
        for (std::size_t r = 0; r < da->n_residues; ++r) {
            for (std::size_t k = 0; k < 4; ++k) {
                if (!da->chiIsDefined(r, k)) continue;
                const double v = da->chiCorrTimeAt(r, k);
                if (!std::isfinite(v)) continue;
                SequenceBarRow row;
                row.residue_index = static_cast<std::int32_t>(r);
                row.value = v;
                row.kind = static_cast<std::uint8_t>(4 + k);  // chi0=4, chi3=7
                rows.push_back(row);
            }
        }
    }
    if (rows.empty()) return nullptr;

    SequenceBarPanel::BindingForRow bindingForRow =
        [rows, descriptorId = descriptor.id](std::size_t row) -> model::SignalBinding {
        model::SignalBinding b;
        b.descriptorId = descriptorId;
        if (row < rows.size())
            b.anchor = model::ResidueAnchor{static_cast<std::size_t>(rows[row].residue_index)};
        return b;
    };

    return std::make_unique<SequenceBarPanel>(
        signal.label.isEmpty() ? descriptor.label : signal.label,
        QStringLiteral("ps"),
        std::move(rows),
        std::move(bindingForRow),
        QColor(135, 211, 124),  // green
        std::nullopt,
        std::nullopt);
}

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildReorientFixedFreqPanel(
        const model::DashboardSignal& signal,
        const model::SignalDescriptor& descriptor) const {
    if (!conformation_) return nullptr;
    const auto* trajectory = conformation_->asTrajectory();
    const auto* h5 = trajectory ? trajectory->h5() : nullptr;
    if (!h5) return nullptr;
    const auto* rd = h5->reorientationalDynamics();
    if (!rd) return nullptr;
    const auto& src = rd->spectral_density_J;
    if (src.n_vectors == 0 || src.n_freqs == 0) return nullptr;

    // Resolve anchor → row (direct BondVectorAnchor or Residue widening).
    std::optional<std::size_t> row;
    if (const auto* vec = std::get_if<model::BondVectorAnchor>(&signal.binding.anchor))
        row = rd->identity.rowFor(vec->residue, vec->kind);
    else if (const auto* res = std::get_if<model::ResidueAnchor>(&signal.binding.anchor))
        row = rd->identity.rowFor(res->residue, /*wildcard=*/0);
    if (!row || *row >= src.n_vectors) return nullptr;

    // Slice the producer-owned (V, n_freqs) buffer down to a one-row
    // copy the panel owns. Keeps lifetime local to the panel — same
    // ownership model as the LagDecayPanel single-row clone path.
    auto view = std::make_unique<model::QtPerBondVectorFixedFreqBlock>();
    view->n_vectors = 1;
    view->n_freqs   = src.n_freqs;
    view->data.assign(src.n_freqs, 0.0);
    for (std::size_t i = 0; i < src.n_freqs; ++i)
        view->data[i] = src.at(*row, i);
    view->freq_values = src.freq_values;
    view->units = src.units;
    view->result_name = src.result_name;

    model::SignalBinding reveal;
    reveal.descriptorId = descriptor.id;
    reveal.anchor = signal.binding.anchor;

    return std::make_unique<FixedFreqPanel>(
        signal.label.isEmpty() ? descriptor.label : signal.label,
        std::move(view),
        /*row=*/0,
        std::move(reveal));
}

std::unique_ptr<AbstractStripPanel>
DashboardDisplayController::buildDihedralLagDecayPanel(
        const model::DashboardSignal& signal,
        const model::SignalDescriptor& descriptor) const {
    if (!conformation_) return nullptr;
    const auto* trajectory = conformation_->asTrajectory();
    const auto* h5 = trajectory ? trajectory->h5() : nullptr;
    if (!h5) return nullptr;
    const auto* da = h5->dihedralAutocorrelation();
    if (!da || da->n_residues == 0) return nullptr;

    const auto* res = std::get_if<model::ResidueAnchor>(&signal.binding.anchor);
    if (!res || res->residue >= da->n_residues) return nullptr;

    const model::QtPerResidueCurve* curve = nullptr;
    if (descriptor.conceptKey == QStringLiteral("dihedral.phi_acf"))      curve = &da->phi_acf;
    else if (descriptor.conceptKey == QStringLiteral("dihedral.psi_acf")) curve = &da->psi_acf;

    auto view = std::make_unique<model::QtPerAtomChannelCurve>();
    view->n_atoms = 1;
    view->axis_label = QStringLiteral("lag");

    if (curve && curve->n_samples >= 2) {
        view->n_channels = 1;
        view->n_samples = curve->n_samples;
        view->data.assign(curve->n_samples, 0.0);
        for (std::size_t s = 0; s < curve->n_samples; ++s)
            view->data[s] = curve->at(res->residue, s);
        view->axis_values = curve->axis_values;
        view->axis_unit = curve->axis_unit;
        view->units = curve->units;
        view->channel_names = QStringList{descriptor.conceptKey};
    } else if (descriptor.conceptKey == QStringLiteral("dihedral.chi_acf")
               && da->n_lags >= 2 && !da->chi_acf.empty()) {
        // L-2a (2026-05-29): pack 4 chi curves as 4 channels — LagDecayPanel
        // already iterates n_channels and draws one polyline per channel.
        view->n_channels = 4;
        view->n_samples = da->n_lags;
        view->data.assign(4 * da->n_lags, 0.0);
        for (std::size_t k = 0; k < 4; ++k) {
            for (std::size_t s = 0; s < da->n_lags; ++s)
                view->data[k * da->n_lags + s] = da->chiAcfAt(res->residue, k, s);
        }
        view->axis_values = da->chi_acf_axis;
        view->axis_unit = QStringLiteral("ps");
        view->units = QStringLiteral("dimensionless");
        view->channel_names = QStringList{QStringLiteral("chi0"),
                                          QStringLiteral("chi1"),
                                          QStringLiteral("chi2"),
                                          QStringLiteral("chi3")};
    } else {
        return nullptr;
    }

    model::SignalBinding reveal;
    reveal.descriptorId = descriptor.id;
    reveal.anchor = model::ResidueAnchor{res->residue};

    return std::make_unique<LagDecayPanel>(
        signal.label.isEmpty() ? descriptor.label : signal.label,
        std::move(view),
        /*atomRow=*/0,
        std::move(reveal));
}

void DashboardDisplayController::refreshPanelVisibility() {
    ASSERT_THREAD(this);
    updateStatusText();
    emit stripTracksChanged();
}

void DashboardDisplayController::updateStatusText() {
    if (!catalog_ || !activeModel_) {
        statusText_ = QStringLiteral("Dashboard signal model is not connected.");
        return;
    }
    if (activeStripSignalCount_ == 0 && activeOwnedPanelCount_ == 0) {
        statusText_ = QStringLiteral("No active strip or panel display modes.");
        return;
    }
    const int displayedSeries = activePanelSeriesCount();
    QString text;
    if (activeStripSignalCount_ > 0) {
        text = QStringLiteral("%1 displayed signal time series from %2 active strip signal%3")
                   .arg(displayedSeries)
                   .arg(activeStripSignalCount_)
                   .arg(activeStripSignalCount_ == 1 ? QString() : QStringLiteral("s"));
    }
    if (activeOwnedPanelCount_ > 0) {
        const QString panelClause = QStringLiteral("%1 panel signal%2")
                                        .arg(activeOwnedPanelCount_)
                                        .arg(activeOwnedPanelCount_ == 1 ? QString() : QStringLiteral("s"));
        text = text.isEmpty()
                   ? panelClause + QStringLiteral(".")
                   : text + QStringLiteral(" + ") + panelClause + QStringLiteral(".");
    } else {
        text += QStringLiteral(".");
    }
    if (activeStripSignalCount_ > 0)
        text += QStringLiteral(" Unimplemented sources append explicit pending gaps.");
    statusText_ = text;
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
    case model::SignalAxis::BondVector:
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
    if (scrubActive_)
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

    long long firstFrameToSample = startFrame;
    if (scrubReleasePending_ && firstFrameToSample < frame) {
        for (long long f = firstFrameToSample; f < frame; ++f) {
            for (ActiveSeries& series : series_) {
                if (series.buffer.lastFrame() < f)
                    series.buffer.append(model::FrameSignalSample::Gap(model::GapReason::Pending));
            }
        }
        firstFrameToSample = frame;
    }

    for (long long f = firstFrameToSample; f <= frame; ++f) {
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

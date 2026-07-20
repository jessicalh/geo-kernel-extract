#include "DashboardSignalModel.h"

#include "../diagnostics/DashboardLogging.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "DisplayModeCapability.h"

#include <QVariantList>
#include <QVariantMap>
#include <QStringList>

#include <algorithm>
#include <utility>

namespace h5reader::model {

namespace {

QStringList normalizedModeList(const QStringList& modes) {
    QStringList result;
    result.reserve(modes.size());
    for (const QString& mode : modes) {
        const QString trimmed = mode.trimmed();
        if (!trimmed.isEmpty() && !result.contains(trimmed))
            result.push_back(trimmed);
    }
    return result;
}

SignalAxis displayAxis(const DashboardSignal& signal) {
    if (signal.requiredAnchor != SignalAxis::None)
        return signal.requiredAnchor;
    if (signal.nativeAxis != SignalAxis::None)
        return signal.nativeAxis;
    return AxisForAnchor(signal.binding.anchor);
}

QList<int> allRoles() {
    return {
        Qt::DisplayRole,
        Qt::EditRole,
        Qt::CheckStateRole,
        DashboardSignalModel::IdRole,
        DashboardSignalModel::UuidRole,
        DashboardSignalModel::LabelRole,
        DashboardSignalModel::SourceRole,
        DashboardSignalModel::SourceKindRole,
        DashboardSignalModel::AxisRole,
        DashboardSignalModel::AxisNameRole,
        DashboardSignalModel::ValueShapeRole,
        DashboardSignalModel::ValueShapeNameRole,
        DashboardSignalModel::DisplayModesRole,
        DashboardSignalModel::EnabledRole,
        DashboardSignalModel::DescriptorIdRole,
        DashboardSignalModel::ConceptKeyRole,
        DashboardSignalModel::ReducerIdRole,
        DashboardSignalModel::DisplayModeRole,
        DashboardSignalModel::FollowsFocusRole,
        DashboardSignalModel::AvailabilityStateRole,
        DashboardSignalModel::AvailabilityRole,
        DashboardSignalModel::AvailabilityReasonRole,
        DashboardSignalModel::RenderableModeCountRole,
        DashboardSignalModel::ModeRenderabilityRole,
    };
}

QString availabilityReasonText(const DashboardSignalModel::SignalAvailability& availability) {
    if (!availability.evaluated)
        return QStringLiteral("Availability has not been evaluated for this run.");

    switch (availability.state) {
    case TrajectoryFieldAvailabilityState::Absent:
        return QStringLiteral("Source field is absent from the loaded run.");
    case TrajectoryFieldAvailabilityState::NoFramePayload:
        return QStringLiteral("Frame-local payload is not loaded for this run.");
    case TrajectoryFieldAvailabilityState::AllMissing:
        return QStringLiteral("Source field has no finite samples in the loaded run.");
    case TrajectoryFieldAvailabilityState::AllZeroStructural:
        return QStringLiteral("Source field is structurally all zero in the loaded run.");
    case TrajectoryFieldAvailabilityState::AllZeroObserved:
        return QStringLiteral("Source field contains only observed zero samples.");
    case TrajectoryFieldAvailabilityState::Available:
        return QStringLiteral("Source field has finite non-zero samples.");
    }
    return QStringLiteral("Availability state is unknown.");
}

QString uuidLogValue(const QUuid& id) {
    return id.toString(QUuid::WithoutBraces);
}

QString modesLogValue(const QStringList& modes) {
    return QStringLiteral("[%1]").arg(modes.join(QStringLiteral(",")));
}

QString selectedIdsLogValue(const QVector<DashboardSignal>& selectedSignals) {
    QStringList ids;
    ids.reserve(selectedSignals.size());
    for (const DashboardSignal& signal : selectedSignals)
        ids.push_back(uuidLogValue(signal.id));
    return modesLogValue(ids);
}

QString anchorLogValue(const SignalAnchor& anchor) {
    auto indexValue = [](const char* kind, std::size_t value) {
        return QStringLiteral("%1:%2").arg(QString::fromLatin1(kind)).arg(static_cast<qulonglong>(value));
    };
    if (std::holds_alternative<NoneAnchor>(anchor))
        return QStringLiteral("none");
    if (const auto* a = std::get_if<AtomAnchor>(&anchor))
        return indexValue("atom", a->atom);
    if (const auto* r = std::get_if<ResidueAnchor>(&anchor))
        return indexValue("residue", r->residue);
    if (const auto* t = std::get_if<AtomTupleAnchor>(&anchor)) {
        QStringList atoms;
        atoms.reserve(static_cast<int>(t->atoms.size()));
        for (std::size_t atom : t->atoms)
            atoms.push_back(QString::number(static_cast<qulonglong>(atom)));
        return QStringLiteral("atom_tuple:%1").arg(atoms.join(QStringLiteral("+")));
    }
    if (const auto* b = std::get_if<BondAnchor>(&anchor))
        return indexValue("bond", b->bond);
    if (const auto* v = std::get_if<BondVectorAnchor>(&anchor)) {
        return QStringLiteral("bond_vector:%1:%2")
            .arg(static_cast<qulonglong>(v->residue))
            .arg(static_cast<unsigned>(v->kind));
    }
    if (const auto* r = std::get_if<RingAnchor>(&anchor))
        return indexValue("ring", r->ring);
    if (const auto* r = std::get_if<AromaticRingAnchor>(&anchor))
        return indexValue("aromatic_ring", r->ring);
    if (const auto* r = std::get_if<SaturatedRingAnchor>(&anchor))
        return indexValue("saturated_ring", r->ring);
    if (const auto* p = std::get_if<RingContributionPairAnchor>(&anchor))
        return indexValue("ring_contribution_pair", p->pair);
    if (const auto* m = std::get_if<RingMembershipAnchor>(&anchor))
        return indexValue("ring_membership", m->membership);
    if (const auto* p = std::get_if<MutationMatchPairAnchor>(&anchor))
        return indexValue("mutation_match_pair", p->pair);
    if (std::holds_alternative<ProteinAnchor>(anchor))
        return QStringLiteral("protein");
    if (std::holds_alternative<SystemAnchor>(anchor))
        return QStringLiteral("system");
    if (std::holds_alternative<EventAnchor>(anchor))
        return QStringLiteral("event");
    return QStringLiteral("unknown");
}

}  // namespace

DashboardSignalModel::DashboardSignalModel(QObject* parent)
    : QAbstractListModel(parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("DashboardSignalModel"));
}

bool DashboardSignalModel::SignalAvailability::isVisible() const {
    return !evaluated || TrajectoryFieldAvailability::isVisibleState(state);
}

bool DashboardSignalModel::SignalAvailability::operator==(const SignalAvailability& other) const {
    return evaluated == other.evaluated
           && state == other.state
           && reason == other.reason
           && finiteSamples == other.finiteSamples
           && nonZeroSamples == other.nonZeroSamples;
}

int DashboardSignalModel::rowCount(const QModelIndex& parent) const {
    if (parent.isValid())
        return 0;
    return static_cast<int>(signals_.size());
}

QVariant DashboardSignalModel::data(const QModelIndex& index, int role) const {
    if (!index.isValid() || index.column() != 0)
        return {};
    const int row = index.row();
    if (row < 0 || row >= signals_.size())
        return {};

    const DashboardSignal& signal = signals_.at(row);
    switch (role) {
    case Qt::DisplayRole:
    case Qt::EditRole:
    case LabelRole:
        return signal.label;
    case Qt::CheckStateRole:
        return signal.enabled ? Qt::Checked : Qt::Unchecked;
    case IdRole:
        return signal.id.toString(QUuid::WithoutBraces);
    case UuidRole:
        return signal.id;
    case SourceRole:
        return ToString(signal.binding.sourceKind);
    case SourceKindRole:
        return static_cast<int>(signal.binding.sourceKind);
    case AxisRole:
        return static_cast<int>(displayAxis(signal));
    case AxisNameRole:
        return ToString(displayAxis(signal));
    case ValueShapeRole:
        return static_cast<int>(signal.valueShape);
    case ValueShapeNameRole:
        return ToString(signal.valueShape);
    case DisplayModesRole:
        return signal.displayModeIds;
    case EnabledRole:
        return signal.enabled;
    case DescriptorIdRole:
        return signal.binding.descriptorId;
    case ConceptKeyRole:
        return signal.binding.conceptKey;
    case ReducerIdRole:
        return signal.binding.reducerId;
    case DisplayModeRole:
        return signal.binding.displayModeId;
    case FollowsFocusRole:
        return signal.binding.followsFocus;
    case AvailabilityStateRole:
        return static_cast<int>(availabilityState(row));
    case AvailabilityRole:
        return availabilityName(row);
    case AvailabilityReasonRole:
        return availabilityReason(row);
    case RenderableModeCountRole:
        return renderableModeCount(row);
    case ModeRenderabilityRole:
        return modeRenderabilityVariantList(row);
    default:
        return {};
    }
}

bool DashboardSignalModel::setData(const QModelIndex& index, const QVariant& value, int role) {
    if (!index.isValid() || index.column() != 0)
        return false;
    const int row = index.row();
    if (row < 0 || row >= signals_.size())
        return false;

    DashboardSignal& signal = signals_[row];
    switch (role) {
    case Qt::EditRole:
    case LabelRole: {
        const QString label = value.toString();
        if (signal.label == label)
            return false;
        signal.label = label;
        emitRowChanged(row, {Qt::DisplayRole, Qt::EditRole, LabelRole});
        return true;
    }
    case Qt::CheckStateRole:
    case EnabledRole: {
        const bool enabled = role == Qt::CheckStateRole ? value.toInt() == Qt::Checked : value.toBool();
        if (signal.enabled == enabled)
            return false;
        signal.enabled = enabled;
        emitRowChanged(row, {Qt::CheckStateRole, EnabledRole});
        return true;
    }
    case DisplayModesRole:
        return setDisplayModes(signal.id, value.toStringList());
    case DisplayModeRole: {
        const QString mode = value.toString().trimmed();
        if (signal.binding.displayModeId == mode)
            return false;
        signal.binding.displayModeId = mode;
        if (!mode.isEmpty() && !signal.displayModeIds.contains(mode))
            signal.displayModeIds.push_back(mode);
        emitRowChanged(row, {DisplayModesRole,
                             DisplayModeRole,
                             RenderableModeCountRole,
                             ModeRenderabilityRole});
        return true;
    }
    default:
        return false;
    }
}

bool DashboardSignalModel::removeRows(int row, int count, const QModelIndex& parent) {
    if (parent.isValid() || row < 0 || count <= 0 || row + count > signals_.size())
        return false;

    QVector<QUuid> removedIds;
    removedIds.reserve(count);
    for (int i = 0; i < count; ++i)
        removedIds.push_back(signals_.at(row + i).id);

    beginRemoveRows(QModelIndex(), row, row + count - 1);
    signals_.erase(signals_.begin() + row, signals_.begin() + row + count);
    availabilityByRow_.erase(availabilityByRow_.begin() + row,
                             availabilityByRow_.begin() + row + count);
    endRemoveRows();

    for (const QUuid& id : removedIds)
        emit signalRemoved(id);
    return true;
}

Qt::ItemFlags DashboardSignalModel::flags(const QModelIndex& index) const {
    if (!index.isValid())
        return Qt::NoItemFlags;
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable | Qt::ItemIsUserCheckable;
}

QHash<int, QByteArray> DashboardSignalModel::roleNames() const {
    QHash<int, QByteArray> names = QAbstractListModel::roleNames();
    names[IdRole] = "id";
    names[UuidRole] = "uuid";
    names[LabelRole] = "label";
    names[SourceRole] = "source";
    names[SourceKindRole] = "sourceKind";
    names[AxisRole] = "axis";
    names[AxisNameRole] = "axisName";
    names[ValueShapeRole] = "valueShape";
    names[ValueShapeNameRole] = "valueShapeName";
    names[DisplayModesRole] = "displayModes";
    names[EnabledRole] = "enabled";
    names[DescriptorIdRole] = "descriptorId";
    names[ConceptKeyRole] = "conceptKey";
    names[ReducerIdRole] = "reducerId";
    names[DisplayModeRole] = "displayMode";
    names[FollowsFocusRole] = "followsFocus";
    names[AvailabilityStateRole] = "availabilityState";
    names[AvailabilityRole] = "availability";
    names[AvailabilityReasonRole] = "availabilityReason";
    names[RenderableModeCountRole] = "renderableModeCount";
    names[ModeRenderabilityRole] = "modeRenderability";
    return names;
}

const DashboardSignal* DashboardSignalModel::signalById(const QUuid& id) const {
    const int row = rowForId(id);
    if (row < 0)
        return nullptr;
    return &signals_.at(row);
}

DashboardSignal* DashboardSignalModel::signalById(const QUuid& id) {
    const int row = rowForId(id);
    if (row < 0)
        return nullptr;
    return &signals_[row];
}

QModelIndex DashboardSignalModel::indexForId(const QUuid& id) const {
    const int row = rowForId(id);
    if (row < 0)
        return {};
    return index(row, 0);
}

int DashboardSignalModel::rowForId(const QUuid& id) const {
    for (int row = 0; row < signals_.size(); ++row) {
        if (signals_.at(row).id == id)
            return row;
    }
    return -1;
}

void DashboardSignalModel::setFieldAvailability(
    std::shared_ptr<const TrajectoryFieldAvailability> availability) {
    ASSERT_THREAD(this);
    availability_ = std::move(availability);
    refreshAvailability();
}

void DashboardSignalModel::refreshAvailability() {
    ASSERT_THREAD(this);
    if (availabilityByRow_.size() != signals_.size())
        availabilityByRow_.resize(signals_.size());

    for (int row = 0; row < signals_.size(); ++row) {
        const SignalAvailability next = resolveAvailability(signals_.at(row));
        if (availabilityByRow_.at(row) == next)
            continue;
        availabilityByRow_[row] = next;
        emitAvailabilityChanged(row);
    }
}

DashboardSignalModel::SignalAvailability DashboardSignalModel::availabilityAt(int row) const {
    if (row < 0 || row >= availabilityByRow_.size())
        return {};
    return availabilityByRow_.at(row);
}

TrajectoryFieldAvailabilityState DashboardSignalModel::availabilityState(int row) const {
    return availabilityAt(row).state;
}

QString DashboardSignalModel::availabilityName(int row) const {
    const SignalAvailability availability = availabilityAt(row);
    if (!availability.evaluated)
        return QStringLiteral("unknown");
    return QString::fromLatin1(ToString(availability.state));
}

QString DashboardSignalModel::availabilityReason(int row) const {
    return availabilityAt(row).reason;
}

bool DashboardSignalModel::isVisibleAvailable(int row) const {
    return availabilityAt(row).isVisible();
}

DashboardSignalModel::ModeRenderability DashboardSignalModel::ModeRenderabilityFor(const QString& mode) {
    const DisplayModeCapability capability = DisplayModeCapabilityFor(mode);
    return ModeRenderability{
        mode,
        capability.hasVisibleSurface,
        capability.buildsPanelWidget,
        capability.emitsPanelRef,
    };
}

DashboardSignalModel::ModeRenderability
DashboardSignalModel::ModeRenderabilityFor(const DashboardSignal& signal,
                                           const QString& mode) {
    ModeRenderability renderability = ModeRenderabilityFor(mode);
    if (signal.binding.sourceKind == SignalSourceKind::ExperimentalShieldingMl
        && signal.binding.descriptorId
               == QStringLiteral("ml:experimental_shielding_t2")
        && mode == QStringLiteral("static.tensor")) {
        renderability.hasVisibleSurface = true;
    }
    return renderability;
}

QVector<DashboardSignalModel::ModeRenderability> DashboardSignalModel::modeRenderability(int row) const {
    QVector<ModeRenderability> out;
    if (row < 0 || row >= signals_.size())
        return out;
    const DashboardSignal& signal = signals_.at(row);
    const QStringList modes = signal.displayModeIds;
    out.reserve(modes.size());
    for (const QString& mode : modes)
        out.push_back(ModeRenderabilityFor(signal, mode));
    return out;
}

int DashboardSignalModel::renderableModeCount(int row) const {
    int count = 0;
    for (const ModeRenderability& mode : modeRenderability(row)) {
        if (mode.isRenderable())
            ++count;
    }
    return count;
}

int DashboardSignalModel::selectedCount() const {
    return static_cast<int>(signals_.size());
}

int DashboardSignalModel::renderableSelectedCount() const {
    int count = 0;
    for (int row = 0; row < signals_.size(); ++row) {
        if (isVisibleAvailable(row) && renderableModeCount(row) > 0)
            ++count;
    }
    return count;
}

int DashboardSignalModel::unavailableCount() const {
    int count = 0;
    for (int row = 0; row < availabilityByRow_.size(); ++row) {
        if (!availabilityByRow_.at(row).isVisible())
            ++count;
    }
    return count;
}

int DashboardSignalModel::noRendererCount() const {
    int count = 0;
    for (int row = 0; row < signals_.size(); ++row) {
        if (renderableModeCount(row) == 0)
            ++count;
    }
    return count;
}

QUuid DashboardSignalModel::addSignal(const DashboardSignal& input) {
    DashboardSignal signal = NormalizeSignal(input);
    if (rowForId(signal.id) >= 0)
        signal.id = QUuid::createUuid();

    const SignalAvailability availability = resolveAvailability(signal);
    const int row = static_cast<int>(signals_.size());
    beginInsertRows(QModelIndex(), row, row);
    signals_.push_back(signal);
    availabilityByRow_.push_back(availability);
    endInsertRows();
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=signal_added id=%1 descriptor_id=%2 concept_key=%3 modes=%4 anchor=%5 selected=%6")
               .arg(uuidLogValue(signal.id))
               .arg(signal.binding.descriptorId)
               .arg(signal.binding.conceptKey)
               .arg(modesLogValue(signal.displayModeIds))
               .arg(anchorLogValue(signal.binding.anchor))
               .arg(selectedIdsLogValue(signals_));
    emit signalAdded(signal.id);
    return signal.id;
}

QUuid DashboardSignalModel::addSignal(const SignalDescriptor& descriptor,
                                      const SignalAnchor& anchor,
                                      const QString& reducerId,
                                      const QStringList& displayModeIds,
                                      bool followsFocus,
                                      const QString& label) {
    DashboardSignal signal;
    signal.binding.sourceKind = descriptor.sourceKind;
    signal.binding.descriptorId = descriptor.id;
    signal.binding.conceptKey = descriptor.conceptKey;
    signal.binding.reducerId = reducerId;
    signal.binding.anchor = anchor;
    signal.binding.followsFocus = followsFocus;
    signal.nativeAxis = descriptor.nativeAxis;
    signal.requiredAnchor = descriptor.requiredAnchor;
    signal.valueShape = descriptor.valueShape;
    signal.label = label.isEmpty() ? descriptor.label : label;

    signal.displayModeIds = normalizedModeList(displayModeIds);
    if (signal.displayModeIds.isEmpty())
        signal.displayModeIds = AllDisplayModes(descriptor);
    if (!signal.displayModeIds.isEmpty())
        signal.binding.displayModeId = signal.displayModeIds.first();

    return addSignal(signal);
}

bool DashboardSignalModel::removeSignal(const QUuid& id) {
    return removeSignalAt(rowForId(id));
}

bool DashboardSignalModel::removeSignalAt(int row) {
    if (row < 0 || row >= signals_.size())
        return false;
    const QUuid id = signals_.at(row).id;
    beginRemoveRows(QModelIndex(), row, row);
    signals_.removeAt(row);
    availabilityByRow_.removeAt(row);
    endRemoveRows();
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=signal_removed id=%1 selected=%2")
               .arg(uuidLogValue(id))
               .arg(selectedIdsLogValue(signals_));
    emit signalRemoved(id);
    return true;
}

void DashboardSignalModel::clear() {
    if (signals_.isEmpty())
        return;
    QVector<QUuid> removedIds;
    removedIds.reserve(signals_.size());
    for (const DashboardSignal& signal : signals_)
        removedIds.push_back(signal.id);

    beginResetModel();
    signals_.clear();
    availabilityByRow_.clear();
    endResetModel();

    for (const QUuid& id : removedIds)
        emit signalRemoved(id);
}

bool DashboardSignalModel::updateSignal(const DashboardSignal& input) {
    if (input.id.isNull())
        return false;
    const int row = rowForId(input.id);
    if (row < 0)
        return false;
    signals_[row] = NormalizeSignal(input);
    availabilityByRow_[row] = resolveAvailability(signals_.at(row));
    emitRowChanged(row, allRoles());
    return true;
}

bool DashboardSignalModel::updateBinding(const QUuid& id, const DisplaySignalBinding& binding) {
    const int row = rowForId(id);
    if (row < 0)
        return false;
    signals_[row].binding = binding;
    if (!binding.displayModeId.isEmpty() && !signals_[row].displayModeIds.contains(binding.displayModeId))
        signals_[row].displayModeIds.push_back(binding.displayModeId);
    availabilityByRow_[row] = resolveAvailability(signals_.at(row));
    emitRowChanged(row,
                   {SourceRole,
                    SourceKindRole,
                    AxisRole,
                    AxisNameRole,
                    DisplayModesRole,
                    DescriptorIdRole,
                    ConceptKeyRole,
                    ReducerIdRole,
                    DisplayModeRole,
                    FollowsFocusRole,
                    AvailabilityStateRole,
                    AvailabilityRole,
                    AvailabilityReasonRole,
                    RenderableModeCountRole,
                    ModeRenderabilityRole});
    return true;
}

bool DashboardSignalModel::setLabel(const QUuid& id, const QString& label) {
    const int row = rowForId(id);
    if (row < 0 || signals_[row].label == label)
        return false;
    signals_[row].label = label;
    emitRowChanged(row, {Qt::DisplayRole, Qt::EditRole, LabelRole});
    return true;
}

bool DashboardSignalModel::setEnabled(const QUuid& id, bool enabled) {
    const int row = rowForId(id);
    if (row < 0 || signals_[row].enabled == enabled)
        return false;
    signals_[row].enabled = enabled;
    emitRowChanged(row, {Qt::CheckStateRole, EnabledRole});
    return true;
}

bool DashboardSignalModel::toggleEnabled(const QUuid& id) {
    const DashboardSignal* signal = signalById(id);
    if (!signal)
        return false;
    return setEnabled(id, !signal->enabled);
}

bool DashboardSignalModel::setDisplayModes(const QUuid& id, const QStringList& displayModeIds) {
    const int row = rowForId(id);
    if (row < 0)
        return false;

    const QStringList modes = normalizedModeList(displayModeIds);
    DashboardSignal& signal = signals_[row];
    if (signal.displayModeIds == modes)
        return false;

    signal.displayModeIds = modes;
    if (signal.displayModeIds.isEmpty()) {
        signal.binding.displayModeId.clear();
    } else if (!signal.displayModeIds.contains(signal.binding.displayModeId)) {
        signal.binding.displayModeId = signal.displayModeIds.first();
    }

    emitRowChanged(row, {DisplayModesRole,
                         DisplayModeRole,
                         RenderableModeCountRole,
                         ModeRenderabilityRole});
    return true;
}

bool DashboardSignalModel::addDisplayMode(const QUuid& id, const QString& displayModeId) {
    const QString mode = displayModeId.trimmed();
    if (mode.isEmpty())
        return false;
    const int row = rowForId(id);
    if (row < 0 || signals_[row].displayModeIds.contains(mode))
        return false;

    DashboardSignal& signal = signals_[row];
    signal.displayModeIds.push_back(mode);
    if (signal.binding.displayModeId.isEmpty())
        signal.binding.displayModeId = mode;
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=display_mode_toggled id=%1 mode=%2 enabled=1")
               .arg(uuidLogValue(id))
               .arg(mode);
    emitRowChanged(row, {DisplayModesRole,
                         DisplayModeRole,
                         RenderableModeCountRole,
                         ModeRenderabilityRole});
    return true;
}

bool DashboardSignalModel::removeDisplayMode(const QUuid& id, const QString& displayModeId) {
    const int row = rowForId(id);
    if (row < 0)
        return false;

    DashboardSignal& signal = signals_[row];
    const int removed = static_cast<int>(signal.displayModeIds.removeAll(displayModeId.trimmed()));
    if (removed == 0)
        return false;
    if (signal.binding.displayModeId == displayModeId)
        signal.binding.displayModeId = signal.displayModeIds.isEmpty() ? QString() : signal.displayModeIds.first();
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=display_mode_toggled id=%1 mode=%2 enabled=0")
               .arg(uuidLogValue(id))
               .arg(displayModeId.trimmed());
    emitRowChanged(row, {DisplayModesRole,
                         DisplayModeRole,
                         RenderableModeCountRole,
                         ModeRenderabilityRole});
    return true;
}

bool DashboardSignalModel::toggleDisplayMode(const QUuid& id, const QString& displayModeId) {
    const DashboardSignal* signal = signalById(id);
    if (!signal)
        return false;
    return toggleDisplayMode(id, displayModeId, !signal->displayModeIds.contains(displayModeId.trimmed()));
}

bool DashboardSignalModel::toggleDisplayMode(const QUuid& id, const QString& displayModeId, bool enabled) {
    return enabled ? addDisplayMode(id, displayModeId) : removeDisplayMode(id, displayModeId);
}

DashboardSignal DashboardSignalModel::NormalizeSignal(DashboardSignal signal) {
    if (signal.id.isNull())
        signal.id = QUuid::createUuid();

    signal.displayModeIds = normalizedModeList(signal.displayModeIds);
    if (signal.displayModeIds.isEmpty() && !signal.binding.displayModeId.trimmed().isEmpty())
        signal.displayModeIds.push_back(signal.binding.displayModeId.trimmed());

    if (!signal.displayModeIds.isEmpty()) {
        if (!signal.displayModeIds.contains(signal.binding.displayModeId))
            signal.binding.displayModeId = signal.displayModeIds.first();
    } else {
        signal.binding.displayModeId.clear();
    }

    if (signal.label.isEmpty()) {
        if (!signal.binding.conceptKey.isEmpty())
            signal.label = signal.binding.conceptKey;
        else
            signal.label = signal.binding.descriptorId;
    }

    return signal;
}

DashboardSignalModel::SignalAvailability DashboardSignalModel::resolveAvailability(
    const DashboardSignal& signal) const {
    SignalAvailability out;
    if (!availability_) {
        out.reason = availabilityReasonText(out);
        return out;
    }

    out.evaluated = true;
    const TrajectoryFieldAvailabilityRecord* record =
        availability_->recordForDescriptor(signal.binding.descriptorId);
    if (record) {
        out.state = record->state;
        out.finiteSamples = record->finiteSamples;
        out.nonZeroSamples = record->nonZeroSamples;
    }
    out.reason = record
        ? availabilityReasonText(out)
        : QStringLiteral("Descriptor was not classified by trajectory availability.");
    return out;
}

QVariantList DashboardSignalModel::modeRenderabilityVariantList(int row) const {
    QVariantList out;
    for (const ModeRenderability& mode : modeRenderability(row)) {
        out.push_back(QVariantMap{
            {QStringLiteral("mode"), mode.mode},
            {QStringLiteral("renderable"), mode.isRenderable()},
            {QStringLiteral("renderable_panel"), mode.buildsPanelWidget},
            {QStringLiteral("has_visible_surface"), mode.hasVisibleSurface},
            {QStringLiteral("is_panel_ref"), mode.emitsPanelRef},
        });
    }
    return out;
}

void DashboardSignalModel::emitRowChanged(int row, const QList<int>& roles) {
    if (row < 0 || row >= signals_.size())
        return;
    const QModelIndex changed = index(row, 0);
    emit dataChanged(changed, changed, roles);
    emit signalChanged(signals_.at(row).id);
}

void DashboardSignalModel::emitAvailabilityChanged(int row) {
    if (row < 0 || row >= signals_.size())
        return;
    const QModelIndex changed = index(row, 0);
    emit dataChanged(changed,
                     changed,
                     {AvailabilityStateRole, AvailabilityRole, AvailabilityReasonRole});
    emit signalChanged(signals_.at(row).id);
}

}  // namespace h5reader::model

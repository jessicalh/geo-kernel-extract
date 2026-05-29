#include "DashboardSignalModel.h"

#include <QStringList>

#include <algorithm>

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
    };
}

}  // namespace

DashboardSignalModel::DashboardSignalModel(QObject* parent)
    : QAbstractListModel(parent) {}

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
        emitRowChanged(row, {DisplayModesRole, DisplayModeRole});
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

QUuid DashboardSignalModel::addSignal(const DashboardSignal& input) {
    DashboardSignal signal = NormalizeSignal(input);
    if (rowForId(signal.id) >= 0)
        signal.id = QUuid::createUuid();

    const int row = static_cast<int>(signals_.size());
    beginInsertRows(QModelIndex(), row, row);
    signals_.push_back(signal);
    endInsertRows();
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
    endRemoveRows();
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
                    FollowsFocusRole});
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

    emitRowChanged(row, {DisplayModesRole, DisplayModeRole});
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
    emitRowChanged(row, {DisplayModesRole, DisplayModeRole});
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
    emitRowChanged(row, {DisplayModesRole, DisplayModeRole});
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

void DashboardSignalModel::emitRowChanged(int row, const QList<int>& roles) {
    if (row < 0 || row >= signals_.size())
        return;
    const QModelIndex changed = index(row, 0);
    emit dataChanged(changed, changed, roles);
    emit signalChanged(signals_.at(row).id);
}

}  // namespace h5reader::model

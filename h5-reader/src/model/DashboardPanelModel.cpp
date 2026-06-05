#include "DashboardPanelModel.h"

#include "../diagnostics/DashboardLogging.h"
#include "DisplayModeCapability.h"
#include "StripVisualization.h"

#include <QStringList>

#include <algorithm>

namespace h5reader::model {

namespace {

QString trimmedOrDefault(QString name, int ordinal) {
    name = name.trimmed();
    if (!name.isEmpty())
        return name;
    return DashboardPanelModel::tr("Panel %1").arg(ordinal);
}

bool modeWantsChannel(const SignalDescriptor& descriptor,
                      const QString& mode,
                      const ChannelDescriptor& channel) {
    return StripModeWantsChannel(descriptor, mode, channel);
}

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

QString uuidLogValue(const QUuid& id) {
    return id.toString(QUuid::WithoutBraces);
}

int signalReferenceCount(const QVector<DashboardPanel>& panels, const QUuid& signalId) {
    if (signalId.isNull())
        return 0;
    int count = 0;
    for (const DashboardPanel& panel : panels) {
        for (const DashboardDisplayRef& ref : panel.displays) {
            if (ref.signalId == signalId)
                ++count;
        }
    }
    return count;
}

void logPanelRefsChanged(const QVector<DashboardPanel>& panels,
                         const QUuid& panelId,
                         const DashboardDisplayRef& ref) {
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=panel_refs_changed signal_id=%1 ref_count=%2 panel_id=%3 mode=%4 channel=%5")
               .arg(uuidLogValue(ref.signalId))
               .arg(signalReferenceCount(panels, ref.signalId))
               .arg(uuidLogValue(panelId))
               .arg(ref.displayModeId)
               .arg(ref.channelId);
}

}  // namespace

QString DashboardDisplayRef::stableKey() const {
    return QStringLiteral("%1|%2|%3")
        .arg(signalId.toString(QUuid::WithoutBraces), displayModeId, channelId);
}

bool IsPanelDisplayMode(const QString& mode) {
    return DisplayModeCapabilityFor(mode).emitsPanelRef;
}

QVector<DashboardDisplayRef> DisplayRefsForSignal(const QUuid& signalId,
                                                  const SignalDescriptor& descriptor,
                                                  const QStringList& displayModeIds) {
    QVector<DashboardDisplayRef> refs;
    if (signalId.isNull())
        return refs;

    const QStringList modes = normalizedModeList(displayModeIds);
    for (const QString& mode : modes) {
        // Panel modes are tracked with a "panel" sentinel channel id so
        // the cleanup cascade (removeDisplayRefsForSignal) drops them
        // alongside any strip-mode refs the signal carries.
        if (DisplayModeCapabilityFor(mode).emitsPanelRef) {
            refs.push_back(DashboardDisplayRef{signalId, mode, QStringLiteral("panel")});
            continue;
        }
        if (!mode.startsWith(QStringLiteral("strip.")))
            continue;

        if (descriptor.channels.isEmpty()) {
            refs.push_back(DashboardDisplayRef{signalId, mode, QStringLiteral("value")});
            continue;
        }

        const int before = static_cast<int>(refs.size());
        for (const ChannelDescriptor& channel : descriptor.channels) {
            if (modeWantsChannel(descriptor, mode, channel))
                refs.push_back(DashboardDisplayRef{signalId, mode, channel.id});
        }
        if (static_cast<int>(refs.size()) == before) {
            for (const ChannelDescriptor& channel : descriptor.channels)
                refs.push_back(DashboardDisplayRef{signalId, mode, channel.id});
        }
    }
    return refs;
}

DashboardPanelModel::DashboardPanelModel(QObject* parent)
    : QAbstractListModel(parent) {
    ensureOnePanel();
}

int DashboardPanelModel::rowCount(const QModelIndex& parent) const {
    if (parent.isValid())
        return 0;
    return static_cast<int>(panels_.size());
}

QVariant DashboardPanelModel::data(const QModelIndex& index, int role) const {
    if (!index.isValid() || index.column() != 0)
        return {};
    const int row = index.row();
    if (row < 0 || row >= panels_.size())
        return {};

    const DashboardPanel& panel = panels_.at(row);
    switch (role) {
    case Qt::DisplayRole:
    case Qt::EditRole:
    case NameRole:
        return panel.name;
    case IdRole:
        return panel.id.toString(QUuid::WithoutBraces);
    case UuidRole:
        return panel.id;
    case DisplayCountRole:
        return panel.displays.size();
    case ActiveRole:
        return panel.id == activePanelId_;
    default:
        return {};
    }
}

bool DashboardPanelModel::setData(const QModelIndex& index, const QVariant& value, int role) {
    if (!index.isValid() || index.column() != 0)
        return false;
    const int row = index.row();
    if (row < 0 || row >= panels_.size())
        return false;

    switch (role) {
    case Qt::EditRole:
    case NameRole:
        return setPanelName(panels_.at(row).id, value.toString());
    default:
        return false;
    }
}

bool DashboardPanelModel::removeRows(int row, int count, const QModelIndex& parent) {
    if (parent.isValid() || row < 0 || count <= 0 || row + count > panels_.size())
        return false;

    bool ok = true;
    for (int i = 0; i < count; ++i)
        ok = removePanelAt(row) && ok;
    return ok;
}

Qt::ItemFlags DashboardPanelModel::flags(const QModelIndex& index) const {
    if (!index.isValid())
        return Qt::NoItemFlags;
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable;
}

QHash<int, QByteArray> DashboardPanelModel::roleNames() const {
    QHash<int, QByteArray> names = QAbstractListModel::roleNames();
    names[IdRole] = "id";
    names[UuidRole] = "uuid";
    names[NameRole] = "name";
    names[DisplayCountRole] = "displayCount";
    names[ActiveRole] = "active";
    return names;
}

const DashboardPanel* DashboardPanelModel::panelById(const QUuid& id) const {
    const int row = rowForId(id);
    if (row < 0)
        return nullptr;
    return &panels_.at(row);
}

DashboardPanel* DashboardPanelModel::panelById(const QUuid& id) {
    const int row = rowForId(id);
    if (row < 0)
        return nullptr;
    return &panels_[row];
}

const DashboardPanel* DashboardPanelModel::activePanel() const {
    return panelById(activePanelId_);
}

int DashboardPanelModel::activePanelRow() const {
    return rowForId(activePanelId_);
}

int DashboardPanelModel::rowForId(const QUuid& id) const {
    for (int row = 0; row < panels_.size(); ++row) {
        if (panels_.at(row).id == id)
            return row;
    }
    return -1;
}

QModelIndex DashboardPanelModel::indexForId(const QUuid& id) const {
    const int row = rowForId(id);
    if (row < 0)
        return {};
    return index(row, 0);
}

QUuid DashboardPanelModel::addPanel(const QString& name) {
    const int row = static_cast<int>(panels_.size());
    DashboardPanel panel = makePanel(name, row + 1);

    beginInsertRows(QModelIndex(), row, row);
    panels_.push_back(panel);
    endInsertRows();

    if (activePanelId_.isNull())
        activePanelId_ = panel.id;
    emit panelAdded(panel.id);
    return panel.id;
}

bool DashboardPanelModel::removePanel(const QUuid& id) {
    return removePanelAt(rowForId(id));
}

bool DashboardPanelModel::removePanelAt(int row) {
    if (row < 0 || row >= panels_.size())
        return false;

    if (panels_.size() == 1) {
        DashboardPanel& panel = panels_[row];
        const QVector<DashboardDisplayRef> removed = panel.displays;
        if (removed.isEmpty())
            return false;
        panel.displays.clear();
        emitPanelRowChanged(row, {DisplayCountRole});
        for (const DashboardDisplayRef& ref : removed)
            emit displayRefRemoved(panel.id, ref);
        emit displayRefsChanged(panel.id);
        emit panelRemoved(panel.id, removed);
        return true;
    }

    const DashboardPanel removedPanel = panels_.at(row);
    const bool removedWasActive = removedPanel.id == activePanelId_;
    beginRemoveRows(QModelIndex(), row, row);
    panels_.removeAt(row);
    endRemoveRows();

    if (removedWasActive) {
        const int nextRow = std::min(row, static_cast<int>(panels_.size()) - 1);
        activePanelId_ = panels_.at(nextRow).id;
        emit activePanelChanged(activePanelId_);
    }
    for (const DashboardDisplayRef& ref : removedPanel.displays)
        emit displayRefRemoved(removedPanel.id, ref);
    emit panelRemoved(removedPanel.id, removedPanel.displays);
    return true;
}

void DashboardPanelModel::clear() {
    const QVector<DashboardPanel> removedPanels = panels_;

    beginResetModel();
    panels_.clear();
    activePanelId_ = {};
    ensureOnePanel();
    endResetModel();

    for (const DashboardPanel& panel : removedPanels) {
        for (const DashboardDisplayRef& ref : panel.displays)
            emit displayRefRemoved(panel.id, ref);
        if (!panel.displays.isEmpty())
            emit displayRefsChanged(panel.id);
        emit panelRemoved(panel.id, panel.displays);
    }
    emit activePanelChanged(activePanelId_);
}

bool DashboardPanelModel::setActivePanel(const QUuid& id) {
    const int row = rowForId(id);
    if (row < 0 || activePanelId_ == id)
        return false;

    const int oldRow = activePanelRow();
    activePanelId_ = id;
    if (oldRow >= 0)
        emitPanelRowChanged(oldRow, {ActiveRole});
    emitPanelRowChanged(row, {ActiveRole});
    emit activePanelChanged(activePanelId_);
    return true;
}

bool DashboardPanelModel::setActivePanelAt(int row) {
    if (row < 0 || row >= panels_.size())
        return false;
    return setActivePanel(panels_.at(row).id);
}

bool DashboardPanelModel::setPanelName(const QUuid& id, const QString& name) {
    const int row = rowForId(id);
    if (row < 0)
        return false;
    const QString normalized = trimmedOrDefault(name, row + 1);
    if (panels_[row].name == normalized)
        return false;
    panels_[row].name = normalized;
    emitPanelRowChanged(row, {Qt::DisplayRole, Qt::EditRole, NameRole});
    emit panelChanged(id);
    return true;
}

bool DashboardPanelModel::addDisplayRef(const QUuid& panelId, const DashboardDisplayRef& ref) {
    if (ref.signalId.isNull() || ref.displayModeId.isEmpty() || ref.channelId.isEmpty())
        return false;
    const int row = rowForId(panelId);
    if (row < 0 || panels_[row].displays.contains(ref))
        return false;
    panels_[row].displays.push_back(ref);
    emitPanelRowChanged(row, {DisplayCountRole});
    logPanelRefsChanged(panels_, panelId, ref);
    emit displayRefAdded(panelId, ref);
    emit displayRefsChanged(panelId);
    return true;
}

bool DashboardPanelModel::addDisplayRefs(const QUuid& panelId, const QVector<DashboardDisplayRef>& refs) {
    bool changed = false;
    for (const DashboardDisplayRef& ref : refs)
        changed = addDisplayRef(panelId, ref) || changed;
    return changed;
}

bool DashboardPanelModel::removeDisplayRef(const QUuid& panelId, const DashboardDisplayRef& ref) {
    const int row = rowForId(panelId);
    if (row < 0)
        return false;
    const int index = static_cast<int>(panels_[row].displays.indexOf(ref));
    if (index < 0)
        return false;
    panels_[row].displays.removeAt(index);
    emitPanelRowChanged(row, {DisplayCountRole});
    logPanelRefsChanged(panels_, panelId, ref);
    emit displayRefRemoved(panelId, ref);
    emit displayRefsChanged(panelId);
    return true;
}

int DashboardPanelModel::removeDisplayRefsForSignal(const QUuid& signalId) {
    if (signalId.isNull())
        return 0;
    int removedCount = 0;
    for (int row = 0; row < panels_.size(); ++row) {
        DashboardPanel& panel = panels_[row];
        bool rowChanged = false;
        for (int i = static_cast<int>(panel.displays.size()) - 1; i >= 0; --i) {
            if (panel.displays.at(i).signalId != signalId)
                continue;
            const DashboardDisplayRef ref = panel.displays.at(i);
            panel.displays.removeAt(i);
            ++removedCount;
            rowChanged = true;
            logPanelRefsChanged(panels_, panel.id, ref);
            emit displayRefRemoved(panel.id, ref);
        }
        if (rowChanged) {
            emitPanelRowChanged(row, {DisplayCountRole});
            emit displayRefsChanged(panel.id);
        }
    }
    return removedCount;
}

int DashboardPanelModel::removeDisplayRefsForSignalMode(const QUuid& signalId, const QString& displayModeId) {
    const QString mode = displayModeId.trimmed();
    if (signalId.isNull() || mode.isEmpty())
        return 0;
    int removedCount = 0;
    for (int row = 0; row < panels_.size(); ++row) {
        DashboardPanel& panel = panels_[row];
        bool rowChanged = false;
        for (int i = static_cast<int>(panel.displays.size()) - 1; i >= 0; --i) {
            const DashboardDisplayRef& candidate = panel.displays.at(i);
            if (candidate.signalId != signalId || candidate.displayModeId != mode)
                continue;
            const DashboardDisplayRef ref = candidate;
            panel.displays.removeAt(i);
            ++removedCount;
            rowChanged = true;
            logPanelRefsChanged(panels_, panel.id, ref);
            emit displayRefRemoved(panel.id, ref);
        }
        if (rowChanged) {
            emitPanelRowChanged(row, {DisplayCountRole});
            emit displayRefsChanged(panel.id);
        }
    }
    return removedCount;
}

int DashboardPanelModel::signalReferenceCount(const QUuid& signalId) const {
    if (signalId.isNull())
        return 0;
    int count = 0;
    for (const DashboardPanel& panel : panels_) {
        for (const DashboardDisplayRef& ref : panel.displays) {
            if (ref.signalId == signalId)
                ++count;
        }
    }
    return count;
}

bool DashboardPanelModel::containsDisplayRef(const QUuid& panelId, const DashboardDisplayRef& ref) const {
    const DashboardPanel* panel = panelById(panelId);
    return panel && panel->displays.contains(ref);
}

QVector<DashboardDisplayRef> DashboardPanelModel::displayRefsForPanel(const QUuid& panelId) const {
    const DashboardPanel* panel = panelById(panelId);
    return panel ? panel->displays : QVector<DashboardDisplayRef>{};
}

QString DashboardPanelModel::defaultPanelName(int ordinal) {
    return QStringLiteral("Panel %1").arg(std::max(1, ordinal));
}

DashboardPanel DashboardPanelModel::makePanel(const QString& name, int ordinal) {
    DashboardPanel panel;
    panel.id = QUuid::createUuid();
    panel.name = name.trimmed().isEmpty() ? defaultPanelName(ordinal) : name.trimmed();
    return panel;
}

void DashboardPanelModel::ensureOnePanel() {
    if (!panels_.isEmpty()) {
        if (activePanelId_.isNull())
            activePanelId_ = panels_.first().id;
        return;
    }
    panels_.push_back(makePanel(QStringLiteral("Dashboard"), 1));
    activePanelId_ = panels_.first().id;
}

void DashboardPanelModel::emitPanelRowChanged(int row, const QList<int>& roles) {
    if (row < 0 || row >= panels_.size())
        return;
    const QModelIndex changed = index(row, 0);
    emit dataChanged(changed, changed, roles);
}

}  // namespace h5reader::model

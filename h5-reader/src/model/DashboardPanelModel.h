#pragma once

#include "DashboardSignal.h"

#include <QAbstractListModel>
#include <QHash>
#include <QList>
#include <QModelIndex>
#include <QUuid>
#include <QVector>

namespace h5reader::model {

struct DashboardDisplayRef {
    QUuid signalId;
    QString displayModeId;
    QString channelId;

    QString stableKey() const;

    friend bool operator==(const DashboardDisplayRef& a, const DashboardDisplayRef& b) {
        return a.signalId == b.signalId
               && a.displayModeId == b.displayModeId
               && a.channelId == b.channelId;
    }
};

struct DashboardPanel {
    QUuid id;
    QString name;
    QVector<DashboardDisplayRef> displays;
};

QVector<DashboardDisplayRef> DisplayRefsForSignal(const QUuid& signalId,
                                                  const SignalDescriptor& descriptor,
                                                  const QStringList& displayModeIds);

class DashboardPanelModel final : public QAbstractListModel {
    Q_OBJECT

public:
    enum Role {
        IdRole = Qt::UserRole + 1,
        UuidRole,
        NameRole,
        DisplayCountRole,
        ActiveRole,
    };
    Q_ENUM(Role)

    explicit DashboardPanelModel(QObject* parent = nullptr);

    int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;
    bool setData(const QModelIndex& index, const QVariant& value, int role = Qt::EditRole) override;
    bool removeRows(int row, int count, const QModelIndex& parent = QModelIndex()) override;
    Qt::ItemFlags flags(const QModelIndex& index) const override;
    QHash<int, QByteArray> roleNames() const override;

    const QVector<DashboardPanel>& panels() const { return panels_; }
    const DashboardPanel* panelById(const QUuid& id) const;
    DashboardPanel* panelById(const QUuid& id);
    const DashboardPanel* activePanel() const;
    QUuid activePanelId() const { return activePanelId_; }
    int activePanelRow() const;
    int rowForId(const QUuid& id) const;
    QModelIndex indexForId(const QUuid& id) const;

    QUuid addPanel(const QString& name = QString());
    bool removePanel(const QUuid& id);
    bool removePanelAt(int row);
    void clear();

    bool setActivePanel(const QUuid& id);
    bool setActivePanelAt(int row);
    bool setPanelName(const QUuid& id, const QString& name);

    bool addDisplayRef(const QUuid& panelId, const DashboardDisplayRef& ref);
    bool addDisplayRefs(const QUuid& panelId, const QVector<DashboardDisplayRef>& refs);
    bool removeDisplayRef(const QUuid& panelId, const DashboardDisplayRef& ref);
    int removeDisplayRefsForSignal(const QUuid& signalId);
    int removeDisplayRefsForSignalMode(const QUuid& signalId, const QString& displayModeId);
    int signalReferenceCount(const QUuid& signalId) const;
    bool containsDisplayRef(const QUuid& panelId, const DashboardDisplayRef& ref) const;
    QVector<DashboardDisplayRef> displayRefsForPanel(const QUuid& panelId) const;

signals:
    void panelAdded(const QUuid& id);
    void panelRemoved(const QUuid& id, const QVector<h5reader::model::DashboardDisplayRef>& removedRefs);
    void panelChanged(const QUuid& id);
    void activePanelChanged(const QUuid& id);
    void displayRefsChanged(const QUuid& panelId);
    void displayRefAdded(const QUuid& panelId, const h5reader::model::DashboardDisplayRef& ref);
    void displayRefRemoved(const QUuid& panelId, const h5reader::model::DashboardDisplayRef& ref);

private:
    static QString defaultPanelName(int ordinal);
    static DashboardPanel makePanel(const QString& name, int ordinal);
    void ensureOnePanel();
    void emitPanelRowChanged(int row, const QList<int>& roles);

    QVector<DashboardPanel> panels_;
    QUuid activePanelId_;
};

}  // namespace h5reader::model

#pragma once

#include "DashboardSignal.h"

#include <QAbstractListModel>
#include <QModelIndex>
#include <QVector>

namespace h5reader::model {

class DashboardSignalModel final : public QAbstractListModel {
    Q_OBJECT

public:
    enum Role {
        IdRole = Qt::UserRole + 1,
        UuidRole,
        LabelRole,
        SourceRole,
        SourceKindRole,
        AxisRole,
        AxisNameRole,
        ValueShapeRole,
        ValueShapeNameRole,
        DisplayModesRole,
        EnabledRole,
        DescriptorIdRole,
        ConceptKeyRole,
        ReducerIdRole,
        DisplayModeRole,
        FollowsFocusRole,
    };
    Q_ENUM(Role)

    explicit DashboardSignalModel(QObject* parent = nullptr);

    int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;
    bool setData(const QModelIndex& index, const QVariant& value, int role = Qt::EditRole) override;
    bool removeRows(int row, int count, const QModelIndex& parent = QModelIndex()) override;
    Qt::ItemFlags flags(const QModelIndex& index) const override;
    QHash<int, QByteArray> roleNames() const override;

    const QVector<DashboardSignal>& activeSignals() const { return signals_; }
    const DashboardSignal* signalById(const QUuid& id) const;
    DashboardSignal* signalById(const QUuid& id);
    QModelIndex indexForId(const QUuid& id) const;
    int rowForId(const QUuid& id) const;

    QUuid addSignal(const DashboardSignal& signal);
    QUuid addSignal(const SignalDescriptor& descriptor,
                    const SignalAnchor& anchor,
                    const QString& reducerId = QString(),
                    const QStringList& displayModeIds = QStringList(),
                    bool followsFocus = false,
                    const QString& label = QString());

    bool removeSignal(const QUuid& id);
    bool removeSignalAt(int row);
    void clear();

    bool updateSignal(const DashboardSignal& signal);
    bool updateBinding(const QUuid& id, const DisplaySignalBinding& binding);
    bool setLabel(const QUuid& id, const QString& label);
    bool setEnabled(const QUuid& id, bool enabled);
    bool toggleEnabled(const QUuid& id);
    bool setDisplayModes(const QUuid& id, const QStringList& displayModeIds);
    bool addDisplayMode(const QUuid& id, const QString& displayModeId);
    bool removeDisplayMode(const QUuid& id, const QString& displayModeId);
    bool toggleDisplayMode(const QUuid& id, const QString& displayModeId);
    bool toggleDisplayMode(const QUuid& id, const QString& displayModeId, bool enabled);

signals:
    void signalAdded(const QUuid& id);
    void signalRemoved(const QUuid& id);
    void signalChanged(const QUuid& id);

private:
    static DashboardSignal NormalizeSignal(DashboardSignal signal);
    void emitRowChanged(int row, const QList<int>& roles);

    QVector<DashboardSignal> signals_;
};

}  // namespace h5reader::model

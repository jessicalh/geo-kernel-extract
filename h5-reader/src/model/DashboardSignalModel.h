#pragma once

#include "DashboardSignal.h"
#include "TrajectoryFieldAvailability.h"

#include <QAbstractListModel>
#include <QModelIndex>
#include <QString>
#include <QVector>

#include <memory>

namespace h5reader::model {

class DashboardSignalModel final : public QAbstractListModel {
    Q_OBJECT

public:
    struct SignalAvailability {
        bool evaluated = false;
        TrajectoryFieldAvailabilityState state = TrajectoryFieldAvailabilityState::Available;
        QString reason;
        qsizetype finiteSamples = 0;
        qsizetype nonZeroSamples = 0;

        bool isVisible() const;
        bool operator==(const SignalAvailability& other) const;
        bool operator!=(const SignalAvailability& other) const { return !(*this == other); }
    };

    struct ModeRenderability {
        QString mode;
        bool hasVisibleSurface = false;
        bool buildsPanelWidget = false;
        bool emitsPanelRef = false;

        bool isRenderable() const { return hasVisibleSurface || buildsPanelWidget; }
    };

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
        AvailabilityStateRole,
        AvailabilityRole,
        AvailabilityReasonRole,
        RenderableModeCountRole,
        ModeRenderabilityRole,
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

    void setFieldAvailability(std::shared_ptr<const TrajectoryFieldAvailability> availability);
    std::shared_ptr<const TrajectoryFieldAvailability> fieldAvailability() const { return availability_; }
    void refreshAvailability();

    SignalAvailability availabilityAt(int row) const;
    TrajectoryFieldAvailabilityState availabilityState(int row) const;
    QString availabilityName(int row) const;
    QString availabilityReason(int row) const;
    bool isVisibleAvailable(int row) const;

    static ModeRenderability ModeRenderabilityFor(const QString& mode);
    QVector<ModeRenderability> modeRenderability(int row) const;
    int renderableModeCount(int row) const;

    int selectedCount() const;
    int renderableSelectedCount() const;
    int unavailableCount() const;
    int noRendererCount() const;

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
    SignalAvailability resolveAvailability(const DashboardSignal& signal) const;
    QVariantList modeRenderabilityVariantList(int row) const;
    void emitRowChanged(int row, const QList<int>& roles);
    void emitAvailabilityChanged(int row);

    QVector<DashboardSignal> signals_;
    QVector<SignalAvailability> availabilityByRow_;
    std::shared_ptr<const TrajectoryFieldAvailability> availability_;
};

}  // namespace h5reader::model

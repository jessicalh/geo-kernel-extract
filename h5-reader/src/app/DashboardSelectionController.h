#pragma once

#include "../model/DashboardSignal.h"

#include <QObject>
#include <QPointer>
#include <QSet>
#include <QString>
#include <QStringList>
#include <QUuid>

namespace h5reader::model {
class DashboardPanelModel;
struct DashboardDisplayRef;
class DashboardSignalModel;
class TrajectorySignalCatalog;
}

namespace h5reader::app {

class DashboardSelectionController final : public QObject {
    Q_OBJECT

public:
    enum class PanelRemovalPolicy {
        RemoveReferencedMetrics,
    };
    Q_ENUM(PanelRemovalPolicy)

    struct PanelTarget {
        QUuid panelId;
        QString newPanelName;
        bool makeActive = false;
    };

    explicit DashboardSelectionController(model::TrajectorySignalCatalog* catalog,
                                          model::DashboardSignalModel* signalModel,
                                          model::DashboardPanelModel* panelModel,
                                          QObject* parent = nullptr);
    ~DashboardSelectionController() override = default;

    QUuid addMetric(const model::SignalDescriptor& descriptor,
                    const model::SignalAnchor& anchor,
                    const QStringList& modes,
                    const PanelTarget& target,
                    bool followsFocus = false,
                    const QString& label = QString(),
                    int* addedRefs = nullptr);
    bool setMetricMode(const QUuid& signalId,
                       const QString& mode,
                       bool enabled);
    bool removeMetric(const QUuid& signalId);
    bool setMetricMode(const QUuid& signalId,
                       const QString& mode,
                       bool enabled,
                       const PanelTarget& target);
    int removePanel(const QUuid& panelId, PanelRemovalPolicy policy);
    void clearAllMetrics();

    int selectedCount() const;

signals:
    void selectedCountChanged(int count);
    void selectionChanged();

private:
    void onSignalAdded(const QUuid& id);
    void onSignalRemoved(const QUuid& id);
    void onDisplayRefRemoved(const model::DashboardDisplayRef& ref);

    QUuid resolvePanelTarget(const PanelTarget& target, bool allowCreate);
    const model::SignalDescriptor* descriptorForSignal(const QUuid& signalId) const;
    int removeAllDisplayRefs();
    void emitSelectedCountIfChanged();

    QPointer<model::TrajectorySignalCatalog> catalog_;
    QPointer<model::DashboardSignalModel> signals_;
    QPointer<model::DashboardPanelModel> panels_;
    QSet<QUuid> signalsBeingRemoved_;
    int lastSelectedCount_ = 0;
};

}  // namespace h5reader::app

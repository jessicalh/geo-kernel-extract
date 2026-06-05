#include "DashboardSelectionController.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/DashboardLogging.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/DashboardPanelModel.h"
#include "../model/DashboardSignalModel.h"
#include "../model/TrajectorySignalCatalog.h"

#include <QAbstractItemModel>

#include <algorithm>

namespace h5reader::app {

namespace {

QString uuidLogValue(const QUuid& id) {
    return id.isNull() ? QStringLiteral("null") : id.toString(QUuid::WithoutBraces);
}

QString modesLogValue(const QStringList& modes) {
    return QStringLiteral("[%1]").arg(modes.join(QStringLiteral(",")));
}

}  // namespace

DashboardSelectionController::DashboardSelectionController(model::TrajectorySignalCatalog* catalog,
                                                           model::DashboardSignalModel* signalModel,
                                                           model::DashboardPanelModel* panelModel,
                                                           QObject* parent)
    : QObject(parent)
    , catalog_(catalog)
    , signals_(signalModel)
    , panels_(panelModel)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("DashboardSelectionController"));

    lastSelectedCount_ = selectedCount();
    if (signals_) {
        ACONNECT(signals_.data(), &model::DashboardSignalModel::signalAdded,
                 this, &DashboardSelectionController::onSignalAdded);
        ACONNECT(signals_.data(), &model::DashboardSignalModel::signalRemoved,
                 this, &DashboardSelectionController::onSignalRemoved);
        ACONNECT(signals_.data(), &QAbstractItemModel::modelReset,
                 this, &DashboardSelectionController::emitSelectedCountIfChanged);
    }
    if (panels_) {
        ACONNECT(panels_.data(), &model::DashboardPanelModel::displayRefRemoved,
                 this, [this](const QUuid&, const model::DashboardDisplayRef& ref) {
                     onDisplayRefRemoved(ref);
                 });
    }
}

QUuid DashboardSelectionController::addMetric(const model::SignalDescriptor& descriptor,
                                              const model::SignalAnchor& anchor,
                                              const QStringList& modes,
                                              const PanelTarget& target,
                                              bool followsFocus,
                                              const QString& label,
                                              int* addedRefs) {
    ASSERT_THREAD(this);
    if (addedRefs)
        *addedRefs = 0;
    if (!signals_)
        return {};

    const QUuid panelId = resolvePanelTarget(target, true);
    const QUuid id = signals_->addSignal(descriptor,
                                         anchor,
                                         QString(),
                                         modes,
                                         followsFocus,
                                         label);

    int refsAdded = 0;
    if (panels_ && !panelId.isNull()) {
        const QVector<model::DashboardDisplayRef> refs =
            model::DisplayRefsForSignal(id, descriptor, modes);
        for (const model::DashboardDisplayRef& ref : refs) {
            if (panels_->addDisplayRef(panelId, ref))
                ++refsAdded;
        }
        if (target.makeActive)
            panels_->setActivePanel(panelId);
    }
    if (addedRefs)
        *addedRefs = refsAdded;

    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=selection_add id=%1 count=%2 panel_id=%3 added_refs=%4 modes=%5")
               .arg(uuidLogValue(id))
               .arg(selectedCount())
               .arg(uuidLogValue(panelId))
               .arg(refsAdded)
               .arg(modesLogValue(modes));
    return id;
}

bool DashboardSelectionController::removeMetric(const QUuid& signalId) {
    ASSERT_THREAD(this);
    const bool removed = signals_ && signals_->removeSignal(signalId);
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=selection_remove id=%1 count=%2 removed=%3")
               .arg(uuidLogValue(signalId))
               .arg(selectedCount())
               .arg(removed ? 1 : 0);
    return removed;
}

bool DashboardSelectionController::setMetricMode(const QUuid& signalId,
                                                 const QString& mode,
                                                 bool enabled) {
    return setMetricMode(signalId, mode, enabled, PanelTarget{});
}

bool DashboardSelectionController::setMetricMode(const QUuid& signalId,
                                                 const QString& mode,
                                                 bool enabled,
                                                 const PanelTarget& target) {
    ASSERT_THREAD(this);
    const QString normalizedMode = mode.trimmed();
    if (!signals_ || signalId.isNull() || normalizedMode.isEmpty())
        return false;

    const model::DashboardSignal* before = signals_->signalById(signalId);
    if (!before)
        return false;

    const bool alreadyDesired = before->displayModeIds.contains(normalizedMode) == enabled;
    if (!alreadyDesired && !signals_->toggleDisplayMode(signalId, normalizedMode, enabled))
        return false;

    int refsChanged = 0;
    const model::SignalDescriptor* descriptor = descriptorForSignal(signalId);
    if (panels_ && descriptor) {
        const QUuid panelId = resolvePanelTarget(target, enabled);
        if (!panelId.isNull()) {
            if (enabled) {
                const QVector<model::DashboardDisplayRef> refs =
                    model::DisplayRefsForSignal(signalId, *descriptor, {normalizedMode});
                for (const model::DashboardDisplayRef& ref : refs) {
                    if (panels_->addDisplayRef(panelId, ref))
                        ++refsChanged;
                }
                if (target.makeActive)
                    panels_->setActivePanel(panelId);
            } else {
                refsChanged = panels_->removeDisplayRefsForSignalMode(signalId, normalizedMode);
            }
        }
    }

    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=selection_set_mode id=%1 count=%2 mode=%3 enabled=%4 refs_changed=%5")
               .arg(uuidLogValue(signalId))
               .arg(selectedCount())
               .arg(normalizedMode)
               .arg(enabled ? 1 : 0)
               .arg(refsChanged);
    return true;
}

int DashboardSelectionController::removePanel(const QUuid& panelId, PanelRemovalPolicy policy) {
    ASSERT_THREAD(this);
    if (policy != PanelRemovalPolicy::RemoveReferencedMetrics || !panels_) {
        qCInfo(diagnostics::cDash).noquote()
            << QStringLiteral("event=selection_remove_panel id=%1 count=%2 removed_metrics=0 removed=0")
                   .arg(uuidLogValue(panelId))
                   .arg(selectedCount());
        return 0;
    }

    const int before = selectedCount();
    const bool removed = panels_->removePanel(panelId);
    const int removedMetrics = removed ? std::max(0, before - selectedCount()) : 0;
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=selection_remove_panel id=%1 count=%2 removed_metrics=%3 removed=%4")
               .arg(uuidLogValue(panelId))
               .arg(selectedCount())
               .arg(removedMetrics)
               .arg(removed ? 1 : 0);
    return removedMetrics;
}

void DashboardSelectionController::clearAllMetrics() {
    ASSERT_THREAD(this);
    const int refsRemoved = removeAllDisplayRefs();
    if (signals_)
        signals_->clear();
    emitSelectedCountIfChanged();
    qCInfo(diagnostics::cDash).noquote()
        << QStringLiteral("event=selection_clear id=all count=%1 refs_removed=%2")
               .arg(selectedCount())
               .arg(refsRemoved);
}

int DashboardSelectionController::selectedCount() const {
    return signals_ ? signals_->selectedCount() : 0;
}

void DashboardSelectionController::onSignalAdded(const QUuid&) {
    ASSERT_THREAD(this);
    emitSelectedCountIfChanged();
}

void DashboardSelectionController::onSignalRemoved(const QUuid& id) {
    ASSERT_THREAD(this);
    if (panels_ && !id.isNull()) {
        signalsBeingRemoved_.insert(id);
        panels_->removeDisplayRefsForSignal(id);
        signalsBeingRemoved_.remove(id);
    }
    emitSelectedCountIfChanged();
}

void DashboardSelectionController::onDisplayRefRemoved(const model::DashboardDisplayRef& ref) {
    ASSERT_THREAD(this);
    if (!signals_ || !panels_ || ref.signalId.isNull())
        return;
    if (signalsBeingRemoved_.contains(ref.signalId) || !signals_->signalById(ref.signalId))
        return;
    if (panels_->signalReferenceCount(ref.signalId) == 0)
        signals_->removeSignal(ref.signalId);
}

QUuid DashboardSelectionController::resolvePanelTarget(const PanelTarget& target, bool allowCreate) {
    if (!panels_)
        return {};
    const QString newPanelName = target.newPanelName.trimmed();
    if (allowCreate && !newPanelName.isEmpty()) {
        const QUuid panelId = panels_->addPanel(newPanelName);
        if (target.makeActive)
            panels_->setActivePanel(panelId);
        return panelId;
    }
    if (!target.panelId.isNull())
        return target.panelId;
    return panels_->activePanelId();
}

const model::SignalDescriptor* DashboardSelectionController::descriptorForSignal(const QUuid& signalId) const {
    if (!catalog_ || !signals_)
        return nullptr;
    const model::DashboardSignal* signal = signals_->signalById(signalId);
    if (!signal)
        return nullptr;
    return catalog_->findDescriptor(signal->binding.descriptorId);
}

int DashboardSelectionController::removeAllDisplayRefs() {
    if (!panels_)
        return 0;
    int removed = 0;
    const QVector<model::DashboardPanel> panels = panels_->panels();
    for (const model::DashboardPanel& panel : panels) {
        for (const model::DashboardDisplayRef& ref : panel.displays) {
            if (panels_->removeDisplayRef(panel.id, ref))
                ++removed;
        }
    }
    return removed;
}

void DashboardSelectionController::emitSelectedCountIfChanged() {
    const int count = selectedCount();
    if (count == lastSelectedCount_)
        return;
    lastSelectedCount_ = count;
    emit selectedCountChanged(count);
    emit selectionChanged();
}

}  // namespace h5reader::app

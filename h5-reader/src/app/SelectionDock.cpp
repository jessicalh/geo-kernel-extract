#include "SelectionDock.h"

#include "../model/AtomSelection.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <QLabel>
#include <QListView>
#include <QVBoxLayout>
#include <QWidget>

namespace h5reader::app {

SelectionDock::SelectionDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Selection"), parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("SelectionDock"));
    setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);

    auto* container = new QWidget(this);
    auto* vbox      = new QVBoxLayout(container);
    vbox->setContentsMargins(4, 4, 4, 4);
    vbox->setSpacing(4);

    header_ = new QLabel(QStringLiteral("No atoms selected"), container);
    header_->setStyleSheet(QStringLiteral("font-weight: bold;"));
    vbox->addWidget(header_);

    list_ = new QListView(container);
    // Display-only in increment 1 — picking happens in the 3-D view, not the
    // panel. A later increment turns row clicks into focus changes.
    list_->setSelectionMode(QAbstractItemView::NoSelection);
    list_->setUniformItemSizes(true);
    vbox->addWidget(list_, 1);

    setWidget(container);
}

void SelectionDock::setModel(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    selection_ = selection;
    if (list_)
        list_->setModel(selection);
    if (selection) {
        ACONNECT(selection, &model::AtomSelection::changed, this, &SelectionDock::refreshHeader);
        ACONNECT(selection, &model::AtomSelection::cleared, this, &SelectionDock::refreshHeader);
    }
    refreshHeader();
}

void SelectionDock::refreshHeader() {
    ASSERT_THREAD(this);
    if (!header_)
        return;
    if (!selection_ || selection_->empty()) {
        header_->setText(QStringLiteral("No atoms selected"));
        return;
    }
    const std::size_t        n = selection_->count();
    const model::GeometryKind k = selection_->geometryKind();
    if (k == model::GeometryKind::None) {
        // One atom: defines no measurement yet — prompt how to grow the set.
        header_->setText(QStringLiteral("%1 atom · Shift-click to add").arg(n));
    } else {
        header_->setText(
            QStringLiteral("%1 atoms · %2").arg(n).arg(QString::fromLatin1(model::NameForGeometryKind(k))));
    }
}

}  // namespace h5reader::app

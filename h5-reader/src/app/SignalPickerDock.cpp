#include "SignalPickerDock.h"

#include "NearbySignalModel.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/AtomSelection.h"
#include "../model/Conformation.h"
#include "../model/DftShieldingStore.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/SignalDictionary.h"

#include <QAbstractItemView>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QItemSelectionModel>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QTableView>
#include <QVBoxLayout>
#include <QWidget>

#include <algorithm>

namespace h5reader::app {

namespace {
bool isNumericSignal(model::SignalValueKind kind) {
    return kind == model::SignalValueKind::Scalar
        || kind == model::SignalValueKind::Angle
        || kind == model::SignalValueKind::Distance
        || kind == model::SignalValueKind::Power;
}
}  // namespace

SignalPickerDock::SignalPickerDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Signal Picker"), parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("SignalPickerDock"));
    setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);

    auto* container = new QWidget(this);
    auto* vbox = new QVBoxLayout(container);
    vbox->setContentsMargins(4, 4, 4, 4);
    vbox->setSpacing(4);

    focusLabel_ = new QLabel(QStringLiteral("Focus: none"), container);
    focusLabel_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    vbox->addWidget(focusLabel_);

    auto* row = new QHBoxLayout;
    liveBox_ = new QCheckBox(QStringLiteral("Live"), container);
    liveBox_->setChecked(true);
    liveBox_->setToolTip(QStringLiteral("On: neighborhood follows the focused atom and current frame. Off: freeze the current candidate set."));
    row->addWidget(liveBox_);
    row->addWidget(new QLabel(QStringLiteral("Radius"), container));
    radiusSpin_ = new QDoubleSpinBox(container);
    radiusSpin_->setRange(1.0, 30.0);
    radiusSpin_->setDecimals(1);
    radiusSpin_->setSingleStep(0.5);
    radiusSpin_->setSuffix(QString::fromUtf8(" Å"));
    radiusSpin_->setValue(5.0);
    row->addWidget(radiusSpin_);
    row->addStretch(1);
    vbox->addLayout(row);

    nearbyModel_ = new NearbySignalModel(this);

    candidatesView_ = new QTableView(container);
    candidatesView_->setModel(nearbyModel_);
    candidatesView_->setSelectionBehavior(QAbstractItemView::SelectRows);
    candidatesView_->setSelectionMode(QAbstractItemView::SingleSelection);
    candidatesView_->setAlternatingRowColors(true);
    candidatesView_->setSortingEnabled(false);
    candidatesView_->verticalHeader()->hide();
    candidatesView_->horizontalHeader()->setStretchLastSection(false);
    candidatesView_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    candidatesView_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    candidatesView_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
    vbox->addWidget(candidatesView_, 3);

    stripList_ = new QListWidget(container);
    stripList_->setSelectionMode(QAbstractItemView::SingleSelection);
    vbox->addWidget(stripList_, 2);

    addButton_ = new QPushButton(QStringLiteral("Add Strip"), container);
    addButton_->setEnabled(false);
    addButton_->setToolTip(QStringLiteral("Dashboard pinning is the next increment; this picker only discovers bindings for now."));
    vbox->addWidget(addButton_);

    setWidget(container);

    ACONNECT(liveBox_.data(), &QCheckBox::toggled, this, &SignalPickerDock::onLiveToggled);
    ACONNECT(radiusSpin_.data(), qOverload<double>(&QDoubleSpinBox::valueChanged),
             this, &SignalPickerDock::onRadiusChanged);
    ACONNECT(candidatesView_->selectionModel(), &QItemSelectionModel::currentRowChanged,
             this, [this](const QModelIndex&, const QModelIndex&) { onCandidateChanged(); });
    ACONNECT(addButton_.data(), &QPushButton::clicked, this, &SignalPickerDock::onAddClicked);
}

void SignalPickerDock::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    ASSERT_THREAD(this);
    protein_ = protein;
    conformation_ = conformation;
    nearbyModel_->setContext(protein, conformation);
    refreshFocusLabel();
}

void SignalPickerDock::setSelection(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    if (selection_)
        disconnect(selection_, nullptr, this, nullptr);
    selection_ = selection;
    if (selection_) {
        ACONNECT(selection_, &model::AtomSelection::focusChanged, this, &SignalPickerDock::onFocusChanged);
        ACONNECT(selection_, &model::AtomSelection::cleared, this, &SignalPickerDock::onSelectionCleared);
        if (selection_->hasFocus())
            onFocusChanged(selection_->focus());
    }
    refreshFocusLabel();
}

void SignalPickerDock::setDftStore(model::DftShieldingStore* store) {
    ASSERT_THREAD(this);
    dftStore_ = store;
    refreshStripChoices();
}

void SignalPickerDock::setFrame(int frame) {
    ASSERT_THREAD(this);
    frame_ = std::max(0, frame);
    if (liveBox_ && liveBox_->isChecked())
        refreshCandidateAnchor();
}

void SignalPickerDock::onFocusChanged(std::size_t atomIdx) {
    ASSERT_THREAD(this);
    latestFocusAtom_ = atomIdx;
    refreshFocusLabel();
    if (liveBox_ && liveBox_->isChecked())
        refreshCandidateAnchor();
}

void SignalPickerDock::onSelectionCleared() {
    ASSERT_THREAD(this);
    latestFocusAtom_.reset();
    refreshFocusLabel();
    if (liveBox_ && liveBox_->isChecked())
        nearbyModel_->clear();
    refreshStripChoices();
}

void SignalPickerDock::onLiveToggled(bool live) {
    ASSERT_THREAD(this);
    if (live)
        refreshCandidateAnchor();
}

void SignalPickerDock::onRadiusChanged(double radius) {
    ASSERT_THREAD(this);
    nearbyModel_->setRadiusAngstrom(radius);
    refreshStripChoices();
}

void SignalPickerDock::onCandidateChanged() {
    ASSERT_THREAD(this);
    refreshStripChoices();
}

void SignalPickerDock::onAddClicked() {
    ASSERT_THREAD(this);
    // Intentionally empty in this increment. The next step is a dashboard strip
    // model that consumes the candidate + chosen strip entry as a pinned binding.
}

void SignalPickerDock::refreshFocusLabel() {
    if (!focusLabel_)
        return;
    if (!latestFocusAtom_) {
        focusLabel_->setText(QStringLiteral("Focus: none"));
        return;
    }
    focusLabel_->setText(QStringLiteral("Focus: %1").arg(atomDisplayLabel(*latestFocusAtom_)));
}

void SignalPickerDock::refreshCandidateAnchor() {
    if (!latestFocusAtom_) {
        nearbyModel_->clear();
        return;
    }
    nearbyModel_->setAnchor(*latestFocusAtom_, frame_);
}

void SignalPickerDock::refreshStripChoices() {
    if (!stripList_)
        return;
    stripList_->clear();

    const QModelIndex current = candidatesView_ ? candidatesView_->currentIndex() : QModelIndex();
    const NearbySignalModel::Candidate* candidate = nearbyModel_->candidateAt(current);
    if (!candidate)
        return;

    const model::SignalAnchorKind anchorKind = candidate->kind == NearbySignalModel::CandidateKind::Atom
                                                   ? model::SignalAnchorKind::Atom
                                                   : model::SignalAnchorKind::Residue;
    for (const model::SignalSpec& spec : model::CoreSignalDictionary()) {
        if (spec.anchorKind != anchorKind)
            continue;
        if (spec.key.family == model::SignalFamily::DftShielding && !dftStore_)
            continue;
        stripList_->addItem(QStringLiteral("Time  %1").arg(spec.label));
        if (isNumericSignal(spec.valueKind))
            stripList_->addItem(QStringLiteral("FFT   %1").arg(spec.label));
    }
}

QString SignalPickerDock::atomDisplayLabel(std::size_t atomIdx) const {
    if (!protein_ || atomIdx >= protein_->atomCount())
        return QStringLiteral("#%1").arg(atomIdx);
    const auto& atom = protein_->atom(atomIdx);
    const QString atomName = protein_->atomNames(atomIdx).iupac;
    if (atom.residueIndex >= 0
        && static_cast<std::size_t>(atom.residueIndex) < protein_->residueCount()) {
        const auto& residue = protein_->residue(static_cast<std::size_t>(atom.residueIndex));
        const QString chain = residue.address.chainId.isEmpty()
                                  ? QString()
                                  : QStringLiteral("%1:").arg(residue.address.chainId);
        const QString resName = protein_->residueLabel(static_cast<std::size_t>(atom.residueIndex),
                                                       model::NamingConvention::Iupac,
                                                       model::NamingSource::Verbatim);
        return QStringLiteral("%1%2%3:%4").arg(chain, resName).arg(residue.address.residueNumber).arg(atomName);
    }
    return QStringLiteral("#%1:%2").arg(atomIdx).arg(atomName);
}

}  // namespace h5reader::app

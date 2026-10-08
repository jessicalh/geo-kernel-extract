#include "SelectionContextWidget.h"

#include "CameraComposer.h"
#include "MoleculeScene.h"
#include "QtPlaybackController.h"
#include "SceneRevealOverlay.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/AtomSelection.h"
#include "../model/ConformationGeometry.h"
#include "../model/QtProtein.h"
#include "../model/QtResidueNames.h"
#include "../model/TrajectorySignalCatalog.h"
#include "../model/TransformedConformation.h"

#include <QHBoxLayout>
#include <QJsonArray>
#include <QLabel>
#include <QLoggingCategory>
#include <QStyle>
#include <QToolButton>
#include <QVBoxLayout>

namespace h5reader::app {
namespace {
Q_LOGGING_CATEGORY(cContext, "h5reader.selection_context")

QLabel* addLabel(QBoxLayout* layout, const char* name, bool bold = false) {
    auto* label = new QLabel;
    label->setObjectName(QString::fromLatin1(name));
    label->setTextFormat(Qt::PlainText);
    label->setWordWrap(true);
    label->setTextInteractionFlags(Qt::TextSelectableByMouse);
    auto font = label->font();
    font.setBold(bold);
    label->setFont(font);
    layout->addWidget(label);
    return label;
}
}

SelectionContextWidget::SelectionContextWidget(
    const model::QtProtein& protein, model::AtomSelection& selection,
    model::TransformedConformation& conformation, MoleculeScene& scene,
    QtPlaybackController& playback, const model::TrajectorySignalCatalog& catalog,
    QWidget* parent)
    : QWidget(parent), protein_(protein), selection_(selection),
      conformation_(conformation), scene_(scene), playback_(playback), catalog_(catalog) {
    setObjectName(QStringLiteral("SelectionContextWidget"));
    auto* layout = new QVBoxLayout(this);
    layout->setContentsMargins(6, 6, 6, 6);
    layout->setSpacing(6);
    auto* heading = new QHBoxLayout;
    kind_ = addLabel(heading, "contextKind", true);
    heading->addStretch();
    value_ = addLabel(heading, "contextValue", true);
    auto* clear = new QToolButton(this);
    clear->setObjectName(QStringLiteral("contextClearSelection"));
    clear->setIcon(style()->standardIcon(QStyle::SP_DialogCloseButton));
    clear->setAutoRaise(true);
    clear->setToolTip(tr("Clear selection"));
    clear->setAccessibleName(tr("Clear selection"));
    heading->addWidget(clear);
    clearSelection_ = clear;
    layout->addLayout(heading);
    atoms_ = addLabel(layout, "contextAtoms");
    meaning_ = addLabel(layout, "contextMeaning");
    const auto addDismissibleRow = [this, layout](QLabel*& label, const char* name,
                                                 const char* buttonName, const QString& action) {
        auto* row = new QHBoxLayout;
        label = addLabel(row, name);
        row->addStretch();
        auto* button = new QToolButton(this);
        button->setObjectName(QString::fromLatin1(buttonName));
        button->setIcon(style()->standardIcon(QStyle::SP_DialogCloseButton));
        button->setAutoRaise(true);
        button->setToolTip(action);
        button->setAccessibleName(action);
        row->addWidget(button, 0, Qt::AlignTop);
        layout->addLayout(row);
        return button;
    };
    clearHighlight_ = addDismissibleRow(highlight_, "contextHighlight", "contextClearHighlight", tr("Clear highlight"));
    releaseCamera_ = addDismissibleRow(camera_, "contextCamera", "contextReleaseCamera", tr("Release camera"));
    connect(clearSelection_, &QAbstractButton::clicked, &selection_, &model::AtomSelection::clear);
    connect(clearHighlight_, &QAbstractButton::clicked, &scene_, &MoleculeScene::clearReveal);
    connect(releaseCamera_, &QAbstractButton::clicked, this, [this] {
        scene_.cameraComposer()->setMode(FreeMode(), FreePolicy(), playback_.currentFrame());
    });
    connect(&selection_, &model::AtomSelection::changed, this, &SelectionContextWidget::refresh);
    connect(&scene_, &MoleculeScene::revealChanged, this, &SelectionContextWidget::refresh);
    connect(scene_.cameraComposer(), &CameraComposer::modeChanged, this, &SelectionContextWidget::refresh);
    connect(&conformation_, &model::TransformedConformation::transformChanged, this, &SelectionContextWidget::refresh);
    connect(&playback_, &QtPlaybackController::frameChanged, this, &SelectionContextWidget::refresh);
    refresh();
}

QString SelectionContextWidget::atomLabel(std::size_t index) const {
    const auto& atom = protein_.atom(index);
    const QString name = protein_.atomNames(index).iupac;
    if (atom.residueIndex < 0)
        return name;
    const auto& residue = protein_.residue(static_cast<std::size_t>(atom.residueIndex));
    const QString chain = residue.address.chainId.isEmpty()
        ? QString() : residue.address.chainId + QLatin1Char(':');
    return QStringLiteral("%1%2%3%4:%5")
        .arg(chain, QString::fromLatin1(model::IupacResidue3LetterFor(residue.aminoAcid)))
        .arg(residue.address.residueNumber).arg(residue.address.insertionCode, name);
}

QString SelectionContextWidget::atomList(const std::vector<std::size_t>& atoms) const {
    QStringList labels;
    for (const auto atom : atoms)
        labels.append(atomLabel(atom));
    return labels.join(QStringLiteral(", "));
}

void SelectionContextWidget::refresh() {
    ASSERT_THREAD(this);
    const auto& atoms = selection_.atoms();
    const int frame = playback_.currentFrame();
    QStringList labels;
    for (std::size_t i = 0; i < atoms.size(); ++i)
        labels.append(QStringLiteral("%1. %2").arg(i + 1).arg(atomLabel(atoms[i])));
    atoms_->setText(labels.join(QLatin1Char('\n')));
    clearSelection_->setVisible(!atoms.empty());
    value_->clear();
    kind_->setText(atoms.empty() ? tr("No atoms selected") : atomLabel(atoms.front()));
    meaning_->clear();
    if (atoms.size() >= 2) {
        const auto measurement = model::Measure(conformation_, frame, atoms);
        kind_->setText(QString::fromLatin1(model::NameForGeometryKind(measurement.kind)));
        if (!measurement.valid)
            value_->setText(tr("Undefined for these positions"));
        else if (measurement.kind == model::GeometryKind::Distance)
            value_->setText(QStringLiteral("%1 %2").arg(measurement.value, 0, 'f', 3).arg(QChar(0x00C5)));
        else
            value_->setText(QStringLiteral("%1%2").arg(measurement.value, 0, 'f', 1).arg(QChar(0x00B0)));
        if (atoms.size() == 2) {
            bool bonded = false;
            for (const auto index : protein_.topology().bondIndicesForAtom(atoms[0])) {
                const auto& bond = protein_.bond(index);
                bonded |= bond.atomIndexA == static_cast<std::int32_t>(atoms[1])
                       || bond.atomIndexB == static_cast<std::int32_t>(atoms[1]);
            }
            meaning_->setText(bonded ? tr("Bond length.") : tr("Atoms are not bonded."));
        } else if (atoms.size() == 3) {
            meaning_->setText(tr("Vertex at atom 2."));
        } else {
            meaning_->setText(tr("Axis through atoms 2 and 3."));
        }
    }
    value_->setVisible(!value_->text().isEmpty());
    atoms_->setVisible(atoms.size() >= 2);
    meaning_->setVisible(!meaning_->text().isEmpty());

    const auto& binding = scene_.activeRevealBinding();
    highlight_->setVisible(binding.has_value());
    clearHighlight_->setVisible(binding.has_value());
    if (binding) {
        const auto* descriptor = catalog_.findDescriptor(binding->descriptorId);
        const auto& highlighted = scene_.revealOverlay()->activeAtoms();
        const QString target = highlighted.size() <= model::AtomSelection::kMaxAtoms
            ? atomList(highlighted) : tr("%1 atoms").arg(highlighted.size());
        highlight_->setText(tr("Highlight: %1\n%2")
            .arg(descriptor ? descriptor->label : model::AnchorLabel(binding->anchor), target));
    } else {
        highlight_->clear();
    }
    const auto mode = scene_.cameraComposer()->mode();
    camera_->setVisible(mode.kind != CameraMode::Kind::Free);
    releaseCamera_->setVisible(mode.kind != CameraMode::Kind::Free);
    switch (mode.kind) {
    case CameraMode::Kind::Free: camera_->setText(tr("Camera: free")); break;
    case CameraMode::Kind::Atom: camera_->setText(tr("Camera follows %1").arg(atomList(mode.atoms))); break;
    case CameraMode::Kind::Bond: camera_->setText(tr("Camera follows atom pair: %1").arg(atomList(mode.atoms))); break;
    case CameraMode::Kind::Plane: camera_->setText(tr("Camera follows plane: %1").arg(atomList(mode.atoms))); break;
    case CameraMode::Kind::Dihedral: camera_->setText(tr("Camera follows dihedral: %1").arg(atomList(mode.atoms))); break;
    case CameraMode::Kind::Subset: camera_->setText(tr("Camera follows %1 atoms").arg(mode.atoms.size())); break;
    }
}

void SelectionContextWidget::setInteractionEnabled(bool enabled) {
    interactionEnabled_ = enabled;
    clearSelection_->setEnabled(enabled);
    clearHighlight_->setEnabled(enabled);
    releaseCamera_->setEnabled(enabled);
}

bool SelectionContextWidget::triggerAction(const QString& action) {
    ASSERT_THREAD(this);
    QAbstractButton* button = nullptr;
    if (action == QStringLiteral("clear_selection")) button = clearSelection_;
    else if (action == QStringLiteral("clear_highlight")) button = clearHighlight_;
    else if (action == QStringLiteral("release_camera")) button = releaseCamera_;
    if (!button || !button->isEnabled() || !button->isVisible())
        return false;
    qCInfo(cContext).noquote() << action;
    button->click();
    return true;
}

QJsonObject SelectionContextWidget::measurementState() const {
    QJsonArray atoms;
    for (const auto atom : selection_.atoms())
        atoms.append(static_cast<qint64>(atom));
    return {{"frame", playback_.currentFrame()}, {"atoms", atoms},
            {"kind", kind_->text()}, {"value", value_->text()}, {"labels", atoms_->text()},
            {"visible", isVisible()}, {"foreground", !visibleRegion().isEmpty()}};
}

QJsonObject SelectionContextWidget::stateJson() const {
    return {{"visible", isVisible()}, {"window", isWindow()},
            {"editable", interactionEnabled_},
            {"measurement", measurementState()}, {"meaning", meaning_->text()},
            {"highlight", highlight_->text()}, {"camera", camera_->text()}};
}
}  // namespace h5reader::app

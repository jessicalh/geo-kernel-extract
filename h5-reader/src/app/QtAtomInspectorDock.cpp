#include "QtAtomInspectorDock.h"
#include "MacWidgetStyle.h"
#include "AtomInspectorGlossary.h"
#include "MetricGlossaryPopup.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/QtConformationSnapshot.h"
#include "../model/QtResidueNames.h"
#include "../model/CsaShape.h"
#include "../model/TrajectoryConformation.h"

// Typed per-frame group views over the snapshot, the inspector's source for
// calculator detail.
#include "../model/QtAimnet2Group.h"
#include "../model/QtApbsGroup.h"
#include "../model/QtBiotSavartGroup.h"
#include "../model/QtBondedGroup.h"
#include "../model/QtCoulombGroup.h"
#include "../model/QtDsspGroup.h"
#include "../model/QtEeqGroup.h"
#include "../model/QtGromacsGroup.h"
#include "../model/QtHaighMallionGroup.h"
#include "../model/QtHBondGroup.h"
#include "../model/QtHydrationGroup.h"
#include "../model/QtLarsenHBondGroup.h"
#include "../model/QtMcConnellGroup.h"
#include "../model/QtMopacCoreGroup.h"
#include "../model/QtMopacCoulombGroup.h"
#include "../model/QtMopacMcConnellGroup.h"
#include "../model/QtOrcaGroup.h"
#include "../model/QtPlanarGeometryGroup.h"
#include "../model/QtSasaGroup.h"
#include "../model/QtWaterFieldGroup.h"
#include "../model/QtWaterPolarizationGroup.h"
#include "../physics/ClassicalSourceMath.h"
#include "../physics/LiteratureAccessors.h"
#include "../physics/RingCurrentScalars.h"
#include "../physics/SphericalBasis.h"
#include "constants/LiteratureConstants.h"

#include <QBrush>
#include <QColor>
#include "TensorGlyphPalette.h"
#include <QHeaderView>
#include <QIcon>
#include <QJsonArray>
#include <QJsonObject>
#include <QLoggingCategory>
#include <QPainter>
#include <QPixmap>
#include <QSizePolicy>
#include <QScrollBar>
#include <QSignalBlocker>
#include <QSplitter>
#include <QString>
#include <QStringList>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <utility>

namespace h5reader::app {

using model::QtEfg;
using model::SphericalTensor;
using model::Vec3;

namespace {
Q_LOGGING_CATEGORY(cDock, "h5reader.inspector")

// Formatting helpers — two-column tree with Field / Value. Keep the
// value text short enough to read at a glance; child nodes carry detail.

QString FmtDouble(double v, int precision = 4) {
    if (!std::isfinite(v))
        return QStringLiteral("nan");
    return QStringLiteral("%1").arg(v, 0, 'g', precision);
}

QString FmtQuantity(double v, const QString& unit = QString(), int precision = 5) {
    const QString text = FmtDouble(v, precision);
    return unit.isEmpty() ? text : text + QStringLiteral(" ") + unit;
}

QString FmtVec3(const Vec3& v, const QString& unit = QString(), int precision = 4) {
    return QStringLiteral("(%1, %2, %3)%4")
        .arg(FmtDouble(v.x(), precision),
             FmtDouble(v.y(), precision),
             FmtDouble(v.z(), precision),
             unit.isEmpty() ? QString() : QStringLiteral(" ") + unit);
}

double T1Magnitude(const SphericalTensor& st) {
    double s = 0.0;
    for (double c : st.T1)
        s += c * c;
    return std::sqrt(s);
}

model::CsaShape ShapeFromSphericalTensor(const SphericalTensor& st) {
    const model::Mat3 sym = h5reader::physics::ReconstructLibraryT2Matrix(st.T0, st.T2);
    return model::ComputeCsaShape(sym);
}

model::CsaShape ShapeFromEfg(const QtEfg& efg) {
    const model::Mat3 sym = h5reader::physics::ReconstructLibraryT2Matrix(efg.t2);
    return model::ComputeCsaShape(sym);
}

model::QtResidue ResidueOrEmpty(const model::QtProtein* protein, int residueIndex) {
    if (!protein || residueIndex < 0)
        return {};
    const auto idx = static_cast<std::size_t>(residueIndex);
    if (idx >= protein->residueCount())
        return {};
    return protein->residue(idx);
}

QString FmtSphericalSummary(const SphericalTensor& st, const QString& unit = QString()) {
    const model::CsaShape shape = ShapeFromSphericalTensor(st);
    QStringList parts;
    parts << QStringLiteral("T0=%1").arg(FmtQuantity(st.T0, unit));
    if (shape.valid) {
        parts << QStringLiteral("span=%1").arg(FmtQuantity(shape.span, unit));
        parts << QStringLiteral("eta=%1").arg(FmtDouble(shape.eta, 4));
    }
    parts << QStringLiteral("|T2|=%1").arg(FmtQuantity(st.T2Magnitude(), unit));
    return parts.join(QStringLiteral("  "));
}

QString FmtEfgSummary(const QtEfg& efg, const QString& unit = QString()) {
    const model::CsaShape shape = ShapeFromEfg(efg);
    QStringList parts;
    if (shape.valid) {
        parts << QStringLiteral("Vzz=%1").arg(FmtQuantity(shape.haeberlen_values[2], unit));
        parts << QStringLiteral("eta=%1").arg(FmtDouble(shape.eta, 4));
    }
    parts << QStringLiteral("|T2|=%1").arg(FmtQuantity(efg.t2Magnitude(), unit));
    return parts.join(QStringLiteral("  "));
}

QTreeWidgetItem* AddKV(QTreeWidgetItem* parent, const QString& field, const QString& value) {
    auto* it = new QTreeWidgetItem(parent);
    it->setText(0, field);
    it->setText(1, value);
    return it;
}

QIcon SwatchIcon(const QColor& color) {
    QPixmap pix(14, 14);
    pix.fill(Qt::transparent);

    QPainter painter(&pix);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen(QPen(color.darker(130), 1.0));
    painter.setBrush(color);
    painter.drawRect(3, 3, 8, 8);
    return QIcon(pix);
}

QTreeWidgetItem* AddSwatchKV(QTreeWidgetItem* parent,
                             const QString& field,
                             const QString& value,
                             const QColor& color) {
    auto* it = AddKV(parent, field, value);
    it->setIcon(0, SwatchIcon(color));
    return it;
}

bool AddScalar(QTreeWidgetItem* parent, const QString& name, double value, const QString& unit = QString()) {
    AddKV(parent, name, unit.isEmpty() ? FmtDouble(value) : FmtDouble(value) + QStringLiteral(" ") + unit);
    return true;
}

// Like AddScalar but returns the row and attaches a provenance/status tooltip
// (the curated metric-inventory label) to both columns, so a primary value
// carries its source + USE/PLACEHOLDER status on hover.
QTreeWidgetItem* AddScalarP(QTreeWidgetItem* parent, const QString& name, double value,
                            const QString& unit, const QString& provenance) {
    auto* it = AddKV(parent, name,
                     unit.isEmpty() ? FmtDouble(value) : FmtDouble(value) + QStringLiteral(" ") + unit);
    it->setToolTip(0, provenance);
    it->setToolTip(1, provenance);
    return it;
}

bool AddVec3(QTreeWidgetItem* parent, const QString& name, const Vec3& v, const QString& unit = QString()) {
    AddKV(parent, name, FmtVec3(v, unit));
    return true;
}

void AddTensorPrincipalRows(QTreeWidgetItem* parent,
                            const model::CsaShape& shape,
                            const QString& valuePrefix,
                            const QString& unit) {
    if (!shape.valid)
        return;

    static constexpr const char* kAxes[] = {"11", "22", "33"};
    const double vals[3] = {shape.principal_values[0], shape.principal_values[1], shape.principal_values[2]};
    for (int i = 0; i < 3; ++i) {
        const auto& rgb = kDefaultTensorColours[i];
        const QColor color = QColor::fromRgbF(rgb[0], rgb[1], rgb[2]);
        AddSwatchKV(parent,
                    QStringLiteral("%1_%2").arg(valuePrefix, QString::fromLatin1(kAxes[i])),
                    FmtQuantity(vals[i], unit),
                    color);
    }
}

void AddSphericalTensorTree(QTreeWidgetItem* it, const SphericalTensor& st, const QString& unit) {
    const model::CsaShape shape = ShapeFromSphericalTensor(st);

    AddScalar(it, QStringLiteral("T0 signed iso"), st.T0, unit);
    AddScalar(it, QStringLiteral("|T2| anisotropy"), st.T2Magnitude(), unit);
    const double t1Mag = T1Magnitude(st);
    if (std::abs(t1Mag) > 1e-12)
        AddScalar(it, QStringLiteral("|T1| antisymmetric"), t1Mag, unit);

    auto* pas = AddKV(it, QStringLiteral("PAS / shape"),
                      shape.valid
                          ? QStringLiteral("span=%1  eta=%2")
                                .arg(FmtQuantity(shape.span, unit), FmtDouble(shape.eta, 4))
                          : QStringLiteral("near-isotropic or unavailable"));
    if (shape.valid) {
        AddScalar(pas, QStringLiteral("span"), shape.span, unit);
        AddScalar(pas, QStringLiteral("eta"), shape.eta);
        AddScalar(pas, QStringLiteral("skew"), shape.skew);
        AddTensorPrincipalRows(pas, shape, QStringLiteral("sigma"), unit);
    }

    auto* raw = AddKV(it, QStringLiteral("raw irreps"), QStringLiteral("[T0, T1, T2]"));
    AddScalar(raw, QStringLiteral("T0"), st.T0, unit);
    AddVec3(raw, QStringLiteral("T1 antisym vector"), Vec3(st.T1[0], st.T1[1], st.T1[2]), unit);
    auto* t2 = AddKV(raw, QStringLiteral("T2 components (library basis)"), QString());
    for (int i = 0; i < 5; ++i)
        AddScalar(t2, QStringLiteral("m=%1").arg(i - 2), st.T2[i], unit);
}

void AddEfgTensorTree(QTreeWidgetItem* it, const QtEfg& efg, const QString& unit) {
    const model::CsaShape shape = ShapeFromEfg(efg);

    AddScalar(it, QStringLiteral("|T2| invariant"), efg.t2Magnitude(), unit);
    auto* pas = AddKV(it, QStringLiteral("PAS / EFG convention"),
                      shape.valid
                          ? QStringLiteral("Vzz=%1  eta=%2")
                                .arg(FmtQuantity(shape.haeberlen_values[2], unit),
                                     FmtDouble(shape.eta, 4))
                          : QStringLiteral("near-isotropic or unavailable"));
    if (shape.valid) {
        AddScalar(pas, QStringLiteral("Vxx"), shape.haeberlen_values[0], unit);
        AddScalar(pas, QStringLiteral("Vyy"), shape.haeberlen_values[1], unit);
        AddScalar(pas, QStringLiteral("Vzz"), shape.haeberlen_values[2], unit);
        AddScalar(pas, QStringLiteral("eta"), shape.eta);
        AddScalar(pas, QStringLiteral("span"), shape.span, unit);
        AddTensorPrincipalRows(pas, shape, QStringLiteral("V"), unit);
    }

    auto* raw = AddKV(it, QStringLiteral("raw T2 components (library basis)"), QString());
    for (int i = 0; i < 5; ++i)
        AddScalar(raw, QStringLiteral("m=%1").arg(i - 2), efg.t2[i], unit);
}

[[maybe_unused]] bool AddSpherical(QTreeWidgetItem* parent, const QString& name, const SphericalTensor& st, const QString& unit = QString()) {
    auto* it = AddKV(parent, name, FmtSphericalSummary(st, unit));
    it->setToolTip(0, QStringLiteral("Tensor tree: signed T0, invariant magnitudes, PAS shape, and raw irreps."));
    it->setToolTip(1, it->toolTip(0));
    AddSphericalTensorTree(it, st, unit);
    return true;
}

void DeleteIfEmpty(QTreeWidgetItem* item) {
    if (item && item->childCount() == 0)
        delete item;
}

// Optional-aware adders — add only real rows. Missing values mean "this
// calculator did not run / this field is absent"; dash-only groups are cut.
bool AddOptScalar(QTreeWidgetItem* p, const QString& name, const std::optional<double>& v, const QString& unit = QString()) {
    if (!v || !std::isfinite(*v))
        return false;
    return AddScalar(p, name, *v, unit);
}
bool AddOptVec3(QTreeWidgetItem* p, const QString& name, const std::optional<Vec3>& v, const QString& unit = QString()) {
    if (!v) return false;
    return AddVec3(p, name, *v, unit);
}
bool AddOptSpherical(QTreeWidgetItem* p, const QString& name, const std::optional<SphericalTensor>& v, const QString& unit = QString()) {
    if (!v) return false;
    return AddSpherical(p, name, *v, unit);
}
bool AddOptEfg(QTreeWidgetItem* p, const QString& name, const std::optional<QtEfg>& v, const QString& unit = QString()) {
    if (!v) return false;
    auto* it = AddKV(p, name, FmtEfgSummary(*v, unit));
    it->setToolTip(0, QStringLiteral("EFG tensor tree: T2 invariant, PAS Vzz/eta, and raw T2 components."));
    it->setToolTip(1, it->toolTip(0));
    AddEfgTensorTree(it, *v, unit);
    return true;
}
bool AddOptInt(QTreeWidgetItem* p, const QString& name, const std::optional<int>& v) {
    if (!v) return false;
    AddKV(p, name, QString::number(*v));
    return true;
}
bool AddOptBool(QTreeWidgetItem* p, const QString& name, const std::optional<bool>& v) {
    if (!v) return false;
    AddKV(p, name, *v ? QStringLiteral("true") : QStringLiteral("false"));
    return true;
}

bool AllowsAny(const std::shared_ptr<const model::TrajectoryFieldAvailability>& availability,
               std::initializer_list<const char*> descriptorIds) {
    if (!availability)
        return true;
    for (const char* id : descriptorIds) {
        const auto state = availability->stateForDescriptor(QString::fromLatin1(id));
        if (model::TrajectoryFieldAvailability::isVisibleState(state))
            return true;
    }
    return false;
}

}  // namespace

QtAtomInspectorDock::QtAtomInspectorDock(QWidget* parent) : QDockWidget(QStringLiteral("Atom Info"), parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtAtomInspectorDock"));
    setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    setMinimumWidth(260);
    setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Expanding);

    auto* splitter = new QSplitter(Qt::Vertical, this);
    splitter->setObjectName(QStringLiteral("atomInfoSplitter"));
    splitter->setChildrenCollapsible(false);
    tree_ = new QTreeWidget(splitter);
    tree_->setObjectName(QStringLiteral("atomFields"));
    auto* lower = new QWidget(splitter);
    tensorLayout_ = new QVBoxLayout(lower);
    tensorLayout_->setContentsMargins(0, 0, 0, 0);
    tensorLayout_->setSpacing(0);
    tensorTree_ = new QTreeWidget(lower);
    tensorTree_->setObjectName(QStringLiteral("tensorFields"));
#ifdef Q_OS_MACOS
    configureMacItemViewChecks(tensorTree_);
#endif
    tensorLayout_->addWidget(tensorTree_);
    for (auto* tree : {tree_.data(), tensorTree_.data()}) {
        tree->setMinimumWidth(0);
        tree->setMinimumHeight(100);
        tree->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Expanding);
        tree->setColumnCount(2);
        tree->setHeaderLabels({QStringLiteral("Field"), QStringLiteral("Value")});
        tree->setAlternatingRowColors(true);
        tree->header()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
        tree->setContextMenuPolicy(Qt::CustomContextMenu);
        connect(tree, &QWidget::customContextMenuRequested, this, [this, tree](const QPoint& position) {
            const auto* item = tree->itemAt(position);
            if (!item) return;
            QStringList path;
            for (const auto* ancestor = item; ancestor; ancestor = ancestor->parent())
                path.prepend(ancestor->text(0));
            const QString source = csa_.sourceDetail.isEmpty()
                ? csa_.sourceLabel : csa_.sourceLabel + QStringLiteral(": ") + csa_.sourceDetail;
            if (const auto help = AtomInspectorGlossary(path, source))
                ShowMetricGlossaryPopup(item->text(0), *help, tree->viewport()->mapToGlobal(position), this);
        });
    }
    connect(tensorTree_, &QTreeWidget::itemChanged, this,
            [this](QTreeWidgetItem* item, int column) {
        if (column != 0) return;
        if (item == csaGroup_)
            setTensorDisplayEnabled(item->checkState(0) == Qt::Checked, orientationEnabled_);
        else if (item == orientationGroup_)
            setTensorDisplayEnabled(shieldingEnabled_, item->checkState(0) == Qt::Checked);
    });
    splitter->setStretchFactor(0, 1);
    splitter->setStretchFactor(1, 1);
    splitter->setSizes({300, 400});
    setWidget(splitter);

    // Starting placeholder.
    auto* hint = new QTreeWidgetItem(tree_);
    hint->setText(0, QStringLiteral("Double-click an atom in the viewport"));
}

void QtAtomInspectorDock::setSelectionContextWidget(QWidget* context) {
    tensorLayout_->insertWidget(0, context);
}

void QtAtomInspectorDock::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    if (conformation_)
        disconnect(conformation_.data(), nullptr, this, nullptr);
    protein_ = protein;
    conformation_ = conformation;
    frame_ = 0;
    if (conformation_) {
        QObject::connect(conformation_.data(), &model::Conformation::snapshotReady,
                 this, &QtAtomInspectorDock::onSnapshotReady);
    } else {
        clearSelection();
    }
}

void QtAtomInspectorDock::setFieldAvailability(
    std::shared_ptr<const model::TrajectoryFieldAvailability> availability) {
    availability_ = std::move(availability);
    rebuild();
}

void QtAtomInspectorDock::setPickedAtom(std::size_t atomIdx) {
    ASSERT_THREAD(this);
    if (!hasSelection_ || atomIdx_ != atomIdx) {
        tree_->verticalScrollBar()->setValue(0);
        tensorTree_->verticalScrollBar()->setValue(0);
    }
    hasSelection_ = true;
    atomIdx_ = atomIdx;
    requestCurrentSnapshot();
}

void QtAtomInspectorDock::requestCurrentSnapshot() {
    ASSERT_THREAD(this);
    if (!hasSelection_)
        return;
    const std::size_t frame = static_cast<std::size_t>(std::max(0, frame_));
    if (conformation_) {
        conformation_->requestSnapshotAsync(frame);
        if (conformation_->snapshot(frame))
            return;  // snapshotReady already rebuilt the tree.
    }
    rebuild();
}

void QtAtomInspectorDock::setFrame(int t) {
    ASSERT_THREAD(this);
    if (t != frame_)
        hasCsa_ = false;  // CSA values are frame-local; never carry old tensors forward.
    frame_ = t;
    if (hasSelection_)
        rebuild();
    else
        rebuildTensors();
}

void QtAtomInspectorDock::onSnapshotReady(std::size_t frame) {
    ASSERT_THREAD(this);
    if (hasSelection_ && static_cast<int>(frame) == frame_)
        rebuild();
}

void QtAtomInspectorDock::clearSelection() {
    ASSERT_THREAD(this);
    hasSelection_ = false;
    hasCsa_ = false;
    hasOrient_ = false;
    csaAtom_ = 0;
    orientAtom_ = 0;
    csaGroup_ = nullptr;
    orientationGroup_ = nullptr;
    shieldingVisible_ = false;
    orientationVisible_ = false;
    tree_->clear();
    tensorTree_->clear();
    auto* hint = new QTreeWidgetItem(tree_);
    hint->setText(0, QStringLiteral("Double-click an atom in the viewport"));
}

void QtAtomInspectorDock::setCsaTensor(std::size_t atom, const CsaTensorInfo& info) {
    csa_ = info;
    csaAtom_ = atom;
    hasCsa_ = true;
    rebuildTensors();
}

void QtAtomInspectorDock::clearCsaTensor() {
    if (!hasCsa_)
        return;
    hasCsa_ = false;
    rebuildTensors();
}

void QtAtomInspectorDock::setOrientationTensor(std::size_t atom, const OrientationTensorInfo& info) {
    orient_ = info;
    orientAtom_ = atom;
    hasOrient_ = true;
    rebuildTensors();
}

void QtAtomInspectorDock::clearOrientationTensor() {
    if (!hasOrient_)
        return;
    hasOrient_ = false;
    rebuildTensors();
}

void QtAtomInspectorDock::setTensorVisibility(bool shielding, bool orientation) {
    ASSERT_THREAD(this);
    if (shieldingVisible_ == shielding && orientationVisible_ == orientation)
        return;
    shieldingVisible_ = shielding;
    orientationVisible_ = orientation;
    refreshTensorHeadings();
}

void QtAtomInspectorDock::refreshTensorHeadings() {
    const QSignalBlocker blocker(tensorTree_);
    const auto markShown = [](QTreeWidgetItem* item, bool enabled, bool shown) {
        if (!item) return;
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(0, enabled ? Qt::Checked : Qt::Unchecked);
        item->setText(1, shown ? QStringLiteral("Shown") : QString());
        auto font = item->font(0);
        font.setBold(shown);
        item->setFont(0, font);
        item->setFont(1, font);
    };
    markShown(csaGroup_, shieldingEnabled_, shieldingVisible_);
    if (csaGroup_ && !csa_.status.isEmpty())
        csaGroup_->setText(1, csa_.status);
    markShown(orientationGroup_, orientationEnabled_, orientationVisible_);
}

void QtAtomInspectorDock::setTensorDisplayEnabled(bool shielding, bool orientation) {
    ASSERT_THREAD(this);
    if (shieldingEnabled_ == shielding && orientationEnabled_ == orientation)
        return;
    const bool requestShielding = shielding && !shieldingEnabled_;
    shieldingEnabled_ = shielding;
    orientationEnabled_ = orientation;
    refreshTensorHeadings();
    qCInfo(cDock) << "tensor display | shielding=" << shielding
                 << "| orientation=" << orientation;
    emit tensorDisplayChanged(shielding, orientation);
    if (requestShielding)
        emit shieldingTensorRequested();
}

void QtAtomInspectorDock::populateCsa(QTreeWidgetItem* root) {
    const QString source = csa_.sourceLabel.isEmpty()
                               ? QStringLiteral("unknown source")
                               : csa_.sourceLabel;
    const QString frame = csa_.frameKind.isEmpty()
                              ? (csa_.framed ? QStringLiteral("molecular frame")
                                             : QStringLiteral("unframed"))
                              : csa_.frameKind;
    auto* group = AddKV(root, QStringLiteral("Shielding tensor (%1)").arg(source), QString());
    csaGroup_ = group;
    group->setExpanded(true);
    group->setToolTip(0, QStringLiteral("Symmetric shielding tensor for the named atom in the current frame."));
    const auto& atom = protein_->atom(csaAtom_);
    const auto residue = ResidueOrEmpty(protein_, atom.residueIndex);
    AddKV(group, QStringLiteral("Atom"),
          QStringLiteral("%1:%2%3%4:%5")
              .arg(residue.address.chainId, QString::fromLatin1(model::IupacResidue3LetterFor(residue.aminoAcid)))
              .arg(residue.address.residueNumber)
              .arg(residue.address.insertionCode, protein_->atomNames(csaAtom_).iupac));
    if (!csa_.status.isEmpty()) {
        AddKV(group, QStringLiteral("Source"), csa_.sourceDetail);
        return;
    }
    AddScalar(group, QStringLiteral("sigma_iso"), csa_.sigmaIso, QStringLiteral("ppm"));
    auto* axes = AddKV(group, QStringLiteral("Principal values"), QStringLiteral("ppm (from mean)"));
    axes->setToolTip(0, QStringLiteral("Shielding along each coloured axis. Parentheses give the difference from sigma_iso, in ppm."));
    axes->setToolTip(1, axes->toolTip(0));
    axes->setExpanded(true);

    // Per-axis principal values. The swatch is the colour key for the scene
    // arrows; no floating labels in the molecule view.
    static constexpr const char* kAxes[] = {"sigma_11", "sigma_22", "sigma_33"};
    const double vals[3] = {csa_.sigma11, csa_.sigma22, csa_.sigma33};
    for (int i = 0; i < 3; ++i) {
        const auto& rgb = kShieldingTensorColours[i];
        const QColor color = QColor::fromRgbF(rgb[0], rgb[1], rgb[2]);
        const double difference = vals[i] - csa_.sigmaIso;
        const QString signedDifference = (difference > 0.0 ? QStringLiteral("+") : QString())
                                         + FmtDouble(difference);
        auto* axis = AddSwatchKV(axes, QString::fromLatin1(kAxes[i]),
                                QStringLiteral("%1 (%2)").arg(FmtDouble(vals[i]), signedDifference), color);
        axis->setToolTip(0, axes->toolTip(0));
        axis->setToolTip(1, axes->toolTip(0));
    }
    AddScalar(group, QStringLiteral("span"), csa_.span, QStringLiteral("ppm"));
    auto* scale = AddKV(group, QStringLiteral("Glyph size"), QStringLiteral("Normalised"));
    scale->setToolTip(0, QStringLiteral("Each tensor is scaled separately. Arrow length follows the size of its deviation from the mean, with a minimum for visibility."));
    scale->setToolTip(1, scale->toolTip(0));
    auto* detail = AddKV(group, QStringLiteral("Details"), QStringLiteral("Source and shape"));
    AddKV(detail, QStringLiteral("Scope"), QStringLiteral("Current frame"));
    if (!csa_.sourceDetail.isEmpty())
        AddKV(detail, QStringLiteral("Source"), csa_.sourceDetail);
    AddKV(detail, QStringLiteral("Frame convention"), frame);
    AddScalar(detail, QStringLiteral("skew"), csa_.skew);
    AddScalar(detail, QStringLiteral("eta"), csa_.eta);
}

void QtAtomInspectorDock::populateOrientation(QTreeWidgetItem* root) {
    auto* group = AddKV(root, QStringLiteral("Bond orientation tensor"), QString());
    orientationGroup_ = group;
    group->setExpanded(true);
    group->setToolTip(0, QStringLiteral("Average of u u^T across aligned trajectory frames, where u is the unit bond vector."));
    AddKV(group, QStringLiteral("Bond"), orient_.bond);
    AddKV(group, QStringLiteral("Average over"), QStringLiteral("Trajectory"));
    auto* order = AddKV(group, QStringLiteral("S^2 (order parameter)"), FmtDouble(orient_.s2));
    order->setToolTip(0, QStringLiteral("Calculated from the three eigenvalues: (3 sum(lambda^2) - 1) / 2. A fixed bond direction gives 1; isotropic directions give 0."));
    order->setToolTip(1, order->toolTip(0));

    // Order-tensor eigenvalues (descending; sum to 1). The swatch is the
    // colour key for the scene arrows; no floating labels in the molecule view.
    static constexpr const char* kAxes[] = {"lambda_1", "lambda_2", "lambda_3"};
    const double vals[3] = {orient_.lambda1, orient_.lambda2, orient_.lambda3};
    for (int i = 0; i < 3; ++i) {
        const auto& rgb = kOrientationTensorColours[i];
        const QColor color = QColor::fromRgbF(rgb[0], rgb[1], rgb[2]);
        auto* axis = AddSwatchKV(group, QString::fromLatin1(kAxes[i]), FmtDouble(vals[i]), color);
        const QString meaning = QStringLiteral("Mean squared projection of the unit bond vector on this principal axis. The three values sum to 1; they are not fractions of frames.");
        axis->setToolTip(0, meaning);
        axis->setToolTip(1, meaning);
    }
    auto* placement = AddKV(group, QStringLiteral("Main axis follows"), QStringLiteral("Current bond"));
    placement->setToolTip(0, QStringLiteral("The largest-eigenvalue axis is rotated onto the current bond. Shape and eigenvalues describe the trajectory average. The arrows do not show the directions of the average tensor in the aligned protein."));
    placement->setToolTip(1, placement->toolTip(0));
    auto* detail = AddKV(group, QStringLiteral("Details"), QString());
    AddKV(detail, QStringLiteral("Source"), QStringLiteral("nmr_extract: reorientational dynamics"));
    AddKV(detail, QStringLiteral("Glyph size"), QStringLiteral("Normalised per tensor"));
}

void QtAtomInspectorDock::rebuild() {
    if (!tree_ || !protein_ || !conformation_)
        return;
    if (!hasSelection_ || atomIdx_ >= protein_->atomCount())
        return;

    // Batch the rebuild into a single repaint: this tree is cleared + fully
    // repopulated on every focus / frame change, so per-item updates would flicker.
    const int scroll = tree_->verticalScrollBar()->value();
    tree_->setUpdatesEnabled(false);
    tree_->clear();

    auto* title = new QTreeWidgetItem(tree_);
    const auto& atom = protein_->atom(atomIdx_);
    const auto res = ResidueOrEmpty(protein_, atom.residueIndex);
    title->setText(
        0,
        QStringLiteral("Atom %1 — %2 %3 #%4")
            .arg(atomIdx_)
            .arg(protein_->atomNames(atomIdx_).amber, QString::fromLatin1(model::IupacResidue3LetterFor(res.aminoAcid)))
            .arg(res.address.residueNumber));
    title->setText(1, QStringLiteral("frame %1 / %2").arg(frame_ + 1).arg(conformation_->frameCount()));
    title->setExpanded(true);

    populateIdentity(title);
    // Raw kernels / diagnostics collapse into ONE drawer at the very bottom: the
    // npy "show your work" stays available but does not compete with the
    // validated metrics. Built as an orphan, filled in populatePerFrame, attached
    // last iff it ended up with content.
    auto* drawer = new QTreeWidgetItem();
    drawer->setText(0, QStringLiteral("Raw kernels & diagnostics"));
    drawer->setText(1, QStringLiteral("raw npy inputs (not validated)"));
    populatePerFrame(title, drawer);
    for (int i = 0; i < title->childCount(); ++i)
        title->child(i)->setExpanded(true);
    if (drawer->childCount() > 0) {
        title->addChild(drawer);     // tree takes ownership; collapse AFTER attach
        drawer->setExpanded(false);
    } else {
        delete drawer;               // orphan, never attached -- we own it
    }

    tree_->verticalScrollBar()->setValue(scroll);
    tree_->setUpdatesEnabled(true);
    rebuildTensors();
}

void QtAtomInspectorDock::rebuildTensors() {
    const int scroll = tensorTree_->verticalScrollBar()->value();
    const QSignalBlocker blocker(tensorTree_);
    tensorTree_->setUpdatesEnabled(false);
    csaGroup_ = nullptr;
    orientationGroup_ = nullptr;
    tensorTree_->clear();
    if (protein_ && hasCsa_)
        populateCsa(tensorTree_->invisibleRootItem());
    if (protein_ && hasSelection_ && hasOrient_ && orientAtom_ == atomIdx_)
        populateOrientation(tensorTree_->invisibleRootItem());
    if (hasSelection_ && tensorTree_->topLevelItemCount() == 0)
        AddKV(tensorTree_->invisibleRootItem(), QStringLiteral("No tensor for this atom"), QString());
    refreshTensorHeadings();
    tensorTree_->verticalScrollBar()->setValue(scroll);
    tensorTree_->setUpdatesEnabled(true);
}

void QtAtomInspectorDock::populateIdentity(QTreeWidgetItem* parent) {
    const auto& atom = protein_->atom(atomIdx_);
    const auto res = ResidueOrEmpty(protein_, atom.residueIndex);

    auto* g = AddKV(parent, QStringLiteral("Identity"), QString());
    g->setExpanded(true);

    AddKV(g, QStringLiteral("Element"), QString::fromLatin1(model::SymbolForElement(atom.element)));
    AddKV(g, QStringLiteral("AMBER name"), protein_->atomNames(atomIdx_).amber);
    AddKV(g, QStringLiteral("IUPAC name"), protein_->atomNames(atomIdx_).iupac);
    AddKV(g, QStringLiteral("BMRB name"), protein_->atomNames(atomIdx_).bmrb);
    AddKV(g, QStringLiteral("Backbone role"), QString::number(static_cast<int>(atom.backboneRole)));
    AddKV(g, QStringLiteral("Locant"), QString::number(static_cast<int>(atom.locant)));
    AddKV(g,
          QStringLiteral("Residue"),
          QStringLiteral("%1 #%2").arg(QString::fromLatin1(model::IupacResidue3LetterFor(res.aminoAcid)),
                                       QString::number(res.address.residueNumber)));
    AddKV(g, QStringLiteral("Chain"), res.address.chainId.isEmpty() ? QStringLiteral("—") : res.address.chainId);
    AddKV(g,
          QStringLiteral("Protonation variant"),
          res.protonationVariantIndex < 0 ? QStringLiteral("default") : QString::number(res.protonationVariantIndex));
    AddScalar(g, QStringLiteral("Covalent radius"), atom.CovalentRadius(), QStringLiteral("Å"));
    AddScalar(g, QStringLiteral("Formal charge"), static_cast<double>(atom.formalCharge), QStringLiteral("e"));

    auto* flags = AddKV(g, QStringLiteral("Substrate flags"), QString());
    AddKV(flags, QStringLiteral("is_backbone"), atom.IsBackbone() ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags,
          QStringLiteral("is_amide_H"),
          (atom.polarH == model::PolarHKind::BackboneAmide) ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_alpha_H"), atom.IsAnyAlphaHydrogen() ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags,
          QStringLiteral("is_methyl"),
          (atom.pseudoatomKind == model::PseudoatomKind::M) ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_aromatic"), atom.aromatic ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_polar_H"), atom.IsPolarH() ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags,
          QStringLiteral("is_hbond_acceptor_elem"),
          atom.IsHBondAcceptorElement() ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_exchangeable"), atom.isExchangeable ? QStringLiteral("true") : QStringLiteral("false"));
}

void QtAtomInspectorDock::populatePerFrame(QTreeWidgetItem* root, QTreeWidgetItem* drawer) {
    const int T = static_cast<int>(conformation_->frameCount());
    const int t = std::clamp(frame_, 0, std::max(0, T - 1));
    const std::size_t st = static_cast<std::size_t>(t);

    const std::size_t a = atomIdx_;

    // Position — the shared seam: the H5 for a trajectory, the snapshot's
    // Pos column for a single pose.
    auto* posG = AddKV(root, QStringLiteral("Position"), QString());
    AddVec3(posG, QStringLiteral("xyz"), conformation_->atomPosition(st, a), QStringLiteral("Å"));

    // Per-frame calculator detail comes from typed views over the snapshot, not
    // from the H5 time series. nullopt means the calculator did not run for this
    // frame, and the panel displays an em dash.
    auto snap = conformation_->snapshot(st);
    if (!snap) {
        auto* g = AddKV(root, QStringLiteral("Per-frame detail"), QStringLiteral("not sampled at this frame"));
        AddKV(g, QStringLiteral("note"),
              QStringLiteral("the pick registered — full per-atom detail is emitted at a frame "
                             "stride, not every frame"));
        if (const auto* traj = conformation_->asTrajectory()) {
            if (auto nf = traj->nearestSampledFrame(st))
                AddKV(g, QStringLiteral("nearest sampled frame"),
                      QStringLiteral("%1  (scrub here for the full pile)").arg(*nf));
        }
        tree_->expandToDepth(1);
        qCDebug(cDock).noquote() << "rebuilt | atom=" << a << "| frame=" << t << "| snapshot= absent";
        return;
    }
    const auto& s = *snap;

    // --- Local classical estimate (VIEWER-DERIVED, TENTATIVE) --- a local
    // estimate computed HERE from loaded npy via the shared physics math
    // (ring, McConnell, Larsen, Buckingham -> ComputeClassicalSigma). Treat the
    // fold like a regression-style explanatory scaffold, not a validated
    // absolute shielding model. The reader NEVER runs the emit; the raw kernels
    // that feed these live in the drawer below. sigma_cl + residual appear only
    // when every mechanism's input is present (and sigma_qm/ORCA, single-pose,
    // for the residual).
    {
        model::QtBiotSavartGroup bsFwd(s);
        model::QtLarsenHBondGroup larsenFwd(s);
        QTreeWidgetItem* fwd = nullptr;
        auto ensureFwd = [&]() -> QTreeWidgetItem* {
            if (!fwd) {
                fwd = AddKV(root, QStringLiteral("Local classical estimate (tentative)"), QString());
                fwd->setExpanded(true);
            }
            return fwd;
        };

        // Per-atom identity for the literature-constant lookups (Buckingham, sigma0).
        const auto& fa = protein_->atom(atomIdx_);
        const auto fres = ResidueOrEmpty(protein_, fa.residueIndex);
        const std::string residueStd = model::IupacResidue3LetterFor(fres.aminoAcid);
        const std::string atomNameStd = protein_->atomNames(atomIdx_).amber.toStdString();
        const std::string frameKindStd =
            (hasCsa_ && csaAtom_ == atomIdx_ && csa_.framed) ? csa_.frameKind.toStdString()
                                                             : std::string();
        h5reader::physics::ClassicalSigmaInputs in;  // each term defaults to {NaN, present=false}

        // ring = sum_t bs_per_type_T0[t] * RingIntensity[t]  (signed T0, ppm)
        if (auto pt = bsFwd.perTypeT0(a)) {
            const double ring =
                h5reader::physics::RingPerTypeT0Ppm(pt->byType.data(), model::kAromaticRingTypeCount);
            if (std::isfinite(ring)) {
                in.ring = {ring, true};
                AddScalarP(ensureFwd(), QStringLiteral("estimated ring contribution"), ring, QStringLiteral("ppm"),
                           QStringLiteral("Biot-Savart per-type signed T0 x Giessner-Prettre ring "
                                          "intensity (ppm); local estimate term, not a measured value. "
                                          "[LOCAL ESTIMATE]"));
            }
        }
        // Larsen = signed T0 + ProCS15 water term  (ppm)
        if (auto sh = larsenFwd.shielding(a)) {
            double lars = sh->T0;
            if (auto wt = larsenFwd.waterTerm(a))
                lars += *wt;
            if (std::isfinite(lars)) {
                in.larsen = {lars, true};
                AddScalarP(ensureFwd(), QStringLiteral("estimated Larsen contribution"), lars, QStringLiteral("ppm"),
                           QStringLiteral("Larsen/ProCS15 signed H-bond T0 + 2.07 ppm ProCS15 water "
                                          "term (ppm); local estimate term, not a measured value. "
                                          "[LOCAL ESTIMATE]"));
            }
        }
        // McConnell = Sum over the 6 forward bond producers of
        // McConnellProducerT0ToPpm(category, producer.T0)  (signed, ppm). The raw
        // per-category _bo kernels feed this; the validated ppm form is this
        // viewer-derived sum. No disulfide producer in the forward set (mirrors
        // the engine's kMcForwardProducerFields).
        {
            static constexpr struct {
                io::FieldKind kind;
                model::BondCategory category;
            } kMcProducers[] = {
                {io::FieldKind::McPeptideCoBo, model::BondCategory::PeptideCO},
                {io::FieldKind::McPeptideCNBo, model::BondCategory::PeptideCN},
                {io::FieldKind::McBackboneOtherBo, model::BondCategory::BackboneOther},
                {io::FieldKind::McSidechainCoBo, model::BondCategory::SidechainCO},
                {io::FieldKind::McSidechainOtherBo, model::BondCategory::SidechainOther},
                {io::FieldKind::McAromaticBo, model::BondCategory::Aromatic},
            };
            model::QtMcConnellGroup mcFwd(s);
            double mc = 0.0;
            bool anyMc = false;
            for (const auto& p : kMcProducers) {
                auto bo = mcFwd.producerBo(p.kind, a);
                if (!bo) continue;
                const double term = h5reader::physics::McConnellProducerT0ToPpm(p.category, bo->T0);
                if (!std::isfinite(term)) continue;
                mc += term;
                anyMc = true;
            }
            if (anyMc && std::isfinite(mc)) {
                in.mcconnell = {mc, true};
                AddScalarP(ensureFwd(), QStringLiteral("estimated McConnell contribution"), mc, QStringLiteral("ppm"),
                           QStringLiteral("Wiberg-weighted (_bo) McConnell signed T0 x DeltaChi x "
                                          "molar prefactor (ppm); local estimate term, not a measured value. "
                                          "[LOCAL ESTIMATE]"));
            }
        }
        // Buckingham inputs: signed E|| = mopac_coulomb_scalars E_bond_proj (the
        // "Buckingham sigma_iso input" = component 1, MOPAC Coulomb field on the
        // parent-to-H bond axis) + element/frame-specific A,B literature constants. The
        // -A*E|| - B*E||^2 fold (and the row) come from ComputeClassicalSigma below.
        if (auto sc = model::QtMopacCoulombGroup(s).scalars(a)) {
            const double ePar = sc->E_bond_proj;
            if (std::isfinite(ePar)) {
                in.e_parallel_mopac = {ePar, true};
                const auto bA = h5reader::physics::BuckinghamA(fa.element, residueStd, atomNameStd, frameKindStd);
                const auto bB = h5reader::physics::BuckinghamB(fa.element, residueStd, atomNameStd, frameKindStd);
                if (std::isfinite(bA.value)) in.buckingham_A = {bA.value, true};
                if (std::isfinite(bB.value)) in.buckingham_B = {bB.value, true};
            }
        }
        // sigma0 is currently only a placeholder in the constants table. Do not
        // feed it into advisor-facing absolute estimates; keep the local
        // contribution terms visible on their own.
        const auto sigma0c = h5reader::physics::Sigma0(fa.element, residueStd, atomNameStd);
        if (sigma0c.status != nmr::constants::LiteratureStatus::Placeholder
            && std::isfinite(sigma0c.value))
            in.sigma0 = {sigma0c.value, true};

        // Fold all terms through the shared engine math (single source of truth);
        // the Buckingham term (-A*E|| - B*E||^2) is computed inside the fold.
        const h5reader::physics::ClassicalSigmaResult folded = h5reader::physics::ComputeClassicalSigma(in);
        if (folded.buckingham.present && std::isfinite(folded.buckingham.value))
            AddScalarP(ensureFwd(), QStringLiteral("estimated Buckingham contribution"),
                       folded.buckingham.value, QStringLiteral("ppm"),
                       QStringLiteral("-A*E_parallel - B*E_parallel^2 with signed MOPAC bond-axis "
                                      "field and literature A,B (ppm); local estimate term, not a measured "
                                      "value. [LOCAL ESTIMATE]"));

        // sigma_cl + residual ONLY when the full classical model is computable
        // (every mechanism's source field present). Otherwise sigma_cl would
        // silently drop a term and mislead -- show just the present contributions.
        const bool fwdComplete = in.ring.present && in.mcconnell.present
                                 && in.larsen.present && in.e_parallel_mopac.present;
        if (fwdComplete && folded.sigma_cl.present && std::isfinite(folded.sigma_cl.value)) {
            if (folded.sigma0.present) {
                AddScalarP(ensureFwd(), QStringLiteral("estimated sigma_cl (tentative)"),
                           folded.sigma_cl.value, QStringLiteral("ppm"),
                           QStringLiteral("Local estimated fold: sigma0 + Buckingham + ring + McConnell + Larsen. "
                                          "Use as a tentative local explanatory/regression-style model. [ESTIMATE]"));
                // residual = sigma_qm - sigma_cl; sigma_qm = ORCA total T0
                // (single-pose DFT only, so absent on a trajectory snapshot).
                if (auto qm = model::QtOrcaGroup(s).total(a)) {
                    if (std::isfinite(qm->T0))
                        AddScalarP(ensureFwd(), QStringLiteral("tentative residual (sigma_qm - estimate)"),
                                   qm->T0 - folded.sigma_cl.value, QStringLiteral("ppm"),
                                   QStringLiteral("sigma_qm (ORCA total signed T0) minus the local tentative "
                                                  "classical estimate (ppm). [ESTIMATE]"));
                }
            }
        }
    }

    // ── Electric field & EFG (PRIMARY descriptors per the curated list) ──
    // signed E|| (the Buckingham linear-term input) + the EFG |T2| invariant
    // (AIMNet2). The raw Coulomb/APBS field vectors and full EFG components stay
    // in the "Electrostatics" drawer below.
    {
        QTreeWidgetItem* ef = nullptr;
        auto ensureEf = [&]() -> QTreeWidgetItem* {
            if (!ef) {
                ef = AddKV(root, QStringLiteral("Electric field & EFG"), QString());
                ef->setExpanded(true);
            }
            return ef;
        };
        if (auto sc = model::QtMopacCoulombGroup(s).scalars(a)) {
            if (std::isfinite(sc->E_bond_proj))
                AddScalarP(ensureEf(), QStringLiteral("signed E_parallel"), sc->E_bond_proj,
                           QStringLiteral("V/Å"),
                           QStringLiteral("MOPAC Coulomb field on the parent-to-H bond axis, "
                                          "signed (V/A); the Buckingham linear-term input. [USE]"));
        }
        if (auto efg = model::QtAimnet2Group(s).efg(a)) {
            if (std::isfinite(efg->t2Magnitude()))
                AddScalarP(ensureEf(), QStringLiteral("EFG |T2| (AIMNet2)"), efg->t2Magnitude(),
                           QStringLiteral("V/Å²"),
                           QStringLiteral("AIMNet2 electric-field-gradient symmetric-traceless invariant "
                                          "|T2| (MOPAC Coulomb is the fallback source). [USE]"));
        }
    }

    // ── Ring current (Biot-Savart / Haigh-Mallion / ring susceptibility) ──
    if (AllowsAny(availability_, {"npy:bs_shielding", "npy:hm_shielding",
                                  "npy:bs_ring_counts"})) {
        model::QtBiotSavartGroup bs(s);
        // raw unit-current kernels -> drawer (the validated ppm form is the
        // viewer-derived ring contribution = bs_per_type_T0 x RingIntensity)
        auto* rk = AddKV(drawer, QStringLiteral("Ring current"), QString());
        bool anyRaw = false;
        anyRaw |= AddOptSpherical(rk, QStringLiteral("bs_shielding"), bs.shielding(a), QStringLiteral("ppm·T/nA"));
        anyRaw |= AddOptSpherical(rk, QStringLiteral("hm_shielding"), model::QtHaighMallionGroup(s).shielding(a), QStringLiteral("Å⁻¹"));
        anyRaw |= AddOptVec3(rk, QStringLiteral("bs_total_B"), bs.totalB(a), QStringLiteral("T"));
        if (!anyRaw) DeleteIfEmpty(rk);
        // ring counts are geometry context -> top
        if (auto rc = bs.ringCounts(a)) {
            auto* g = AddKV(root, QStringLiteral("Ring geometry"), QString());
            AddKV(g, QStringLiteral("ring counts (3/5/8/12 Å)"),
                  QStringLiteral("%1 / %2 / %3 / %4").arg(rc->within3A).arg(rc->within5A).arg(rc->within8A).arg(rc->within12A));
        }
    }


    // ── Bond anisotropy (McConnell) ──
    if (AllowsAny(availability_, {"npy:mc_peptide_co_fixed",
                                  "npy:mc_peptide_cn_fixed",
                                  "npy:mc_backbone_other_fixed",
                                  "npy:mc_sidechain_co_fixed",
                                  "npy:mc_sidechain_other_fixed",
                                  "npy:mc_disulfide_fixed",
                                  "npy:mc_aromatic_fixed",
                                  "npy:mc_backbone_xh_fixed",
                                  "npy:mc_sidechain_xh_fixed",
                                  "npy:mc_s_h_fixed"})) {
        model::QtMcConnellGroup mc(s);
        auto* mk = AddKV(drawer, QStringLiteral("Bond anisotropy (McConnell)"), QString());
        bool anyMk = false;
        static constexpr struct {
            io::FieldKind kind;
            const char* label;
        } kFixedSources[] = {
            {io::FieldKind::McPeptideCoFixed, "peptide C=O"},
            {io::FieldKind::McPeptideCNFixed, "peptide C-N"},
            {io::FieldKind::McBackboneOtherFixed, "backbone other"},
            {io::FieldKind::McSidechainCoFixed, "sidechain C=O"},
            {io::FieldKind::McSidechainOtherFixed, "sidechain other"},
            {io::FieldKind::McDisulfideFixed, "disulfide"},
            {io::FieldKind::McAromaticFixed, "aromatic"},
            {io::FieldKind::McBackboneXhFixed, "backbone X-H"},
            {io::FieldKind::McSidechainXhFixed, "sidechain X-H"},
            {io::FieldKind::McSHFixed, "S-H"},
        };
        for (const auto& source : kFixedSources) {
            anyMk |= AddOptSpherical(
                mk,
                QString::fromLatin1(source.label),
                mc.tensor(source.kind, a),
                QStringLiteral("Angstrom^-3"));
        }
        if (!anyMk) DeleteIfEmpty(mk);
    }

    // ── Electrostatics (Coulomb / APBS / AIMNet2 EFG) ──
    if (AllowsAny(availability_, {"npy:coulomb_shielding", "npy:coulomb_E",
                                  "npy:apbs_E", "npy:apbs_efg",
                                  "npy:aimnet2_efg"})) {
        auto* g = AddKV(drawer, QStringLiteral("Electrostatics"), QString());
        model::QtCoulombGroup coul(s);
        model::QtApbsGroup apbs(s);
        bool any = false;
        any |= AddOptSpherical(g, QStringLiteral("coulomb_shielding"), coul.shielding(a), QStringLiteral("V/Å²"));
        any |= AddOptVec3(g, QStringLiteral("coulomb_E"), coul.E(a), QStringLiteral("V/Å"));
        any |= AddOptVec3(g, QStringLiteral("apbs_E (APBS diagnostic)"), apbs.E(a), QStringLiteral("V/Å"));
        any |= AddOptEfg(g, QStringLiteral("apbs_efg (APBS diagnostic)"), apbs.efg(a), QStringLiteral("V/Å²"));
        any |= AddOptEfg(g, QStringLiteral("aimnet2_efg"), model::QtAimnet2Group(s).efg(a), QStringLiteral("V/Å²"));
        if (!any) DeleteIfEmpty(g);
    }

    // ── H-bond (kernel form) ──
    if (AllowsAny(availability_, {"npy:hbond_scalars"})) {
        auto* g = AddKV(root, QStringLiteral("H-bond"), QString());
        model::QtHBondGroup hb(s);
        bool any = false;
        if (auto sc = hb.scalars(a)) {
            AddScalar(g, QStringLiteral("nearest dist"), sc->nearest_dist, QStringLiteral("Å"));
            AddScalar(g, QStringLiteral("1/r³"), sc->inv_d3);
            AddScalar(g, QStringLiteral("count ≤ 3.5 Å"), sc->count_3_5A);
            any = true;
        }
        if (!any) DeleteIfEmpty(g);
    }

    // ── SASA ──
    if (AllowsAny(availability_, {"npy:atom_sasa", "npy:sasa_normal"})) {
        auto* g = AddKV(root, QStringLiteral("SASA"), QString());
        model::QtSasaGroup sasa(s);
        bool any = false;
        any |= AddOptScalar(g, QStringLiteral("atom_sasa"), sasa.sasa(a), QStringLiteral("Å²"));
        any |= AddOptVec3(g, QStringLiteral("surface normal"), sasa.normal(a));
        if (!any) DeleteIfEmpty(g);
    }

    // ── Water environment ──
    if (AllowsAny(availability_, {"npy:water_efield", "npy:water_efg",
                                  "npy:water_shell_counts", "npy:hydration_shell",
                                  "npy:water_polarization"})) {
        auto* g = AddKV(root, QStringLiteral("Water"), QString());
        model::QtWaterFieldGroup wf(s);
        bool any = false;
        // raw field vectors -> drawer; shell / hydration / polarization (the
        // environmental descriptors) stay at top.
        auto* wk = AddKV(drawer, QStringLiteral("Water field"), QString());
        bool anyWk = false;
        anyWk |= AddOptVec3(wk, QStringLiteral("water_efield"), wf.efield(a), QStringLiteral("V/Å"));
        anyWk |= AddOptEfg(wk, QStringLiteral("water_efg"), wf.efg(a), QStringLiteral("V/Å²"));
        if (!anyWk) DeleteIfEmpty(wk);
        if (auto wc = wf.shellCounts(a)) {
            AddKV(g, QStringLiteral("shell counts (1st/2nd)"), QStringLiteral("%1 / %2").arg(wc->nFirst).arg(wc->nSecond));
            any = true;
        }
        if (auto h = model::QtHydrationGroup(s).shell(a)) {
            AddScalar(g, QStringLiteral("half-shell asymmetry"), h->halfShellAsymmetry);
            AddScalar(g, QStringLiteral("mean water dipole cos"), h->meanWaterDipoleCos);
            AddKV(g, QStringLiteral("nearest ion"),
                  h->hasNearestIon()
                      ? QStringLiteral("%1 Å, q=%2 e").arg(FmtDouble(h->nearestIonDist), FmtDouble(h->nearestIonCharge))
                      : QStringLiteral("none in cutoff"));
            any = true;
        }
        if (auto p = model::QtWaterPolarizationGroup(s).polarization(a)) {
            AddScalar(g, QStringLiteral("dipole alignment"), p->alignment);
            AddScalar(g, QStringLiteral("dipole coherence"), p->coherence);
            any = true;
        }
        if (!any) DeleteIfEmpty(g);
    }

    // ── Charges ──
    if (AllowsAny(availability_, {"npy:aimnet2_charges", "npy:eeq_charges",
                                  "npy:eeq_cn",
                                  "npy:aimnet2_charge_response_gradient"})) {
        auto* g = AddKV(drawer, QStringLiteral("Charges"), QString());
        model::QtAimnet2Group aim(s);
        model::QtEeqGroup eeq(s);
        bool any = false;
        any |= AddOptScalar(g, QStringLiteral("AIMNet2 (Hirshfeld)"), aim.charge(a), QStringLiteral("e"));
        any |= AddOptScalar(g, QStringLiteral("EEQ"), eeq.charge(a), QStringLiteral("e"));
        any |= AddOptScalar(g, QStringLiteral("EEQ coord. number"), eeq.coordinationNumber(a));
        any |= AddOptScalar(g, QStringLiteral("|charge-response grad|"), aim.chargeResponseGradientNorm(a), QStringLiteral("e²/Å"));
        if (!any) DeleteIfEmpty(g);
    }

    // ── DSSP (per-residue, broadcast to the atom) ──
    if (AllowsAny(availability_, {"npy:dssp_backbone", "npy:dssp_ss8"})) {
        if (auto bb = model::QtDsspGroup(s).backbone(a)) {
            auto* g = AddKV(root, QStringLiteral("DSSP (secondary structure)"), QString());
            AddScalar(g, QStringLiteral("φ (neg-IUPAC)"), bb->phi, QStringLiteral("rad"));
            AddScalar(g, QStringLiteral("ψ (neg-IUPAC)"), bb->psi, QStringLiteral("rad"));
            AddScalar(g, QStringLiteral("residue SASA"), bb->sasa, QStringLiteral("Å²"));
            if (auto ss = model::QtDsspGroup(s).ss8(a))
                AddKV(g, QStringLiteral("SS8 class (ordinal)"), QString::number(static_cast<int>(ss->dominant())));
        }
    }

    // ── Planar geometry ──
    if (AllowsAny(availability_, {"npy:pyramidalization", "npy:omega_actual",
                                  "npy:omega_deviation", "npy:omega_is_xpro"})) {
        model::QtPlanarGeometryGroup pg(s);
        const int resIdx = protein_->atom(a).residueIndex;
        auto* g = AddKV(root, QStringLiteral("Planar geometry"), QString());
        bool any = false;
        any |= AddOptScalar(g, QStringLiteral("pyramidalization"), pg.pyramidalization(a), QStringLiteral("Å"));
        if (resIdx >= 0) {
            const std::size_t r = static_cast<std::size_t>(resIdx);
            any |= AddOptScalar(g, QStringLiteral("ω (peptide)"), pg.omegaActual(r), QStringLiteral("rad"));
            any |= AddOptScalar(g, QStringLiteral("ω deviation"), pg.omegaDeviation(r), QStringLiteral("rad"));
            any |= AddOptBool(g, QStringLiteral("X→Pro context"), pg.omegaIsXpro(r));
        }
        if (!any) DeleteIfEmpty(g);
    }

    // ── Energy (per-atom bonded share + whole-frame GROMACS) ──
    {
        if (auto be = model::QtBondedGroup(s).energy(a)) {
            auto* g = AddKV(root, QStringLiteral("Bonded energy (per-atom share)"), QString());
            AddScalar(g, QStringLiteral("total"), be->total, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("bond"), be->bond, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("angle"), be->angle, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("Urey-Bradley"), be->ureyBradley, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("proper dih"), be->proper, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("harmonic improper dih"), be->harmonicImproper, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("periodic improper dih"), be->periodicImproper, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("CMAP"), be->cmap, QStringLiteral("kJ/mol"));
        }
        if (auto ge = model::QtGromacsGroup(s).energy()) {
            auto* g = AddKV(root, QStringLiteral("Frame energy (GROMACS)"), QString());
            AddScalar(g, QStringLiteral("potential"), ge->potential(), QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("temperature"), ge->temperature(), QStringLiteral("K"));
            AddScalar(g, QStringLiteral("pressure"), ge->pressure(), QStringLiteral("bar"));
        }
    }

    // ── MOPAC (PM7+MOZYME; FullFat --mopac runs only) ──
    if (AllowsAny(availability_, {"npy:mopac_charges", "npy:mopac_scalars",
                                  "npy:mopac_coulomb_shielding",
                                  "npy:mopac_mc_shielding", "npy:mopac_global"})) {
        model::QtMopacCoreGroup mopac(s);
        if (auto ms = mopac.scalars(a)) {
            auto* g = AddKV(root, QStringLiteral("MOPAC (PM7+MOZYME)"), QString());
            AddScalar(g, QStringLiteral("charge"), ms->charge, QStringLiteral("e"));
            AddScalar(g, QStringLiteral("s population"), ms->sPop, QStringLiteral("e"));
            AddScalar(g, QStringLiteral("p population"), ms->pPop, QStringLiteral("e"));
            AddScalar(g, QStringLiteral("Wiberg valency"), ms->valency);
            auto* mk = AddKV(drawer, QStringLiteral("MOPAC kernels"), QString());
            bool anyMk = false;
            anyMk |= AddOptSpherical(mk, QStringLiteral("mopac_coulomb_shielding"), model::QtMopacCoulombGroup(s).shielding(a), QStringLiteral("V/Å²"));
            anyMk |= AddOptSpherical(mk, QStringLiteral("mopac_mc_shielding"), model::QtMopacMcConnellGroup(s).shielding(a), QStringLiteral("Å⁻³"));
            if (!anyMk) DeleteIfEmpty(mk);
            if (auto mg = mopac.global())
                AddScalar(g, QStringLiteral("ΔHf (frame)"), mg->heatOfFormation, QStringLiteral("kcal/mol"));
        }
    }

    // ── ORCA reference shielding and Larsen H-bond contribution ──
    {
        model::QtOrcaGroup orca(s);
        if (AllowsAny(availability_, {"orca_dft:total", "orca_dft:diamagnetic",
                                      "orca_dft:paramagnetic"})) {
            if (auto tot = orca.total(a)) {
                auto* g = AddKV(root, QStringLiteral("DFT reference (ORCA)"), QString());
                AddOptSpherical(g, QStringLiteral("σ total"), tot, QStringLiteral("ppm"));
                auto* dk = AddKV(drawer, QStringLiteral("Shielding components (gauge-dependent)"), QString());
                bool anyDk = false;
                anyDk |= AddOptSpherical(dk, QStringLiteral("σ diamagnetic"), orca.diamagnetic(a), QStringLiteral("ppm"));
                anyDk |= AddOptSpherical(dk, QStringLiteral("σ paramagnetic"), orca.paramagnetic(a), QStringLiteral("ppm"));
                if (!anyDk) DeleteIfEmpty(dk);
            }
        }
        model::QtLarsenHBondGroup larsen(s);
        if (AllowsAny(availability_, {"npy:larsen_hbond_shielding",
                                      "npy:larsen_hbond_water_term",
                                      "npy:larsen_hbond_count"})
            && larsen.hasContribution(a)) {
            auto* g = AddKV(root, QStringLiteral("Larsen H-bond (ProCS15 grid)"), QString());
            AddOptSpherical(g, QStringLiteral("Δσ total"), larsen.shielding(a), QStringLiteral("ppm"));
            AddOptScalar(g, QStringLiteral("water term"), larsen.waterTerm(a), QStringLiteral("ppm"));
            AddOptInt(g, QStringLiteral("H-bond pair count"), larsen.count(a));
        }
    }

    qCDebug(cDock).noquote() << "rebuilt | atom=" << a << "| frame=" << t << "| snapshot= resident";
}

namespace {
// Recursively serialize a tree item's children to JSON (field / value / tooltip).
QJsonArray serializeInspectorChildren(const QTreeWidgetItem* item) {
    QJsonArray out;
    for (int i = 0; i < item->childCount(); ++i) {
        const QTreeWidgetItem* c = item->child(i);
        QJsonObject o{{QStringLiteral("field"), c->text(0)},
                      {QStringLiteral("value"), c->text(1)}};
        if (c->data(0, Qt::CheckStateRole).isValid())
            o.insert(QStringLiteral("checked"), c->checkState(0) == Qt::Checked);
        const QString tip = c->toolTip(0);
        if (!tip.isEmpty())
            o.insert(QStringLiteral("tooltip"), tip);
        if (c->childCount() > 0)
            o.insert(QStringLiteral("children"), serializeInspectorChildren(c));
        out.append(o);
    }
    return out;
}
}  // namespace

// Serialize the focused atom's panel tree for the REST harness so the curated
// display + provenance tooltips are programmatically assertable. Read-only.
QJsonArray QtAtomInspectorDock::dumpTree() const {
    QJsonArray out = serializeInspectorChildren(tree_->invisibleRootItem());
    for (const auto& tensor : serializeInspectorChildren(tensorTree_->invisibleRootItem()))
        out.append(tensor);
    return out;
}

}  // namespace h5reader::app

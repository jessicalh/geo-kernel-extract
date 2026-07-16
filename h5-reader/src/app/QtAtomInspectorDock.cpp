#include "QtAtomInspectorDock.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/QtConformationSnapshot.h"
#include "../model/QtResidueNames.h"
#include "../model/CsaShape.h"
#include "../model/TrajectoryConformation.h"

// Typed per-frame group views over the snapshot — the inspector's single
// source for per-frame calculator detail (build 2; tier-mirror memory).
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
#include "../model/QtTripeptideGroup.h"
#include "../model/QtWaterFieldGroup.h"
#include "../model/QtWaterPolarizationGroup.h"
#include "../physics/ClassicalSourceMath.h"
#include "../physics/LiteratureAccessors.h"
#include "../physics/RingCurrentScalars.h"
#include "../physics/SphericalBasis.h"
#include "constants/LiteratureConstants.h"

#include <QBrush>
#include <QColor>
#include <QFont>
#include <QHeaderView>
#include <QIcon>
#include <QJsonArray>
#include <QJsonObject>
#include <QLoggingCategory>
#include <QPainter>
#include <QPixmap>
#include <QSizePolicy>
#include <QString>
#include <QStringList>
#include <QTreeWidget>
#include <QTreeWidgetItem>

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
// value text short enough to read at a glance; expand on a child
// node if the user wants details (future).

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

    static constexpr struct { const char* name; double r, g, b; } kAxes[3] = {
        {"11", 0.96, 0.66, 0.16},  // amber
        {"22", 0.18, 0.74, 0.74},  // teal
        {"33", 0.74, 0.36, 0.86},  // violet
    };
    const double vals[3] = {shape.principal_values[0], shape.principal_values[1], shape.principal_values[2]};
    for (int i = 0; i < 3; ++i) {
        const QColor color = QColor::fromRgbF(kAxes[i].r, kAxes[i].g, kAxes[i].b);
        AddSwatchKV(parent,
                    QStringLiteral("%1_%2").arg(valuePrefix, QString::fromLatin1(kAxes[i].name)),
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
    if (!v) return false;
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
    QFont compactFont = font();
    if (compactFont.pointSize() > 8)
        compactFont.setPointSize(compactFont.pointSize() - 1);
    else if (compactFont.pixelSize() > 10)
        compactFont.setPixelSize(compactFont.pixelSize() - 1);
    setFont(compactFont);

    tree_ = new QTreeWidget(this);
    tree_->setMinimumWidth(0);
    tree_->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Expanding);
    tree_->setColumnCount(2);
    tree_->setHeaderLabels({QStringLiteral("Field"), QStringLiteral("Value")});
    tree_->setAlternatingRowColors(true);
    tree_->header()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    setWidget(tree_);

    // Starting placeholder.
    auto* hint = new QTreeWidgetItem(tree_);
    hint->setText(0, QStringLiteral("Double-click an atom in the viewport"));
}

void QtAtomInspectorDock::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    if (conformation_)
        disconnect(conformation_.data(), nullptr, this, nullptr);
    protein_ = protein;
    conformation_ = conformation;
    if (conformation_) {
        ACONNECT(conformation_.data(), &model::Conformation::snapshotReady,
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
    hasSelection_ = true;
    atomIdx_ = atomIdx;
    const std::size_t frame = static_cast<std::size_t>(std::max(0, frame_));
    if (conformation_) {
        conformation_->requestSnapshot(frame);
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
    tree_->clear();
    auto* hint = new QTreeWidgetItem(tree_);
    hint->setText(0, QStringLiteral("Double-click an atom in the viewport"));
}

void QtAtomInspectorDock::setCsaTensor(std::size_t atom, const CsaTensorInfo& info) {
    csa_ = info;
    csaAtom_ = atom;
    hasCsa_ = true;
    if (hasSelection_)
        rebuild();
}

void QtAtomInspectorDock::clearCsaTensor() {
    if (!hasCsa_)
        return;
    hasCsa_ = false;
    if (hasSelection_)
        rebuild();
}

void QtAtomInspectorDock::setOrientationTensor(std::size_t atom, const OrientationTensorInfo& info) {
    orient_ = info;
    orientAtom_ = atom;
    hasOrient_ = true;
    if (hasSelection_)
        rebuild();
}

void QtAtomInspectorDock::clearOrientationTensor() {
    if (!hasOrient_)
        return;
    hasOrient_ = false;
    if (hasSelection_)
        rebuild();
}

void QtAtomInspectorDock::populateCsa(QTreeWidgetItem* root) {
    auto* group = AddKV(root, QStringLiteral("CSA shielding tensor (DFT)"),
                        csa_.framed ? csa_.frameKind : QStringLiteral("unframed"));
    group->setExpanded(true);
    AddScalar(group, QStringLiteral("sigma_iso"), csa_.sigmaIso, QStringLiteral("ppm"));
    AddScalar(group, QStringLiteral("span"), csa_.span, QStringLiteral("ppm"));
    AddScalar(group, QStringLiteral("skew"), csa_.skew);
    AddScalar(group, QStringLiteral("eta"), csa_.eta);

    // Per-axis principal values. The swatch is the colour key for the scene
    // arrows; no floating labels in the molecule view.
    static constexpr struct { const char* name; double r, g, b; } kAxes[3] = {
        {"sigma_11", 0.96, 0.66, 0.16},  // amber
        {"sigma_22", 0.18, 0.74, 0.74},  // teal
        {"sigma_33", 0.74, 0.36, 0.86},  // violet
    };
    const double vals[3] = {csa_.sigma11, csa_.sigma22, csa_.sigma33};
    for (int i = 0; i < 3; ++i) {
        const QColor color = QColor::fromRgbF(kAxes[i].r, kAxes[i].g, kAxes[i].b);
        AddSwatchKV(group, QString::fromLatin1(kAxes[i].name),
                    FmtDouble(vals[i]) + QStringLiteral(" ppm"), color);
    }
}

void QtAtomInspectorDock::populateOrientation(QTreeWidgetItem* root) {
    auto* group = AddKV(root, QStringLiteral("Bond orientation tensor"), orient_.bond);
    group->setExpanded(true);
    AddScalar(group, QStringLiteral("S^2 (order parameter)"), orient_.s2);

    // Order-tensor eigenvalues (descending; sum to 1). The swatch is the
    // colour key for the scene arrows; no floating labels in the molecule view.
    static constexpr struct { const char* name; double r, g, b; } kAxes[3] = {
        {"lambda_1", 0.96, 0.66, 0.16},  // amber
        {"lambda_2", 0.18, 0.74, 0.74},  // teal
        {"lambda_3", 0.74, 0.36, 0.86},  // violet
    };
    const double vals[3] = {orient_.lambda1, orient_.lambda2, orient_.lambda3};
    for (int i = 0; i < 3; ++i) {
        const QColor color = QColor::fromRgbF(kAxes[i].r, kAxes[i].g, kAxes[i].b);
        AddSwatchKV(group, QString::fromLatin1(kAxes[i].name), FmtDouble(vals[i]), color);
    }
}

void QtAtomInspectorDock::rebuild() {
    if (!tree_ || !protein_ || !conformation_)
        return;
    if (!hasSelection_ || atomIdx_ >= protein_->atomCount())
        return;

    // Batch the rebuild into a single repaint: this tree is cleared + fully
    // repopulated on every focus / frame change, so per-item updates would flicker.
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
    // last (below CSA / orientation) iff it ended up with content.
    auto* drawer = new QTreeWidgetItem();
    drawer->setText(0, QStringLiteral("Raw kernels & diagnostics"));
    drawer->setText(1, QStringLiteral("raw npy inputs (not validated)"));
    populatePerFrame(title, drawer);
    if (hasCsa_ && csaAtom_ == atomIdx_)
        populateCsa(title);
    if (hasOrient_ && orientAtom_ == atomIdx_)
        populateOrientation(title);
    if (drawer->childCount() > 0) {
        title->addChild(drawer);     // tree takes ownership; collapse AFTER attach
        drawer->setExpanded(false);
    } else {
        delete drawer;               // orphan, never attached -- we own it
    }

    tree_->setUpdatesEnabled(true);
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

    // The full per-frame calculator pile comes from the snapshot via the typed
    // group views. ONE SOURCE PER ROLE: this panel reads the snapshot only
    // (tier-mirror memory), never the H5 time series. A group view's nullopt is
    // "this calculator did not run this frame" (absent, not faked) → em-dash.
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
                {io::FieldKind::McAromaticZeroedBo, model::BondCategory::Aromatic},
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
        // primary bond axis) + element/frame-specific A,B literature constants. The
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
                           QStringLiteral("MOPAC Coulomb field on the local bond / molecular-z axis, "
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
    if (AllowsAny(availability_, {"npy:mc_shielding", "npy:mc_category_T2",
                                  "npy:mc_scalars"})) {
        model::QtMcConnellGroup mc(s);
        // raw geometry kernel + angular sums -> drawer (the validated ppm form is
        // the viewer-derived McConnell contribution, _bo signed T0 x DeltaChi)
        auto* mk = AddKV(drawer, QStringLiteral("Bond anisotropy (McConnell)"), QString());
        bool anyMk = false;
        anyMk |= AddOptSpherical(mk, QStringLiteral("mc_shielding"), mc.shielding(a), QStringLiteral("Å⁻³"));
        if (auto sc = mc.scalars(a)) {
            AddScalar(mk, QStringLiteral("Σ C=O angular"), sc->co_sum);
            AddScalar(mk, QStringLiteral("Σ C–N angular"), sc->cn_sum);
            AddScalar(mk, QStringLiteral("Σ sidechain"), sc->sidechain_sum);
            AddScalar(mk, QStringLiteral("Σ aromatic"), sc->aromatic_sum);
            anyMk = true;
            // nearest carbonyl / amide are geometry -> top
            auto* g = AddKV(root, QStringLiteral("Nearest backbone partners"), QString());
            AddKV(g, QStringLiteral("nearest C=O"),
                  sc->hasNearestCO() ? FmtDouble(sc->nearest_CO_dist) + QStringLiteral(" Å") : QStringLiteral("none"));
            AddKV(g, QStringLiteral("nearest C–N"),
                  sc->hasNearestCN() ? FmtDouble(sc->nearest_CN_dist) + QStringLiteral(" Å") : QStringLiteral("none"));
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
        if (AllowsAny(availability_, {"npy:bonded_energy"})) {
            if (auto be = model::QtBondedGroup(s).energy(a)) {
                auto* g = AddKV(root, QStringLiteral("Bonded energy (per-atom share)"), QString());
                AddScalar(g, QStringLiteral("total"), be->total, QStringLiteral("kJ/mol"));
                AddScalar(g, QStringLiteral("bond"), be->bond, QStringLiteral("kJ/mol"));
                AddScalar(g, QStringLiteral("angle"), be->angle, QStringLiteral("kJ/mol"));
                AddScalar(g, QStringLiteral("proper dih"), be->proper, QStringLiteral("kJ/mol"));
                AddScalar(g, QStringLiteral("improper dih"), be->improper, QStringLiteral("kJ/mol"));
            }
        }
        if (AllowsAny(availability_, {"npy:gromacs_energy"})) {
            if (auto ge = model::QtGromacsGroup(s).energy()) {
                auto* g = AddKV(root, QStringLiteral("Frame energy (GROMACS)"), QString());
                AddScalar(g, QStringLiteral("potential"), ge->potential(), QStringLiteral("kJ/mol"));
                AddScalar(g, QStringLiteral("temperature"), ge->temperature(), QStringLiteral("K"));
                AddScalar(g, QStringLiteral("pressure"), ge->pressure(), QStringLiteral("bar"));
            }
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

    // ── DFT / ProCS15 reference shielding (single-pose --orca + tripeptide / Larsen) ──
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
        model::QtTripeptideGroup trip(s);
        if (AllowsAny(availability_, {"npy:tripeptide_bb_shielding",
                                      "npy:tripeptide_bb_match_distance",
                                      "npy:tripeptide_neighbor_shielding"})
            && trip.hasBackboneMatch(a)) {
            auto* g = AddKV(root, QStringLiteral("Tripeptide reference (ProCS15)"), QString());
            AddOptSpherical(g, QStringLiteral("backbone σ"), trip.backboneShielding(a), QStringLiteral("ppm"));
            AddOptScalar(g, QStringLiteral("match distance"), trip.backboneMatchDistance(a), QStringLiteral("Å"));
            AddOptSpherical(g, QStringLiteral("neighbor σ (i±1)"), trip.neighborShielding(a), QStringLiteral("ppm"));
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

    tree_->expandToDepth(1);
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
    QJsonArray out;
    if (!tree_)
        return out;
    for (int i = 0; i < tree_->topLevelItemCount(); ++i) {
        const QTreeWidgetItem* top = tree_->topLevelItem(i);
        QJsonObject o{{QStringLiteral("field"), top->text(0)},
                      {QStringLiteral("value"), top->text(1)}};
        const QString tip = top->toolTip(0);
        if (!tip.isEmpty())
            o.insert(QStringLiteral("tooltip"), tip);
        if (top->childCount() > 0)
            o.insert(QStringLiteral("children"), serializeInspectorChildren(top));
        out.append(o);
    }
    return out;
}

}  // namespace h5reader::app

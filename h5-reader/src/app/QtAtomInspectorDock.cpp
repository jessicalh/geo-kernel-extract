#include "QtAtomInspectorDock.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../io/QtNamingRegistry.h"
#include "../model/QtFrame.h"

#include <QHeaderView>
#include <QLoggingCategory>
#include <QString>
#include <QTreeWidget>
#include <QTreeWidgetItem>

#include <algorithm>
#include <cmath>

namespace h5reader::app {

using model::SphericalTensor;
using model::Vec3;

namespace {
Q_LOGGING_CATEGORY(cDock, "h5reader.inspector")

// Formatting helpers — two-column tree with Field / Value. Keep the
// value text short enough to read at a glance; expand on a child
// node if the user wants details (future).

QString FmtDouble(double v, int precision = 4) {
    if (!std::isfinite(v)) return QStringLiteral("nan");
    return QStringLiteral("%1").arg(v, 0, 'g', precision);
}

QString FmtVec3(const Vec3& v, const QString& unit = QString(),
                int precision = 4) {
    return QStringLiteral("(%1, %2, %3)%4")
        .arg(FmtDouble(v.x(), precision),
             FmtDouble(v.y(), precision),
             FmtDouble(v.z(), precision),
             unit.isEmpty() ? QString() : QStringLiteral(" ") + unit);
}

QString FmtSphericalSummary(const SphericalTensor& st) {
    double t2Mag = 0.0;
    for (double c : st.T2) t2Mag += c * c;
    t2Mag = std::sqrt(t2Mag);
    return QStringLiteral("T0=%1  |T2|=%2")
        .arg(FmtDouble(st.T0, 5), FmtDouble(t2Mag, 5));
}

QTreeWidgetItem* AddKV(QTreeWidgetItem* parent,
                        const QString& field,
                        const QString& value) {
    auto* it = new QTreeWidgetItem(parent);
    it->setText(0, field);
    it->setText(1, value);
    return it;
}

void AddScalar(QTreeWidgetItem* parent, const QString& name,
                double value, const QString& unit = QString()) {
    AddKV(parent, name,
          unit.isEmpty() ? FmtDouble(value)
                         : FmtDouble(value) + QStringLiteral(" ") + unit);
}

void AddVec3(QTreeWidgetItem* parent, const QString& name,
              const Vec3& v, const QString& unit = QString()) {
    AddKV(parent, name, FmtVec3(v, unit));
}

void AddSpherical(QTreeWidgetItem* parent, const QString& name,
                   const SphericalTensor& st, const QString& unit = QString()) {
    auto* it = AddKV(parent, name,
        FmtSphericalSummary(st) + (unit.isEmpty()
                                     ? QString()
                                     : QStringLiteral(" ") + unit));
    // T2 component breakdown as a child row each — useful for
    // verifying the angular decomposition frame-to-frame.
    auto* t2 = AddKV(it, QStringLiteral("T2 components"), QString());
    for (int i = 0; i < 5; ++i) {
        AddScalar(t2, QStringLiteral("m=%1").arg(i - 2), st.T2[i]);
    }
    AddVec3(it, QStringLiteral("T1 (antisym)"),
            Vec3(st.T1[0], st.T1[1], st.T1[2]));
}

const char* AtomRoleName(model::AtomRole r) { return model::NameForAtomRole(r); }
const char* HybName(model::Hybridisation h) { return model::NameForHybridisation(h); }

}  // namespace

QtAtomInspectorDock::QtAtomInspectorDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Atom Inspector"), parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtAtomInspectorDock"));
    setFeatures(QDockWidget::DockWidgetMovable
                | QDockWidget::DockWidgetFloatable);

    tree_ = new QTreeWidget(this);
    tree_->setColumnCount(2);
    tree_->setHeaderLabels({QStringLiteral("Field"), QStringLiteral("Value")});
    tree_->setAlternatingRowColors(true);
    tree_->header()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    setWidget(tree_);

    // Starting placeholder.
    auto* hint = new QTreeWidgetItem(tree_);
    hint->setText(0, QStringLiteral("Double-click an atom in the viewport"));
}

void QtAtomInspectorDock::setContext(const model::QtProtein*      protein,
                                       const model::QtConformation* conformation) {
    protein_      = protein;
    conformation_ = conformation;
}

void QtAtomInspectorDock::setPickedAtom(std::size_t atomIdx) {
    ASSERT_THREAD(this);
    hasSelection_ = true;
    atomIdx_      = atomIdx;
    rebuild();
}

void QtAtomInspectorDock::setFrame(int t) {
    ASSERT_THREAD(this);
    frame_ = t;
    if (hasSelection_) rebuild();
}

void QtAtomInspectorDock::clearSelection() {
    ASSERT_THREAD(this);
    hasSelection_ = false;
    tree_->clear();
    auto* hint = new QTreeWidgetItem(tree_);
    hint->setText(0, QStringLiteral("Double-click an atom in the viewport"));
}

void QtAtomInspectorDock::rebuild() {
    if (!tree_ || !protein_ || !conformation_) return;
    if (!hasSelection_ || atomIdx_ >= protein_->atomCount()) return;

    tree_->clear();

    auto* title = new QTreeWidgetItem(tree_);
    const auto& atom = protein_->atom(atomIdx_);
    const auto& res  = atom.residueIndex >= 0
                       ? protein_->residue(atom.residueIndex)
                       : model::QtResidue{};
    title->setText(0,
        QStringLiteral("Atom %1 — %2 %3 #%4")
            .arg(atomIdx_)
            .arg(atom.h5AtomName,
                 QString::fromLatin1(io::ThreeLetterCodeForAminoAcid(res.aminoAcid)))
            .arg(res.residueNumber));
    title->setText(1, QStringLiteral("frame %1 / %2")
        .arg(frame_ + 1)
        .arg(conformation_->frameCount()));
    title->setExpanded(true);

    populateIdentity(title);
    populatePerFrame(title);
}

void QtAtomInspectorDock::populateIdentity(QTreeWidgetItem* parent) {
    const auto& atom = protein_->atom(atomIdx_);
    const auto& res  = atom.residueIndex >= 0
                       ? protein_->residue(atom.residueIndex)
                       : model::QtResidue{};

    auto* g = AddKV(parent, QStringLiteral("Identity"), QString());
    g->setExpanded(true);

    AddKV(g, QStringLiteral("Element"),
          QString::fromLatin1(model::SymbolForElement(atom.element)));
    AddKV(g, QStringLiteral("H5 name"), atom.h5AtomName);
    AddKV(g, QStringLiteral("Role"),
          QString::fromLatin1(AtomRoleName(atom.role)));
    AddKV(g, QStringLiteral("Hybridisation"),
          QString::fromLatin1(HybName(atom.hybridisation)));
    AddKV(g, QStringLiteral("Residue"),
          QStringLiteral("%1 #%2")
              .arg(QString::fromLatin1(
                  io::ThreeLetterCodeForAminoAcid(res.aminoAcid)),
                   QString::number(res.residueNumber)));
    AddKV(g, QStringLiteral("Chain"),
          res.chainId.isEmpty() ? QStringLiteral("—") : res.chainId);
    AddKV(g, QStringLiteral("Variant"),
          QString::fromLatin1(io::ProtonationVariantName(res.variant)));
    AddScalar(g, QStringLiteral("Partial charge"), atom.partialCharge, QStringLiteral("e"));
    AddScalar(g, QStringLiteral("vdW radius"),     atom.vdwRadius,     QStringLiteral("Å"));

    auto* flags = AddKV(g, QStringLiteral("Flags"), QString());
    AddKV(flags, QStringLiteral("is_backbone"),
          atom.isBackbone ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_amide_H"),
          atom.isAmideH ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_alpha_H"),
          atom.isAlphaH ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_methyl"),
          atom.isMethyl ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_aromatic_H"),
          atom.isAromaticH ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_hbond_donor"),
          atom.isHBondDonor ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("is_hbond_acceptor"),
          atom.isHBondAcceptor ? QStringLiteral("true") : QStringLiteral("false"));
    AddKV(flags, QStringLiteral("parent_is_sp2"),
          atom.parentIsSp2 ? QStringLiteral("true") : QStringLiteral("false"));
}

void QtAtomInspectorDock::populatePerFrame(QTreeWidgetItem* root) {
    const int T = static_cast<int>(conformation_->frameCount());
    const int t = std::clamp(frame_, 0, T - 1);
    const auto& frame = conformation_->frame(static_cast<size_t>(t));

    // Position at this frame.
    auto* posG = AddKV(root, QStringLiteral("Position"), QString());
    AddVec3(posG, QStringLiteral("xyz"), frame.position(atomIdx_),
            QStringLiteral("Å"));

    // Ring-current shielding.
    {
        auto* g = AddKV(root, QStringLiteral("Ring current"), QString());
        AddSpherical(g, QStringLiteral("bs_shielding"), frame.bsShielding(atomIdx_),
                     QStringLiteral("ppm"));
        AddSpherical(g, QStringLiteral("hm_shielding"), frame.hmShielding(atomIdx_),
                     QStringLiteral("ppm"));
        AddSpherical(g, QStringLiteral("rs_shielding"), frame.rsShielding(atomIdx_),
                     QStringLiteral("ppm"));
        AddVec3(g, QStringLiteral("total_B_field"), frame.totalBField(atomIdx_));
        AddScalar(g, QStringLiteral("n_rings ≤ 3 Å"), frame.nRings3A(atomIdx_));
        AddScalar(g, QStringLiteral("n_rings ≤ 5 Å"), frame.nRings5A(atomIdx_));
        AddScalar(g, QStringLiteral("n_rings ≤ 8 Å"), frame.nRings8A(atomIdx_));
        AddScalar(g, QStringLiteral("mean ring dist"),   frame.meanRingDist(atomIdx_),
                  QStringLiteral("Å"));
        AddScalar(g, QStringLiteral("nearest ring atom"), frame.nearestRingAtom(atomIdx_),
                  QStringLiteral("Å"));
    }

    // Bond anisotropy.
    {
        auto* g = AddKV(root, QStringLiteral("Bond anisotropy (McConnell)"), QString());
        AddSpherical(g, QStringLiteral("mc_shielding"), frame.mcShielding(atomIdx_),
                     QStringLiteral("Å⁻³"));
        AddScalar(g, QStringLiteral("co_sum"), frame.mcCOSum(atomIdx_));
        AddScalar(g, QStringLiteral("nearest C=O dist"),
                  frame.mcNearestCODist(atomIdx_), QStringLiteral("Å"));
        AddVec3(g, QStringLiteral("nearest C=O direction"),
                frame.mcNearestCODir(atomIdx_));
    }

    // Quadrupole + dispersion.
    {
        auto* g = AddKV(root, QStringLiteral("Quadrupole / Dispersion"), QString());
        AddSpherical(g, QStringLiteral("pq_shielding"),   frame.pqShielding(atomIdx_));
        AddSpherical(g, QStringLiteral("disp_shielding"), frame.dispShielding(atomIdx_));
    }

    // Electrostatics.
    {
        auto* g = AddKV(root, QStringLiteral("Electrostatics"), QString());
        AddSpherical(g, QStringLiteral("coulomb_shielding"),
                     frame.coulombShielding(atomIdx_), QStringLiteral("ppm"));
        AddVec3(g, QStringLiteral("coulomb_E_total"),
                frame.coulombETotal(atomIdx_), QStringLiteral("V/Å"));
        AddScalar(g, QStringLiteral("coulomb_E_magnitude"),
                  frame.coulombEMagnitude(atomIdx_), QStringLiteral("V/Å"));
        AddSpherical(g, QStringLiteral("apbs_efg"), frame.apbsEfg(atomIdx_),
                     QStringLiteral("V/Å²"));
        AddVec3(g, QStringLiteral("apbs_efield"), frame.apbsEfield(atomIdx_),
                QStringLiteral("V/Å"));
        AddSpherical(g, QStringLiteral("aimnet2_shielding"),
                     frame.aimnet2Shielding(atomIdx_), QStringLiteral("ppm"));
    }

    // H-bond.
    {
        auto* g = AddKV(root, QStringLiteral("H-bond"), QString());
        AddSpherical(g, QStringLiteral("hbond_shielding"),
                     frame.hbondShielding(atomIdx_));
        AddScalar(g, QStringLiteral("nearest dist"),
                  frame.hbondNearestDist(atomIdx_), QStringLiteral("Å"));
        AddVec3(g, QStringLiteral("nearest direction"),
                frame.hbondNearestDir(atomIdx_));
        AddScalar(g, QStringLiteral("count ≤ 3.5 Å"),
                  frame.hbondCount35A(atomIdx_));
        AddKV(g, QStringLiteral("is_donor (this frame)"),
              frame.hbondIsDonor(atomIdx_) ? QStringLiteral("true") : QStringLiteral("false"));
        AddKV(g, QStringLiteral("is_acceptor (this frame)"),
              frame.hbondIsAcceptor(atomIdx_) ? QStringLiteral("true") : QStringLiteral("false"));
    }

    // SASA.
    {
        auto* g = AddKV(root, QStringLiteral("SASA"), QString());
        AddScalar(g, QStringLiteral("sasa"), frame.sasa(atomIdx_), QStringLiteral("Å²"));
        AddVec3(g, QStringLiteral("normal"), frame.sasaNormal(atomIdx_));
    }

    // Water environment.
    {
        auto* g = AddKV(root, QStringLiteral("Water"), QString());
        AddVec3(g, QStringLiteral("efield"),
                frame.waterEfield(atomIdx_), QStringLiteral("V/Å"));
        AddScalar(g, QStringLiteral("n_first shell"),  frame.waterNFirst(atomIdx_));
        AddScalar(g, QStringLiteral("n_second shell"), frame.waterNSecond(atomIdx_));
        AddScalar(g, QStringLiteral("half-shell asymmetry"),
                  frame.waterHalfShellAsymmetry(atomIdx_));
        AddScalar(g, QStringLiteral("dipole cos"), frame.waterDipoleCos(atomIdx_));
    }

    // Charges.
    {
        auto* g = AddKV(root, QStringLiteral("Charges"), QString());
        AddScalar(g, QStringLiteral("AIMNet2"), frame.aimnet2Charge(atomIdx_),
                  QStringLiteral("e"));
        AddScalar(g, QStringLiteral("EEQ"),     frame.eeqCharge(atomIdx_),
                  QStringLiteral("e"));
        AddScalar(g, QStringLiteral("EEQ CN"),  frame.eeqCoordinationNumber(atomIdx_));
    }

    tree_->expandToDepth(1);
    qCDebug(cDock).noquote()
        << "rebuilt | atom=" << atomIdx_ << "| frame=" << t;
}

}  // namespace h5reader::app

#include "QtAtomInspectorDock.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/QtConformationSnapshot.h"
#include "../model/QtResidueNames.h"
#include "../model/TrajectoryConformation.h"

// Typed per-frame group views over the snapshot — the inspector's single
// source for per-frame calculator detail (build 2; tier-mirror memory).
#include "../model/QtAimnet2Group.h"
#include "../model/QtApbsGroup.h"
#include "../model/QtBiotSavartGroup.h"
#include "../model/QtBondedGroup.h"
#include "../model/QtCoulombGroup.h"
#include "../model/QtDispersionGroup.h"
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
#include "../model/QtPiQuadrupoleGroup.h"
#include "../model/QtPlanarGeometryGroup.h"
#include "../model/QtRingSusceptibilityGroup.h"
#include "../model/QtSasaGroup.h"
#include "../model/QtTripeptideGroup.h"
#include "../model/QtWaterFieldGroup.h"
#include "../model/QtWaterPolarizationGroup.h"

#include <QFont>
#include <QHeaderView>
#include <QLoggingCategory>
#include <QSizePolicy>
#include <QString>
#include <QTreeWidget>
#include <QTreeWidgetItem>

#include <algorithm>
#include <cmath>

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

QString FmtVec3(const Vec3& v, const QString& unit = QString(), int precision = 4) {
    return QStringLiteral("(%1, %2, %3)%4")
        .arg(FmtDouble(v.x(), precision),
             FmtDouble(v.y(), precision),
             FmtDouble(v.z(), precision),
             unit.isEmpty() ? QString() : QStringLiteral(" ") + unit);
}

[[maybe_unused]] QString FmtSphericalSummary(const SphericalTensor& st) {
    double t2Mag = 0.0;
    for (double c : st.T2)
        t2Mag += c * c;
    t2Mag = std::sqrt(t2Mag);
    return QStringLiteral("T0=%1  |T2|=%2").arg(FmtDouble(st.T0, 5), FmtDouble(t2Mag, 5));
}

QTreeWidgetItem* AddKV(QTreeWidgetItem* parent, const QString& field, const QString& value) {
    auto* it = new QTreeWidgetItem(parent);
    it->setText(0, field);
    it->setText(1, value);
    return it;
}

void AddScalar(QTreeWidgetItem* parent, const QString& name, double value, const QString& unit = QString()) {
    AddKV(parent, name, unit.isEmpty() ? FmtDouble(value) : FmtDouble(value) + QStringLiteral(" ") + unit);
}

void AddVec3(QTreeWidgetItem* parent, const QString& name, const Vec3& v, const QString& unit = QString()) {
    AddKV(parent, name, FmtVec3(v, unit));
}

[[maybe_unused]] void AddSpherical(QTreeWidgetItem* parent, const QString& name, const SphericalTensor& st, const QString& unit = QString()) {
    auto* it = AddKV(parent, name, FmtSphericalSummary(st) + (unit.isEmpty() ? QString() : QStringLiteral(" ") + unit));
    // T2 component breakdown as a child row each — useful for
    // verifying the angular decomposition frame-to-frame.
    auto* t2 = AddKV(it, QStringLiteral("T2 components"), QString());
    for (int i = 0; i < 5; ++i) {
        AddScalar(t2, QStringLiteral("m=%1").arg(i - 2), st.T2[i]);
    }
    AddVec3(it, QStringLiteral("T1 (antisym)"), Vec3(st.T1[0], st.T1[1], st.T1[2]));
}

// Optional-aware adders — show the value, or an em-dash when the group view
// returned nullopt ("this calculator did not run this frame"; absent, not faked).
void AddOptScalar(QTreeWidgetItem* p, const QString& name, const std::optional<double>& v, const QString& unit = QString()) {
    if (v) AddScalar(p, name, *v, unit);
    else AddKV(p, name, QStringLiteral("—"));
}
void AddOptVec3(QTreeWidgetItem* p, const QString& name, const std::optional<Vec3>& v, const QString& unit = QString()) {
    if (v) AddVec3(p, name, *v, unit);
    else AddKV(p, name, QStringLiteral("—"));
}
void AddOptSpherical(QTreeWidgetItem* p, const QString& name, const std::optional<SphericalTensor>& v, const QString& unit = QString()) {
    if (v) AddSpherical(p, name, *v, unit);
    else AddKV(p, name, QStringLiteral("—"));
}
void AddOptEfg(QTreeWidgetItem* p, const QString& name, const std::optional<QtEfg>& v, const QString& unit = QString()) {
    if (!v) {
        AddKV(p, name, QStringLiteral("—"));
        return;
    }
    auto* it = AddKV(p, name,
                     QStringLiteral("|T2|=%1%2").arg(FmtDouble(v->t2Magnitude(), 5),
                                                     unit.isEmpty() ? QString() : QStringLiteral(" ") + unit));
    auto* t2 = AddKV(it, QStringLiteral("T2 components"), QString());
    for (int i = 0; i < 5; ++i)
        AddScalar(t2, QStringLiteral("m=%1").arg(i - 2), v->t2[i]);
}
void AddOptInt(QTreeWidgetItem* p, const QString& name, const std::optional<int>& v) {
    AddKV(p, name, v ? QString::number(*v) : QStringLiteral("—"));
}
void AddOptBool(QTreeWidgetItem* p, const QString& name, const std::optional<bool>& v) {
    AddKV(p, name, v ? (*v ? QStringLiteral("true") : QStringLiteral("false")) : QStringLiteral("—"));
}

}  // namespace

QtAtomInspectorDock::QtAtomInspectorDock(QWidget* parent) : QDockWidget(QStringLiteral("Atom Info"), parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtAtomInspectorDock"));
    setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    setMinimumWidth(48);
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
    protein_ = protein;
    conformation_ = conformation;
    if (conformation_) {
        ACONNECT(conformation_.data(), &model::Conformation::snapshotReady,
                 this, &QtAtomInspectorDock::onSnapshotReady);
    }
}

void QtAtomInspectorDock::setPickedAtom(std::size_t atomIdx) {
    ASSERT_THREAD(this);
    hasSelection_ = true;
    atomIdx_ = atomIdx;
    if (conformation_)
        conformation_->requestSnapshot(static_cast<std::size_t>(std::max(0, frame_)));
    rebuild();
}

void QtAtomInspectorDock::setFrame(int t) {
    ASSERT_THREAD(this);
    frame_ = t;
    if (hasSelection_) {
        if (conformation_)
            conformation_->requestSnapshot(static_cast<std::size_t>(std::max(0, t)));
        rebuild();
    }
}

void QtAtomInspectorDock::onSnapshotReady(std::size_t frame) {
    ASSERT_THREAD(this);
    if (hasSelection_ && static_cast<int>(frame) == frame_)
        rebuild();
}

void QtAtomInspectorDock::clearSelection() {
    ASSERT_THREAD(this);
    hasSelection_ = false;
    tree_->clear();
    auto* hint = new QTreeWidgetItem(tree_);
    hint->setText(0, QStringLiteral("Double-click an atom in the viewport"));
}

void QtAtomInspectorDock::rebuild() {
    if (!tree_ || !protein_ || !conformation_)
        return;
    if (!hasSelection_ || atomIdx_ >= protein_->atomCount())
        return;

    tree_->clear();

    auto* title = new QTreeWidgetItem(tree_);
    const auto& atom = protein_->atom(atomIdx_);
    const auto& res = atom.residueIndex >= 0 ? protein_->residue(atom.residueIndex) : model::QtResidue{};
    title->setText(
        0,
        QStringLiteral("Atom %1 — %2 %3 #%4")
            .arg(atomIdx_)
            .arg(protein_->atomNames(atomIdx_).amber, QString::fromLatin1(model::IupacResidue3LetterFor(res.aminoAcid)))
            .arg(res.address.residueNumber));
    title->setText(1, QStringLiteral("frame %1 / %2").arg(frame_ + 1).arg(conformation_->frameCount()));
    title->setExpanded(true);

    populateIdentity(title);
    populatePerFrame(title);
}

void QtAtomInspectorDock::populateIdentity(QTreeWidgetItem* parent) {
    const auto& atom = protein_->atom(atomIdx_);
    const auto& res = atom.residueIndex >= 0 ? protein_->residue(atom.residueIndex) : model::QtResidue{};

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

void QtAtomInspectorDock::populatePerFrame(QTreeWidgetItem* root) {
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

    // ── Ring current (Biot-Savart / Haigh-Mallion / ring susceptibility) ──
    {
        auto* g = AddKV(root, QStringLiteral("Ring current"), QString());
        model::QtBiotSavartGroup bs(s);
        AddOptSpherical(g, QStringLiteral("bs_shielding"), bs.shielding(a), QStringLiteral("ppm·T/nA"));
        AddOptSpherical(g, QStringLiteral("hm_shielding"), model::QtHaighMallionGroup(s).shielding(a), QStringLiteral("Å⁻¹"));
        AddOptSpherical(g, QStringLiteral("ringchi_shielding"), model::QtRingSusceptibilityGroup(s).shielding(a), QStringLiteral("Å⁻³"));
        AddOptVec3(g, QStringLiteral("bs_total_B"), bs.totalB(a), QStringLiteral("T"));
        if (auto rc = bs.ringCounts(a))
            AddKV(g, QStringLiteral("ring counts (3/5/8/12 Å)"),
                  QStringLiteral("%1 / %2 / %3 / %4").arg(rc->within3A).arg(rc->within5A).arg(rc->within8A).arg(rc->within12A));
    }

    // ── π-quadrupole / dispersion ──
    {
        auto* g = AddKV(root, QStringLiteral("Quadrupole / Dispersion"), QString());
        AddOptSpherical(g, QStringLiteral("pq_shielding"), model::QtPiQuadrupoleGroup(s).shielding(a), QStringLiteral("Å⁻⁵"));
        AddOptSpherical(g, QStringLiteral("disp_shielding"), model::QtDispersionGroup(s).shielding(a), QStringLiteral("Å⁻⁶"));
    }

    // ── Bond anisotropy (McConnell) ──
    {
        auto* g = AddKV(root, QStringLiteral("Bond anisotropy (McConnell)"), QString());
        model::QtMcConnellGroup mc(s);
        AddOptSpherical(g, QStringLiteral("mc_shielding"), mc.shielding(a), QStringLiteral("Å⁻³"));
        if (auto sc = mc.scalars(a)) {
            AddScalar(g, QStringLiteral("Σ C=O angular"), sc->co_sum);
            AddScalar(g, QStringLiteral("Σ C–N angular"), sc->cn_sum);
            AddScalar(g, QStringLiteral("Σ sidechain"), sc->sidechain_sum);
            AddScalar(g, QStringLiteral("Σ aromatic"), sc->aromatic_sum);
            AddKV(g, QStringLiteral("nearest C=O"),
                  sc->hasNearestCO() ? FmtDouble(sc->nearest_CO_dist) + QStringLiteral(" Å") : QStringLiteral("none"));
            AddKV(g, QStringLiteral("nearest C–N"),
                  sc->hasNearestCN() ? FmtDouble(sc->nearest_CN_dist) + QStringLiteral(" Å") : QStringLiteral("none"));
        }
    }

    // ── Electrostatics (Coulomb / APBS / AIMNet2 EFG) ──
    {
        auto* g = AddKV(root, QStringLiteral("Electrostatics"), QString());
        model::QtCoulombGroup coul(s);
        model::QtApbsGroup apbs(s);
        AddOptSpherical(g, QStringLiteral("coulomb_shielding"), coul.shielding(a), QStringLiteral("V/Å²"));
        AddOptVec3(g, QStringLiteral("coulomb_E"), coul.E(a), QStringLiteral("V/Å"));
        AddOptVec3(g, QStringLiteral("apbs_E"), apbs.E(a), QStringLiteral("V/Å"));
        AddOptEfg(g, QStringLiteral("apbs_efg"), apbs.efg(a), QStringLiteral("V/Å²"));
        AddOptEfg(g, QStringLiteral("aimnet2_efg"), model::QtAimnet2Group(s).efg(a), QStringLiteral("V/Å²"));
    }

    // ── H-bond (kernel form) ──
    {
        auto* g = AddKV(root, QStringLiteral("H-bond"), QString());
        model::QtHBondGroup hb(s);
        AddOptSpherical(g, QStringLiteral("hbond_shielding"), hb.shielding(a), QStringLiteral("Å⁻³"));
        if (auto sc = hb.scalars(a)) {
            AddScalar(g, QStringLiteral("nearest dist"), sc->nearest_dist, QStringLiteral("Å"));
            AddScalar(g, QStringLiteral("1/r³"), sc->inv_d3);
            AddScalar(g, QStringLiteral("count ≤ 3.5 Å"), sc->count_3_5A);
        }
    }

    // ── SASA ──
    {
        auto* g = AddKV(root, QStringLiteral("SASA"), QString());
        model::QtSasaGroup sasa(s);
        AddOptScalar(g, QStringLiteral("atom_sasa"), sasa.sasa(a), QStringLiteral("Å²"));
        AddOptVec3(g, QStringLiteral("surface normal"), sasa.normal(a));
    }

    // ── Water environment ──
    {
        auto* g = AddKV(root, QStringLiteral("Water"), QString());
        model::QtWaterFieldGroup wf(s);
        AddOptVec3(g, QStringLiteral("water_efield"), wf.efield(a), QStringLiteral("V/Å"));
        AddOptEfg(g, QStringLiteral("water_efg"), wf.efg(a), QStringLiteral("V/Å²"));
        if (auto wc = wf.shellCounts(a))
            AddKV(g, QStringLiteral("shell counts (1st/2nd)"), QStringLiteral("%1 / %2").arg(wc->nFirst).arg(wc->nSecond));
        if (auto h = model::QtHydrationGroup(s).shell(a)) {
            AddScalar(g, QStringLiteral("half-shell asymmetry"), h->halfShellAsymmetry);
            AddScalar(g, QStringLiteral("mean water dipole cos"), h->meanWaterDipoleCos);
            AddKV(g, QStringLiteral("nearest ion"),
                  h->hasNearestIon()
                      ? QStringLiteral("%1 Å, q=%2 e").arg(FmtDouble(h->nearestIonDist), FmtDouble(h->nearestIonCharge))
                      : QStringLiteral("none in cutoff"));
        }
        if (auto p = model::QtWaterPolarizationGroup(s).polarization(a)) {
            AddScalar(g, QStringLiteral("dipole alignment"), p->alignment);
            AddScalar(g, QStringLiteral("dipole coherence"), p->coherence);
        }
    }

    // ── Charges ──
    {
        auto* g = AddKV(root, QStringLiteral("Charges"), QString());
        model::QtAimnet2Group aim(s);
        model::QtEeqGroup eeq(s);
        AddOptScalar(g, QStringLiteral("AIMNet2 (Hirshfeld)"), aim.charge(a), QStringLiteral("e"));
        AddOptScalar(g, QStringLiteral("EEQ"), eeq.charge(a), QStringLiteral("e"));
        AddOptScalar(g, QStringLiteral("EEQ coord. number"), eeq.coordinationNumber(a));
        AddOptScalar(g, QStringLiteral("|charge-response grad|"), aim.chargeResponseGradientNorm(a), QStringLiteral("e²/Å"));
    }

    // ── DSSP (per-residue, broadcast to the atom) ──
    if (auto bb = model::QtDsspGroup(s).backbone(a)) {
        auto* g = AddKV(root, QStringLiteral("DSSP (secondary structure)"), QString());
        AddScalar(g, QStringLiteral("φ (neg-IUPAC)"), bb->phi, QStringLiteral("rad"));
        AddScalar(g, QStringLiteral("ψ (neg-IUPAC)"), bb->psi, QStringLiteral("rad"));
        AddScalar(g, QStringLiteral("residue SASA"), bb->sasa, QStringLiteral("Å²"));
        if (auto ss = model::QtDsspGroup(s).ss8(a))
            AddKV(g, QStringLiteral("SS8 class (ordinal)"), QString::number(static_cast<int>(ss->dominant())));
    }

    // ── Planar geometry ──
    {
        model::QtPlanarGeometryGroup pg(s);
        const int resIdx = protein_->atom(a).residueIndex;
        auto* g = AddKV(root, QStringLiteral("Planar geometry"), QString());
        AddOptScalar(g, QStringLiteral("pyramidalization"), pg.pyramidalization(a), QStringLiteral("Å"));
        if (resIdx >= 0) {
            const std::size_t r = static_cast<std::size_t>(resIdx);
            AddOptScalar(g, QStringLiteral("ω (peptide)"), pg.omegaActual(r), QStringLiteral("rad"));
            AddOptScalar(g, QStringLiteral("ω deviation"), pg.omegaDeviation(r), QStringLiteral("rad"));
            AddOptBool(g, QStringLiteral("X→Pro context"), pg.omegaIsXpro(r));
        }
    }

    // ── Energy (per-atom bonded share + whole-frame GROMACS) ──
    {
        if (auto be = model::QtBondedGroup(s).energy(a)) {
            auto* g = AddKV(root, QStringLiteral("Bonded energy (per-atom share)"), QString());
            AddScalar(g, QStringLiteral("total"), be->total, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("bond"), be->bond, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("angle"), be->angle, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("proper dih"), be->proper, QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("improper dih"), be->improper, QStringLiteral("kJ/mol"));
        }
        if (auto ge = model::QtGromacsGroup(s).energy()) {
            auto* g = AddKV(root, QStringLiteral("Frame energy (GROMACS)"), QString());
            AddScalar(g, QStringLiteral("potential"), ge->potential(), QStringLiteral("kJ/mol"));
            AddScalar(g, QStringLiteral("temperature"), ge->temperature(), QStringLiteral("K"));
            AddScalar(g, QStringLiteral("pressure"), ge->pressure(), QStringLiteral("bar"));
        }
    }

    // ── MOPAC (PM7+MOZYME; FullFat --mopac runs only) ──
    {
        model::QtMopacCoreGroup mopac(s);
        if (auto ms = mopac.scalars(a)) {
            auto* g = AddKV(root, QStringLiteral("MOPAC (PM7+MOZYME)"), QString());
            AddScalar(g, QStringLiteral("charge"), ms->charge, QStringLiteral("e"));
            AddScalar(g, QStringLiteral("s population"), ms->sPop, QStringLiteral("e"));
            AddScalar(g, QStringLiteral("p population"), ms->pPop, QStringLiteral("e"));
            AddScalar(g, QStringLiteral("Wiberg valency"), ms->valency);
            AddOptSpherical(g, QStringLiteral("mopac_coulomb_shielding"), model::QtMopacCoulombGroup(s).shielding(a), QStringLiteral("V/Å²"));
            AddOptSpherical(g, QStringLiteral("mopac_mc_shielding"), model::QtMopacMcConnellGroup(s).shielding(a), QStringLiteral("Å⁻³"));
            if (auto mg = mopac.global())
                AddScalar(g, QStringLiteral("ΔHf (frame)"), mg->heatOfFormation, QStringLiteral("kcal/mol"));
        }
    }

    // ── DFT / ProCS15 reference shielding (single-pose --orca + tripeptide / Larsen) ──
    {
        model::QtOrcaGroup orca(s);
        if (auto tot = orca.total(a)) {
            auto* g = AddKV(root, QStringLiteral("DFT reference (ORCA)"), QString());
            AddOptSpherical(g, QStringLiteral("σ total"), tot, QStringLiteral("ppm"));
            AddOptSpherical(g, QStringLiteral("σ diamagnetic"), orca.diamagnetic(a), QStringLiteral("ppm"));
            AddOptSpherical(g, QStringLiteral("σ paramagnetic"), orca.paramagnetic(a), QStringLiteral("ppm"));
        }
        model::QtTripeptideGroup trip(s);
        if (trip.hasBackboneMatch(a)) {
            auto* g = AddKV(root, QStringLiteral("Tripeptide reference (ProCS15)"), QString());
            AddOptSpherical(g, QStringLiteral("backbone σ"), trip.backboneShielding(a), QStringLiteral("ppm"));
            AddOptScalar(g, QStringLiteral("match distance"), trip.backboneMatchDistance(a), QStringLiteral("Å"));
            AddOptSpherical(g, QStringLiteral("neighbor σ (i±1)"), trip.neighborShielding(a), QStringLiteral("ppm"));
        }
        model::QtLarsenHBondGroup larsen(s);
        if (larsen.hasContribution(a)) {
            auto* g = AddKV(root, QStringLiteral("Larsen H-bond (ProCS15 grid)"), QString());
            AddOptSpherical(g, QStringLiteral("Δσ total"), larsen.shielding(a), QStringLiteral("ppm"));
            AddOptScalar(g, QStringLiteral("water term"), larsen.waterTerm(a), QStringLiteral("ppm"));
            AddOptInt(g, QStringLiteral("H-bond pair count"), larsen.count(a));
        }
    }

    tree_->expandToDepth(1);
    qCDebug(cDock).noquote() << "rebuilt | atom=" << a << "| frame=" << t << "| snapshot= resident";
}

}  // namespace h5reader::app

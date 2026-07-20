// metric_taxonomy_tests -- headless QtTest for the PURE mechanism->concept->form
// classification (model/MetricTaxonomy.h) that organizes the flat catalog into a
// navigable, hypothesis-first tree. Header-only / no I/O, so links Qt6::Test only
// (mirrors csa_math_tests). Verifies the family->group map, the cross-family
// concept overrides (AIMNet2 / explicit-water straddle electrostatic-EFG vs
// charges/solvation), the charge-model sub-tag, fold-suffix form derivation, that
// namespace dots are NOT folded, and that no descriptor lands in Other.

#include "model/MetricTaxonomy.h"

#include <QObject>
#include <QtTest>

using namespace h5reader::model;

namespace {
SignalDescriptor desc(const char* id, const char* concept, const char* family,
                      SignalSourceKind src) {
    SignalDescriptor d;
    d.id = QString::fromLatin1(id);
    d.conceptKey = QString::fromLatin1(concept);
    d.family = QString::fromLatin1(family);
    d.sourceKind = src;
    d.label = QString::fromLatin1(concept);
    return d;
}
MetricClass cls(const char* id, const char* concept, const char* family, SignalSourceKind src) {
    return ClassifyMetric(desc(id, concept, family, src));
}
}  // namespace

class MetricTaxonomyTests : public QObject {
    Q_OBJECT
private slots:
    void shieldingKernelsGroup();
    void electrostaticByChargeModel();
    void crossFamilyOverrides();
    void referenceAndForms();
    void experimentalShieldingMlIsNotReference();
    void namespaceDotsNotFolded();
    void formSuffixesFold();
    void formCollapseGathersForms();
};

void MetricTaxonomyTests::shieldingKernelsGroup() {
    QVERIFY(cls("npy:bs_shielding", "bs_shielding", "biot_savart", SignalSourceKind::FrameNpySnapshot).group == MetricGroup::RingCurrent);
    QVERIFY(cls("h5:hm", "hm_shielding", "haigh_mallion", SignalSourceKind::DenseH5Trajectory).group == MetricGroup::RingCurrent);
    QVERIFY(cls("npy:mc", "mc_shielding", "mcconnell", SignalSourceKind::FrameNpySnapshot).group == MetricGroup::BondAnisotropy);
    QVERIFY(cls("npy:l", "larsen_hbond_1pHB_shielding", "larsen_hbond", SignalSourceKind::FrameNpySnapshot).group == MetricGroup::HBond);
    // roles: the four kernels are HYPOTHESES.
    QVERIFY(RoleForGroup(MetricGroup::RingCurrent) == MetricRole::Hypothesis);
    QVERIFY(RoleForGroup(MetricGroup::Electrostatic) == MetricRole::Hypothesis);
    QVERIFY(RoleForGroup(MetricGroup::DftReference) == MetricRole::Reference);
}

void MetricTaxonomyTests::electrostaticByChargeModel() {
    const auto a = cls("npy:coulomb_E", "coulomb_E", "coulomb", SignalSourceKind::FrameNpySnapshot);
    QVERIFY(a.group == MetricGroup::Electrostatic);
    QCOMPARE(a.chargeModel, QStringLiteral("Coulomb (point charge)"));
    QCOMPARE(cls("npy:mc_E", "mopac_coulomb_E", "mopac_coulomb", SignalSourceKind::FrameNpySnapshot).chargeModel, QStringLiteral("MOPAC charge"));
    QCOMPARE(cls("npy:apbs", "apbs_E", "apbs", SignalSourceKind::FrameNpySnapshot).chargeModel, QStringLiteral("APBS continuum"));
}

void MetricTaxonomyTests::crossFamilyOverrides() {
    // AIMNet2 straddles: EFG -> Electrostatic, charges -> Charges.
    const auto efg = cls("npy:a_efg", "aimnet2_efg", "aimnet2", SignalSourceKind::FrameNpySnapshot);
    QVERIFY(efg.group == MetricGroup::Electrostatic);
    QCOMPARE(efg.chargeModel, QStringLiteral("AIMNet2"));
    QVERIFY(cls("npy:a_q", "aimnet2_charges", "aimnet2", SignalSourceKind::FrameNpySnapshot).group == MetricGroup::Charges);
    // explicit water straddles: EFG -> Electrostatic, shell counts -> Solvation.
    QVERIFY(cls("npy:w_efg", "water_efg", "water_field", SignalSourceKind::FrameNpySnapshot).group == MetricGroup::Electrostatic);
    QVERIFY(cls("npy:w_sc", "water_shell_counts", "water_field", SignalSourceKind::FrameNpySnapshot).group == MetricGroup::Solvation);
    // ring geometry sits in the identity family but is ring-current.
    QVERIFY(cls("npy:rc", "ring_contributions", "identity", SignalSourceKind::FrameNpySnapshot).group == MetricGroup::RingCurrent);
    QVERIFY(cls("npy:mq", "mopac_charges", "mopac_core", SignalSourceKind::FrameNpySnapshot).group == MetricGroup::Charges);
    QVERIFY(cls("npy:d", "delta_shielding", "mutation_delta", SignalSourceKind::FrameNpySnapshot).group == MetricGroup::Mutation);
}

void MetricTaxonomyTests::referenceAndForms() {
    const auto o = cls("orca_dft:total", "orca_total", "orca", SignalSourceKind::OrcaDftFrame);
    QVERIFY(o.group == MetricGroup::DftReference);
    QVERIFY(o.form == MetricForm::Reference);
    QVERIFY(cls("topology:res", "topology.residues", "topology", SignalSourceKind::Topology).form == MetricForm::Spine);
    QVERIFY(cls("geo:d", "geometry.distance", "geometry", SignalSourceKind::DerivedGeometry).form == MetricForm::Derived);
}

void MetricTaxonomyTests::experimentalShieldingMlIsNotReference() {
    const auto e = cls("ml:t2", "experimental_shielding_ml.t2",
                       "experimental_shielding_ml",
                       SignalSourceKind::ExperimentalShieldingMl);
    QVERIFY(e.group == MetricGroup::Experimental);
    QVERIFY(e.form == MetricForm::Derived);
    QVERIFY(RoleForGroup(e.group) == MetricRole::Experimental);
}

void MetricTaxonomyTests::namespaceDotsNotFolded() {
    // topology.residues / kernel_dynamics.acf are NAMESPACE dots -> base unchanged.
    const auto t = cls("topology:res", "topology.residues", "topology", SignalSourceKind::Topology);
    QCOMPARE(t.baseConcept, QStringLiteral("topology.residues"));
    QVERIFY(t.group == MetricGroup::Scaffold);
    const auto k = cls("h5:kd", "kernel_dynamics.acf", "kernel_dynamics", SignalSourceKind::DenseH5Trajectory);
    QCOMPARE(k.baseConcept, QStringLiteral("kernel_dynamics.acf"));
    QVERIFY(k.group == MetricGroup::Dynamics);
}

void MetricTaxonomyTests::formSuffixesFold() {
    // .stats / .autocorrelation / .transition fold into the base concept.
    const auto s = cls("h5:bsw", "bs_shielding.stats", "biot_savart", SignalSourceKind::DenseH5Trajectory);
    QCOMPARE(s.baseConcept, QStringLiteral("bs_shielding"));
    QVERIFY(s.form == MetricForm::Rollup);
    QVERIFY(s.group == MetricGroup::RingCurrent);
    const auto ac = cls("h5:bsa", "bs_shielding.autocorrelation", "biot_savart", SignalSourceKind::DenseH5Trajectory);
    QCOMPARE(ac.baseConcept, QStringLiteral("bs_shielding"));
    QVERIFY(ac.form == MetricForm::Dynamics);
    const auto tr = cls("h5:dt", "dssp_ss8.transition", "dssp", SignalSourceKind::DenseH5Trajectory);
    QCOMPARE(tr.baseConcept, QStringLiteral("dssp_ss8"));
    QVERIFY(tr.form == MetricForm::Transition);
}

void MetricTaxonomyTests::formCollapseGathersForms() {
    // The SAME quantity across snapshot/series/rollup/dynamics becomes ONE concept
    // with several FORMS (related items shown together, not a flattened wall).
    QVector<SignalDescriptor> cat{
        desc("npy:bs_shielding", "bs_shielding", "biot_savart", SignalSourceKind::FrameNpySnapshot),
        desc("h5:bs_shielding_time_series", "bs_shielding", "biot_savart", SignalSourceKind::DenseH5Trajectory),
        desc("bs_welford", "bs_shielding.stats", "biot_savart", SignalSourceKind::DenseH5Trajectory),
        desc("bs_t0_autocorrelation", "bs_shielding.autocorrelation", "biot_savart", SignalSourceKind::DenseH5Trajectory),
        desc("orca_dft:total", "orca_total", "orca", SignalSourceKind::OrcaDftFrame),
    };
    const QVector<MetricGroupNode> tree = GroupCatalog(cat);
    // RingCurrent group present with ONE concept (bs_shielding) carrying 4 forms.
    bool foundRC = false, foundRef = false;
    for (const MetricGroupNode& g : tree) {
        if (g.group == MetricGroup::RingCurrent) {
            foundRC = true;
            QCOMPARE(g.concepts.size(), qsizetype(1));
            QCOMPARE(g.concepts[0].baseConcept, QStringLiteral("bs_shielding"));
            QCOMPARE(g.concepts[0].forms.size(), qsizetype(4));
        }
        if (g.group == MetricGroup::DftReference) foundRef = true;
    }
    QVERIFY(foundRC);
    QVERIFY(foundRef);
}

QTEST_APPLESS_MAIN(MetricTaxonomyTests)
#include "metric_taxonomy_tests.moc"

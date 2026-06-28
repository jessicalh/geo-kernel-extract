// MetricTaxonomy -- the mechanism -> concept -> form classification that turns
// the flat ~188-descriptor catalog into a navigable, hypothesis-first tree.
//
// Pure functions over SignalDescriptor (DashboardSignal.h); no Qt widgets / I/O,
// so the SAME classification is shared by the picker UI, REST (/catalog/tree),
// and the headless test -- the catalog can never disagree with itself about how
// it groups.
//
// Three axes (see notes/metrics-tree-2026-06-22):
//   * GROUP   -- the mechanism / hypothesis a quantity belongs to. The four
//                shielding-contribution kernels (ring-current, bond-anisotropy,
//                electrostatic, H-bond) are the competing HYPOTHESES; DFT/ProCS15
//                is the REFERENCE they are judged against; charges/solvation/
//                structure are conditioning INPUTS; plus dynamics and scaffold.
//   * CONCEPT -- the base quantity with the form suffix stripped, so the SAME
//                quantity is ONE concept across its snapshot/series/rollup forms.
//   * FORM    -- which representation a descriptor is (snapshot = per-frame npy,
//                series = h5 time-series, rollup = .stats welford, dynamics =
//                .autocorrelation, transition = .transition, ...).
// Electrostatic carries a chargeModel sub-tag (Coulomb / MOPAC / APBS / explicit
// water / AIMNet2) -- the user navigates electrostatic shielding BY charge model.

#pragma once

#include "DashboardSignal.h"

#include <QHash>
#include <QLatin1String>
#include <QString>
#include <QVector>

#include <cstdint>

namespace h5reader::model {

enum class MetricGroup : std::uint8_t {
    RingCurrent,     // A1  aromatic ring-current (Biot-Savart / Haigh-Mallion / pi-quad / disp)
    BondAnisotropy,  // A2  McConnell point-dipole
    Electrostatic,   // A3  E-field / EFG, sub-grouped by charge model
    HBond,           // A4  Larsen / ProCS15 H-bond
    DftReference,    // B   ORCA + tripeptide/ProCS15 -- the targets
    Experimental,    // B'  local experimental model estimates
    Charges,         // C   MOPAC / EEQ / AIMNet2 charges + electronic structure
    Solvation,       // C   SASA / hydration / water environment
    Structure,       // C   secondary structure + geometry
    Dynamics,        // D   autocorrelation / relaxation / order parameters
    Scaffold,        // E   identity / topology / energies / QC
    Mutation,        // E   WT/mut/delta (mutant-pair runs only)
    Other,
};

enum class MetricRole : std::uint8_t {
    Hypothesis,  // the competing shielding-contribution kernels (A1-A4)
    Reference,   // the DFT / empirical truth they are compared against (B)
    Experimental,  // local experimental estimates, not QM references
    Input,       // conditioning features (C)
    Dynamics,    // time-evolution summaries (D)
    Scaffold,    // identity / bookkeeping / perturbation (E)
};

enum class MetricForm : std::uint8_t {
    Snapshot,    // per-frame NPY
    Series,      // dense-H5 time series
    Rollup,      // .stats welford mean/var
    Dynamics,    // .autocorrelation (a dynamics form of a base concept)
    Transition,  // .transition event matrix
    Reference,   // live .out DFT (orca_dft:*)
    Spine,       // startup-loaded topology sidecar
    Derived,     // computed on demand (geometry)
    Other,
};

inline const char* ToString(MetricGroup g) {
    switch (g) {
        case MetricGroup::RingCurrent:    return "ring_current";
        case MetricGroup::BondAnisotropy: return "bond_anisotropy";
        case MetricGroup::Electrostatic:  return "electrostatic";
        case MetricGroup::HBond:          return "hbond";
        case MetricGroup::DftReference:   return "dft_reference";
        case MetricGroup::Experimental:   return "experimental";
        case MetricGroup::Charges:        return "charges";
        case MetricGroup::Solvation:      return "solvation";
        case MetricGroup::Structure:      return "structure";
        case MetricGroup::Dynamics:       return "dynamics";
        case MetricGroup::Scaffold:       return "scaffold";
        case MetricGroup::Mutation:       return "mutation";
        case MetricGroup::Other:          return "other";
    }
    return "other";
}

inline const char* ToString(MetricRole r) {
    switch (r) {
        case MetricRole::Hypothesis: return "hypothesis";
        case MetricRole::Reference:  return "reference";
        case MetricRole::Experimental: return "experimental";
        case MetricRole::Input:      return "input";
        case MetricRole::Dynamics:   return "dynamics";
        case MetricRole::Scaffold:   return "scaffold";
    }
    return "scaffold";
}

inline const char* ToString(MetricForm f) {
    switch (f) {
        case MetricForm::Snapshot:   return "snapshot";
        case MetricForm::Series:     return "series";
        case MetricForm::Rollup:     return "rollup";
        case MetricForm::Dynamics:   return "dynamics";
        case MetricForm::Transition: return "transition";
        case MetricForm::Reference:  return "reference";
        case MetricForm::Spine:      return "spine";
        case MetricForm::Derived:    return "derived";
        case MetricForm::Other:      return "other";
    }
    return "other";
}

inline MetricRole RoleForGroup(MetricGroup g) {
    switch (g) {
        case MetricGroup::RingCurrent:
        case MetricGroup::BondAnisotropy:
        case MetricGroup::Electrostatic:
        case MetricGroup::HBond:
            return MetricRole::Hypothesis;
        case MetricGroup::DftReference:
            return MetricRole::Reference;
        case MetricGroup::Experimental:
            return MetricRole::Experimental;
        case MetricGroup::Charges:
        case MetricGroup::Solvation:
        case MetricGroup::Structure:
            return MetricRole::Input;
        case MetricGroup::Dynamics:
            return MetricRole::Dynamics;
        case MetricGroup::Scaffold:
        case MetricGroup::Mutation:
        case MetricGroup::Other:
            return MetricRole::Scaffold;
    }
    return MetricRole::Scaffold;
}

// Strip a fold-form suffix (.stats / .autocorrelation / .transition) from a
// concept key, returning the base concept and (via *form) the form that suffix
// implies. NAMESPACE dots (topology.residues, kernel_dynamics.acf, geometry.x)
// are NOT fold suffixes -- they are left intact, so those stay their own concept.
inline QString StripFormSuffix(const QString& conceptKey, MetricForm* form) {
    struct Fold { QLatin1String suffix; MetricForm form; };
    const Fold folds[] = {
        {QLatin1String(".stats"), MetricForm::Rollup},
        {QLatin1String(".autocorrelation"), MetricForm::Dynamics},
        {QLatin1String(".transition"), MetricForm::Transition},
    };
    for (const Fold& f : folds) {
        if (conceptKey.endsWith(f.suffix)) {
            if (form) *form = f.form;
            return conceptKey.chopped(f.suffix.size());
        }
    }
    return conceptKey;
}

// The electrostatic charge model (Electrostatic group only); empty otherwise.
inline QString ChargeModelFor(const QString& family) {
    if (family == QLatin1String("coulomb"))       return QStringLiteral("Coulomb (point charge)");
    if (family == QLatin1String("mopac_coulomb")) return QStringLiteral("MOPAC charge");
    if (family == QLatin1String("apbs"))          return QStringLiteral("APBS continuum");
    if (family == QLatin1String("water_field"))   return QStringLiteral("Explicit water");
    if (family == QLatin1String("aimnet2"))       return QStringLiteral("AIMNet2");
    return {};
}

// Mechanism group for a descriptor, from its family + base concept. Concept-level
// overrides come FIRST (the cross-family cases: AIMNet2 and explicit-water each
// straddle electrostatic-EFG and a non-electrostatic role; ring geometry sits in
// the identity family; the MOPAC/FF reconciliation is a QC scalar, not a charge).
inline MetricGroup GroupForDescriptor(const QString& family, const QString& baseConcept) {
    if (baseConcept.startsWith(QLatin1String("aimnet2_efg")))    return MetricGroup::Electrostatic;
    if (baseConcept.startsWith(QLatin1String("water_efg"))
        || baseConcept.startsWith(QLatin1String("water_efield"))) return MetricGroup::Electrostatic;
    if (baseConcept == QLatin1String("water_shell_counts"))      return MetricGroup::Solvation;
    if (baseConcept.startsWith(QLatin1String("ring_contributions"))
        || baseConcept.startsWith(QLatin1String("ring_geometry"))) return MetricGroup::RingCurrent;
    if (baseConcept == QLatin1String("mopac_vs_ff14sb_reconciliation")) return MetricGroup::Scaffold;
    if (baseConcept == QLatin1String("positions") || baseConcept == QLatin1String("element")
        || baseConcept.startsWith(QLatin1String("residue_"))
        || baseConcept == QLatin1String("atoms_category_info")) return MetricGroup::Scaffold;

    if (family == QLatin1String("biot_savart") || family == QLatin1String("haigh_mallion")
        || family == QLatin1String("pi_quadrupole") || family == QLatin1String("dispersion")
        || family == QLatin1String("ring_current")) return MetricGroup::RingCurrent;
    if (family == QLatin1String("mcconnell") || family == QLatin1String("mopac_mcconnell"))
        return MetricGroup::BondAnisotropy;
    if (family == QLatin1String("coulomb") || family == QLatin1String("mopac_coulomb")
        || family == QLatin1String("apbs")) return MetricGroup::Electrostatic;
    if (family == QLatin1String("larsen_hbond") || family == QLatin1String("hbond"))
        return MetricGroup::HBond;
    if (family == QLatin1String("orca") || family == QLatin1String("tripeptide"))
        return MetricGroup::DftReference;
    if (family == QLatin1String("experimental_shielding_ml"))
        return MetricGroup::Experimental;
    if (family == QLatin1String("eeq") || family == QLatin1String("mopac_core")
        || family == QLatin1String("aimnet2")) return MetricGroup::Charges;
    if (family == QLatin1String("sasa") || family == QLatin1String("hydration")
        || family == QLatin1String("water_polarization") || family == QLatin1String("water_field"))
        return MetricGroup::Solvation;
    if (family == QLatin1String("dssp") || family == QLatin1String("planar_geometry")
        || family == QLatin1String("geometry") || family == QLatin1String("j_coupling"))
        return MetricGroup::Structure;
    if (family == QLatin1String("kernel_dynamics") || family == QLatin1String("kernel_coherence")
        || family == QLatin1String("dihedral_autocorrelation")
        || family == QLatin1String("reorientational_dynamics")
        || family == QLatin1String("lipari_szabo") || family == QLatin1String("rmsd"))
        return MetricGroup::Dynamics;
    if (family == QLatin1String("identity") || family == QLatin1String("topology")
        || family == QLatin1String("selections") || family == QLatin1String("gromacs")
        || family == QLatin1String("bonded")) return MetricGroup::Scaffold;
    if (family == QLatin1String("mutation_delta")) return MetricGroup::Mutation;
    return MetricGroup::Other;
}

inline MetricForm FormForDescriptor(const SignalDescriptor& d, MetricForm suffixForm,
                                    bool hadSuffix) {
    if (hadSuffix) return suffixForm;
    switch (d.sourceKind) {
        case SignalSourceKind::FrameNpySnapshot:  return MetricForm::Snapshot;
        case SignalSourceKind::DenseH5Trajectory: return MetricForm::Series;
        case SignalSourceKind::OrcaDftFrame:      return MetricForm::Reference;
        case SignalSourceKind::Topology:          return MetricForm::Spine;
        case SignalSourceKind::DerivedGeometry:   return MetricForm::Derived;
        case SignalSourceKind::SelectionEvents:   return MetricForm::Other;
        case SignalSourceKind::ExperimentalShieldingMl: return MetricForm::Derived;
    }
    return MetricForm::Other;
}

struct MetricClass {
    MetricGroup group = MetricGroup::Other;
    QString baseConcept;     // form-suffix stripped
    MetricForm form = MetricForm::Other;
    QString chargeModel;     // Electrostatic only; empty otherwise
};

inline MetricClass ClassifyMetric(const SignalDescriptor& d) {
    MetricClass out;
    MetricForm suffixForm = MetricForm::Other;
    out.baseConcept = StripFormSuffix(d.conceptKey, &suffixForm);
    const bool hadSuffix = (out.baseConcept != d.conceptKey);
    out.group = GroupForDescriptor(d.family, out.baseConcept);
    out.form = FormForDescriptor(d, suffixForm, hadSuffix);
    if (out.group == MetricGroup::Electrostatic)
        out.chargeModel = ChargeModelFor(d.family);
    return out;
}

// ---- The grouped tree (mechanism -> concept -> forms) -------------------------

struct MetricFormEntry {
    MetricForm form = MetricForm::Other;
    QString descriptorId;
};

struct MetricConceptNode {
    QString baseConcept;
    QString label;        // human label (from the first descriptor seen)
    QString chargeModel;  // Electrostatic only
    MetricGroup group = MetricGroup::Other;
    QVector<MetricFormEntry> forms;
};

struct MetricGroupNode {
    MetricGroup group = MetricGroup::Other;
    MetricRole role = MetricRole::Scaffold;
    QVector<MetricConceptNode> concepts;
};

// Collapse a flat descriptor list into the mechanism -> concept -> form tree.
// The 188 catalog rows fold to ~35 base concepts in ~11 groups: the same quantity
// across snapshot/series/rollup/... becomes ONE concept with several forms.
inline QVector<MetricGroupNode> GroupCatalog(const QVector<SignalDescriptor>& descriptors) {
    QHash<QString, MetricConceptNode> byConcept;  // key = "<group>|<base>" (avoid cross-group base collisions)
    QVector<QString> conceptOrder;
    for (const SignalDescriptor& d : descriptors) {
        const MetricClass c = ClassifyMetric(d);
        const QString key = QString::number(static_cast<int>(c.group)) + QLatin1Char('|') + c.baseConcept;
        auto it = byConcept.find(key);
        if (it == byConcept.end()) {
            MetricConceptNode node;
            node.baseConcept = c.baseConcept;
            node.label = d.label.isEmpty() ? c.baseConcept : d.label;
            node.chargeModel = c.chargeModel;
            node.group = c.group;
            it = byConcept.insert(key, node);
            conceptOrder.push_back(key);
        }
        it->forms.push_back({c.form, d.id});
    }

    QVector<MetricGroupNode> groups;
    auto groupIndex = [&](MetricGroup g) -> int {
        for (int i = 0; i < groups.size(); ++i)
            if (groups[i].group == g) return i;
        groups.push_back({g, RoleForGroup(g), {}});
        return static_cast<int>(groups.size()) - 1;
    };
    for (const QString& key : conceptOrder) {
        const MetricConceptNode& node = byConcept[key];
        groups[groupIndex(node.group)].concepts.push_back(node);
    }
    return groups;
}

}  // namespace h5reader::model

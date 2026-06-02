// main_extract — headless `h5reader-rediscover` CLI. QCoreApplication only
// (no widgets, no VTK). Loads one 1P9J calcset into the all-frames-resident
// RunData and runs the requested relationships, each writing its two CSV row
// kinds + sidecar NPYs.
//
// Usage:
//   h5reader_extract --run <calcset_dir_or_LGS> --out <dir> [--case ring|mc|all]
//                    [--engine composed|procedural]
//
// ring_current / mcconnell run through the COMPOSED functional engine by
// default (Layer-1 verbs + Layer-2 curried closures + the Layer-3
// RunRelationship loop; SURFACE_DESIGN.md), or the PROCEDURAL reference cells
// under `--engine procedural` (the oracle the composed path is diffed against).
// The driver builds a plain vector of RunnableCases and runs them in a loop —
// no scheduler, no dependency graph (DESIGN.md). Progress flows through the
// StructuredLogger (UDP 9997 + stderr).

#include "ExtractionSupport.h"
#include "AllAtomEquivariant.h"
#include "AllAtomEquivariantSink.h"
#include "Aimnet2Feature.h"
#include "Aimnet2FeatureSink.h"
#include "BroadBackbone.h"
#include "BroadBackboneSink.h"
#include "BuckinghamEfield.h"
#include "BuckinghamEfieldSink.h"
#include "Catalog.h"
#include "ChargeDipoleNeighborhood.h"
#include "ComposedRelationships.h"
#include "EfgFeature.h"
#include "EfgFeatureSink.h"
#include "McConnellNeighborhood.h"
#include "OutputManifest.h"
#include "RecordSink.h"
#include "RediscoveryExtraction.h"
#include "Relationship.h"
#include "RelationshipEngine.h"
#include "ResidentIndexes.h"
#include "RingCurrentNeighborhood.h"
#include "RunData.h"

#include "../diagnostics/StructuredLogger.h"

#include <QCommandLineParser>
#include <QCoreApplication>
#include <QLoggingCategory>
#include <QStringList>

#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
Q_LOGGING_CATEGORY(cMain, "h5reader.rediscover.main")

QString relationshipKindName(h5reader::rediscover::RelationshipKind kind) {
    switch (kind) {
    case h5reader::rediscover::RelationshipKind::SourceSum:
        return QStringLiteral("source_sum");
    case h5reader::rediscover::RelationshipKind::PerAtomFeature:
        return QStringLiteral("per_atom_feature");
    }
    return QStringLiteral("unknown");
}

bool isFailLoudStub(const QString& which) {
    return which == QStringLiteral("charge_quadrupole")
           || which == QStringLiteral("larsen_hbond");
}

QString stubMessage(const QString& which, const QString& chargeSource,
                    const h5reader::rediscover::Catalog& catalog) {
    if (which == QStringLiteral("charge_quadrupole")) {
        if (chargeSource == QStringLiteral("mopac")
            && !catalog.has(h5reader::rediscover::ArrayId::MopacCharge)) {
            return QStringLiteral("%1 requires charge_source=mopac, but per-frame MOPAC charges are absent")
                .arg(which);
        }
        if (chargeSource == QStringLiteral("ff14sb")
            && !catalog.has(h5reader::rediscover::ArrayId::Ff14sbCharge)) {
            return QStringLiteral("%1 requires FF14SB charges, but topol.top charges were not loaded")
                .arg(which);
        }
        if (chargeSource == QStringLiteral("aimnet2")
            && !catalog.has(h5reader::rediscover::ArrayId::Aimnet2Charge)) {
            return QStringLiteral("%1 requires AIMNet2 charges, but aimnet2_charge is absent")
                .arg(which);
        }
        return QStringLiteral("%1 is registered but not runnable yet: charge multipole reducer/schema is fail-loud, no zero output")
            .arg(which);
    }
    if (which == QStringLiteral("efg") && !catalog.has(h5reader::rediscover::ArrayId::ApbsEfg))
        return QStringLiteral("efg requires APBS EFG, but apbs_efg is absent");
    if (which == QStringLiteral("buckingham_efield") && !catalog.has(h5reader::rediscover::ArrayId::ApbsEfield))
        return QStringLiteral("buckingham_efield requires APBS E-field, but apbs_efield is absent");
    return QStringLiteral("%1 is registered as a fail-loud stub; data/decision is not ready, no zeros emitted")
        .arg(which);
}

h5reader::rediscover::ArrayId chargeArrayForSource(const QString& chargeSource) {
    if (chargeSource == QStringLiteral("ff14sb")) return h5reader::rediscover::ArrayId::Ff14sbCharge;
    if (chargeSource == QStringLiteral("aimnet2")) return h5reader::rediscover::ArrayId::Aimnet2Charge;
    if (chargeSource == QStringLiteral("mopac")) return h5reader::rediscover::ArrayId::MopacCharge;
    return h5reader::rediscover::ArrayId::MopacCharge;
}

QString chargeSourceMissingMessage(const QString& relationship, const QString& chargeSource) {
    if (chargeSource == QStringLiteral("ff14sb"))
        return QStringLiteral("%1 requires FF14SB charges, but topol.top charges were not loaded")
            .arg(relationship);
    if (chargeSource == QStringLiteral("aimnet2"))
        return QStringLiteral("%1 requires AIMNet2 Hirshfeld charges, but aimnet2_charge is absent")
            .arg(relationship);
    if (chargeSource == QStringLiteral("mopac"))
        return QStringLiteral("%1 requires charge_source=mopac, but per-frame MOPAC charges are absent")
            .arg(relationship);
    return QStringLiteral("%1 received unknown charge_source=%2").arg(relationship, chargeSource);
}

bool isAimnet2FeatureCase(const QString& which) {
    return which == QStringLiteral("aimnet2_features")
           || which == QStringLiteral("aimnet2_embedding")
           || which == QStringLiteral("charge_response_gradient")
           || which == QStringLiteral("crg");
}

QString aimnet2FeatureMissingMessage(const h5reader::rediscover::Catalog& catalog) {
    QStringList missing;
    if (!catalog.has(h5reader::rediscover::ArrayId::Aimnet2Charge))
        missing << QStringLiteral("aimnet2_charge");
    if (!catalog.has(h5reader::rediscover::ArrayId::Aimnet2ChargeRespScalar))
        missing << QStringLiteral("aimnet2_charge_response_gradient_scalar");
    if (!catalog.has(h5reader::rediscover::ArrayId::Aimnet2ChargeRespVector))
        missing << QStringLiteral("aimnet2_charge_response_gradient_vector");
    if (!catalog.has(h5reader::rediscover::ArrayId::Aimnet2Embedding))
        missing << QStringLiteral("aimnet2_embedding");
    if (missing.isEmpty()) return {};
    return QStringLiteral("aimnet2_features requires AIMNet2 charge, CRG scalar/vector, and 256-d embedding; missing: %1")
        .arg(missing.join(QStringLiteral(", ")));
}
}

int main(int argc, char** argv) {
    QCoreApplication app(argc, argv);
    QCoreApplication::setApplicationName(QStringLiteral("h5reader-rediscover"));
    QCoreApplication::setApplicationVersion(QStringLiteral(H5READER_VERSION));

    // Structured logger first, before anything else (reader discipline).
    h5reader::diagnostics::StructuredLogger::Install();

    QCommandLineParser parser;
    parser.setApplicationDescription(
        QStringLiteral("Rediscover substrate extractor: per-(atom,frame) "
                       "feature/target CSVs for ring-current and McConnell strata."));
    parser.addHelpOption();
    parser.addVersionOption();
    QCommandLineOption runOpt(QStringLiteral("run"),
                              QStringLiteral("Calcset directory holding the single .LGS, or the .LGS path."),
                              QStringLiteral("path"));
    QCommandLineOption outOpt(QStringLiteral("out"),
                              QStringLiteral("Output directory for the CSV files."),
                              QStringLiteral("dir"));
    QCommandLineOption caseOpt(QStringLiteral("case"),
                               QStringLiteral("Which extraction(s): ring_current | mcconnell | charge_dipole | broad_backbone | all_atom_equivariant | efg | buckingham_efield | aimnet2_features | ring | mc | all, or a registered fail-loud stub."),
                               QStringLiteral("case"), QStringLiteral("all"));
    // McConnell source-discovery cutoff (Å). Surfaced + recorded per the
    // substrate conventions' no-hidden-cutoffs rule; 10.0 Å matches the
    // producer McConnell bond-anisotropy default.
    QCommandLineOption mcCutoffOpt(QStringLiteral("mc-cutoff"),
                                   QStringLiteral("McConnell bond-discovery cutoff in Angstrom (default 10.0)."),
                                   QStringLiteral("angstrom"), QStringLiteral("10.0"));
    QCommandLineOption chargeSourceOpt(QStringLiteral("charge-source"),
                                       QStringLiteral("Charge source for charge cases: ff14sb | aimnet2 | mopac."),
                                       QStringLiteral("source"), QStringLiteral("ff14sb"));
    QCommandLineOption chargeCutoffOpt(QStringLiteral("charge-cutoff"),
                                       QStringLiteral("Charge-site source cutoff in Angstrom for charge_dipole / broad_backbone charge field. Sweep 6/10/12 (6 truncates the long-range 1/r^2 field). Default 6.0."),
                                       QStringLiteral("angstrom"), QStringLiteral("6.0"));
    // broad_backbone ring/bond source-discovery cutoffs (Å), recorded per the
    // no-hidden-cutoffs rule. Ring 8.0 = conventions' aromatic-neighbourhood
    // recommendation; bond 10.0 matches the producer McConnell reach.
    QCommandLineOption ringCutoffOpt(QStringLiteral("ring-cutoff"),
                                     QStringLiteral("broad_backbone ring-centre cutoff in Angstrom (default 8.0)."),
                                     QStringLiteral("angstrom"), QStringLiteral("8.0"));
    QCommandLineOption bondCutoffOpt(QStringLiteral("bond-cutoff"),
                                     QStringLiteral("broad_backbone anisotropic-bond cutoff in Angstrom (default 10.0)."),
                                     QStringLiteral("angstrom"), QStringLiteral("10.0"));
    QCommandLineOption mcNearFieldRatioOpt(
        QStringLiteral("mc-near-field-ratio"),
        QStringLiteral("broad_backbone McConnell near-field exclusion ratio (default 0.5)."),
        QStringLiteral("ratio"), QStringLiteral("0.5"));
    // Which traversal stands ring_current / mcconnell up: the composed
    // functional API (Layer 1 verbs + Layer 2 curried closures + the Layer 3
    // engine — the default, what we validate) or the original procedural cells
    // (the reference oracle the composed path is checked against). SURFACE_DESIGN.
    QCommandLineOption engineOpt(QStringLiteral("engine"),
                                 QStringLiteral("ring/mc traversal: composed (functional API, default) | procedural (reference oracle cells)."),
                                 QStringLiteral("engine"), QStringLiteral("composed"));
    parser.addOption(runOpt);
    parser.addOption(outOpt);
    parser.addOption(caseOpt);
    parser.addOption(mcCutoffOpt);
    parser.addOption(chargeSourceOpt);
    parser.addOption(chargeCutoffOpt);
    parser.addOption(ringCutoffOpt);
    parser.addOption(bondCutoffOpt);
    parser.addOption(mcNearFieldRatioOpt);
    parser.addOption(engineOpt);
    parser.process(app);

    if (!parser.isSet(runOpt) || !parser.isSet(outOpt)) {
        qCCritical(cMain) << "both --run and --out are required";
        parser.showHelp(2);
    }
    const QString runPath = parser.value(runOpt);
    const QString outDir = parser.value(outOpt);
    const QString which = parser.value(caseOpt);
    const QString chargeSource = parser.value(chargeSourceOpt);
    const QString engine = parser.value(engineOpt);
    if (engine != QStringLiteral("composed") && engine != QStringLiteral("procedural")) {
        qCCritical(cMain).noquote() << "invalid --engine" << engine
                                    << "(expected composed|procedural)";
        return 2;
    }
    if (chargeSource != QStringLiteral("ff14sb") && chargeSource != QStringLiteral("aimnet2")
        && chargeSource != QStringLiteral("mopac")) {
        qCCritical(cMain).noquote() << "invalid --charge-source" << chargeSource
                                    << "(expected ff14sb|aimnet2|mopac)";
        return 2;
    }
    bool cutoffOk = false;
    const double mcCutoff = parser.value(mcCutoffOpt).toDouble(&cutoffOk);
    if (!cutoffOk || !(mcCutoff > 0.0)) {
        qCCritical(cMain).noquote() << "invalid --mc-cutoff" << parser.value(mcCutoffOpt);
        return 2;
    }
    bool chargeCutoffOk = false;
    const double chargeCutoff = parser.value(chargeCutoffOpt).toDouble(&chargeCutoffOk);
    if (!chargeCutoffOk || !(chargeCutoff > 0.0)) {
        qCCritical(cMain).noquote() << "invalid --charge-cutoff" << parser.value(chargeCutoffOpt);
        return 2;
    }
    bool ringCutoffOk = false;
    const double ringCutoff = parser.value(ringCutoffOpt).toDouble(&ringCutoffOk);
    if (!ringCutoffOk || !(ringCutoff > 0.0)) {
        qCCritical(cMain).noquote() << "invalid --ring-cutoff" << parser.value(ringCutoffOpt);
        return 2;
    }
    bool bondCutoffOk = false;
    const double bondCutoff = parser.value(bondCutoffOpt).toDouble(&bondCutoffOk);
    if (!bondCutoffOk || !(bondCutoff > 0.0)) {
        qCCritical(cMain).noquote() << "invalid --bond-cutoff" << parser.value(bondCutoffOpt);
        return 2;
    }
    bool mcNearFieldRatioOk = false;
    const double mcNearFieldRatio = parser.value(mcNearFieldRatioOpt).toDouble(&mcNearFieldRatioOk);
    if (!mcNearFieldRatioOk || !(mcNearFieldRatio >= 0.0)) {
        qCCritical(cMain).noquote()
            << "invalid --mc-near-field-ratio" << parser.value(mcNearFieldRatioOpt);
        return 2;
    }

    qCInfo(cMain).noquote() << "loading run" << runPath;
    QString err;
    auto run = h5reader::rediscover::RunLoader::Load(runPath, &err);
    if (!run) {
        qCCritical(cMain).noquote() << "load failed:" << err;
        return 1;
    }

    qCInfo(cMain) << "building rediscover catalog and resident indexes";
    h5reader::rediscover::Catalog catalog(*run);
    h5reader::rediscover::ResidentIndexes indexes = h5reader::rediscover::BuildResidentIndexes(*run);
    const h5reader::rediscover::Body body{*run, indexes, catalog};

    // Resolve the T2 Cartesian-frame caveat: is the ORCA-input geometry (the
    // frame the DFT tensors live in) the same orientation as the H5 positions?
    const auto align = h5reader::rediscover::CheckDftFrameAlignment(*run);
    qCInfo(cMain).noquote()
        << "DFT frame check | frames=" << align.n_frames << "| atoms=" << align.n_atoms_used
        << "| rotation mean=" << align.mean_angle_deg << "deg max=" << align.max_angle_deg
        << "deg | RMSD mean=" << align.mean_rmsd_A << "A max=" << align.max_rmsd_A << "A | T2 components"
        << (align.max_angle_deg < 1.0 ? "FRAME-ALIGNED (comparable as emitted)"
                                      : "ROTATED (need tensor rotation; T0/|T2| invariants safe)");

    if (isFailLoudStub(which)) {
        qCCritical(cMain).noquote() << "ValidateScenario failed:" << stubMessage(which, chargeSource, catalog);
        return 2;
    }

    if (which == QStringLiteral("charge_dipole")) {
        const h5reader::rediscover::ArrayId chargeArray = chargeArrayForSource(chargeSource);
        if (!catalog.has(chargeArray)) {
            qCCritical(cMain).noquote()
                << "ValidateScenario failed:" << chargeSourceMissingMessage(QStringLiteral("charge_dipole"), chargeSource);
            return 2;
        }
    }

    // ── efg — focused per_atom_feature sibling carrier. This is #29's second
    // non-source_sum data point after broad_backbone; keep it direct here and
    // do not widen RunRelationship / RecordSink in this spike. ──────────────
    if (which == QStringLiteral("efg")) {
        if (!catalog.has(h5reader::rediscover::ArrayId::ApbsEfg)) {
            qCCritical(cMain).noquote()
                << "ValidateScenario failed: efg requires APBS EFG, but apbs_efg is absent";
            return 2;
        }
        h5reader::rediscover::EfgFeatureSink sink(outDir, QStringLiteral("efg"));
        if (!sink.Ok()) {
            qCCritical(cMain).noquote() << "efg sink open failed";
            return 3;
        }
        h5reader::rediscover::EfgFeatureStats stats;
        try {
            stats = h5reader::rediscover::RunEfgPerAtomFeature(body, sink);
        } catch (const std::exception& e) {
            qCCritical(cMain).noquote() << "efg failed:" << e.what();
            return 1;
        }
        const bool committed = sink.Commit();
        qCInfo(cMain).noquote() << "efg | rows=" << sink.rowsWritten()
                                << "| dft_present=" << stats.dft_present
                                << "| apbs_efg_present=" << stats.apbs_efg_present
                                << "| finite_efg=" << stats.finite_efg
                                << "| committed=" << committed;
        if (!committed) return 4;
        std::vector<h5reader::rediscover::OutputEntry> outputs = {
            {QStringLiteral("efg"), QStringLiteral("per_atom_feature"), QString(),
             QStringLiteral("efg_aggregated.csv"), sink.sidecarFiles(), stats.rows,
             0, sink.rowsWritten()}};
        QString manifestErr;
        if (!h5reader::rediscover::WriteOutputManifest(outDir, outputs, align, 0, &manifestErr)) {
            qCCritical(cMain).noquote() << "manifest write failed:" << manifestErr;
            return 4;
        }
        return 0;
    }

    // ── buckingham_efield — APBS solvated-PB E-field projected in the local
    // backbone frame, T0-only fit target. T1/T2 target payloads are emitted for
    // audit completeness; T1 remains convention-unverified and is not fitted. ─
    if (which == QStringLiteral("buckingham_efield")) {
        if (!catalog.has(h5reader::rediscover::ArrayId::ApbsEfield)) {
            qCCritical(cMain).noquote()
                << "ValidateScenario failed: buckingham_efield requires APBS E-field, but apbs_efield is absent";
            return 2;
        }
        h5reader::rediscover::BuckinghamEfieldSink sink(outDir, QStringLiteral("buckingham_efield"));
        if (!sink.Ok()) {
            qCCritical(cMain).noquote() << "buckingham_efield sink open failed";
            return 3;
        }
        h5reader::rediscover::BuckinghamEfieldStats stats;
        try {
            stats = h5reader::rediscover::RunBuckinghamEfieldPerAtomFeature(body, sink);
        } catch (const std::exception& e) {
            qCCritical(cMain).noquote() << "buckingham_efield failed:" << e.what();
            return 1;
        }
        const bool committed = sink.Commit();
        qCInfo(cMain).noquote() << "buckingham_efield | rows=" << sink.rowsWritten()
                                << "| dft_present=" << stats.dft_present
                                << "| frame_valid=" << stats.frame_valid
                                << "| apbs_efield_present=" << stats.apbs_efield_present
                                << "| finite_efield=" << stats.finite_efield
                                << "| committed=" << committed;
        if (!committed) return 4;
        std::vector<h5reader::rediscover::OutputEntry> outputs = {
            {QStringLiteral("buckingham_efield"), QStringLiteral("per_atom_feature"), QString(),
             QStringLiteral("buckingham_efield_aggregated.csv"), sink.sidecarFiles(), stats.rows,
             0, sink.rowsWritten()}};
        QString manifestErr;
        if (!h5reader::rediscover::WriteOutputManifest(outDir, outputs, align, 0, &manifestErr)) {
            qCCritical(cMain).noquote() << "manifest write failed:" << manifestErr;
            return 4;
        }
        return 0;
    }

    // -- AIMNet2 features: charge source, charge-response-gradient (CRG, not
    // polarizability), and 256-d embedding. One row-aligned per-atom substrate.
    if (isAimnet2FeatureCase(which)) {
        const QString missing = aimnet2FeatureMissingMessage(catalog);
        if (!missing.isEmpty()) {
            qCCritical(cMain).noquote() << "ValidateScenario failed:" << missing;
            return 2;
        }
        h5reader::rediscover::Aimnet2FeatureSink sink(outDir, QStringLiteral("aimnet2_features"), 256);
        if (!sink.Ok()) {
            qCCritical(cMain).noquote() << "aimnet2_features sink open failed";
            return 3;
        }
        h5reader::rediscover::Aimnet2FeatureStats stats;
        try {
            stats = h5reader::rediscover::RunAimnet2PerAtomFeature(body, sink);
        } catch (const std::exception& e) {
            qCCritical(cMain).noquote() << "aimnet2_features failed:" << e.what();
            return 1;
        }
        const bool committed = sink.Commit();
        qCInfo(cMain).noquote() << "aimnet2_features | rows=" << sink.rowsWritten()
                                << "| dft_present=" << stats.dft_present
                                << "| charge_present=" << stats.charge_present
                                << "| crg_present=" << stats.crg_present
                                << "| embedding_present=" << stats.embedding_present
                                << "| embedding_dims=" << stats.embedding_dims
                                << "| committed=" << committed;
        if (!committed) return 4;
        std::vector<h5reader::rediscover::OutputEntry> outputs = {
            {QStringLiteral("aimnet2_features"), QStringLiteral("per_atom_feature"), QString(),
             QStringLiteral("aimnet2_features_aggregated.csv"), sink.sidecarFiles(), stats.rows,
             0, sink.rowsWritten()}};
        QString manifestErr;
        if (!h5reader::rediscover::WriteOutputManifest(outDir, outputs, align, 0, &manifestErr)) {
            qCCritical(cMain).noquote() << "manifest write failed:" << manifestErr;
            return 4;
        }
        return 0;
    }

    // -- all_atom_equivariant: corrected e3nn substrate. Every atom, all ring
    // types, all producer bond categories, FF14SB/AIMNet2 charge-site rows when
    // available, and per-target APBS/AIMNet2 source rows. Everything is emitted
    // in the H5/ORCA-aligned molecular lab frame; no local frame is imposed.
    if (which == QStringLiteral("all_atom_equivariant")
        || which == QStringLiteral("all_atom_equiv")) {
        h5reader::rediscover::AllAtomEquivariantConfig cfg;
        cfg.ring_cutoff_A = ringCutoff;
        cfg.bond_cutoff_A = bondCutoff;
        cfg.charge_cutoff_A = chargeCutoff;
        cfg.mc_near_field_ratio = mcNearFieldRatio;
        h5reader::rediscover::AllAtomEquivariantSink sink(outDir,
                                                          QStringLiteral("all_atom_equivariant"),
                                                          256);
        if (!sink.Ok()) {
            qCCritical(cMain).noquote() << "all_atom_equivariant sink open failed";
            return 3;
        }
        h5reader::rediscover::AllAtomEquivariantStats stats;
        try {
            stats = h5reader::rediscover::RunAllAtomEquivariantEmit(body, sink, cfg);
        } catch (const std::exception& e) {
            qCCritical(cMain).noquote() << "all_atom_equivariant failed:" << e.what();
            return 1;
        }
        const bool committed = sink.Commit();
        qCInfo(cMain).noquote()
            << "all_atom_equivariant | atoms=" << stats.atom_count
            << "| dft_rows=" << stats.dft_rows
            << "| target_rows=" << sink.targetRowsWritten()
            << "| dft_present=" << stats.dft_present
            << "| source_rows=" << sink.sourceRowsWritten()
            << "| ring_rows=" << stats.ring_rows
            << "| bond_rows=" << stats.bond_rows
            << "| charge_ff14sb_rows=" << stats.charge_ff14sb_rows
            << "| charge_aimnet2_rows=" << stats.charge_aimnet2_rows
            << "| apbs_efield_rows=" << stats.apbs_efield_rows
            << "| apbs_efg_rows=" << stats.apbs_efg_rows
            << "| aimnet2_atom_rows=" << stats.aimnet2_atom_rows
            << "| committed=" << committed;
        if (!committed) return 4;
        std::vector<h5reader::rediscover::OutputEntry> outputs = {
            {QStringLiteral("all_atom_equivariant"), QStringLiteral("source_sum"),
             QStringLiteral("all_atom_equivariant_sources.csv"),
             QStringLiteral("all_atom_equivariant_targets.csv"), sink.sidecarFiles(),
             stats.target_rows, sink.sourceRowsWritten(), sink.targetRowsWritten()}};
        QString manifestErr;
        if (!h5reader::rediscover::WriteOutputManifest(outDir, outputs, align, 0, &manifestErr)) {
            qCCritical(cMain).noquote() << "manifest write failed:" << manifestErr;
            return 4;
        }
        return 0;
    }

    // ── broad_backbone — the composed heterogeneous relationship (its own
    // two-kind carrier; not a RunnableCase/RecordSink). ValidateScenario, run,
    // commit, manifest, return. ────────────────────────────────────────────
    if (which == QStringLiteral("broad_backbone")) {
        const h5reader::rediscover::ArrayId chargeArray = chargeArrayForSource(chargeSource);
        if (!catalog.has(chargeArray)) {
            qCCritical(cMain).noquote()
                << "ValidateScenario failed:" << chargeSourceMissingMessage(QStringLiteral("broad_backbone"), chargeSource);
            return 2;
        }
        const h5reader::rediscover::BroadRelationship brel =
            h5reader::rediscover::MakeBroadBackboneRelationship(ringCutoff, bondCutoff, chargeCutoff,
                                                                chargeSource, /*exclude_residue=*/true,
                                                                mcNearFieldRatio);
        h5reader::rediscover::BroadBackboneSink sink(outDir, QStringLiteral("broad_backbone"));
        if (!sink.Ok()) {
            qCCritical(cMain).noquote() << "broad_backbone sink open failed";
            return 3;
        }
        std::size_t cases = 0;
        try {
            cases = h5reader::rediscover::RunBroadBackbone(brel, body, sink);
        } catch (const std::exception& e) {
            qCCritical(cMain).noquote() << "broad_backbone failed:" << e.what();
            return 1;
        }
        const bool committed = sink.Commit();
        qCInfo(cMain).noquote() << "broad_backbone | cases=" << cases
                                << "| source_rows=" << sink.sourceRowsWritten()
                                << "| agg_rows=" << sink.aggregatedRowsWritten()
                                << "| committed=" << committed;
        if (!committed) return 4;
        std::vector<h5reader::rediscover::OutputEntry> outputs = {
            {QStringLiteral("broad_backbone"), QStringLiteral("source_sum"),
             QStringLiteral("broad_backbone_sources.csv"),
             QStringLiteral("broad_backbone_aggregated.csv"), sink.sidecarFiles(), cases,
             sink.sourceRowsWritten(), sink.aggregatedRowsWritten()}};
        QString manifestErr;
        if (!h5reader::rediscover::WriteOutputManifest(outDir, outputs, align, 0, &manifestErr)) {
            qCCritical(cMain).noquote() << "manifest write failed:" << manifestErr;
            return 4;
        }
        return 0;
    }

    // A runnable case: its name, its declared schema, and a runner that walks
    // the body into the sink. ring_current / mcconnell run either through the
    // COMPOSED functional engine (default — RunRelationship over a Relationship
    // built from curried closures) or the PROCEDURAL reference cells
    // (--engine procedural), so the two can be diffed for the oracle gate.
    // charge_dipole stays on its procedural cell (out of scope for the API).
    struct RunnableCase {
        QString name;
        h5reader::rediscover::FeatureSchema schema;
        std::function<std::size_t(const h5reader::rediscover::Body&,
                                  h5reader::rediscover::RecordSink&)>
            run;
    };
    const bool composed = (engine == QStringLiteral("composed"));
    std::vector<RunnableCase> cases_to_run;

    const bool wantRing = which == QStringLiteral("ring")
                          || which == QStringLiteral("ring_current") || which == QStringLiteral("all");
    const bool wantMc = which == QStringLiteral("mc") || which == QStringLiteral("mcconnell")
                        || which == QStringLiteral("all");

    if (wantRing) {
        if (composed) {
            auto rel = std::make_shared<h5reader::rediscover::Relationship>(
                h5reader::rediscover::MakeRingCurrentRelationship());
            cases_to_run.push_back(
                {rel->name, rel->schema,
                 [rel](const h5reader::rediscover::Body& b, h5reader::rediscover::RecordSink& s) {
                     return h5reader::rediscover::RunRelationship(*rel, b, s);
                 }});
        } else {
            auto cell = std::make_shared<h5reader::rediscover::RingCurrentNeighborhood>();
            cases_to_run.push_back(
                {cell->name(), cell->schema(),
                 [cell](const h5reader::rediscover::Body& b, h5reader::rediscover::RecordSink& s) {
                     return cell->extract(b, s);
                 }});
        }
    }
    if (wantMc) {
        if (composed) {
            auto rel = std::make_shared<h5reader::rediscover::Relationship>(
                h5reader::rediscover::MakeMcConnellRelationship(mcCutoff));
            cases_to_run.push_back(
                {rel->name, rel->schema,
                 [rel](const h5reader::rediscover::Body& b, h5reader::rediscover::RecordSink& s) {
                     return h5reader::rediscover::RunRelationship(*rel, b, s);
                 }});
        } else {
            auto cell = std::make_shared<h5reader::rediscover::McConnellNeighborhood>();
            cell->cutoff_A = mcCutoff;
            cases_to_run.push_back(
                {cell->name(), cell->schema(),
                 [cell](const h5reader::rediscover::Body& b, h5reader::rediscover::RecordSink& s) {
                     return cell->extract(b, s);
                 }});
        }
    }
    if (which == QStringLiteral("charge_dipole")) {
        auto cell = std::make_shared<h5reader::rediscover::ChargeDipoleNeighborhood>();
        cell->charge_source = chargeSource;
        cell->cutoff_A = chargeCutoff;
        cases_to_run.push_back(
            {cell->name(), cell->schema(),
             [cell](const h5reader::rediscover::Body& b, h5reader::rediscover::RecordSink& s) {
                 return cell->extract(b, s);
             }});
    }
    if (cases_to_run.empty()) {
        qCCritical(cMain).noquote() << "unknown --case" << which
                                    << "(expected ring|mc|charge_dipole|broad_backbone|all_atom_equivariant|efg|buckingham_efield|aimnet2_features|all)";
        return 2;
    }
    qCInfo(cMain).noquote() << "engine =" << engine << "| cases =" << cases_to_run.size();

    int rc = 0;
    std::vector<h5reader::rediscover::OutputEntry> outputs;
    for (const auto& rcase : cases_to_run) {
        const h5reader::rediscover::FeatureSchema schema = rcase.schema;
        h5reader::rediscover::RecordSink sink(outDir, schema);
        if (!sink.Ok()) {
            qCCritical(cMain).noquote() << "sink open failed for" << rcase.name;
            rc = 3;
            continue;
        }
        std::size_t cases = 0;
        try {
            cases = rcase.run(body, sink);
        } catch (const std::exception& e) {
            qCCritical(cMain).noquote() << rcase.name << "failed:" << e.what();
            return 1;
        }
        const bool committed = sink.Commit();
        qCInfo(cMain).noquote() << rcase.name << "| cases=" << cases
                                << "| source_rows=" << sink.sourceRowsWritten()
                                << "| agg_rows=" << sink.aggregatedRowsWritten()
                                << "| committed=" << committed;
        if (!committed) rc = 4;
        if (committed) {
            outputs.push_back({rcase.name,
                               relationshipKindName(schema.relationshipKind),
                               QStringLiteral("%1_sources.csv").arg(rcase.name),
                               QStringLiteral("%1_aggregated.csv").arg(rcase.name),
                               sink.sidecarFiles(),
                               cases,
                               sink.sourceRowsWritten(),
                               sink.aggregatedRowsWritten()});
        }
    }
    if (rc == 0) {
        QString manifestErr;
        if (!h5reader::rediscover::WriteOutputManifest(outDir, outputs, align, rc, &manifestErr)) {
            qCCritical(cMain).noquote() << "manifest write failed:" << manifestErr;
            return 4;
        }
    }
    return rc;
}

// main_extract — headless `h5reader-rediscover` CLI. QCoreApplication only
// (no widgets, no VTK). Loads one 1P9J calcset into the all-frames-resident
// RunData and runs the two extractions, each writing its two CSV row kinds.
//
// Usage:
//   h5reader_extract --run <calcset_dir_or_LGS> --out <dir> [--case ring|mc|all]
//
// The driver runs a plain vector<unique_ptr<RediscoveryExtraction>> loop — no
// scheduler, no dependency graph (DESIGN.md). Progress flows through the
// StructuredLogger (UDP 9997 + stderr).

#include "ExtractionSupport.h"
#include "Catalog.h"
#include "McConnellNeighborhood.h"
#include "OutputManifest.h"
#include "RecordSink.h"
#include "RediscoveryExtraction.h"
#include "ResidentIndexes.h"
#include "RingCurrentNeighborhood.h"
#include "RunData.h"

#include "../diagnostics/StructuredLogger.h"

#include <QCommandLineParser>
#include <QCoreApplication>
#include <QLoggingCategory>

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
    return which == QStringLiteral("buckingham_efield")
           || which == QStringLiteral("efg")
           || which == QStringLiteral("charge_dipole")
           || which == QStringLiteral("charge_quadrupole")
           || which == QStringLiteral("larsen_hbond")
           || which == QStringLiteral("charge_response_gradient")
           || which == QStringLiteral("crg")
           || which == QStringLiteral("aimnet2_embedding");
}

QString stubMessage(const QString& which, const QString& chargeSource,
                    const h5reader::rediscover::Catalog& catalog) {
    if (which == QStringLiteral("charge_dipole") || which == QStringLiteral("charge_quadrupole")) {
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
    if ((which == QStringLiteral("charge_response_gradient") || which == QStringLiteral("crg"))
        && (!catalog.has(h5reader::rediscover::ArrayId::Aimnet2ChargeRespScalar)
            || !catalog.has(h5reader::rediscover::ArrayId::Aimnet2ChargeRespVector)))
        return QStringLiteral("charge_response_gradient requires AIMNet2 CRG scalar/vector arrays, but they are absent");
    if (which == QStringLiteral("aimnet2_embedding")
        && !catalog.has(h5reader::rediscover::ArrayId::Aimnet2Embedding))
        return QStringLiteral("aimnet2_embedding requires AIMNet2 embedding, but aimnet2_embedding is absent");
    return QStringLiteral("%1 is registered as a fail-loud stub; data/decision is not ready, no zeros emitted")
        .arg(which);
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
                               QStringLiteral("Which extraction(s): ring_current | mcconnell | ring | mc | all, or a registered fail-loud stub."),
                               QStringLiteral("case"), QStringLiteral("all"));
    // McConnell source-discovery cutoff (Å). Surfaced + recorded per the
    // substrate conventions' no-hidden-cutoffs rule; 8.0 Å is the conventions'
    // aromatic-neighbourhood recommendation, reused as the anisotropy reach.
    QCommandLineOption mcCutoffOpt(QStringLiteral("mc-cutoff"),
                                   QStringLiteral("McConnell bond-discovery cutoff in Angstrom (default 8.0)."),
                                   QStringLiteral("angstrom"), QStringLiteral("8.0"));
    QCommandLineOption chargeSourceOpt(QStringLiteral("charge-source"),
                                       QStringLiteral("Charge source for charge cases: ff14sb | aimnet2 | mopac."),
                                       QStringLiteral("source"), QStringLiteral("ff14sb"));
    parser.addOption(runOpt);
    parser.addOption(outOpt);
    parser.addOption(caseOpt);
    parser.addOption(mcCutoffOpt);
    parser.addOption(chargeSourceOpt);
    parser.process(app);

    if (!parser.isSet(runOpt) || !parser.isSet(outOpt)) {
        qCCritical(cMain) << "both --run and --out are required";
        parser.showHelp(2);
    }
    const QString runPath = parser.value(runOpt);
    const QString outDir = parser.value(outOpt);
    const QString which = parser.value(caseOpt);
    const QString chargeSource = parser.value(chargeSourceOpt);
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

    std::vector<std::unique_ptr<h5reader::rediscover::RediscoveryExtraction>> extractions;
    if (which == QStringLiteral("ring") || which == QStringLiteral("ring_current") || which == QStringLiteral("all"))
        extractions.push_back(std::make_unique<h5reader::rediscover::RingCurrentNeighborhood>());
    if (which == QStringLiteral("mc") || which == QStringLiteral("mcconnell") || which == QStringLiteral("all")) {
        auto mc = std::make_unique<h5reader::rediscover::McConnellNeighborhood>();
        mc->cutoff_A = mcCutoff;
        extractions.push_back(std::move(mc));
    }
    if (extractions.empty()) {
        qCCritical(cMain).noquote() << "unknown --case" << which << "(expected ring|mc|all)";
        return 2;
    }

    int rc = 0;
    std::vector<h5reader::rediscover::OutputEntry> outputs;
    for (const auto& ex : extractions) {
        const h5reader::rediscover::FeatureSchema schema = ex->schema();
        h5reader::rediscover::RecordSink sink(outDir, schema);
        if (!sink.Ok()) {
            qCCritical(cMain).noquote() << "sink open failed for" << ex->name();
            rc = 3;
            continue;
        }
        std::size_t cases = 0;
        try {
            cases = ex->extract(body, sink);
        } catch (const std::exception& e) {
            qCCritical(cMain).noquote() << ex->name() << "failed:" << e.what();
            return 1;
        }
        const bool committed = sink.Commit();
        qCInfo(cMain).noquote() << ex->name() << "| cases=" << cases
                                << "| source_rows=" << sink.sourceRowsWritten()
                                << "| agg_rows=" << sink.aggregatedRowsWritten()
                                << "| committed=" << committed;
        if (!committed) rc = 4;
        if (committed) {
            outputs.push_back({ex->name(),
                               relationshipKindName(schema.relationshipKind),
                               QStringLiteral("%1_sources.csv").arg(ex->name()),
                               QStringLiteral("%1_aggregated.csv").arg(ex->name()),
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

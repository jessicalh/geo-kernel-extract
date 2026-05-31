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
#include "McConnellNeighborhood.h"
#include "RecordSink.h"
#include "RediscoveryExtraction.h"
#include "RingCurrentNeighborhood.h"
#include "RunData.h"

#include "../diagnostics/StructuredLogger.h"

#include <QCommandLineParser>
#include <QCoreApplication>
#include <QLoggingCategory>

#include <memory>
#include <vector>

namespace {
Q_LOGGING_CATEGORY(cMain, "h5reader.rediscover.main")
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
                               QStringLiteral("Which extraction(s): ring | mc | all (default all)."),
                               QStringLiteral("case"), QStringLiteral("all"));
    // McConnell source-discovery cutoff (Å). Surfaced + recorded per the
    // substrate conventions' no-hidden-cutoffs rule; 8.0 Å is the conventions'
    // aromatic-neighbourhood recommendation, reused as the anisotropy reach.
    QCommandLineOption mcCutoffOpt(QStringLiteral("mc-cutoff"),
                                   QStringLiteral("McConnell bond-discovery cutoff in Angstrom (default 8.0)."),
                                   QStringLiteral("angstrom"), QStringLiteral("8.0"));
    parser.addOption(runOpt);
    parser.addOption(outOpt);
    parser.addOption(caseOpt);
    parser.addOption(mcCutoffOpt);
    parser.process(app);

    if (!parser.isSet(runOpt) || !parser.isSet(outOpt)) {
        qCCritical(cMain) << "both --run and --out are required";
        parser.showHelp(2);
    }
    const QString runPath = parser.value(runOpt);
    const QString outDir = parser.value(outOpt);
    const QString which = parser.value(caseOpt);
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

    // Resolve the T2 Cartesian-frame caveat: is the ORCA-input geometry (the
    // frame the DFT tensors live in) the same orientation as the H5 positions?
    const auto align = h5reader::rediscover::CheckDftFrameAlignment(*run);
    qCInfo(cMain).noquote()
        << "DFT frame check | frames=" << align.n_frames << "| atoms=" << align.n_atoms_used
        << "| rotation mean=" << align.mean_angle_deg << "deg max=" << align.max_angle_deg
        << "deg | RMSD mean=" << align.mean_rmsd_A << "A max=" << align.max_rmsd_A << "A | T2 components"
        << (align.max_angle_deg < 1.0 ? "FRAME-ALIGNED (comparable as emitted)"
                                      : "ROTATED (need tensor rotation; T0/|T2| invariants safe)");

    std::vector<std::unique_ptr<h5reader::rediscover::RediscoveryExtraction>> extractions;
    if (which == QStringLiteral("ring") || which == QStringLiteral("all"))
        extractions.push_back(std::make_unique<h5reader::rediscover::RingCurrentNeighborhood>());
    if (which == QStringLiteral("mc") || which == QStringLiteral("all")) {
        auto mc = std::make_unique<h5reader::rediscover::McConnellNeighborhood>();
        mc->cutoff_A = mcCutoff;
        extractions.push_back(std::move(mc));
    }
    if (extractions.empty()) {
        qCCritical(cMain).noquote() << "unknown --case" << which << "(expected ring|mc|all)";
        return 2;
    }

    int rc = 0;
    for (const auto& ex : extractions) {
        h5reader::rediscover::RecordSink sink(outDir, ex->schema());
        if (!sink.Ok()) {
            qCCritical(cMain).noquote() << "sink open failed for" << ex->name();
            rc = 3;
            continue;
        }
        const std::size_t cases = ex->extract(*run, sink);
        const bool committed = sink.Commit();
        qCInfo(cMain).noquote() << ex->name() << "| cases=" << cases
                                << "| source_rows=" << sink.sourceRowsWritten()
                                << "| agg_rows=" << sink.aggregatedRowsWritten()
                                << "| committed=" << committed;
        if (!committed) rc = 4;
    }
    return rc;
}

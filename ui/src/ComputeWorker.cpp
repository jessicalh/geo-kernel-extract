// ComputeWorker: calls the library API, holds the result as shared_ptr<Protein>.
//
// Data flow:
//   BuildFromPdb → OperationRunner::Run → protein is fully const → emit
//
// The viewer reads the library objects directly. No copy into flat structs.

#include "ComputeWorker.h"

#include "Protein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "PdbFileReader.h"
#include "BuildResult.h"
#include "OperationRunner.h"
#include "OrcaRunLoader.h"
#include "BiotSavartResult.h"
#include "OperationLog.h"
#include "Atom.h"
#include "Residue.h"
#include "Ring.h"
#include "Bond.h"

#include <QElapsedTimer>
#include <cmath>

using namespace nmr;

ComputeWorker::ComputeWorker(QObject* parent) : QObject(parent) {}

void ComputeWorker::cancel() {
    cancelled_.store(true, std::memory_order_relaxed);
}

void ComputeWorker::computeAll(ViewerLoadRequest request) {
    cancelled_.store(false);
    QElapsedTimer timer;
    timer.start();

    ComputeResult result;
    result.proteinName = request.pdbPath;

    // ================================================================
    // Phase 1: Build protein (protonates, detects bonds, assigns charges)
    // ================================================================
    emit progress(0, 100, QStringLiteral("Building protein (protonation + bond perception)..."));

    OperationLog::Info(LogViewer, "ComputeWorker", "Phase1: loading " + request.pdbPath);

    BuildResult buildResult = BuildFromPdb(request.pdbPath);

    OperationLog::Info(LogViewer, "ComputeWorker",
        std::string("Phase1: Build returned ok=") + (buildResult.Ok() ? "true" : "false"));
    if (!buildResult.Ok()) {
        emit progress(0, 100, QString("Build failed: %1")
            .arg(QString::fromStdString(buildResult.error)));
        emit finished(result);
        return;
    }

    // Transfer ownership: unique_ptr → shared_ptr.
    // Protein is non-movable, non-copyable. shared_ptr takes the raw pointer.
    result.protein = std::shared_ptr<Protein>(buildResult.protein.release());

    auto& protein = *result.protein;
    auto& conf = protein.Conformation();

    OperationLog::Info(LogViewer, "ComputeWorker",
        "BuildFromPdb: " + std::to_string(protein.AtomCount()) + " atoms, " +
        std::to_string(protein.ResidueCount()) + " residues, " +
        std::to_string(protein.RingCount()) + " rings");

    if (cancelled_.load()) { emit finished(result); return; }

    // ================================================================
    // Phase 2: Run all calculators via OperationRunner
    // ================================================================
    emit progress(10, 100, QStringLiteral("Computing shielding tensors..."));

    RunOptions opts;
    if (buildResult.charges) {
        opts.charge_source = buildResult.charges.get();
        opts.net_charge = buildResult.net_charge;
    }

    OperationLog::Info(LogViewer, "ComputeWorker", "Phase2: calling OperationRunner::Run");
    auto runResult = OperationRunner::Run(conf, opts);
    OperationLog::Info(LogViewer, "ComputeWorker",
        std::string("Phase2: Run returned ok=") + (runResult.Ok() ? "true" : "false") +
        " attached=" + std::to_string(runResult.attached.size()));
    if (!runResult.Ok()) {
        emit progress(15, 100, QString("Run error: %1")
            .arg(QString::fromStdString(runResult.error)));
        // Continue — partial results are still useful
    }

    if (cancelled_.load()) { emit finished(result); return; }

    // ================================================================
    // Phase 3: T0 field grid computation (for isosurfaces around rings)
    // ================================================================
    if (conf.HasResult<BiotSavartResult>()) {
        const auto& bs = conf.Result<BiotSavartResult>();
        const int G = 20;
        const double extent = 7.0;
        int nRings = static_cast<int>(protein.RingCount());

        OperationLog::Info(LogViewer, "ComputeWorker::Phase3",
            "Computing T0 field grids: " + std::to_string(nRings) +
            " rings, " + std::to_string(G) + "^3 grid");

        for (size_t ri = 0; ri < protein.RingCount() && !cancelled_.load(); ++ri) {
            emit progress(60 + 25 * static_cast<int>(ri) / std::max(1, nRings), 100,
                QString("Field grid %1/%2 (%3)...")
                    .arg(ri + 1).arg(nRings)
                    .arg(QString::fromStdString(protein.RingAt(ri).TypeName())));

            const auto& geo = conf.ring_geometries[ri];
            ViewerFieldGrid grid;
            grid.ringCenter = geo.center;
            grid.ringType = protein.RingAt(ri).TypeName();

            grid.origin[0] = geo.center.x() - extent;
            grid.origin[1] = geo.center.y() - extent;
            grid.origin[2] = geo.center.z() - extent;
            grid.spacing[0] = grid.spacing[1] = grid.spacing[2] = 2.0 * extent / (G - 1);
            grid.dims[0] = grid.dims[1] = grid.dims[2] = G;

            int nPoints = G * G * G;
            grid.T0.resize(nPoints, 0.0);
            grid.bsT0.resize(nPoints, 0.0);

            for (int iz = 0; iz < G; ++iz) {
                for (int iy = 0; iy < G; ++iy) {
                    for (int ix = 0; ix < G; ++ix) {
                        if (cancelled_.load()) break;

                        Vec3 pt(grid.origin[0] + ix * grid.spacing[0],
                                grid.origin[1] + iy * grid.spacing[1],
                                grid.origin[2] + iz * grid.spacing[2]);

                        double dist = (pt - geo.center).norm();
                        if (dist > extent * 1.2 || dist < 0.5) continue;

                        auto stResult = bs.SampleShieldingAt(pt);
                        double bsT0val = stResult.T0;

                        if (!std::isfinite(bsT0val)) continue;

                        int idx = ix + iy * G + iz * G * G;
                        grid.bsT0[idx] = bsT0val;
                        grid.T0[idx] = bsT0val;
                    }
                }
            }
            result.fieldGrids.push_back(std::move(grid));
        }
    }

    // ================================================================
    // Phase 4: Butterfly field computation (optional, expensive)
    // ================================================================
    if (request.computeButterfly && protein.RingCount() > 0 &&
        conf.HasResult<BiotSavartResult>()) {
        emit progress(90, 100, QStringLiteral("Computing B-field grid..."));
        const auto& bs = conf.Result<BiotSavartResult>();

        const int G = 15;
        const double extent = 6.0;

        for (size_t ri = 0; ri < protein.RingCount() && !cancelled_.load(); ++ri) {
            const auto& geo = conf.ring_geometries[ri];
            ViewerButterflyData bfd;
            bfd.ringCenter = geo.center;
            bfd.ringNormal = geo.normal;
            bfd.ringRadius = geo.radius;
            bfd.ringType = protein.RingAt(ri).TypeName();
            bfd.gridDims[0] = bfd.gridDims[1] = bfd.gridDims[2] = G;

            Vec3 n = geo.normal;
            Vec3 arbitrary = (std::abs(n.x()) < 0.9) ? Vec3(1,0,0) : Vec3(0,1,0);
            Vec3 u = n.cross(arbitrary).normalized();
            Vec3 v = n.cross(u);

            bfd.positions.reserve(G * G * G);
            bfd.fields.reserve(G * G * G);

            for (int iz = 0; iz < G; ++iz) {
                for (int iy = 0; iy < G; ++iy) {
                    for (int ix = 0; ix < G; ++ix) {
                        double fx = -extent + 2.0 * extent * ix / (G - 1);
                        double fy = -extent + 2.0 * extent * iy / (G - 1);
                        double fz = -extent + 2.0 * extent * iz / (G - 1);
                        Vec3 pt = geo.center + fx * u + fy * v + fz * n;
                        Vec3 B = bs.SampleBFieldAt(pt);
                        bfd.positions.push_back(pt);
                        bfd.fields.push_back(B);
                    }
                }
            }
            result.butterflyFields.push_back(std::move(bfd));
        }
    }

    int N = static_cast<int>(conf.AtomCount());
    OperationLog::Info(LogViewer, "ComputeWorker",
        "Done: " + std::to_string(N) + " atoms, " +
        std::to_string(protein.RingCount()) + " rings, " +
        std::to_string(result.fieldGrids.size()) + " grids");

    emit progress(100, 100,
        QString("Done: %1 atoms, %2 rings, %3 ms")
            .arg(N).arg(protein.RingCount()).arg(timer.elapsed()));
    emit finished(result);
}

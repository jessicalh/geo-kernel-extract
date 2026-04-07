// ComputeWorker: dispatches on JobSpec.mode to the correct builder,
// then runs OperationRunner and computes viewer grids.
//
// Data flow:
//   JobSpec → builder → OperationRunner::Run → protein is fully const → emit
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
#include "GromacsEnsembleLoader.h"
#include "CalculatorConfig.h"
#include "BiotSavartResult.h"
#include "ConformationResult.h"
#include "OperationLog.h"
#include "Atom.h"
#include "Residue.h"
#include "Ring.h"
#include "Bond.h"

#include <QElapsedTimer>
#include <cmath>
#include <filesystem>

namespace fs = std::filesystem;
using namespace nmr;

ComputeWorker::ComputeWorker(QObject* parent) : QObject(parent) {}

void ComputeWorker::cancel() {
    cancelled_.store(true, std::memory_order_relaxed);
}

void ComputeWorker::computeAll(nmr::JobSpec spec) {
    cancelled_.store(false);
    QElapsedTimer timer;
    timer.start();

    ComputeResult result;

    // ================================================================
    // Phase 1: Build protein — dispatch on JobSpec mode
    // ================================================================
    emit progress(0, 100, QStringLiteral("Building protein..."));

    BuildResult buildResult;
    BuildResult alaBuild;  // only for mutant mode

    switch (spec.mode) {

    case JobMode::Pdb:
        OperationLog::Info(LogViewer, "ComputeWorker",
            "Phase1: BuildFromPdb " + spec.pdb_path + " pH=" + std::to_string(spec.pH));
        result.proteinName = spec.pdb_path;
        buildResult = BuildFromPdb(spec.pdb_path, spec.pH);
        break;

    case JobMode::ProtonatedPdb:
        OperationLog::Info(LogViewer, "ComputeWorker",
            "Phase1: BuildFromProtonatedPdb " + spec.pdb_path);
        result.proteinName = spec.pdb_path;
        buildResult = BuildFromProtonatedPdb(spec.pdb_path);
        break;

    case JobMode::Orca:
        OperationLog::Info(LogViewer, "ComputeWorker",
            "Phase1: BuildFromOrca xyz=" + spec.orca_files.xyz_path);
        result.proteinName = spec.orca_files.xyz_path;
        buildResult = BuildFromOrca(spec.orca_files);
        break;

    case JobMode::Mutant:
        OperationLog::Info(LogViewer, "ComputeWorker",
            "Phase1: BuildFromOrca mutant wt=" + spec.wt_files.xyz_path +
            " ala=" + spec.ala_files.xyz_path);
        result.proteinName = spec.wt_files.xyz_path;
        buildResult = BuildFromOrca(spec.wt_files);
        alaBuild = BuildFromOrca(spec.ala_files);
        if (!alaBuild.Ok()) {
            emit progress(0, 100, QString("ALA build failed: %1")
                .arg(QString::fromStdString(alaBuild.error)));
            emit finished(result);
            return;
        }
        break;

    case JobMode::Fleet:
        OperationLog::Info(LogViewer, "ComputeWorker",
            "Phase1: BuildFromGromacs poses=" + spec.fleet_paths.sampled_poses_dir);
        result.proteinName = spec.fleet_paths.sampled_poses_dir;
        buildResult = BuildFromGromacs(spec.fleet_paths);
        break;

    case JobMode::None:
        emit finished(result);
        return;
    }

    if (!buildResult.Ok()) {
        emit progress(0, 100, QString("Build failed: %1")
            .arg(QString::fromStdString(buildResult.error)));
        OperationLog::Error("ComputeWorker", buildResult.error);
        emit finished(result);
        return;
    }

    // Transfer ownership: unique_ptr → shared_ptr.
    result.protein = std::shared_ptr<Protein>(buildResult.protein.release());
    auto& protein = *result.protein;

    OperationLog::Info(LogViewer, "ComputeWorker",
        "Build: " + std::to_string(protein.AtomCount()) + " atoms, " +
        std::to_string(protein.ResidueCount()) + " residues, " +
        std::to_string(protein.RingCount()) + " rings, " +
        std::to_string(protein.ConformationCount()) + " conformations");

    if (cancelled_.load()) { emit finished(result); return; }

    // ================================================================
    // Phase 2: Run calculators
    // ================================================================
    emit progress(10, 100, QStringLiteral("Computing shielding tensors..."));

    RunOptions opts;
    if (buildResult.charges) {
        opts.charge_source = buildResult.charges.get();
        opts.net_charge = buildResult.net_charge;
    }

    if (spec.mode == JobMode::Mutant) {
        // Mutant: run both WT and ALA, compute delta
        auto& wt_conf = protein.Conformation();
        auto& ala_conf = alaBuild.protein->Conformation();

        RunOptions wt_opts = opts;
        if (!spec.wt_files.nmr_out_path.empty())
            wt_opts.orca_nmr_path = spec.wt_files.nmr_out_path;

        RunOptions ala_opts;
        ala_opts.charge_source = alaBuild.charges.get();
        ala_opts.net_charge = alaBuild.net_charge;
        if (!spec.ala_files.nmr_out_path.empty())
            ala_opts.orca_nmr_path = spec.ala_files.nmr_out_path;

        auto runResult = OperationRunner::RunMutantComparison(
            wt_conf, wt_opts, ala_conf, ala_opts);
        OperationLog::Info(LogViewer, "ComputeWorker",
            std::string("RunMutantComparison ok=") + (runResult.Ok() ? "true" : "false") +
            " attached=" + std::to_string(runResult.attached.size()));

    } else if (spec.mode == JobMode::Fleet) {
        // Fleet: run all conformations
        auto results = OperationRunner::RunEnsemble(protein, opts);
        OperationLog::Info(LogViewer, "ComputeWorker",
            "RunEnsemble: " + std::to_string(results.size()) + " frames");

    } else {
        // Single conformation: PDB, ProtonatedPdb, Orca
        auto& conf = protein.Conformation();

        if (spec.mode == JobMode::Orca && !spec.orca_files.nmr_out_path.empty())
            opts.orca_nmr_path = spec.orca_files.nmr_out_path;

        auto runResult = OperationRunner::Run(conf, opts);
        OperationLog::Info(LogViewer, "ComputeWorker",
            std::string("Run ok=") + (runResult.Ok() ? "true" : "false") +
            " attached=" + std::to_string(runResult.attached.size()));
        if (!runResult.Ok()) {
            emit progress(15, 100, QString("Run error: %1")
                .arg(QString::fromStdString(runResult.error)));
            // Continue — partial results are still useful
        }
    }

    // Write features if output dir specified
    if (!spec.output_dir.empty()) {
        fs::create_directories(spec.output_dir);
        if (spec.mode == JobMode::Fleet) {
            for (size_t i = 0; i < protein.ConformationCount(); ++i) {
                auto& conf = protein.ConformationAt(i);
                std::string frame_dir = spec.output_dir + "/frame_" + std::to_string(i + 1);
                fs::create_directories(frame_dir);
                ConformationResult::WriteAllFeatures(conf, frame_dir);
            }
        } else {
            auto& conf = protein.Conformation();
            ConformationResult::WriteAllFeatures(conf, spec.output_dir);
        }
        OperationLog::Info(LogViewer, "ComputeWorker",
            "features written to " + spec.output_dir);
    }

    if (cancelled_.load()) { emit finished(result); return; }

    // ================================================================
    // Phase 3: T0 field grid computation (for isosurfaces around rings)
    // Use first conformation for viewer grids.
    // ================================================================
    auto& conf = protein.Conformation();

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
    // Phase 4: Butterfly field computation (always for viewer)
    // ================================================================
    if (protein.RingCount() > 0 && conf.HasResult<BiotSavartResult>()) {
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

// ComputeWorker: dispatches on JobSpec.mode to the correct builder,
// then runs OperationRunner and computes viewer grids.
//
// Data flow:
//   JobSpec → builder → OperationRunner::Run → protein is fully const → emit
//
// The viewer reads the library objects directly. No copy into flat structs.

#include "AIMNet2Result.h"   // must precede any GROMACS headers (DIM macro conflict)

#include "ComputeWorker.h"

#include "Atom.h"
#include "BiotSavartResult.h"
#include "Bond.h"
#include "BuildResult.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "ConformationResult.h"
#include "MolecularGraphResult.h"
#include "OperationLog.h"
#include "OperationRunner.h"
#include "OrcaRunLoader.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "Ring.h"
#include "Session.h"
#include "SpatialIndexResult.h"

#include <QElapsedTimer>
#include <cmath>
#include <filesystem>
#include <stdexcept>  // TrajectoryH5 ctor throws on structural failure

namespace fs = std::filesystem;
using namespace nmr;

ComputeWorker::ComputeWorker(nmr::Session& session, QObject* parent) : QObject(parent), session_(session) {}

void ComputeWorker::cancel() {
    cancelled_.store(true, std::memory_order_relaxed);
}

void ComputeWorker::computeAll(const nmr::JobSpec& spec) {
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

    case JobMode::Trajectory:
        // Trajectory mode is CLI-only (two-pass batch over full-system XTC).
        // Use nmr_extract --trajectory. Not loadable interactively.
        emit finished(result);
        return;

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

    // Session owns the AIMNet2 model, TripeptideDftTable (libpq), and
    // LarsenHBondGrid for the entire process. Pre-loaded in
    // main_viewer.cpp before MainWindow construction. Each per-run
    // RunOptions just borrows the pointers. Null tripeptide/larsen
    // pointers (DSN unconfigured, grid dir empty) silently skip those
    // calculator families — same contract nmr_extract uses.
    RunOptions opts;
    if (buildResult.charges) {
        opts.charge_source = buildResult.charges.get();
        opts.net_charge = buildResult.net_charge;
    }
    opts.skip_mopac    = true;                // MOPAC is batch/calibration-only
    opts.skip_coulomb  = true;                // APBS preferred (6x faster, better physics)
    opts.skip_apbs     = spec.skip_apbs;      // APBS runs unless --no-apbs
    opts.aimnet2_model = session_.Aimnet2Model();
    opts.tripeptide_dft_table = session_.TripeptideDftTablePtr();
    opts.larsen_hbond_grid = session_.LarsenHBondGridPtr();

    if (spec.mode == JobMode::Mutant) {
        // Mutant: run both WT and ALA, compute delta
        auto& wt_conf = protein.Conformation();
        auto& ala_conf = alaBuild.protein->Conformation();

        RunOptions wt_opts = opts;
        if (!spec.wt_files.nmr_out_path.empty()) {
            wt_opts.orca_nmr_path = spec.wt_files.nmr_out_path;
        }

        RunOptions ala_opts;
        ala_opts.charge_source = alaBuild.charges.get();
        ala_opts.net_charge = alaBuild.net_charge;
        ala_opts.skip_mopac    = true;
        ala_opts.skip_coulomb  = true;
        ala_opts.skip_apbs     = spec.skip_apbs;
        ala_opts.aimnet2_model = session_.Aimnet2Model();
        ala_opts.tripeptide_dft_table = session_.TripeptideDftTablePtr();
        ala_opts.larsen_hbond_grid = session_.LarsenHBondGridPtr();
        if (!spec.ala_files.nmr_out_path.empty()) {
            ala_opts.orca_nmr_path = spec.ala_files.nmr_out_path;
        }

        auto runResult = OperationRunner::RunMutantComparison(
            wt_conf, wt_opts, ala_conf, ala_opts);
        OperationLog::Info(LogViewer, "ComputeWorker",
            std::string("RunMutantComparison ok=") + (runResult.Ok() ? "true" : "false") +
            " attached=" + std::to_string(runResult.attached.size()));

        // Attach MolecularGraphResult on both sides so the Inspector
        // surfaces graph topology for the displayed (WT) conformation
        // and the ALA companion stays symmetric. Same rationale as the
        // single-conformation branch below — OperationRunner does not
        // attach it; the viewer does.
        auto attachGraph = [](ProteinConformation& c) {
            if (c.HasResult<SpatialIndexResult>() && !c.HasResult<MolecularGraphResult>()) {
                auto g = MolecularGraphResult::Compute(c);
                if (g) {
                    c.AttachResult(std::move(g));
                }
            }
        };
        attachGraph(wt_conf);
        attachGraph(ala_conf);

    } else {
        // Single conformation: PDB, ProtonatedPdb, Orca
        auto& conf = protein.Conformation();

        if (spec.mode == JobMode::Orca && !spec.orca_files.nmr_out_path.empty()) {
            opts.orca_nmr_path = spec.orca_files.nmr_out_path;
        }

        auto runResult = OperationRunner::Run(conf, opts);
        OperationLog::Info(LogViewer, "ComputeWorker",
            std::string("Run ok=") + (runResult.Ok() ? "true" : "false") +
            " attached=" + std::to_string(runResult.attached.size()));
        if (!runResult.Ok()) {
            emit progress(15, 100, QString("Run error: %1")
                .arg(QString::fromStdString(runResult.error)));
            // Continue — partial results are still useful
        }

        // MolecularGraphResult: BFS through-bond features (graph_dist_*,
        // eneg_sum_*, n_pi_bonds_3, is_conjugated, bfs_decay). Library's
        // OperationRunner::Run does not attach this result; only
        // MutationDeltaResult consumes it. The viewer attaches it here
        // so the per-atom Inspector "Graph topology" section reads real
        // values instead of default sentinels. Dependency is just
        // SpatialIndexResult (already attached by OperationRunner).
        if (conf.HasResult<SpatialIndexResult>()) {
            auto graph = MolecularGraphResult::Compute(conf);
            if (graph) {
                conf.AttachResult(std::move(graph));
                OperationLog::Info(LogViewer, "ComputeWorker", "MolecularGraphResult: attached");
            }
        }
    }

    // Write features if output dir specified
    if (!spec.output_dir.empty()) {
        fs::create_directories(spec.output_dir);
        {
            auto& conf = protein.Conformation();
            ConformationResult::WriteAllFeatures(conf, spec.output_dir);
        }
        OperationLog::Info(LogViewer, "ComputeWorker",
            "features written to " + spec.output_dir);
    }

    if (cancelled_.load()) { emit finished(result); return; }

    // ================================================================
    // Phase 2b: Load companion trajectory H5 (read-only).
    //
    // The viewer never writes H5 or triggers new extractions. If the
    // spec carries a path, the file was proven to exist by
    // ValidateJobSpec. Here we construct TrajectoryH5 (which validates
    // the structural minimum — /atoms, /trajectory/frames, root attrs
    // — and eagerly loads the bounded slices the inspector consumes),
    // verify atom-count + per-index element agreement against the
    // library Protein, then publish the result as
    // shared_ptr<const TrajectoryH5>. Identity mismatch = leave null
    // and log the specific discrepancy; the viewer renders the protein
    // normally, the Time Series tab is empty.
    //
    // HDF5/HighFive is an external library boundary: exceptions are
    // caught here and converted to log lines (same pattern as DSSP,
    // cif++, UDP socket in the library).
    // ================================================================
    if (!spec.analysis_h5_path.empty()) {
        emit progress(40, 100, QStringLiteral("Loading trajectory H5..."));
        OperationLog::Info(LogViewer, "ComputeWorker", "Phase2b: TrajectoryH5 " + spec.analysis_h5_path);
        std::shared_ptr<TrajectoryH5> h5;
        try {
            h5 = std::make_shared<TrajectoryH5>(spec.analysis_h5_path);
        } catch (const std::exception& e) {
            OperationLog::Error("ComputeWorker",
                                std::string("TrajectoryH5 read failed for ") + spec.analysis_h5_path + ": " + e.what());
        }

        if (h5) {
            // Identity check: atom count + per-index element. The writer
            // emits `static_cast<int>(a.element)` (Element enum ordinal)
            // in /atoms/element; see TrajectoryProtein::WriteH5.
            // Compare the same way — Element enum ordinal, NOT atomic
            // number. Mismatch = the H5 describes a different protein
            // and every time-series value would land on the wrong atom.
            const size_t n_lib = protein.AtomCount();
            if (h5->AtomCount() != n_lib) {
                // A common cause when the library is short by ~3 H atoms:
                // `--pdb` runs reduce, which re-evaluates HIS/LYS protonation
                // from H-bond geometry and can drop HIP→HID/HIE H atoms that
                // the H5 (built via FullSystemReader from a TPR with HIP)
                // expects to be present. Use `--protonated-pdb` to trust the
                // PDB's existing protonation and skip reduce.
                std::string hint;
                if (h5->AtomCount() > n_lib && (h5->AtomCount() - n_lib) <= 10) {
                    hint = "  (library is short; if running --pdb on an MD"
                           " snapshot, try --protonated-pdb to skip reduce"
                           " and preserve HIP/LYS protonation)";
                }
                OperationLog::Error("ComputeWorker",
                                    "trajectory H5 atom-count mismatch: h5.n_atoms=" + std::to_string(h5->AtomCount())
                                        + " library.AtomCount=" + std::to_string(n_lib) + " — time series will not be attached"
                                        + hint);
            } else {
                size_t first_bad = SIZE_MAX;
                int bad_h5 = 0;
                int bad_lib = 0;
                for (size_t i = 0; i < n_lib; ++i) {
                    int const lib_e = static_cast<int>(protein.AtomAt(i).element);
                    if (h5->ElementAt(i) != lib_e) {
                        first_bad = i;
                        bad_h5 = h5->ElementAt(i);
                        bad_lib = lib_e;
                        break;
                    }
                }
                if (first_bad != SIZE_MAX) {
                    OperationLog::Error("ComputeWorker",
                                        "trajectory H5 element mismatch at atom " + std::to_string(first_bad)
                                            + ": h5=" + std::to_string(bad_h5) + " library=" + std::to_string(bad_lib)
                                            + " (Element enum ordinals)" + " — time series will not be attached");
                } else {
                    // Identity check passed. Log sparse-set diagnostic
                    // so future readers know which TR groups are
                    // available in this file.
                    std::string groups_text;
                    for (const auto& gn : h5->GroupsPresent()) {
                        if (!groups_text.empty())
                            groups_text += ", ";
                        groups_text += gn;
                    }
                    OperationLog::Info(LogViewer,
                                       "ComputeWorker",
                                       "trajectory H5 identity check ok: n_atoms=" + std::to_string(n_lib)
                                           + " n_frames=" + std::to_string(h5->FrameCount()));
                    OperationLog::Info(LogViewer,
                                       "ComputeWorker",
                                       "  groups present (" + std::to_string(h5->GroupsPresent().size())
                                           + "): " + (groups_text.empty() ? std::string("(none)") : groups_text));

                    AnalysisBinding binding;
                    binding.h5 = h5;
                    binding.libToH5.resize(n_lib);
                    for (size_t i = 0; i < n_lib; ++i) {
                        binding.libToH5[i] = i;   // identity — documented contract
                    }

                    // Collect atom-name deltas. The H5 atom name comes
                    // from the wrapped Protein at extraction time
                    // (TrajectoryProtein::WriteH5 emits a.pdb_atom_name);
                    // mismatches surface convention drift between
                    // extraction-time and viewer-time canonicalisation.
                    // Informational only — element check above already
                    // proved alignment.
                    for (size_t i = 0; i < n_lib; ++i) {
                        const std::string& libName = protein.AtomAt(i).pdb_atom_name;
                        const std::string& h5Name = h5->AtomNameAt(i);
                        if (libName != h5Name) {
                            binding.nameMismatches.push_back({i, libName, h5Name});
                        }
                    }
                    OperationLog::Info(LogViewer,
                                       "ComputeWorker",
                                       "AnalysisBinding: identity libToH5, " + std::to_string(binding.nameMismatches.size())
                                           + " atom-name mismatch(es)");
                    const size_t show = std::min<size_t>(5,
                        binding.nameMismatches.size());
                    for (size_t k = 0; k < show; ++k) {
                        const auto& nm = binding.nameMismatches[k];
                        OperationLog::Info(LogViewer, "ComputeWorker",
                            "  atom " + std::to_string(nm.idx) +
                            ": lib='" + nm.lib + "' h5='" + nm.h5 + "'");
                    }

                    result.analysisBinding = std::move(binding);
                }
            }
        }
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
        int const nRings = static_cast<int>(protein.RingCount());

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

            int const nPoints = G * G * G;
            grid.T0.resize(nPoints, 0.0);
            grid.bsT0.resize(nPoints, 0.0);

            for (int iz = 0; iz < G; ++iz) {
                for (int iy = 0; iy < G; ++iy) {
                    for (int ix = 0; ix < G; ++ix) {
                        if (cancelled_.load()) break;

                        Vec3 const pt(grid.origin[0] + ix * grid.spacing[0],
                                      grid.origin[1] + iy * grid.spacing[1],
                                      grid.origin[2] + iz * grid.spacing[2]);

                        double const dist = (pt - geo.center).norm();
                        if (dist > extent * 1.2 || dist < 0.5) continue;

                        auto stResult = bs.SampleShieldingAt(pt);
                        double const bsT0val = stResult.T0;

                        if (!std::isfinite(bsT0val)) continue;

                        int const idx = ix + iy * G + iz * G * G;
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
            Vec3 const arbitrary = (std::abs(n.x()) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
            Vec3 const u = n.cross(arbitrary).normalized();
            Vec3 const v = n.cross(u);

            bfd.positions.reserve(G * G * G);
            bfd.fields.reserve(G * G * G);

            for (int iz = 0; iz < G; ++iz) {
                for (int iy = 0; iy < G; ++iy) {
                    for (int ix = 0; ix < G; ++ix) {
                        double const fx = -extent + 2.0 * extent * ix / (G - 1);
                        double const fy = -extent + 2.0 * extent * iy / (G - 1);
                        double const fz = -extent + 2.0 * extent * iz / (G - 1);
                        Vec3 const pt = geo.center + fx * u + fy * v + fz * n;
                        Vec3 const B = bs.SampleBFieldAt(pt);
                        bfd.positions.push_back(pt);
                        bfd.fields.push_back(B);
                    }
                }
            }
            result.butterflyFields.push_back(std::move(bfd));
        }
    }

    int const N = static_cast<int>(conf.AtomCount());
    OperationLog::Info(LogViewer, "ComputeWorker",
        "Done: " + std::to_string(N) + " atoms, " +
        std::to_string(protein.RingCount()) + " rings, " +
        std::to_string(result.fieldGrids.size()) + " grids");

    emit progress(100, 100,
        QString("Done: %1 atoms, %2 rings, %3 ms")
            .arg(N).arg(protein.RingCount()).arg(timer.elapsed()));
    emit finished(result);
}

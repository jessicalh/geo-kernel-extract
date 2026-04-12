//
// trajectory_scout: scan a full-system trajectory, accumulate per-atom
// water statistics, write a compact CSV catalog for frame selection.
//
// Reads: full-system .xtc + .tpr (topology) + COLVAR + .edr
// Writes: one CSV with per-atom statistics across all frames.
//         No NPY, no full extraction — just the catalog.
//
// Usage:
//   trajectory_scout --tpr FILE --xtc FILE --colvar FILE --edr FILE --output FILE
//
// Output CSV columns (one row per protein atom):
//   atom_idx, element, resname, resnum, atom_name,
//   water_n_first_{mean,std,min,min_frame,max,max_frame},
//   water_emag_{mean,std,min,min_frame,max,max_frame},
//   sasa_{mean,std,min,min_frame,max,max_frame},
//   n_frames_dry, n_frames_exposed
//

#include "FullSystemReader.h"
#include "SolventEnvironment.h"
#include "OperationLog.h"
#include "xtc_reader.h"
#include "PdbFileReader.h"
#include "RuntimeEnvironment.h"

// For SASA we use SpatialIndex + SasaResult directly
#include "ProteinConformation.h"
#include "Protein.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "SasaResult.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace nmr;

// ── Welford accumulator ──────────────────────────────────────────

struct Welford {
    int    count = 0;
    double mean  = 0.0;
    double M2    = 0.0;
    double min_val = 1e30;
    double max_val = -1e30;
    int    min_frame = -1;
    int    max_frame = -1;

    void Update(double x, int frame) {
        ++count;
        double delta = x - mean;
        mean += delta / count;
        double delta2 = x - mean;
        M2 += delta * delta2;

        if (x < min_val) { min_val = x; min_frame = frame; }
        if (x > max_val) { max_val = x; max_frame = frame; }
    }

    double Std() const {
        return (count > 1) ? std::sqrt(M2 / (count - 1)) : 0.0;
    }
};

// ── Per-atom catalog entry ───────────────────────────────────────

struct AtomCatalog {
    Welford water_n_first;
    Welford water_emag;
    Welford sasa;
    int n_frames_dry = 0;      // frames where water_n_first == 0
};

// ── Water shell + E-field magnitude (no tensor, no EFG) ──────────

static constexpr double FIRST_SHELL = 3.5;    // A
static constexpr double EFIELD_CUTOFF = 15.0;  // A
static constexpr double KE = 14.3996;          // V*A/e

static void ComputeWaterScalars(
    const Vec3& atom_pos,
    const SolventEnvironment& solvent,
    int& n_first,
    double& efield_mag)
{
    double first_sq = FIRST_SHELL * FIRST_SHELL;
    double cutoff_sq = EFIELD_CUTOFF * EFIELD_CUTOFF;
    n_first = 0;
    Vec3 E = Vec3::Zero();

    for (size_t wi = 0; wi < solvent.WaterCount(); ++wi) {
        const auto& water = solvent.waters[wi];
        Vec3 r_O = water.O_pos - atom_pos;
        double d_sq = r_O.squaredNorm();

        if (d_sq < first_sq)
            ++n_first;

        if (d_sq > cutoff_sq) continue;

        // E-field from all 3 charge sites (scalar magnitude only)
        auto add_charge = [&](const Vec3& q_pos, double q) {
            Vec3 r = q_pos - atom_pos;
            double r2 = r.squaredNorm();
            if (r2 < 0.01) return;
            double r3 = r2 * std::sqrt(r2);
            E += (KE * q / r3) * r;
        };

        add_charge(water.O_pos,  water.O_charge);
        add_charge(water.H1_pos, water.H_charge);
        add_charge(water.H2_pos, water.H_charge);
    }

    efield_mag = E.norm();
}

// ── CLI ──────────────────────────────────────────────────────────

static std::string GetArg(int argc, char* argv[], const char* flag) {
    for (int i = 1; i < argc - 1; ++i)
        if (std::string(argv[i]) == flag) return argv[i + 1];
    return "";
}

int main(int argc, char* argv[]) {
    std::string tpr_path    = GetArg(argc, argv, "--tpr");
    std::string xtc_path    = GetArg(argc, argv, "--xtc");
    std::string ref_pdb     = GetArg(argc, argv, "--ref");
    std::string output_path = GetArg(argc, argv, "--output");

    if (tpr_path.empty() || xtc_path.empty() || ref_pdb.empty() || output_path.empty()) {
        fprintf(stderr,
            "Usage: trajectory_scout --tpr FILE --xtc FILE --ref FILE --output FILE\n"
            "\n"
            "  --tpr     Full-system TPR (topology + charges)\n"
            "  --xtc     Full-system XTC (protein + water + ions)\n"
            "  --ref     Reference PDB (protein only, for topology)\n"
            "  --output  Output CSV path\n");
        return 1;
    }

    RuntimeEnvironment::Load();
    OperationLog::SetChannelMask(0xFFFFFFFF);

    // 1. Read topology
    FullSystemReader sys_reader;
    if (!sys_reader.ReadTopology(tpr_path)) {
        fprintf(stderr, "ERROR: %s\n", sys_reader.error().c_str());
        return 1;
    }
    auto& topo = sys_reader.Topology();
    fprintf(stderr, "Topology: %zu protein, %zu water, %zu ions\n",
            topo.protein_count, topo.water_count, topo.ion_count);

    // 2. Build protein from reference PDB (for atom names, residues)
    auto build = BuildFromProtonatedPdb(ref_pdb);
    if (!build.Ok()) {
        fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        return 1;
    }
    const size_t N = build.protein->Conformation().AtomCount();
    fprintf(stderr, "Protein: %zu atoms\n", N);

    // 3. Read all XTC frames
    fprintf(stderr, "Reading %s ...\n", xtc_path.c_str());
    auto frames = read_all_xtc_frames(xtc_path);
    fprintf(stderr, "Frames: %zu\n", frames.size());

    // 4. Allocate per-atom catalogs
    std::vector<AtomCatalog> catalog(N);

    // 5. Stream through frames
    for (size_t fi = 0; fi < frames.size(); ++fi) {
        if (fi % 100 == 0)
            fprintf(stderr, "  frame %zu / %zu\n", fi, frames.size());

        // Split frame into protein + solvent
        std::vector<Vec3> protein_pos;
        SolventEnvironment solvent;
        if (!sys_reader.ExtractFrame(frames[fi].x, protein_pos, solvent)) {
            fprintf(stderr, "WARNING: frame %zu extract failed, skipping\n", fi);
            continue;
        }

        // Create temporary conformation for SASA
        auto conf = std::make_unique<ProteinConformation>(
            build.protein.get(), protein_pos);

        // Compute SASA (needs Geometry + SpatialIndex)
        auto geom = GeometryResult::Compute(*conf);
        if (geom) conf->AttachResult(std::move(geom));
        auto spatial = SpatialIndexResult::Compute(*conf);
        if (spatial) conf->AttachResult(std::move(spatial));
        auto sasa = SasaResult::Compute(*conf);
        if (sasa) conf->AttachResult(std::move(sasa));

        // Per-atom: water scalars + SASA
        for (size_t ai = 0; ai < N; ++ai) {
            int n_first = 0;
            double emag = 0.0;
            ComputeWaterScalars(protein_pos[ai], solvent, n_first, emag);

            catalog[ai].water_n_first.Update(
                static_cast<double>(n_first), static_cast<int>(fi));
            catalog[ai].water_emag.Update(emag, static_cast<int>(fi));

            if (n_first == 0) catalog[ai].n_frames_dry++;

            // SASA from the result
            if (sasa) {
                double atom_sasa = conf->AtomAt(ai).atom_sasa;
                catalog[ai].sasa.Update(atom_sasa, static_cast<int>(fi));
            }
        }

        // Conformation dies here — no disk, no NPY
    }

    // 6. Write CSV catalog
    fprintf(stderr, "Writing %s\n", output_path.c_str());
    std::ofstream out(output_path);
    out << "atom_idx,element,resname,resnum,atom_name,"
        << "water_n_first_mean,water_n_first_std,"
        << "water_n_first_min,water_n_first_min_frame,"
        << "water_n_first_max,water_n_first_max_frame,"
        << "water_emag_mean,water_emag_std,"
        << "water_emag_min,water_emag_min_frame,"
        << "water_emag_max,water_emag_max_frame,"
        << "sasa_mean,sasa_std,"
        << "sasa_min,sasa_min_frame,"
        << "sasa_max,sasa_max_frame,"
        << "n_frames_dry,n_frames_total\n";

    auto& ref_conf = build.protein->Conformation();
    for (size_t ai = 0; ai < N; ++ai) {
        const auto& atom = ref_conf.AtomAt(ai);
        const auto& cat = catalog[ai];

        // Element name
        const char* elem = "?";
        switch (atom.Element()) {
            case Element::H: elem = "H"; break;
            case Element::C: elem = "C"; break;
            case Element::N: elem = "N"; break;
            case Element::O: elem = "O"; break;
            case Element::S: elem = "S"; break;
            default: break;
        }

        out << ai << ","
            << elem << ","
            << atom.ResidueName() << ","
            << atom.ResidueNumber() << ","
            << atom.Name() << ","
            << cat.water_n_first.mean << ","
            << cat.water_n_first.Std() << ","
            << cat.water_n_first.min_val << ","
            << cat.water_n_first.min_frame << ","
            << cat.water_n_first.max_val << ","
            << cat.water_n_first.max_frame << ","
            << cat.water_emag.mean << ","
            << cat.water_emag.Std() << ","
            << cat.water_emag.min_val << ","
            << cat.water_emag.min_frame << ","
            << cat.water_emag.max_val << ","
            << cat.water_emag.max_frame << ","
            << cat.sasa.mean << ","
            << cat.sasa.Std() << ","
            << cat.sasa.min_val << ","
            << cat.sasa.min_frame << ","
            << cat.sasa.max_val << ","
            << cat.sasa.max_frame << ","
            << cat.n_frames_dry << ","
            << cat.water_n_first.count << "\n";
    }

    fprintf(stderr, "Done: %zu atoms x %zu frames -> %s\n",
            N, frames.size(), output_path.c_str());
    return 0;
}

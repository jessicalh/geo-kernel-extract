#pragma once
//
// ComputeWorker: bridge between the NMR field library and the Qt/VTK viewer.
//
// Runs all physics calculations on a background QThread:
//   1. Build protein via the builder matching JobSpec.mode
//   2. Run OperationRunner (DSSP, geometry, all calculators)
//   3. Compute viewer-specific grids (field grids, butterfly fields)
//   4. Signal completion — MainWindow reads directly from library objects
//
// After Run completes, the Protein and its ProteinConformation are fully
// const. shared_ptr<Protein> crosses the thread boundary safely.
// The viewer reads the object model directly. No adapter layer.
//

#include <QObject>
#include <QMetaType>
#include <atomic>
#include <vector>
#include <memory>
#include <string>

#include "JobSpec.h"
#include "Types.h"
#include "analysis_file.h"   // standalone H5 reader (read-only)

namespace nmr { class Protein; }

// ---- Viewer-specific grid computations (NOT library data copies) ----
// These sample the BiotSavartResult at arbitrary 3D points for isosurfaces
// and streamline visualisation. They are viewer operations, not adapters.

// Biot-Savart butterfly: magnetic field on a grid around one aromatic ring
struct ViewerButterflyData {
    nmr::Vec3 ringCenter;
    nmr::Vec3 ringNormal;
    double ringRadius;
    std::string ringType;
    int gridDims[3];
    std::vector<nmr::Vec3> positions;
    std::vector<nmr::Vec3> fields;
};

// Shielding T0 field on a structured 3D grid around one aromatic ring.
struct ViewerFieldGrid {
    double origin[3];
    double spacing[3];
    int dims[3];
    std::vector<double> T0;
    std::vector<double> bsT0;
    nmr::Vec3 ringCenter;
    std::string ringType;
};

// ---- Output: the protein model + viewer-specific grid data ----

// AnalysisBinding: explicit per-atom correspondence between the library
// Protein (built from the ns0 PDB, or whatever the mode loaded) and the
// analysis H5.  For every H5 produced by the frozen extraction pipeline,
// the atom ordering matches BuildFromProtonatedPdb on the matching ns0
// snapshot — both derive from the same GROMACS system.  `libToH5` is
// identity today and is the ONE place that fact is declared.  If a
// future pipeline emits a non-identity ordering, the mapping is learned
// here and every read site (`H5IndexFor`) stays correct without change.
//
// Read-side contract.  The viewer never writes H5 or re-runs extraction.
struct AnalysisBinding {
    std::shared_ptr<const AnalysisFile> h5;   // null if no H5 loaded
    std::vector<size_t>                 libToH5;  // currently identity [0..N-1]

    struct NameMismatch {
        size_t       idx;
        std::string  lib;   // protein.AtomAt(idx).pdb_atom_name (ff14SB)
        std::string  h5;    // h5.atoms.atom_name[idx]          (CHARMM)
    };
    std::vector<NameMismatch> nameMismatches;  // informational; logged at bind

    // libToH5 is populated together with h5 or not at all, so the null
    // check alone expresses validity.  No defensive two-part test.
    bool   Valid()              const { return h5 != nullptr; }
    size_t H5IndexFor(size_t i) const { return libToH5[i]; }
};

struct ComputeResult {
    std::shared_ptr<nmr::Protein> protein;
    std::vector<ViewerButterflyData> butterflyFields;
    std::vector<ViewerFieldGrid> fieldGrids;
    std::string proteinName;
    // Optional read-only companion time-series binding.  Valid() iff the
    // H5 was supplied AND passed the atom-count + per-index element
    // identity checks against the library Protein.  The viewer never
    // writes H5 or triggers a new extraction run.
    AnalysisBinding analysisBinding;
};

Q_DECLARE_METATYPE(nmr::JobSpec)
Q_DECLARE_METATYPE(ComputeResult)

class ComputeWorker : public QObject {
    Q_OBJECT
public:
    explicit ComputeWorker(QObject* parent = nullptr);

public slots:
    void computeAll(nmr::JobSpec spec);
    void cancel();

signals:
    void progress(int current, int total, QString phase);
    void finished(ComputeResult result);

private:
    std::atomic<bool> cancelled_{false};
};

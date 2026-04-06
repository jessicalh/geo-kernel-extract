#pragma once
//
// ComputeWorker: bridge between the NMR field library and the Qt/VTK viewer.
//
// Runs all physics calculations on a background QThread:
//   1. Build protein via BuildFromPdb (protonates, assigns charges)
//   2. Run OperationRunner::Run (DSSP, geometry, all 8 calculators)
//   3. Signal completion — MainWindow reads directly from library objects
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

#include "Types.h"

namespace nmr { class Protein; }

// ---- Input: what to compute ----

struct ViewerLoadRequest {
    std::string pdbPath;

    // Comparison mode: WT vs ALA ORCA outputs
    bool comparisonMode = false;
    std::string wtOrcaOut;
    std::string wtXyz;
    std::string alaOrcaOut;
    std::string alaXyz;

    // Butterfly field computation (expensive — only when requested)
    bool computeButterfly = false;
};

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

struct ComputeResult {
    std::shared_ptr<nmr::Protein> protein;
    std::vector<ViewerButterflyData> butterflyFields;
    std::vector<ViewerFieldGrid> fieldGrids;
    std::string proteinName;
};

Q_DECLARE_METATYPE(ViewerLoadRequest)
Q_DECLARE_METATYPE(ComputeResult)

class ComputeWorker : public QObject {
    Q_OBJECT
public:
    explicit ComputeWorker(QObject* parent = nullptr);

public slots:
    void computeAll(ViewerLoadRequest request);
    void cancel();

signals:
    void progress(int current, int total, QString phase);
    void finished(ComputeResult result);

private:
    std::atomic<bool> cancelled_{false};
};

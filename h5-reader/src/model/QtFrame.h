// QtFrame — a per-frame slab view into the AnalysisFile.
//
// One QtFrame per sampled XTC frame. Light-weight: holds the frame
// index and back-pointers to the conformation and the underlying
// AnalysisFile. All per-atom and per-residue accessors compute
// slab offsets into the AnalysisFile's flat vectors on demand — no
// copying of tensor data per frame.
//
// Lifetime: owned by QtConformation. Never outlives it. The raw
// pointers never dangle because QtConformation holds a shared_ptr
// to the AnalysisFile and QtFrames are destroyed with it.

#pragma once

#include "QtRing.h"   // for RingGeometry struct
#include "Types.h"

#include "analysis_file.h"

#include <cstddef>
#include <vector>

namespace h5reader::model {

class QtConformation;

class QtFrame {
public:
    QtFrame(const QtConformation* conformation,
            const AnalysisFile*   h5,
            size_t                tIndex);

    ~QtFrame() = default;
    QtFrame(const QtFrame&)            = default;
    QtFrame& operator=(const QtFrame&) = default;

    // ----- Frame identity -----
    size_t tIndex()            const { return tIndex_; }
    double timePicoseconds()   const;
    int    xtcFrameIndex()     const;
    const QtConformation* conformation() const { return conformation_; }

    // ----- Per-atom positions -----
    // Reads positions/xyz[t, atomIdx, :] and returns a Vec3 in Angstroms.
    Vec3 position(size_t atomIdx) const;

    // Convenience — atom count shortcut.
    size_t atomCount() const;

    // ----- Per-ring geometry -----
    // Reads ring_geometry/data[t, ringIdx, :] and returns a RingGeometry
    // with center, normal, and radius. Returns zero geometry if the
    // layout is invalid (loader logs this case via ErrorBus).
    RingGeometry ringGeometry(size_t ringIdx) const;

    // Returns the ring's vertex positions at this frame, assembled
    // from the ring's atomIndices and position() slabs. The overlays
    // that sample ring-local geometry (polygons, BS/HM evaluators,
    // B-field streamline grids) all need this; centralising here
    // keeps the three call sites in sync when a future extension
    // (e.g. hydrogen inclusion, alternative vertex ordering) appears.
    std::vector<Vec3> ringVertices(size_t ringIdx) const;

    // ----- Per-residue DSSP -----
    // Reads dssp/ss8[t, residueIdx] and returns a typed DsspCode.
    DsspCode dsspCode(size_t residueIdx) const;

    // ========================================================================
    // Per-atom slab accessors for the atom inspector. Each reads one H5
    // dataset at (tIndex_, atomIdx). Zero / default on out-of-range.
    // ========================================================================

    // Ring-current
    SphericalTensor bsShielding(size_t atomIdx) const;
    SphericalTensor hmShielding(size_t atomIdx) const;
    SphericalTensor rsShielding(size_t atomIdx) const;
    Vec3            totalBField(size_t atomIdx) const;
    int             nRings3A(size_t atomIdx) const;
    int             nRings5A(size_t atomIdx) const;
    int             nRings8A(size_t atomIdx) const;
    double          meanRingDist(size_t atomIdx) const;
    double          nearestRingAtom(size_t atomIdx) const;

    // Bond anisotropy (McConnell)
    SphericalTensor mcShielding(size_t atomIdx) const;
    double          mcCOSum(size_t atomIdx) const;
    double          mcNearestCODist(size_t atomIdx) const;
    Vec3            mcNearestCODir(size_t atomIdx) const;

    // Quadrupole / dispersion
    SphericalTensor pqShielding(size_t atomIdx) const;
    SphericalTensor dispShielding(size_t atomIdx) const;

    // Electrostatics (Coulomb ff14SB, APBS, AIMNet2)
    SphericalTensor coulombShielding(size_t atomIdx) const;
    Vec3            coulombETotal(size_t atomIdx) const;
    double          coulombEMagnitude(size_t atomIdx) const;
    SphericalTensor apbsEfg(size_t atomIdx) const;
    Vec3            apbsEfield(size_t atomIdx) const;
    SphericalTensor aimnet2Shielding(size_t atomIdx) const;

    // H-bond
    SphericalTensor hbondShielding(size_t atomIdx) const;
    double          hbondNearestDist(size_t atomIdx) const;
    Vec3            hbondNearestDir(size_t atomIdx) const;
    int             hbondCount35A(size_t atomIdx) const;
    bool            hbondIsDonor(size_t atomIdx) const;
    bool            hbondIsAcceptor(size_t atomIdx) const;

    // SASA
    double          sasa(size_t atomIdx) const;
    Vec3            sasaNormal(size_t atomIdx) const;

    // Water
    Vec3            waterEfield(size_t atomIdx) const;
    int             waterNFirst(size_t atomIdx) const;
    int             waterNSecond(size_t atomIdx) const;
    double          waterHalfShellAsymmetry(size_t atomIdx) const;
    double          waterDipoleCos(size_t atomIdx) const;

    // Charges
    double          aimnet2Charge(size_t atomIdx) const;
    double          eeqCharge(size_t atomIdx) const;
    double          eeqCoordinationNumber(size_t atomIdx) const;

private:
    const QtConformation* conformation_;
    const AnalysisFile*   h5_;
    size_t                tIndex_;
};

}  // namespace h5reader::model

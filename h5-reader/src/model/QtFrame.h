// QtFrame — per-frame typed view into TrajectoryConformation's per-TR buffers.
//
// Light-weight: holds (conformation*, tIndex_) — copyable, value type.
// All accessors delegate to the right per-TR buffer via
// conformation_->h5(). Absent TRs return default-initialised values
// (Vec3::Zero, SphericalTensor{}, 0.0, etc.) so consumer code that
// asks for `frame.bsShielding(atom)` against an H5 without bs_shielding
// reads as "no contribution" naturally.
//
// API preserved from the pre-2026-05-23 QtFrame so existing overlays /
// inspectors compile without changes; the implementation reads from
// QtTrajectoryH5 typed buffers. Consumers wanting explicit "absent vs
// zero" semantics should consult TrajectoryConformation::h5()->xxx() pointer-
// or-null directly, or use QtFrameAtomView's typed-optional accessors.

#pragma once

#include "QtRing.h"  // RingGeometry, OrthoBasisFromNormal
#include "Types.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace h5reader::model {

class TrajectoryConformation;

class QtFrame {
public:
    QtFrame(const TrajectoryConformation* conformation, std::size_t tIndex);

    ~QtFrame() = default;
    QtFrame(const QtFrame&) = default;
    QtFrame& operator=(const QtFrame&) = default;

    // ----- Frame identity -----
    std::size_t tIndex() const { return tIndex_; }
    double timePicoseconds() const;
    int xtcFrameIndex() const;
    const TrajectoryConformation* conformation() const { return conformation_; }

    // ----- Per-atom positions -----
    Vec3 position(std::size_t atomIdx) const;

    // Convenience
    std::size_t atomCount() const;

    // ----- Per-ring geometry (computed from vertex positions) -----
    RingGeometry ringGeometry(std::size_t ringIdx) const;
    std::vector<Vec3> ringVertices(std::size_t ringIdx) const;

    // ----- Per-residue DSSP -----
    DsspCode dsspCode(std::size_t residueIdx) const;

    // ============================================================
    // Per-atom slab accessors (API preserved from pre-rewrite).
    // Default-initialised values when the source TR is absent.
    // ============================================================

    // Ring-current
    SphericalTensor bsShielding(std::size_t atomIdx) const;
    SphericalTensor hmShielding(std::size_t atomIdx) const;
    SphericalTensor rsShielding(std::size_t atomIdx) const;  // ringchi
    Vec3 totalBField(std::size_t atomIdx) const;
    int nRings3A(std::size_t atomIdx) const;  // derived; ring_neighbourhood
    int nRings5A(std::size_t atomIdx) const;
    int nRings8A(std::size_t atomIdx) const;
    double meanRingDist(std::size_t atomIdx) const;
    double nearestRingAtom(std::size_t atomIdx) const;

    // Bond anisotropy (McConnell)
    SphericalTensor mcShielding(std::size_t atomIdx) const;
    double mcCOSum(std::size_t atomIdx) const;
    double mcNearestCODist(std::size_t atomIdx) const;
    Vec3 mcNearestCODir(std::size_t atomIdx) const;

    // Quadrupole / dispersion
    SphericalTensor pqShielding(std::size_t atomIdx) const;
    SphericalTensor dispShielding(std::size_t atomIdx) const;

    // Electrostatics
    SphericalTensor coulombShielding(std::size_t atomIdx) const;
    Vec3 coulombETotal(std::size_t atomIdx) const;
    double coulombEMagnitude(std::size_t atomIdx) const;
    SphericalTensor apbsEfg(std::size_t atomIdx) const;  // T2-only; T0/T1 zero
    Vec3 apbsEfield(std::size_t atomIdx) const;
    SphericalTensor aimnet2Shielding(std::size_t atomIdx) const;

    // H-bond (kernel form)
    SphericalTensor hbondShielding(std::size_t atomIdx) const;
    double hbondNearestDist(std::size_t atomIdx) const;
    Vec3 hbondNearestDir(std::size_t atomIdx) const;
    int hbondCount35A(std::size_t atomIdx) const;
    bool hbondIsDonor(std::size_t atomIdx) const;
    bool hbondIsAcceptor(std::size_t atomIdx) const;

    // SASA
    double sasa(std::size_t atomIdx) const;
    Vec3 sasaNormal(std::size_t atomIdx) const;

    // Water
    Vec3 waterEfield(std::size_t atomIdx) const;
    int waterNFirst(std::size_t atomIdx) const;
    int waterNSecond(std::size_t atomIdx) const;
    double waterHalfShellAsymmetry(std::size_t atomIdx) const;
    double waterDipoleCos(std::size_t atomIdx) const;

    // Charges
    double aimnet2Charge(std::size_t atomIdx) const;
    double eeqCharge(std::size_t atomIdx) const;
    double eeqCoordinationNumber(std::size_t atomIdx) const;

private:
    const TrajectoryConformation* conformation_ = nullptr;
    std::size_t tIndex_ = 0;
};

}  // namespace h5reader::model

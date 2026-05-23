// QtFrame implementation — accessors delegate to QtTrajectoryH5
// per-TR buffers. Absent TRs return default values (zero / empty)
// per the "absent-is-zero" backward-compat policy.
//
// Ring geometry is computed from the ring's vertex positions per call
// (the new format dropped the per-frame ring_geometry dataset; ring
// vertices come from ring_membership.npy and per-frame positions).

#include "QtFrame.h"

#include "../io/QtTrajectoryH5.h"
#include "QtConformation.h"
#include "QtProtein.h"
#include "QtTopology.h"

#include <Eigen/SVD>
#include <cstdint>

namespace h5reader::model {

namespace {

// Compute ring geometry (center, normal, radius) from vertex positions
// via SVD. Mirrors nmr::Ring::ComputeGeometry's least-squares plane
// fit shape.
RingGeometry FitRingGeometry(const std::vector<Vec3>& verts) {
    RingGeometry g;
    if (verts.empty())
        return g;

    // Center = centroid
    g.center = Vec3::Zero();
    for (const auto& v : verts)
        g.center += v;
    g.center /= static_cast<double>(verts.size());

    // Normal = smallest singular vector of the centered point matrix
    if (verts.size() >= 3) {
        Eigen::Matrix<double, Eigen::Dynamic, 3> M(verts.size(), 3);
        for (std::size_t i = 0; i < verts.size(); ++i)
            M.row(static_cast<Eigen::Index>(i)) = (verts[i] - g.center).transpose();
        Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, 3>> svd(M, Eigen::ComputeFullV);
        g.normal = svd.matrixV().col(2);
        if (g.normal.norm() > 1e-12)
            g.normal.normalize();
    }

    // Radius = mean distance from center
    double r_sum = 0.0;
    for (const auto& v : verts)
        r_sum += (v - g.center).norm();
    g.radius = r_sum / static_cast<double>(verts.size());
    return g;
}

}  // namespace


QtFrame::QtFrame(const QtConformation* conformation, std::size_t tIndex) : conformation_(conformation), tIndex_(tIndex) {}


// ── Frame identity ──────────────────────────────────────────────────

double QtFrame::timePicoseconds() const {
    if (!conformation_ || !conformation_->h5())
        return 0.0;
    const auto& ft = conformation_->h5()->frameTimes();
    return tIndex_ < ft.size() ? ft[tIndex_] : 0.0;
}

int QtFrame::xtcFrameIndex() const {
    if (!conformation_ || !conformation_->h5())
        return -1;
    const auto& fi = conformation_->h5()->frameIndices();
    return tIndex_ < fi.size() ? static_cast<int>(fi[tIndex_]) : -1;
}

std::size_t QtFrame::atomCount() const {
    return (conformation_ && conformation_->h5()) ? conformation_->h5()->atomCount() : 0;
}


// ── Positions ────────────────────────────────────────────────────────

Vec3 QtFrame::position(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return Vec3::Zero();
    const auto* pos = conformation_->h5()->positions();
    if (!pos)
        return Vec3::Zero();
    return pos->at(atomIdx, tIndex_);
}


// ── Ring geometry ────────────────────────────────────────────────────

RingGeometry QtFrame::ringGeometry(std::size_t ringIdx) const {
    if (!conformation_)
        return {};
    const QtProtein* p = conformation_->protein();
    if (!p || ringIdx >= p->ringCount())
        return {};
    return FitRingGeometry(ringVertices(ringIdx));
}

std::vector<Vec3> QtFrame::ringVertices(std::size_t ringIdx) const {
    std::vector<Vec3> verts;
    if (!conformation_)
        return verts;
    const QtProtein* p = conformation_->protein();
    if (!p || ringIdx >= p->ringCount())
        return verts;
    const QtRing& ring = p->ring(ringIdx);
    verts.reserve(ring.atomIndices.size());
    for (int32_t a : ring.atomIndices)
        verts.push_back(position(static_cast<std::size_t>(a)));
    return verts;
}


// ── DSSP ────────────────────────────────────────────────────────────

DsspCode QtFrame::dsspCode(std::size_t residueIdx) const {
    if (!conformation_ || !conformation_->h5())
        return DsspCode::Unknown;
    const auto* ds = conformation_->h5()->dssp8();
    if (!ds || !ds->sourceAttachedAt(tIndex_))
        return DsspCode::Unknown;
    return ds->codeAt(residueIdx, tIndex_);
}


// ── Ring-current accessors ─────────────────────────────────────────

SphericalTensor QtFrame::bsShielding(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return {};
    const auto* ts = conformation_->h5()->bsShielding();
    return ts ? ts->at(atomIdx, tIndex_) : SphericalTensor{};
}
SphericalTensor QtFrame::hmShielding(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return {};
    const auto* ts = conformation_->h5()->hmShielding();
    return ts ? ts->at(atomIdx, tIndex_) : SphericalTensor{};
}
SphericalTensor QtFrame::rsShielding(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return {};
    const auto* ts = conformation_->h5()->ringChiShielding();
    return ts ? ts->at(atomIdx, tIndex_) : SphericalTensor{};
}
Vec3 QtFrame::totalBField(std::size_t /*atomIdx*/) const {
    // Not in new format directly; would need to sum BS B-fields per ring
    // via QtBiotSavartCalc. Returning zero in v1; Session 5 wires the
    // recompute on demand.
    return Vec3::Zero();
}
int QtFrame::nRings3A(std::size_t atomIdx) const {
    // Derive from ring_neighbourhood_trajectory_stats: count ring slots
    // with distance <= 3 Å at this frame.
    if (!conformation_ || !conformation_->h5())
        return 0;
    const auto* rn = conformation_->h5()->ringNeighbourhood();
    if (!rn)
        return 0;
    int n = 0;
    for (std::size_t s = 0; s < rn->n_slots; ++s) {
        if (rn->ringIndexAt(atomIdx, s) < 0)
            continue;
        const auto ch = rn->at(atomIdx, tIndex_, s);
        if (ch[0] > 0 && ch[0] <= 3.0)
            ++n;
    }
    return n;
}
int QtFrame::nRings5A(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return 0;
    const auto* rn = conformation_->h5()->ringNeighbourhood();
    if (!rn)
        return 0;
    int n = 0;
    for (std::size_t s = 0; s < rn->n_slots; ++s) {
        if (rn->ringIndexAt(atomIdx, s) < 0)
            continue;
        const auto ch = rn->at(atomIdx, tIndex_, s);
        if (ch[0] > 0 && ch[0] <= 5.0)
            ++n;
    }
    return n;
}
int QtFrame::nRings8A(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return 0;
    const auto* rn = conformation_->h5()->ringNeighbourhood();
    if (!rn)
        return 0;
    int n = 0;
    for (std::size_t s = 0; s < rn->n_slots; ++s) {
        if (rn->ringIndexAt(atomIdx, s) < 0)
            continue;
        const auto ch = rn->at(atomIdx, tIndex_, s);
        if (ch[0] > 0 && ch[0] <= 8.0)
            ++n;
    }
    return n;
}
double QtFrame::meanRingDist(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return 0.0;
    const auto* rn = conformation_->h5()->ringNeighbourhood();
    if (!rn)
        return 0.0;
    double sum = 0.0;
    int n = 0;
    for (std::size_t s = 0; s < rn->n_slots; ++s) {
        if (rn->ringIndexAt(atomIdx, s) < 0)
            continue;
        sum += rn->at(atomIdx, tIndex_, s)[0];
        ++n;
    }
    return n > 0 ? sum / n : 0.0;
}
double QtFrame::nearestRingAtom(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return 0.0;
    const auto* rn = conformation_->h5()->ringNeighbourhood();
    if (!rn)
        return 0.0;
    double mn = -1.0;
    for (std::size_t s = 0; s < rn->n_slots; ++s) {
        if (rn->ringIndexAt(atomIdx, s) < 0)
            continue;
        const double d = rn->at(atomIdx, tIndex_, s)[0];
        if (d > 0 && (mn < 0 || d < mn))
            mn = d;
    }
    return mn < 0 ? 0.0 : mn;
}


// ── McConnell ─────────────────────────────────────────────────────

SphericalTensor QtFrame::mcShielding(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return {};
    const auto* ts = conformation_->h5()->mcShielding();
    return ts ? ts->at(atomIdx, tIndex_) : SphericalTensor{};
}
double QtFrame::mcCOSum(std::size_t /*atomIdx*/) const {
    return 0.0;
}  // not in new format
double QtFrame::mcNearestCODist(std::size_t /*atomIdx*/) const {
    return 0.0;
}
Vec3 QtFrame::mcNearestCODir(std::size_t /*atomIdx*/) const {
    return Vec3::Zero();
}


// ── Quadrupole / Dispersion ──────────────────────────────────────

SphericalTensor QtFrame::pqShielding(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return {};
    const auto* ts = conformation_->h5()->piQuadShielding();
    return ts ? ts->at(atomIdx, tIndex_) : SphericalTensor{};
}
SphericalTensor QtFrame::dispShielding(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return {};
    const auto* ts = conformation_->h5()->dispShielding();
    return ts ? ts->at(atomIdx, tIndex_) : SphericalTensor{};
}


// ── Electrostatics ───────────────────────────────────────────────

SphericalTensor QtFrame::coulombShielding(std::size_t /*atomIdx*/) const {
    // The current format does not emit a coulomb shielding TR — only
    // the per-frame APBS efield/efg + AIMNet2 charges. Return zero so
    // legacy overlays compile; consumer code uses what's available.
    return {};
}
Vec3 QtFrame::coulombETotal(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return Vec3::Zero();
    const auto* ts = conformation_->h5()->apbsEfield();
    return ts ? ts->at(atomIdx, tIndex_) : Vec3::Zero();
}
double QtFrame::coulombEMagnitude(std::size_t atomIdx) const {
    return coulombETotal(atomIdx).norm();
}
SphericalTensor QtFrame::apbsEfg(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return {};
    const auto* ts = conformation_->h5()->apbsEfg();
    if (!ts)
        return {};
    SphericalTensor st;
    const auto t2 = ts->at(atomIdx, tIndex_);
    for (std::size_t k = 0; k < 5; ++k)
        st.T2[k] = t2[k];
    return st;
}
Vec3 QtFrame::apbsEfield(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return Vec3::Zero();
    const auto* ts = conformation_->h5()->apbsEfield();
    return ts ? ts->at(atomIdx, tIndex_) : Vec3::Zero();
}
SphericalTensor QtFrame::aimnet2Shielding(std::size_t /*atomIdx*/) const {
    return {};  // Not exposed; AIMNet2 in new format = embedding + charge + crg
}


// ── H-bond ───────────────────────────────────────────────────────

SphericalTensor QtFrame::hbondShielding(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return {};
    const auto* ts = conformation_->h5()->hbondShielding();
    return ts ? ts->at(atomIdx, tIndex_) : SphericalTensor{};
}
double QtFrame::hbondNearestDist(std::size_t /*atomIdx*/) const {
    return 0.0;
}
Vec3 QtFrame::hbondNearestDir(std::size_t /*atomIdx*/) const {
    return Vec3::Zero();
}
int QtFrame::hbondCount35A(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return 0;
    const auto* ts = conformation_->h5()->larsenHBondCount();
    return ts ? static_cast<int>(ts->at(atomIdx, tIndex_)) : 0;
}
bool QtFrame::hbondIsDonor(std::size_t /*atomIdx*/) const {
    return false;
}
bool QtFrame::hbondIsAcceptor(std::size_t /*atomIdx*/) const {
    return false;
}


// ── SASA ─────────────────────────────────────────────────────────

double QtFrame::sasa(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return 0.0;
    const auto* ts = conformation_->h5()->sasa();
    return ts ? ts->at(atomIdx, tIndex_) : 0.0;
}
Vec3 QtFrame::sasaNormal(std::size_t /*atomIdx*/) const {
    return Vec3::Zero();
}


// ── Water ────────────────────────────────────────────────────────

Vec3 QtFrame::waterEfield(std::size_t /*atomIdx*/) const {
    return Vec3::Zero();
}
int QtFrame::waterNFirst(std::size_t /*atomIdx*/) const {
    return 0;
}
int QtFrame::waterNSecond(std::size_t /*atomIdx*/) const {
    return 0;
}
double QtFrame::waterHalfShellAsymmetry(std::size_t /*atomIdx*/) const {
    return 0.0;
}
double QtFrame::waterDipoleCos(std::size_t /*atomIdx*/) const {
    return 0.0;
}


// ── Charges ──────────────────────────────────────────────────────

double QtFrame::aimnet2Charge(std::size_t atomIdx) const {
    if (!conformation_ || !conformation_->h5())
        return 0.0;
    const auto* ts = conformation_->h5()->aimnet2Charge();
    return ts ? ts->at(atomIdx, tIndex_) : 0.0;
}
double QtFrame::eeqCharge(std::size_t /*atomIdx*/) const {
    // EEQ is in the eeq_welford TR (no per-frame TS); return Welford mean
    // as a "current value" proxy. Or 0.0 if absent.
    if (!conformation_ || !conformation_->h5())
        return 0.0;
    const auto* w = conformation_->h5()->eeqWelford();
    return w ? 0.0 : 0.0;  // intentional placeholder until Session 5
}
double QtFrame::eeqCoordinationNumber(std::size_t /*atomIdx*/) const {
    return 0.0;
}


}  // namespace h5reader::model

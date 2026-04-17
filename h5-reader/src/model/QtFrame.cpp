#include "QtFrame.h"
#include "QtConformation.h"
#include "QtProtein.h"

#include <cstdint>

namespace h5reader::model {

QtFrame::QtFrame(const QtConformation* conformation,
                 const AnalysisFile*   h5,
                 size_t                tIndex)
    : conformation_(conformation), h5_(h5), tIndex_(tIndex) {}

double QtFrame::timePicoseconds() const {
    if (tIndex_ < h5_->meta.frame_times.size())
        return h5_->meta.frame_times[tIndex_];
    return 0.0;
}

int QtFrame::xtcFrameIndex() const {
    if (tIndex_ < h5_->meta.frame_indices.size())
        return h5_->meta.frame_indices[tIndex_];
    return -1;
}

size_t QtFrame::atomCount() const {
    return h5_->n_atoms;
}

Vec3 QtFrame::position(size_t atomIdx) const {
    // positions/xyz is (T, N, 3) row-major.
    const size_t N = h5_->n_atoms;
    const size_t base = tIndex_ * N * 3 + atomIdx * 3;
    return Vec3(h5_->positions.xyz[base + 0],
                h5_->positions.xyz[base + 1],
                h5_->positions.xyz[base + 2]);
}

RingGeometry QtFrame::ringGeometry(size_t ringIdx) const {
    // ring_geometry/data is (T, n_rings, 7) row-major. Column layout
    // comes from QtConformation's parsed RingGeometryLayout.
    RingGeometry g;
    const auto& layout = conformation_->ringGeometryLayout();
    if (!layout.IsValid()) return g;

    const size_t R = conformation_->ringCount();
    if (ringIdx >= R) return g;
    if (h5_->ring_geometry.data.empty()) return g;

    const size_t base = tIndex_ * R * 7 + ringIdx * 7;
    const auto& d = h5_->ring_geometry.data;
    g.center = Vec3(d[base + layout.centerX],
                    d[base + layout.centerY],
                    d[base + layout.centerZ]);
    g.normal = Vec3(d[base + layout.normalX],
                    d[base + layout.normalY],
                    d[base + layout.normalZ]);
    g.radius = d[base + layout.radius];
    return g;
}

std::vector<Vec3> QtFrame::ringVertices(size_t ringIdx) const {
    std::vector<Vec3> verts;
    if (!conformation_) return verts;
    const auto* p = conformation_->protein();
    if (!p || ringIdx >= p->ringCount()) return verts;
    const auto& ring = p->ring(ringIdx);
    verts.reserve(ring.atomIndices.size());
    for (size_t a : ring.atomIndices) verts.push_back(position(a));
    return verts;
}

// ============================================================================
// Per-atom slab helpers. The H5 stores per-atom tensors/vectors/scalars as
// flat vectors — we index into them with (t, atom, component). One helper
// per data shape (SphericalTensor 9-float, Vec3 3-float, scalar double,
// scalar int16, scalar int8-as-bool).
// ============================================================================

namespace {

SphericalTensor ReadSpherical(const std::vector<double>& v,
                               size_t N, size_t t, size_t a) {
    SphericalTensor st;
    if (v.empty()) return st;
    const size_t base = t * N * 9 + a * 9;
    if (base + 9 > v.size()) return st;
    st.T0 = v[base + 0];
    st.T1[0] = v[base + 1]; st.T1[1] = v[base + 2]; st.T1[2] = v[base + 3];
    st.T2[0] = v[base + 4]; st.T2[1] = v[base + 5]; st.T2[2] = v[base + 6];
    st.T2[3] = v[base + 7]; st.T2[4] = v[base + 8];
    return st;
}

Vec3 ReadVec3(const std::vector<double>& v,
               size_t N, size_t t, size_t a) {
    if (v.empty()) return Vec3::Zero();
    const size_t base = t * N * 3 + a * 3;
    if (base + 3 > v.size()) return Vec3::Zero();
    return Vec3(v[base + 0], v[base + 1], v[base + 2]);
}

double ReadScalarD(const std::vector<double>& v,
                    size_t N, size_t t, size_t a) {
    if (v.empty()) return 0.0;
    const size_t idx = t * N + a;
    return idx < v.size() ? v[idx] : 0.0;
}

int ReadScalarI16(const std::vector<int16_t>& v,
                   size_t N, size_t t, size_t a) {
    if (v.empty()) return 0;
    const size_t idx = t * N + a;
    return idx < v.size() ? static_cast<int>(v[idx]) : 0;
}

bool ReadBoolI8(const std::vector<int8_t>& v,
                 size_t N, size_t t, size_t a) {
    if (v.empty()) return false;
    const size_t idx = t * N + a;
    return idx < v.size() && v[idx] != 0;
}

}  // namespace

// ---- Ring-current ------------------------------------------------------

SphericalTensor QtFrame::bsShielding(size_t a) const {
    return ReadSpherical(h5_->ring_current.bs_shielding, h5_->n_atoms, tIndex_, a);
}
SphericalTensor QtFrame::hmShielding(size_t a) const {
    return ReadSpherical(h5_->ring_current.hm_shielding, h5_->n_atoms, tIndex_, a);
}
SphericalTensor QtFrame::rsShielding(size_t a) const {
    return ReadSpherical(h5_->ring_current.rs_shielding, h5_->n_atoms, tIndex_, a);
}
Vec3 QtFrame::totalBField(size_t a) const {
    return ReadVec3(h5_->ring_current.total_B_field, h5_->n_atoms, tIndex_, a);
}
int QtFrame::nRings3A(size_t a) const {
    return ReadScalarI16(h5_->ring_current.n_rings_3A, h5_->n_atoms, tIndex_, a);
}
int QtFrame::nRings5A(size_t a) const {
    return ReadScalarI16(h5_->ring_current.n_rings_5A, h5_->n_atoms, tIndex_, a);
}
int QtFrame::nRings8A(size_t a) const {
    return ReadScalarI16(h5_->ring_current.n_rings_8A, h5_->n_atoms, tIndex_, a);
}
double QtFrame::meanRingDist(size_t a) const {
    return ReadScalarD(h5_->ring_current.mean_ring_dist, h5_->n_atoms, tIndex_, a);
}
double QtFrame::nearestRingAtom(size_t a) const {
    return ReadScalarD(h5_->ring_current.nearest_ring_atom, h5_->n_atoms, tIndex_, a);
}

// ---- Bond anisotropy (McConnell) ---------------------------------------

SphericalTensor QtFrame::mcShielding(size_t a) const {
    return ReadSpherical(h5_->bond_aniso.mc_shielding, h5_->n_atoms, tIndex_, a);
}
double QtFrame::mcCOSum(size_t a) const {
    return ReadScalarD(h5_->bond_aniso.co_sum, h5_->n_atoms, tIndex_, a);
}
double QtFrame::mcNearestCODist(size_t a) const {
    return ReadScalarD(h5_->bond_aniso.nearest_CO_dist, h5_->n_atoms, tIndex_, a);
}
Vec3 QtFrame::mcNearestCODir(size_t a) const {
    return ReadVec3(h5_->bond_aniso.dir_nearest_CO, h5_->n_atoms, tIndex_, a);
}

// ---- Quadrupole / Dispersion -------------------------------------------

SphericalTensor QtFrame::pqShielding(size_t a) const {
    return ReadSpherical(h5_->quadrupole.pq_shielding, h5_->n_atoms, tIndex_, a);
}
SphericalTensor QtFrame::dispShielding(size_t a) const {
    return ReadSpherical(h5_->dispersion.disp_shielding, h5_->n_atoms, tIndex_, a);
}

// ---- Electrostatics (Coulomb ff14SB, APBS, AIMNet2) -------------------

SphericalTensor QtFrame::coulombShielding(size_t a) const {
    return ReadSpherical(h5_->efg.coulomb_shielding, h5_->n_atoms, tIndex_, a);
}
Vec3 QtFrame::coulombETotal(size_t a) const {
    return ReadVec3(h5_->efg.E_total, h5_->n_atoms, tIndex_, a);
}
double QtFrame::coulombEMagnitude(size_t a) const {
    return ReadScalarD(h5_->efg.E_magnitude, h5_->n_atoms, tIndex_, a);
}
SphericalTensor QtFrame::apbsEfg(size_t a) const {
    return ReadSpherical(h5_->efg.apbs_efg, h5_->n_atoms, tIndex_, a);
}
Vec3 QtFrame::apbsEfield(size_t a) const {
    return ReadVec3(h5_->efg.apbs_efield, h5_->n_atoms, tIndex_, a);
}
SphericalTensor QtFrame::aimnet2Shielding(size_t a) const {
    return ReadSpherical(h5_->efg.aimnet2_shielding, h5_->n_atoms, tIndex_, a);
}

// ---- H-bond -------------------------------------------------------------

SphericalTensor QtFrame::hbondShielding(size_t a) const {
    return ReadSpherical(h5_->hbond.hbond_shielding, h5_->n_atoms, tIndex_, a);
}
double QtFrame::hbondNearestDist(size_t a) const {
    return ReadScalarD(h5_->hbond.nearest_dist, h5_->n_atoms, tIndex_, a);
}
Vec3 QtFrame::hbondNearestDir(size_t a) const {
    return ReadVec3(h5_->hbond.nearest_dir, h5_->n_atoms, tIndex_, a);
}
int QtFrame::hbondCount35A(size_t a) const {
    return ReadScalarI16(h5_->hbond.count_3_5A, h5_->n_atoms, tIndex_, a);
}
bool QtFrame::hbondIsDonor(size_t a) const {
    return ReadBoolI8(h5_->hbond.is_donor, h5_->n_atoms, tIndex_, a);
}
bool QtFrame::hbondIsAcceptor(size_t a) const {
    return ReadBoolI8(h5_->hbond.is_acceptor, h5_->n_atoms, tIndex_, a);
}

// ---- SASA ---------------------------------------------------------------

double QtFrame::sasa(size_t a) const {
    return ReadScalarD(h5_->sasa.sasa, h5_->n_atoms, tIndex_, a);
}
Vec3 QtFrame::sasaNormal(size_t a) const {
    return ReadVec3(h5_->sasa.normal, h5_->n_atoms, tIndex_, a);
}

// ---- Water --------------------------------------------------------------

Vec3 QtFrame::waterEfield(size_t a) const {
    return ReadVec3(h5_->water.efield, h5_->n_atoms, tIndex_, a);
}
int QtFrame::waterNFirst(size_t a) const {
    return ReadScalarI16(h5_->water.n_first, h5_->n_atoms, tIndex_, a);
}
int QtFrame::waterNSecond(size_t a) const {
    return ReadScalarI16(h5_->water.n_second, h5_->n_atoms, tIndex_, a);
}
double QtFrame::waterHalfShellAsymmetry(size_t a) const {
    return ReadScalarD(h5_->water.half_shell_asymmetry, h5_->n_atoms, tIndex_, a);
}
double QtFrame::waterDipoleCos(size_t a) const {
    return ReadScalarD(h5_->water.dipole_cos, h5_->n_atoms, tIndex_, a);
}

// ---- Charges ------------------------------------------------------------

double QtFrame::aimnet2Charge(size_t a) const {
    return ReadScalarD(h5_->charges.aimnet2_charge, h5_->n_atoms, tIndex_, a);
}
double QtFrame::eeqCharge(size_t a) const {
    return ReadScalarD(h5_->charges.eeq_charge, h5_->n_atoms, tIndex_, a);
}
double QtFrame::eeqCoordinationNumber(size_t a) const {
    return ReadScalarD(h5_->charges.eeq_cn, h5_->n_atoms, tIndex_, a);
}

DsspCode QtFrame::dsspCode(size_t residueIdx) const {
    // dssp/ss8 is (T, R) row-major int8.
    if (h5_->dssp.ss8.empty()) return DsspCode::Unknown;
    const size_t R = h5_->n_residues;
    if (residueIdx >= R) return DsspCode::Unknown;
    const int8_t v = h5_->dssp.ss8[tIndex_ * R + residueIdx];
    if (v < 0 || v > static_cast<int8_t>(DsspCode::Unknown))
        return DsspCode::Unknown;
    return static_cast<DsspCode>(v);
}

}  // namespace h5reader::model

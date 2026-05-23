#include "TrajectoryH5.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// HighFive boundary lives here, not in the header. Eager-load means
// the HighFive::File handle's lifetime ends in the constructor.
#include "highfive/H5Attribute.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5Group.hpp"

// OperationLog: file-scoped warnings when a TR group is present but
// has a shape that disagrees with the file's declared n_atoms.
// Group-absent is silent (sparse-set normal); group-malformed is loud.
#include "OperationLog.h"

namespace {

// Try to read a 1D dataset into out, return false on absence rather
// than throwing. (HighFive's getDataSet throws when missing; this
// wrapper turns that into a tolerant sparse-set probe.)
template <typename T>
bool TryReadDataset(const HighFive::Group& grp, const std::string& name, std::vector<T>& out) {
    if (!grp.exist(name))
        return false;
    grp.getDataSet(name).read(out);
    return true;
}

// Same shape for attribute reads.
template <typename T>
bool TryReadAttribute(const HighFive::Group& grp, const std::string& name, T& out) {
    if (!grp.hasAttribute(name))
        return false;
    grp.getAttribute(name).read(out);
    return true;
}

template <typename T>
bool TryReadFileAttribute(const HighFive::File& file, const std::string& name, T& out) {
    if (!file.hasAttribute(name))
        return false;
    file.getAttribute(name).read(out);
    return true;
}

// Loud-on-shape-mismatch: when a TR group is present in the file but
// its dataset shapes disagree with the declared n_atoms (or the
// expected per-atom-rank), the file is structurally broken — we leave
// the reader silently absent, but log a Warn so the failure is
// surfaced via UDP. Group-absent stays silent (sparse-set is normal);
// group-malformed is loud.
void WarnShapeMismatch(const char* group_path, const std::string& detail) {
    nmr::OperationLog::Warn("TrajectoryH5", std::string(group_path) + " present but malformed: " + detail);
}

// Tensor-source Welford reader (BS / HM / McConnell). The schema for
// these is identical: t0_{mean,std} + t2mag_{mean,std} as N-shaped
// float64 datasets. The producer ALSO emits t2magnitude_{mean,std}
// (new canonical name) — accept either, prefer t2mag for compat with
// May-8 fixtures.
void ReadShieldingWelford(const HighFive::File& file,
                          const char* group_path,
                          std::size_t n_atoms,
                          std::optional<std::vector<TrajectoryH5::ShieldingWelfordRow>>& out) {
    if (!file.exist(group_path))
        return;
    auto grp = file.getGroup(group_path);
    std::vector<double> t0_mean;
    std::vector<double> t0_std;
    std::vector<double> t2_mean;
    std::vector<double> t2_std;
    const bool has = TryReadDataset(grp, "t0_mean", t0_mean) && TryReadDataset(grp, "t0_std", t0_std)
                     && (TryReadDataset(grp, "t2mag_mean", t2_mean) || TryReadDataset(grp, "t2magnitude_mean", t2_mean))
                     && (TryReadDataset(grp, "t2mag_std", t2_std) || TryReadDataset(grp, "t2magnitude_std", t2_std));
    if (!has) {
        WarnShapeMismatch(group_path, "missing one of t0_mean/t0_std/t2mag_{mean,std}");
        return;
    }
    if (t0_mean.size() != n_atoms || t0_std.size() != n_atoms || t2_mean.size() != n_atoms || t2_std.size() != n_atoms) {
        WarnShapeMismatch(group_path, "dataset length disagrees with n_atoms=" + std::to_string(n_atoms));
        return;
    }
    std::vector<TrajectoryH5::ShieldingWelfordRow> rows(n_atoms);
    for (std::size_t i = 0; i < n_atoms; ++i) {
        rows[i].t0.mean = t0_mean[i];
        rows[i].t0.std = t0_std[i];
        rows[i].t2magnitude.mean = t2_mean[i];
        rows[i].t2magnitude.std = t2_std[i];
    }
    out = std::move(rows);
}

// Shielding time-series frame-0 reader. Reads xyz[:, 0, :] from a
// [N, T, 9] dataset; returns T0 (idx 0) and |T2| (rms of idx 4..8).
void ReadShieldingFrame0(const HighFive::File& file,
                         const char* group_path,
                         std::size_t n_atoms,
                         std::optional<std::vector<TrajectoryH5::ShieldingFrame0Row>>& out) {
    if (!file.exist(group_path))
        return;
    auto grp = file.getGroup(group_path);
    if (!grp.exist("xyz")) {
        WarnShapeMismatch(group_path, "missing 'xyz' dataset");
        return;
    }
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getDimensions();
    if (dims.size() != 3 || dims[0] != n_atoms || dims[2] != 9) {
        WarnShapeMismatch(group_path, "xyz shape != [n_atoms=" + std::to_string(n_atoms) + ", T, 9]");
        return;
    }
    std::vector<double> slab(n_atoms * 9);
    ds.select({0, 0, 0}, {n_atoms, 1, 9}).read(slab.data());
    std::vector<TrajectoryH5::ShieldingFrame0Row> rows(n_atoms);
    for (std::size_t i = 0; i < n_atoms; ++i) {
        const std::size_t base = i * 9;
        rows[i].T0 = slab[base + 0];
        double t2_sq = 0.0;
        for (std::size_t k = 4; k < 9; ++k) {
            const double v = slab[base + k];
            t2_sq += v * v;
        }
        rows[i].T2_magnitude = std::sqrt(t2_sq);
    }
    out = std::move(rows);
}

// Scalar [N, T] time-series frame-0 reader (SASA, AIMNet2 charge).
void ReadScalarFrame0(const HighFive::File& file,
                      const char* group_path,
                      const char* dataset_name,
                      std::size_t n_atoms,
                      std::optional<std::vector<double>>& out) {
    if (!file.exist(group_path))
        return;
    auto grp = file.getGroup(group_path);
    if (!grp.exist(dataset_name)) {
        WarnShapeMismatch(group_path, std::string("missing '") + dataset_name + "' dataset");
        return;
    }
    auto ds = grp.getDataSet(dataset_name);
    const auto dims = ds.getDimensions();
    if (dims.size() != 2 || dims[0] != n_atoms) {
        WarnShapeMismatch(group_path, std::string(dataset_name) + " shape != [n_atoms=" + std::to_string(n_atoms) + ", T]");
        return;
    }
    std::vector<double> slab(n_atoms);
    ds.select({0, 0}, {n_atoms, 1}).read(slab.data());
    out = std::move(slab);
}

// Vec3 [N, T, 3] time-series frame-0 reader (APBS efield).
void ReadVec3Frame0(const HighFive::File& file,
                    const char* group_path,
                    const char* dataset_name,
                    std::size_t n_atoms,
                    std::optional<std::vector<TrajectoryH5::Vec3F0>>& out) {
    if (!file.exist(group_path))
        return;
    auto grp = file.getGroup(group_path);
    if (!grp.exist(dataset_name)) {
        WarnShapeMismatch(group_path, std::string("missing '") + dataset_name + "' dataset");
        return;
    }
    auto ds = grp.getDataSet(dataset_name);
    const auto dims = ds.getDimensions();
    if (dims.size() != 3 || dims[0] != n_atoms || dims[2] != 3) {
        WarnShapeMismatch(group_path, std::string(dataset_name) + " shape != [n_atoms=" + std::to_string(n_atoms) + ", T, 3]");
        return;
    }
    std::vector<double> slab(n_atoms * 3);
    ds.select({0, 0, 0}, {n_atoms, 1, 3}).read(slab.data());
    std::vector<TrajectoryH5::Vec3F0> rows(n_atoms);
    for (std::size_t i = 0; i < n_atoms; ++i) {
        rows[i].x = slab[i * 3 + 0];
        rows[i].y = slab[i * 3 + 1];
        rows[i].z = slab[i * 3 + 2];
    }
    out = std::move(rows);
}

}  // namespace


TrajectoryH5::TrajectoryH5(const std::string& path) {
    using namespace HighFive;
    File const file(path, File::ReadOnly);

    // ── Root attrs (n_atoms required; protein_id best-effort) ────────
    if (!TryReadFileAttribute(file, "n_atoms", n_atoms_)) {
        throw std::runtime_error("TrajectoryH5: root attribute 'n_atoms' missing");
    }
    TryReadFileAttribute(file, "protein_id", protein_id_);

    // ── /atoms identity (required) ──────────────────────────────────
    if (!file.exist("/atoms")) {
        throw std::runtime_error("TrajectoryH5: '/atoms' group missing");
    }
    {
        auto atoms = file.getGroup("/atoms");
        atoms.getDataSet("element").read(atom_element_);
        atoms.getDataSet("residue_index").read(residue_index_);
        atoms.getDataSet("pdb_atom_name").read(atom_name_);
        if (atom_element_.size() != n_atoms_ || residue_index_.size() != n_atoms_ || atom_name_.size() != n_atoms_) {
            throw std::runtime_error("TrajectoryH5: /atoms dataset sizes inconsistent with n_atoms");
        }
    }

    // ── /trajectory/frames (required) ───────────────────────────────
    if (!file.exist("/trajectory/frames")) {
        throw std::runtime_error("TrajectoryH5: '/trajectory/frames' missing");
    }
    {
        auto frames = file.getGroup("/trajectory/frames");
        frames.getDataSet("time_ps").read(frame_times_);
        n_frames_ = frame_times_.size();
    }

    // ── /trajectory attrs (best-effort) ─────────────────────────────
    if (file.exist("/trajectory")) {
        auto traj = file.getGroup("/trajectory");
        TryReadAttribute(traj, "xtc_path", xtc_path_);
        TryReadAttribute(traj, "tpr_path", tpr_path_);
        TryReadAttribute(traj, "edr_path", edr_path_);
        // listObjectNames includes everything directly under /trajectory,
        // so the inspector would see /trajectory/frames (frame metadata)
        // and /trajectory/selections (run-scope event bag) alongside the
        // per-TR groups. Filter both: groups_present_ is for "which TR
        // groups produced data on this atom axis" and the metadata
        // siblings would just clutter the listing.
        const std::vector<std::string> kMetadataGroups = {"frames", "selections"};
        for (const auto& g : traj.listObjectNames()) {
            const bool is_metadata = std::find(kMetadataGroups.begin(), kMetadataGroups.end(), g) != kMetadataGroups.end();
            if (!is_metadata) {
                groups_present_.push_back(g);
            }
        }
    }

    // ── Shielding-source Welford rollups (uniform schema) ──────────
    ReadShieldingWelford(file, "/trajectory/bs_welford", n_atoms_, bs_welford_);
    ReadShieldingWelford(file, "/trajectory/hm_welford", n_atoms_, hm_welford_);
    ReadShieldingWelford(file, "/trajectory/mc_welford", n_atoms_, mc_welford_);

    // ── Scalar-source Welford rollups (per-physics channel name) ───
    // Producer emits both canonical (`<channel>_mean`, `<channel>_std`)
    // and legacy (`mean`, `std`) datasets; accept either to stay
    // compatible with older fixtures.
    if (file.exist("/trajectory/sasa_welford")) {
        auto grp = file.getGroup("/trajectory/sasa_welford");
        std::vector<double> m;
        std::vector<double> s;
        if ((TryReadDataset(grp, "sasa_mean", m) || TryReadDataset(grp, "mean", m))
            && (TryReadDataset(grp, "sasa_std", s) || TryReadDataset(grp, "std", s)) && m.size() == n_atoms_
            && s.size() == n_atoms_) {
            std::vector<SasaWelfordRow> rows(n_atoms_);
            for (std::size_t i = 0; i < n_atoms_; ++i) {
                rows[i].sasa.mean = m[i];
                rows[i].sasa.std = s[i];
            }
            sasa_welford_ = std::move(rows);
        }
    }
    if (file.exist("/trajectory/eeq_welford")) {
        auto grp = file.getGroup("/trajectory/eeq_welford");
        std::vector<double> m;
        std::vector<double> s;
        if ((TryReadDataset(grp, "charge_mean", m) || TryReadDataset(grp, "mean", m))
            && (TryReadDataset(grp, "charge_std", s) || TryReadDataset(grp, "std", s)) && m.size() == n_atoms_
            && s.size() == n_atoms_) {
            std::vector<EeqWelfordRow> rows(n_atoms_);
            for (std::size_t i = 0; i < n_atoms_; ++i) {
                rows[i].charge.mean = m[i];
                rows[i].charge.std = s[i];
            }
            eeq_welford_ = std::move(rows);
        }
    }
    if (file.exist("/trajectory/hbond_count_welford")) {
        auto grp = file.getGroup("/trajectory/hbond_count_welford");
        std::vector<double> m;
        std::vector<double> s;
        if ((TryReadDataset(grp, "count_mean", m) || TryReadDataset(grp, "mean", m))
            && (TryReadDataset(grp, "count_std", s) || TryReadDataset(grp, "std", s)) && m.size() == n_atoms_
            && s.size() == n_atoms_) {
            std::vector<HBondCountWelfordRow> rows(n_atoms_);
            for (std::size_t i = 0; i < n_atoms_; ++i) {
                rows[i].count.mean = m[i];
                rows[i].count.std = s[i];
            }
            hbond_count_welford_ = std::move(rows);
        }
    }

    // ── Tensor shielding time-series frame-0 slabs ─────────────────
    // 7 calculators share xyz[N, T, 9] schema. Irrep layout per the
    // file's `irrep_layout` attr: idx 0 = T0, idx 4..8 = T2 m-basis.
    // T1 deferred (Cartesian/m basis mismatch with runtime
    // SphericalTensor — see header doc).
    ReadShieldingFrame0(file, "/trajectory/bs_shielding_time_series", n_atoms_, bs_shielding_f0_);
    ReadShieldingFrame0(file, "/trajectory/hm_shielding_time_series", n_atoms_, hm_shielding_f0_);
    ReadShieldingFrame0(file, "/trajectory/mc_shielding_time_series", n_atoms_, mc_shielding_f0_);
    ReadShieldingFrame0(file, "/trajectory/piquad_shielding_time_series", n_atoms_, piquad_shielding_f0_);
    ReadShieldingFrame0(file, "/trajectory/ringchi_shielding_time_series", n_atoms_, ringchi_shielding_f0_);
    ReadShieldingFrame0(file, "/trajectory/disp_shielding_time_series", n_atoms_, disp_shielding_f0_);
    ReadShieldingFrame0(file, "/trajectory/hbond_shielding_time_series", n_atoms_, hbond_shielding_f0_);

    // ── Scalar + Vec3 time-series frame-0 slabs ────────────────────
    ReadScalarFrame0(file, "/trajectory/sasa_time_series", "sasa", n_atoms_, sasa_f0_);
    ReadScalarFrame0(file, "/trajectory/aimnet2_charge_time_series", "charge", n_atoms_, aimnet2_charge_f0_);
    ReadVec3Frame0(file, "/trajectory/apbs_efield_time_series", "xyz", n_atoms_, apbs_efield_f0_);
}


// All per-atom accessors gate on (group-present AND atom_idx in range).
// The range check is defensive: today libToH5 is identity and the inner
// vectors are sized n_atoms_, so atom_idx >= n_atoms_ can only happen
// if a future caller wires non-identity libToH5 incorrectly. Returning
// nullopt on OOB is graceful — silent UB would be a debugging
// nightmare.
#define TRAJH5_RETURN_IF_ABSENT_OR_OOB(buffer)                                                                                 \
    do {                                                                                                                       \
        if (!(buffer) || atom_idx >= n_atoms_)                                                                                 \
            return std::nullopt;                                                                                               \
    } while (0)

std::optional<TrajectoryH5::ShieldingWelfordRow> TrajectoryH5::BsWelford(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(bs_welford_);
    return (*bs_welford_)[atom_idx];
}

std::optional<TrajectoryH5::ShieldingWelfordRow> TrajectoryH5::HmWelford(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(hm_welford_);
    return (*hm_welford_)[atom_idx];
}

std::optional<TrajectoryH5::ShieldingWelfordRow> TrajectoryH5::McWelford(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(mc_welford_);
    return (*mc_welford_)[atom_idx];
}

std::optional<TrajectoryH5::SasaWelfordRow> TrajectoryH5::SasaWelford(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(sasa_welford_);
    return (*sasa_welford_)[atom_idx];
}

std::optional<TrajectoryH5::EeqWelfordRow> TrajectoryH5::EeqWelford(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(eeq_welford_);
    return (*eeq_welford_)[atom_idx];
}

std::optional<TrajectoryH5::HBondCountWelfordRow> TrajectoryH5::HBondCountWelford(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(hbond_count_welford_);
    return (*hbond_count_welford_)[atom_idx];
}

std::optional<TrajectoryH5::ShieldingFrame0Row> TrajectoryH5::BsShieldingFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(bs_shielding_f0_);
    return (*bs_shielding_f0_)[atom_idx];
}

std::optional<TrajectoryH5::ShieldingFrame0Row> TrajectoryH5::HmShieldingFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(hm_shielding_f0_);
    return (*hm_shielding_f0_)[atom_idx];
}

std::optional<TrajectoryH5::ShieldingFrame0Row> TrajectoryH5::McShieldingFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(mc_shielding_f0_);
    return (*mc_shielding_f0_)[atom_idx];
}

std::optional<TrajectoryH5::ShieldingFrame0Row> TrajectoryH5::PiQuadShieldingFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(piquad_shielding_f0_);
    return (*piquad_shielding_f0_)[atom_idx];
}

std::optional<TrajectoryH5::ShieldingFrame0Row> TrajectoryH5::RingChiShieldingFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(ringchi_shielding_f0_);
    return (*ringchi_shielding_f0_)[atom_idx];
}

std::optional<TrajectoryH5::ShieldingFrame0Row> TrajectoryH5::DispShieldingFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(disp_shielding_f0_);
    return (*disp_shielding_f0_)[atom_idx];
}

std::optional<TrajectoryH5::ShieldingFrame0Row> TrajectoryH5::HBondShieldingFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(hbond_shielding_f0_);
    return (*hbond_shielding_f0_)[atom_idx];
}

std::optional<double> TrajectoryH5::SasaFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(sasa_f0_);
    return (*sasa_f0_)[atom_idx];
}

std::optional<double> TrajectoryH5::Aimnet2ChargeFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(aimnet2_charge_f0_);
    return (*aimnet2_charge_f0_)[atom_idx];
}

std::optional<TrajectoryH5::Vec3F0> TrajectoryH5::ApbsEfieldFrame0(std::size_t atom_idx) const {
    TRAJH5_RETURN_IF_ABSENT_OR_OOB(apbs_efield_f0_);
    return (*apbs_efield_f0_)[atom_idx];
}

#undef TRAJH5_RETURN_IF_ABSENT_OR_OOB

bool TrajectoryH5::HasGroup(const std::string& name) const {
    return std::find(groups_present_.begin(), groups_present_.end(), name) != groups_present_.end();
}

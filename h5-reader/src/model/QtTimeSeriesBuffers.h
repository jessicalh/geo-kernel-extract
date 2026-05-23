// QtTimeSeriesBuffers.h — atom-axis per-frame typed buffer family.
//
// Each struct represents the eager-loaded contents of one
// /trajectory/<name>/ subgroup in trajectory.h5 for an atom-axis TR
// (shape (N, T, K) for K-component per-atom-per-frame data). Loaded
// once at QtTrajectoryH5 construction; held for the session.
//
// Storage strategy: flat std::vector<double> indexed (atom, frame,
// component) row-major. Per-atom slabs are contiguous in memory, which
// matches how the writer emitted them — see Agent 3's fixture deep-
// dive: shapes are (N, T, K) with N leading.
//
// Per-frame metadata (frame_indices, frame_times) duplicated from
// /trajectory/frames; carried per-TR because some TRs are
// conditional-attach and may align to a subset. source_attached_per_frame
// is the (T,) uint8 mask per the "absent, not faked" discipline —
// empty vector means always-attached.
//
// Display-only string attrs (irrep_layout, units, parity,
// normalization, result_name) are kept on the buffer for inspector /
// glossary display. NEVER dispatched on by any rendering or
// calculator code.

#pragma once

#include "Types.h"

#include <QString>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace h5reader::model {

// ──────────────────────────────────────────────────────────────────
// Common per-TR per-frame metadata (shared shape across all TS types).
// ──────────────────────────────────────────────────────────────────

struct QtTimeSeriesFrameMeta {
    std::vector<uint64_t> frame_indices;   // (T,)
    std::vector<double> frame_times;       // (T,) ps
    std::vector<uint8_t> source_attached;  // (T,) — empty == always-attached

    std::size_t frameCount() const { return frame_times.size(); }

    bool sourceAttachedAt(std::size_t t) const {
        return source_attached.empty() || (t < source_attached.size() && source_attached[t] != 0);
    }
};


// ──────────────────────────────────────────────────────────────────
// QtShieldingTimeSeries — (N, T, 9) tensor TS (T0+T1+T2 spherical).
// Used by ~13 shielding TRs (bs/hm/mc/piquad/ringchi/disp/hbond +
// mopac_coulomb/mopac_mc/mopac_vs_ff14sb + tripeptide_bb/neighbor +
// 4× larsen_hbond + water_field).
// ──────────────────────────────────────────────────────────────────

struct QtShieldingTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<double> xyz;  // (N*T*9,) row-major
    QtTimeSeriesFrameMeta meta;

    // Attrs — display only
    QString irrep_layout;  // "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,..."
    QString units;
    QString parity;         // e.g. "0e+1o+2e"
    QString normalization;  // "isometric_real_sph"
    QString result_name;

    SphericalTensor at(std::size_t atomIdx, std::size_t t) const {
        SphericalTensor st;
        if (atomIdx >= n_atoms || t >= n_frames)
            return st;
        const std::size_t base = (atomIdx * n_frames + t) * 9;
        st.T0 = xyz[base + 0];
        st.T1[0] = xyz[base + 1];
        st.T1[1] = xyz[base + 2];
        st.T1[2] = xyz[base + 3];
        st.T2[0] = xyz[base + 4];
        st.T2[1] = xyz[base + 5];
        st.T2[2] = xyz[base + 6];
        st.T2[3] = xyz[base + 7];
        st.T2[4] = xyz[base + 8];
        return st;
    }

    bool sourceAttachedAt(std::size_t t) const { return meta.sourceAttachedAt(t); }
};


// ──────────────────────────────────────────────────────────────────
// QtScalarTimeSeries — (N, T) per-atom scalar TS.
// Used by sasa, aimnet2_charge, j_coupling, larsen_hbond_count,
// larsen_hbond_water_term, bonded_energy_*.
// ──────────────────────────────────────────────────────────────────

struct QtScalarTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<double> data;  // (N*T,) row-major
    QtTimeSeriesFrameMeta meta;

    QString units;
    QString result_name;
    QString dataset_name;  // the per-TR dataset key (e.g. "sasa", "charge")

    double at(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return 0.0;
        return data[atomIdx * n_frames + t];
    }

    bool sourceAttachedAt(std::size_t t) const { return meta.sourceAttachedAt(t); }
};


// ──────────────────────────────────────────────────────────────────
// QtVec3TimeSeries — (N, T, 3) per-atom Cartesian Vec3 TS.
// Used by apbs_efield, tripeptide_*_residual_vec_*,
// aimnet2_charge_response_gradient (vec component).
// ──────────────────────────────────────────────────────────────────

struct QtVec3TimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<double> xyz;  // (N*T*3,) row-major
    QtTimeSeriesFrameMeta meta;

    QString units;
    QString result_name;
    QString irrep_layout;  // typically "x,y,z" or "v_x,v_y,v_z"

    Vec3 at(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return Vec3::Zero();
        const std::size_t base = (atomIdx * n_frames + t) * 3;
        return Vec3(xyz[base + 0], xyz[base + 1], xyz[base + 2]);
    }

    bool sourceAttachedAt(std::size_t t) const { return meta.sourceAttachedAt(t); }
};


// ──────────────────────────────────────────────────────────────────
// QtT2TimeSeries — (N, T, 5) per-atom T2-only tensor TS.
// Used by apbs_efg (EFG is symmetric-traceless, so T0/T1 are
// structural zeros and the writer projects only the 5 T2 components).
// ──────────────────────────────────────────────────────────────────

struct QtT2TimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<double> t2;  // (N*T*5,) row-major
    QtTimeSeriesFrameMeta meta;

    QString units;
    QString result_name;
    QString irrep_layout;  // "m-2,m-1,m0,m+1,m+2"

    std::array<double, 5> at(std::size_t atomIdx, std::size_t t) const {
        std::array<double, 5> out{};
        if (atomIdx >= n_atoms || t >= n_frames)
            return out;
        const std::size_t base = (atomIdx * n_frames + t) * 5;
        for (std::size_t k = 0; k < 5; ++k)
            out[k] = t2[base + k];
        return out;
    }

    bool sourceAttachedAt(std::size_t t) const { return meta.sourceAttachedAt(t); }
};


// ──────────────────────────────────────────────────────────────────
// QtTagTimeSeries — (N, T) per-atom uint8 tag.
// Used by tripeptide_bb_method_tag (encoding method-source per frame).
// ──────────────────────────────────────────────────────────────────

struct QtTagTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<uint8_t> tag;  // (N*T,) row-major
    QtTimeSeriesFrameMeta meta;

    QString units;  // typically "tag"
    QString result_name;
    QString dataset_name;

    uint8_t at(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return 0;
        return tag[atomIdx * n_frames + t];
    }
};


// ──────────────────────────────────────────────────────────────────
// QtEmbeddingTimeSeries — (N, T, 256) per-atom float32 embedding.
// Used by aimnet2_embedding only. Float32 storage to match the writer
// and avoid bloat — 256 doubles × 846 atoms × 751 frames = 1.3 GB,
// float32 halves that.
// ──────────────────────────────────────────────────────────────────

struct QtEmbeddingTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;
    std::size_t n_dims = 0;  // typically 256

    std::vector<float> data;  // (N*T*n_dims,) row-major
    QtTimeSeriesFrameMeta meta;

    QString result_name;

    // Returns a span into the embedding for one (atom, frame). View
    // is read-only; pointer valid for the buffer's lifetime.
    const float* dataAt(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return nullptr;
        return data.data() + (atomIdx * n_frames + t) * n_dims;
    }
};


// ──────────────────────────────────────────────────────────────────
// QtPositionsTimeSeries — (N, T, 3) per-atom Cartesian position TS.
// Always-attached (positions are present at tp.Seed time). The one
// TR that drives playback animation; every frame transition reads
// from this.
// ──────────────────────────────────────────────────────────────────

struct QtPositionsTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<double> xyz;  // (N*T*3,) row-major; Angstroms
    QtTimeSeriesFrameMeta meta;

    Vec3 at(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return Vec3::Zero();
        const std::size_t base = (atomIdx * n_frames + t) * 3;
        return Vec3(xyz[base + 0], xyz[base + 1], xyz[base + 2]);
    }
};


// ──────────────────────────────────────────────────────────────────
// QtChargeResponseGradientTimeSeries — (N, T) scalar +
// (N, T, 3) vec for the AIMNet2 charge-response gradient TR. The
// writer emits both `charge_response_gradient_scalar` and
// `charge_response_gradient_vector` in one group; we hold them
// alongside for the single-TR-equals-one-Qt-buffer discipline.
// ──────────────────────────────────────────────────────────────────

struct QtAimnet2ChargeResponseGradientTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<double> scalar;  // (N*T,)
    std::vector<double> vec;     // (N*T*3,)
    QtTimeSeriesFrameMeta meta;

    QString units;
    QString result_name;

    double scalarAt(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return 0.0;
        return scalar[atomIdx * n_frames + t];
    }
    Vec3 vecAt(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return Vec3::Zero();
        const std::size_t base = (atomIdx * n_frames + t) * 3;
        return Vec3(vec[base + 0], vec[base + 1], vec[base + 2]);
    }
};


}  // namespace h5reader::model

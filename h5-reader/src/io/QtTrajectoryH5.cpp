// QtTrajectoryH5 implementation — eager loader for the per-TR-group
// trajectory.h5 format. HighFive boundary lives entirely in this .cpp;
// callers see only the typed buffer accessors on the header.
//
// Pattern (mirrors ui/src/TrajectoryH5.cpp):
//   - File-local helpers: TryReadDataset / TryReadAttribute /
//     WarnShapeMismatch / TryReadFrameMeta.
//   - Per-shape buffer readers: one function per buffer struct,
//     parameterised by group path + dataset names.
//   - Per-Welford readers: one per Welford struct (channel handling).
//   - Special-TR readers: one per unique-shape TR (DSSP, dihedral,
//     gromacs_energy, ring_neighbourhood, rmsd_tracking, etc.).
//   - Constructor: walks /trajectory child groups, dispatches by name.
//
// Sparse-tolerant: absent groups leave the unique_ptr null. Group-
// present-but-malformed logs Warn and also leaves null.

#include "QtTrajectoryH5.h"

#include "../diagnostics/ErrorBus.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonParseError>

#include "highfive/H5Attribute.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5Group.hpp"

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <string>

namespace h5reader::io {

using namespace h5reader::model;

namespace {

// ── Diagnostics ────────────────────────────────────────────────────

void WarnShapeMismatch(const char* group_path, const QString& detail) {
    h5reader::diagnostics::ErrorBus::Report(
        h5reader::diagnostics::Severity::Warning,
        QStringLiteral("QtTrajectoryH5"),
        QStringLiteral("%1 present but malformed: %2").arg(QString::fromUtf8(group_path), detail),
        QString());
}

void WarnGroupAbsent(const char* /*group_path*/) {
    // Sparse-set normal — do not log per absent group (would spam at
    // every load for the 30+ TRs not present in any given run).
}

// Qt + HighFive exception boundary. Per Qt rules, exceptions must
// never escape back into the event-loop-managed caller — the
// QtTrajectoryH5 constructor runs from QtProteinLoader, which runs
// from the loader worker thread. Each Read* function calls this
// helper from its catch handlers so a single malformed TR group
// degrades to a logged warning + nullptr buffer (the existing
// "absent, not faked" semantics) instead of unwinding the whole
// load. `noexcept` because the handler itself must not throw — if
// WarnShapeMismatch ever did, the whole std::terminate would fire.
void LogReadException(const char* group_path, const char* kind, const std::exception& e) noexcept {
    try {
        WarnShapeMismatch(group_path,
            QStringLiteral("%1: %2").arg(QString::fromLatin1(kind),
                                         QString::fromUtf8(e.what())));
    } catch (...) {
        // Truly defensive; never reachable in practice.
    }
}

void LogReadException(const char* group_path, const char* kind) noexcept {
    try {
        WarnShapeMismatch(group_path,
            QStringLiteral("%1: unknown exception").arg(QString::fromLatin1(kind)));
    } catch (...) {
    }
}

void ValidateRequiredPositions(const QtPositionsTimeSeries* positions,
                               std::size_t expected_frames,
                               std::size_t expected_atoms) {
    if (!positions) {
        throw std::runtime_error(
            "QtTrajectoryH5: /trajectory/positions missing or malformed; positions buffer is null");
    }
    if (positions->n_frames != expected_frames) {
        throw std::runtime_error(
            QStringLiteral("QtTrajectoryH5: /trajectory/positions frame count mismatch: "
                           "positions n_frames=%1, /trajectory/frames count=%2")
                .arg(positions->n_frames)
                .arg(expected_frames)
                .toStdString());
    }
    if (positions->n_atoms != expected_atoms) {
        throw std::runtime_error(
            QStringLiteral("QtTrajectoryH5: /trajectory/positions atom count mismatch: "
                           "positions n_atoms=%1, /atoms count=%2")
                .arg(positions->n_atoms)
                .arg(expected_atoms)
                .toStdString());
    }
}

// ── Generic readers ────────────────────────────────────────────────

template <typename T>
bool TryReadDataset(const HighFive::Group& grp, const std::string& name, std::vector<T>& out) {
    if (!grp.exist(name))
        return false;
    grp.getDataSet(name).read(out);
    return true;
}

template <typename T>
bool TryReadAttribute(const HighFive::Group& grp, const std::string& name, T& out) {
    if (!grp.hasAttribute(name))
        return false;
    grp.getAttribute(name).read(out);
    return true;
}

bool TryReadAttributeQ(const HighFive::Group& grp, const std::string& name, QString& out) {
    if (!grp.hasAttribute(name))
        return false;
    std::string s;
    grp.getAttribute(name).read(s);
    out = QString::fromStdString(s);
    return true;
}

template <typename T>
bool TryReadFileAttribute(const HighFive::File& file, const std::string& name, T& out) {
    if (!file.hasAttribute(name))
        return false;
    file.getAttribute(name).read(out);
    return true;
}

bool TryReadFileAttributeQ(const HighFive::File& file, const std::string& name, QString& out) {
    if (!file.hasAttribute(name))
        return false;
    std::string s;
    file.getAttribute(name).read(s);
    out = QString::fromStdString(s);
    return true;
}

// Read a (T,) sized raw-byte dataset into a vector with size assumed.
// For HighFive's T* read overload to be safe, the buffer must be
// presized to the dataset's element count.
template <typename T>
bool ReadFlat(const HighFive::DataSet& ds, std::vector<T>& out, std::size_t expected_count) {
    if (out.size() != expected_count)
        out.resize(expected_count);
    ds.read(out.data());
    return true;
}

// Read /trajectory/<grp>/frame_indices + frame_times + optional
// source_attached_per_frame into a QtTimeSeriesFrameMeta or
// QtPerResidueFrameMeta. Returns false if frame_times absent or wrong size.
template <typename FrameMetaT>
bool TryReadFrameMeta(const HighFive::Group& grp, FrameMetaT& meta, std::size_t expected_T, const char* group_path) {
    bool ok = true;
    if (grp.exist("frame_times")) {
        meta.frame_times.clear();
        grp.getDataSet("frame_times").read(meta.frame_times);
        if (meta.frame_times.size() != expected_T) {
            WarnShapeMismatch(
                group_path,
                QStringLiteral("frame_times size %1 != expected_T %2").arg(meta.frame_times.size()).arg(expected_T));
            ok = false;
        }
    } else {
        WarnShapeMismatch(group_path, QStringLiteral("missing frame_times"));
        ok = false;
    }
    if (grp.exist("frame_indices")) {
        meta.frame_indices.clear();
        grp.getDataSet("frame_indices").read(meta.frame_indices);
    }
    if (grp.exist("source_attached_per_frame")) {
        meta.source_attached.clear();
        grp.getDataSet("source_attached_per_frame").read(meta.source_attached);
    }
    return ok;
}

// ── Per-shape buffer readers ──────────────────────────────────────

void ReadShieldingTimeSeries(HighFive::File& file,
                             const char* group_path,
                             std::size_t n_atoms,
                             std::unique_ptr<QtShieldingTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("xyz")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing xyz"));
        return;
    }
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getDimensions();
    if (dims.size() != 3 || dims[0] != n_atoms || dims[2] != 9) {
        WarnShapeMismatch(group_path, QStringLiteral("xyz shape != [n_atoms=%1, T, 9]").arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtShieldingTimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    ReadFlat<double>(ds, buf->xyz, n_atoms * buf->n_frames * 9);
    if (!TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path)) {
        // Frame meta missing is a soft fail; the data is still usable
        // for per-atom-slice display but per-frame indexing is degraded.
    }
    TryReadAttributeQ(grp, "irrep_layout", buf->irrep_layout);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "parity", buf->parity);
    TryReadAttributeQ(grp, "normalization", buf->normalization);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadScalarTimeSeries(HighFive::File& file,
                          const char* group_path,
                          const char* dataset_name,
                          std::size_t n_atoms,
                          std::unique_ptr<QtScalarTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist(dataset_name)) {
        WarnShapeMismatch(group_path, QStringLiteral("missing dataset %1").arg(QString::fromUtf8(dataset_name)));
        return;
    }
    auto ds = grp.getDataSet(dataset_name);
    const auto dims = ds.getDimensions();
    if (dims.size() != 2 || dims[0] != n_atoms) {
        WarnShapeMismatch(group_path,
                          QStringLiteral("%1 shape != [n_atoms=%2, T]").arg(QString::fromUtf8(dataset_name)).arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtScalarTimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    buf->dataset_name = QString::fromUtf8(dataset_name);
    ReadFlat<double>(ds, buf->data, n_atoms * buf->n_frames);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadVec3TimeSeries(HighFive::File& file,
                        const char* group_path,
                        const char* dataset_name,
                        std::size_t n_atoms,
                        std::unique_ptr<QtVec3TimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist(dataset_name)) {
        WarnShapeMismatch(group_path, QStringLiteral("missing dataset %1").arg(QString::fromUtf8(dataset_name)));
        return;
    }
    auto ds = grp.getDataSet(dataset_name);
    const auto dims = ds.getDimensions();
    if (dims.size() != 3 || dims[0] != n_atoms || dims[2] != 3) {
        WarnShapeMismatch(group_path,
                          QStringLiteral("%1 shape != [n_atoms=%2, T, 3]").arg(QString::fromUtf8(dataset_name)).arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtVec3TimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    ReadFlat<double>(ds, buf->xyz, n_atoms * buf->n_frames * 3);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttributeQ(grp, "irrep_layout", buf->irrep_layout);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadT2TimeSeries(HighFive::File& file,
                      const char* group_path,
                      const char* dataset_name,
                      std::size_t n_atoms,
                      std::unique_ptr<QtT2TimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist(dataset_name)) {
        WarnShapeMismatch(group_path, QStringLiteral("missing dataset %1").arg(QString::fromUtf8(dataset_name)));
        return;
    }
    auto ds = grp.getDataSet(dataset_name);
    const auto dims = ds.getDimensions();
    if (dims.size() != 3 || dims[0] != n_atoms || dims[2] != 5) {
        WarnShapeMismatch(group_path,
                          QStringLiteral("%1 shape != [n_atoms=%2, T, 5]").arg(QString::fromUtf8(dataset_name)).arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtT2TimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    ReadFlat<double>(ds, buf->t2, n_atoms * buf->n_frames * 5);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttributeQ(grp, "irrep_layout", buf->irrep_layout);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadTagTimeSeries(HighFive::File& file,
                       const char* group_path,
                       const char* dataset_name,
                       std::size_t n_atoms,
                       std::unique_ptr<QtTagTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist(dataset_name)) {
        WarnShapeMismatch(group_path, QStringLiteral("missing dataset %1").arg(QString::fromUtf8(dataset_name)));
        return;
    }
    auto ds = grp.getDataSet(dataset_name);
    const auto dims = ds.getDimensions();
    if (dims.size() != 2 || dims[0] != n_atoms) {
        WarnShapeMismatch(group_path,
                          QStringLiteral("%1 shape != [n_atoms=%2, T]").arg(QString::fromUtf8(dataset_name)).arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtTagTimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    buf->dataset_name = QString::fromUtf8(dataset_name);
    ReadFlat<uint8_t>(ds, buf->tag, n_atoms * buf->n_frames);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadEmbeddingTimeSeries(HighFive::File& file,
                             const char* group_path,
                             const char* dataset_name,
                             std::size_t n_atoms,
                             std::unique_ptr<QtEmbeddingTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist(dataset_name)) {
        WarnShapeMismatch(group_path, QStringLiteral("missing dataset %1").arg(QString::fromUtf8(dataset_name)));
        return;
    }
    auto ds = grp.getDataSet(dataset_name);
    const auto dims = ds.getDimensions();
    if (dims.size() != 3 || dims[0] != n_atoms) {
        WarnShapeMismatch(
            group_path,
            QStringLiteral("%1 shape != [n_atoms=%2, T, n_dims]").arg(QString::fromUtf8(dataset_name)).arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtEmbeddingTimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    buf->n_dims = dims[2];
    ReadFlat<float>(ds, buf->data, n_atoms * buf->n_frames * buf->n_dims);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadPositionsTimeSeries(HighFive::File& file,
                             const char* group_path,
                             std::size_t n_atoms,
                             std::unique_ptr<QtPositionsTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("xyz")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing xyz"));
        return;
    }
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getDimensions();
    if (dims.size() != 3 || dims[2] != 3) {
        WarnShapeMismatch(group_path, QStringLiteral("xyz shape != [n_atoms=%1, T, 3]").arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtPositionsTimeSeries>();
    buf->n_atoms = dims[0];
    buf->n_frames = dims[1];
    ReadFlat<double>(ds, buf->xyz, buf->n_atoms * buf->n_frames * 3);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadChargeResponseGradientTimeSeries(HighFive::File& file,
                                          const char* group_path,
                                          std::size_t n_atoms,
                                          std::unique_ptr<QtAimnet2ChargeResponseGradientTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    // Two datasets: charge_response_gradient_scalar (N, T), charge_response_gradient_vector (N, T, 3).
    if (!grp.exist("charge_response_gradient_scalar") || !grp.exist("charge_response_gradient_vector")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing charge_response_gradient_{scalar,vector}"));
        return;
    }
    auto ds_s = grp.getDataSet("charge_response_gradient_scalar");
    auto ds_v = grp.getDataSet("charge_response_gradient_vector");
    const auto dims_s = ds_s.getDimensions();
    const auto dims_v = ds_v.getDimensions();
    if (dims_s.size() != 2 || dims_s[0] != n_atoms || dims_v.size() != 3 || dims_v[0] != n_atoms || dims_v[2] != 3
        || dims_s[1] != dims_v[1]) {
        WarnShapeMismatch(group_path, QStringLiteral("shape mismatch on scalar/vector pair"));
        return;
    }
    auto buf = std::make_unique<QtAimnet2ChargeResponseGradientTimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims_s[1];
    ReadFlat<double>(ds_s, buf->scalar, n_atoms * buf->n_frames);
    ReadFlat<double>(ds_v, buf->vec, n_atoms * buf->n_frames * 3);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

// ── Welford readers ───────────────────────────────────────────────

// Reads a 7-stat block (mean/std/m2/min/max/min_frame/max_frame) for one
// channel into a vector<QtWelfordMoments>. Channel name prefix lets us
// share the function across the 4+ channels per Welford TR.
//
// Producer convention: scalar-source Welfords use `<channel>_<stat>`
// names (canonical) and some emit legacy `<stat>` aliases. The reader
// accepts either.
bool ReadWelfordChannelScalar(const HighFive::Group& grp,
                              const std::string& channel_prefix,
                              std::size_t n_atoms,
                              std::vector<QtWelfordMoments>& out,
                              bool accept_legacy_no_prefix = false) {
    auto try_read = [&](const std::string& name, std::vector<double>& v) -> bool {
        if (grp.exist(name)) {
            grp.getDataSet(name).read(v);
            return v.size() == n_atoms;
        }
        return false;
    };
    auto try_read_u64 = [&](const std::string& name, std::vector<uint64_t>& v) -> bool {
        if (grp.exist(name)) {
            grp.getDataSet(name).read(v);
            return v.size() == n_atoms;
        }
        return false;
    };

    std::vector<double> mean, std_, m2, mn, mx;
    std::vector<uint64_t> mn_f, mx_f;

    const std::string p = channel_prefix.empty() ? std::string() : channel_prefix + "_";

    const bool ok_mean = try_read(p + "mean", mean) || (accept_legacy_no_prefix && try_read("mean", mean));
    const bool ok_std = try_read(p + "std", std_) || (accept_legacy_no_prefix && try_read("std", std_));
    if (!ok_mean || !ok_std)
        return false;

    try_read(p + "m2", m2);
    try_read(p + "min", mn);
    try_read(p + "max", mx);
    try_read_u64(p + "min_frame", mn_f);
    try_read_u64(p + "max_frame", mx_f);

    out.resize(n_atoms);
    for (std::size_t i = 0; i < n_atoms; ++i) {
        out[i].mean = mean[i];
        out[i].std = std_[i];
        out[i].m2 = i < m2.size() ? m2[i] : 0.0;
        out[i].min = i < mn.size() ? mn[i] : 0.0;
        out[i].max = i < mx.size() ? mx[i] : 0.0;
        out[i].min_frame = i < mn_f.size() ? mn_f[i] : 0;
        out[i].max_frame = i < mx_f.size() ? mx_f[i] : 0;
    }
    return true;
}

// Per-component variant: reads (N, K) shape into vector<array<Moments, K>>.
template <std::size_t K>
bool ReadWelfordChannelMultiComp(const HighFive::Group& grp,
                                 const std::string& channel_prefix,
                                 std::size_t n_atoms,
                                 std::vector<std::array<QtWelfordMoments, K>>& out) {
    // HighFive's read(vector<T>) infers a 1D dataset; (N, K) 2D
    // datasets must use the raw-pointer overload with the buffer
    // pre-sized. See ReadFlat's pattern.
    auto try_read_2d = [&](const std::string& name, std::vector<double>& v) -> bool {
        if (!grp.exist(name))
            return false;
        auto ds = grp.getDataSet(name);
        const auto dims = ds.getDimensions();
        if (dims.size() != 2 || dims[0] != n_atoms || dims[1] != K)
            return false;
        v.resize(n_atoms * K);
        ds.read(v.data());
        return true;
    };
    auto try_read_2d_u64 = [&](const std::string& name, std::vector<uint64_t>& v) -> bool {
        if (!grp.exist(name))
            return false;
        auto ds = grp.getDataSet(name);
        const auto dims = ds.getDimensions();
        if (dims.size() != 2 || dims[0] != n_atoms || dims[1] != K)
            return false;
        v.resize(n_atoms * K);
        ds.read(v.data());
        return true;
    };

    std::vector<double> mean, std_, m2, mn, mx;
    std::vector<uint64_t> mn_f, mx_f;

    const std::string p = channel_prefix + "_";
    if (!try_read_2d(p + "mean", mean) || !try_read_2d(p + "std", std_))
        return false;
    try_read_2d(p + "m2", m2);
    try_read_2d(p + "min", mn);
    try_read_2d(p + "max", mx);
    try_read_2d_u64(p + "min_frame", mn_f);
    try_read_2d_u64(p + "max_frame", mx_f);

    out.resize(n_atoms);
    for (std::size_t i = 0; i < n_atoms; ++i) {
        for (std::size_t k = 0; k < K; ++k) {
            const std::size_t idx = i * K + k;
            auto& m = out[i][k];
            m.mean = mean[idx];
            m.std = std_[idx];
            m.m2 = idx < m2.size() ? m2[idx] : 0.0;
            m.min = idx < mn.size() ? mn[idx] : 0.0;
            m.max = idx < mx.size() ? mx[idx] : 0.0;
            m.min_frame = idx < mn_f.size() ? mn_f[idx] : 0;
            m.max_frame = idx < mx_f.size() ? mx_f[idx] : 0;
        }
    }
    return true;
}

void ReadShieldingWelford(HighFive::File& file,
                          const char* group_path,
                          std::size_t n_atoms,
                          std::unique_ptr<QtShieldingWelford>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    auto buf = std::make_unique<QtShieldingWelford>();
    buf->n_atoms = n_atoms;

    if (!ReadWelfordChannelScalar(grp, "t0", n_atoms, buf->t0)) {
        WarnShapeMismatch(group_path, QStringLiteral("t0 channel missing/malformed"));
        return;
    }
    ReadWelfordChannelMultiComp<3>(grp, "t1", n_atoms, buf->t1);
    ReadWelfordChannelMultiComp<5>(grp, "t2", n_atoms, buf->t2);
    // Accept t2magnitude (canonical) and legacy t2mag aliases.
    if (!ReadWelfordChannelScalar(grp, "t2magnitude", n_atoms, buf->t2magnitude)) {
        ReadWelfordChannelScalar(grp, "t2mag", n_atoms, buf->t2magnitude);
    }
    ReadWelfordChannelScalar(grp, "t0_delta", n_atoms, buf->t0_delta);
    ReadWelfordChannelScalar(grp, "t0_abs_delta", n_atoms, buf->t0_abs_delta);
    ReadWelfordChannelScalar(grp, "t0_delta_squared", n_atoms, buf->t0_delta_squared);
    ReadWelfordChannelScalar(grp, "t0_dxdt", n_atoms, buf->t0_dxdt);
    if (grp.exist("t0_rms_delta"))
        grp.getDataSet("t0_rms_delta").read(buf->t0_rms_delta);
    if (grp.exist("n_frames_per_atom"))
        grp.getDataSet("n_frames_per_atom").read(buf->n_frames_per_atom);
    if (grp.exist("delta_n_per_atom"))
        grp.getDataSet("delta_n_per_atom").read(buf->delta_n_per_atom);
    if (grp.exist("dxdt_n_per_atom"))
        grp.getDataSet("dxdt_n_per_atom").read(buf->dxdt_n_per_atom);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttributeQ(grp, "irrep_layout_t1", buf->irrep_layout_t1);
    TryReadAttributeQ(grp, "irrep_layout_t2", buf->irrep_layout_t2);
    int ddof_attr = 1;
    TryReadAttribute(grp, "ddof", ddof_attr);
    buf->ddof = ddof_attr;
    TryReadAttribute(grp, "mean_dt_ps", buf->mean_dt_ps);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadScalarWelford(HighFive::File& file,
                       const char* group_path,
                       const std::string& channel_name,
                       std::size_t n_atoms,
                       std::unique_ptr<QtScalarWelford>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    auto buf = std::make_unique<QtScalarWelford>();
    buf->n_atoms = n_atoms;
    buf->channel_name = QString::fromStdString(channel_name);

    // Try canonical <channel>_mean first; fall back to legacy `mean`.
    if (!ReadWelfordChannelScalar(grp,
                                  channel_name,
                                  n_atoms,
                                  buf->value,
                                  /*accept_legacy_no_prefix=*/true)) {
        WarnShapeMismatch(group_path, QStringLiteral("scalar Welford value channel missing"));
        return;
    }
    ReadWelfordChannelScalar(grp, channel_name + "_delta", n_atoms, buf->delta);
    ReadWelfordChannelScalar(grp, channel_name + "_abs_delta", n_atoms, buf->abs_delta);
    ReadWelfordChannelScalar(grp, channel_name + "_delta_squared", n_atoms, buf->delta_squared);
    ReadWelfordChannelScalar(grp, channel_name + "_dxdt", n_atoms, buf->dxdt);
    if (grp.exist("rms_delta"))
        grp.getDataSet("rms_delta").read(buf->rms_delta);
    if (grp.exist("n_frames_per_atom"))
        grp.getDataSet("n_frames_per_atom").read(buf->n_frames_per_atom);
    if (grp.exist("delta_n_per_atom"))
        grp.getDataSet("delta_n_per_atom").read(buf->delta_n_per_atom);
    if (grp.exist("dxdt_n_per_atom"))
        grp.getDataSet("dxdt_n_per_atom").read(buf->dxdt_n_per_atom);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    int ddof_attr = 1;
    TryReadAttribute(grp, "ddof", ddof_attr);
    buf->ddof = ddof_attr;
    TryReadAttribute(grp, "mean_dt_ps", buf->mean_dt_ps);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadVec3Welford(HighFive::File& file, const char* group_path, std::size_t n_atoms, std::unique_ptr<QtVec3Welford>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    auto buf = std::make_unique<QtVec3Welford>();
    buf->n_atoms = n_atoms;
    // Try common channel names: "value" / generic / "vector"
    if (!ReadWelfordChannelMultiComp<3>(grp, "value", n_atoms, buf->components)) {
        ReadWelfordChannelMultiComp<3>(grp, "vector", n_atoms, buf->components);
    }
    ReadWelfordChannelScalar(grp, "magnitude", n_atoms, buf->magnitude);
    if (grp.exist("n_frames_per_atom"))
        grp.getDataSet("n_frames_per_atom").read(buf->n_frames_per_atom);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttributeQ(grp, "irrep_layout", buf->irrep_layout);
    int ddof_attr = 1;
    TryReadAttribute(grp, "ddof", ddof_attr);
    buf->ddof = ddof_attr;
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadBondOrderWelford(HighFive::File& file, const char* group_path, std::unique_ptr<QtBondOrderWelford>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    // Bond-axis size discovered from the data shape (not n_atoms).
    // Writer emits `order_*` dataset names (not `bond_order_*`).
    if (!grp.exist("order_mean")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing order_mean"));
        return;
    }
    std::vector<double> tmp;
    grp.getDataSet("order_mean").read(tmp);
    const std::size_t n_bonds = tmp.size();
    auto buf = std::make_unique<QtBondOrderWelford>();
    buf->n_bonds = n_bonds;
    ReadWelfordChannelScalar(grp, "order", n_bonds, buf->bond_order, true);
    // Writer emits `n_per_bond` (canonical) + `n_total_per_bond` (legacy);
    // prefer canonical, fall back to legacy.
    if (grp.exist("n_per_bond"))
        grp.getDataSet("n_per_bond").read(buf->n_frames_per_bond);
    else if (grp.exist("n_frames_per_bond"))
        grp.getDataSet("n_frames_per_bond").read(buf->n_frames_per_bond);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    int ddof_attr = 1;
    TryReadAttribute(grp, "ddof", ddof_attr);
    buf->ddof = ddof_attr;
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadHydrationWelford(HighFive::File& file,
                          const char* group_path,
                          std::size_t n_atoms,
                          std::unique_ptr<QtHydrationWelford>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    auto buf = std::make_unique<QtHydrationWelford>();
    buf->n_atoms = n_atoms;
    // Hydration TRs vary in channel set; walk dataset names ending in
    // "_mean" and pull the corresponding 7-stat block per channel.
    const auto names = grp.listObjectNames();
    for (const auto& name : names) {
        if (name.size() > 5 && name.substr(name.size() - 5) == "_mean") {
            const std::string prefix = name.substr(0, name.size() - 5);
            QtHydrationChannel ch;
            ch.name = QString::fromStdString(prefix);
            if (ReadWelfordChannelScalar(grp, prefix, n_atoms, ch.moments)) {
                buf->channels.push_back(std::move(ch));
            }
        }
    }
    if (grp.exist("n_frames_per_atom"))
        grp.getDataSet("n_frames_per_atom").read(buf->n_frames_per_atom);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    int ddof_attr = 1;
    TryReadAttribute(grp, "ddof", ddof_attr);
    buf->ddof = ddof_attr;
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadIRedOrderParameters(HighFive::File& file,
                             const char* group_path,
                             std::unique_ptr<QtIRedOrderParameters>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("s2_ired") || !grp.exist("residue_index")
        || !grp.exist("n_atom") || !grp.exist("h_atom")) {
        WarnShapeMismatch(group_path,
            QStringLiteral("missing one of s2_ired / residue_index / n_atom / h_atom"));
        return;
    }
    auto s2_ds = grp.getDataSet("s2_ired");
    const auto s2_dims = s2_ds.getDimensions();
    if (s2_dims.size() != 1) {
        WarnShapeMismatch(group_path, QStringLiteral("s2_ired shape != (M,)"));
        return;
    }
    const std::size_t M = s2_dims[0];

    // Every identity dataset must match s2_ired's M; mismatched lengths
    // mean malformed/partial-write H5 and downstream lookups (rowFor,
    // SceneRevealOverlay) would index OOB. Bail loudly rather than
    // silently truncate. This guard is the template Phases D-G clone.
    auto same_M = [&](const char* name) -> bool {
        if (!grp.exist(name))
            return false;
        const auto d = grp.getDataSet(name).getDimensions();
        if (d.size() != 1 || d[0] != M) {
            WarnShapeMismatch(group_path,
                QStringLiteral("%1 shape != (%2,)")
                    .arg(QString::fromLatin1(name)).arg(M));
            return false;
        }
        return true;
    };
    if (!same_M("residue_index") || !same_M("n_atom") || !same_M("h_atom"))
        return;
    if (grp.exist("eigenvalues") && grp.getDataSet("eigenvalues").getDimensions().size() == 1
        && grp.getDataSet("eigenvalues").getDimensions()[0] != M) {
        WarnShapeMismatch(group_path,
            QStringLiteral("eigenvalues shape != (%1,)").arg(M));
        return;
    }

    auto buf = std::make_unique<QtIRedOrderParameters>();
    buf->identity.n_vectors = M;
    buf->identity.kind.assign(M, 1);  // all amide N-H per the producer
    s2_ds.read(buf->s2_ired);
    if (grp.exist("eigenvalues"))
        grp.getDataSet("eigenvalues").read(buf->eigenvalues);
    grp.getDataSet("residue_index").read(buf->identity.residue_index);
    grp.getDataSet("n_atom").read(buf->identity.tail_atom);
    grp.getDataSet("h_atom").read(buf->identity.head_atom);
    buf->identity.owning_atom = buf->identity.head_atom;   // highlight = H

    if (grp.hasAttribute("separability_gap"))
        grp.getAttribute("separability_gap").read(buf->separability_gap);
    std::size_t n_frames = 0;
    if (grp.hasAttribute("n_frames")) {
        grp.getAttribute("n_frames").read(n_frames);
        buf->n_frames = n_frames;
    }
    TryReadAttributeQ(grp, "reference", buf->reference);
    TryReadAttributeQ(grp, "frame", buf->frame);
    TryReadAttributeQ(grp, "vector_set", buf->vector_set);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadKernelCoherence(HighFive::File& file,
                         const char* group_path,
                         std::size_t n_atoms,
                         std::unique_ptr<QtKernelCoherence>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("correlation_matrix") || !grp.exist("channel_names")) {
        WarnShapeMismatch(group_path,
            QStringLiteral("missing correlation_matrix or channel_names"));
        return;
    }
    const auto dims = grp.getDataSet("correlation_matrix").getDimensions();
    if (dims.size() != 3 || dims[0] != n_atoms || dims[1] != dims[2]) {
        WarnShapeMismatch(group_path,
            QStringLiteral("correlation_matrix shape != (n_atoms=%1, C, C)").arg(n_atoms));
        return;
    }
    const std::size_t C = dims[1];

    auto buf = std::make_unique<QtKernelCoherence>();
    buf->matrix.n_atoms = n_atoms;
    buf->matrix.n_channels = C;
    ReadFlat<double>(grp.getDataSet("correlation_matrix"),
                     buf->matrix.data, n_atoms * C * C);

    std::vector<std::string> names, units;
    grp.getDataSet("channel_names").read(names);
    if (grp.exist("channel_units"))
        grp.getDataSet("channel_units").read(units);
    if (names.size() != C) {
        WarnShapeMismatch(group_path,
            QStringLiteral("channel_names length != C=%1").arg(C));
        return;
    }
    for (const auto& s : names) buf->matrix.channel_names.push_back(QString::fromStdString(s));
    for (const auto& s : units) buf->matrix.channel_units.push_back(QString::fromStdString(s));
    buf->matrix.units = QStringLiteral("Pearson r (dimensionless)");
    TryReadAttributeQ(grp, "result_name", buf->matrix.result_name);

    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadDihedralAutocorrelation(HighFive::File& file,
                                 const char* group_path,
                                 std::unique_ptr<QtDihedralAutocorrelation>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    for (const char* req : {"phi_acf", "psi_acf",
                            "phi_corr_time_ps", "psi_corr_time_ps",
                            "lag_times_ps"}) {
        if (!grp.exist(req)) {
            WarnShapeMismatch(group_path,
                QStringLiteral("missing %1").arg(QString::fromLatin1(req)));
            return;
        }
    }
    const auto phi_dims = grp.getDataSet("phi_acf").getDimensions();
    if (phi_dims.size() != 2 || phi_dims[0] == 0) {
        WarnShapeMismatch(group_path, QStringLiteral("phi_acf shape != (R, L)"));
        return;
    }
    const std::size_t R = phi_dims[0];
    const std::size_t L = phi_dims[1];

    const auto psi_dims = grp.getDataSet("psi_acf").getDimensions();
    if (psi_dims.size() != 2 || psi_dims[0] != R || psi_dims[1] != L) {
        WarnShapeMismatch(group_path,
            QStringLiteral("psi_acf shape != (R=%1, L=%2)").arg(R).arg(L));
        return;
    }
    auto check_R = [&](const char* name) {
        const auto d = grp.getDataSet(name).getDimensions();
        if (d.size() != 1 || d[0] != R) {
            WarnShapeMismatch(group_path,
                QStringLiteral("%1 shape != (R=%2,)")
                    .arg(QString::fromLatin1(name)).arg(R));
            return false;
        }
        return true;
    };
    if (!check_R("phi_corr_time_ps") || !check_R("psi_corr_time_ps"))
        return;

    auto buf = std::make_unique<QtDihedralAutocorrelation>();
    buf->n_residues = R;
    buf->n_lags = L;

    auto fill_curve = [&](const char* name, QtPerResidueCurve& dst) {
        dst.n_residues = R;
        dst.n_samples = L;
        ReadFlat<double>(grp.getDataSet(name), dst.data, R * L);
        grp.getDataSet("lag_times_ps").read(dst.axis_values);
        dst.axis_unit = QStringLiteral("ps");
        dst.units = QStringLiteral("dimensionless");
    };
    fill_curve("phi_acf", buf->phi_acf);
    fill_curve("psi_acf", buf->psi_acf);

    auto fill_scalar = [&](const char* values_name, const char* defined_name,
                           QtPerResidueScalar& dst) {
        dst.n_residues = R;
        grp.getDataSet(values_name).read(dst.values);
        if (grp.exist(defined_name))
            grp.getDataSet(defined_name).read(dst.defined);
        dst.units = QStringLiteral("ps");
    };
    fill_scalar("phi_corr_time_ps", "phi_defined", buf->phi_corr_time);
    fill_scalar("psi_corr_time_ps", "psi_defined", buf->psi_corr_time);

    // Chi[0..3] composite payload — L-2a (2026-05-29). Producer emits
    // chi_acf (R, 4, L) + chi_corr_time_ps (R, 4) + chi_defined (R, 4).
    // Gracefully absent on older runs that predate the chi expansion;
    // we just leave the chi_* vectors empty in that case.
    if (grp.exist("chi_acf") && grp.exist("chi_corr_time_ps")) {
        const auto chi_dims = grp.getDataSet("chi_acf").getDimensions();
        if (chi_dims.size() == 3 && chi_dims[0] == R && chi_dims[1] == 4 && chi_dims[2] == L) {
            ReadFlat<double>(grp.getDataSet("chi_acf"), buf->chi_acf, R * 4 * L);
            buf->chi_acf_axis = buf->phi_acf.axis_values;  // shared lag grid
            const auto ct_dims = grp.getDataSet("chi_corr_time_ps").getDimensions();
            if (ct_dims.size() == 2 && ct_dims[0] == R && ct_dims[1] == 4) {
                ReadFlat<double>(grp.getDataSet("chi_corr_time_ps"), buf->chi_corr_time, R * 4);
            } else {
                WarnShapeMismatch(group_path,
                    QStringLiteral("chi_corr_time_ps shape != (R=%1, 4)").arg(R));
            }
            if (grp.exist("chi_defined")) {
                const auto d_dims = grp.getDataSet("chi_defined").getDimensions();
                if (d_dims.size() == 2 && d_dims[0] == R && d_dims[1] == 4) {
                    ReadFlat<uint8_t>(grp.getDataSet("chi_defined"), buf->chi_defined, R * 4);
                }
            }
        } else {
            WarnShapeMismatch(group_path,
                QStringLiteral("chi_acf shape != (R=%1, 4, L=%2)").arg(R).arg(L));
        }
    }

    if (grp.hasAttribute("sample_interval_ps"))
        grp.getAttribute("sample_interval_ps").read(buf->sample_interval_ps);

    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadReorientationalDynamics(HighFive::File& file,
                                 const char* group_path,
                                 std::unique_ptr<QtReorientationalDynamics>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    for (const char* req : {"bond_vector_autocorrelation",
                            "bond_vector_autocorrelation_lab",
                            "order_parameter_S2",
                            "lipari_szabo_tau_e",
                            "bond_orientation_tensor",
                            "vector_kind",
                            "owning_atom", "tail_atom", "head_atom",
                            "residue_index",
                            "lag_times_ps",
                            "spectral_density_j",
                            "relaxation_R1", "relaxation_R2",
                            "relaxation_NOE",
                            "relaxation_larmor_freqs_rad_per_s"}) {
        if (!grp.exist(req)) {
            WarnShapeMismatch(group_path,
                QStringLiteral("missing %1").arg(QString::fromLatin1(req)));
            return;
        }
    }

    // V = number of vectors; L = lag count.
    const auto s2_dims = grp.getDataSet("order_parameter_S2").getDimensions();
    if (s2_dims.size() != 1 || s2_dims[0] == 0) {
        WarnShapeMismatch(group_path, QStringLiteral("order_parameter_S2 shape != (V,)"));
        return;
    }
    const std::size_t V = s2_dims[0];

    const auto tcf_dims = grp.getDataSet("bond_vector_autocorrelation").getDimensions();
    if (tcf_dims.size() != 2 || tcf_dims[0] != V) {
        WarnShapeMismatch(group_path,
            QStringLiteral("bond_vector_autocorrelation shape != (V=%1, L)").arg(V));
        return;
    }
    const std::size_t L = tcf_dims[1];

    auto check_V = [&](const char* name) {
        const auto d = grp.getDataSet(name).getDimensions();
        if (d.size() != 1 || d[0] != V) {
            WarnShapeMismatch(group_path,
                QStringLiteral("%1 shape != (V=%2,)")
                    .arg(QString::fromLatin1(name)).arg(V));
            return false;
        }
        return true;
    };
    if (!check_V("lipari_szabo_tau_e") || !check_V("vector_kind")
        || !check_V("owning_atom") || !check_V("tail_atom")
        || !check_V("head_atom") || !check_V("residue_index")
        || !check_V("relaxation_R1") || !check_V("relaxation_R2")
        || !check_V("relaxation_NOE")) {
        return;
    }

    const auto tens_dims = grp.getDataSet("bond_orientation_tensor").getDimensions();
    if (tens_dims.size() != 3 || tens_dims[0] != V || tens_dims[1] != 3 || tens_dims[2] != 3) {
        WarnShapeMismatch(group_path,
            QStringLiteral("bond_orientation_tensor shape != (V=%1, 3, 3)").arg(V));
        return;
    }

    const auto j_dims = grp.getDataSet("spectral_density_j").getDimensions();
    if (j_dims.size() != 2 || j_dims[0] != V) {
        WarnShapeMismatch(group_path,
            QStringLiteral("spectral_density_j shape != (V=%1, K)").arg(V));
        return;
    }
    const std::size_t K = j_dims[1];

    auto buf = std::make_unique<QtReorientationalDynamics>();

    // Identity table
    buf->identity.n_vectors = V;
    grp.getDataSet("vector_kind").read(buf->identity.kind);
    grp.getDataSet("residue_index").read(buf->identity.residue_index);
    grp.getDataSet("owning_atom").read(buf->identity.owning_atom);
    grp.getDataSet("tail_atom").read(buf->identity.tail_atom);
    grp.getDataSet("head_atom").read(buf->identity.head_atom);

    // Curves (body + lab)
    auto fill_curve = [&](const char* name, QtPerBondVectorCurve& dst) {
        dst.n_vectors = V;
        dst.n_samples = L;
        ReadFlat<double>(grp.getDataSet(name), dst.data, V * L);
        grp.getDataSet("lag_times_ps").read(dst.axis_values);
        dst.axis_unit = QStringLiteral("ps");
        dst.units = QStringLiteral("dimensionless");
    };
    fill_curve("bond_vector_autocorrelation",     buf->acf_internal);
    fill_curve("bond_vector_autocorrelation_lab", buf->acf_lab);

    // Per-vector scalars (S², τ_e, R1, R2, NOE)
    auto fill_scalar = [&](const char* name, QtPerBondVectorScalar& dst, const char* unit) {
        dst.n_vectors = V;
        grp.getDataSet(name).read(dst.values);
        dst.units = QString::fromLatin1(unit);
    };
    fill_scalar("order_parameter_S2",  buf->s2,    "dimensionless");
    fill_scalar("lipari_szabo_tau_e",  buf->tau_e, "ps");
    fill_scalar("relaxation_R1",       buf->r1,    "1/s");
    fill_scalar("relaxation_R2",       buf->r2,    "1/s");
    fill_scalar("relaxation_NOE",      buf->noe,   "dimensionless");

    // Mat3 orientation tensor
    buf->orientation_tensor.n_vectors = V;
    ReadFlat<double>(grp.getDataSet("bond_orientation_tensor"),
                     buf->orientation_tensor.data, V * 9);

    // Spectral-density J at the K KTB Larmor frequencies
    buf->spectral_density_J.n_vectors = V;
    buf->spectral_density_J.n_freqs = K;
    ReadFlat<double>(grp.getDataSet("spectral_density_j"),
                     buf->spectral_density_J.data, V * K);
    grp.getDataSet("relaxation_larmor_freqs_rad_per_s").read(buf->spectral_density_J.freq_values);
    buf->spectral_density_J.units = QStringLiteral("s");

    // Per-trajectory attrs
    if (grp.hasAttribute("tau_m_ps")) grp.getAttribute("tau_m_ps").read(buf->tau_m_ps);
    if (grp.hasAttribute("trajectory_length_over_tau_m"))
        grp.getAttribute("trajectory_length_over_tau_m").read(buf->trajectory_length_over_tau_m);
    if (grp.hasAttribute("tau_m_converged"))
        grp.getAttribute("tau_m_converged").read(buf->tau_m_converged);
    if (grp.hasAttribute("relaxation_field_tesla"))
        grp.getAttribute("relaxation_field_tesla").read(buf->relaxation_field_tesla);
    TryReadAttributeQ(grp, "result_name", buf->result_name);

    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadKernelDynamics(HighFive::File& file,
                        const char* group_path,
                        std::size_t n_atoms,
                        std::unique_ptr<QtKernelDynamics>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    for (const char* req : {"acf", "power_spectrum", "decay_time_ps",
                            "peak_freq_per_ps", "spectral_centroid_per_ps",
                            "channel_names", "channel_units",
                            "lag_times_ps", "frequencies_per_ps"}) {
        if (!grp.exist(req)) {
            WarnShapeMismatch(group_path,
                QStringLiteral("missing %1").arg(QString::fromLatin1(req)));
            return;
        }
    }

    const auto acf_dims = grp.getDataSet("acf").getDimensions();
    if (acf_dims.size() != 3 || acf_dims[0] != n_atoms) {
        WarnShapeMismatch(group_path,
            QStringLiteral("acf shape != (n_atoms=%1, n_channels, n_lags)").arg(n_atoms));
        return;
    }
    const std::size_t N = acf_dims[0];
    const std::size_t C = acf_dims[1];
    const std::size_t L = acf_dims[2];

    const auto psd_dims = grp.getDataSet("power_spectrum").getDimensions();
    if (psd_dims.size() != 3 || psd_dims[0] != N || psd_dims[1] != C) {
        WarnShapeMismatch(group_path,
            QStringLiteral("power_spectrum shape != (N=%1, C=%2, F)").arg(N).arg(C));
        return;
    }
    const std::size_t F = psd_dims[2];

    auto check_NC = [&](const char* name) {
        const auto d = grp.getDataSet(name).getDimensions();
        if (d.size() != 2 || d[0] != N || d[1] != C) {
            WarnShapeMismatch(group_path,
                QStringLiteral("%1 shape != (N=%2, C=%3)")
                    .arg(QString::fromLatin1(name)).arg(N).arg(C));
            return false;
        }
        return true;
    };
    if (!check_NC("decay_time_ps") || !check_NC("peak_freq_per_ps")
        || !check_NC("spectral_centroid_per_ps")) {
        return;
    }

    auto buf = std::make_unique<QtKernelDynamics>();
    buf->n_atoms = N;
    buf->n_channels = C;

    // Channel metadata
    {
        std::vector<std::string> names, units;
        grp.getDataSet("channel_names").read(names);
        grp.getDataSet("channel_units").read(units);
        if (names.size() != C || units.size() != C) {
            WarnShapeMismatch(group_path,
                QStringLiteral("channel_names/channel_units length != C=%1").arg(C));
            return;
        }
        for (const auto& s : names) buf->channel_names.push_back(QString::fromStdString(s));
        for (const auto& s : units) buf->channel_units.push_back(QString::fromStdString(s));
    }
    if (grp.hasAttribute("sample_interval_ps"))
        grp.getAttribute("sample_interval_ps").read(buf->sample_interval_ps);

    // ACF curves
    buf->acf.n_atoms = N;
    buf->acf.n_channels = C;
    buf->acf.n_samples = L;
    ReadFlat<double>(grp.getDataSet("acf"), buf->acf.data, N * C * L);
    grp.getDataSet("lag_times_ps").read(buf->acf.axis_values);
    buf->acf.axis_unit = QStringLiteral("ps");
    buf->acf.axis_label = QStringLiteral("lag");
    buf->acf.units = QStringLiteral("dimensionless");
    buf->acf.channel_names = buf->channel_names;
    buf->acf.channel_units = buf->channel_units;

    // PSD curves
    buf->power_spectrum.n_atoms = N;
    buf->power_spectrum.n_channels = C;
    buf->power_spectrum.n_samples = F;
    ReadFlat<double>(grp.getDataSet("power_spectrum"), buf->power_spectrum.data, N * C * F);
    grp.getDataSet("frequencies_per_ps").read(buf->power_spectrum.axis_values);
    buf->power_spectrum.axis_unit = QStringLiteral("1/ps");
    buf->power_spectrum.axis_label = QStringLiteral("frequency");
    buf->power_spectrum.units = QStringLiteral("channel_units^2*ps");
    buf->power_spectrum.channel_names = buf->channel_names;
    buf->power_spectrum.channel_units = buf->channel_units;

    // Scalar reductions
    auto fill_scalar = [&](const char* name, QtPerAtomChannelScalar& dst, const char* unit) {
        dst.n_atoms = N;
        dst.n_channels = C;
        ReadFlat<double>(grp.getDataSet(name), dst.data, N * C);
        dst.units = QString::fromLatin1(unit);
        dst.channel_names = buf->channel_names;
        dst.channel_units = buf->channel_units;
    };
    fill_scalar("decay_time_ps",            buf->decay_time_ps,            "ps");
    fill_scalar("peak_freq_per_ps",         buf->peak_freq_per_ps,         "1/ps");
    fill_scalar("spectral_centroid_per_ps", buf->spectral_centroid_per_ps, "1/ps");

    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadAutocorrelation(HighFive::File& file,
                         const char* group_path,
                         std::size_t n_atoms,
                         std::unique_ptr<QtAutocorrelation>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("rho")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing rho"));
        return;
    }
    auto ds = grp.getDataSet("rho");
    const auto dims = ds.getDimensions();
    if (dims.size() != 2 || dims[0] != n_atoms) {
        WarnShapeMismatch(group_path, QStringLiteral("rho shape != [n_atoms=%1, n_lags]").arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtAutocorrelation>();
    buf->n_atoms = n_atoms;
    buf->n_lags = dims[1];
    ReadFlat<double>(ds, buf->rho, n_atoms * buf->n_lags);
    if (grp.exist("lag_times_ps"))
        grp.getDataSet("lag_times_ps").read(buf->lag_times_ps);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

// ── Composite-shape TR readers ────────────────────────────────────

void ReadWaterFieldTimeSeries(HighFive::File& file,
                              const char* group_path,
                              std::size_t n_atoms,
                              std::unique_ptr<QtWaterFieldTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("efg") || !grp.exist("efield")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing efg or efield"));
        return;
    }
    auto ds_efg = grp.getDataSet("efg");
    const auto dims = ds_efg.getDimensions();
    if (dims.size() != 3 || dims[0] != n_atoms || dims[2] != 5) {
        WarnShapeMismatch(group_path, QStringLiteral("efg shape != [n_atoms=%1, T, 5]").arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtWaterFieldTimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    ReadFlat<double>(ds_efg, buf->efg, n_atoms * buf->n_frames * 5);
    if (grp.exist("efg_first"))
        ReadFlat<double>(grp.getDataSet("efg_first"), buf->efg_first, n_atoms * buf->n_frames * 5);
    ReadFlat<double>(grp.getDataSet("efield"), buf->efield, n_atoms * buf->n_frames * 3);
    if (grp.exist("efield_first"))
        ReadFlat<double>(grp.getDataSet("efield_first"), buf->efield_first, n_atoms * buf->n_frames * 3);
    if (grp.exist("n_first"))
        ReadFlat<uint32_t>(grp.getDataSet("n_first"), buf->n_first, n_atoms * buf->n_frames);
    if (grp.exist("n_second"))
        ReadFlat<uint32_t>(grp.getDataSet("n_second"), buf->n_second, n_atoms * buf->n_frames);
    if (grp.exist("frame_indices"))
        grp.getDataSet("frame_indices").read(buf->frame_indices);
    if (grp.exist("frame_times"))
        grp.getDataSet("frame_times").read(buf->frame_times);
    if (grp.exist("source_attached_per_frame"))
        grp.getDataSet("source_attached_per_frame").read(buf->source_attached);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadHydrationShellTimeSeries(HighFive::File& file,
                                  const char* group_path,
                                  std::size_t n_atoms,
                                  std::unique_ptr<QtHydrationShellTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("half_shell_asymmetry")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing half_shell_asymmetry"));
        return;
    }
    auto ds = grp.getDataSet("half_shell_asymmetry");
    const auto dims = ds.getDimensions();
    if (dims.size() != 2 || dims[0] != n_atoms) {
        WarnShapeMismatch(group_path, QStringLiteral("half_shell_asymmetry shape != [n_atoms=%1, T]").arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtHydrationShellTimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    ReadFlat<double>(ds, buf->half_shell_asymmetry, n_atoms * buf->n_frames);
    if (grp.exist("mean_water_dipole_cos"))
        ReadFlat<double>(grp.getDataSet("mean_water_dipole_cos"), buf->mean_water_dipole_cos, n_atoms * buf->n_frames);
    if (grp.exist("nearest_ion_charge"))
        ReadFlat<double>(grp.getDataSet("nearest_ion_charge"), buf->nearest_ion_charge, n_atoms * buf->n_frames);
    if (grp.exist("nearest_ion_distance"))
        ReadFlat<double>(grp.getDataSet("nearest_ion_distance"), buf->nearest_ion_distance, n_atoms * buf->n_frames);
    if (grp.exist("frame_indices"))
        grp.getDataSet("frame_indices").read(buf->frame_indices);
    if (grp.exist("frame_times"))
        grp.getDataSet("frame_times").read(buf->frame_times);
    if (grp.exist("source_attached_per_frame"))
        grp.getDataSet("source_attached_per_frame").read(buf->source_attached);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadHydrationGeometryTimeSeries(HighFive::File& file,
                                     const char* group_path,
                                     std::size_t n_atoms,
                                     std::unique_ptr<QtHydrationGeometryTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("dipole_alignment")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing dipole_alignment"));
        return;
    }
    auto ds = grp.getDataSet("dipole_alignment");
    const auto dims = ds.getDimensions();
    if (dims.size() != 2 || dims[0] != n_atoms) {
        WarnShapeMismatch(group_path, QStringLiteral("dipole_alignment shape != [n_atoms=%1, T]").arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtHydrationGeometryTimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    ReadFlat<double>(ds, buf->dipole_alignment, n_atoms * buf->n_frames);
    if (grp.exist("dipole_coherence"))
        ReadFlat<double>(grp.getDataSet("dipole_coherence"), buf->dipole_coherence, n_atoms * buf->n_frames);
    if (grp.exist("dipole_vector"))
        ReadFlat<double>(grp.getDataSet("dipole_vector"), buf->dipole_vector, n_atoms * buf->n_frames * 3);
    if (grp.exist("first_shell_count"))
        ReadFlat<uint32_t>(grp.getDataSet("first_shell_count"), buf->first_shell_count, n_atoms * buf->n_frames);
    if (grp.exist("half_shell_asymmetry"))
        ReadFlat<double>(grp.getDataSet("half_shell_asymmetry"), buf->half_shell_asymmetry, n_atoms * buf->n_frames);
    if (grp.exist("surface_normal"))
        ReadFlat<double>(grp.getDataSet("surface_normal"), buf->surface_normal, n_atoms * buf->n_frames * 3);
    if (grp.exist("frame_indices"))
        grp.getDataSet("frame_indices").read(buf->frame_indices);
    if (grp.exist("frame_times"))
        grp.getDataSet("frame_times").read(buf->frame_times);
    if (grp.exist("source_attached_per_frame"))
        grp.getDataSet("source_attached_per_frame").read(buf->source_attached);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadRingPuckerTimeSeries(HighFive::File& file, const char* group_path, std::unique_ptr<QtRingPuckerTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("aromatic_chi2")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing aromatic_chi2"));
        return;
    }
    auto ds_a = grp.getDataSet("aromatic_chi2");
    const auto dims_a = ds_a.getDimensions();
    if (dims_a.size() != 2) {
        WarnShapeMismatch(group_path, QStringLiteral("aromatic_chi2 rank != 2"));
        return;
    }
    auto buf = std::make_unique<QtRingPuckerTimeSeries>();
    buf->n_aromatic_rings = dims_a[0];
    buf->n_frames = dims_a[1];
    ReadFlat<double>(ds_a, buf->aromatic_chi2, buf->n_aromatic_rings * buf->n_frames);
    if (grp.exist("aromatic_parent_residue_index"))
        grp.getDataSet("aromatic_parent_residue_index").read(buf->aromatic_parent_residue_index);
    if (grp.exist("pucker_Q")) {
        auto ds_q = grp.getDataSet("pucker_Q");
        const auto dims_q = ds_q.getDimensions();
        if (dims_q.size() == 2 && dims_q[1] == buf->n_frames) {
            buf->n_saturated_rings = dims_q[0];
            ReadFlat<double>(ds_q, buf->pucker_Q, buf->n_saturated_rings * buf->n_frames);
            if (grp.exist("pucker_theta"))
                ReadFlat<double>(grp.getDataSet("pucker_theta"), buf->pucker_theta, buf->n_saturated_rings * buf->n_frames);
            if (grp.exist("saturated_parent_residue_index"))
                grp.getDataSet("saturated_parent_residue_index").read(buf->saturated_parent_residue_index);
        }
    }
    if (grp.exist("frame_indices"))
        grp.getDataSet("frame_indices").read(buf->frame_indices);
    if (grp.exist("frame_times"))
        grp.getDataSet("frame_times").read(buf->frame_times);
    if (grp.exist("source_attached_per_frame"))
        grp.getDataSet("source_attached_per_frame").read(buf->source_attached);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadJCouplingTimeSeries(HighFive::File& file,
                             const char* group_path,
                             std::size_t n_residues,
                             std::size_t n_atoms,
                             std::unique_ptr<QtJCouplingTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    auto buf = std::make_unique<QtJCouplingTimeSeries>();
    buf->n_residues = n_residues;
    auto try_J = [&](const char* name, std::vector<double>& v) {
        if (!grp.exist(name))
            return;
        auto ds = grp.getDataSet(name);
        const auto dims = ds.getDimensions();
        if (dims.size() != 2 || dims[0] != n_residues)
            return;
        if (buf->n_frames == 0)
            buf->n_frames = dims[1];
        ReadFlat<double>(ds, v, n_residues * buf->n_frames);
    };
    auto try_exist = [&](const char* name, std::vector<uint8_t>& v) {
        if (grp.exist(name))
            grp.getDataSet(name).read(v);
    };
    try_J("J_Cprime_Cgamma", buf->J_Cprime_Cgamma);
    try_J("J_HN_Cbeta", buf->J_HN_Cbeta);
    try_J("J_HN_Cprime", buf->J_HN_Cprime);
    try_J("J_HN_Halpha", buf->J_HN_Halpha);
    try_J("J_HN_Halpha_Vogeli", buf->J_HN_Halpha_Vogeli);
    try_J("J_Halpha_Cprime", buf->J_Halpha_Cprime);
    try_J("J_Halpha_Hbeta2", buf->J_Halpha_Hbeta2);
    try_J("J_Halpha_Hbeta3", buf->J_Halpha_Hbeta3);
    try_J("J_N_Cgamma", buf->J_N_Cgamma);
    try_exist("J_Cprime_Cgamma_exists", buf->J_Cprime_Cgamma_exists);
    try_exist("J_HN_Cbeta_exists", buf->J_HN_Cbeta_exists);
    try_exist("J_HN_Cprime_exists", buf->J_HN_Cprime_exists);
    try_exist("J_HN_Halpha_exists", buf->J_HN_Halpha_exists);
    try_exist("J_Halpha_Cprime_exists", buf->J_Halpha_Cprime_exists);
    try_exist("J_Halpha_Hbeta_exists", buf->J_Halpha_Hbeta_exists);
    try_exist("J_N_Cgamma_exists", buf->J_N_Cgamma_exists);
    try_exist("J_chi1_exists", buf->J_chi1_exists);
    if (grp.exist("residue_index_per_atom"))
        ReadFlat<int32_t>(grp.getDataSet("residue_index_per_atom"), buf->residue_index_per_atom, n_atoms);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttributeQ(grp, "units", buf->units);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

// ── Special-shape TR readers ──────────────────────────────────────

void ReadDihedralTimeSeries(HighFive::File& file,
                            const char* group_path,
                            std::size_t n_residues,
                            std::size_t n_atoms,
                            std::unique_ptr<QtDihedralTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("phi") || !grp.exist("psi") || !grp.exist("omega")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing phi/psi/omega"));
        return;
    }
    auto ds_phi = grp.getDataSet("phi");
    const auto dims_phi = ds_phi.getDimensions();
    if (dims_phi.size() != 2 || dims_phi[0] != n_residues) {
        WarnShapeMismatch(group_path, QStringLiteral("phi shape != [n_residues=%1, T]").arg(n_residues));
        return;
    }
    auto buf = std::make_unique<QtDihedralTimeSeries>();
    buf->n_residues = n_residues;
    buf->n_frames = dims_phi[1];
    ReadFlat<double>(ds_phi, buf->phi, n_residues * buf->n_frames);
    ReadFlat<double>(grp.getDataSet("psi"), buf->psi, n_residues * buf->n_frames);
    ReadFlat<double>(grp.getDataSet("omega"), buf->omega, n_residues * buf->n_frames);
    if (grp.exist("omega_deviation"))
        ReadFlat<double>(grp.getDataSet("omega_deviation"), buf->omega_deviation, n_residues * buf->n_frames);
    if (grp.exist("chi"))
        ReadFlat<double>(grp.getDataSet("chi"), buf->chi, n_residues * buf->n_frames * 4);
    if (grp.exist("rama_region"))
        ReadFlat<uint8_t>(grp.getDataSet("rama_region"), buf->rama_region, n_residues * buf->n_frames);
    if (grp.exist("chi_exists"))
        ReadFlat<uint8_t>(grp.getDataSet("chi_exists"), buf->chi_exists, n_residues * 4);
    if (grp.exist("is_proline"))
        ReadFlat<uint8_t>(grp.getDataSet("is_proline"), buf->is_proline, n_residues);
    if (grp.exist("is_glycine"))
        ReadFlat<uint8_t>(grp.getDataSet("is_glycine"), buf->is_glycine, n_residues);
    if (grp.exist("is_pre_proline"))
        ReadFlat<uint8_t>(grp.getDataSet("is_pre_proline"), buf->is_pre_proline, n_residues);
    if (grp.exist("omega_is_xpro"))
        ReadFlat<uint8_t>(grp.getDataSet("omega_is_xpro"), buf->omega_is_xpro, n_residues);
    if (grp.exist("residue_terminal_state"))
        ReadFlat<uint8_t>(grp.getDataSet("residue_terminal_state"), buf->residue_terminal_state, n_residues);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    TryReadAttributeQ(grp, "angle_units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttributeQ(grp, "angle_convention", buf->angle_convention);
    (void)n_atoms;  // chain_id_per_residue is intentionally not read
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadDssp8TimeSeries(HighFive::File& file,
                         const char* group_path,
                         std::size_t n_residues,
                         std::size_t n_atoms,
                         std::unique_ptr<QtDssp8TimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("ss8_code")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing ss8_code"));
        return;
    }
    auto ds_ss = grp.getDataSet("ss8_code");
    const auto dims_ss = ds_ss.getDimensions();
    if (dims_ss.size() != 2 || dims_ss[0] != n_residues) {
        WarnShapeMismatch(group_path, QStringLiteral("ss8_code shape != [n_residues=%1, T]").arg(n_residues));
        return;
    }
    auto buf = std::make_unique<QtDssp8TimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_residues = n_residues;
    buf->n_frames = dims_ss[1];
    ReadFlat<uint8_t>(ds_ss, buf->ss8_code, n_residues * buf->n_frames);
    if (grp.exist("hbond_acceptor_partner"))
        ReadFlat<int32_t>(grp.getDataSet("hbond_acceptor_partner"),
                          buf->hbond_acceptor_partner,
                          n_residues * buf->n_frames * 2);
    if (grp.exist("hbond_acceptor_energy"))
        ReadFlat<double>(grp.getDataSet("hbond_acceptor_energy"), buf->hbond_acceptor_energy, n_residues * buf->n_frames * 2);
    if (grp.exist("hbond_donor_partner"))
        ReadFlat<int32_t>(grp.getDataSet("hbond_donor_partner"), buf->hbond_donor_partner, n_residues * buf->n_frames * 2);
    if (grp.exist("hbond_donor_energy"))
        ReadFlat<double>(grp.getDataSet("hbond_donor_energy"), buf->hbond_donor_energy, n_residues * buf->n_frames * 2);
    if (grp.exist("residue_index_per_atom"))
        ReadFlat<int32_t>(grp.getDataSet("residue_index_per_atom"), buf->residue_index_per_atom, n_atoms);
    TryReadFrameMeta(grp, buf->meta, buf->n_frames, group_path);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttributeQ(grp, "ss8_legend", buf->ss8_legend);
    TryReadAttributeQ(grp, "hbond_threshold", buf->hbond_threshold);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadRingNeighbourhoodTimeSeries(HighFive::File& file,
                                     const char* group_path,
                                     std::size_t n_atoms,
                                     std::unique_ptr<QtRingNeighbourhoodTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("geometry") || !grp.exist("ring_membership_per_atom")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing geometry/ring_membership_per_atom"));
        return;
    }
    auto ds = grp.getDataSet("geometry");
    const auto dims = ds.getDimensions();
    if (dims.size() != 4 || dims[0] != n_atoms || dims[3] != 4) {
        WarnShapeMismatch(group_path, QStringLiteral("geometry shape != [n_atoms=%1, T, n_slots, 4]").arg(n_atoms));
        return;
    }
    auto buf = std::make_unique<QtRingNeighbourhoodTimeSeries>();
    buf->n_atoms = n_atoms;
    buf->n_frames = dims[1];
    buf->n_slots = dims[2];
    buf->n_channels = dims[3];
    ReadFlat<double>(ds, buf->geometry, n_atoms * buf->n_frames * buf->n_slots * buf->n_channels);
    ReadFlat<int32_t>(grp.getDataSet("ring_membership_per_atom"), buf->ring_membership_per_atom, n_atoms * buf->n_slots);
    if (grp.exist("frame_indices"))
        grp.getDataSet("frame_indices").read(buf->frame_indices);
    if (grp.exist("frame_times"))
        grp.getDataSet("frame_times").read(buf->frame_times);
    if (grp.exist("source_attached_per_frame"))
        grp.getDataSet("source_attached_per_frame").read(buf->source_attached);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "channel_layout", buf->channel_layout);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttribute(grp, "ring_current_spatial_cutoff_A", buf->ring_current_spatial_cutoff_A);
    TryReadAttributeQ(grp, "z_sign_convention", buf->z_sign_convention);
    TryReadAttributeQ(grp, "in_plane_angle_range", buf->in_plane_angle_range);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadBondLengthStats(HighFive::File& file, const char* group_path, std::unique_ptr<QtBondLengthStats>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("length_mean")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing length_mean"));
        return;
    }
    auto ds = grp.getDataSet("length_mean");
    std::vector<double> tmp;
    ds.read(tmp);
    const std::size_t n_bonds = tmp.size();
    auto buf = std::make_unique<QtBondLengthStats>();
    buf->n_bonds = n_bonds;
    buf->length_mean = std::move(tmp);
    if (grp.exist("length_std"))
        grp.getDataSet("length_std").read(buf->length_std);
    if (grp.exist("length_min"))
        grp.getDataSet("length_min").read(buf->length_min);
    if (grp.exist("length_max"))
        grp.getDataSet("length_max").read(buf->length_max);
    if (grp.exist("length_delta_mean"))
        grp.getDataSet("length_delta_mean").read(buf->length_delta_mean);
    if (grp.exist("length_delta_std"))
        grp.getDataSet("length_delta_std").read(buf->length_delta_std);
    if (grp.exist("atom_a"))
        grp.getDataSet("atom_a").read(buf->atom_a);
    if (grp.exist("atom_b"))
        grp.getDataSet("atom_b").read(buf->atom_b);
    if (grp.exist("bond_order"))
        grp.getDataSet("bond_order").read(buf->bond_order);
    if (grp.exist("bond_category"))
        grp.getDataSet("bond_category").read(buf->bond_category);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadGromacsEnergyTimeSeries(HighFive::File& file, const char* group_path, std::unique_ptr<QtSystemEnergyTimeSeries>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    auto try_1d = [&](const char* name, std::vector<double>& v) {
        if (grp.exist(name))
            grp.getDataSet(name).read(v);
    };
    // 2D (T, 9) tensor reads need the raw-pointer overload.
    auto try_2d_flat = [&](const char* name, std::vector<double>& v) {
        if (!grp.exist(name))
            return;
        auto ds = grp.getDataSet(name);
        const auto dims = ds.getDimensions();
        if (dims.size() != 2)
            return;
        v.resize(dims[0] * dims[1]);
        ds.read(v.data());
    };
    auto buf = std::make_unique<QtSystemEnergyTimeSeries>();
    try_1d("kinetic", buf->kinetic);
    try_1d("potential", buf->potential);
    try_1d("total_energy", buf->total_energy);
    try_1d("enthalpy", buf->enthalpy);
    try_1d("temperature", buf->temperature);
    try_1d("T_protein", buf->T_protein);
    try_1d("T_non_protein", buf->T_non_protein);
    try_1d("pressure", buf->pressure);
    try_1d("density", buf->density);
    try_1d("volume", buf->volume);
    try_1d("box_x", buf->box_x);
    try_1d("box_y", buf->box_y);
    try_1d("box_z", buf->box_z);
    try_1d("bond", buf->bond);
    try_1d("angle", buf->angle);
    try_1d("proper_dih", buf->proper_dih);
    try_1d("improper_dih", buf->improper_dih);
    try_1d("urey_bradley", buf->urey_bradley);
    try_1d("cmap_dih", buf->cmap_dih);
    try_1d("coulomb_sr", buf->coulomb_sr);
    try_1d("coulomb_recip", buf->coulomb_recip);
    try_1d("coulomb_14", buf->coulomb_14);
    try_1d("lj_sr", buf->lj_sr);
    try_1d("lj_14", buf->lj_14);
    try_1d("disper_corr", buf->disper_corr);
    try_2d_flat("pressure_tensor", buf->pressure_tensor);
    try_2d_flat("virial", buf->virial);
    try_1d("energy_frame_times_ps", buf->energy_frame_times_ps);
    if (grp.exist("frame_indices"))
        grp.getDataSet("frame_indices").read(buf->frame_indices);
    if (grp.exist("frame_times"))
        grp.getDataSet("frame_times").read(buf->frame_times);
    if (grp.exist("source_attached_per_frame"))
        grp.getDataSet("source_attached_per_frame").read(buf->source_attached);
    buf->n_frames = buf->frame_times.size();
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttributeQ(grp, "tensor_layout", buf->tensor_layout);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadRmsdTracking(HighFive::File& file, const char* group_path, std::unique_ptr<QtRmsdTracking>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    if (!grp.exist("rmsd")) {
        WarnShapeMismatch(group_path, QStringLiteral("missing rmsd"));
        return;
    }
    auto buf = std::make_unique<QtRmsdTracking>();
    grp.getDataSet("rmsd").read(buf->rmsd);
    buf->n_frames = buf->rmsd.size();
    if (grp.exist("atom_indices"))
        grp.getDataSet("atom_indices").read(buf->atom_indices);
    if (grp.exist("frame_indices"))
        grp.getDataSet("frame_indices").read(buf->frame_indices);
    if (grp.exist("frame_times"))
        grp.getDataSet("frame_times").read(buf->frame_times);
    if (grp.exist("source_attached_per_frame"))
        grp.getDataSet("source_attached_per_frame").read(buf->source_attached);
    TryReadAttributeQ(grp, "units", buf->units);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    TryReadAttributeQ(grp, "alignment_method", buf->alignment_method);
    TryReadAttributeQ(grp, "atom_selection", buf->atom_selection);
    TryReadAttributeQ(grp, "reference_frame_origin", buf->reference_frame_origin);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadDssp8Transitions(HighFive::File& file,
                          const char* group_path,
                          std::size_t n_residues,
                          std::unique_ptr<QtDssp8Transitions>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    auto buf = std::make_unique<QtDssp8Transitions>();
    buf->n_residues = n_residues;
    if (grp.exist("ss8_occupancy"))
        ReadFlat<uint32_t>(grp.getDataSet("ss8_occupancy"), buf->ss8_occupancy, n_residues * 8);
    if (grp.exist("ss8_transition_matrix"))
        ReadFlat<uint32_t>(grp.getDataSet("ss8_transition_matrix"), buf->ss8_transition_matrix, n_residues * 8 * 8);
    if (grp.exist("ss8_transition_count"))
        ReadFlat<uint32_t>(grp.getDataSet("ss8_transition_count"), buf->ss8_transition_count, n_residues);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

void ReadDihedralBinTransitions(HighFive::File& file,
                                const char* group_path,
                                std::size_t n_residues,
                                std::unique_ptr<QtDihedralBinTransitions>& out) {
    try {
    if (!file.exist(group_path)) {
        WarnGroupAbsent(group_path);
        return;
    }
    auto grp = file.getGroup(group_path);
    auto buf = std::make_unique<QtDihedralBinTransitions>();
    buf->n_residues = n_residues;
    if (grp.exist("backbone_bin_occupancy"))
        ReadFlat<uint32_t>(grp.getDataSet("backbone_bin_occupancy"), buf->backbone_bin_occupancy, n_residues * 6);
    if (grp.exist("backbone_transition_count"))
        ReadFlat<uint32_t>(grp.getDataSet("backbone_transition_count"), buf->backbone_transition_count, n_residues);
    if (grp.exist("backbone_dominant_region"))
        ReadFlat<uint8_t>(grp.getDataSet("backbone_dominant_region"), buf->backbone_dominant_region, n_residues);
    if (grp.exist("chi_rotamer_occupancy"))
        ReadFlat<uint32_t>(grp.getDataSet("chi_rotamer_occupancy"), buf->chi_rotamer_occupancy, n_residues * 4 * 3);
    if (grp.exist("chi_transition_count"))
        ReadFlat<uint32_t>(grp.getDataSet("chi_transition_count"), buf->chi_transition_count, n_residues * 4);
    if (grp.exist("chi_dominant_rotamer"))
        ReadFlat<uint8_t>(grp.getDataSet("chi_dominant_rotamer"), buf->chi_dominant_rotamer, n_residues * 4);
    TryReadAttributeQ(grp, "result_name", buf->result_name);
    out = std::move(buf);
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
        out.reset();
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
        out.reset();
    } catch (...) {
        LogReadException(group_path, "unknown");
        out.reset();
    }
}

// ── Selections reader ─────────────────────────────────────────────

void ParseSelectionMeta(QtSelectionEvent& ev) {
    if (ev.metadata_json_raw.isEmpty())
        return;
    QJsonParseError perr{};
    const QJsonDocument doc = QJsonDocument::fromJson(ev.metadata_json_raw.toUtf8(), &perr);
    if (perr.error != QJsonParseError::NoError || !doc.isObject())
        return;
    const QJsonObject obj = doc.object();
    switch (ev.kind) {
    case QtSelectionKind::DftPoseCoordinator: {
        QtDftPoseMeta m;
        if (obj.contains(QStringLiteral("score")))
            m.score = obj.value(QStringLiteral("score")).toDouble();
        if (obj.contains(QStringLiteral("method")))
            m.method = obj.value(QStringLiteral("method")).toString();
        ev.meta = m;
        break;
    }
    case QtSelectionKind::RmsdSpikeSelection: {
        QtRmsdSpikeMeta m;
        if (obj.contains(QStringLiteral("rmsd_angstroms")))
            m.rmsd_angstroms = obj.value(QStringLiteral("rmsd_angstroms")).toDouble();
        if (obj.contains(QStringLiteral("threshold")))
            m.threshold = obj.value(QStringLiteral("threshold")).toDouble();
        ev.meta = m;
        break;
    }
    case QtSelectionKind::ChiRotamerSelection: {
        QtChiRotamerMeta m;
        if (obj.contains(QStringLiteral("residue_index")))
            m.residue_index = static_cast<int32_t>(obj.value(QStringLiteral("residue_index")).toInt());
        if (obj.contains(QStringLiteral("chi_axis")))
            m.chi_axis = static_cast<int8_t>(obj.value(QStringLiteral("chi_axis")).toInt());
        if (obj.contains(QStringLiteral("from_rotamer")))
            m.from_rotamer = static_cast<int8_t>(obj.value(QStringLiteral("from_rotamer")).toInt());
        if (obj.contains(QStringLiteral("to_rotamer")))
            m.to_rotamer = static_cast<int8_t>(obj.value(QStringLiteral("to_rotamer")).toInt());
        ev.meta = m;
        break;
    }
    case QtSelectionKind::Unknown:
        break;
    }
}

void ReadSelections(HighFive::File& file, QtSelectionBag& bag) {
    constexpr const char* group_path = "/trajectory/selections";
    try {
    if (!file.exist("/trajectory/selections"))
        return;
    auto sel_grp = file.getGroup("/trajectory/selections");
    for (const auto& mangled : sel_grp.listObjectNames()) {
        const QString qmangled = QString::fromStdString(mangled);
        const QtSelectionKind kind = ParseSelectionGroupName(qmangled);
        if (kind == QtSelectionKind::Unknown) {
            h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                    QStringLiteral("QtTrajectoryH5"),
                                                    QStringLiteral("unknown selection kind mangled name: %1").arg(qmangled),
                                                    QString());
            continue;
        }
        auto kind_grp = sel_grp.getGroup(mangled);
        std::vector<uint64_t> frame_idx;
        std::vector<double> time_ps;
        std::vector<std::string> reason, meta_json;
        if (kind_grp.exist("frame_idx"))
            kind_grp.getDataSet("frame_idx").read(frame_idx);
        if (kind_grp.exist("time_ps"))
            kind_grp.getDataSet("time_ps").read(time_ps);
        if (kind_grp.exist("reason"))
            kind_grp.getDataSet("reason").read(reason);
        if (kind_grp.exist("metadata_json"))
            kind_grp.getDataSet("metadata_json").read(meta_json);
        const std::size_t n = std::min({frame_idx.size(), time_ps.size(), reason.size(), meta_json.size()});
        for (std::size_t i = 0; i < n; ++i) {
            QtSelectionEvent ev;
            ev.kind = kind;
            ev.frame_idx = frame_idx[i];
            ev.time_ps = time_ps[i];
            ev.reason = QString::fromStdString(reason[i]);
            ev.metadata_json_raw = QString::fromStdString(meta_json[i]);
            ParseSelectionMeta(ev);
            bag.push(std::move(ev));
        }
    }
    } catch (const HighFive::Exception& e) {
        LogReadException(group_path, "HighFive", e);
    } catch (const std::exception& e) {
        LogReadException(group_path, "std::exception", e);
    } catch (...) {
        LogReadException(group_path, "unknown");
    }
}

}  // namespace


// ─────────────────────────────────────────────────────────────────────
// Constructor — eager-load everything.
// ─────────────────────────────────────────────────────────────────────

QtTrajectoryH5::QtTrajectoryH5(const QString& h5_path) {
    using namespace HighFive;
    File const file(h5_path.toStdString(), File::ReadOnly);

    // ── Root attrs ────────────────────────────────────────────────
    uint64_t n_atoms_u64 = 0;
    if (!TryReadFileAttribute(file, "n_atoms", n_atoms_u64)) {
        throw std::runtime_error("QtTrajectoryH5: root attribute 'n_atoms' missing");
    }
    n_atoms_ = static_cast<std::size_t>(n_atoms_u64);
    TryReadFileAttributeQ(file, "protein_id", protein_id_);

    // ── /atoms identity (required) ───────────────────────────────
    if (!file.exist("/atoms"))
        throw std::runtime_error("QtTrajectoryH5: '/atoms' group missing");
    {
        auto atoms = file.getGroup("/atoms");
        atoms.getDataSet("element").read(atom_element_);
        atoms.getDataSet("residue_index").read(atom_residue_index_);
        std::vector<std::string> names_std;
        atoms.getDataSet("pdb_atom_name").read(names_std);
        atom_pdb_name_.reserve(names_std.size());
        for (const auto& s : names_std)
            atom_pdb_name_.push_back(QString::fromStdString(s));
        if (atom_element_.size() != n_atoms_ || atom_residue_index_.size() != n_atoms_ || atom_pdb_name_.size() != n_atoms_) {
            throw std::runtime_error("QtTrajectoryH5: /atoms dataset sizes inconsistent with n_atoms");
        }
    }

    // ── /trajectory/frames (required) ────────────────────────────
    if (!file.exist("/trajectory/frames"))
        throw std::runtime_error("QtTrajectoryH5: '/trajectory/frames' missing");
    {
        auto frames = file.getGroup("/trajectory/frames");
        frames.getDataSet("time_ps").read(frame_times_);
        if (frames.exist("original_index"))
            frames.getDataSet("original_index").read(frame_indices_);
        n_frames_ = frame_times_.size();
    }

    // ── /trajectory attrs ────────────────────────────────────────
    if (file.exist("/trajectory")) {
        auto traj = file.getGroup("/trajectory");
        TryReadAttributeQ(traj, "xtc_path", xtc_path_);
        TryReadAttributeQ(traj, "tpr_path", tpr_path_);
        TryReadAttributeQ(traj, "edr_path", edr_path_);
        TryReadAttributeQ(traj, "configuration", configuration_);
        // Group inventory for sparse-set introspection. Skip metadata
        // siblings (frames, selections, source — these aren't TRs).
        static const std::vector<std::string> kMetadataGroups = {
            "frames",
            "selections",
            "source",
        };
        for (const auto& g : traj.listObjectNames()) {
            const bool is_metadata = std::find(kMetadataGroups.begin(), kMetadataGroups.end(), g) != kMetadataGroups.end();
            if (!is_metadata)
                groups_present_.append(QString::fromStdString(g));
        }
    }

    // ── Positions (required for the animator) ────────────────────
    ReadPositionsTimeSeries(const_cast<File&>(file), "/trajectory/positions", n_atoms_, positions_);
    ValidateRequiredPositions(positions_.get(), n_frames_, n_atoms_);

    // ── Shielding TS family ──────────────────────────────────────
    ReadShieldingTimeSeries(const_cast<File&>(file), "/trajectory/bs_shielding_time_series", n_atoms_, bs_shielding_);
    ReadShieldingTimeSeries(const_cast<File&>(file), "/trajectory/hm_shielding_time_series", n_atoms_, hm_shielding_);
    ReadShieldingTimeSeries(const_cast<File&>(file), "/trajectory/mc_shielding_time_series", n_atoms_, mc_shielding_);
    ReadShieldingTimeSeries(const_cast<File&>(file), "/trajectory/piquad_shielding_time_series", n_atoms_, piquad_shielding_);
    ReadShieldingTimeSeries(const_cast<File&>(file), "/trajectory/ringchi_shielding_time_series", n_atoms_, ringchi_shielding_);
    ReadShieldingTimeSeries(const_cast<File&>(file), "/trajectory/disp_shielding_time_series", n_atoms_, disp_shielding_);
    ReadShieldingTimeSeries(const_cast<File&>(file), "/trajectory/hbond_shielding_time_series", n_atoms_, hbond_shielding_);
    // mopac_coulomb_shielding is T2-only (no T0 / T1); read as QtT2TimeSeries.
    ReadT2TimeSeries(const_cast<File&>(file),
                     "/trajectory/mopac_coulomb_shielding_time_series",
                     "t2",
                     n_atoms_,
                     mopac_coulomb_shielding_);
    ReadShieldingTimeSeries(const_cast<File&>(file),
                            "/trajectory/mopac_mc_shielding_time_series",
                            n_atoms_,
                            mopac_mc_shielding_);
    // mopac_vs_ff14sb is a single scalar per atom per frame (cos similarity
    // of MOPAC vs FF14SB EFG T2 vectors).
    ReadScalarTimeSeries(const_cast<File&>(file),
                         "/trajectory/mopac_vs_ff14sb_reconciliation",
                         "cos_t2",
                         n_atoms_,
                         mopac_vs_ff14sb_reconciliation_);
    ReadShieldingTimeSeries(const_cast<File&>(file),
                            "/trajectory/tripeptide_bb_shielding_time_series",
                            n_atoms_,
                            tripeptide_bb_shielding_);
    ReadShieldingTimeSeries(const_cast<File&>(file),
                            "/trajectory/tripeptide_neighbor_shielding_time_series",
                            n_atoms_,
                            tripeptide_neighbor_shielding_);
    ReadShieldingTimeSeries(const_cast<File&>(file),
                            "/trajectory/larsen_hbond_1pHB_shielding_time_series",
                            n_atoms_,
                            larsen_1pHB_);
    ReadShieldingTimeSeries(const_cast<File&>(file),
                            "/trajectory/larsen_hbond_1pHaB_shielding_time_series",
                            n_atoms_,
                            larsen_1pHaB_);
    ReadShieldingTimeSeries(const_cast<File&>(file),
                            "/trajectory/larsen_hbond_2pHB_shielding_time_series",
                            n_atoms_,
                            larsen_2pHB_);
    ReadShieldingTimeSeries(const_cast<File&>(file),
                            "/trajectory/larsen_hbond_2pHaB_shielding_time_series",
                            n_atoms_,
                            larsen_2pHaB_);
    // Water field is a composite (efg + efg_first + efield + efield_first + n_first + n_second).
    ReadWaterFieldTimeSeries(const_cast<File&>(file), "/trajectory/water_field_time_series", n_atoms_, water_field_ts_);

    // ── Scalar TS family ─────────────────────────────────────────
    ReadScalarTimeSeries(const_cast<File&>(file), "/trajectory/sasa_time_series", "sasa", n_atoms_, sasa_);
    ReadScalarTimeSeries(const_cast<File&>(file),
                         "/trajectory/aimnet2_charge_time_series",
                         "charge",
                         n_atoms_,
                         aimnet2_charge_);
    // Larsen H-bond count: the dataset is just "count", and the larsen
    // term is "water_term". Other Larsen TRs (1pHB/2pHB/etc.) are
    // shielding shape, handled above.
    ReadScalarTimeSeries(const_cast<File&>(file),
                         "/trajectory/larsen_hbond_count_time_series",
                         "count",
                         n_atoms_,
                         larsen_hbond_count_);
    ReadScalarTimeSeries(const_cast<File&>(file),
                         "/trajectory/larsen_hbond_water_term_time_series",
                         "water_term",
                         n_atoms_,
                         larsen_hbond_water_term_);
    ReadScalarTimeSeries(const_cast<File&>(file),
                         "/trajectory/bonded_energy_time_series",
                         "total",
                         n_atoms_,
                         bonded_energy_total_);
    // Hydration is composite (4 channels for shell, 5 for geometry).
    ReadHydrationShellTimeSeries(const_cast<File&>(file),
                                 "/trajectory/hydration_shell_time_series",
                                 n_atoms_,
                                 hydration_shell_ts_);
    ReadHydrationGeometryTimeSeries(const_cast<File&>(file),
                                    "/trajectory/hydration_geometry_time_series",
                                    n_atoms_,
                                    hydration_geometry_ts_);

    // ── Vec3 TS family ───────────────────────────────────────────
    ReadVec3TimeSeries(const_cast<File&>(file), "/trajectory/apbs_efield_time_series", "xyz", n_atoms_, apbs_efield_);
    ReadVec3TimeSeries(const_cast<File&>(file),
                       "/trajectory/tripeptide_bb_residual_vec_time_series",
                       "xyz",
                       n_atoms_,
                       tripeptide_bb_residual_vec_);
    ReadVec3TimeSeries(const_cast<File&>(file),
                       "/trajectory/tripeptide_neighbor_residual_vec_prev_time_series",
                       "xyz",
                       n_atoms_,
                       tripeptide_neighbor_residual_vec_prev_);
    ReadVec3TimeSeries(const_cast<File&>(file),
                       "/trajectory/tripeptide_neighbor_residual_vec_next_time_series",
                       "xyz",
                       n_atoms_,
                       tripeptide_neighbor_residual_vec_next_);

    // ── Special TS shapes ────────────────────────────────────────
    ReadT2TimeSeries(const_cast<File&>(file), "/trajectory/apbs_efg_time_series", "t2", n_atoms_, apbs_efg_);
    ReadTagTimeSeries(const_cast<File&>(file),
                      "/trajectory/tripeptide_bb_method_tag_time_series",
                      "method_tag",
                      n_atoms_,
                      tripeptide_bb_method_tag_);
    ReadEmbeddingTimeSeries(const_cast<File&>(file),
                            "/trajectory/aimnet2_embedding_time_series",
                            "embedding",
                            n_atoms_,
                            aimnet2_embedding_);
    ReadChargeResponseGradientTimeSeries(const_cast<File&>(file),
                                         "/trajectory/aimnet2_charge_response_gradient_time_series",
                                         n_atoms_,
                                         aimnet2_crg_);

    // ── Per-residue TRs ─────────────────────────────────────────
    // We need n_residues; pull from one of the per-residue groups'
    // dataset shapes. If neither dssp8 nor dihedral is present, the
    // per-residue readers no-op gracefully (n_residues=0 path).
    std::size_t n_residues = 0;
    if (file.exist("/trajectory/dihedral_time_series")) {
        auto grp = file.getGroup("/trajectory/dihedral_time_series");
        if (grp.exist("phi")) {
            const auto dims = grp.getDataSet("phi").getDimensions();
            if (!dims.empty())
                n_residues = dims[0];
        }
    } else if (file.exist("/trajectory/dssp8_time_series")) {
        auto grp = file.getGroup("/trajectory/dssp8_time_series");
        if (grp.exist("ss8_code")) {
            const auto dims = grp.getDataSet("ss8_code").getDimensions();
            if (!dims.empty())
                n_residues = dims[0];
        }
    }

    if (n_residues > 0) {
        ReadDihedralTimeSeries(const_cast<File&>(file), "/trajectory/dihedral_time_series", n_residues, n_atoms_, dihedrals_);
        ReadDssp8TimeSeries(const_cast<File&>(file), "/trajectory/dssp8_time_series", n_residues, n_atoms_, dssp8_);
        ReadJCouplingTimeSeries(const_cast<File&>(file),
                                "/trajectory/j_coupling_time_series",
                                n_residues,
                                n_atoms_,
                                j_coupling_full_);
        ReadDssp8Transitions(const_cast<File&>(file), "/trajectory/dssp8_transition", n_residues, dssp8_transitions_);
        ReadDihedralBinTransitions(const_cast<File&>(file),
                                   "/trajectory/dihedral_bin_transition",
                                   n_residues,
                                   dihedral_bin_transitions_);
    }
    // Ring pucker is per-ring (aromatic + saturated axes), not per-residue —
    // doesn't depend on n_residues discovery.
    ReadRingPuckerTimeSeries(const_cast<File&>(file), "/trajectory/ring_pucker_time_series", ring_pucker_);

    // ── Special-shape TRs ────────────────────────────────────────
    ReadRingNeighbourhoodTimeSeries(const_cast<File&>(file),
                                    "/trajectory/ring_neighbourhood_trajectory_stats",
                                    n_atoms_,
                                    ring_neighbourhood_);
    ReadBondLengthStats(const_cast<File&>(file), "/trajectory/bond_length_stats", bond_length_stats_);
    ReadGromacsEnergyTimeSeries(const_cast<File&>(file), "/trajectory/gromacs_energy_time_series", gromacs_energy_);
    ReadRmsdTracking(const_cast<File&>(file), "/trajectory/rmsd_tracking", rmsd_tracking_);

    // ── Welford rollups ──────────────────────────────────────────
    ReadShieldingWelford(const_cast<File&>(file), "/trajectory/bs_welford", n_atoms_, bs_welford_);
    ReadShieldingWelford(const_cast<File&>(file), "/trajectory/hm_welford", n_atoms_, hm_welford_);
    ReadShieldingWelford(const_cast<File&>(file), "/trajectory/mc_welford", n_atoms_, mc_welford_);
    ReadScalarWelford(const_cast<File&>(file), "/trajectory/sasa_welford", "sasa", n_atoms_, sasa_welford_);
    ReadScalarWelford(const_cast<File&>(file), "/trajectory/eeq_welford", "charge", n_atoms_, eeq_welford_);
    ReadScalarWelford(const_cast<File&>(file), "/trajectory/hbond_count_welford", "count", n_atoms_, hbond_count_welford_);
    ReadScalarWelford(const_cast<File&>(file), "/trajectory/mopac_charge_welford", "charge", n_atoms_, mopac_charge_welford_);
    // mopac_bond_order_welford uses channel prefix "order" (not "bond_order").
    ReadBondOrderWelford(const_cast<File&>(file), "/trajectory/mopac_bond_order_welford", mopac_bond_order_welford_);
    ReadVec3Welford(const_cast<File&>(file), "/trajectory/water_field_welford", n_atoms_, water_field_welford_);
    ReadVec3Welford(const_cast<File&>(file),
                    "/trajectory/aimnet2_charge_response_gradient_welford",
                    n_atoms_,
                    aimnet2_crg_welford_);
    ReadHydrationWelford(const_cast<File&>(file), "/trajectory/hydration_shell_welford", n_atoms_, hydration_shell_welford_);
    ReadHydrationWelford(const_cast<File&>(file),
                         "/trajectory/hydration_geometry_welford",
                         n_atoms_,
                         hydration_geometry_welford_);
    ReadAutocorrelation(const_cast<File&>(file), "/trajectory/bs_t0_autocorrelation", n_atoms_, bs_t0_autocorrelation_);

    // ── Bond-vector axis ─────────────────────────────────────────
    ReadIRedOrderParameters(const_cast<File&>(file), "/trajectory/ired_order_parameters", ired_order_parameters_);

    // ── Per-atom × per-channel composites ────────────────────────
    ReadKernelDynamics(const_cast<File&>(file), "/trajectory/kernel_dynamics", n_atoms_, kernel_dynamics_);
    ReadReorientationalDynamics(const_cast<File&>(file), "/trajectory/reorientational_dynamics", reorientational_dynamics_);
    ReadDihedralAutocorrelation(const_cast<File&>(file), "/trajectory/dihedral_autocorrelation", dihedral_autocorrelation_);
    ReadKernelCoherence(const_cast<File&>(file), "/trajectory/kernel_coherence", n_atoms_, kernel_coherence_);

    // ── Selections ───────────────────────────────────────────────
    ReadSelections(const_cast<File&>(file), selections_);
}

}  // namespace h5reader::io

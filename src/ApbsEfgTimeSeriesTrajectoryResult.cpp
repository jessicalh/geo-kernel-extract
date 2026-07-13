#include "ApbsEfgTimeSeriesTrajectoryResult.h"
#include "ApbsFieldResult.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <limits>
#include <typeinfo>

namespace nmr {


std::unique_ptr<ApbsEfgTimeSeriesTrajectoryResult>
ApbsEfgTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<ApbsEfgTimeSeriesTrajectoryResult>();
    r->per_atom_efg_.assign(tp.AtomCount(), std::vector<SphericalTensor>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Source-attached gate: in production ApbsFieldResult is
// RequireConformationResult'd in PerFrameExtractionSet, so
// HasResult<ApbsFieldResult>() is always true here. The gate is
// defensive — if a non-standard
// RunConfiguration omits the Require, this TR NaN-fills the absent
// frames and records mask=0 instead of capturing the zero-default
// apbs_efg_spherical (canonical "absent, not faked" — see
// feedback_capture_at_the_boundary).

void ApbsEfgTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const bool source_present = conf.HasResult<ApbsFieldResult>();
    if (!source_present) {
        OperationLog::Warn(
            "ApbsEfgTimeSeriesTrajectoryResult::Compute",
            "ApbsFieldResult not attached at frame " +
            std::to_string(frame_idx) +
            " — NaN-fill emitted + mask=0. Always-attached policy "
            "requires ApbsFieldResult; the run's PerFrameExtractionSet "
            "must RequireConformationResult(ApbsFieldResult).");
    }

    SphericalTensor nan_tensor{};
    if (!source_present) {
        const double nan_d = std::numeric_limits<double>::quiet_NaN();
        nan_tensor.T0 = nan_d;
        nan_tensor.T1 = {nan_d, nan_d, nan_d};
        nan_tensor.T2 = {nan_d, nan_d, nan_d, nan_d, nan_d};
    }

    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_efg_[i].push_back(
            source_present ? conf.AtomAt(i).apbs_efg_spherical : nan_tensor);
    }
    if (source_present) {
        const auto& apbs = conf.Result<ApbsFieldResult>();
        grid_dims_per_frame_.push_back(apbs.GridDims());
        grid_lengths_A_per_frame_.push_back(apbs.GridLengthsA());
        grid_origin_A_per_frame_.push_back(apbs.GridOriginA());
        grid_spacing_A_per_frame_.push_back(apbs.GridSpacingA());
        apbs_manual_grid_padding_A_ = apbs.ManualGridPaddingA();
        apbs_manual_grid_min_dim_A_ = apbs.ManualGridMinDimA();
        apbs_temperature_K_ = apbs.TemperatureK();
        apbs_thermal_voltage_V_ = apbs.ThermalVoltageV();
    } else {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        grid_dims_per_frame_.push_back({0, 0, 0});
        grid_lengths_A_per_frame_.push_back({nan, nan, nan});
        grid_origin_A_per_frame_.push_back({nan, nan, nan});
        grid_spacing_A_per_frame_.push_back({nan, nan, nan});
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_present ? 1u : 0u);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Transfer growing per-atom buffers into a contiguous atom-major
// DenseBuffer<SphericalTensor> owned by TrajectoryProtein.
//
// Idempotent after the first successful ownership transfer. History lengths
// are validated before any per-atom storage is released, so a malformed
// partial buffer is never adopted.

void ApbsEfgTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                  Trajectory& traj) {
    (void)traj;
    if (finalized_) return;
    const std::size_t N = tp.AtomCount();

    for (std::size_t i = 0; i < N; ++i) {
        if (per_atom_efg_[i].size() != n_frames_) {
            OperationLog::Error(
                "ApbsEfgTimeSeriesTrajectoryResult::Finalize",
                "per-atom APBS EFG history length mismatch at atom " +
                std::to_string(i));
            return;
        }
    }

    auto buffer = std::make_unique<DenseBuffer<SphericalTensor>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_efg_[i];
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<SphericalTensor>().swap(per_atom_efg_[i]);
        ++atoms_written;
    }
    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<SphericalTensor>(
        std::move(buffer),
        std::type_index(typeid(ApbsEfgTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "ApbsEfgTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) APBS EFG T2 SphericalTensor to tp dense buffer");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// Flat (N · T · 5) double array via explicit T2[k] component access.
// T2-only emission: APBS EFG is rank-2 symmetric-traceless after the
// symmetrization and trace projection applied at source, so T0 and T1
// are structurally zero.
//
// The 5-component trailing axis is the project-native real-tesseral
// (m=-2 ... m=+2) ordering matching SphericalTensor::T2 layout. The
// primary t2_* attrs pin that convention; legacy irrep attrs remain only
// as deprecated compatibility aliases.

void ApbsEfgTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(
            typeid(ApbsEfgTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "ApbsEfgTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();
    if (source_attached_per_frame_.size() != T ||
        grid_dims_per_frame_.size() != T ||
        grid_lengths_A_per_frame_.size() != T ||
        grid_origin_A_per_frame_.size() != T ||
        grid_spacing_A_per_frame_.size() != T) {
        OperationLog::Error(
            "ApbsEfgTimeSeriesTrajectoryResult::WriteH5Group",
            "APBS EFG history/metadata shape mismatch; group not emitted");
        return;
    }

    auto grp = file.createGroup("/trajectory/apbs_efg_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("t2_basis",
        std::string("project_native_t2_isometric_real_tesseral_v1"));
    grp.createAttribute("t2_component_order",
        std::string("T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("t2_frame",      std::string("cartesian_xyz_emitted_frame"));
    grp.createAttribute("t2_parity",     std::string("even"));
    grp.createAttribute("efg_t0_structural_zero", true);
    grp.createAttribute("efg_t1_structural_zero", true);
    grp.createAttribute("e3nn_export", std::string(
        "raw project tensor; call to_e3nn()/to_e3nn_T2() or "
        "project_t2_to_e3nn() before using e3nn Irreps"));
    grp.createAttribute("irrep_layout",
        std::string("T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("parity",        std::string("2e"));
    grp.createAttribute("legacy_irrep_attrs_deprecated", true);
    grp.createAttribute("units",         std::string("V/Å^2"));
    grp.createAttribute("directional_metadata_scope", std::string(
        "t2_basis,t2_component_order,t2_frame,t2_parity,e3nn_export,"
        "irrep_layout,normalization,parity,units describe t2 only; "
        "apbs_grid_*_per_frame datasets carry per-dataset lab-axis contracts"));
    grp.createAttribute("source_result", std::string("ApbsFieldResult"));
    grp.createAttribute("source_field",  std::string("apbs_efg_spherical"));
    grp.createAttribute("operation", std::string(
        "linearized_poisson_boltzmann_grid_hessian_traceless_t2"));
    grp.createAttribute("source", std::string(
        "ApbsFieldResult.apbs_efg_spherical; canonical APBS reaction "
        "EFG (total PB minus homogeneous-vacuum reference), sign-aligned "
        "potential Hessian, symmetrized and trace-projected at source; "
        "T2 components 0..4 only."));
    grp.createAttribute("source_attached_policy", std::string(
        "required_conformation_result -- production RunConfiguration "
        "requires ApbsFieldResult for this trajectory result; Compute "
        "defensively writes NaN-fill + mask=0 when the source is absent; "
        "source_attached_per_frame records that gate."));
    grp.createAttribute("field_quantity", std::string(
        "reaction_field_total_minus_homogeneous_vacuum_reference"));
    grp.createAttribute("reference_dielectric", 1.0);
    grp.createAttribute("reference_ionic_strength_M", 0.0);
    grp.createAttribute("reference_mobile_ion_count", 0);
    grp.createAttribute("reference_subtracts",
                        std::string("all_solute_charges"));
    grp.createAttribute("diagnostic_total_unclamped", true);
    grp.createAttribute("apbs_grid_mode", std::string("single_manual"));
    grp.createAttribute("apbs_manual_grid_padding_A",
                        apbs_manual_grid_padding_A_);
    grp.createAttribute("apbs_manual_grid_min_dim_A",
                        apbs_manual_grid_min_dim_A_);
    grp.createAttribute("apbs_temperature_K", apbs_temperature_K_);
    grp.createAttribute("apbs_thermal_voltage_V", apbs_thermal_voltage_V_);
    grp.createAttribute("max_potential_derivative_rank", 2);
    grp.createAttribute("higher_derivatives_present", false);
    grp.createAttribute("rank3_policy",
                        std::string("not_emitted_no_local_frame"));

    // (N, T, 5) via explicit T2[k] component access. No reinterpret,
    // no struct-packing assumption. Atom-major:
    // [atom_0_frame_0_T2(-2..+2), atom_0_frame_1_T2(-2..+2), ...,
    //  atom_1_frame_0_T2(-2..+2), ...].
    std::vector<double> flat(N * T * 5);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const SphericalTensor& st = buffer->At(i, t);
            const std::size_t base = (i * T + t) * 5;
            st.PackT2(&flat[base]);
        }
    }
    std::vector<std::size_t> dims = {N, T, std::size_t(5)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("t2", space);
    ds.write_raw(flat.data());
    ds.createAttribute("coordinate_frame",
        std::string("conformation_cartesian_xyz"));
    ds.createAttribute("transformation", std::string(
        "even_rank2: T'=R T R^T; emitted values are project-native T2 "
        "coefficients"));
    ds.createAttribute("parity", std::string("2e"));

    auto write_grid_u64 = [&](const std::string& name,
            const std::vector<std::array<std::uint64_t, 3>>& rows) {
        std::vector<std::uint64_t> values(T * 3);
        for (std::size_t t = 0; t < T; ++t)
            for (std::size_t d = 0; d < 3; ++d)
                values[t * 3 + d] = rows[t][d];
        HighFive::DataSpace grid_space({T, std::size_t(3)});
        auto grid_ds = grp.createDataSet<std::uint64_t>(name, grid_space);
        grid_ds.write_raw(values.data());
        grid_ds.createAttribute("component_order", std::string("x,y,z"));
        grid_ds.createAttribute("units", std::string("grid_points"));
        grid_ds.createAttribute("coordinate_frame",
            std::string("apbs_lab_axis_aligned_grid_xyz"));
        grid_ds.createAttribute("parity", std::string("even"));
        grid_ds.createAttribute("transformation", std::string(
            "configured cubic grid point counts (n,n,n); invariant under "
            "proper/improper orthogonal transforms and translations"));
    };
    auto write_grid_f64 = [&](const std::string& name,
            const std::vector<std::array<double, 3>>& rows,
            const std::string& transformation) {
        std::vector<double> values(T * 3);
        for (std::size_t t = 0; t < T; ++t)
            for (std::size_t d = 0; d < 3; ++d)
                values[t * 3 + d] = rows[t][d];
        HighFive::DataSpace grid_space({T, std::size_t(3)});
        auto grid_ds = grp.createDataSet<double>(name, grid_space);
        grid_ds.write_raw(values.data());
        grid_ds.createAttribute("component_order", std::string("x,y,z"));
        grid_ds.createAttribute("units", std::string("Angstrom"));
        grid_ds.createAttribute("coordinate_frame",
            std::string("apbs_lab_axis_aligned_grid_xyz"));
        grid_ds.createAttribute("parity", std::string("mixed"));
        grid_ds.createAttribute("transformation", transformation);
    };
    write_grid_u64("apbs_grid_dims_per_frame", grid_dims_per_frame_);
    write_grid_f64("apbs_grid_lengths_A_per_frame",
                   grid_lengths_A_per_frame_,
                   "lab-axis grid extents; translation invariant; no closed "
                   "O(3) transformation law under arbitrary orthogonal transforms");
    write_grid_f64("apbs_grid_origin_A_per_frame", grid_origin_A_per_frame_,
                   "lab-axis grid origin; o'=o+t under pure translations; no "
                   "affine-position/O(3) law under arbitrary orthogonal "
                   "transforms because grid axes remain lab-fixed");
    write_grid_f64("apbs_grid_spacing_A_per_frame",
                   grid_spacing_A_per_frame_,
                   "lab-axis grid step lengths; translation invariant; no "
                   "closed O(3) transformation law under arbitrary orthogonal "
                   "transforms");

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);
}

}  // namespace nmr

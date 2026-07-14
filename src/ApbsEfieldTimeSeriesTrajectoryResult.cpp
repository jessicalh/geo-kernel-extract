#include "ApbsEfieldTimeSeriesTrajectoryResult.h"
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


std::unique_ptr<ApbsEfieldTimeSeriesTrajectoryResult>
ApbsEfieldTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<ApbsEfieldTimeSeriesTrajectoryResult>();
    r->per_atom_efield_.assign(tp.AtomCount(), std::vector<Vec3>{});
    r->per_atom_clamp_mask_.assign(
        tp.AtomCount(), std::vector<std::uint8_t>{});
    r->per_atom_clamp_scale_.assign(
        tp.AtomCount(), std::vector<double>{});
    return r;
}


void ApbsEfieldTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const bool source_present = conf.HasResult<ApbsFieldResult>();
    if (!source_present) {
        OperationLog::Warn(
            "ApbsEfieldTimeSeriesTrajectoryResult::Compute",
            "ApbsFieldResult not attached at frame " +
            std::to_string(frame_idx) +
            " — NaN-fill emitted + source mask=0; clamp mask=0 and "
            "clamp scale=NaN");
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const Vec3 nan_vector(nan, nan, nan);
    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_efield_[i].push_back(
            source_present ? conf.AtomAt(i).apbs_efield : nan_vector);
        per_atom_clamp_mask_[i].push_back(
            source_present ? conf.AtomAt(i).apbs_efield_clamp_mask : 0u);
        per_atom_clamp_scale_[i].push_back(
            source_present ? conf.AtomAt(i).apbs_efield_clamp_scale : nan);
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
        efield_clamp_threshold_ = apbs.EfieldClampThreshold();
    } else {
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


void ApbsEfieldTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                    Trajectory& traj) {
    (void)traj;
    if (finalized_) return;
    const std::size_t N = tp.AtomCount();

    // Pre-validate every per-atom history length BEFORE any per-atom
    // storage is released, so a malformed partial buffer is never adopted
    // and no earlier atom's history is destroyed on a later mismatch.
    // Mirrors ApbsEfgTimeSeriesTrajectoryResult::Finalize's prevalidate
    // order (this loop previously validated and swap-released in one pass,
    // so a late mismatch returned after earlier histories were freed).
    for (std::size_t i = 0; i < N; ++i) {
        if (per_atom_efield_[i].size() != n_frames_ ||
            per_atom_clamp_mask_[i].size() != n_frames_ ||
            per_atom_clamp_scale_[i].size() != n_frames_) {
            OperationLog::Error(
                "ApbsEfieldTimeSeriesTrajectoryResult::Finalize",
                "per-atom APBS E/clamp history length mismatch at atom " +
                std::to_string(i));
            return;
        }
    }

    auto buffer = std::make_unique<DenseBuffer<Vec3>>(N, n_frames_);
    clamp_mask_flat_.assign(N * n_frames_, 0u);
    clamp_scale_flat_.assign(
        N * n_frames_, std::numeric_limits<double>::quiet_NaN());
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_efield_[i];
        const auto& mask_src = per_atom_clamp_mask_[i];
        const auto& scale_src = per_atom_clamp_scale_[i];
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
            const std::size_t offset = i * n_frames_ + f;
            clamp_mask_flat_[offset] = mask_src[f];
            clamp_scale_flat_[offset] = scale_src[f];
        }
        std::vector<Vec3>().swap(per_atom_efield_[i]);
        std::vector<std::uint8_t>().swap(per_atom_clamp_mask_[i]);
        std::vector<double>().swap(per_atom_clamp_scale_[i]);
        ++atoms_written;
    }
    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<Vec3>(
        std::move(buffer),
        std::type_index(typeid(ApbsEfieldTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "ApbsEfieldTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) + " frames) APBS E-field Vec3 to tp dense buffer");
}


void ApbsEfieldTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<Vec3>(std::type_index(
            typeid(ApbsEfieldTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "ApbsEfieldTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();
    if (clamp_mask_flat_.size() != N * T ||
        clamp_scale_flat_.size() != N * T ||
        source_attached_per_frame_.size() != T ||
        grid_dims_per_frame_.size() != T ||
        grid_lengths_A_per_frame_.size() != T ||
        grid_origin_A_per_frame_.size() != T ||
        grid_spacing_A_per_frame_.size() != T) {
        OperationLog::Error(
            "ApbsEfieldTimeSeriesTrajectoryResult::WriteH5Group",
            "APBS E-field history/metadata shape mismatch; group not emitted");
        return;
    }

    auto grp = file.createGroup("/trajectory/apbs_efield_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("irrep_layout",  std::string("x,y,z"));
    grp.createAttribute("normalization", std::string("cartesian"));
    grp.createAttribute("parity",        std::string("1o"));
    grp.createAttribute("units",         std::string("V/Angstrom"));
    grp.createAttribute("directional_metadata_scope", std::string(
        "irrep_layout,normalization,parity,units describe xyz only; "
        "clamp_* and apbs_grid_*_per_frame datasets carry per-dataset "
        "scalar/lab-axis contracts"));
    grp.createAttribute("source_result", std::string("ApbsFieldResult"));
    grp.createAttribute("source_field", std::string("apbs_efield"));
    grp.createAttribute("source", std::string(
        "ApbsFieldResult.apbs_efield; canonical APBS reaction field "
        "(total PB minus homogeneous-vacuum reference)"));
    grp.createAttribute("source_attached_policy", std::string(
        "required_conformation_result -- production requires "
        "ApbsFieldResult; absent frames are NaN-filled with source mask=0"));

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
    grp.createAttribute("efield_clamp_units", std::string("V/Angstrom"));
    grp.createAttribute("efield_clamp_threshold", efield_clamp_threshold_);
    grp.createAttribute("max_potential_derivative_rank", 2);
    grp.createAttribute("higher_derivatives_present", false);
    grp.createAttribute("rank3_policy",
                        std::string("not_emitted_no_local_frame"));

    // (N, T, 3) via explicit .x()/.y()/.z() access.
    std::vector<double> flat(N * T * 3);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const Vec3& v = buffer->At(i, t);
            const std::size_t base = (i * T + t) * 3;
            flat[base + 0] = v.x();
            flat[base + 1] = v.y();
            flat[base + 2] = v.z();
        }
    }
    std::vector<std::size_t> dims = {N, T, std::size_t(3)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());
    ds.createAttribute("coordinate_frame",
        std::string("conformation_cartesian_xyz"));
    ds.createAttribute("transformation",
        std::string(
            "continuum polar_vector: v'=R v; translation invariant. The "
            "live axis-aligned finite-difference APBS solve has no exact "
            "O(3) law; transformed production reruns use the recorded "
            "1.8e-2 V/Angstrom absolute + 5e-2 relative finite-grid envelope"));
    ds.createAttribute("parity", std::string("1o"));

    HighFive::DataSpace scalar_space({N, T});
    auto clamp_mask_ds =
        grp.createDataSet<std::uint8_t>("clamp_mask", scalar_space);
    clamp_mask_ds.write_raw(clamp_mask_flat_.data());
    clamp_mask_ds.createAttribute("units", std::string("dimensionless"));
    clamp_mask_ds.createAttribute("coordinate_frame", std::string(
        "lab_fixed_apbs_finite_difference_grid"));
    clamp_mask_ds.createAttribute("parity", std::string("mixed"));
    clamp_mask_ds.createAttribute("transformation", std::string(
        "continuum rotation/translation/reflection-invariant scalar "
        "threshold diagnostic derived from |E|; the live axis-aligned "
        "finite-difference APBS solve has no exact O(3) law"));
    clamp_mask_ds.createAttribute("validity", std::string(
        "uint8 0/1; 1 iff the canonical reaction E-field row was "
        "magnitude-clamped; absent-source frames use 0 and must be "
        "interpreted with source_attached_per_frame"));
    auto clamp_scale_ds =
        grp.createDataSet<double>("clamp_scale", scalar_space);
    clamp_scale_ds.write_raw(clamp_scale_flat_.data());
    clamp_scale_ds.createAttribute("units", std::string("dimensionless"));
    clamp_scale_ds.createAttribute("coordinate_frame", std::string(
        "lab_fixed_apbs_finite_difference_grid"));
    clamp_scale_ds.createAttribute("parity", std::string("mixed"));
    clamp_scale_ds.createAttribute("transformation", std::string(
        "continuum rotation/translation/reflection-invariant scalar "
        "derived from |E| and the configured clamp threshold; the live "
        "axis-aligned finite-difference APBS solve has no exact O(3) law"));
    clamp_scale_ds.createAttribute("validity", std::string(
        "finite in (0,1] when source is attached; 1 when unclamped; "
        "absent-source frames are NaN and must be interpreted with "
        "source_attached_per_frame"));

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
    grp.createDataSet("frame_times", frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame",
                      source_attached_per_frame_);
}

}  // namespace nmr

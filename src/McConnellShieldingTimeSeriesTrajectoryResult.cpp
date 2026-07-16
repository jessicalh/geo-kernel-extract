#include "McConnellShieldingTimeSeriesTrajectoryResult.h"
#include "McConnellResult.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <typeinfo>

namespace nmr {


std::unique_ptr<McConnellShieldingTimeSeriesTrajectoryResult>
McConnellShieldingTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<McConnellShieldingTimeSeriesTrajectoryResult>();
    r->per_atom_shielding_.assign(tp.AtomCount(),
                                  std::vector<SphericalTensor>{});
    return r;
}


void McConnellShieldingTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_shielding_[i].push_back(
            conf.AtomAt(i).mc_shielding_contribution);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void McConnellShieldingTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer =
        std::make_unique<DenseBuffer<SphericalTensor>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_shielding_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<SphericalTensor>().swap(per_atom_shielding_[i]);
        ++atoms_written;
    }

    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<SphericalTensor>(
        std::move(buffer),
        std::type_index(typeid(McConnellShieldingTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "McConnellShieldingTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) SphericalTensor time-series to tp dense buffer");
}


void McConnellShieldingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(
            typeid(McConnellShieldingTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "McConnellShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup("/trajectory/mc_shielding_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("tensor_basis",
        std::string(kMcConnellPackFull9TensorBasis));
    grp.createAttribute("tensor_component_order",
        std::string(kMcConnellPackFull9ComponentOrder));
    grp.createAttribute("tensor_frame",
        std::string(kMcConnellPackFull9TensorFrame));
    grp.createAttribute("coordinate_frame",
        std::string(kMcConnellPackFull9TensorFrame));
    grp.createAttribute("parity", std::string("0e+1e+2e"));
    grp.createAttribute("tensor_parity", std::string("even"));
    grp.createAttribute("tensor_transformation",
        std::string("even_rank2: T'=R T R^T"));
    grp.createAttribute("tensor_t1_semantics", std::string(
        "Cartesian Levi-Civita dual x,y,z (not real-Y1m): "
        "a=((T_yz-T_zy)/2,(T_zx-T_xz)/2,(T_xy-T_yx)/2); "
        "axial a'=det(R) R a; generically nonzero"));
    grp.createAttribute("tensor_t1_structural_zero", false);
    grp.createAttribute("tensor_structural_zero_components",
        std::string("none"));
    grp.createAttribute("e3nn_export",
        std::string(kMcConnellPackFull9E3nnExport));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("normalization_scope", std::string(
        "xyz tensor payload: T2 uses isometric real-tesseral "
        "normalization; T1 uses the tensor_t1_semantics convention"));
    grp.createAttribute("units",         std::string("Angstrom^-3"));

    std::vector<double> flat(N * T * 9);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const SphericalTensor& st = buffer->At(i, t);
            const std::size_t base = (i * T + t) * 9;
            st.PackFull9(&flat[base]);
        }
    }

    std::vector<std::size_t> dims = {N, T, std::size_t(9)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);
}

}  // namespace nmr

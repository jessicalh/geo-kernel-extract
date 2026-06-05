#include "ApbsFieldResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "ForceFieldChargeTable.h"
#include "GeometryResult.h"
#include "NpyWriter.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "SasaResult.h"
#include "SpatialIndexResult.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Args {
    std::string trajectory_dir;
    std::string h5_path;
    std::string params_path;
    std::string out_dir;
    std::string emit_rows_csv;
    std::vector<std::size_t> sample_rows{0, 250, 501};
};

struct H5Positions {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;
    std::vector<double> xyz;  // atom-major: (N, T, 3)
    std::vector<std::uint64_t> frame_indices;
    std::vector<double> frame_times;

    std::vector<nmr::Vec3> frame(std::size_t h5_row) const {
        if (h5_row >= n_frames) {
            throw std::runtime_error("h5_row out of range");
        }
        std::vector<nmr::Vec3> out(n_atoms);
        for (std::size_t atom = 0; atom < n_atoms; ++atom) {
            const std::size_t base = (atom * n_frames + h5_row) * 3;
            out[atom] = nmr::Vec3(xyz[base + 0], xyz[base + 1], xyz[base + 2]);
        }
        return out;
    }
};

struct SolveOutput {
    std::vector<double> e;      // (N,3)
    std::vector<double> efg_t2; // (N,5)
    std::vector<double> sasa;   // optional, size N for sample rows
};

struct Metrics {
    std::size_t h5_row = 0;
    std::uint64_t original_index = 0;
    double time_ps = 0.0;
    std::string stratum;
    std::size_t n = 0;
    double e_cos_median = std::numeric_limits<double>::quiet_NaN();
    double e_cos_p05 = std::numeric_limits<double>::quiet_NaN();
    double e_mag_ratio_median = std::numeric_limits<double>::quiet_NaN();
    double e_mag_ratio_p05 = std::numeric_limits<double>::quiet_NaN();
    double e_mag_ratio_p95 = std::numeric_limits<double>::quiet_NaN();
    double efg_t2_component_r = std::numeric_limits<double>::quiet_NaN();
    double efg_abs_t2_r = std::numeric_limits<double>::quiet_NaN();
    double efg_mag_ratio_median = std::numeric_limits<double>::quiet_NaN();
    double efg_mag_ratio_p05 = std::numeric_limits<double>::quiet_NaN();
    double efg_mag_ratio_p95 = std::numeric_limits<double>::quiet_NaN();
    double real_E_median = std::numeric_limits<double>::quiet_NaN();
    double placeholder_E_median = std::numeric_limits<double>::quiet_NaN();
    double real_T2_median = std::numeric_limits<double>::quiet_NaN();
    double placeholder_T2_median = std::numeric_limits<double>::quiet_NaN();
};

void usage(const char* argv0) {
    std::cerr
        << "Usage: " << argv0 << " --trajectory-dir DIR --h5 trajectory.h5 "
        << "--params ff14sb_params.dat --out DIR [--frames 0,250,501] "
        << "[--emit-rows-csv buckingham_efield_aggregated.csv]\n";
}

std::vector<std::size_t> parse_rows(const std::string& s) {
    std::vector<std::size_t> rows;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        rows.push_back(static_cast<std::size_t>(std::stoull(item)));
    }
    return rows;
}

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto need_value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string(name) + " requires a value");
            }
            return argv[++i];
        };
        if (key == "--trajectory-dir") args.trajectory_dir = need_value("--trajectory-dir");
        else if (key == "--h5") args.h5_path = need_value("--h5");
        else if (key == "--params") args.params_path = need_value("--params");
        else if (key == "--out") args.out_dir = need_value("--out");
        else if (key == "--frames") args.sample_rows = parse_rows(need_value("--frames"));
        else if (key == "--emit-rows-csv") args.emit_rows_csv = need_value("--emit-rows-csv");
        else if (key == "--help" || key == "-h") {
            usage(argv[0]);
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + key);
        }
    }
    if (args.trajectory_dir.empty() || args.h5_path.empty() ||
        args.params_path.empty() || args.out_dir.empty()) {
        throw std::runtime_error("missing required argument");
    }
    return args;
}

H5Positions read_positions(const std::string& path) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    auto ds = file.getDataSet("/trajectory/positions/xyz");
    const auto dims = ds.getDimensions();
    if (dims.size() != 3 || dims[2] != 3) {
        throw std::runtime_error("/trajectory/positions/xyz shape is not (N,T,3)");
    }
    H5Positions out;
    out.n_atoms = dims[0];
    out.n_frames = dims[1];
    out.xyz.resize(out.n_atoms * out.n_frames * 3);
    ds.read(out.xyz.data());
    file.getDataSet("/trajectory/positions/frame_indices").read(out.frame_indices);
    file.getDataSet("/trajectory/positions/frame_times").read(out.frame_times);
    if (out.frame_indices.size() != out.n_frames || out.frame_times.size() != out.n_frames) {
        throw std::runtime_error("position frame metadata length mismatch");
    }
    return out;
}

std::vector<std::size_t> read_h5_rows_csv(const std::string& path) {
    if (path.empty()) return {};
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open emit rows CSV: " + path);
    std::string header;
    if (!std::getline(in, header)) throw std::runtime_error("empty CSV: " + path);

    std::vector<std::string> cols;
    std::stringstream hs(header);
    std::string col;
    while (std::getline(hs, col, ',')) cols.push_back(col);
    auto it = std::find(cols.begin(), cols.end(), "h5_row");
    if (it == cols.end()) throw std::runtime_error("CSV has no h5_row column: " + path);
    const std::size_t h5_col = static_cast<std::size_t>(std::distance(cols.begin(), it));

    std::set<std::size_t> unique;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ls(line);
        std::string token;
        std::size_t idx = 0;
        while (std::getline(ls, token, ',')) {
            if (idx == h5_col) {
                unique.insert(static_cast<std::size_t>(std::stoull(token)));
                break;
            }
            ++idx;
        }
    }
    return std::vector<std::size_t>(unique.begin(), unique.end());
}

std::unique_ptr<nmr::TrajectoryProtein> build_seeded_tp(
        const std::string& trajectory_dir,
        const std::vector<nmr::Vec3>& first_positions,
        double first_time_ps) {
    auto tp = std::make_unique<nmr::TrajectoryProtein>();
    if (!tp->BuildFromTrajectory(trajectory_dir)) {
        throw std::runtime_error("BuildFromTrajectory failed: " + tp->Error());
    }
    tp->Seed(first_positions, first_time_ps);
    return tp;
}

std::vector<nmr::AtomChargeRadius> make_charge_rows(
        const nmr::TrajectoryProtein& tp,
        const std::string& params_path,
        bool real_static_radii) {
    const auto& conf = tp.CanonicalConformation();
    std::string err;
    auto tpr_rows = tp.Charges()->LoadCharges(tp.ProteinRef(), conf, err);
    if (tpr_rows.empty()) {
        throw std::runtime_error("TPR charge load failed: " + err);
    }

    if (!real_static_radii) {
        for (auto& row : tpr_rows) {
            row.pb_radius = nmr::kCompatibilityPlaceholderPbRadiusAngstrom;
            row.status = nmr::ChargeAssignmentStatus::PlaceholderPbRadius;
        }
        return tpr_rows;
    }

    nmr::ParamFileChargeSource params(params_path);
    auto param_rows = params.LoadCharges(tp.ProteinRef(), conf, err);
    if (param_rows.empty()) {
        throw std::runtime_error("ff14SB param radii load failed: " + err);
    }
    if (param_rows.size() != tpr_rows.size()) {
        throw std::runtime_error("TPR/param charge row count mismatch");
    }
    for (std::size_t i = 0; i < tpr_rows.size(); ++i) {
        tpr_rows[i].pb_radius = param_rows[i].pb_radius;
        tpr_rows[i].status = nmr::ChargeAssignmentStatus::Matched;
    }
    return tpr_rows;
}

void install_charge_table(
        nmr::TrajectoryProtein& tp,
        const std::vector<nmr::AtomChargeRadius>& rows,
        nmr::ChargeModelKind kind,
        const std::string& source_description) {
    std::string err;
    auto rows_copy = rows;
    auto table = nmr::ForceFieldChargeTable::FromValues(
        nmr::ForceField::Amber_ff14SB, kind, source_description,
        std::move(rows_copy), tp.AtomCount(), err);
    if (!table) {
        throw std::runtime_error("charge-table install failed: " + err);
    }
    auto& protein = const_cast<nmr::Protein&>(tp.ProteinRef());
    protein.SetForceFieldCharges(std::move(table));
}

bool attach_basic_geometry_for_sasa(nmr::ProteinConformation& conf) {
    auto geom = nmr::GeometryResult::Compute(conf);
    if (!geom || !conf.AttachResult(std::move(geom))) return false;
    auto spatial = nmr::SpatialIndexResult::Compute(conf);
    if (!spatial || !conf.AttachResult(std::move(spatial))) return false;
    auto sasa = nmr::SasaResult::Compute(conf);
    if (!sasa || !conf.AttachResult(std::move(sasa))) return false;
    return true;
}

SolveOutput solve_frame(const nmr::TrajectoryProtein& tp,
                        const nmr::ChargeSource& charge_source,
                        const std::vector<nmr::Vec3>& positions,
                        bool compute_sasa) {
    auto conf = tp.TickConformation(positions);
    if (compute_sasa && !attach_basic_geometry_for_sasa(*conf)) {
        throw std::runtime_error("failed to attach geometry/spatial/SASA");
    }
    auto charge = nmr::ChargeAssignmentResult::Compute(*conf, charge_source);
    if (!charge || !conf->AttachResult(std::move(charge))) {
        throw std::runtime_error("ChargeAssignmentResult failed");
    }
    auto apbs = nmr::ApbsFieldResult::Compute(*conf);
    if (!apbs || !conf->AttachResult(std::move(apbs))) {
        throw std::runtime_error("ApbsFieldResult failed");
    }

    const std::size_t n = conf->AtomCount();
    SolveOutput out;
    out.e.resize(n * 3);
    out.efg_t2.resize(n * 5);
    if (compute_sasa) out.sasa.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& atom = conf->AtomAt(i);
        out.e[i * 3 + 0] = atom.apbs_efield.x();
        out.e[i * 3 + 1] = atom.apbs_efield.y();
        out.e[i * 3 + 2] = atom.apbs_efield.z();
        atom.apbs_efg_spherical.PackT2(&out.efg_t2[i * 5]);
        if (compute_sasa) out.sasa[i] = atom.atom_sasa;
    }
    return out;
}

double percentile(std::vector<double> values, double p) {
    values.erase(std::remove_if(values.begin(), values.end(),
        [](double v) { return !std::isfinite(v); }), values.end());
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const double pos = p * static_cast<double>(values.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(pos));
    const auto hi = static_cast<std::size_t>(std::ceil(pos));
    if (lo == hi) return values[lo];
    const double w = pos - static_cast<double>(lo);
    return values[lo] * (1.0 - w) + values[hi] * w;
}

double pearson(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.size() < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double sx = 0.0, sy = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
        sx += a[i];
        sy += b[i];
        ++n;
    }
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    const double mx = sx / static_cast<double>(n);
    const double my = sy / static_cast<double>(n);
    double num = 0.0, xx = 0.0, yy = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
        const double dx = a[i] - mx;
        const double dy = b[i] - my;
        num += dx * dy;
        xx += dx * dx;
        yy += dy * dy;
    }
    if (xx <= 0.0 || yy <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return num / std::sqrt(xx * yy);
}

Metrics compare_outputs(const SolveOutput& placeholder,
                        const SolveOutput& real,
                        const std::vector<std::size_t>& atoms,
                        std::size_t h5_row,
                        std::uint64_t original_index,
                        double time_ps,
                        std::string stratum) {
    Metrics m;
    m.h5_row = h5_row;
    m.original_index = original_index;
    m.time_ps = time_ps;
    m.stratum = std::move(stratum);
    m.n = atoms.size();

    std::vector<double> e_cos;
    std::vector<double> e_ratio;
    std::vector<double> e_real_mag;
    std::vector<double> e_placeholder_mag;
    std::vector<double> t2_real_mag;
    std::vector<double> t2_placeholder_mag;
    std::vector<double> t2_ratio;
    std::vector<double> t2_components_real;
    std::vector<double> t2_components_placeholder;

    constexpr double eps = 1e-12;
    for (std::size_t atom : atoms) {
        const double* ep = &placeholder.e[atom * 3];
        const double* er = &real.e[atom * 3];
        double dot = 0.0, np = 0.0, nr = 0.0;
        for (int d = 0; d < 3; ++d) {
            dot += ep[d] * er[d];
            np += ep[d] * ep[d];
            nr += er[d] * er[d];
        }
        np = std::sqrt(np);
        nr = std::sqrt(nr);
        if (np > eps && nr > eps) {
            e_cos.push_back(dot / (np * nr));
            e_ratio.push_back(nr / np);
        }
        e_real_mag.push_back(nr);
        e_placeholder_mag.push_back(np);

        const double* tp = &placeholder.efg_t2[atom * 5];
        const double* tr = &real.efg_t2[atom * 5];
        double mp = 0.0, mr = 0.0;
        for (int k = 0; k < 5; ++k) {
            mp += tp[k] * tp[k];
            mr += tr[k] * tr[k];
            t2_components_placeholder.push_back(tp[k]);
            t2_components_real.push_back(tr[k]);
        }
        mp = std::sqrt(mp);
        mr = std::sqrt(mr);
        if (mp > eps) t2_ratio.push_back(mr / mp);
        t2_real_mag.push_back(mr);
        t2_placeholder_mag.push_back(mp);
    }

    m.e_cos_median = percentile(e_cos, 0.50);
    m.e_cos_p05 = percentile(e_cos, 0.05);
    m.e_mag_ratio_median = percentile(e_ratio, 0.50);
    m.e_mag_ratio_p05 = percentile(e_ratio, 0.05);
    m.e_mag_ratio_p95 = percentile(e_ratio, 0.95);
    m.efg_t2_component_r = pearson(t2_components_placeholder, t2_components_real);
    m.efg_abs_t2_r = pearson(t2_placeholder_mag, t2_real_mag);
    m.efg_mag_ratio_median = percentile(t2_ratio, 0.50);
    m.efg_mag_ratio_p05 = percentile(t2_ratio, 0.05);
    m.efg_mag_ratio_p95 = percentile(t2_ratio, 0.95);
    m.real_E_median = percentile(e_real_mag, 0.50);
    m.placeholder_E_median = percentile(e_placeholder_mag, 0.50);
    m.real_T2_median = percentile(t2_real_mag, 0.50);
    m.placeholder_T2_median = percentile(t2_placeholder_mag, 0.50);
    return m;
}

std::vector<std::size_t> all_atoms(std::size_t n) {
    std::vector<std::size_t> out(n);
    std::iota(out.begin(), out.end(), std::size_t{0});
    return out;
}

std::vector<std::size_t> atoms_by_sasa_quantile(
        const std::vector<double>& sasa,
        bool surface) {
    std::vector<double> vals = sasa;
    const double cutoff = percentile(vals, surface ? 0.75 : 0.25);
    std::vector<std::size_t> atoms;
    for (std::size_t i = 0; i < sasa.size(); ++i) {
        if (surface) {
            if (sasa[i] >= cutoff) atoms.push_back(i);
        } else {
            if (sasa[i] <= cutoff) atoms.push_back(i);
        }
    }
    return atoms;
}

void write_sample_outputs(const fs::path& out_dir,
                          std::size_t h5_row,
                          const std::string& label,
                          const SolveOutput& out,
                          std::size_t n_atoms) {
    const fs::path frame_dir =
        out_dir / label / ("frame_" + std::to_string(h5_row));
    fs::create_directories(frame_dir);
    nmr::NpyWriter::WriteFloat64((frame_dir / "apbs_E.npy").string(),
                                 out.e.data(), n_atoms, 3);
    nmr::NpyWriter::WriteFloat64((frame_dir / "apbs_efg.npy").string(),
                                 out.efg_t2.data(), n_atoms, 5);
    if (!out.sasa.empty()) {
        nmr::NpyWriter::WriteFloat64((frame_dir / "atom_sasa.npy").string(),
                                     out.sasa.data(), n_atoms);
    }
}

void write_metrics_csv(const fs::path& path, const std::vector<Metrics>& rows) {
    std::ofstream out(path);
    out << "h5_row,original_index,time_ps,stratum,n,"
        << "E_cos_median,E_cos_p05,E_mag_ratio_median,E_mag_ratio_p05,E_mag_ratio_p95,"
        << "EFG_T2_component_r,EFG_absT2_r,EFG_mag_ratio_median,EFG_mag_ratio_p05,EFG_mag_ratio_p95,"
        << "real_E_median,placeholder_E_median,real_T2_median,placeholder_T2_median\n";
    for (const auto& m : rows) {
        out << m.h5_row << ',' << m.original_index << ',' << m.time_ps << ','
            << m.stratum << ',' << m.n << ','
            << m.e_cos_median << ',' << m.e_cos_p05 << ','
            << m.e_mag_ratio_median << ',' << m.e_mag_ratio_p05 << ','
            << m.e_mag_ratio_p95 << ',' << m.efg_t2_component_r << ','
            << m.efg_abs_t2_r << ',' << m.efg_mag_ratio_median << ','
            << m.efg_mag_ratio_p05 << ',' << m.efg_mag_ratio_p95 << ','
            << m.real_E_median << ',' << m.placeholder_E_median << ','
            << m.real_T2_median << ',' << m.placeholder_T2_median << '\n';
    }
}

void write_radii_summary(const fs::path& path,
                         const std::vector<nmr::AtomChargeRadius>& real_rows,
                         const std::vector<nmr::AtomChargeRadius>& placeholder_rows) {
    std::map<double, std::size_t> counts;
    std::vector<double> abs_delta;
    double min_r = std::numeric_limits<double>::infinity();
    double max_r = -std::numeric_limits<double>::infinity();
    double sum = 0.0;
    for (std::size_t i = 0; i < real_rows.size(); ++i) {
        const double r = real_rows[i].pb_radius;
        counts[r]++;
        min_r = std::min(min_r, r);
        max_r = std::max(max_r, r);
        sum += r;
        abs_delta.push_back(std::abs(r - placeholder_rows[i].pb_radius));
    }
    std::ofstream out(path);
    out << "source,n,min_radius,max_radius,mean_radius,median_abs_delta_from_placeholder,counts\n";
    out << "ff14SB flat pb_radius column," << real_rows.size() << ','
        << min_r << ',' << max_r << ','
        << (sum / static_cast<double>(real_rows.size())) << ','
        << percentile(abs_delta, 0.50) << ',';
    bool first = true;
    for (const auto& [radius, n] : counts) {
        if (!first) out << ';';
        first = false;
        out << radius << ':' << n;
    }
    out << '\n';
}

void write_fixed_h5(const fs::path& path,
                    std::size_t n_atoms,
                    const std::vector<std::size_t>& h5_rows,
                    const H5Positions& positions,
                    const std::vector<double>& e,
                    const std::vector<double>& efg_t2,
                    const std::string& params_path) {
    HighFive::File file(path.string(), HighFive::File::Overwrite);
    file.createAttribute("provenance", std::string(
        "APBS static-radii workaround: same 1P9J H5 positions and TPR charges; "
        "PB radii replaced with ff14SB flat pb_radius values; trajectory TPR "
        "producer path had used the uniform 1.5 A placeholder."));
    file.createAttribute("radii_source", params_path);
    file.createAttribute("charges_source", std::string("production.tpr stored charges"));
    file.createAttribute("units_E", std::string("V/Angstrom"));
    file.createAttribute("units_EFG", std::string("V/Angstrom^2"));

    std::vector<std::uint64_t> h5_rows_u64(h5_rows.begin(), h5_rows.end());
    std::vector<std::uint64_t> original_indices;
    std::vector<double> times;
    original_indices.reserve(h5_rows.size());
    times.reserve(h5_rows.size());
    for (std::size_t row : h5_rows) {
        original_indices.push_back(positions.frame_indices[row]);
        times.push_back(positions.frame_times[row]);
    }
    file.createDataSet("h5_rows", h5_rows_u64);
    file.createDataSet("original_indices", original_indices);
    file.createDataSet("time_ps", times);

    HighFive::DataSpace e_space({n_atoms, h5_rows.size(), std::size_t{3}});
    auto e_ds = file.createDataSet<double>("apbs_E", e_space);
    e_ds.write_raw(e.data());
    e_ds.createAttribute("units", std::string("V/Angstrom"));
    e_ds.createAttribute("layout", std::string("(atom, emitted_row, xyz)"));

    HighFive::DataSpace t2_space({n_atoms, h5_rows.size(), std::size_t{5}});
    auto t2_ds = file.createDataSet<double>("apbs_efg", t2_space);
    t2_ds.write_raw(efg_t2.data());
    t2_ds.createAttribute("units", std::string("V/Angstrom^2"));
    t2_ds.createAttribute("layout", std::string("(atom, emitted_row, T2_m-2..T2_m+2)"));
    t2_ds.createAttribute("normalization", std::string("isometric_real_sph"));
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Args args = parse_args(argc, argv);
        fs::create_directories(args.out_dir);

        const H5Positions positions = read_positions(args.h5_path);
        auto placeholder_tp = build_seeded_tp(args.trajectory_dir,
                                              positions.frame(0),
                                              positions.frame_times[0]);
        auto real_tp = build_seeded_tp(args.trajectory_dir,
                                       positions.frame(0),
                                       positions.frame_times[0]);

        auto placeholder_rows = make_charge_rows(*placeholder_tp, args.params_path, false);
        auto real_rows = make_charge_rows(*real_tp, args.params_path, true);
        install_charge_table(
            *placeholder_tp, placeholder_rows, nmr::ChargeModelKind::GromacsTpr,
            "ff14SB:tpr_charges+placeholder_1p5A");
        install_charge_table(
            *real_tp, real_rows, nmr::ChargeModelKind::Ff14SBParamFile,
            "ff14SB:tpr_charges+flat_pb_radius:" + args.params_path);
        write_radii_summary(fs::path(args.out_dir) / "radii_summary.csv",
                            real_rows, placeholder_rows);

        nmr::PreloadedChargeSource placeholder_source(
            placeholder_rows, nmr::ForceField::Amber_ff14SB);
        nmr::PreloadedChargeSource real_source(
            real_rows, nmr::ForceField::Amber_ff14SB);

        std::vector<Metrics> metrics;
        for (std::size_t row : args.sample_rows) {
            if (row >= positions.n_frames) {
                std::cerr << "skip sample row out of range: " << row << "\n";
                continue;
            }
            const auto frame = positions.frame(row);
            auto placeholder = solve_frame(*placeholder_tp, placeholder_source, frame,
                                           true);
            auto real = solve_frame(*real_tp, real_source, frame, true);
            write_sample_outputs(args.out_dir, row, "placeholder_1p5A",
                                 placeholder, positions.n_atoms);
            write_sample_outputs(args.out_dir, row, "real_static_radii",
                                 real, positions.n_atoms);

            const auto all = all_atoms(positions.n_atoms);
            const auto buried = atoms_by_sasa_quantile(real.sasa, false);
            const auto surface = atoms_by_sasa_quantile(real.sasa, true);
            metrics.push_back(compare_outputs(
                placeholder, real, all, row, positions.frame_indices[row],
                positions.frame_times[row], "all"));
            metrics.push_back(compare_outputs(
                placeholder, real, buried, row, positions.frame_indices[row],
                positions.frame_times[row], "buried_sasa_q1"));
            metrics.push_back(compare_outputs(
                placeholder, real, surface, row, positions.frame_indices[row],
                positions.frame_times[row], "surface_sasa_q4"));
            std::cerr << "A/B solved h5_row=" << row << "\n";
        }
        write_metrics_csv(fs::path(args.out_dir) / "ab_metrics.csv", metrics);

        const auto emit_rows = read_h5_rows_csv(args.emit_rows_csv);
        if (!emit_rows.empty()) {
            std::vector<double> e_fixed(positions.n_atoms * emit_rows.size() * 3);
            std::vector<double> efg_fixed(positions.n_atoms * emit_rows.size() * 5);

            for (std::size_t t = 0; t < emit_rows.size(); ++t) {
                const std::size_t row = emit_rows[t];
                if (row >= positions.n_frames) {
                    throw std::runtime_error("emit h5_row out of range: " +
                                             std::to_string(row));
                }
                auto solved = solve_frame(*real_tp, real_source,
                                          positions.frame(row), false);
                for (std::size_t atom = 0; atom < positions.n_atoms; ++atom) {
                    for (int d = 0; d < 3; ++d) {
                        e_fixed[(atom * emit_rows.size() + t) * 3 + d] =
                            solved.e[atom * 3 + d];
                    }
                    for (int k = 0; k < 5; ++k) {
                        efg_fixed[(atom * emit_rows.size() + t) * 5 + k] =
                            solved.efg_t2[atom * 5 + k];
                    }
                }
                std::cerr << "fixed real-radii solved " << (t + 1) << "/"
                          << emit_rows.size() << " h5_row=" << row << "\n";
            }
            write_fixed_h5(fs::path(args.out_dir) / "real_static_radii_apbs.h5",
                           positions.n_atoms, emit_rows, positions,
                           e_fixed, efg_fixed, args.params_path);
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        usage(argv[0]);
        return 1;
    }
}

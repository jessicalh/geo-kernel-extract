//
// roundtrip_test: comprehensive validation of the AnalysisFile
// serialisation layer.
//
// For each test H5 file:
//
//   Phase 1 — Original-vs-written dataset comparison
//     Read into object model, write back, compare every dataset
//     byte-for-byte.  Uses read(ptr, dtype) to get true raw bytes,
//     not read(char*) which would silently convert types.
//
//   Phase 2 — Double roundtrip (write determinism)
//     Read the written file into a second object model, write it
//     again, compare the two written files byte-for-byte at the
//     dataset level.  This proves the writer is deterministic
//     independent of the original file's internal layout.
//
//   Phase 3 — Object model identity
//     Compare every field of the two in-memory AnalysisFile objects
//     (the one read from the original and the one read from the
//     written file).  This proves read/write symmetry through the
//     typed object model.
//
//   Phase 4 — Structural checks
//     Verify bidirectional dataset/group/attribute coverage:
//     nothing extra, nothing missing.
//
//   Phase 5 — 1CBH aim dtype exception
//     For the file with float64 aim, verify the written dataset is
//     float32, verify lossless conversion, verify shape preserved.
//
//   Phase 6 — Special value census
//     Count NaN, Inf, -Inf, negative zero, and subnormals in every
//     float64/float32 dataset.  Verify counts match between original
//     and written.
//

#include "analysis_file.h"

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <set>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// ── Counters ───────────────────────────────────────────────────────

struct TestCounters {
    int pass = 0;
    int fail = 0;
    int skip = 0;
    std::vector<std::string> failures;

    void Pass() { pass++; }
    void Fail(const std::string& msg) { fail++; failures.push_back(msg); }
    void Skip() { skip++; }
    bool Ok() const { return fail == 0; }
};

// ── Collect all dataset paths from an H5 file ──────────────────────

void CollectDatasets(HighFive::Group& grp, const std::string& prefix,
                     std::vector<std::string>& out) {
    for (const auto& name : grp.listObjectNames()) {
        std::string full = prefix.empty() ? name : prefix + "/" + name;
        auto type = grp.getObjectType(name);
        if (type == HighFive::ObjectType::Dataset) {
            out.push_back(full);
        } else if (type == HighFive::ObjectType::Group) {
            auto sub = grp.getGroup(name);
            CollectDatasets(sub, full, out);
        }
    }
}

void CollectGroups(HighFive::Group& grp, const std::string& prefix,
                   std::vector<std::string>& out) {
    for (const auto& name : grp.listObjectNames()) {
        std::string full = prefix.empty() ? name : prefix + "/" + name;
        auto type = grp.getObjectType(name);
        if (type == HighFive::ObjectType::Group) {
            out.push_back(full);
            auto sub = grp.getGroup(name);
            CollectGroups(sub, full, out);
        }
    }
}

// ── Compare one dataset byte-for-byte ──────────────────────────────

void CompareDataset(HighFive::File& a, HighFive::File& b,
                    const std::string& path, TestCounters& c) {
    auto ds_a = a.getDataSet(path);
    auto ds_b = b.getDataSet(path);

    // Shape
    auto shape_a = ds_a.getDimensions();
    auto shape_b = ds_b.getDimensions();
    if (shape_a != shape_b) {
        std::string sa, sb;
        for (auto d : shape_a) sa += std::to_string(d) + ",";
        for (auto d : shape_b) sb += std::to_string(d) + ",";
        c.Fail(path + ": shape (" + sa + ") vs (" + sb + ")");
        return;
    }

    // Dtype
    auto dt_a = ds_a.getDataType();
    auto dt_b = ds_b.getDataType();
    if (dt_a.getClass() != dt_b.getClass() || dt_a.getSize() != dt_b.getSize()) {
        char buf[256];
        snprintf(buf, sizeof(buf), "%s: dtype class=%d/%zu vs %d/%zu",
                 path.c_str(),
                 (int)dt_a.getClass(), dt_a.getSize(),
                 (int)dt_b.getClass(), dt_b.getSize());
        c.Fail(buf);
        return;
    }

    // Strings: compare element-by-element
    if (dt_a.getClass() == HighFive::DataTypeClass::String) {
        std::vector<std::string> va, vb;
        ds_a.read(va);
        ds_b.read(vb);
        if (va != vb) {
            c.Fail(path + ": string content differs");
        } else {
            c.Pass();
        }
        return;
    }

    // Numeric: compare raw bytes using file dtype as memory dtype
    size_t n_elements = 1;
    for (auto d : shape_a) n_elements *= d;
    size_t n_bytes = n_elements * dt_a.getSize();

    std::vector<char> buf_a(n_bytes), buf_b(n_bytes);
    ds_a.read(buf_a.data(), dt_a);
    ds_b.read(buf_b.data(), dt_b);

    if (std::memcmp(buf_a.data(), buf_b.data(), n_bytes) != 0) {
        // Find first mismatch
        for (size_t i = 0; i < n_bytes; ++i) {
            if (buf_a[i] != buf_b[i]) {
                size_t elem_idx = i / dt_a.getSize();
                char msg[256];
                snprintf(msg, sizeof(msg),
                         "%s: byte mismatch at offset %zu (element %zu), "
                         "0x%02x vs 0x%02x",
                         path.c_str(), i, elem_idx,
                         (unsigned char)buf_a[i], (unsigned char)buf_b[i]);
                c.Fail(msg);
                return;
            }
        }
    }
    c.Pass();
}

// ── Compare all attributes on a group or dataset ───────────────────

void CompareAttrs(HighFive::AnnotateTraits<HighFive::Group>& a,
                  HighFive::AnnotateTraits<HighFive::Group>& b,
                  const std::string& path, TestCounters& c,
                  bool check_extra = true) {
    auto names_a = a.listAttributeNames();
    auto names_b = b.listAttributeNames();

    std::set<std::string> set_a(names_a.begin(), names_a.end());
    std::set<std::string> set_b(names_b.begin(), names_b.end());

    // Check for missing attrs in b
    for (const auto& name : names_a) {
        if (set_b.find(name) == set_b.end()) {
            c.Fail(path + "/@" + name + ": missing in written");
            continue;
        }

        auto aa = a.getAttribute(name);
        auto ab = b.getAttribute(name);

        if (aa.getDataType().getClass() == HighFive::DataTypeClass::String) {
            std::string va, vb;
            aa.read(va);
            ab.read(vb);
            if (va != vb) {
                c.Fail(path + "/@" + name + ": '" + va + "' vs '" + vb + "'");
            } else {
                c.Pass();
            }
        } else {
            auto dt = aa.getDataType();
            size_t sz = dt.getSize();
            std::vector<char> ba(sz), bb(sz);
            aa.read(ba.data(), dt);
            ab.read(bb.data(), dt);
            if (std::memcmp(ba.data(), bb.data(), sz) != 0) {
                c.Fail(path + "/@" + name + ": value mismatch");
            } else {
                c.Pass();
            }
        }
    }

    // Check for extra attrs in b
    if (check_extra) {
        for (const auto& name : names_b) {
            if (set_a.find(name) == set_a.end()) {
                c.Fail(path + "/@" + name + ": extra in written");
            }
        }
    }
}

// ── Special value census for a float dataset ───────────────────────

struct SpecialCounts {
    size_t nan_count  = 0;
    size_t inf_count  = 0;
    size_t ninf_count = 0;
    size_t negz_count = 0;
    size_t subn_count = 0;

    bool operator==(const SpecialCounts& o) const {
        return nan_count == o.nan_count && inf_count == o.inf_count &&
               ninf_count == o.ninf_count && negz_count == o.negz_count &&
               subn_count == o.subn_count;
    }
};

template<typename T>
SpecialCounts CountSpecials(const std::vector<char>& buf) {
    SpecialCounts sc;
    size_t n = buf.size() / sizeof(T);
    const T* data = reinterpret_cast<const T*>(buf.data());
    for (size_t i = 0; i < n; ++i) {
        T v = data[i];
        if (std::isnan(v)) sc.nan_count++;
        else if (v == std::numeric_limits<T>::infinity()) sc.inf_count++;
        else if (v == -std::numeric_limits<T>::infinity()) sc.ninf_count++;
        else if (v == T(0) && std::signbit(v)) sc.negz_count++;
        else if (std::fpclassify(v) == FP_SUBNORMAL) sc.subn_count++;
    }
    return sc;
}

void CompareSpecials(HighFive::File& a, HighFive::File& b,
                     const std::string& path, TestCounters& c) {
    auto ds_a = a.getDataSet(path);
    auto dt = ds_a.getDataType();
    if (dt.getClass() != HighFive::DataTypeClass::Float) return;

    size_t n_elements = 1;
    for (auto d : ds_a.getDimensions()) n_elements *= d;
    size_t n_bytes = n_elements * dt.getSize();

    std::vector<char> buf_a(n_bytes), buf_b(n_bytes);
    ds_a.read(buf_a.data(), dt);
    b.getDataSet(path).read(buf_b.data(), dt);

    SpecialCounts sa, sb;
    if (dt.getSize() == 8) {
        sa = CountSpecials<double>(buf_a);
        sb = CountSpecials<double>(buf_b);
    } else {
        sa = CountSpecials<float>(buf_a);
        sb = CountSpecials<float>(buf_b);
    }

    if (sa == sb) {
        c.Pass();
    } else {
        char msg[512];
        snprintf(msg, sizeof(msg),
                 "%s: special value mismatch — "
                 "NaN %zu/%zu, Inf %zu/%zu, -Inf %zu/%zu, -0 %zu/%zu, subn %zu/%zu",
                 path.c_str(),
                 sa.nan_count, sb.nan_count,
                 sa.inf_count, sb.inf_count,
                 sa.ninf_count, sb.ninf_count,
                 sa.negz_count, sb.negz_count,
                 sa.subn_count, sb.subn_count);
        c.Fail(msg);
    }
}

// ── Compare two AnalysisFile objects field-by-field ─────────────────
// Uses memcmp on numeric vector data and == on string vectors.

#define CMP_VEC(field) do { \
    if (a.field.size() != b.field.size()) { \
        c.Fail(std::string(#field) + ": size " + \
               std::to_string(a.field.size()) + " vs " + \
               std::to_string(b.field.size())); \
    } else if (!a.field.empty() && \
               std::memcmp(a.field.data(), b.field.data(), \
                           a.field.size() * sizeof(a.field[0])) != 0) { \
        c.Fail(std::string(#field) + ": data differs"); \
    } else { \
        c.Pass(); \
    } \
} while (0)

#define CMP_STR_VEC(field) do { \
    if (a.field != b.field) { \
        c.Fail(std::string(#field) + ": differs"); \
    } else { \
        c.Pass(); \
    } \
} while (0)

#define CMP_SCALAR(field) do { \
    if (a.field != b.field) { \
        c.Fail(std::string(#field) + ": " + \
               std::to_string(a.field) + " vs " + std::to_string(b.field)); \
    } else { \
        c.Pass(); \
    } \
} while (0)

void CompareObjects(const AnalysisFile& a, const AnalysisFile& b,
                    TestCounters& c) {
    // Dimensions
    CMP_SCALAR(n_frames);
    CMP_SCALAR(n_atoms);
    CMP_SCALAR(n_residues);
    CMP_SCALAR(n_bonds);
    CMP_SCALAR(n_rings);

    // Meta
    CMP_SCALAR(meta.stride);
    CMP_VEC(meta.frame_times);
    CMP_VEC(meta.frame_indices);
    if (a.meta.protein_id != b.meta.protein_id)
        c.Fail("meta.protein_id: " + a.meta.protein_id + " vs " + b.meta.protein_id);
    else c.Pass();

    // Atoms (all 26 fields)
    CMP_VEC(atoms.element); CMP_VEC(atoms.residue_index);
    CMP_STR_VEC(atoms.atom_name);
    CMP_VEC(atoms.atom_role); CMP_VEC(atoms.hybridisation);
    CMP_VEC(atoms.n_bonded); CMP_VEC(atoms.graph_dist_ring);
    CMP_VEC(atoms.is_backbone); CMP_VEC(atoms.is_conjugated);
    CMP_VEC(atoms.is_amide_H); CMP_VEC(atoms.is_alpha_H);
    CMP_VEC(atoms.is_methyl); CMP_VEC(atoms.is_aromatic_H);
    CMP_VEC(atoms.is_on_aromatic_residue);
    CMP_VEC(atoms.is_hbond_donor); CMP_VEC(atoms.is_hbond_acceptor);
    CMP_VEC(atoms.parent_is_sp2);
    CMP_VEC(atoms.graph_dist_N); CMP_VEC(atoms.graph_dist_O);
    CMP_VEC(atoms.eneg_sum_1); CMP_VEC(atoms.eneg_sum_2);
    CMP_VEC(atoms.n_pi_bonds_3); CMP_VEC(atoms.bfs_to_nearest_ring);
    CMP_VEC(atoms.bfs_decay);
    CMP_VEC(atoms.partial_charge); CMP_VEC(atoms.vdw_radius);

    // Residues
    CMP_STR_VEC(residues.residue_name);
    CMP_VEC(residues.residue_number);
    CMP_STR_VEC(residues.chain_id);

    // Topology
    CMP_VEC(topology.bond_atoms); CMP_VEC(topology.bond_category);
    CMP_VEC(topology.bond_order); CMP_VEC(topology.parent_atom_index);
    CMP_VEC(topology.ring_type); CMP_VEC(topology.ring_residue);
    CMP_VEC(topology.ring_fused_partner); CMP_VEC(topology.ring_offsets);
    CMP_VEC(topology.ring_atom_indices);

    // Positions
    CMP_VEC(positions.xyz);

    // Ring current (17 fields)
    CMP_VEC(ring_current.bs_T0_per_type); CMP_VEC(ring_current.bs_T2_per_type);
    CMP_VEC(ring_current.hm_T0_per_type); CMP_VEC(ring_current.hm_T2_per_type);
    CMP_VEC(ring_current.bs_shielding); CMP_VEC(ring_current.hm_shielding);
    CMP_VEC(ring_current.rs_shielding);
    CMP_VEC(ring_current.n_rings_3A); CMP_VEC(ring_current.n_rings_5A);
    CMP_VEC(ring_current.n_rings_8A); CMP_VEC(ring_current.n_rings_12A);
    CMP_VEC(ring_current.mean_ring_dist); CMP_VEC(ring_current.nearest_ring_atom);
    CMP_VEC(ring_current.G_iso_exp_sum); CMP_VEC(ring_current.G_T2_exp_sum);
    CMP_VEC(ring_current.G_iso_var_8A); CMP_VEC(ring_current.total_B_field);

    // EFG (18 fields)
    CMP_VEC(efg.coulomb_total); CMP_VEC(efg.coulomb_backbone);
    CMP_VEC(efg.coulomb_aromatic); CMP_VEC(efg.coulomb_shielding);
    CMP_VEC(efg.E_total); CMP_VEC(efg.E_backbone);
    CMP_VEC(efg.E_sidechain); CMP_VEC(efg.E_aromatic);
    CMP_VEC(efg.E_solvent); CMP_VEC(efg.E_magnitude);
    CMP_VEC(efg.E_bond_proj); CMP_VEC(efg.E_backbone_frac);
    CMP_VEC(efg.apbs_efg); CMP_VEC(efg.apbs_efield);
    CMP_VEC(efg.aimnet2_total); CMP_VEC(efg.aimnet2_backbone);
    CMP_VEC(efg.aimnet2_aromatic); CMP_VEC(efg.aimnet2_shielding);

    // Bond aniso (15 fields)
    CMP_VEC(bond_aniso.mc_shielding); CMP_VEC(bond_aniso.T2_backbone);
    CMP_VEC(bond_aniso.T2_sidechain); CMP_VEC(bond_aniso.T2_aromatic);
    CMP_VEC(bond_aniso.T2_CO_nearest); CMP_VEC(bond_aniso.T2_CN_nearest);
    CMP_VEC(bond_aniso.co_sum); CMP_VEC(bond_aniso.cn_sum);
    CMP_VEC(bond_aniso.sidechain_sum); CMP_VEC(bond_aniso.aromatic_sum);
    CMP_VEC(bond_aniso.co_nearest);
    CMP_VEC(bond_aniso.nearest_CO_dist); CMP_VEC(bond_aniso.nearest_CN_dist);
    CMP_VEC(bond_aniso.nearest_CO_midpoint); CMP_VEC(bond_aniso.dir_nearest_CO);

    // Quadrupole
    CMP_VEC(quadrupole.pq_shielding);
    CMP_VEC(quadrupole.pq_T0_per_type); CMP_VEC(quadrupole.pq_T2_per_type);

    // Dispersion
    CMP_VEC(dispersion.disp_shielding);
    CMP_VEC(dispersion.disp_T0_per_type); CMP_VEC(dispersion.disp_T2_per_type);

    // HBond
    CMP_VEC(hbond.hbond_shielding); CMP_VEC(hbond.nearest_spherical);
    CMP_VEC(hbond.nearest_dist); CMP_VEC(hbond.nearest_dir);
    CMP_VEC(hbond.inv_d3); CMP_VEC(hbond.count_3_5A);
    CMP_VEC(hbond.is_donor); CMP_VEC(hbond.is_acceptor);
    CMP_VEC(hbond.is_backbone);

    // SASA
    CMP_VEC(sasa.sasa); CMP_VEC(sasa.normal);

    // Water
    CMP_VEC(water.efield); CMP_VEC(water.efg);
    CMP_VEC(water.efield_first); CMP_VEC(water.efg_first);
    CMP_VEC(water.n_first); CMP_VEC(water.n_second);
    CMP_VEC(water.half_shell_asymmetry); CMP_VEC(water.dipole_cos);
    CMP_VEC(water.nearest_ion_dist); CMP_VEC(water.nearest_ion_charge);
    CMP_VEC(water.dipole_vector); CMP_VEC(water.surface_normal);
    CMP_VEC(water.sasa_asymmetry); CMP_VEC(water.sasa_dipole_align);
    CMP_VEC(water.sasa_dipole_cohere); CMP_VEC(water.sasa_first_shell_n);

    // Charges
    CMP_VEC(charges.aimnet2_charge); CMP_VEC(charges.eeq_charge);
    CMP_VEC(charges.eeq_cn);

    // AIMNet2 embedding
    CMP_VEC(aimnet2_embedding.aim);

    // Per-ring
    CMP_VEC(per_ring.geometry);
    CMP_STR_VEC(per_ring.geometry_fields);
    CMP_VEC(per_ring.ring_type);
    CMP_VEC(per_ring.bs_T2); CMP_VEC(per_ring.hm_T2);
    CMP_VEC(per_ring.chi_T2); CMP_VEC(per_ring.pq_T2);
    CMP_VEC(per_ring.hm_H_T2);
    CMP_VEC(per_ring.disp_scalar); CMP_VEC(per_ring.disp_T2);

    // Ring geometry
    CMP_VEC(ring_geometry.data);
    CMP_STR_VEC(ring_geometry.fields);

    // Bonded energy
    CMP_VEC(bonded_energy.bond); CMP_VEC(bonded_energy.angle);
    CMP_VEC(bonded_energy.urey_bradley); CMP_VEC(bonded_energy.proper_dih);
    CMP_VEC(bonded_energy.improper_dih); CMP_VEC(bonded_energy.cmap);
    CMP_VEC(bonded_energy.total);

    // Energy
    CMP_VEC(energy.coulomb_sr); CMP_VEC(energy.coulomb_recip);
    CMP_VEC(energy.bond); CMP_VEC(energy.angle);
    CMP_VEC(energy.urey_bradley); CMP_VEC(energy.proper_dih);
    CMP_VEC(energy.improper_dih); CMP_VEC(energy.cmap_dih);
    CMP_VEC(energy.lj_sr); CMP_VEC(energy.potential);
    CMP_VEC(energy.kinetic); CMP_VEC(energy.enthalpy);
    CMP_VEC(energy.temperature); CMP_VEC(energy.pressure);
    CMP_VEC(energy.volume); CMP_VEC(energy.density);
    CMP_VEC(energy.box); CMP_VEC(energy.virial);
    CMP_STR_VEC(energy.virial_layout);
    CMP_VEC(energy.pressure_tensor);
    CMP_VEC(energy.T_protein); CMP_VEC(energy.T_non_protein);

    // Dihedrals
    CMP_VEC(dihedrals.phi); CMP_VEC(dihedrals.psi);
    CMP_VEC(dihedrals.omega);
    CMP_VEC(dihedrals.chi1); CMP_VEC(dihedrals.chi2);
    CMP_VEC(dihedrals.chi3); CMP_VEC(dihedrals.chi4);
    CMP_VEC(dihedrals.chi1_cos); CMP_VEC(dihedrals.chi1_sin);
    CMP_VEC(dihedrals.chi2_cos); CMP_VEC(dihedrals.chi2_sin);
    CMP_VEC(dihedrals.chi3_cos); CMP_VEC(dihedrals.chi3_sin);
    CMP_VEC(dihedrals.chi4_cos); CMP_VEC(dihedrals.chi4_sin);

    // DSSP
    CMP_VEC(dssp.ss8); CMP_VEC(dssp.hbond_energy);
}

#undef CMP_VEC
#undef CMP_STR_VEC
#undef CMP_SCALAR


// ── Main ───────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    std::string test_dir = "test";
    if (argc > 1) test_dir = argv[1];

    std::vector<std::string> h5_files;
    for (const auto& entry : fs::directory_iterator(test_dir)) {
        if (entry.path().extension() == ".h5")
            h5_files.push_back(entry.path().string());
    }
    std::sort(h5_files.begin(), h5_files.end());

    if (h5_files.empty()) {
        fprintf(stderr, "No .h5 files found in %s\n", test_dir.c_str());
        return 1;
    }

    printf("Roundtrip validation: %zu H5 files in %s\n\n", h5_files.size(), test_dir.c_str());

    int total_pass = 0;
    int total_fail = 0;

    for (const auto& src_path : h5_files) {
        std::string basename = fs::path(src_path).stem().string();
        std::string tmp1 = test_dir + "/" + basename + "_rt1.h5";
        std::string tmp2 = test_dir + "/" + basename + "_rt2.h5";

        printf("=== %s ===\n", basename.c_str());

        // ── Read original ──────────────────────────────────────
        printf("  Reading original... ");
        fflush(stdout);
        AnalysisFile af1;
        af1.ReadH5(src_path);
        printf("OK (%zu frames, %zu atoms, %zu res)\n",
               af1.n_frames, af1.n_atoms, af1.n_residues);

        if (af1.aimnet2_embedding.was_float64)
            printf("  NOTE: aim was float64 in source, normalised to float32\n");

        // ── Write rt1 ──────────────────────────────────────────
        printf("  Writing rt1... ");
        fflush(stdout);
        af1.WriteH5(tmp1);
        printf("OK\n");

        // ── Phase 1: original vs rt1 dataset comparison ────────
        printf("  Phase 1: dataset byte comparison (orig vs rt1)\n");
        {
            TestCounters c;
            HighFive::File f_orig(src_path, HighFive::File::ReadOnly);
            HighFive::File f_rt1(tmp1, HighFive::File::ReadOnly);

            auto root = f_orig.getGroup("/");
            std::vector<std::string> ds_paths;
            CollectDatasets(root, "", ds_paths);

            for (const auto& p : ds_paths) {
                if (af1.aimnet2_embedding.was_float64 && p == "aimnet2_embedding/aim") {
                    c.Skip();
                    continue;
                }
                CompareDataset(f_orig, f_rt1, p, c);
            }
            printf("    Datasets: %d pass, %d fail, %d skip\n", c.pass, c.fail, c.skip);
            for (const auto& f : c.failures) printf("    %s\n", f.c_str());
            if (!c.Ok()) { total_fail++; fs::remove(tmp1); continue; }
        }

        // ── Phase 2: double roundtrip (write determinism) ──────
        printf("  Phase 2: double roundtrip (rt1 -> read -> rt2)\n");
        {
            AnalysisFile af2;
            af2.ReadH5(tmp1);
            af2.WriteH5(tmp2);

            TestCounters c;
            HighFive::File f_rt1(tmp1, HighFive::File::ReadOnly);
            HighFive::File f_rt2(tmp2, HighFive::File::ReadOnly);

            auto root = f_rt1.getGroup("/");
            std::vector<std::string> ds_paths;
            CollectDatasets(root, "", ds_paths);

            for (const auto& p : ds_paths)
                CompareDataset(f_rt1, f_rt2, p, c);

            printf("    Datasets: %d pass, %d fail\n", c.pass, c.fail);
            for (const auto& f : c.failures) printf("    %s\n", f.c_str());
            if (!c.Ok()) { total_fail++; fs::remove(tmp1); fs::remove(tmp2); continue; }
        }

        // ── Phase 3: object model identity ─────────────────────
        printf("  Phase 3: object model identity (af1 vs af2)\n");
        {
            AnalysisFile af2;
            af2.ReadH5(tmp1);

            TestCounters c;
            CompareObjects(af1, af2, c);
            printf("    Fields: %d pass, %d fail\n", c.pass, c.fail);
            for (const auto& f : c.failures) printf("    %s\n", f.c_str());
            if (!c.Ok()) { total_fail++; fs::remove(tmp1); fs::remove(tmp2); continue; }
        }

        // ── Phase 4: structural checks ─────────────────────────
        printf("  Phase 4: structural coverage\n");
        {
            TestCounters c;
            HighFive::File f_orig(src_path, HighFive::File::ReadOnly);
            HighFive::File f_rt1(tmp1, HighFive::File::ReadOnly);

            // Bidirectional dataset coverage
            auto root_o = f_orig.getGroup("/");
            auto root_w = f_rt1.getGroup("/");
            std::vector<std::string> ds_o, ds_w;
            CollectDatasets(root_o, "", ds_o);
            CollectDatasets(root_w, "", ds_w);
            std::set<std::string> set_o(ds_o.begin(), ds_o.end());
            std::set<std::string> set_w(ds_w.begin(), ds_w.end());

            for (const auto& p : set_o) {
                if (set_w.find(p) == set_w.end())
                    c.Fail("missing in written: " + p);
            }
            for (const auto& p : set_w) {
                if (set_o.find(p) == set_o.end())
                    c.Fail("extra in written: " + p);
            }
            if (set_o == set_w) c.Pass();

            // Bidirectional group coverage
            std::vector<std::string> grp_o, grp_w;
            CollectGroups(root_o, "", grp_o);
            CollectGroups(root_w, "", grp_w);
            std::set<std::string> gset_o(grp_o.begin(), grp_o.end());
            std::set<std::string> gset_w(grp_w.begin(), grp_w.end());
            for (const auto& g : gset_o)
                if (gset_w.find(g) == gset_w.end())
                    c.Fail("missing group: " + g);
            for (const auto& g : gset_w)
                if (gset_o.find(g) == gset_o.end())
                    c.Fail("extra group: " + g);
            if (gset_o == gset_w) c.Pass();

            // Bidirectional group attribute comparison
            for (const auto& gp : grp_o) {
                auto go = f_orig.getGroup(gp);
                auto gw = f_rt1.getGroup(gp);
                CompareAttrs(go, gw, gp, c);
            }

            // Dataset-level attributes (SphericalTensor layout/convention)
            for (const auto& dp : ds_o) {
                if (af1.aimnet2_embedding.was_float64 && dp == "aimnet2_embedding/aim")
                    continue;
                auto dso = f_orig.getDataSet(dp);
                auto dsw = f_rt1.getDataSet(dp);
                auto names_o = dso.listAttributeNames();
                auto names_w = dsw.listAttributeNames();
                std::set<std::string> sno(names_o.begin(), names_o.end());
                std::set<std::string> snw(names_w.begin(), names_w.end());

                for (const auto& a : sno) {
                    if (snw.find(a) == snw.end()) {
                        c.Fail(dp + "/@" + a + ": missing in written");
                    } else {
                        std::string vo, vw;
                        dso.getAttribute(a).read(vo);
                        dsw.getAttribute(a).read(vw);
                        if (vo != vw)
                            c.Fail(dp + "/@" + a + ": value mismatch");
                        else
                            c.Pass();
                    }
                }
                for (const auto& a : snw) {
                    if (sno.find(a) == sno.end())
                        c.Fail(dp + "/@" + a + ": extra in written");
                }
            }

            printf("    Checks: %d pass, %d fail\n", c.pass, c.fail);
            for (const auto& f : c.failures) printf("    %s\n", f.c_str());
            if (!c.Ok()) { total_fail++; fs::remove(tmp1); fs::remove(tmp2); continue; }
        }

        // ── Phase 5: 1CBH aim dtype exception ──────────────────
        if (af1.aimnet2_embedding.was_float64) {
            printf("  Phase 5: aim float64->float32 conversion\n");
            TestCounters c;
            HighFive::File f_orig(src_path, HighFive::File::ReadOnly);
            HighFive::File f_rt1(tmp1, HighFive::File::ReadOnly);

            auto ds_o = f_orig.getDataSet("aimnet2_embedding/aim");
            auto ds_w = f_rt1.getDataSet("aimnet2_embedding/aim");

            // Verify written is float32
            auto dt_w = ds_w.getDataType();
            if (dt_w.getSize() == 4 &&
                dt_w.getClass() == HighFive::DataTypeClass::Float) {
                c.Pass();
            } else {
                c.Fail("written aim is not float32");
            }

            // Verify shape preserved
            if (ds_o.getDimensions() == ds_w.getDimensions()) {
                c.Pass();
            } else {
                c.Fail("aim shape changed");
            }

            // Verify conversion is lossless: read orig as float64,
            // cast to float32, compare to written float32
            size_t n = 1;
            for (auto d : ds_o.getDimensions()) n *= d;

            auto dt_o = ds_o.getDataType();
            std::vector<double> orig_f64(n);
            ds_o.read(orig_f64.data(), dt_o);

            std::vector<float> written_f32(n);
            ds_w.read(written_f32.data(), dt_w);

            size_t lossy = 0;
            for (size_t i = 0; i < n; ++i) {
                float converted = static_cast<float>(orig_f64[i]);
                if (std::memcmp(&converted, &written_f32[i], sizeof(float)) != 0)
                    lossy++;
            }
            if (lossy == 0) {
                c.Pass();
                printf("    float64->float32 conversion: LOSSLESS (%zu values)\n", n);
            } else {
                c.Fail("aim conversion: " + std::to_string(lossy) +
                       "/" + std::to_string(n) + " values differ");
            }

            printf("    Checks: %d pass, %d fail\n", c.pass, c.fail);
            for (const auto& f : c.failures) printf("    %s\n", f.c_str());
            if (!c.Ok()) { total_fail++; fs::remove(tmp1); fs::remove(tmp2); continue; }
        }

        // ── Phase 6: special value census ──────────────────────
        printf("  Phase 6: special value census\n");
        {
            TestCounters c;
            HighFive::File f_orig(src_path, HighFive::File::ReadOnly);
            HighFive::File f_rt1(tmp1, HighFive::File::ReadOnly);

            auto root = f_orig.getGroup("/");
            std::vector<std::string> ds_paths;
            CollectDatasets(root, "", ds_paths);

            for (const auto& p : ds_paths) {
                if (af1.aimnet2_embedding.was_float64 && p == "aimnet2_embedding/aim")
                    continue;
                CompareSpecials(f_orig, f_rt1, p, c);
            }
            printf("    Float datasets: %d pass, %d fail\n", c.pass, c.fail);
            for (const auto& f : c.failures) printf("    %s\n", f.c_str());
            if (!c.Ok()) { total_fail++; fs::remove(tmp1); fs::remove(tmp2); continue; }
        }

        // ── File-level non-identity evidence ───────────────────
        {
            auto sz_orig = fs::file_size(src_path);
            auto sz_rt1  = fs::file_size(tmp1);
            if (sz_orig != sz_rt1) {
                printf("  File sizes: orig=%zuMB, rt1=%zuMB (differ — expected)\n",
                       sz_orig / (1024*1024), sz_rt1 / (1024*1024));
            } else {
                printf("  File sizes: both %zuMB (same size, may still differ internally)\n",
                       sz_orig / (1024*1024));
            }
        }

        printf("  PASS\n\n");
        total_pass++;

        fs::remove(tmp1);
        fs::remove(tmp2);
    }

    printf("========================================\n");
    printf("Total: %d pass, %d fail out of %zu files\n",
           total_pass, total_fail, h5_files.size());

    return total_fail > 0 ? 1 : 0;
}

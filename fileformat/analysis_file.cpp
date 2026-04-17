#include "analysis_file.h"

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Group.hpp>

#include <cassert>
#include <cstring>
#include <stdexcept>

// ── Read helpers ────────────────────────────────────────────────────
// These read a dataset into a pre-sized flat vector using read_raw,
// which preserves exact binary content (no type conversion).

namespace {

void ReadDouble(const HighFive::Group& grp, const std::string& name,
                std::vector<double>& out, size_t expected) {
    auto ds = grp.getDataSet(name);
    out.resize(expected);
    ds.read(out.data());
}

void ReadFloat(const HighFive::Group& grp, const std::string& name,
               std::vector<float>& out, size_t expected) {
    auto ds = grp.getDataSet(name);
    out.resize(expected);
    ds.read(out.data());
}

void ReadInt32(const HighFive::Group& grp, const std::string& name,
               std::vector<int32_t>& out, size_t expected) {
    auto ds = grp.getDataSet(name);
    out.resize(expected);
    ds.read(out.data());
}

void ReadInt16(const HighFive::Group& grp, const std::string& name,
               std::vector<int16_t>& out, size_t expected) {
    auto ds = grp.getDataSet(name);
    out.resize(expected);
    ds.read(out.data());
}

void ReadInt8(const HighFive::Group& grp, const std::string& name,
              std::vector<int8_t>& out, size_t expected) {
    auto ds = grp.getDataSet(name);
    out.resize(expected);
    ds.read(out.data());
}

void ReadStrings(const HighFive::Group& grp, const std::string& name,
                 std::vector<std::string>& out) {
    auto ds = grp.getDataSet(name);
    ds.read(out);
}

// ── Write helpers ───────────────────────────────────────────────────

void WriteDouble(HighFive::Group& grp, const std::string& name,
                 const double* data, std::vector<size_t> shape) {
    HighFive::DataSpace space(shape);
    auto ds = grp.createDataSet<double>(name, space);
    ds.write_raw(data);
}

void WriteFloat(HighFive::Group& grp, const std::string& name,
                const float* data, std::vector<size_t> shape) {
    HighFive::DataSpace space(shape);
    auto ds = grp.createDataSet<float>(name, space);
    ds.write_raw(data);
}

void WriteInt32(HighFive::Group& grp, const std::string& name,
                const int32_t* data, std::vector<size_t> shape) {
    HighFive::DataSpace space(shape);
    auto ds = grp.createDataSet<int32_t>(name, space);
    ds.write_raw(data);
}

void WriteInt16(HighFive::Group& grp, const std::string& name,
                const int16_t* data, std::vector<size_t> shape) {
    HighFive::DataSpace space(shape);
    auto ds = grp.createDataSet<int16_t>(name, space);
    ds.write_raw(data);
}

void WriteInt8(HighFive::Group& grp, const std::string& name,
               const int8_t* data, std::vector<size_t> shape) {
    HighFive::DataSpace space(shape);
    auto ds = grp.createDataSet<int8_t>(name, space);
    ds.write_raw(data);
}

// SphericalTensor: (T, N, 9) double with layout + convention attrs.
void WriteSphericalTensor(HighFive::Group& grp, const std::string& name,
                          const double* data, size_t T, size_t N) {
    HighFive::DataSpace space({T, N, size_t{9}});
    auto ds = grp.createDataSet<double>(name, space);
    ds.write_raw(data);
    ds.createAttribute("layout",
        std::string("T0,T1[3],T2[5]"));
    ds.createAttribute("convention",
        std::string("T0=isotropic, T1[0..2]=antisymmetric, "
                     "T2[0..4]=traceless symmetric (sphericart)"));
}

}  // anonymous namespace


// ════════════════════════════════════════════════════════════════════
//  ReadH5
// ════════════════════════════════════════════════════════════════════

void AnalysisFile::ReadH5(const std::string& path) {
    using namespace HighFive;
    File file(path, File::ReadOnly);

    // ── /meta ──────────────────────────────────────────────────
    {
        auto grp = file.getGroup("meta");
        grp.getAttribute("protein_id").read(meta.protein_id);

        // Read dimension attrs as size_t-compatible types.
        // HighFive stores these as the type they were created with;
        // the writer uses size_t, but h5py shows them as int64.
        // We read into size_t via a temp to handle either.
        grp.getAttribute("n_atoms").read(n_atoms);
        grp.getAttribute("n_frames").read(n_frames);
        grp.getAttribute("n_residues").read(n_residues);
        grp.getAttribute("stride").read(meta.stride);

        ReadDouble(grp, "frame_times",   meta.frame_times,   n_frames);
        ReadInt32 (grp, "frame_indices", meta.frame_indices, n_frames);
    }

    const size_t T = n_frames;
    const size_t N = n_atoms;
    const size_t R = n_residues;

    // ── /atoms ─────────────────────────────────────────────────
    {
        auto grp = file.getGroup("atoms");
        ReadInt32  (grp, "element",         atoms.element,         N);
        ReadInt32  (grp, "residue_index",   atoms.residue_index,   N);
        ReadStrings(grp, "atom_name",       atoms.atom_name);
        ReadInt32  (grp, "atom_role",       atoms.atom_role,       N);
        ReadInt32  (grp, "hybridisation",   atoms.hybridisation,   N);
        ReadInt32  (grp, "n_bonded",        atoms.n_bonded,        N);
        ReadInt32  (grp, "graph_dist_ring", atoms.graph_dist_ring, N);
        ReadInt8   (grp, "is_backbone",     atoms.is_backbone,     N);
        ReadInt8   (grp, "is_conjugated",   atoms.is_conjugated,   N);
        // Enrichment
        ReadInt8   (grp, "is_amide_H",           atoms.is_amide_H,           N);
        ReadInt8   (grp, "is_alpha_H",            atoms.is_alpha_H,            N);
        ReadInt8   (grp, "is_methyl",             atoms.is_methyl,             N);
        ReadInt8   (grp, "is_aromatic_H",         atoms.is_aromatic_H,         N);
        ReadInt8   (grp, "is_on_aromatic_residue",atoms.is_on_aromatic_residue,N);
        ReadInt8   (grp, "is_hbond_donor",        atoms.is_hbond_donor,        N);
        ReadInt8   (grp, "is_hbond_acceptor",     atoms.is_hbond_acceptor,     N);
        ReadInt8   (grp, "parent_is_sp2",         atoms.parent_is_sp2,         N);
        ReadInt32  (grp, "graph_dist_N",          atoms.graph_dist_N,          N);
        ReadInt32  (grp, "graph_dist_O",          atoms.graph_dist_O,          N);
        ReadDouble (grp, "eneg_sum_1",            atoms.eneg_sum_1,            N);
        ReadDouble (grp, "eneg_sum_2",            atoms.eneg_sum_2,            N);
        ReadInt32  (grp, "n_pi_bonds_3",          atoms.n_pi_bonds_3,          N);
        ReadInt32  (grp, "bfs_to_nearest_ring",   atoms.bfs_to_nearest_ring,   N);
        ReadDouble (grp, "bfs_decay",             atoms.bfs_decay,             N);
        ReadDouble (grp, "partial_charge",        atoms.partial_charge,        N);
        ReadDouble (grp, "vdw_radius",            atoms.vdw_radius,            N);
    }

    // ── /residues ──────────────────────────────────────────────
    {
        auto grp = file.getGroup("residues");
        ReadStrings(grp, "residue_name",   residues.residue_name);
        ReadInt32  (grp, "residue_number", residues.residue_number, R);
        ReadStrings(grp, "chain_id",       residues.chain_id);
    }

    // ── /topology ──────────────────────────────────────────────
    {
        auto grp = file.getGroup("topology");
        grp.getAttribute("n_bonds").read(n_bonds);
        grp.getAttribute("n_rings").read(n_rings);

        ReadInt32(grp, "bond_atoms",        topology.bond_atoms,        n_bonds * 2);
        ReadInt32(grp, "bond_category",     topology.bond_category,     n_bonds);
        ReadInt32(grp, "bond_order",        topology.bond_order,        n_bonds);
        ReadInt32(grp, "parent_atom_index", topology.parent_atom_index, N);
        ReadInt32(grp, "ring_type",         topology.ring_type,         n_rings);
        ReadInt32(grp, "ring_residue",      topology.ring_residue,      n_rings);
        ReadInt32(grp, "ring_fused_partner",topology.ring_fused_partner,n_rings);
        ReadInt32(grp, "ring_offsets",      topology.ring_offsets,      n_rings + 1);
        // ring_atom_indices: size from the last offset
        size_t n_ring_atoms = 0;
        if (!topology.ring_offsets.empty())
            n_ring_atoms = topology.ring_offsets.back();
        ReadInt32(grp, "ring_atom_indices", topology.ring_atom_indices, n_ring_atoms);
    }

    // ── /positions ─────────────────────────────────────────────
    {
        auto grp = file.getGroup("positions");
        ReadDouble(grp, "xyz", positions.xyz, T * N * 3);
    }

    // ── /ring_current ──────────────────────────────────────────
    {
        auto grp = file.getGroup("ring_current");
        ReadDouble(grp, "bs_T0_per_type",  ring_current.bs_T0_per_type,  T*N*8);
        ReadDouble(grp, "bs_T2_per_type",  ring_current.bs_T2_per_type,  T*N*8*5);
        ReadDouble(grp, "hm_T0_per_type",  ring_current.hm_T0_per_type,  T*N*8);
        ReadDouble(grp, "hm_T2_per_type",  ring_current.hm_T2_per_type,  T*N*8*5);
        ReadDouble(grp, "bs_shielding",    ring_current.bs_shielding,    T*N*9);
        ReadDouble(grp, "hm_shielding",    ring_current.hm_shielding,    T*N*9);
        ReadDouble(grp, "rs_shielding",    ring_current.rs_shielding,    T*N*9);
        ReadInt16 (grp, "n_rings_3A",      ring_current.n_rings_3A,      T*N);
        ReadInt16 (grp, "n_rings_5A",      ring_current.n_rings_5A,      T*N);
        ReadInt16 (grp, "n_rings_8A",      ring_current.n_rings_8A,      T*N);
        ReadInt16 (grp, "n_rings_12A",     ring_current.n_rings_12A,     T*N);
        ReadDouble(grp, "mean_ring_dist",  ring_current.mean_ring_dist,  T*N);
        ReadDouble(grp, "nearest_ring_atom",ring_current.nearest_ring_atom,T*N);
        ReadDouble(grp, "G_iso_exp_sum",   ring_current.G_iso_exp_sum,   T*N);
        ReadDouble(grp, "G_T2_exp_sum",    ring_current.G_T2_exp_sum,    T*N*5);
        ReadDouble(grp, "G_iso_var_8A",    ring_current.G_iso_var_8A,    T*N);
        ReadDouble(grp, "total_B_field",   ring_current.total_B_field,   T*N*3);
    }

    // ── /efg ───────────────────────────────────────────────────
    {
        auto grp = file.getGroup("efg");
        ReadDouble(grp, "coulomb_total",     efg.coulomb_total,     T*N*9);
        ReadDouble(grp, "coulomb_backbone",  efg.coulomb_backbone,  T*N*9);
        ReadDouble(grp, "coulomb_aromatic",  efg.coulomb_aromatic,  T*N*9);
        ReadDouble(grp, "coulomb_shielding", efg.coulomb_shielding, T*N*9);
        ReadDouble(grp, "E_total",           efg.E_total,           T*N*3);
        ReadDouble(grp, "E_backbone",        efg.E_backbone,        T*N*3);
        ReadDouble(grp, "E_sidechain",       efg.E_sidechain,       T*N*3);
        ReadDouble(grp, "E_aromatic",        efg.E_aromatic,        T*N*3);
        ReadDouble(grp, "E_solvent",         efg.E_solvent,         T*N*3);
        ReadDouble(grp, "E_magnitude",       efg.E_magnitude,       T*N);
        ReadDouble(grp, "E_bond_proj",       efg.E_bond_proj,       T*N);
        ReadDouble(grp, "E_backbone_frac",   efg.E_backbone_frac,   T*N);
        ReadDouble(grp, "apbs_efg",          efg.apbs_efg,          T*N*9);
        ReadDouble(grp, "apbs_efield",       efg.apbs_efield,       T*N*3);
        ReadDouble(grp, "aimnet2_total",     efg.aimnet2_total,     T*N*9);
        ReadDouble(grp, "aimnet2_backbone",  efg.aimnet2_backbone,  T*N*9);
        ReadDouble(grp, "aimnet2_aromatic",  efg.aimnet2_aromatic,  T*N*9);
        ReadDouble(grp, "aimnet2_shielding", efg.aimnet2_shielding, T*N*9);
    }

    // ── /bond_aniso ────────────────────────────────────────────
    {
        auto grp = file.getGroup("bond_aniso");
        ReadDouble(grp, "mc_shielding",       bond_aniso.mc_shielding,       T*N*9);
        ReadDouble(grp, "T2_backbone",        bond_aniso.T2_backbone,        T*N*9);
        ReadDouble(grp, "T2_sidechain",       bond_aniso.T2_sidechain,       T*N*9);
        ReadDouble(grp, "T2_aromatic",        bond_aniso.T2_aromatic,        T*N*9);
        ReadDouble(grp, "T2_CO_nearest",      bond_aniso.T2_CO_nearest,      T*N*9);
        ReadDouble(grp, "T2_CN_nearest",      bond_aniso.T2_CN_nearest,      T*N*9);
        ReadDouble(grp, "co_sum",             bond_aniso.co_sum,             T*N);
        ReadDouble(grp, "cn_sum",             bond_aniso.cn_sum,             T*N);
        ReadDouble(grp, "sidechain_sum",      bond_aniso.sidechain_sum,      T*N);
        ReadDouble(grp, "aromatic_sum",       bond_aniso.aromatic_sum,       T*N);
        ReadDouble(grp, "co_nearest",         bond_aniso.co_nearest,         T*N);
        ReadDouble(grp, "nearest_CO_dist",    bond_aniso.nearest_CO_dist,    T*N);
        ReadDouble(grp, "nearest_CN_dist",    bond_aniso.nearest_CN_dist,    T*N);
        ReadDouble(grp, "nearest_CO_midpoint",bond_aniso.nearest_CO_midpoint,T*N*3);
        ReadDouble(grp, "dir_nearest_CO",     bond_aniso.dir_nearest_CO,     T*N*3);
    }

    // ── /quadrupole ────────────────────────────────────────────
    {
        auto grp = file.getGroup("quadrupole");
        ReadDouble(grp, "pq_shielding",    quadrupole.pq_shielding,    T*N*9);
        ReadDouble(grp, "pq_T0_per_type",  quadrupole.pq_T0_per_type,  T*N*8);
        ReadDouble(grp, "pq_T2_per_type",  quadrupole.pq_T2_per_type,  T*N*8*5);
    }

    // ── /dispersion ────────────────────────────────────────────
    {
        auto grp = file.getGroup("dispersion");
        ReadDouble(grp, "disp_shielding",   dispersion.disp_shielding,   T*N*9);
        ReadDouble(grp, "disp_T0_per_type", dispersion.disp_T0_per_type, T*N*8);
        ReadDouble(grp, "disp_T2_per_type", dispersion.disp_T2_per_type, T*N*8*5);
    }

    // ── /hbond ─────────────────────────────────────────────────
    {
        auto grp = file.getGroup("hbond");
        ReadDouble(grp, "hbond_shielding",   hbond.hbond_shielding,   T*N*9);
        ReadDouble(grp, "nearest_spherical", hbond.nearest_spherical, T*N*9);
        ReadDouble(grp, "nearest_dist",      hbond.nearest_dist,      T*N);
        ReadDouble(grp, "nearest_dir",       hbond.nearest_dir,       T*N*3);
        ReadDouble(grp, "inv_d3",            hbond.inv_d3,            T*N);
        ReadInt16 (grp, "count_3_5A",        hbond.count_3_5A,        T*N);
        ReadInt8  (grp, "is_donor",          hbond.is_donor,          T*N);
        ReadInt8  (grp, "is_acceptor",       hbond.is_acceptor,       T*N);
        ReadInt8  (grp, "is_backbone",       hbond.is_backbone,       T*N);
    }

    // ── /sasa ──────────────────────────────────────────────────
    {
        auto grp = file.getGroup("sasa");
        ReadDouble(grp, "sasa",   sasa.sasa,   T*N);
        ReadDouble(grp, "normal", sasa.normal, T*N*3);
    }

    // ── /water ─────────────────────────────────────────────────
    {
        auto grp = file.getGroup("water");
        ReadDouble(grp, "efield",            water.efield,            T*N*3);
        ReadDouble(grp, "efg",               water.efg,               T*N*9);
        ReadDouble(grp, "efield_first",      water.efield_first,      T*N*3);
        ReadDouble(grp, "efg_first",         water.efg_first,         T*N*9);
        ReadInt16 (grp, "n_first",           water.n_first,           T*N);
        ReadInt16 (grp, "n_second",          water.n_second,          T*N);
        ReadDouble(grp, "half_shell_asymmetry", water.half_shell_asymmetry, T*N);
        ReadDouble(grp, "dipole_cos",        water.dipole_cos,        T*N);
        ReadDouble(grp, "nearest_ion_dist",  water.nearest_ion_dist,  T*N);
        ReadDouble(grp, "nearest_ion_charge",water.nearest_ion_charge,T*N);
        ReadDouble(grp, "dipole_vector",     water.dipole_vector,     T*N*3);
        ReadDouble(grp, "surface_normal",    water.surface_normal,    T*N*3);
        ReadDouble(grp, "sasa_asymmetry",    water.sasa_asymmetry,    T*N);
        ReadDouble(grp, "sasa_dipole_align", water.sasa_dipole_align, T*N);
        ReadDouble(grp, "sasa_dipole_cohere",water.sasa_dipole_cohere,T*N);
        ReadInt16 (grp, "sasa_first_shell_n",water.sasa_first_shell_n,T*N);
    }

    // ── /charges ───────────────────────────────────────────────
    {
        auto grp = file.getGroup("charges");
        ReadDouble(grp, "aimnet2_charge", charges.aimnet2_charge, T*N);
        ReadDouble(grp, "eeq_charge",    charges.eeq_charge,    T*N);
        ReadDouble(grp, "eeq_cn",        charges.eeq_cn,        T*N);
    }

    // ── /aimnet2_embedding ─────────────────────────────────────
    {
        auto grp = file.getGroup("aimnet2_embedding");
        auto ds = grp.getDataSet("aim");
        auto dtype = ds.getDataType();
        const size_t total = T * N * AIM_DIMS;

        if (dtype.getSize() == 8) {
            // Source file has float64 (pre-fix).  Read as double,
            // convert to float32 for the canonical representation.
            aimnet2_embedding.was_float64 = true;
            std::vector<double> tmp(total);
            ds.read(tmp.data());
            aimnet2_embedding.aim.resize(total);
            for (size_t i = 0; i < total; ++i)
                aimnet2_embedding.aim[i] = static_cast<float>(tmp[i]);
        } else {
            aimnet2_embedding.was_float64 = false;
            aimnet2_embedding.aim.resize(total);
            ds.read(aimnet2_embedding.aim.data());
        }
    }

    // ── /per_ring ──────────────────────────────────────────────
    {
        auto grp = file.getGroup("per_ring");
        ReadDouble (grp, "geometry",    per_ring.geometry,    T*N*K_RINGS*6);
        ReadStrings(grp, "geometry_fields", per_ring.geometry_fields);
        ReadInt8   (grp, "ring_type",   per_ring.ring_type,   T*N*K_RINGS);
        ReadDouble (grp, "bs_T2",       per_ring.bs_T2,       T*N*K_RINGS*5);
        ReadDouble (grp, "hm_T2",       per_ring.hm_T2,       T*N*K_RINGS*5);
        ReadDouble (grp, "chi_T2",      per_ring.chi_T2,      T*N*K_RINGS*5);
        ReadDouble (grp, "pq_T2",       per_ring.pq_T2,       T*N*K_RINGS*5);
        ReadDouble (grp, "hm_H_T2",     per_ring.hm_H_T2,     T*N*K_RINGS*5);
        ReadDouble (grp, "disp_scalar", per_ring.disp_scalar, T*N*K_RINGS);
        ReadDouble (grp, "disp_T2",     per_ring.disp_T2,     T*N*K_RINGS*5);
    }

    // ── /ring_geometry ─────────────────────────────────────────
    {
        auto grp = file.getGroup("ring_geometry");
        ReadDouble (grp, "data",   ring_geometry.data,   T * n_rings * 7);
        ReadStrings(grp, "fields", ring_geometry.fields);
    }

    // ── /bonded_energy ─────────────────────────────────────────
    {
        auto grp = file.getGroup("bonded_energy");
        ReadDouble(grp, "bond",         bonded_energy.bond,         T*N);
        ReadDouble(grp, "angle",        bonded_energy.angle,        T*N);
        ReadDouble(grp, "urey_bradley", bonded_energy.urey_bradley, T*N);
        ReadDouble(grp, "proper_dih",   bonded_energy.proper_dih,   T*N);
        ReadDouble(grp, "improper_dih", bonded_energy.improper_dih, T*N);
        ReadDouble(grp, "cmap",         bonded_energy.cmap,         T*N);
        ReadDouble(grp, "total",        bonded_energy.total,        T*N);
    }

    // ── /energy ────────────────────────────────────────────────
    {
        auto grp = file.getGroup("energy");
        ReadDouble (grp, "coulomb_sr",     energy.coulomb_sr,     T);
        ReadDouble (grp, "coulomb_recip",  energy.coulomb_recip,  T);
        ReadDouble (grp, "bond",           energy.bond,           T);
        ReadDouble (grp, "angle",          energy.angle,          T);
        ReadDouble (grp, "urey_bradley",   energy.urey_bradley,   T);
        ReadDouble (grp, "proper_dih",     energy.proper_dih,     T);
        ReadDouble (grp, "improper_dih",   energy.improper_dih,   T);
        ReadDouble (grp, "cmap_dih",       energy.cmap_dih,       T);
        ReadDouble (grp, "lj_sr",          energy.lj_sr,          T);
        ReadDouble (grp, "potential",      energy.potential,      T);
        ReadDouble (grp, "kinetic",        energy.kinetic,        T);
        ReadDouble (grp, "enthalpy",       energy.enthalpy,       T);
        ReadDouble (grp, "temperature",    energy.temperature,    T);
        ReadDouble (grp, "pressure",       energy.pressure,       T);
        ReadDouble (grp, "volume",         energy.volume,         T);
        ReadDouble (grp, "density",        energy.density,        T);
        ReadDouble (grp, "box",            energy.box,            T*3);
        ReadDouble (grp, "virial",         energy.virial,         T*9);
        ReadStrings(grp, "virial_layout",  energy.virial_layout);
        ReadDouble (grp, "pressure_tensor",energy.pressure_tensor,T*9);
        ReadDouble (grp, "T_protein",      energy.T_protein,      T);
        ReadDouble (grp, "T_non_protein",  energy.T_non_protein,  T);
    }

    // ── /dihedrals ─────────────────────────────────────────────
    {
        auto grp = file.getGroup("dihedrals");
        ReadDouble(grp, "phi",      dihedrals.phi,      T*R);
        ReadDouble(grp, "psi",      dihedrals.psi,      T*R);
        ReadDouble(grp, "omega",    dihedrals.omega,    T*R);
        ReadDouble(grp, "chi1",     dihedrals.chi1,     T*R);
        ReadDouble(grp, "chi2",     dihedrals.chi2,     T*R);
        ReadDouble(grp, "chi3",     dihedrals.chi3,     T*R);
        ReadDouble(grp, "chi4",     dihedrals.chi4,     T*R);
        ReadDouble(grp, "chi1_cos", dihedrals.chi1_cos, T*R);
        ReadDouble(grp, "chi1_sin", dihedrals.chi1_sin, T*R);
        ReadDouble(grp, "chi2_cos", dihedrals.chi2_cos, T*R);
        ReadDouble(grp, "chi2_sin", dihedrals.chi2_sin, T*R);
        ReadDouble(grp, "chi3_cos", dihedrals.chi3_cos, T*R);
        ReadDouble(grp, "chi3_sin", dihedrals.chi3_sin, T*R);
        ReadDouble(grp, "chi4_cos", dihedrals.chi4_cos, T*R);
        ReadDouble(grp, "chi4_sin", dihedrals.chi4_sin, T*R);
    }

    // ── /dssp ──────────────────────────────────────────────────
    {
        auto grp = file.getGroup("dssp");
        ReadInt8  (grp, "ss8",          dssp.ss8,          T*R);
        ReadDouble(grp, "hbond_energy", dssp.hbond_energy, T*R);
    }
}


// ════════════════════════════════════════════════════════════════════
//  WriteH5
// ════════════════════════════════════════════════════════════════════

void AnalysisFile::WriteH5(const std::string& path) const {
    using namespace HighFive;
    File file(path, File::Truncate);

    const size_t T = n_frames;
    const size_t N = n_atoms;
    const size_t R = n_residues;

    // ── /meta ──────────────────────────────────────────────────
    {
        auto grp = file.createGroup("meta");
        grp.createAttribute("protein_id", meta.protein_id);
        grp.createAttribute("n_atoms",    n_atoms);
        grp.createAttribute("n_frames",   n_frames);
        grp.createAttribute("n_residues", n_residues);
        grp.createAttribute("stride",     meta.stride);
        grp.createDataSet("frame_times",   meta.frame_times);
        grp.createDataSet("frame_indices", meta.frame_indices);
    }

    // ── /atoms ─────────────────────────────────────────────────
    {
        auto grp = file.createGroup("atoms");
        grp.createDataSet("element",         atoms.element);
        grp.createDataSet("residue_index",   atoms.residue_index);
        grp.createDataSet("atom_name",       atoms.atom_name);
        grp.createDataSet("atom_role",       atoms.atom_role);
        grp.createDataSet("hybridisation",   atoms.hybridisation);
        grp.createDataSet("n_bonded",        atoms.n_bonded);
        grp.createDataSet("graph_dist_ring", atoms.graph_dist_ring);
        grp.createDataSet("is_backbone",     atoms.is_backbone);
        grp.createDataSet("is_conjugated",   atoms.is_conjugated);
        grp.createDataSet("is_amide_H",            atoms.is_amide_H);
        grp.createDataSet("is_alpha_H",             atoms.is_alpha_H);
        grp.createDataSet("is_methyl",              atoms.is_methyl);
        grp.createDataSet("is_aromatic_H",          atoms.is_aromatic_H);
        grp.createDataSet("is_on_aromatic_residue", atoms.is_on_aromatic_residue);
        grp.createDataSet("is_hbond_donor",         atoms.is_hbond_donor);
        grp.createDataSet("is_hbond_acceptor",      atoms.is_hbond_acceptor);
        grp.createDataSet("parent_is_sp2",          atoms.parent_is_sp2);
        grp.createDataSet("graph_dist_N",           atoms.graph_dist_N);
        grp.createDataSet("graph_dist_O",           atoms.graph_dist_O);
        grp.createDataSet("eneg_sum_1",             atoms.eneg_sum_1);
        grp.createDataSet("eneg_sum_2",             atoms.eneg_sum_2);
        grp.createDataSet("n_pi_bonds_3",           atoms.n_pi_bonds_3);
        grp.createDataSet("bfs_to_nearest_ring",    atoms.bfs_to_nearest_ring);
        grp.createDataSet("bfs_decay",              atoms.bfs_decay);
        grp.createDataSet("partial_charge",         atoms.partial_charge);
        grp.createDataSet("vdw_radius",             atoms.vdw_radius);
    }

    // ── /residues ──────────────────────────────────────────────
    {
        auto grp = file.createGroup("residues");
        grp.createDataSet("residue_name",   residues.residue_name);
        grp.createDataSet("residue_number", residues.residue_number);
        grp.createDataSet("chain_id",       residues.chain_id);
    }

    // ── /topology ──────────────────────────────────────────────
    {
        auto grp = file.createGroup("topology");
        grp.createAttribute("n_bonds", n_bonds);
        grp.createAttribute("n_rings", n_rings);
        WriteInt32(grp, "bond_atoms", topology.bond_atoms.data(), {n_bonds, size_t{2}});
        grp.createDataSet("bond_category",     topology.bond_category);
        grp.createDataSet("bond_order",        topology.bond_order);
        grp.createDataSet("parent_atom_index", topology.parent_atom_index);
        grp.createDataSet("ring_type",         topology.ring_type);
        grp.createDataSet("ring_residue",      topology.ring_residue);
        grp.createDataSet("ring_fused_partner",topology.ring_fused_partner);
        grp.createDataSet("ring_offsets",      topology.ring_offsets);
        grp.createDataSet("ring_atom_indices", topology.ring_atom_indices);
    }

    // ── /positions ─────────────────────────────────────────────
    {
        auto grp = file.createGroup("positions");
        WriteDouble(grp, "xyz", positions.xyz.data(), {T, N, size_t{3}});
    }

    // ── /ring_current ──────────────────────────────────────────
    {
        auto grp = file.createGroup("ring_current");
        WriteDouble(grp, "bs_T0_per_type", ring_current.bs_T0_per_type.data(), {T,N,size_t{8}});
        WriteDouble(grp, "bs_T2_per_type", ring_current.bs_T2_per_type.data(), {T,N,size_t{8},size_t{5}});
        WriteDouble(grp, "hm_T0_per_type", ring_current.hm_T0_per_type.data(), {T,N,size_t{8}});
        WriteDouble(grp, "hm_T2_per_type", ring_current.hm_T2_per_type.data(), {T,N,size_t{8},size_t{5}});
        WriteSphericalTensor(grp, "bs_shielding", ring_current.bs_shielding.data(), T, N);
        WriteSphericalTensor(grp, "hm_shielding", ring_current.hm_shielding.data(), T, N);
        WriteSphericalTensor(grp, "rs_shielding", ring_current.rs_shielding.data(), T, N);
        WriteInt16(grp, "n_rings_3A",  ring_current.n_rings_3A.data(),  {T,N});
        WriteInt16(grp, "n_rings_5A",  ring_current.n_rings_5A.data(),  {T,N});
        WriteInt16(grp, "n_rings_8A",  ring_current.n_rings_8A.data(),  {T,N});
        WriteInt16(grp, "n_rings_12A", ring_current.n_rings_12A.data(), {T,N});
        WriteDouble(grp, "mean_ring_dist",   ring_current.mean_ring_dist.data(),   {T,N});
        WriteDouble(grp, "nearest_ring_atom",ring_current.nearest_ring_atom.data(),{T,N});
        WriteDouble(grp, "G_iso_exp_sum",    ring_current.G_iso_exp_sum.data(),    {T,N});
        WriteDouble(grp, "G_T2_exp_sum",     ring_current.G_T2_exp_sum.data(),     {T,N,size_t{5}});
        WriteDouble(grp, "G_iso_var_8A",     ring_current.G_iso_var_8A.data(),     {T,N});
        WriteDouble(grp, "total_B_field",    ring_current.total_B_field.data(),    {T,N,size_t{3}});
    }

    // ── /efg ───────────────────────────────────────────────────
    {
        auto grp = file.createGroup("efg");
        WriteSphericalTensor(grp, "coulomb_total",     efg.coulomb_total.data(),     T, N);
        WriteSphericalTensor(grp, "coulomb_backbone",  efg.coulomb_backbone.data(),  T, N);
        WriteSphericalTensor(grp, "coulomb_aromatic",  efg.coulomb_aromatic.data(),  T, N);
        WriteSphericalTensor(grp, "coulomb_shielding", efg.coulomb_shielding.data(), T, N);
        WriteDouble(grp, "E_total",         efg.E_total.data(),         {T,N,size_t{3}});
        WriteDouble(grp, "E_backbone",      efg.E_backbone.data(),      {T,N,size_t{3}});
        WriteDouble(grp, "E_sidechain",     efg.E_sidechain.data(),     {T,N,size_t{3}});
        WriteDouble(grp, "E_aromatic",      efg.E_aromatic.data(),      {T,N,size_t{3}});
        WriteDouble(grp, "E_solvent",       efg.E_solvent.data(),       {T,N,size_t{3}});
        WriteDouble(grp, "E_magnitude",     efg.E_magnitude.data(),     {T,N});
        WriteDouble(grp, "E_bond_proj",     efg.E_bond_proj.data(),     {T,N});
        WriteDouble(grp, "E_backbone_frac", efg.E_backbone_frac.data(), {T,N});
        WriteSphericalTensor(grp, "apbs_efg",          efg.apbs_efg.data(),          T, N);
        WriteDouble(grp, "apbs_efield",     efg.apbs_efield.data(),     {T,N,size_t{3}});
        WriteSphericalTensor(grp, "aimnet2_total",     efg.aimnet2_total.data(),     T, N);
        WriteSphericalTensor(grp, "aimnet2_backbone",  efg.aimnet2_backbone.data(),  T, N);
        WriteSphericalTensor(grp, "aimnet2_aromatic",  efg.aimnet2_aromatic.data(),  T, N);
        WriteSphericalTensor(grp, "aimnet2_shielding", efg.aimnet2_shielding.data(), T, N);
    }

    // ── /bond_aniso ────────────────────────────────────────────
    {
        auto grp = file.createGroup("bond_aniso");
        WriteSphericalTensor(grp, "mc_shielding",  bond_aniso.mc_shielding.data(),  T, N);
        WriteSphericalTensor(grp, "T2_backbone",   bond_aniso.T2_backbone.data(),   T, N);
        WriteSphericalTensor(grp, "T2_sidechain",  bond_aniso.T2_sidechain.data(),  T, N);
        WriteSphericalTensor(grp, "T2_aromatic",   bond_aniso.T2_aromatic.data(),   T, N);
        WriteSphericalTensor(grp, "T2_CO_nearest", bond_aniso.T2_CO_nearest.data(), T, N);
        WriteSphericalTensor(grp, "T2_CN_nearest", bond_aniso.T2_CN_nearest.data(), T, N);
        WriteDouble(grp, "co_sum",        bond_aniso.co_sum.data(),        {T,N});
        WriteDouble(grp, "cn_sum",        bond_aniso.cn_sum.data(),        {T,N});
        WriteDouble(grp, "sidechain_sum", bond_aniso.sidechain_sum.data(), {T,N});
        WriteDouble(grp, "aromatic_sum",  bond_aniso.aromatic_sum.data(),  {T,N});
        WriteDouble(grp, "co_nearest",    bond_aniso.co_nearest.data(),    {T,N});
        WriteDouble(grp, "nearest_CO_dist",    bond_aniso.nearest_CO_dist.data(),    {T,N});
        WriteDouble(grp, "nearest_CN_dist",    bond_aniso.nearest_CN_dist.data(),    {T,N});
        WriteDouble(grp, "nearest_CO_midpoint",bond_aniso.nearest_CO_midpoint.data(),{T,N,size_t{3}});
        WriteDouble(grp, "dir_nearest_CO",     bond_aniso.dir_nearest_CO.data(),     {T,N,size_t{3}});
    }

    // ── /quadrupole ────────────────────────────────────────────
    {
        auto grp = file.createGroup("quadrupole");
        WriteSphericalTensor(grp, "pq_shielding", quadrupole.pq_shielding.data(), T, N);
        WriteDouble(grp, "pq_T0_per_type", quadrupole.pq_T0_per_type.data(), {T,N,size_t{8}});
        WriteDouble(grp, "pq_T2_per_type", quadrupole.pq_T2_per_type.data(), {T,N,size_t{8},size_t{5}});
    }

    // ── /dispersion ────────────────────────────────────────────
    {
        auto grp = file.createGroup("dispersion");
        WriteSphericalTensor(grp, "disp_shielding", dispersion.disp_shielding.data(), T, N);
        WriteDouble(grp, "disp_T0_per_type", dispersion.disp_T0_per_type.data(), {T,N,size_t{8}});
        WriteDouble(grp, "disp_T2_per_type", dispersion.disp_T2_per_type.data(), {T,N,size_t{8},size_t{5}});
    }

    // ── /hbond ─────────────────────────────────────────────────
    {
        auto grp = file.createGroup("hbond");
        WriteSphericalTensor(grp, "hbond_shielding",   hbond.hbond_shielding.data(),   T, N);
        WriteSphericalTensor(grp, "nearest_spherical",  hbond.nearest_spherical.data(),  T, N);
        WriteDouble(grp, "nearest_dist", hbond.nearest_dist.data(), {T,N});
        WriteDouble(grp, "nearest_dir",  hbond.nearest_dir.data(),  {T,N,size_t{3}});
        WriteDouble(grp, "inv_d3",       hbond.inv_d3.data(),       {T,N});
        WriteInt16 (grp, "count_3_5A",   hbond.count_3_5A.data(),   {T,N});
        WriteInt8  (grp, "is_donor",     hbond.is_donor.data(),     {T,N});
        WriteInt8  (grp, "is_acceptor",  hbond.is_acceptor.data(),  {T,N});
        WriteInt8  (grp, "is_backbone",  hbond.is_backbone.data(),  {T,N});
    }

    // ── /sasa ──────────────────────────────────────────────────
    {
        auto grp = file.createGroup("sasa");
        WriteDouble(grp, "sasa",   sasa.sasa.data(),   {T,N});
        WriteDouble(grp, "normal", sasa.normal.data(), {T,N,size_t{3}});
    }

    // ── /water ─────────────────────────────────────────────────
    {
        auto grp = file.createGroup("water");
        WriteDouble(grp, "efield",            water.efield.data(),            {T,N,size_t{3}});
        WriteSphericalTensor(grp, "efg",      water.efg.data(),               T, N);
        WriteDouble(grp, "efield_first",      water.efield_first.data(),      {T,N,size_t{3}});
        WriteSphericalTensor(grp, "efg_first",water.efg_first.data(),         T, N);
        WriteInt16 (grp, "n_first",           water.n_first.data(),           {T,N});
        WriteInt16 (grp, "n_second",          water.n_second.data(),          {T,N});
        WriteDouble(grp, "half_shell_asymmetry", water.half_shell_asymmetry.data(), {T,N});
        WriteDouble(grp, "dipole_cos",        water.dipole_cos.data(),        {T,N});
        WriteDouble(grp, "nearest_ion_dist",  water.nearest_ion_dist.data(),  {T,N});
        WriteDouble(grp, "nearest_ion_charge",water.nearest_ion_charge.data(),{T,N});
        WriteDouble(grp, "dipole_vector",     water.dipole_vector.data(),     {T,N,size_t{3}});
        WriteDouble(grp, "surface_normal",    water.surface_normal.data(),    {T,N,size_t{3}});
        WriteDouble(grp, "sasa_asymmetry",    water.sasa_asymmetry.data(),    {T,N});
        WriteDouble(grp, "sasa_dipole_align", water.sasa_dipole_align.data(), {T,N});
        WriteDouble(grp, "sasa_dipole_cohere",water.sasa_dipole_cohere.data(),{T,N});
        WriteInt16 (grp, "sasa_first_shell_n",water.sasa_first_shell_n.data(),{T,N});
    }

    // ── /charges ───────────────────────────────────────────────
    {
        auto grp = file.createGroup("charges");
        WriteDouble(grp, "aimnet2_charge", charges.aimnet2_charge.data(), {T,N});
        WriteDouble(grp, "eeq_charge",    charges.eeq_charge.data(),    {T,N});
        WriteDouble(grp, "eeq_cn",        charges.eeq_cn.data(),        {T,N});
    }

    // ── /aimnet2_embedding ─────────────────────────────────────
    {
        auto grp = file.createGroup("aimnet2_embedding");
        WriteFloat(grp, "aim", aimnet2_embedding.aim.data(), {T, N, AIM_DIMS});
    }

    // ── /per_ring ──────────────────────────────────────────────
    {
        auto grp = file.createGroup("per_ring");
        grp.createAttribute("K", K_RINGS);
        WriteDouble(grp, "geometry", per_ring.geometry.data(), {T,N,K_RINGS,size_t{6}});
        grp.createDataSet("geometry_fields", per_ring.geometry_fields);
        WriteInt8  (grp, "ring_type", per_ring.ring_type.data(), {T,N,K_RINGS});
        WriteDouble(grp, "bs_T2",    per_ring.bs_T2.data(),    {T,N,K_RINGS,size_t{5}});
        WriteDouble(grp, "hm_T2",    per_ring.hm_T2.data(),    {T,N,K_RINGS,size_t{5}});
        WriteDouble(grp, "chi_T2",   per_ring.chi_T2.data(),   {T,N,K_RINGS,size_t{5}});
        WriteDouble(grp, "pq_T2",    per_ring.pq_T2.data(),    {T,N,K_RINGS,size_t{5}});
        WriteDouble(grp, "hm_H_T2",  per_ring.hm_H_T2.data(),  {T,N,K_RINGS,size_t{5}});
        WriteDouble(grp, "disp_scalar",per_ring.disp_scalar.data(),{T,N,K_RINGS});
        WriteDouble(grp, "disp_T2",  per_ring.disp_T2.data(),  {T,N,K_RINGS,size_t{5}});
    }

    // ── /ring_geometry ─────────────────────────────────────────
    {
        auto grp = file.createGroup("ring_geometry");
        WriteDouble(grp, "data", ring_geometry.data.data(), {T, n_rings, size_t{7}});
        grp.createDataSet("fields", ring_geometry.fields);
    }

    // ── /bonded_energy ─────────────────────────────────────────
    {
        auto grp = file.createGroup("bonded_energy");
        grp.createAttribute("units",
            std::string("kJ/mol (split evenly among participating atoms)"));
        WriteDouble(grp, "bond",         bonded_energy.bond.data(),         {T,N});
        WriteDouble(grp, "angle",        bonded_energy.angle.data(),        {T,N});
        WriteDouble(grp, "urey_bradley", bonded_energy.urey_bradley.data(), {T,N});
        WriteDouble(grp, "proper_dih",   bonded_energy.proper_dih.data(),   {T,N});
        WriteDouble(grp, "improper_dih", bonded_energy.improper_dih.data(), {T,N});
        WriteDouble(grp, "cmap",         bonded_energy.cmap.data(),         {T,N});
        WriteDouble(grp, "total",        bonded_energy.total.data(),        {T,N});
    }

    // ── /energy ────────────────────────────────────────────────
    {
        auto grp = file.createGroup("energy");
        grp.createAttribute("units_energy",      std::string("kJ/mol"));
        grp.createAttribute("units_pressure",    std::string("bar"));
        grp.createAttribute("units_temperature", std::string("K"));
        grp.createAttribute("units_volume",      std::string("nm^3"));
        grp.createAttribute("units_density",     std::string("kg/m^3"));
        grp.createAttribute("units_box",         std::string("nm"));
        grp.createDataSet("coulomb_sr",     energy.coulomb_sr);
        grp.createDataSet("coulomb_recip",  energy.coulomb_recip);
        grp.createDataSet("bond",           energy.bond);
        grp.createDataSet("angle",          energy.angle);
        grp.createDataSet("urey_bradley",   energy.urey_bradley);
        grp.createDataSet("proper_dih",     energy.proper_dih);
        grp.createDataSet("improper_dih",   energy.improper_dih);
        grp.createDataSet("cmap_dih",       energy.cmap_dih);
        grp.createDataSet("lj_sr",          energy.lj_sr);
        grp.createDataSet("potential",      energy.potential);
        grp.createDataSet("kinetic",        energy.kinetic);
        grp.createDataSet("enthalpy",       energy.enthalpy);
        grp.createDataSet("temperature",    energy.temperature);
        grp.createDataSet("pressure",       energy.pressure);
        grp.createDataSet("volume",         energy.volume);
        grp.createDataSet("density",        energy.density);
        WriteDouble(grp, "box",    energy.box.data(),    {T,size_t{3}});
        WriteDouble(grp, "virial", energy.virial.data(), {T,size_t{9}});
        grp.createDataSet("virial_layout", energy.virial_layout);
        WriteDouble(grp, "pressure_tensor", energy.pressure_tensor.data(), {T,size_t{9}});
        grp.createDataSet("T_protein",     energy.T_protein);
        grp.createDataSet("T_non_protein", energy.T_non_protein);
    }

    // ── /dihedrals ─────────────────────────────────────────────
    {
        auto grp = file.createGroup("dihedrals");
        grp.createAttribute("units", std::string("radians"));
        grp.createAttribute("omega_convention",
            std::string("CA(i)-C(i)-N(i+1)-CA(i+1)"));
        WriteDouble(grp, "phi",      dihedrals.phi.data(),      {T,R});
        WriteDouble(grp, "psi",      dihedrals.psi.data(),      {T,R});
        WriteDouble(grp, "omega",    dihedrals.omega.data(),    {T,R});
        WriteDouble(grp, "chi1",     dihedrals.chi1.data(),     {T,R});
        WriteDouble(grp, "chi2",     dihedrals.chi2.data(),     {T,R});
        WriteDouble(grp, "chi3",     dihedrals.chi3.data(),     {T,R});
        WriteDouble(grp, "chi4",     dihedrals.chi4.data(),     {T,R});
        WriteDouble(grp, "chi1_cos", dihedrals.chi1_cos.data(), {T,R});
        WriteDouble(grp, "chi1_sin", dihedrals.chi1_sin.data(), {T,R});
        WriteDouble(grp, "chi2_cos", dihedrals.chi2_cos.data(), {T,R});
        WriteDouble(grp, "chi2_sin", dihedrals.chi2_sin.data(), {T,R});
        WriteDouble(grp, "chi3_cos", dihedrals.chi3_cos.data(), {T,R});
        WriteDouble(grp, "chi3_sin", dihedrals.chi3_sin.data(), {T,R});
        WriteDouble(grp, "chi4_cos", dihedrals.chi4_cos.data(), {T,R});
        WriteDouble(grp, "chi4_sin", dihedrals.chi4_sin.data(), {T,R});
    }

    // ── /dssp ──────────────────────────────────────────────────
    {
        auto grp = file.createGroup("dssp");
        grp.createAttribute("ss8_encoding",
            std::string("0=H(alpha),1=G(310),2=I(pi),3=E(strand),"
                         "4=B(bridge),5=T(turn),6=S(bend),7=C(coil)"));
        WriteInt8  (grp, "ss8",          dssp.ss8.data(),          {T,R});
        WriteDouble(grp, "hbond_energy", dssp.hbond_energy.data(), {T,R});
    }
}

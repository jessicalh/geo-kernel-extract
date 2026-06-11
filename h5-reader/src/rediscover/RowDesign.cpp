#include "RowDesign.h"

#include "RowDesignCatalogCoverage.h"

#include <QStringList>

namespace h5reader::rediscover {

namespace {

RowColumnSpec col(const char* name,
                  RowColType type = RowColType::Double,
                  const char* unit = "",
                  const char* irrep = "",
                  RowNativeAxis axis = RowNativeAxis::Row,
                  bool timeVarying = false,
                  bool isFeature = true) {
    RowColumnSpec c;
    c.name = QString::fromLatin1(name);
    c.type = type;
    c.unit = QString::fromLatin1(unit);
    c.irrep = QString::fromLatin1(irrep);
    c.nativeAxis = axis;
    c.timeVarying = timeVarying;
    c.isFeature = isFeature;
    return c;
}

std::vector<RowColumnSpec> buildSchema() {
    std::vector<RowColumnSpec> s = {
        col("dataset_id", RowColType::String, "", "", RowNativeAxis::Protein, false, false),
        col("protein_id", RowColType::String, "", "", RowNativeAxis::Protein, false, false),
        col("case_id", RowColType::String, "", "", RowNativeAxis::Row, false, false),
        col("pose_kind", RowColType::String, "", "", RowNativeAxis::Frame, false, false),
        col("frame_slot", RowColType::Int, "", "", RowNativeAxis::Frame, true, false),
        col("split_group_id", RowColType::String, "", "", RowNativeAxis::Protein, false, false),
        col("atom_index", RowColType::Int, "", "", RowNativeAxis::Atom, false, false),
        col("h5_row", RowColType::Int, "", "", RowNativeAxis::Frame, true, false),
        col("original_index", RowColType::Int, "", "", RowNativeAxis::Frame, true, false),
        col("time_ps", RowColType::Double, "ps", "0e", RowNativeAxis::Frame, true, false),
        col("element", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("residue_index", RowColType::Int, "", "", RowNativeAxis::Residue),
        col("residue_number", RowColType::Int, "", "", RowNativeAxis::Residue),
        col("residue_type", RowColType::Int, "", "", RowNativeAxis::Residue),
        col("atom_name", RowColType::String, "", "", RowNativeAxis::Atom),
        col("iupac_role", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("backbone_role", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("locant", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("branch_outer", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("branch_inner", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("di_index", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("prochiral", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("ring_position_primary", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("ring_position_secondary", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("equivalence_class", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("formal_charge", RowColType::Int, "e", "0e", RowNativeAxis::Atom),
        col("is_exchangeable", RowColType::Bool, "", "", RowNativeAxis::Atom),
        col("polar_h_kind", RowColType::Int, "", "", RowNativeAxis::Atom),
        col("aromatic", RowColType::Bool, "", "", RowNativeAxis::Atom),
        col("residue_class", RowColType::Int, "", "", RowNativeAxis::Residue),
        col("prev_class", RowColType::Int, "", "", RowNativeAxis::Residue),
        col("next_class", RowColType::Int, "", "", RowNativeAxis::Residue),
        col("prev_restype", RowColType::Int, "", "", RowNativeAxis::Residue),
        col("next_restype", RowColType::Int, "", "", RowNativeAxis::Residue),
        col("pre_proline", RowColType::Bool, "", "", RowNativeAxis::Residue),
        col("is_pro", RowColType::Bool, "", "", RowNativeAxis::Residue),
        col("is_gly", RowColType::Bool, "", "", RowNativeAxis::Residue),
        col("terminal_state", RowColType::Int, "", "", RowNativeAxis::Residue),
        col("phi", RowColType::Double, "rad", "0e", RowNativeAxis::Residue, true),
        col("psi", RowColType::Double, "rad", "0e", RowNativeAxis::Residue, true),
        col("sin_phi", RowColType::Double, "", "0e", RowNativeAxis::Residue, true),
        col("cos_phi", RowColType::Double, "", "0e", RowNativeAxis::Residue, true),
        col("sin_psi", RowColType::Double, "", "0e", RowNativeAxis::Residue, true),
        col("cos_psi", RowColType::Double, "", "0e", RowNativeAxis::Residue, true),
        col("rama_region", RowColType::Int, "", "", RowNativeAxis::Residue, true),
    };
    for (int k = 1; k <= 4; ++k) {
        s.push_back(col(QStringLiteral("sin_chi%1").arg(k).toLatin1().constData(), RowColType::Double, "", "0e", RowNativeAxis::Residue, true));
        s.push_back(col(QStringLiteral("cos_chi%1").arg(k).toLatin1().constData(), RowColType::Double, "", "0e", RowNativeAxis::Residue, true));
        s.push_back(col(QStringLiteral("chi%1_exists").arg(k).toLatin1().constData(), RowColType::Bool, "", "", RowNativeAxis::Residue, true));
    }
    s.push_back(col("omega", RowColType::Double, "rad", "0e", RowNativeAxis::Residue, true));
    s.push_back(col("sin_omega", RowColType::Double, "", "0e", RowNativeAxis::Residue, true));
    s.push_back(col("cos_omega", RowColType::Double, "", "0e", RowNativeAxis::Residue, true));
    s.push_back(col("pucker_Q", RowColType::Double, "", "0e", RowNativeAxis::Ring));
    s.push_back(col("pucker_theta", RowColType::Double, "degrees", "0e", RowNativeAxis::Ring));
    s.push_back(col("pyramidalization", RowColType::Double, "", "0e", RowNativeAxis::Atom));
    s.push_back(col("aromatic_ringflip_state", RowColType::Double, "rad", "0e", RowNativeAxis::Ring));
    s.push_back(col("dssp_ss8", RowColType::Int, "", "", RowNativeAxis::Residue, true));
    s.push_back(col("dssp_ss3", RowColType::Int, "", "", RowNativeAxis::Residue, true));
    s.push_back(col("dssp_hbond_energy", RowColType::Double, "kcal/mol", "0e", RowNativeAxis::Residue, true));
    s.push_back(col("hbond_donor", RowColType::Bool, "", "", RowNativeAxis::Residue, true));
    s.push_back(col("hbond_acceptor", RowColType::Bool, "", "", RowNativeAxis::Residue, true));
    for (int cutoff : {4, 6, 8, 10}) {
        s.push_back(col(QStringLiteral("n_atoms_%1A").arg(cutoff).toLatin1().constData(), RowColType::Int, "", "", RowNativeAxis::Atom, true));
        s.push_back(col(QStringLiteral("n_rings_%1A").arg(cutoff).toLatin1().constData(), RowColType::Int, "", "", RowNativeAxis::Ring, true));
        s.push_back(col(QStringLiteral("n_charges_%1A").arg(cutoff).toLatin1().constData(), RowColType::Int, "", "", RowNativeAxis::Atom, true));
        s.push_back(col(QStringLiteral("n_bonds_%1A").arg(cutoff).toLatin1().constData(), RowColType::Int, "", "", RowNativeAxis::Bond, true));
    }
    const char* scalarFeaturesBeforeRing[] = {
        "nearest_ring_dist", "nearest_ring_identity_ord", "nearest_charge_dist",
        "nearest_charge_identity_ord", "nearest_charge_sign", "nearest_bond_dist",
        "nearest_bond_identity_ord", "nearest_atom_dist", "nearest_atom_identity_ord",
        "ring_cyl_z", "ring_cyl_rho", "ring_angle_to_normal", "ring_cos_phi", "ring_sin_phi",
        "self_or_bonded_atom_count", "self_or_bonded_bond_count",
        "ring_bs_T0", "ring_bs_absT2", "ring_hm_T0", "ring_hm_absT2",
    };
    for (const char* name : scalarFeaturesBeforeRing)
        s.push_back(col(name, RowColType::Double, "", "0e", RowNativeAxis::Atom, true));

    s.push_back(col("ring_jb_T0", RowColType::Double, "ppm", "0e", RowNativeAxis::Atom, true));
    for (int i = 0; i < 3; ++i)
        s.push_back(col(QStringLiteral("ring_jb_T1_%1").arg(i).toLatin1().constData(),
                        RowColType::Double, "ppm", "1e", RowNativeAxis::Atom, true));
    for (int i = 0; i < 5; ++i)
        s.push_back(col(QStringLiteral("ring_jb_T2_%1").arg(i).toLatin1().constData(),
                        RowColType::Double, "ppm", "2e", RowNativeAxis::Atom, true));
    s.push_back(col("ring_jb_absT2", RowColType::Double, "ppm", "0e", RowNativeAxis::Atom, true));

    s.push_back(col("ring_hm_jb_T0", RowColType::Double, "ppm", "0e", RowNativeAxis::Atom, true));
    for (int i = 0; i < 3; ++i)
        s.push_back(col(QStringLiteral("ring_hm_jb_T1_%1").arg(i).toLatin1().constData(),
                        RowColType::Double, "ppm", "1e", RowNativeAxis::Atom, true));
    for (int i = 0; i < 5; ++i)
        s.push_back(col(QStringLiteral("ring_hm_jb_T2_%1").arg(i).toLatin1().constData(),
                        RowColType::Double, "ppm", "2e", RowNativeAxis::Atom, true));
    s.push_back(col("ring_hm_jb_absT2", RowColType::Double, "ppm", "0e", RowNativeAxis::Atom, true));

    const char* attribution[] = {
        "ring_jb_T0_phe", "ring_jb_T0_tyr", "ring_jb_T0_trp6", "ring_jb_T0_trp5",
        "ring_jb_T0_trp9", "ring_jb_T0_his", "ring_jb_T0_hid", "ring_jb_T0_hie",
    };
    for (const char* name : attribution)
        s.push_back(col(name, RowColType::Double, "ppm", "0e", RowNativeAxis::Atom, true));

    const char* scalarFeaturesAfterRing[] = {
        "mc_lit_T0", "mc_lit_absT2",
        "mopac_bare_efg_kernel_absT2", "apbs_bare_efg_kernel_absT2",
        "aimnet2_bare_efg_kernel_absT2", "ff14sb_efield_mag",
        "mopac_efield_mag", "apbs_efield_mag", "aimnet2_charge",
        "ff14sb_charge", "mopac_charge", "eeq_charge", "eeq_cn", "sasa",
        "larsen_hbond_count", "larsen_hbond_absT2", "target_T0"
    };
    for (const char* name : scalarFeaturesAfterRing)
        s.push_back(col(name, RowColType::Double, "", "0e", RowNativeAxis::Atom, true));
    s.push_back(col("tensor_frame", RowColType::String, "", "", RowNativeAxis::Row, false, false));
    s.push_back(col("valid_for_T2_model", RowColType::Bool, "", "", RowNativeAxis::Row, false, false));
    s.push_back(col("region_def_id", RowColType::String, "", "", RowNativeAxis::Row, false, false));
    s.push_back(col("rama_region_hdr", RowColType::String, "", "", RowNativeAxis::Residue, true, false));
    s.push_back(col("rotamer_id", RowColType::String, "", "", RowNativeAxis::Residue, true, false));
    for (const CatalogRowColumn& c : CatalogRowColumns()) s.push_back(c.spec);
    return s;
}

}  // namespace

ResidueClass ClassifyResidue(model::AminoAcid aa) {
    using model::AminoAcid;
    switch (aa) {
    case AminoAcid::GLY: return ResidueClass::Glycine;
    case AminoAcid::PRO: return ResidueClass::Proline;
    case AminoAcid::CYS: return ResidueClass::Cys;
    case AminoAcid::PHE:
    case AminoAcid::TYR:
    case AminoAcid::TRP:
    case AminoAcid::HIS: return ResidueClass::Aromatic;
    case AminoAcid::ARG:
    case AminoAcid::LYS: return ResidueClass::PosCharged;
    case AminoAcid::ASP:
    case AminoAcid::GLU: return ResidueClass::NegCharged;
    case AminoAcid::SER:
    case AminoAcid::THR:
    case AminoAcid::ASN:
    case AminoAcid::GLN: return ResidueClass::Polar;
    case AminoAcid::ALA:
    case AminoAcid::VAL:
    case AminoAcid::LEU:
    case AminoAcid::ILE:
    case AminoAcid::MET: return ResidueClass::Hydrophobic;
    default: return ResidueClass::Other;
    }
}

const char* NameForResidueClass(ResidueClass c) {
    switch (c) {
    case ResidueClass::Hydrophobic: return "hydrophobic";
    case ResidueClass::Polar: return "polar";
    case ResidueClass::Aromatic: return "aromatic";
    case ResidueClass::PosCharged: return "pos_charged";
    case ResidueClass::NegCharged: return "neg_charged";
    case ResidueClass::Glycine: return "glycine";
    case ResidueClass::Proline: return "proline";
    case ResidueClass::Cys: return "cys";
    case ResidueClass::Other: return "other";
    }
    return "other";
}

const std::vector<RowColumnSpec>& RowDesignSchema() {
    static const std::vector<RowColumnSpec> schema = buildSchema();
    return schema;
}

QStringList RowDesignHeader() {
    QStringList out;
    for (const RowColumnSpec& c : RowDesignSchema()) out << c.name;
    return out;
}

int RowDesignColumnIndex(const QString& name) {
    const auto& schema = RowDesignSchema();
    for (std::size_t i = 0; i < schema.size(); ++i) {
        if (schema[i].name == name) return static_cast<int>(i);
    }
    return -1;
}

}  // namespace h5reader::rediscover

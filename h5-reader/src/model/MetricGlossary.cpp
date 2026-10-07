#include "MetricGlossary.h"

#include "MetricTaxonomy.h"

#include <QLatin1String>
#include <QStringList>

namespace h5reader::model {

namespace {

QString withFullStop(QString text) {
    text = text.trimmed();
    if (!text.isEmpty() && !text.endsWith(QLatin1Char('.')))
        text += QLatin1Char('.');
    return text;
}

QString baseLabel(const SignalDescriptor& descriptor) {
    QString label = descriptor.label.trimmed();
    for (const QString& suffix : {QStringLiteral(" time series"),
                                  QStringLiteral(" trajectory stats"),
                                  QStringLiteral(" Welford rollup"),
                                  QStringLiteral(" statistics")}) {
        if (label.endsWith(suffix, Qt::CaseInsensitive)) {
            label.chop(suffix.size());
            break;
        }
    }
    return label;
}

QString axisName(SignalAxis axis) {
    switch (axis) {
    case SignalAxis::None:                 return QStringLiteral("item");
    case SignalAxis::Atom:                 return QStringLiteral("atom");
    case SignalAxis::Residue:              return QStringLiteral("residue");
    case SignalAxis::AtomTuple:            return QStringLiteral("atom set");
    case SignalAxis::Bond:                 return QStringLiteral("bond");
    case SignalAxis::BondVector:           return QStringLiteral("bond vector");
    case SignalAxis::Ring:                 return QStringLiteral("ring");
    case SignalAxis::AromaticRing:         return QStringLiteral("aromatic ring");
    case SignalAxis::SaturatedRing:        return QStringLiteral("saturated ring");
    case SignalAxis::RingContributionPair: return QStringLiteral("atom-ring pair");
    case SignalAxis::RingMembership:       return QStringLiteral("ring membership");
    case SignalAxis::MutationMatchPair:    return QStringLiteral("matched mutation site");
    case SignalAxis::Protein:              return QStringLiteral("protein");
    case SignalAxis::System:               return QStringLiteral("trajectory");
    case SignalAxis::Event:                return QStringLiteral("event");
    }
    return QStringLiteral("item");
}

QString meaningFor(const SignalDescriptor& descriptor) {
    const QString concept = descriptor.conceptKey;

    if (concept == QLatin1String("selections")) {
        if (descriptor.storagePath == QLatin1String("/trajectory/selections"))
            return QStringLiteral("These events record the objects and frames selected during extraction.");
        return QStringLiteral("This reports the number of atoms currently selected in Reader.");
    }

    static const struct {
        const char* key;
        const char* text;
    } meanings[] = {
        {"positions", "Coordinates locate each atom in three dimensions."},
        {"element", "This names the atom's chemical element."},
        {"topology.atoms", "The topology identifies and names each atom."},
        {"topology.residues", "The topology identifies each residue and its place in the sequence."},
        {"topology.bonds", "These records identify which atoms are covalently bonded."},
        {"topology.rings", "These records identify rings in the protein."},
        {"topology.ring_membership", "These records identify the atoms in each ring."},
        {"geometry.bond_length", "Bond length is the distance between two bonded atoms."},
        {"geometry.distance", "This measures the distance between the two selected atoms."},
        {"geometry.angle", "This measures the angle at the middle of three selected atoms."},
        {"geometry.dihedral", "This measures the signed angle between the planes formed by four selected atoms."},
        {"geometry.atom_displacement", "This is how far the atom has moved from its position in the first loaded frame."},
        {"selection_counts", "This is the number of atoms currently selected in Reader."},
        {"orca_total", "The tensor describes magnetic shielding at the atom for different directions of the applied field."},
        {"orca_diamagnetic", "The diamagnetic term contributes to the atom's shielding tensor."},
        {"orca_paramagnetic", "The paramagnetic term contributes to the atom's shielding tensor."},
        {"experimental_shielding_ml.iso", "This scalar is the orientation-averaged part of the predicted shielding tensor."},
        {"experimental_shielding_ml.t2", "These five values form the symmetric traceless part of the predicted shielding tensor. They describe its variation with orientation."},
        {"experimental_shielding_ml.t2_norm", "The norm of the five T2 values measures how much the predicted shielding varies with orientation."},
        {"bs_shielding", "This tensor describes the response to unit currents in aromatic rings at the atom."},
        {"hm_shielding", "This tensor describes the response from aromatic rings before ring-current strengths are applied."},
        {"mc_shielding", "This tensor describes the geometric response to magnetic anisotropy in surrounding bonds."},
        {"mopac_mc_shielding", "This tensor describes the geometric response to magnetic anisotropy in surrounding bonds, weighted by their MOPAC bond orders."},
        {"mc_peptide_co_rhombic", "This tensor is the change in the peptide C=O response when the axial susceptibility tensor is replaced by the rhombic tensor."},
        {"larsen_hbond_shielding", "The total sums the four Larsen hydrogen-bond tensor terms. The separate water correction is not included."},
        {"larsen_hbond_1pHB_shielding", "This primary term describes shielding changes on the donor side of a hydrogen bond donated by an amide H."},
        {"larsen_hbond_2pHB_shielding", "This secondary term describes shielding changes on the acceptor side of a hydrogen bond donated by an amide H."},
        {"larsen_hbond_1pHaB_shielding", "This primary term describes shielding changes on the donor side of a hydrogen bond donated by an alpha H."},
        {"larsen_hbond_2pHaB_shielding", "This secondary term describes shielding changes on the acceptor side of a hydrogen bond donated by an alpha H."},
        {"larsen_hbond_count", "The count includes hydrogen bonds that contribute to the Larsen calculation."},
        {"larsen_hbond_water_term", "The Larsen water correction adds isotropic shielding at an amide H when no geometric hydrogen bond is found."},
        {"bs_total_B", "This vector gives the total magnetic field from the rings included in the Biot-Savart calculation."},
        {"bs_ring_counts", "The counts group rings included in the Biot-Savart calculation by type."},
        {"pq_per_type_T0", "Each value is the scalar pi-quadrupole term for one ring type."},
        {"ringchi_per_type_T0", "These factors describe the atom's position relative to each ring type for the ring-susceptibility calculation."},
        {"disp_per_type_T0", "This measures proximity to aromatic rings by ring type, with closer ring vertices weighted more strongly."},
        {"tau_N_CA_C", "N, CA and C define this angle, with CA at the centre."},
        {"angle_N_CA_CB", "N, CA and CB define this angle, with CA at the centre."},
        {"angle_CB_CA_C", "CB, CA and C define this angle, with CA at the centre."},
        {"angle_Cprev_N_CA", "The angle at N runs from the preceding residue's carbonyl C to this residue's CA."},
        {"angle_CA_C_Nnext", "The angle at C runs from this residue's CA to the following residue's N."},
        {"cb_deviation", "This measures how far CB lies from its ideal position."},
        {"cb_residual_vector", "This vector runs from CB's ideal position to its observed position."},
        {"atom_sasa", "SASA measures the area of the atom's surface that a solvent probe can reach."},
        {"atom_sasa_fraction", "This fraction reports how much of the atom's probe-expanded surface solvent can reach."},
        {"sasa_normal", "This vector gives the mean outward direction of the atom's solvent-exposed surface."},
        {"ff_partial_charge", "The force field assigns this partial charge to the atom."},
        {"ff_pb_radius", "The Poisson-Boltzmann calculation uses this atomic radius."},
        {"eeq_charges", "Electronegativity equilibration assigns this partial charge to the atom."},
        {"eeq_cn", "This measures the atom's coordination from its neighbours and their distances."},
        {"eeq_chi_eff", "Coordination adjusts the atom's electronegativity to give this value."},
        {"eeq_hardness", "The two values are the hardness parameter and matrix diagonal used to solve the atom's charge."},
        {"apbs_phi", "The potential comes from the solvent's response to the protein's charges."},
        {"aimnet2_charges", "AIMNet2 predicts this partial charge for the atom."},
        {"aimnet2_embedding", "These numbers encode the atom's chemical environment as learned by AIMNet2."},
        {"aimnet2_charge_response_gradient", "This vector describes how the sum of squared AIMNet2 charges changes when this atom moves."},
        {"aimnet2_charge_response_gradient_scalar", "The vector length measures how strongly the sum of squared AIMNet2 charges changes when this atom moves."},
        {"aimnet2_energy_mlp", "AIMNet2 assigns this energy before applying the atom's energy offset."},
        {"aimnet2_energy_shifted_local", "AIMNet2 assigns this energy after applying the atom's energy offset."},
        {"aimnet2_d3_e_disp_atom", "This value is the atom's share of the D3 dispersion energy."},
        {"aimnet2_d3_cn", "The D3 dispersion calculation uses this coordination number for the atom."},
        {"aimnet2_d3_c6_stats", "The three values are the sum, mean and maximum of the atom's D3 C6 coefficients across its neighbours."},
        {"pyramidalization", "This measures how far an sp2 atom lies out of its local plane."},
        {"omega_actual", "Omega measures rotation about the peptide C-N bond."},
        {"omega_deviation", "The signed difference from 180 degrees measures omega's departure from trans."},
        {"aromatic_chi2", "Chi2 describes rotation of the aromatic side chain about its second torsion bond."},
        {"aromatic_ring_chi2", "Chi2 describes rotation of the aromatic side chain about its second torsion bond."},
        {"pucker_Q", "The Cremer-Pople amplitude measures how far the five-membered ring departs from a plane."},
        {"saturated_ring_pucker_amplitude", "The puckering amplitude measures how far the ring departs from a plane."},
        {"pucker_theta", "The Cremer-Pople phase describes the pattern of puckering around the five-membered ring."},
        {"saturated_ring_pucker_phase", "The phase describes the pattern of puckering around the ring."},
        {"mopac_charges_full_precision", "MOPAC assigns this Coulson partial charge to the atom."},
        {"mopac_bond_valencies_full_precision", "The diagonal of MOPAC's Wiberg bond-order matrix gives this atomic valency."},
        {"mopac_atom_s_population", "MOPAC assigns this electron population to the atom's s orbital."},
        {"mopac_atom_p_population", "MOPAC assigns this total electron population to the atom's p orbitals."},
        {"mopac_atom_d_population", "MOPAC assigns this total electron population to the atom's d orbitals."},
        {"mopac_lewis_bond_count", "The count gives the atom's bonds in MOZYME's Lewis structure."},
        {"mopac_mc_co_sum", "This is the sum of scalar McConnell terms from C=O bonds, each weighted by its MOPAC bond order."},
        {"mopac_mc_cn_sum", "This is the sum of scalar McConnell terms from C-N bonds, each weighted by its MOPAC bond order."},
        {"mopac_mc_sidechain_sum", "This is the sum of scalar McConnell terms from side chains, each weighted by its MOPAC bond order."},
        {"mopac_mc_aromatic_sum", "This is the sum of scalar McConnell terms from aromatic groups, each weighted by its MOPAC bond order."},
        {"mopac_mc_nearest_co_dist", "This is the distance from the atom to the midpoint of the nearest accepted C=O source."},
        {"mopac_mc_nearest_cn_dist", "This is the distance from the atom to the midpoint of the nearest accepted C-N source."},
        {"hbond_nearest_dir", "This unit vector points from the nearest accepted hydrogen-bond donor H to the atom."},
        {"mc_nearest_co_dir", "This unit vector points from the midpoint of the nearest accepted peptide C=O bond to the atom."},
        {"mc_nearfield_counts", "These counts distinguish included and rejected McConnell sources less than 3 angstroms away."},
        {"hbond_scalars", "These values give the distance to the nearest accepted donor H, its inverse cube, the number of nearby hydrogen bonds and the sum of their angular terms."},
        {"dssp_ss8", "DSSP assigns the residue to one of eight secondary-structure classes."},
        {"dssp_observed", "This indicates whether DSSP returned a result for the atom's residue."},
        {"dssp_ppii", "This indicates whether DSSP assigned a polyproline-II conformation."},
        {"dssp_hbond_energy", "These are DSSP's electrostatic scores for the residue's hydrogen bonds as donor and acceptor."},
        {"dssp_torsion_angle", "The phi, psi and four chi angles describe the residue's backbone and side-chain torsions."},
        {"water_shell_counts", "The counts report water molecules in the first and second shells around the atom."},
        {"hydration_shell", "These values describe the distribution and orientation of nearby water, and the distance and charge of the nearest ion."},
        {"hydration_geometry", "These values describe the number, net dipole and orientation of first-shell water molecules relative to the exposed surface."},
        {"residue_dihedral", "The phi, psi, omega and chi angles describe rotation about the residue's backbone and side-chain bonds."},
        {"j_coupling", "The residue's torsion angles provide estimates of its scalar spin-spin couplings."},
        {"ring_neighbourhood", "The distances and angles locate the atom relative to nearby rings."},
        {"bonded_energy", "Bonded force-field terms assign this energy to the atom."},
        {"gromacs_energy", "These records report the simulation's energy, temperature, pressure and volume."},
        {"rmsd_tracking", "RMSD measures how far the fitted structure differs from its reference."},
        {"bs_shielding.autocorrelation", "This measures how the Biot-Savart T0 value correlates with itself at later times."},
        {"kernel_dynamics.acf", "Each curve measures how a calculated signal correlates with itself at later times."},
        {"kernel_dynamics.psd", "Each spectrum shows how the power in a signal is distributed over frequency."},
        {"kernel_dynamics.decay_time", "This estimates how long correlations persist in each signal."},
        {"kernel_dynamics.peak_freq", "The peak is the nonzero frequency with the greatest power in each signal's spectrum."},
        {"kernel_dynamics.spectral_centroid", "The centroid is the mean frequency weighted by spectral power."},
        {"kernel_coherence.matrix", "Each matrix entry is the Pearson correlation between two of the seven calculated signals."},
        {"dihedral.phi_corr_time", "This estimates how long correlations persist in the residue's phi angle."},
        {"dihedral.psi_corr_time", "This estimates how long correlations persist in the residue's psi angle."},
        {"dihedral.chi_corr_time", "These times estimate how long correlations persist in the residue's chi angles."},
        {"dihedral.phi_acf", "This curve measures how the phi angle correlates with itself at later times."},
        {"dihedral.psi_acf", "This curve measures how the psi angle correlates with itself at later times."},
        {"dihedral.chi_acf", "These curves measure how each chi angle correlates with itself at later times."},
        {"reorient.s2", "S2 describes how strongly a bond's orientation is restricted."},
        {"reorient.tau_e", "This estimates the timescale of the bond's internal reorientation."},
        {"reorient.r1", "R1 estimates the rate at which nitrogen-15 magnetisation recovers along the applied field."},
        {"reorient.r2", "R2 estimates the rate at which nitrogen-15 magnetisation decays perpendicular to the applied field."},
        {"reorient.noe", "This estimates the nitrogen-15 NOE produced by irradiation of its bonded proton."},
        {"reorient.acf_internal", "This measures how a bond's orientation correlates with itself at later times after removing the protein's tumbling."},
        {"reorient.acf_lab", "This measures how a bond's orientation correlates with itself at later times, including the protein's tumbling."},
        {"reorient.orientation_tensor", "This tensor describes the directions a bond occupies after the protein's overall rotation has been removed."},
        {"reorient.spectral_density", "These values describe motion at five frequency combinations used to calculate N-H relaxation."},
        {"ired.s2", "This order parameter describes how strongly an N-H bond's orientation is restricted, using iRED."},
        {"dssp_ss8.transition", "This is the number of DSSP class changes between successive observed frames."},
        {"residue_dihedral.transition", "This is the number of backbone dihedral-bin changes between successive observed frames."},
    };
    for (const auto& entry : meanings) {
        if (concept == QLatin1String(entry.key))
            return QString::fromLatin1(entry.text);
    }

    if (concept.endsWith(QLatin1String(".stats"))) {
        QString quantity = baseLabel(descriptor);
        if (concept == QLatin1String("bond_length.stats"))
            quantity = QStringLiteral("bond length");
        else if (concept == QLatin1String("bs_shielding.stats"))
            quantity = QStringLiteral("the Biot-Savart response to unit ring currents");
        else if (concept == QLatin1String("hm_shielding.stats"))
            quantity = QStringLiteral("the Haigh-Mallion response before ring-current strengths are applied");
        else if (concept == QLatin1String("mc_shielding.stats"))
            quantity = QStringLiteral("the McConnell geometric response");
        else if (concept == QLatin1String("atom_sasa.stats"))
            quantity = QStringLiteral("solvent-accessible area");
        else if (concept == QLatin1String("eeq_charges.stats"))
            quantity = QStringLiteral("the EEQ charge");
        else if (concept == QLatin1String("hbond_count.stats"))
            quantity = QStringLiteral("the hydrogen-bond count");
        else if (concept == QLatin1String("mopac_charges.stats"))
            quantity = QStringLiteral("the MOPAC charge");
        else if (concept == QLatin1String("mopac_bond_orders.stats"))
            quantity = QStringLiteral("the MOPAC bond order");
        else if (concept == QLatin1String("water_field.stats"))
            quantity = QStringLiteral("the water electric field, EFG and shell counts");
        else if (concept == QLatin1String("aimnet2_charge_response_gradient.stats"))
            quantity = QStringLiteral("the AIMNet2 charge-response gradient");
        else if (concept == QLatin1String("hydration_shell.stats"))
            quantity = QStringLiteral("the water and ion measurements around the atom");
        else if (concept == QLatin1String("hydration_geometry.stats"))
            quantity = QStringLiteral("the count, net dipole and orientation of first-shell waters");
        return QStringLiteral("Mean and variance summarise %1 across the trajectory.").arg(quantity);
    }

    if (concept.startsWith(QLatin1String("bs_per_type_"))
        || concept.startsWith(QLatin1String("hm_per_type_"))) {
        const QString method = concept.startsWith(QLatin1String("bs_"))
                                   ? QStringLiteral("Biot-Savart") : QStringLiteral("Haigh-Mallion");
        if (concept.endsWith(QLatin1String("T0")))
            return QStringLiteral("T0 is the scalar part of the %1 response for each ring type.").arg(method);
        if (concept.endsWith(QLatin1String("T1")))
            return QStringLiteral("The three T1 values form the antisymmetric part of the %1 response for each ring type.").arg(method);
        return QStringLiteral("The five T2 values form the symmetric traceless part of the %1 response for each ring type.").arg(method);
    }

    if (concept.contains(QLatin1String("_efg")) || concept.contains(QLatin1String("_efield"))
        || concept.endsWith(QLatin1String("_E")) || concept.contains(QLatin1String("_E_"))) {
        QString source;
        if (descriptor.family == QLatin1String("apbs"))
            source = QStringLiteral("the solvent's response to the protein's charges");
        else if (descriptor.family == QLatin1String("water_field"))
            source = concept.endsWith(QLatin1String("_first"))
                         ? QStringLiteral("water in the first shell") : QStringLiteral("surrounding water");
        else if (concept.endsWith(QLatin1String("_backbone")))
            source = QStringLiteral("backbone charges");
        else if (concept.endsWith(QLatin1String("_sidechain")))
            source = QStringLiteral("non-aromatic side-chain charges");
        else if (concept.endsWith(QLatin1String("_aromatic")))
            source = QStringLiteral("aromatic-group charges");
        else
            source = QStringLiteral("surrounding atomic charges");
        if (concept.contains(QLatin1String("_efg")))
            return QStringLiteral("This traceless tensor describes spatial variation in the electric field from %1. The reported EFG is the potential Hessian, the negative of the field gradient.").arg(source);
        return QStringLiteral("This vector gives the electric field from %1.").arg(source);
    }

    if (concept.endsWith(QLatin1String("coulomb_scalars")))
        return QStringLiteral("The four values report the field magnitude, its projection along the primary bond, the backbone-to-total magnitude ratio and the aromatic field magnitude.");

    if (descriptor.family == QLatin1String("mcconnell")
        || descriptor.family == QLatin1String("mopac_mcconnell")) {
        QString bonds;
        if (concept.contains(QLatin1String("nearest_co")))
            bonds = QStringLiteral("the nearest included C=O bond");
        else if (concept.contains(QLatin1String("nearest_cn")))
            bonds = QStringLiteral("the nearest included C-N bond");
        else if (concept.contains(QLatin1String("peptide_co")))
            bonds = QStringLiteral("peptide C=O bonds");
        else if (concept.contains(QLatin1String("peptide_cn")))
            bonds = QStringLiteral("peptide C-N bonds");
        else if (concept.contains(QLatin1String("backbone_other")))
            bonds = QStringLiteral("other backbone bonds");
        else if (concept.contains(QLatin1String("sidechain_co")))
            bonds = QStringLiteral("side-chain C=O bonds");
        else if (concept.contains(QLatin1String("sidechain_other")))
            bonds = QStringLiteral("other side-chain bonds");
        else if (concept.contains(QLatin1String("disulfide")))
            bonds = QStringLiteral("disulfide bonds");
        else if (concept.contains(QLatin1String("aromatic")))
            bonds = QStringLiteral("aromatic groups");
        else if (concept.contains(QLatin1String("backbone_xh")))
            bonds = QStringLiteral("backbone bonds to hydrogen");
        else if (concept.contains(QLatin1String("sidechain_xh")))
            bonds = QStringLiteral("side-chain bonds to hydrogen");
        else if (concept.contains(QLatin1String("s_h")))
            bonds = QStringLiteral("S-H bonds");
        else if (concept.contains(QLatin1String("backbone")))
            bonds = QStringLiteral("backbone bonds");
        else if (concept.contains(QLatin1String("sidechain")))
            bonds = QStringLiteral("side-chain bonds");
        if (!bonds.isEmpty())
            return QStringLiteral("This tensor describes the geometric response to magnetic anisotropy in %1.").arg(bonds);
    }

    return withFullStop(QStringLiteral("%1 for the selected %2")
                            .arg(baseLabel(descriptor), axisName(descriptor.nativeAxis)));
}

QString formSentence(const SignalDescriptor& descriptor) {
    const MetricClass metricClass = ClassifyMetric(descriptor);
    switch (metricClass.form) {
    case MetricForm::Snapshot:
    case MetricForm::Series:
        return {};
    case MetricForm::Rollup:
        return QStringLiteral("Welford's method accumulates the mean and variance across frames.");
    case MetricForm::Dynamics:
    case MetricForm::Transition:
    case MetricForm::Reference:
    case MetricForm::Spine:
    case MetricForm::Derived:
    case MetricForm::Other:
        return {};
    }
    return {};
}

struct FamilyText {
    QString calculation;
    QString origin;
};

std::optional<FamilyText> familyText(const SignalDescriptor& descriptor) {
    const QString& family = descriptor.family;
    const QString& concept = descriptor.conceptKey;

    if (family == QLatin1String("biot_savart")) {
        const QString origin = QStringLiteral("Biot-Savart calculator in nmr_extract.");
        if (concept == QLatin1String("bs_ring_counts"))
            return FamilyText{QStringLiteral("The calculation counts accepted atom-ring pairs separately for each ring type."), origin};
        if (concept == QLatin1String("bs_total_B"))
            return FamilyText{QStringLiteral("The calculation sums magnetic-field vectors from the unit-current double-loop model of each accepted aromatic ring."), origin};
        if (concept == QLatin1String("bs_shielding.autocorrelation"))
            return FamilyText{QStringLiteral("The calculation subtracts the mean Biot-Savart T0, correlates the fluctuations at each lag and divides by the zero-lag covariance."), origin};
        if (concept.startsWith(QLatin1String("bs_per_type_")))
            return FamilyText{QStringLiteral("The calculation forms a response tensor from each ring's unit-current field and normal, separates its T0, T1 and T2 parts, and sums the named part by ring type."), origin};
        return FamilyText{QStringLiteral("The calculation forms a response tensor from each ring's unit-current double-loop field and normal, then sums the accepted rings."), origin};
    }
    if (family == QLatin1String("haigh_mallion")) {
        const QString origin = QStringLiteral("Haigh-Mallion calculator in nmr_extract.");
        if (concept.startsWith(QLatin1String("hm_per_type_")))
            return FamilyText{QStringLiteral("The calculation integrates over each ring's triangulated surface, forms a response tensor using the ring normal, and sums the named tensor part by ring type."), origin};
        return FamilyText{QStringLiteral("The calculation integrates over each ring's triangulated surface and uses the ring normal to form its response tensor. It then sums the accepted ring contributions."), origin};
    }
    if (family == QLatin1String("pi_quadrupole"))
        return FamilyText{QStringLiteral("The calculation sums (3 cos^2(theta) - 1) / r^4 by ring type. Here r is the distance from the ring centre to the atom, and theta is the angle between that displacement and the ring normal."),
                          QStringLiteral("Pi-quadrupole calculator in nmr_extract.")};
    if (family == QLatin1String("ring_susceptibility"))
        return FamilyText{QStringLiteral("The calculation sums (3 cos^2(theta) - 1) / r^3 by ring type. Here r is the distance from the ring centre to the atom, and theta is the angle between that displacement and the ring normal."),
                          QStringLiteral("Ring-susceptibility calculator in nmr_extract.")};
    if (family == QLatin1String("dispersion"))
        return FamilyText{QStringLiteral("The calculation sums 1 / r^6 over accepted ring vertices, with r measured from the atom to each vertex. The configured switching function modifies each term."),
                          QStringLiteral("Dispersion calculator in nmr_extract.")};
    if (family == QLatin1String("ring_current"))
        return FamilyText{QStringLiteral("Ring centres, normals, radii and atom-ring positions are calculated from the conformation."),
                          QStringLiteral("Ring and atom records from nmr_extract.")};
    if (family == QLatin1String("identity"))
        return FamilyText{QStringLiteral("Reader loads the Cartesian coordinates stored for each atom and frame."),
                          QStringLiteral("Atomic coordinates from nmr_extract.")};

    if (family == QLatin1String("mcconnell")) {
        if (concept == QLatin1String("mc_nearest_co_dir"))
            return FamilyText{QStringLiteral("The calculation subtracts the nearest accepted peptide C=O midpoint from the atom's position and normalises the resulting vector."),
                              QStringLiteral("McConnell calculator in nmr_extract.")};
        if (concept == QLatin1String("mc_nearfield_counts"))
            return FamilyText{QStringLiteral("The calculation counts accepted and rejected McConnell sources within 3 angstroms of the atom."),
                              QStringLiteral("McConnell calculator in nmr_extract.")};
        if (concept == QLatin1String("mc_peptide_co_rhombic"))
            return FamilyText{QStringLiteral("The calculation subtracts the axial response from the rhombic response for each accepted peptide C=O bond and sums the differences."),
                              QStringLiteral("McConnell calculator in nmr_extract.")};
        if (concept.contains(QLatin1String("nearest_")))
            return FamilyText{QStringLiteral("The calculation selects the nearest accepted peptide bond of this type and combines its axial susceptibility tensor with the 1 / r^3 dipolar tensor."),
                              QStringLiteral("McConnell calculator in nmr_extract.")};
        return FamilyText{QStringLiteral("The calculation combines each accepted bond's unit susceptibility tensor with the 1 / r^3 dipolar tensor and sums the responses."),
                          QStringLiteral("McConnell calculator in nmr_extract.")};
    }
    if (family == QLatin1String("mopac_mcconnell")) {
        if (concept.contains(QLatin1String("nearest_co_dist"))
            || concept.contains(QLatin1String("nearest_cn_dist"))) {
            return FamilyText{QStringLiteral("The calculation measures the distance from the atom to the midpoint of the nearest accepted peptide bond of this type."),
                              QStringLiteral("MOPAC-weighted McConnell calculator in nmr_extract.")};
        }
        if (concept.endsWith(QLatin1String("_sum"))) {
            return FamilyText{QStringLiteral("The calculation weights each scalar contribution from the named source type by its MOPAC bond order, then sums the contributions."),
                              QStringLiteral("MOPAC-weighted McConnell calculator in nmr_extract.")};
        }
        if (concept.contains(QLatin1String("nearest_")))
            return FamilyText{QStringLiteral("The calculation selects the nearest accepted peptide bond of this type and multiplies its axial McConnell response by its MOPAC bond order."),
                              QStringLiteral("MOPAC-weighted McConnell calculator in nmr_extract.")};
        return FamilyText{QStringLiteral("The calculation combines each bond's unit susceptibility tensor with the 1 / r^3 dipolar tensor, weights the response by its MOPAC bond order and sums the accepted bonds."),
                          QStringLiteral("MOPAC-weighted McConnell calculator in nmr_extract.")};
    }

    if (family == QLatin1String("coulomb") || family == QLatin1String("eeq_coulomb")
        || family == QLatin1String("mopac_coulomb")
        || (family == QLatin1String("aimnet2")
            && (concept.startsWith(QLatin1String("aimnet2_E"))
                || concept.startsWith(QLatin1String("aimnet2_efg"))))) {
        const QString chargeSource = family == QLatin1String("coulomb")
                                         ? QStringLiteral("force-field")
                                         : (family == QLatin1String("eeq_coulomb")
                                                ? QStringLiteral("EEQ")
                                                : (family == QLatin1String("aimnet2")
                                                       ? QStringLiteral("AIMNet2") : QStringLiteral("MOPAC")));
        const QString origin = QStringLiteral("Electrostatics calculator in nmr_extract, using %1 charges.").arg(chargeSource);
        if (concept.endsWith(QLatin1String("coulomb_scalars")))
            return FamilyText{QStringLiteral("The calculation takes vector lengths and a bond-direction dot product from the charge-derived fields. The backbone ratio is |E_backbone| / |E_total|."), origin};
        QString sources = QStringLiteral("accepted atoms");
        if (concept.endsWith(QLatin1String("_aromatic")))
            sources = QStringLiteral("accepted aromatic atoms");
        else if (concept.endsWith(QLatin1String("_backbone")))
            sources = QStringLiteral("accepted backbone atoms");
        else if (concept.endsWith(QLatin1String("_sidechain")))
            sources = QStringLiteral("accepted non-aromatic side-chain atoms");
        if (concept.contains(QLatin1String("_efg")))
            return FamilyText{QStringLiteral("The calculation sums second spatial derivatives of the Coulomb potential from %1, using %2 charges, and removes the tensor trace.").arg(sources, chargeSource), origin};
        return FamilyText{QStringLiteral("The calculation sums Coulomb electric-field vectors from %1, using their positions and %2 charges.").arg(sources, chargeSource), origin};
    }
    if (family == QLatin1String("apbs")) {
        const QString origin = QStringLiteral("APBS continuum electrostatics, processed by nmr_extract.");
        if (concept == QLatin1String("apbs_phi"))
            return FamilyText{QStringLiteral("The calculation subtracts the homogeneous-vacuum reference potential from the solvated APBS potential and interpolates the difference at the atom."), origin};
        if (concept == QLatin1String("apbs_E"))
            return FamilyText{QStringLiteral("The calculation takes the negative gradient of the APBS reaction potential at the atom, using central differences on the grid."), origin};
        return FamilyText{QStringLiteral("The calculation takes second spatial derivatives of the APBS reaction potential at the atom using grid differences, symmetrises the tensor and removes its trace."), origin};
    }

    if (family == QLatin1String("larsen_hbond")) {
        if (concept == QLatin1String("larsen_hbond_water_term"))
            return FamilyText{QStringLiteral("The calculation assigns 2.07 ppm to an eligible amide H when no geometric hydrogen bond is found."),
                              QStringLiteral("Larsen hydrogen-bond calculator in nmr_extract.")};
        if (concept == QLatin1String("larsen_hbond_count"))
            return FamilyText{QStringLiteral("The calculation counts hydrogen bonds that contributed in any Larsen class."),
                              QStringLiteral("Larsen hydrogen-bond calculator in nmr_extract.")};
        return FamilyText{QStringLiteral("For each accepted hydrogen bond, the calculation uses its distance and angles to query the Larsen tensor grid, rotates the tensor into the conformation frame and sums it at the atom."),
                          QStringLiteral("Larsen hydrogen-bond calculator in nmr_extract.")};
    }
    if (family == QLatin1String("hbond")) {
        const QString origin = QStringLiteral("Hydrogen-bond calculator in nmr_extract, using DSSP backbone hydrogen bonds.");
        if (concept == QLatin1String("hbond_nearest_dir"))
            return FamilyText{QStringLiteral("The calculation selects the nearest accepted donor H, subtracts its position from the atom's position and normalises the vector."), origin};
        if (concept == QLatin1String("hbond_count.stats"))
            return FamilyText{QStringLiteral("The calculation counts accepted hydrogen bonds whose donor H lies within the configured counting radius of the atom."), origin};
        return FamilyText{QStringLiteral("Distances run from donor H to the atom. The angular sum uses (3 cos^2(theta) - 1) / r^3, where theta is the angle to the H-bond axis; nearby sources are counted within the configured radius."), origin};
    }

    if (family == QLatin1String("aimnet2")) {
        QString calculation;
        if (concept == QLatin1String("aimnet2_embedding")) {
            calculation = QStringLiteral("AIMNet2 produces this per-atom vector while evaluating the conformation.");
        } else if (concept == QLatin1String("aimnet2_charge_response_gradient_scalar")) {
            calculation = QStringLiteral("Automatic differentiation gives the gradient of the sum of squared AIMNet2 charges with respect to the atom's three coordinates. The reported scalar is its Euclidean length.");
        } else if (concept.startsWith(QLatin1String("aimnet2_charge_response_gradient"))) {
            calculation = QStringLiteral("Automatic differentiation gives the gradient of the sum of squared AIMNet2 charges with respect to the atom's three position coordinates.");
        } else if (concept == QLatin1String("aimnet2_energy_mlp")) {
            calculation = QStringLiteral("AIMNet2 produces a local neural-network energy for each atom before adding the atomic energy shift.");
        } else if (concept == QLatin1String("aimnet2_energy_shifted_local")) {
            calculation = QStringLiteral("AIMNet2 adds the atomic energy shift to the local neural-network energy before summing over atoms.");
        } else if (concept == QLatin1String("aimnet2_d3_c6_stats")) {
            calculation = QStringLiteral("The calculation takes the sum, mean and maximum of the D3 module's pairwise C6 coefficients over the atom's valid neighbours.");
        } else if (concept == QLatin1String("aimnet2_d3_cn")) {
            calculation = QStringLiteral("The D3 calculation sums smooth distance-dependent neighbour weights based on covalent radii, then caps the sum at the element's maximum coordination number.");
        } else if (concept == QLatin1String("aimnet2_d3_e_disp_atom")) {
            calculation = QStringLiteral("The calculation sums the D3 module's damped C6 / r^6 and C8 / r^8 pair terms for the atom and converts the energy to electronvolts.");
        } else {
            calculation = QStringLiteral("AIMNet2 predicts atomic partial charges from the elements and positions in the conformation.");
        }
        return FamilyText{calculation, QStringLiteral("AIMNet2, run by nmr_extract.")};
    }
    if (family == QLatin1String("eeq")) {
        const QString origin = QStringLiteral("The charge-equilibration model implemented in nmr_extract.");
        if (concept == QLatin1String("eeq_cn"))
            return FamilyText{QStringLiteral("The calculation sums smooth neighbour weights based on each separation relative to the sum of the two covalent radii."), origin};
        if (concept == QLatin1String("eeq_chi_eff"))
            return FamilyText{QStringLiteral("The calculation adds an element-specific coefficient times the square root of the coordination number to the element's base electronegativity."), origin};
        if (concept == QLatin1String("eeq_hardness"))
            return FamilyText{QStringLiteral("The first value is the element's hardness parameter. The second adds the Gaussian charge self-interaction term to form the charge-equation diagonal."), origin};
        return FamilyText{QStringLiteral("The calculation solves coupled charge equations using coordination-adjusted electronegativities, atomic hardnesses and interatomic distances, subject to the total molecular charge."), origin};
    }
    if (family == QLatin1String("mopac") || family == QLatin1String("mopac_core")) {
        const QString origin = QStringLiteral("MOPAC PM7 with MOZYME and 1SCF, run by nmr_extract.");
        if (concept.contains(QLatin1String("population")))
            return FamilyText{QStringLiteral("The calculation sums diagonal entries of MOPAC's atomic density matrix over the named orbital shell."), origin};
        if (concept == QLatin1String("mopac_bond_valencies_full_precision"))
            return FamilyText{QStringLiteral("nmr_extract reads the atomic valency from the diagonal of MOPAC's Wiberg bond-order matrix."), origin};
        if (concept == QLatin1String("mopac_lewis_bond_count"))
            return FamilyText{QStringLiteral("nmr_extract reads the number of bonds assigned to the atom in MOZYME's Lewis structure."), origin};
        if (concept == QLatin1String("mopac_bond_orders.stats"))
            return FamilyText{QStringLiteral("MOPAC calculates Wiberg bond orders for the conformation."), origin};
        return FamilyText{QStringLiteral("MOPAC's population analysis assigns Coulson partial charges for the conformation."), origin};
    }
    if (family == QLatin1String("force_field") && concept == QLatin1String("ff_pb_radius"))
        return FamilyText{QStringLiteral("nmr_extract uses the supplied Poisson-Boltzmann radius or derives it with the mbondi2 radius rules."),
                          QStringLiteral("Prepared force-field data and radius assignment in nmr_extract.")};
    if (family == QLatin1String("force_field"))
        return FamilyText{QStringLiteral("The force-field topology assigns this value to the atom."),
                          QStringLiteral("Force-field topology read by nmr_extract.")};

    if (family == QLatin1String("sasa")) {
        const QString origin = QStringLiteral("SASA calculator in nmr_extract.");
        if (concept == QLatin1String("sasa_normal"))
            return FamilyText{QStringLiteral("The calculation sums directions to unoccluded sample points on the probe-expanded atomic sphere and normalises the sum. A zero sum gives a zero vector."), origin};
        if (concept == QLatin1String("atom_sasa_fraction"))
            return FamilyText{QStringLiteral("The calculation divides the accessible area by 4 pi times the square of the sum of the atom's Bondi radius and the solvent-probe radius."), origin};
        return FamilyText{QStringLiteral("The calculation samples points on the probe-expanded atomic sphere. The unoccluded fraction times the sphere's area gives SASA."), origin};
    }
    if (family == QLatin1String("hydration")) {
        const QString origin = QStringLiteral("Hydration calculator in nmr_extract.");
        if (concept.startsWith(QLatin1String("hydration_geometry")))
            return FamilyText{QStringLiteral("The calculation sums first-shell water dipoles and compares their direction and positions with the SASA surface normal."), origin};
        return FamilyText{QStringLiteral("The calculation measures the fraction of waters on the side away from the protein centre, averages dipole alignment with atom-to-water directions and finds the nearest ion within the cutoff."), origin};
    }
    if (family == QLatin1String("water_field")) {
        const QString origin = QStringLiteral("Water-field calculator in nmr_extract.");
        if (concept == QLatin1String("water_shell_counts"))
            return FamilyText{QStringLiteral("The calculation counts water oxygens in the first shell and in the region between the first and second shell boundaries, using the configured distances."), origin};
        if (concept == QLatin1String("water_field.stats"))
            return FamilyText{QStringLiteral("The calculation sums Coulomb fields and potential Hessians from explicit water charges and counts water oxygens in each shell."), origin};
        const QString waters = concept.endsWith(QLatin1String("_first"))
                                   ? QStringLiteral("water molecules whose oxygens lie in the first shell")
                                   : QStringLiteral("water molecules whose oxygens lie within the field cutoff");
        if (concept.contains(QLatin1String("_efg")))
            return FamilyText{QStringLiteral("The calculation sums Coulomb potential Hessians from the charge sites of %1 and removes the tensor trace.").arg(waters), origin};
        return FamilyText{QStringLiteral("The calculation sums Coulomb electric-field vectors from the charge sites of %1.").arg(waters), origin};
    }
    if (family == QLatin1String("water_polarization"))
        return FamilyText{QStringLiteral("Charges and dipoles represent the water around the atom."),
                          QStringLiteral("Water-polarization calculator in nmr_extract.")};

    if (family == QLatin1String("dssp")) {
        if (concept == QLatin1String("dssp_torsion_angle"))
            return FamilyText{QStringLiteral("nmr_extract stores phi and psi from DSSP and calculates the four chi angles from their defining atom coordinates."),
                              QStringLiteral("DSSP and the torsion calculation in nmr_extract.")};
        const QString origin = QStringLiteral("DSSP, run by nmr_extract.");
        if (concept == QLatin1String("dssp_hbond_energy"))
            return FamilyText{QStringLiteral("DSSP evaluates backbone hydrogen bonds with its electrostatic distance formula. nmr_extract records two acceptor and two donor scores for each residue."), origin};
        if (concept == QLatin1String("dssp_observed"))
            return FamilyText{QStringLiteral("nmr_extract marks an atom as observed when its parent residue has a DSSP result."), origin};
        if (concept == QLatin1String("dssp_ppii"))
            return FamilyText{QStringLiteral("nmr_extract tests whether the residue's DSSP assignment is polyproline II and records the result for each atom in that residue."), origin};
        if (concept == QLatin1String("dssp_ss8.transition"))
            return FamilyText{QStringLiteral("nmr_extract compares successive observed DSSP assignments and counts transitions between different classes for each residue."), origin};
        if (descriptor.id == QLatin1String("npy:dssp_ss8"))
            return FamilyText{QStringLiteral("DSSP assigns secondary structure from backbone geometry and hydrogen bonds. nmr_extract marks the assigned class in an eight-column array for each atom."), origin};
        return FamilyText{QStringLiteral("DSSP assigns secondary structure from backbone geometry and hydrogen bonds. nmr_extract stores the class for each residue and frame."), origin};
    }
    if (family == QLatin1String("local_geometry")) {
        const QString origin = QStringLiteral("Local backbone geometry calculator in nmr_extract.");
        if (concept == QLatin1String("cb_deviation"))
            return FamilyText{QStringLiteral("The calculation constructs an ideal CB position from N, CA and C using fixed coefficients, then measures its distance from the observed CB."), origin};
        if (concept == QLatin1String("cb_residual_vector"))
            return FamilyText{QStringLiteral("The calculation constructs an ideal CB position from N, CA and C using fixed coefficients, then subtracts it from the observed CB position."), origin};
        return FamilyText{QStringLiteral("The calculation forms two vectors from the central atom to the named endpoints and takes the arccosine of their normalised dot product."), origin};
    }
    if (family == QLatin1String("planar_geometry")) {
        const QString origin = QStringLiteral("Geometry calculator in nmr_extract.");
        if (concept == QLatin1String("pyramidalization"))
            return FamilyText{QStringLiteral("The calculation measures the absolute perpendicular distance from the atom to the plane through its three bonded neighbours."), origin};
        if (concept == QLatin1String("omega_deviation"))
            return FamilyText{QStringLiteral("The calculation subtracts 180 degrees from the peptide omega torsion and wraps the difference to the signed half-circle range."), origin};
        if (concept == QLatin1String("omega_actual"))
            return FamilyText{QStringLiteral("The calculation measures the signed CA-C-N-CA dihedral across the peptide bond."), origin};
        if (concept.contains(QLatin1String("pucker"))) {
            if (concept.contains(QLatin1String("phase")) || concept == QLatin1String("pucker_theta"))
                return FamilyText{QStringLiteral("The calculation projects five ring atoms' displacements from the Cremer-Pople mean plane onto sine and cosine modes, then takes the phase of the two coefficients."), origin};
            return FamilyText{QStringLiteral("The calculation projects five ring atoms' displacements from the Cremer-Pople mean plane onto sine and cosine modes. The amplitude is the square root of the sum of the squared coefficients."), origin};
        }
        if (concept == QLatin1String("residue_dihedral.transition"))
            return FamilyText{QStringLiteral("nmr_extract assigns backbone phi and psi to Ramachandran bins and counts transitions between different bins in successive observed frames."), origin};
        return FamilyText{QStringLiteral("Each torsion is the signed angle between the two planes defined by its four ordered atoms, measured about the central bond."), origin};
    }
    if (family == QLatin1String("geometry")) {
        const QString origin = QStringLiteral("Reader geometry tools.");
        if (concept == QLatin1String("geometry.atom_displacement"))
            return FamilyText{QStringLiteral("Reader subtracts the atom's position in the first loaded frame from its current position and plots the resulting vector's length."), origin};
        if (concept == QLatin1String("geometry.distance"))
            return FamilyText{QStringLiteral("Reader takes the Euclidean distance between the two selected atom coordinates."), origin};
        if (concept == QLatin1String("geometry.angle"))
            return FamilyText{QStringLiteral("Reader measures the angle between vectors from the middle selected atom to the other two atoms."), origin};
        return FamilyText{QStringLiteral("Reader measures the signed angle between the first and last three-atom planes, about the bond between the middle two selected atoms."), origin};
    }
    if (family == QLatin1String("j_coupling"))
        return FamilyText{QStringLiteral("The coupling is calculated from the relevant torsion angle using the configured Karplus relation."),
                          QStringLiteral("J-coupling calculator in nmr_extract.")};

    if (family == QLatin1String("kernel_dynamics")) {
        const QString origin = QStringLiteral("Kernel dynamics calculator in nmr_extract.");
        if (concept == QLatin1String("kernel_dynamics.acf"))
            return FamilyText{QStringLiteral("The calculation correlates each signal's fluctuations about its mean at increasing lags and divides by the zero-lag covariance."), origin};
        if (concept == QLatin1String("kernel_dynamics.decay_time"))
            return FamilyText{QStringLiteral("The calculation sums the normalised autocorrelation up to its first nonpositive value and multiplies by the sampling interval. If there is no crossing, the sum uses the full lag window."), origin};
        if (concept == QLatin1String("kernel_dynamics.peak_freq"))
            return FamilyText{QStringLiteral("The calculation selects the strongest nonzero-frequency bin in the Parzen-windowed power spectrum."), origin};
        if (concept == QLatin1String("kernel_dynamics.spectral_centroid"))
            return FamilyText{QStringLiteral("The calculation sums frequency times spectral power over all bins of the Parzen-windowed spectrum, then divides by the total power."), origin};
        return FamilyText{QStringLiteral("The calculation applies a Parzen window to the signal's autocovariance and transforms it into a power spectrum."), origin};
    }
    if (family == QLatin1String("kernel_coherence"))
        return FamilyText{QStringLiteral("Pearson correlation compares each pair of signals across the trajectory."),
                          QStringLiteral("Kernel correlation calculator in nmr_extract.")};
    if (family == QLatin1String("dihedral_autocorrelation")) {
        const QString origin = QStringLiteral("Dihedral dynamics calculator in nmr_extract.");
        if (concept.endsWith(QLatin1String("_corr_time")))
            return FamilyText{QStringLiteral("The calculation finds the first lag where the mean cosine of the angle change falls to 1/e and interpolates the crossing time. If there is no crossing, it reports the full lag window."), origin};
        return FamilyText{QStringLiteral("At each time lag, the calculation averages the cosine of the difference between the two torsion angles over all available frame pairs."), origin};
    }
    if (family == QLatin1String("reorientational_dynamics")) {
        const QString origin = QStringLiteral("Reorientational dynamics calculator in nmr_extract.");
        if (concept == QLatin1String("reorient.orientation_tensor"))
            return FamilyText{QStringLiteral("The calculation aligns each frame to the reference and averages u u^T, where u is the unit bond vector in the aligned frame."), origin};
        if (concept == QLatin1String("reorient.spectral_density"))
            return FamilyText{QStringLiteral("The Lipari-Szabo model combines S2, the internal timescale and the trajectory-wide correlation time to evaluate J at five frequencies."), origin};
        if (concept == QLatin1String("reorient.r1")
            || concept == QLatin1String("reorient.r2")
            || concept == QLatin1String("reorient.noe")) {
            return FamilyText{QStringLiteral("Standard nitrogen-15 relaxation equations combine the calculated spectral densities at the required frequencies."), origin};
        }
        if (concept == QLatin1String("reorient.s2"))
            return FamilyText{QStringLiteral("After aligning frames, the calculation averages products of unit bond-vector components. S2 is 1.5 times the sum of the squared averages, counting off-diagonal terms twice, minus 0.5."), origin};
        if (concept == QLatin1String("reorient.tau_e"))
            return FamilyText{QStringLiteral("The calculation sums (C_internal - S2) / (1 - S2) up to the first crossing of S2 and multiplies by the sampling interval. It uses the full lag window if there is no crossing."), origin};
        if (concept == QLatin1String("reorient.acf_internal"))
            return FamilyText{QStringLiteral("After aligning frames to the reference, the calculation averages (3 cos^2(theta) - 1) / 2 at each lag, where theta is the angle between the bond directions."), origin};
        return FamilyText{QStringLiteral("Without removing overall rotation, the calculation averages (3 cos^2(theta) - 1) / 2 at each lag, where theta is the angle between the bond directions."), origin};
    }
    if (family == QLatin1String("lipari_szabo"))
        return FamilyText{QStringLiteral("The matrix of average P2 correlations between N-H directions is diagonalised. Each bond's order parameter sums its squared projections onto the five largest modes, weighted by their eigenvalues."),
                          QStringLiteral("iRED calculator in nmr_extract.")};
    if (family == QLatin1String("rmsd"))
        return FamilyText{QStringLiteral("Each frame is fitted to the reference before calculating the root-mean-square atomic displacement."),
                          QStringLiteral("Trajectory alignment in nmr_extract.")};

    if (family == QLatin1String("topology")) {
        if (concept == QLatin1String("bond_length.stats"))
            return FamilyText{QStringLiteral("nmr_extract measures the distance between each bonded atom pair in every frame."),
                              QStringLiteral("Bond-length statistics calculator in nmr_extract.")};
        if (concept == QLatin1String("geometry.bond_length")) {
            return FamilyText{QStringLiteral("Reader measures the distance between the two atoms named by each bond record."),
                              QStringLiteral("Atomic coordinates and bond records from nmr_extract.")};
        }
        return FamilyText{QStringLiteral("Reader loads these records with the run."),
                          QStringLiteral("Topology sidecar from nmr_extract.")};
    }
    if (family == QLatin1String("selections")) {
        if (descriptor.storagePath == QLatin1String("/trajectory/selections"))
            return FamilyText{QStringLiteral("Reader loads the extraction events and counts those recorded at each trajectory frame."),
                              QStringLiteral("Selection events stored by nmr_extract.")};
        return FamilyText{QStringLiteral("Reader counts the atoms in its current selection. This is a current count, not a history of selection changes."),
                          QStringLiteral("Reader atom selection.")};
    }
    if (family == QLatin1String("gromacs"))
        return FamilyText{QStringLiteral("The values come from the energy and runtime records for the trajectory."),
                          QStringLiteral("GROMACS, recorded by nmr_extract.")};
    if (family == QLatin1String("bonded"))
        return FamilyText{QStringLiteral("The calculation evaluates the force-field bonded interactions and divides each interaction's energy equally among its participating atoms."),
                          QStringLiteral("Bonded-energy calculator in nmr_extract.")};
    if (family == QLatin1String("mutation_delta"))
        return FamilyText{QStringLiteral("Wild-type and mutant values are matched at the same site and subtracted in the stated order."),
                          QStringLiteral("nmr_extract mutation-pair comparison.")};

    if (family == QLatin1String("orca"))
        return FamilyText{QStringLiteral("Reader separates the Cartesian shielding tensor into its isotropic, antisymmetric and symmetric traceless parts."),
                          QStringLiteral("ORCA shielding output, read by Reader.")};
    if (family == QLatin1String("experimental_shielding_ml")) {
        const QString origin = QStringLiteral("Reader's bundled equivariant model for predicting ORCA total shielding.");
        if (concept == QLatin1String("experimental_shielding_ml.t2_norm"))
            return FamilyText{QStringLiteral("Reader supplies the extracted atomic features to the model, then takes the Euclidean norm of its five predicted T2 components."), origin};
        if (concept == QLatin1String("experimental_shielding_ml.iso"))
            return FamilyText{QStringLiteral("Reader supplies the extracted atomic features to the model and reads its predicted isotropic shielding output."), origin};
        return FamilyText{QStringLiteral("Reader supplies the extracted atomic features to the model and reads its five symmetric traceless shielding components."), origin};
    }

    return std::nullopt;
}

}  // namespace

std::optional<MetricGlossaryEntry> MetricGlossaryFor(const SignalDescriptor& descriptor) {
    const std::optional<FamilyText> family = familyText(descriptor);
    if (!family)
        return std::nullopt;

    QString calculation = family->calculation;
    const QString form = formSentence(descriptor);
    if (!form.isEmpty())
        calculation += QLatin1Char(' ') + form;

    MetricGlossaryEntry entry;
    entry.meaning = meaningFor(descriptor);
    entry.calculation = calculation;
    entry.origin = family->origin;
    if (entry.meaning.isEmpty() || entry.calculation.isEmpty() || entry.origin.isEmpty())
        return std::nullopt;
    return entry;
}

}  // namespace h5reader::model

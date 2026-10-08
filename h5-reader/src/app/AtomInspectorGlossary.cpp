#include "AtomInspectorGlossary.h"

namespace h5reader::app {
namespace {
using model::MetricGlossaryEntry;

std::optional<MetricGlossaryEntry> catalogHelp(const char* family, const char* concept,
                                              const QString& label) {
    model::SignalDescriptor descriptor;
    descriptor.family = QString::fromLatin1(family);
    descriptor.conceptKey = QString::fromLatin1(concept);
    descriptor.label = label;
    return model::MetricGlossaryFor(descriptor);
}

std::optional<MetricGlossaryEntry> identityHelp(const QString& field) {
    static const struct { const char* field; const char* meaning; const char* calculation; } entries[] = {
        {"Element", "Chemical element of the atom.", "Reader displays the element assigned in the topology."},
        {"AMBER name", "Atom name in the AMBER naming scheme.", "Reader reads the atom's AMBER name from the name sidecar."},
        {"IUPAC name", "Atom name in the IUPAC naming scheme.", "Reader reads the atom's IUPAC name from the name sidecar."},
        {"BMRB name", "Atom name used for BMRB assignments.", "Reader reads the atom's BMRB name from the name sidecar."},
        {"Backbone role", "Numeric code identifying the atom's backbone role.", "Codes: 0 none, 1 N, 2 CA, 3 carbonyl C, 4 carbonyl O, 5 amide H, 6 alpha H."},
        {"Locant", "Position within the residue, such as alpha or beta.", "Codes: 0 none, 1 alpha, 2 beta, 3 gamma, 4 delta, 5 epsilon, 6 zeta, 7 eta."},
        {"Residue", "Residue type and sequence number containing this atom.", "Reader follows the atom's residue reference and displays the residue's name and number."},
        {"Chain", "Chain containing this atom's residue.", "Reader displays the chain identifier from the residue record."},
        {"Protonation variant", "Variant of the residue's protonation state.", "Reader displays the stored variant index, or default when no variant is assigned."},
        {"Covalent radius", "Element-based covalent radius, in angstroms.", "Reader looks up the radius for this element; it is not measured from the current bonds."},
        {"Formal charge", "Integer charge assigned to the atom in the topology, not a partial charge.", "Reader displays the stored formal charge in elementary-charge units."},
        {"is_backbone", "Whether the atom has a backbone role.", "True when BackboneRole is not None."},
        {"is_amide_H", "Whether the atom is classified as a backbone amide hydrogen.", "True when PolarHKind is BackboneAmide."},
        {"is_alpha_H", "Whether this is an alpha hydrogen, including glycine HA2 and HA3.", "Reader tests the alpha-hydrogen role or the glycine-compatible hydrogen and alpha-locant combination."},
        {"is_methyl", "Whether the atom belongs to the stored methyl grouping.", "True when PseudoatomKind is M."},
        {"is_aromatic", "Whether the atom is marked aromatic.", "Reader displays the topology's aromatic flag."},
        {"is_polar_H", "Whether the atom is classified as a polar hydrogen.", "True when PolarHKind is not NotPolar."},
        {"is_hbond_acceptor_elem", "Whether the element is N or O. This is not a complete test of hydrogen-bond acceptance.", "Reader tests the element only, not protonation or bonding."},
        {"is_exchangeable", "Whether the atom is marked as exchangeable.", "Reader displays the topology's exchangeability flag."},
    };
    for (const auto& entry : entries) {
        if (field == QLatin1String(entry.field))
            return MetricGlossaryEntry{QString::fromLatin1(entry.meaning), QString::fromLatin1(entry.calculation),
                QStringLiteral("nmr_extract topology and name records; Reader's typed atom model.")};
    }
    return std::nullopt;
}

std::optional<MetricGlossaryEntry> tensorHelp(const QString& field, bool bond,
                                              const QString& source) {
    const QString origin = bond
        ? QStringLiteral("Reorientational dynamics in nmr_extract; Reader calculates the principal values and draws the glyph.")
        : QStringLiteral("%1; Reader calculates the principal values and draws the glyph.").arg(source);
    const auto entry = [&origin](const char* meaning, const char* calculation) {
        return MetricGlossaryEntry{QString::fromLatin1(meaning), QString::fromLatin1(calculation), origin};
    };
    if (field == QLatin1String("Atom"))
        return MetricGlossaryEntry{
            QStringLiteral("Atom represented by this shielding tensor."),
            QStringLiteral("Reader displays its chain, residue and IUPAC atom name."),
            QStringLiteral("nmr_extract topology and name records.")};
    if (field == QLatin1String("sigma_iso"))
        return entry("Mean shielding over field directions, in ppm. This is shielding, not chemical shift.",
                     "One third of the tensor trace, or the mean of its three principal values.");
    if (field.startsWith(QLatin1String("sigma_")) || field == QLatin1String("Principal values"))
        return entry("Shielding along the three coloured axes, in ppm. Parentheses show each value minus the mean.",
                     "Reader diagonalises the symmetric part of the tensor and orders the eigenvalues from smallest to largest. Each swatch identifies the corresponding axis.");
    if (field == QLatin1String("span"))
        return entry("Difference between the largest and smallest principal values.", "sigma_33 - sigma_11.");
    if (field == QLatin1String("skew"))
        return entry("Position of the middle principal value relative to the mean and span.", "3 (sigma_22 - sigma_iso) / span.");
    if (field == QLatin1String("eta"))
        return entry("Asymmetry of the tensor, using Haeberlen ordering.",
                     "Order axes by distance from the mean: zz farthest, xx next, yy nearest. Then eta = (sigma_yy - sigma_xx) / (sigma_zz - sigma_iso).");
    if (field.startsWith(QLatin1String("lambda_")))
        return entry("Mean squared projection of the unit bond vector on this axis. The three values sum to one; they are not fractions of frames.",
                     "Reader diagonalises the stored average of u u^T and orders the eigenvalues from largest to smallest. Each swatch identifies the corresponding axis.");
    if (field == QLatin1String("S^2 (order parameter)"))
        return entry("Order of the bond direction after removing overall rotation. A fixed direction gives one; isotropic directions give zero.",
                     "S^2 = (3 (lambda_1^2 + lambda_2^2 + lambda_3^2) - 1) / 2.");
    if (field == QLatin1String("Main axis follows"))
        return entry("The glyph follows the current bond, while its shape describes the trajectory average.",
                     "Reader rotates the largest-eigenvalue axis onto the current bond. The displayed arrows are not the average tensor's axes in the aligned protein.");
    if (field == QLatin1String("Glyph size"))
        return entry("Glyphs are scaled separately. Sizes cannot be compared as absolute magnitudes between atoms.",
                     "Arrow lengths use absolute deviations of the principal values from their mean. Reader normalises them separately for each tensor and keeps a minimum length for visibility.");
    if (field == QLatin1String("Bond"))
        return entry("Bond used to calculate this orientation tensor.", "Reader displays the stored bond endpoints using chain, residue and IUPAC atom names.");
    if (field == QLatin1String("Average over"))
        return entry("These values describe the trajectory, not just the current frame.", "nmr_extract aligns frames and averages u u^T, where u is the unit bond vector.");
    if (field == QLatin1String("Scope"))
        return entry("Frame represented by the shielding values.", "Reader uses the shielding tensor for the named atom at the displayed frame.");
    if (field == QLatin1String("Source"))
        return entry("Calculation or model supplying the tensor.", "Reader displays the method or model identifier attached to the tensor.");
    if (field == QLatin1String("Frame convention"))
        return entry("Coordinate frame in which the supplied tensor is expressed.", "Reader uses this convention when rotating the tensor into the displayed molecule.");
    if (bond)
        return entry("Variation in a bond's direction across the aligned trajectory.", "nmr_extract averages u u^T for unit bond vectors. Reader diagonalises this matrix to obtain its three principal values.");
    return entry("Magnetic shielding at the atom for different directions of the applied field.", "Reader diagonalises (sigma + sigma^T) / 2. Its mean and deviations describe the isotropic and directional parts.");
}
MetricGlossaryEntry tensorComponentHelp(const QString& field, MetricGlossaryEntry help, bool efg) {
    QString meaning;
    QString calculation;
    if (field == QLatin1String("T0") || field == QLatin1String("T0 signed iso")) {
        meaning = QStringLiteral("Isotropic scalar part of this tensor, in the displayed units.");
        calculation = QStringLiteral("T0 = trace(T) / 3.");
    } else if (field.startsWith(QLatin1String("|T2|"))) {
        meaning = QStringLiteral("Magnitude of the symmetric traceless part, independent of coordinate orientation.");
        calculation = QStringLiteral("Square root of the sum of the five squared T2 components. The library basis preserves the Cartesian Frobenius norm.");
    } else if (field == QLatin1String("|T1| antisymmetric") || field == QLatin1String("T1 antisym vector")) {
        meaning = QStringLiteral("Antisymmetric part of this tensor, represented by three components.");
        calculation = QStringLiteral("T1 = (Tyz-Tzy, Tzx-Txz, Txy-Tyx)/2. The magnitude is the square root of the three squared components' sum.");
    } else if (field.startsWith(QLatin1String("m=")) || field.contains(QLatin1String("T2 components"))) {
        meaning = QStringLiteral("Components of the symmetric traceless tensor in the extractor's library basis.");
        calculation = QStringLiteral("Remove the trace from the symmetric matrix S. In order m=-2 to 2: sqrt(2) Sxy, sqrt(2) Syz, sqrt(3/2) Szz, sqrt(2) Sxz, (Sxx-Syy)/sqrt(2).");
    } else if (field == QLatin1String("span")) {
        meaning = QStringLiteral("Difference between the largest and smallest principal values.");
        calculation = QStringLiteral("Largest minus smallest eigenvalue of the symmetric matrix.");
    } else if (field == QLatin1String("skew")) {
        meaning = QStringLiteral("Position of the middle principal value relative to the mean and span.");
        calculation = QStringLiteral("3 (middle eigenvalue - trace(T)/3) / span.");
    } else if (field == QLatin1String("eta")) {
        meaning = QStringLiteral("Asymmetry of the tensor, using Haeberlen ordering.");
        calculation = efg
            ? QStringLiteral("Order the traceless tensor's axes by absolute eigenvalue: zz largest, xx next, yy smallest. Then eta = (Vyy - Vxx) / Vzz.")
            : QStringLiteral("Order axes by distance from the mean: zz farthest, xx next, yy nearest. Then eta = (Tyy - Txx) / (Tzz - trace(T)/3).");
    } else if (field.startsWith(QLatin1String("sigma_")) || field.startsWith(QLatin1String("V_"))) {
        meaning = QStringLiteral("Principal value of the symmetric tensor, in the displayed units.");
        calculation = QStringLiteral("Reader diagonalises the symmetric matrix and orders its eigenvalues from smallest (11) to largest (33).");
    } else if (field == QLatin1String("Vxx") || field == QLatin1String("Vyy") || field == QLatin1String("Vzz")) {
        meaning = QStringLiteral("Principal value of the traceless electric-field-gradient tensor.");
        calculation = QStringLiteral("Reader orders the eigenvalues by absolute size: zz largest, xx next, yy smallest.");
    } else if (field.startsWith(QLatin1String("PAS /"))) {
        meaning = QStringLiteral("Principal-axis values and measures of tensor shape.");
        calculation = QStringLiteral("Reader reconstructs the symmetric matrix from T0 and T2, then diagonalises it.");
    } else if (field == QLatin1String("raw irreps")) {
        meaning = QStringLiteral("The tensor separated into an isotropic scalar, three antisymmetric components and five symmetric traceless components.");
        calculation = QStringLiteral("T0 is trace(T)/3. T1 represents (T-T^T)/2; T2 represents (T+T^T)/2 - T0 I.");
    }
    if (!meaning.isEmpty()) {
        help.meaning = meaning + QLatin1Char(' ') + help.meaning;
        help.calculation = calculation;
        help.origin += QStringLiteral(" Tensor decomposition and shape values are displayed by Reader.");
    }
    return help;
}
}  // namespace

std::optional<MetricGlossaryEntry> AtomInspectorGlossary(
    const QStringList& path, const QString& shieldingSource) {
    if (path.isEmpty()) return std::nullopt;
    const QString& field = path.back();
    for (const QString& section : path) {
        if (section.startsWith(QLatin1String("Shielding tensor (")))
            return tensorHelp(field, false, shieldingSource);
        if (section == QLatin1String("Bond orientation tensor"))
            return tensorHelp(field, true, {});
    }
    if (path.contains(QStringLiteral("Identity"))) {
        if (auto help = identityHelp(field)) return help;
        return MetricGlossaryEntry{QStringLiteral("Names, chemical roles and flags identifying the atom."),
            QStringLiteral("Reader reads the topology and resolves the atom's residue and naming records."),
            QStringLiteral("nmr_extract topology and name sidecars.")};
    }
    if (path.contains(QStringLiteral("Position")))
        return MetricGlossaryEntry{QStringLiteral("Position of the atom in the displayed conformation, in angstroms."),
            QStringLiteral("Reader reads the frame coordinates and applies the active molecular alignment."),
            QStringLiteral("Extracted coordinates and Reader's conformation transform.")};

    if (path.contains(QStringLiteral("Water"))) {
        static const struct { const char* field; const char* meaning; const char* calculation; const char* origin; } water[] = {
            {"half-shell asymmetry", "Fraction of first-shell waters on the side away from the protein centre.",
             "Count waters in the outward half-shell and divide by the first-shell count.", "Hydration calculator in nmr_extract."},
            {"mean water dipole cos", "Mean alignment of water dipoles with the atom-to-water direction.",
             "Average the cosine between each first-shell water dipole and its atom-to-water vector.", "Hydration calculator in nmr_extract."},
            {"nearest ion", "Distance and charge of the nearest ion within the search cutoff.",
             "Find the nearest accepted ion. Distance is in angstroms; charge is in elementary-charge units.", "Hydration calculator in nmr_extract."},
            {"dipole alignment", "Alignment of the net first-shell water dipole with the exposed surface normal.",
             "Cosine of the angle between the summed water dipoles and the SASA surface normal.", "Water-polarization result from nmr_extract."},
            {"dipole coherence", "Magnitude of the mean first-shell water dipole. Opposing dipoles cancel.",
             "Take the length of the summed dipole vector and divide by the number of first-shell waters.", "Water-polarization result from nmr_extract."},
        };
        for (const auto& entry : water)
            if (field == QLatin1String(entry.field))
                return MetricGlossaryEntry{QString::fromLatin1(entry.meaning), QString::fromLatin1(entry.calculation), QString::fromLatin1(entry.origin)};
    }

    // These rows are alternate presentations of the existing metric glossary.
    static const struct { const char* field; const char* family; const char* concept; } aliases[] = {
        {"bs_shielding", "biot_savart", "bs_shielding"},
        {"hm_shielding", "haigh_mallion", "hm_shielding"},
        {"bs_total_B", "biot_savart", "bs_total_B"},
        {"coulomb_shielding", "coulomb", "coulomb_efg"},
        {"coulomb_E", "coulomb", "coulomb_E"},
        {"apbs_E (APBS diagnostic)", "apbs", "apbs_E"},
        {"apbs_efg (APBS diagnostic)", "apbs", "apbs_efg"},
        {"aimnet2_efg", "aimnet2", "aimnet2_efg"},
        {"atom_sasa", "sasa", "atom_sasa"},
        {"surface normal", "sasa", "sasa_normal"},
        {"water_efield", "water_field", "water_efield"},
        {"water_efg", "water_field", "water_efg"},
        {"shell counts (1st/2nd)", "water_field", "water_shell_counts"},
        {"AIMNet2 (Hirshfeld)", "aimnet2", "aimnet2_charges"},
        {"EEQ", "eeq", "eeq_charges"},
        {"EEQ coord. number", "eeq", "eeq_cn"},
        {"|charge-response grad|", "aimnet2", "aimnet2_charge_response_gradient_scalar"},
        {"pyramidalization", "planar_geometry", "pyramidalization"},
        {"charge", "mopac", "mopac_charges_full_precision"},
        {"s population", "mopac", "mopac_atom_s_population"},
        {"p population", "mopac", "mopac_atom_p_population"},
        {"Wiberg valency", "mopac", "mopac_bond_valencies_full_precision"},
        {"mopac_coulomb_shielding", "mopac_coulomb", "mopac_coulomb_efg"},
        {"mopac_mc_shielding", "mopac_mcconnell", "mopac_mc_shielding"},
        {"water term", "larsen_hbond", "larsen_hbond_water_term"},
        {"H-bond pair count", "larsen_hbond", "larsen_hbond_count"},
        {"\xCF\x83 total", "orca", "orca_total"},
        {"\xCF\x83 diamagnetic", "orca", "orca_diamagnetic"},
        {"\xCF\x83 paramagnetic", "orca", "orca_paramagnetic"},
        {"\xCE\x94\xCF\x83 total", "larsen_hbond", "larsen_hbond_shielding"},
    };
    for (const auto& alias : aliases) {
        if (field == QString::fromUtf8(alias.field))
            return catalogHelp(alias.family, alias.concept, field);
    }

    if (path.contains(QStringLiteral("H-bond"))) {
        auto help = catalogHelp("hbond", "hbond_scalars", field);
        if (help && field == QLatin1String("nearest dist")) {
            help->meaning = QStringLiteral("Distance from this atom to the nearest accepted donor H, in angstroms.");
            help->calculation = QStringLiteral("Measure distances to the accepted donor hydrogens and select the shortest.");
        } else if (help && field == QStringLiteral("1/r\u00b3")) {
            help->meaning = QStringLiteral("Inverse cube of the distance to the nearest accepted donor H, in inverse cubic angstroms.");
            help->calculation = QStringLiteral("1 / r^3, using the nearest donor-H distance.");
        } else if (help && field == QStringLiteral("count \u2264 3.5 \u00c5")) {
            help->meaning = QStringLiteral("Number of accepted hydrogen-bond donor hydrogens within 3.5 angstroms of this atom.");
            help->calculation = QStringLiteral("Count accepted sources at or below the stated distance.");
        }
        return help;
    }
    if (path.contains(QStringLiteral("SASA")))
        return catalogHelp("sasa", "atom_sasa", field);
    if (path.contains(QStringLiteral("Bonded energy (per-atom share)"))) {
        auto help = catalogHelp("bonded", "bonded_energy", field);
        if (help) help->meaning = field == QLatin1String("total")
            ? QStringLiteral("Sum of this atom's shares of the bonded interaction energies, in kJ/mol.")
            : QStringLiteral("This atom's share of the %1 interaction energy, in kJ/mol.").arg(field);
        return help;
    }
    if (path.contains(QStringLiteral("Frame energy (GROMACS)"))) {
        auto help = catalogHelp("gromacs", "gromacs_energy", field);
        if (help) help->meaning = QStringLiteral("%1 of the whole simulated system, not an atomic contribution.").arg(field);
        return help;
    }
    if (path.contains(QStringLiteral("DSSP (secondary structure)"))) {
        if (field == QStringLiteral("\u03c6 (neg-IUPAC)") || field == QStringLiteral("\u03c8 (neg-IUPAC)"))
            return MetricGlossaryEntry{QStringLiteral("Backbone torsion for this atom's residue, in radians."),
                QStringLiteral("nmr_extract reads DSSP phi and psi and stores them in its negative-IUPAC convention."),
                QStringLiteral("DSSP backbone results from nmr_extract.")};
        if (field == QLatin1String("residue SASA"))
            return MetricGlossaryEntry{QStringLiteral("Solvent-accessible area of the whole residue, not this atom alone."),
                QStringLiteral("Reader displays DSSP's residue area for each atom in that residue."), QStringLiteral("DSSP, run by nmr_extract.")};
        auto help = catalogHelp("dssp", "dssp_ss8", field);
        if (help) help->meaning += QStringLiteral(" Codes: 0 alpha helix, 1 3-10 helix, 2 pi helix, 3 extended strand, 4 beta bridge, 5 turn, 6 bend, 7 coil, 255 unknown.");
        return help;
    }
    if (field == QStringLiteral("\u03c9 (peptide)")) return catalogHelp("planar_geometry", "omega_actual", field);
    if (field == QStringLiteral("\u03c9 deviation")) return catalogHelp("planar_geometry", "omega_deviation", field);
    if (field == QStringLiteral("X\u2192Pro context"))
        return MetricGlossaryEntry{QStringLiteral("Whether this peptide bond leads into proline."),
            QStringLiteral("nmr_extract tests the following residue's type."), QStringLiteral("Planar geometry in nmr_extract.")};
    if (field == QStringLiteral("\u0394Hf (frame)"))
        return MetricGlossaryEntry{QStringLiteral("Heat of formation of the whole conformation, in kcal/mol. It is not an atomic contribution."),
            QStringLiteral("Reader displays the heat of formation reported by MOPAC."), QStringLiteral("MOPAC PM7/MOZYME output from nmr_extract.")};
    if (path.contains(QStringLiteral("Ring geometry")))
        return MetricGlossaryEntry{QStringLiteral("Numbers of nearby rings within 3, 5, 8 and 12 angstroms."),
            QStringLiteral("Reader displays the four distance counts from the Biot-Savart result."), QStringLiteral("Biot-Savart calculator in nmr_extract.")};

    if (path.contains(QStringLiteral("Local classical estimate (tentative)"))) {
        static const struct { const char* field; const char* calculation; } estimates[] = {
            {"estimated ring contribution", "Sum the Biot-Savart scalar for each ring type times its Giessner-Prettre ring intensity."},
            {"estimated Larsen contribution", "Add the Larsen hydrogen-bond T0 and the separate ProCS15 water correction when present."},
            {"estimated McConnell contribution", "Sum the available bond-order-weighted McConnell T0 terms times their susceptibility anisotropies and molar conversion factor."},
            {"estimated Buckingham contribution", "Evaluate -A E_parallel - B E_parallel^2 using the signed MOPAC bond-axis field and the literature A and B constants."},
            {"estimated sigma_cl (tentative)", "Add sigma0 and the ring, McConnell, Larsen and Buckingham contributions."},
            {"tentative residual (sigma_qm - estimate)", "Subtract the local classical estimate from ORCA's total isotropic shielding."},
        };
        for (const auto& estimate : estimates)
            if (field == QLatin1String(estimate.field))
                return MetricGlossaryEntry{QStringLiteral("Local shielding estimate in ppm, not a measured value or a validated absolute shielding model."),
                    QString::fromLatin1(estimate.calculation), QStringLiteral("Reader's classical calculation using extracted values and literature constants.")};
    }
    if (field == QLatin1String("signed E_parallel"))
        return MetricGlossaryEntry{QStringLiteral("Electric field along the parent-to-H bond, in V/angstrom. The sign records its direction."),
            QStringLiteral("Dot product of the MOPAC-charge Coulomb field and the parent-to-H unit vector."),
            QStringLiteral("MOPAC Coulomb calculator in nmr_extract.")};
    if (field == QLatin1String("EFG |T2| (AIMNet2)")) {
        if (auto help = catalogHelp("aimnet2", "aimnet2_efg", field))
            return tensorComponentHelp(QStringLiteral("|T2| invariant"), *help, true);
    }

    // The parent identifies the source; a tensor child identifies the operation.
    for (const QString& ancestor : path) {
        for (const auto& alias : aliases) {
            if (ancestor == QString::fromUtf8(alias.field)) {
                auto help = catalogHelp(alias.family, alias.concept, ancestor);
                if (!help) return std::nullopt;
                return tensorComponentHelp(field, *help, path.contains(QStringLiteral("PAS / EFG convention")));
            }
        }
    }
    if (path.contains(QStringLiteral("Bond anisotropy (McConnell)"))) {
        if (auto help = catalogHelp("mcconnell", "mc_shielding", field))
            return tensorComponentHelp(field, *help, false);
    }
    return std::nullopt;
}

}  // namespace h5reader::app

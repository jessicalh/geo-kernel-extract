# DCHI Compute Pilot

Date: 2026-06-02

## Outcome

This pilot did **not** produce replacement `Delta chi` values. ORCA is installed, but ORCA 6.1.1 does not expose a documented public input keyword/block for the requested closed-shell **full molecular magnetizability / magnetic susceptibility tensor** at the SCF/DFT GIAO level.

I therefore did not run the production formaldehyde/N-methylacetamide optimizations, did not calculate final anisotropies, did not touch any C++ kernel, and did not rerun any McConnell correlation.

## ORCA Availability

Installed executable:

```text
/opt/orca/orca
```

Version banner:

```text
Program Version 6.1.1  -  RELEASE   -
(GIT: $487d211c$)
($2025-11-21 10:33:24 +0100$)
```

## Manual Check

Manual checked: ORCA 6.1.1 manual, hosted at:

- https://orca-manual.mpi-muelheim.mpg.de/
- PDF used for text search: https://orca-manual.mpi-muelheim.mpg.de/_static/pdf/orca.pdf

Relevant manual findings:

- The NMR section documents GIAO magnetic-field response for shielding tensors through `%eprnmr` and states that shieldings use Gauge Including Atomic Orbitals.
- The CP-SCF section states that magnetic response solver controls are in `%eprnmr`.
- The electric-property section documents `%elprop` for dipole/quadrupole/polarizability/hyperpolarizability.
- The QDPT magnetic-property sections document `DoMagnetization` and `DoSusceptibility` for MRCI/CASSCF/SINGLE_ANISO-style magnetic susceptibility workflows, not the closed-shell molecular GIAO magnetizability tensor requested here.
- A full-text search of the 6.1.1 manual for `magnetizability` found no documented public keyword.

## Syntax Probes

ORCA input/output probes are kept in:

```text
src/rediscover/dchi_compute_pilot_orca/
```

Rejected public input attempts:

```text
%elprop
  Magnetizability true
end
```

```text
%eprnmr
  Magnetizability true
end
```

```text
! HF STO-3G TightSCF Magnetizability
```

```text
%prop
  Magnetizability true
end
```

ORCA rejected these with unknown-identifier/input errors. Example outputs:

- `syntax_elprop_mag.out`
- `syntax_eprnmr_mag.out`
- `syntax_simple_mag.out`
- `syntax_prop_mag.out`

## Hidden `orca_prop` Check

`orca_prop` contains an internal `Magnetizability` property slot. I verified this by capturing ORCA's generated property-driver input from a small GIAO NMR job. The generated `propinp` line labels the slot:

```text
... NMRSpinSpin, Magnetizability OptDSOC OptDSS SOCEnergy DBOC
```

Using a GIAO NMR response job and manually setting only that internal `Magnetizability` flag allowed:

```text
/opt/orca/orca_prop syntax_nmr_for_mag_cap_magnetizability.propinp syntax_nmr_for_mag_cap 1
```

This printed:

```text
MAGNETIZABILITY TENSOR
Paramagnetic magnetizability tensor (au):
...
Isotropic (paramagnetic) magnetizability :   3.12106
```

This is **not acceptable for the requested protocol**:

- It is not a documented ORCA input syntax.
- It requires preserving and editing ORCA temporary driver files.
- It prints only the **paramagnetic** magnetizability tensor in atomic units.
- It does not print the full molecular tensor, the diamagnetic contribution, or a total tensor.

Because diamagnetic molecular susceptibilities are negative and the requested McConnell quantity needs the full susceptibility anisotropy, using this paramagnetic-only tensor would be a units/sign failure.

## Units

The hidden ORCA output unit is `au` for the **paramagnetic** magnetizability tensor. No validated ORCA output unit was found for a full molecular susceptibility tensor because no documented full-tensor protocol was found.

Useful conversion if a future ORCA/full-tensor route reports magnetizability in atomic units:

```text
1 a.u. magnetizability = 7.891036607e-29 J T^-2
```

Source for the atomic unit value: NIST/CODATA constants listing for `atomic unit of magnetizability` (https://physics.nist.gov/cgi-bin/cuu/Results?search_for=atomic+unit). Molar cgs conversion:

```text
chi_m(SI, m^3/mol) = mu_0 * N_A * xi(SI, J/T^2)
chi_m(cgs, cm^3/mol) = chi_m(SI) / (4*pi*1e-6)
```

Therefore:

```text
1 a.u. = 4.7520 x 10^-6 cm^3/mol
```

This conversion was **not applied to production `Delta chi` values**, because no full ORCA tensor was obtained.

## Formaldehyde Sanity Anchor

Known experimental formaldehyde molecular-beam anisotropy data are available from Kukolich, J. Chem. Phys. 54, 105-110 (1971), DOI: https://doi.org/10.1063/1.1674644. The abstract reports susceptibility anisotropy combinations in `10^-6 erg/G^2 mol`:

```text
2 chi_aa - chi_bb - chi_cc = 25.4 +/- 0.2
2 chi_bb - chi_aa - chi_cc = -4.1 +/- 0.15
```

Those are the correct kind of cgs molar susceptibility-anisotropy quantities for the requested sanity check. The hidden ORCA paramagnetic-only output cannot be compared to them as a total susceptibility.

## Sign Convention

For the McConnell kernel, the relevant `q` must be a **full molar magnetic susceptibility anisotropy** in `10^-6 cm^3/mol`. Diamagnetic susceptibilities are negative. The existing kernel form

```text
sigma_ppm = -0.5535 * q * f
```

therefore depends on the sign of `q` after conversion to the same cgs molar convention. A paramagnetic-only ORCA tensor would have the wrong physical content and could flip or distort the McConnell sign/magnitude.

## Basis And Molecule Runs

Not run:

- formaldehyde r2SCAN/def2-SVP/CPCM(water)
- formaldehyde r2SCAN/def2-TZVP/CPCM(water)
- trans N-methylacetamide r2SCAN/def2-SVP/CPCM(water)

Reason: Step 0 protocol verification failed for the requested **full** molecular magnetizability tensor. Running optimized structures and reporting paramagnetic-only tensors would create exactly the unit/sign ambiguity this pilot is meant to prevent.

## Recommendation

Do not wire computed `Delta chi` values from ORCA 6.1.1 yet. Route the calculation to a package/protocol that explicitly reports the full molecular magnetizability tensor, including diamagnetic and paramagnetic terms, or get confirmation from ORCA/FACCTs for the supported 6.1.1 input syntax and whether `orca_prop` can output the total tensor.

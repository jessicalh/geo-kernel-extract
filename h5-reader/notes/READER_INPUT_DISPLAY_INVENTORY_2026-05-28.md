# Reader Input and Display Inventory

Date: 2026-05-28

> **Snapshot only (trued 2026-06-04).** This inventory predates later source,
> catalog, and dashboard evolution. Verify all counts and drift notes against
> current source before acting; keep the file as historical inventory.

Pass: 1 of 3.

Purpose: collect the full h5-reader import surface aggressively before
designing the `StripCalculation` class tree. This is not the consumer design
yet. It is the map of what the consumer design must cover.

Implementation guard: this inventory is the catalog target, not background
reading. Early rendering can be staged, but the first catalog/dialog pass should
be shaped around this complete source surface so that demo-first UI assumptions
do not fossilize.

## Rules

- Do not turn storage format into ontology. NPY, dense H5, ORCA/DFT, and
  geometry are different import/export sets and may be read differently, but
  they often represent overlapping measurement families.
- Do not force unification where the export sets are genuinely different.
  Instead, track overlap explicitly.
- Do unify the consumer side. The working data surface is trajectories plus
  frame snapshots, with calculated values derived from those sources. Every
  source is either startup-loaded, frame-loaded, or derived, and every source
  exposes data sets that can be described by the same axes, value shapes, gap
  semantics, and display modes. The class tree should encode those
  commonalities without becoming a plugin framework.
- Every imported signal family needs two display-form answers:
  temporal strip recorder candidates and static/current-frame display
  candidates.
- A "strip" means persistent sampled display history owned by a
  `StripCalculation`. Source readers provide frame-local or imported data;
  they are not the user's trail.
- Rollups and precomputed statistics are imported signals too, but many are
  static/statistical rather than frame-time strips.

## Authoritative Reader Surfaces

### SDK / NPY Catalog

Source: `python/nmr_extract/_catalog.py`.

This is the source of truth for every NPY file the extractor can emit. Current
counts:

- 109 fields
- 28 groups
- 10 native axes
- Generated C++ mirror: `h5-reader/src/io/QtFieldCatalog.gen.h`
- Drift found: SDK has `gromacs_energy` with 43 columns; generated C++ catalog
  currently has 42 columns. `FrameNpyLoader` already trusts the NPY header and
  warns on this mismatch.

Native axes:

- atom: 84
- mutation_match_pair: 10
- residue: 4
- aromatic_ring: 2
- saturated_ring: 2
- protein: 2
- bond: 2
- ring_contribution_pair: 1
- ring: 1
- ring_membership: 1

Mechanism counts:

- electrostatic_efg: 15
- topology: 10
- ring_current: 10
- charges: 10
- quantum_reference: 10
- mutation_delta: 10
- solvation: 9
- hbond_grid: 8
- geometry: 7
- bond_anisotropy: 6
- secondary_structure: 4
- ring_efg: 3
- ring_dispersion: 3
- hbond_kernel: 2
- gromacs_runtime: 2

### h5-reader Snapshot Reader

Source: `FrameNpyLoader`, `QtConformationSnapshot`, and `Qt*Group` views.

Contract:

- One frame directory is loaded as frame-local source data.
- Filename stem resolves once through `QtFieldCatalog.gen.h`.
- Topology sidecar fields are not stored as numeric calculator columns.
- Present columns are stored by `FieldKind`.
- Typed group views decode columns semantically and generally return
  `std::optional`.
- Short or drifted writer output is handled at the reader boundary where known
  (`hbond_scalars` and `gromacs_energy` are current examples).

Snapshot group views currently cover:

- ring-current kernels: `QtBiotSavartGroup`, `QtHaighMallionGroup`,
  `QtRingSusceptibilityGroup`
- ring EFG/dispersion: `QtPiQuadrupoleGroup`, `QtDispersionGroup`
- bond anisotropy: `QtMcConnellGroup`, `QtMopacMcConnellGroup`
- electrostatics: `QtCoulombGroup`, `QtMopacCoulombGroup`, `QtApbsGroup`
- H-bond kernels/grids: `QtHBondGroup`, `QtLarsenHBondGroup`
- solvation: `QtSasaGroup`, `QtWaterFieldGroup`, `QtHydrationGroup`,
  `QtWaterPolarizationGroup`
- charges/quantum: `QtEeqGroup`, `QtAimnet2Group`, `QtOrcaGroup`,
  `QtMopacCoreGroup`, `QtTripeptideGroup`
- runtime/system: `QtGromacsGroup`, `QtBondedGroup`
- secondary structure/geometry: `QtDsspGroup`, `QtPlanarGeometryGroup`

### Dense H5 Trajectory Reader

Source: `QtTrajectoryH5`, `QtTimeSeriesBuffers`, `QtPerResidueBuffers`,
`QtSpecialBuffers`, and `QtWelfordBuffers`.

This is a different import set, not a separate conceptual universe. Some H5
groups overlap the NPY catalog; others are trajectory products that do not have
one-frame NPY equivalents.

Current dense reader groups:

- positions
- shielding time series: BS, HM, MC, piquad, ringchi, dispersion, H-bond,
  MOPAC Coulomb T2, MOPAC McConnell, tripeptide backbone, tripeptide neighbor,
  Larsen 1pHB/1pHaB/2pHB/2pHaB
- scalar time series: SASA, AIMNet2 charge, Larsen count, Larsen water term,
  bonded energy total, MOPAC-vs-FF14SB reconciliation
- composite atom time series: water field, hydration shell, hydration geometry,
  AIMNet2 charge-response gradient
- vector time series: APBS E-field, tripeptide residual vectors
- T2/tag/embedding time series: APBS EFG, tripeptide method tag, AIMNet2
  embedding
- residue time series: dihedrals, DSSP8, J-couplings
- ring/special time series: ring pucker, ring neighbourhood trajectory stats,
  RMSD tracking, Gromacs energy
- static/statistical H5 products: bond length stats, Welford rollups,
  autocorrelation, DSSP transitions, dihedral-bin transitions, selections

### ORCA / DFT Reader

Source: `DftShieldingStore`, `OrcaShieldingParser`, `DftShielding`.

Import set:

- ORCA shielding output per original trajectory frame.
- One parsed frame contains all atom total/diamagnetic/paramagnetic shielding
  tensors.
- It overlaps conceptually with `orca_total`, `orca_diamagnetic`, and
  `orca_paramagnetic` NPY fields, but it is a separate local reader path and
  should remain an orthogonal frame source.

### Geometry Reader

Source: `ConformationGeometry::Measure`, positions from `Conformation`.

Import set:

- derived from positions and topology/selection, not emitted as SDK NPY.
- temporal display is direct distance/angle/dihedral sampling.
- static display is measurement overlay and selected atom tuple labels.

## Display Primitive Codes

## Unifying Consumer Model

All fixed trajectory/snapshot sources should feed one consumer-side hierarchy.
The source is metadata, not the type of display object.

Common source residency:

- `StartupLoadedSource`: data is loaded when the run opens and can be sampled
  without per-frame source loading. Current examples: dense H5 time series,
  H5 rollups/stat products, topology sidecar, manifest/identity data.
- `FrameLoadedSource`: data for a single frame is loaded on demand, observers
  sample it, and source payload is released. Current examples: per-frame NPY
  snapshots and ORCA/DFT frames.
- `DerivedSource`: data is computed from already available source state, not
  imported as its own array. Current examples: selected-atom geometry
  distance/angle/dihedral and derived ring counts from neighbourhood data.

Common axes:

- atom
- residue
- atom tuple
- bond
- ring / aromatic ring / saturated ring
- ring contribution pair
- mutation match pair
- protein/system
- event/time marker

Common value shapes:

- scalar
- category/tag
- count
- vector3
- EFG/T2-only tensor
- spherical shielding tensor T0/T1/T2
- matrix/tensor component set
- per-class block
- embedding/high-dimensional vector
- rollup moments
- event record

Common temporal semantics:

- frame-time sample
- conditional/gapped frame-time sample
- system frame-time sample
- precomputed rollup/statistic with no frame axis
- event timeline
- derived transform of a sampled strip, such as FFT/spectrum

Therefore the class-tree pass should define base concepts like:

- `DisplaySignalDescriptor`
- `TrajectorySignalCatalog`
- fixed source helpers for dense H5, frame NPY snapshots, ORCA/DFT frames,
  topology, selections/events, and derived geometry
- `StripCalculation`
- `StaticDisplayCalculation`
- value-specific calculation bases for scalar, vector, tensor, category,
  per-class, rollup, event, and embedding signals

The design goal is not to hide the difference between NPY/H5/ORCA. The design
goal is that once a signal is described, the picker, strip recorder, static
panel, and overlay layers consume it through the same vocabulary.

Temporal strip recorder forms:

- `strip.scalar`: y = scalar over frame/time.
- `strip.count`: count/category count over frame/time.
- `strip.category`: categorical state over frame/time; display as step strip,
  occupancy strip, event/rug strip, or colored band.
- `strip.vector.component`: x/y/z component over frame/time.
- `strip.vector.magnitude`: vector norm over frame/time.
- `strip.tensor.T0`: isotropic scalar over frame/time.
- `strip.tensor.T1`: T1 components or T1 magnitude over frame/time.
- `strip.tensor.T2`: T2 components or T2 magnitude over frame/time.
- `strip.tensor.component`: explicit tensor component channel.
- `strip.per-class`: one channel per ring type, bond category, shell, residue
  class, coupling type, etc.
- `strip.system`: system-wide scalar over frame/time.
- `strip.rollup`: statistic channel such as mean/std/min/max/rms_delta; not a
  live frame trail unless the stat itself has a frame axis.
- `strip.event`: selection/transition events over time.
- `strip.spectrum`: FFT/power spectrum derived from a strip-compatible scalar
  signal. This is a consumer transform, not a new source signal.

Static/current-frame forms:

- `static.scalar`: inspector readout, table cell, atom/residue color map.
- `static.category`: label, color class, legend, residue/atom band.
- `static.vector`: arrow/glyph/field line; optional component table.
- `static.tensor`: T0 color, T1/T2 table, tensor glyph/ellipsoid/principal
  axes, decomposition panel.
- `static.efg`: T2 glyph/table, magnitude color, principal-axis display.
- `static.per-class`: stacked table/bar by ring type, bond category, shell, or
  coupling type.
- `static.geometry`: overlay measurement, bond/ring/residue labels, distance
  shells.
- `static.topology`: topology labels, graph/bond/ring membership overlay.
- `static.system`: dashboard readout/table; no atom anchor required.
- `static.rollup`: statistics table, min/max frame markers, histogram.
- `static.event`: event list, timeline marks, jump-to-frame action.
- `static.embedding`: vector table/projection/nearest-neighbor display; not a
  default strip for all 256 dimensions.

## Import Families and Display Modes

The following table is the main inventory. "NPY fields" are the SDK/import
surface. "H5 overlap/products" lists dense trajectory representations where
present. H5 entries are not forced to unify with NPY entries; overlap means the
consumer design should know they are the same or related measurement family.

| Family | Axis/anchor | NPY fields | H5 overlap/products | Temporal strip modes | Static/current-frame modes |
| --- | --- | --- | --- | --- | --- |
| Identity / positions | atom | `pos`, `element`, `residue_index`, `residue_type`, `atoms_category_info` | `positions`; `/atoms` identity | `strip.vector.component`, `strip.vector.magnitude` for position/displacement; derived geometry strips | `static.topology`, `static.geometry`, atom labels, color by element/residue/category |
| Topology sidecar | residue, bond, ring, ring_membership | `residues`, `bonds`, `rings`, `ring_membership`, `ring_geometry` | topology loaded through sidecar/protein; ring pucker/neighbourhood use related ring axes | generally static; possible `strip.count` for selected neighbourhood over time via derived sources | `static.topology`, `static.geometry`, ring/bond overlays, residue/atom picker context |
| Ring contributions | ring_contribution_pair | `ring_contributions` | no direct dense TS listed | possible `strip.tensor.*` per selected pair if exposed; likely advanced | `static.tensor`, `static.per-class`, ring-pair table/heatmap |
| Biot-Savart ring current | atom | `bs_shielding`, `bs_per_type_T0`, `bs_per_type_T2`, `bs_total_B`, `bs_ring_counts` | `bs_shielding_time_series`, `bs_welford`, `bs_t0_autocorrelation`, ring neighbourhood stats overlap | `strip.tensor.T0/T1/T2`, `strip.vector.*` for B, `strip.count`, `strip.per-class`, `strip.rollup`, `strip.spectrum` | `static.tensor`, `static.vector`, `static.per-class`, atom heatmap, ring count table, autocorr panel |
| Haigh-Mallion ring current | atom | `hm_shielding`, `hm_per_type_T0`, `hm_per_type_T2` | `hm_shielding_time_series`, `hm_welford` | `strip.tensor.T0/T1/T2`, `strip.per-class`, `strip.rollup`, `strip.spectrum` | `static.tensor`, `static.per-class`, atom heatmap |
| Ring susceptibility | atom | `ringchi_shielding` | `ringchi_shielding_time_series` | `strip.tensor.T0/T1/T2`, `strip.spectrum` | `static.tensor`, atom heatmap |
| Pi quadrupole | atom | `pq_shielding`, `pq_per_type_T0`, `pq_per_type_T2` | `piquad_shielding_time_series` | `strip.tensor.T0/T1/T2`, `strip.per-class`, `strip.spectrum` | `static.tensor`, `static.per-class`, atom heatmap |
| Dispersion | atom | `disp_shielding`, `disp_per_type_T0`, `disp_per_type_T2` | `disp_shielding_time_series` | `strip.tensor.T0/T1/T2`, `strip.per-class`, `strip.spectrum` | `static.tensor`, `static.per-class`, atom heatmap |
| McConnell bond anisotropy | atom | `mc_shielding`, `mc_category_T2`, `mc_scalars` | `mc_shielding_time_series`, `mc_welford` | `strip.tensor.T0/T1/T2`, `strip.per-class` by McConnell category, scalar distance/sum strips, `strip.rollup` | `static.tensor`, `static.per-class`, nearest-bond scalar table |
| Coulomb electrostatics | atom | `coulomb_shielding`, `coulomb_E`, `coulomb_efg_backbone`, `coulomb_efg_aromatic`, `coulomb_scalars` | mostly per-frame NPY; dense H5 uses APBS/AIMNet2/MOPAC variants instead | `strip.tensor.*`, `strip.vector.*`, `strip.efg/T2`, scalar magnitude/projection strips | `static.tensor`, `static.vector`, `static.efg`, electrostatic table |
| APBS electrostatics | atom | `apbs_E`, `apbs_efg` | `apbs_efield_time_series`, `apbs_efg_time_series` | `strip.vector.component/magnitude`, `strip.tensor.T2`, `strip.spectrum` | `static.vector`, `static.efg`, field/glyph overlay |
| AIMNet2 | atom | `aimnet2_charges`, `aimnet2_aim`, `aimnet2_efg`, `aimnet2_efg_aromatic`, `aimnet2_efg_backbone`, `aimnet2_charge_response_gradient`, `aimnet2_charge_response_gradient_scalar` | `aimnet2_charge_time_series`, `aimnet2_embedding_time_series`, `aimnet2_charge_response_gradient_time_series`, `aimnet2_charge_response_gradient_welford` | charge scalar, CRG vector component/magnitude, EFG T2/magnitude if available per time source, embedding projections/components only by explicit request, rollup stats | charge heatmap, `static.efg`, CRG vector glyph, embedding projection/table, nearest-neighbor/projection view |
| EEQ | atom | `eeq_charges`, `eeq_cn` | `eeq_welford` | if frame-local NPY sampled: charge/CN scalar strips; H5 rollup as `strip.rollup` only | charge/CN heatmap, stats table |
| MOPAC core | atom, bond, protein | `mopac_charges`, `mopac_scalars`, `mopac_bond_orders`, `mopac_global` | `mopac_charge_welford`, `mopac_bond_order_welford` | atom charge/scalar strips if frame source present; bond order strips by selected bond; global scalar strips; rollup stats | atom/bond heatmaps, bond labels, global table |
| MOPAC Coulomb | atom | `mopac_coulomb_shielding`, `mopac_coulomb_E`, `mopac_coulomb_efg_backbone`, `mopac_coulomb_efg_aromatic`, `mopac_coulomb_scalars` | `mopac_coulomb_shielding_time_series` is T2-only dense product | `strip.tensor.T2`, vector/scalar strips where NPY source present, category comparisons | `static.tensor`, `static.vector`, `static.efg`, MOPAC-vs-force-field comparison |
| MOPAC McConnell | atom | `mopac_mc_shielding`, `mopac_mc_category_T2`, `mopac_mc_scalars` | `mopac_mc_shielding_time_series` | `strip.tensor.*`, category T2 strips, scalar distance/sum strips | `static.tensor`, `static.per-class`, comparison table |
| H-bond kernel | atom | `hbond_shielding`, `hbond_scalars` | `hbond_shielding_time_series`, `hbond_count_welford` | tensor T0/T1/T2, nearest distance, inverse distance, count, McConnell scalar, rollups | `static.tensor`, scalar table, H-bond participant overlay |
| Larsen H-bond grid | atom | `larsen_hbond_shielding`, `larsen_hbond_1pHB_shielding`, `larsen_hbond_2pHB_shielding`, `larsen_hbond_1pHaB_shielding`, `larsen_hbond_2pHaB_shielding`, `larsen_hbond_diagnostic_CB_shielding`, `larsen_hbond_water_term`, `larsen_hbond_count` | corresponding four shielding time series plus count/water term time series | tensor T0/T1/T2 per component/source, count strip, water term strip, source comparison strips | `static.tensor`, per-component H-bond grid table, count/water heatmap |
| DSSP / secondary structure | atom/residue | `dssp_backbone`, `dssp_ss8`, `dssp_hbond_energy`, `dssp_chi` | `dssp8_time_series`, `dssp8_transition` | `strip.category` SS8, phi/psi-like scalars from backbone, H-bond energy strips, transition/event strips | residue color bands, SS labels, H-bond partner table, transition matrix |
| Residue dihedral geometry | residue | `omega_actual`, `omega_deviation`, `omega_is_xpro`, `dssp_chi`; plus runtime geometry | `dihedral_time_series`, `dihedral_bin_transition` | phi/psi/omega/chi strips, omega deviation, category/bin strips, transition counts, FFT/spectrum | Ramachandran/rotamer static view, residue labels, current angle table |
| Planar/ring geometry | atom, residue, aromatic_ring, saturated_ring | `pyramidalization`, `aromatic_chi2`, `pucker_Q`, `pucker_theta` | `ring_pucker_time_series` | scalar angle/amplitude strips by atom/residue/ring, category/event strips for puckers | ring geometry overlay, pucker/chi table, atom heatmap |
| User-selected geometry | atom tuple | derived, not SDK NPY | positions-backed runtime measurement | distance, angle, dihedral strips; FFT/spectrum | measurement overlay, tuple labels, reveal/jump target |
| SASA | atom | `atom_sasa`, `sasa_normal` | `sasa_time_series`, `sasa_welford` | SASA scalar, normal vector components/magnitude if source exists, rollups | atom SASA heatmap, normal glyph overlay, stats table |
| Water field | atom | `water_efield`, `water_efield_first`, `water_efg`, `water_efg_first`, `water_shell_counts` | `water_field_time_series`, `water_field_welford` | E-field vector strips, EFG T2 strips, first/second shell count strips, rollups | field glyphs, EFG glyph/table, hydration shell count heatmap |
| Hydration shell/geometry | atom | `hydration_shell`, `water_polarization` | `hydration_shell_time_series`, `hydration_geometry_time_series`, `hydration_shell_welford`, `hydration_geometry_welford` | half-shell asymmetry, dipole alignment/coherence, ion charge/distance, first-shell count, vector strips, rollups | hydration overlays, dipole/surface normal glyphs, multi-channel table |
| Gromacs energy/runtime | protein/system | `gromacs_energy` | `gromacs_energy_time_series` | `strip.system` for energy, temperature, pressure, density, box, bonded/nonbonded terms; tensor component strips for pressure tensor/virial | system dashboard table, tensor table, timeline selection |
| Bonded energy | atom | `bonded_energy` | `bonded_energy_time_series` currently total only | atom scalar strip for total and NPY split terms when available | atom heatmap, per-term table |
| Bond length stats | bond | no direct NPY field except topology bonds | `bond_length_stats` | `strip.rollup` only unless frame-time bond length derived from positions | bond overlay, min/max/stat table |
| RMSD tracking | system/subset | no NPY catalog field | `rmsd_tracking` | system RMSD strip, event/rug for spikes | RMSD readout, subset highlight, jump-to-frame |
| Ring neighbourhood stats | atom/ring slot | overlaps ring proximity concepts, not direct SDK field | `ring_neighbourhood_trajectory_stats` | distance/rho/z/in-plane-angle strips by selected atom and slot, ring count derived strips | ring neighbourhood overlay, nearest-ring table |
| J-couplings | residue | not in current NPY catalog | `j_coupling_time_series` | per-coupling scalar strips, existence-mask-aware gaps | coupling table by residue, existence labels |
| ORCA DFT shielding | atom | `orca_total`, `orca_diamagnetic`, `orca_paramagnetic` | local ORCA reader via `DftShieldingStore`; may also exist as NPY fields | tensor T0/T1/T2 strips for total/dia/para, component comparison strips, FFT/spectrum | `static.tensor`, total/dia/para table, atom heatmap |
| Tripeptide reference | atom | `tripeptide_bb_shielding`, `tripeptide_bb_residual_vec`, `tripeptide_bb_match_distance`, `tripeptide_bb_method_tag`, `tripeptide_neighbor_shielding`, `tripeptide_neighbor_residual_vec_prev`, `tripeptide_neighbor_residual_vec_next` | corresponding shielding/residual/method tag time series | tensor strips, residual vector strips, match distance strip, method category strip | reference tensor table, residual vector overlay, method label |
| Mutation delta | mutation_match_pair | `delta_shielding`, `delta_scalars`, `delta_apbs`, `delta_ring_proximity`, wild-type/mutant dia/para and delta dia/para shieldings | no dense H5 trajectory reader listed | mostly static/comparison; temporal only if mutation frames are sequenced as a run | delta tensor table, APBS/ring proximity comparison, pair heatmap |
| Selections/events | event/system | no SDK NPY field | `/trajectory/selections` | event/rug strip, jump-to-frame marks, selection counts | event table, reveal/jump actions |

## Complete SDK / NPY Groups

Notation: `field[axis, columns, wrapper, units]`; `?` means optional in the SDK
catalog; `var` means catalog `cols=None`.

- `identity`:
  `pos[atom,3,VectorField,A]`,
  `element[atom,var,np.ndarray,-]`,
  `residue_index[atom,var,np.ndarray,-]`,
  `residue_type[atom,var,np.ndarray,-]`,
  `ring_contributions[ring_contribution_pair,58,RingContributions,-]`,
  `ring_geometry[aromatic_ring,10,RingGeometry,A]`,
  `atoms_category_info[atom,var?,np.ndarray,-]`.
- `biot_savart`:
  `bs_shielding[atom,9,ShieldingTensor,ppm_T_per_nA]`,
  `bs_per_type_T0[atom,8,PerRingTypeT0,ppm_T_per_nA]`,
  `bs_per_type_T2[atom,40,PerRingTypeT2,ppm_T_per_nA]`,
  `bs_total_B[atom,3,VectorField,T]`,
  `bs_ring_counts[atom,4,RingCounts,-]`.
- `haigh_mallion`:
  `hm_shielding`, `hm_per_type_T0`, `hm_per_type_T2`.
- `pi_quadrupole`:
  `pq_shielding`, `pq_per_type_T0`, `pq_per_type_T2`.
- `dispersion`:
  `disp_shielding`, `disp_per_type_T0`, `disp_per_type_T2`.
- `ring_susceptibility`:
  `ringchi_shielding`.
- `mcconnell`:
  `mc_shielding`, `mc_category_T2`, `mc_scalars`.
- `coulomb`:
  `coulomb_shielding`, `coulomb_E`, `coulomb_efg_backbone`,
  `coulomb_efg_aromatic`, `coulomb_scalars`.
- `hbond`:
  `hbond_shielding`, `hbond_scalars`.
- `dssp`:
  `dssp_backbone`, `dssp_ss8`, `dssp_hbond_energy`, `dssp_chi`.
- `sasa`:
  `atom_sasa`, `sasa_normal`.
- `water_field`:
  `water_efield`, `water_efield_first`, `water_efg`,
  `water_efg_first`, `water_shell_counts`.
- `hydration`:
  `hydration_shell`.
- `water_polarization`:
  `water_polarization`.
- `eeq`:
  `eeq_charges`, `eeq_cn`.
- `gromacs`:
  `gromacs_energy` (SDK 43 cols; C++ generated catalog currently 42).
- `bonded`:
  `bonded_energy`.
- `mopac_core`:
  `mopac_charges`, `mopac_scalars`, `mopac_bond_orders`,
  `mopac_global`.
- `mopac_coulomb`:
  `mopac_coulomb_shielding`, `mopac_coulomb_E`,
  `mopac_coulomb_efg_backbone`, `mopac_coulomb_efg_aromatic`,
  `mopac_coulomb_scalars`.
- `mopac_mcconnell`:
  `mopac_mc_shielding`, `mopac_mc_category_T2`, `mopac_mc_scalars`.
- `apbs`:
  `apbs_E`, `apbs_efg`.
- `orca`:
  `orca_total`, `orca_diamagnetic`, `orca_paramagnetic`.
- `delta`:
  `delta_shielding`, `delta_scalars`, `delta_apbs`,
  `delta_ring_proximity`, `wt_shielding_diamagnetic`,
  `wt_shielding_paramagnetic`, `mut_shielding_diamagnetic`,
  `mut_shielding_paramagnetic`, `delta_shielding_diamagnetic`,
  `delta_shielding_paramagnetic`.
- `aimnet2`:
  `aimnet2_charges`, `aimnet2_aim`, `aimnet2_efg`,
  `aimnet2_efg_aromatic`, `aimnet2_efg_backbone`,
  `aimnet2_charge_response_gradient`,
  `aimnet2_charge_response_gradient_scalar`.
- `planar_geometry`:
  `pyramidalization`, `omega_actual`, `omega_deviation`,
  `aromatic_chi2`, `pucker_Q`, `pucker_theta`, `omega_is_xpro`.
- `tripeptide`:
  `tripeptide_bb_shielding`, `tripeptide_bb_residual_vec`,
  `tripeptide_bb_match_distance`, `tripeptide_bb_method_tag`,
  `tripeptide_neighbor_shielding`,
  `tripeptide_neighbor_residual_vec_prev`,
  `tripeptide_neighbor_residual_vec_next`.
- `larsen_hbond`:
  `larsen_hbond_shielding`, `larsen_hbond_1pHB_shielding`,
  `larsen_hbond_2pHB_shielding`,
  `larsen_hbond_1pHaB_shielding`,
  `larsen_hbond_2pHaB_shielding`,
  `larsen_hbond_diagnostic_CB_shielding`,
  `larsen_hbond_water_term`, `larsen_hbond_count`.
- `topology`:
  `residues`, `bonds`, `rings`, `ring_membership`.

## Dense H5 Product Groups

These are import products currently decoded by `QtTrajectoryH5`.

- Frame identity/root: `/atoms`, `/trajectory/frames`, trajectory attrs.
- Atom time series:
  `/trajectory/positions`,
  `/trajectory/bs_shielding_time_series`,
  `/trajectory/hm_shielding_time_series`,
  `/trajectory/mc_shielding_time_series`,
  `/trajectory/piquad_shielding_time_series`,
  `/trajectory/ringchi_shielding_time_series`,
  `/trajectory/disp_shielding_time_series`,
  `/trajectory/hbond_shielding_time_series`,
  `/trajectory/mopac_coulomb_shielding_time_series`,
  `/trajectory/mopac_mc_shielding_time_series`,
  `/trajectory/tripeptide_bb_shielding_time_series`,
  `/trajectory/tripeptide_neighbor_shielding_time_series`,
  `/trajectory/larsen_hbond_1pHB_shielding_time_series`,
  `/trajectory/larsen_hbond_1pHaB_shielding_time_series`,
  `/trajectory/larsen_hbond_2pHB_shielding_time_series`,
  `/trajectory/larsen_hbond_2pHaB_shielding_time_series`.
- Atom scalar/vector/composite time series:
  `/trajectory/mopac_vs_ff14sb_reconciliation`,
  `/trajectory/sasa_time_series`,
  `/trajectory/aimnet2_charge_time_series`,
  `/trajectory/larsen_hbond_count_time_series`,
  `/trajectory/larsen_hbond_water_term_time_series`,
  `/trajectory/bonded_energy_time_series`,
  `/trajectory/water_field_time_series`,
  `/trajectory/hydration_shell_time_series`,
  `/trajectory/hydration_geometry_time_series`,
  `/trajectory/apbs_efield_time_series`,
  `/trajectory/apbs_efg_time_series`,
  `/trajectory/aimnet2_embedding_time_series`,
  `/trajectory/aimnet2_charge_response_gradient_time_series`,
  `/trajectory/tripeptide_bb_residual_vec_time_series`,
  `/trajectory/tripeptide_neighbor_residual_vec_prev_time_series`,
  `/trajectory/tripeptide_neighbor_residual_vec_next_time_series`,
  `/trajectory/tripeptide_bb_method_tag_time_series`.
- Residue/ring/system/special time series:
  `/trajectory/dihedral_time_series`,
  `/trajectory/dssp8_time_series`,
  `/trajectory/j_coupling_time_series`,
  `/trajectory/ring_pucker_time_series`,
  `/trajectory/ring_neighbourhood_trajectory_stats`,
  `/trajectory/gromacs_energy_time_series`,
  `/trajectory/rmsd_tracking`.
- Statistical/static products:
  `/trajectory/bond_length_stats`,
  `/trajectory/bs_welford`,
  `/trajectory/hm_welford`,
  `/trajectory/mc_welford`,
  `/trajectory/sasa_welford`,
  `/trajectory/eeq_welford`,
  `/trajectory/hbond_count_welford`,
  `/trajectory/mopac_charge_welford`,
  `/trajectory/mopac_bond_order_welford`,
  `/trajectory/water_field_welford`,
  `/trajectory/aimnet2_charge_response_gradient_welford`,
  `/trajectory/hydration_shell_welford`,
  `/trajectory/hydration_geometry_welford`,
  `/trajectory/bs_t0_autocorrelation`,
  `/trajectory/dssp8_transition`,
  `/trajectory/dihedral_bin_transition`,
  `/trajectory/selections`.

## Immediate Design Pressure

The class-tree pass should not start with the current UI. It should start from
these display primitives:

- atom scalar calculation
- atom vector calculation
- atom tensor calculation
- atom per-class calculation
- residue scalar/category calculation
- ring/bond calculation
- system calculation
- rollup/statistical calculation
- event/selection calculation
- derived geometry calculation
- static overlay/panel calculation

Then each import family above maps onto those bases. That is the point of this
inventory: the picker and dashboard should be driven by this map, not by the
first few strips that happened to work.

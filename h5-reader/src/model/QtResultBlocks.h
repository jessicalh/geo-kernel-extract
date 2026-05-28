// QtResultBlocks.h — small typed value blocks shared by the per-frame
// result-group views (QtBiotSavartGroup, QtMcConnellGroup, ...).
//
// These mirror the SDK's _tensors.py wrappers. This header holds the
// ring-current blocks (PerRingTypeT0, PerRingTypeT2, RingCounts);
// sibling blocks (PerBondCategoryT2, the scalar blocks, RingContributions,
// ...) arrive with the groups that use them. "Objects answer questions
// about themselves": a block exposes the same semantic accessors the SDK
// group fields offer (forType(RingTypeIndex), total(), ...) over the
// raw NPY columns, so consuming code never indexes a bare double[8].
//
// Library/contract provenance: the per-ring-type axis is ordered by
// RingTypeIndex (8 aromatic types; ProPyrrolidine excluded — saturated,
// no ring current). Per GEOMETRIC_KERNEL_CATALOGUE + _catalog.py
// (bs_per_type_T0 cols=8, bs_per_type_T2 cols=40 = 8x5).

#pragma once

#include "Types.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace h5reader::model {

// bs/hm/pq/disp _per_type_T0 — one isotropic (T0) scalar per aromatic
// ring type. 8 columns, ordered by RingTypeIndex (0..7).
struct PerRingTypeT0 {
    std::array<double, kAromaticRingTypeCount> byType = {};

    double forType(RingTypeIndex t) const { return byType[static_cast<std::size_t>(t)]; }
    double total() const {
        double s = 0.0;
        for (double v : byType)
            s += v;
        return s;
    }
};

// bs/hm/pq/disp _per_type_T2 — a 5-component T2 per aromatic ring type.
// 40 columns = 8 types x 5 (m = -2..+2), ordered by RingTypeIndex.
struct PerRingTypeT2 {
    std::array<std::array<double, 5>, kAromaticRingTypeCount> byType = {};

    const std::array<double, 5>& forType(RingTypeIndex t) const {
        return byType[static_cast<std::size_t>(t)];
    }
    // Sum over types — reconstructs the summed-shielding T2 these decompose.
    std::array<double, 5> total() const {
        std::array<double, 5> s = {};
        for (const auto& comp : byType)
            for (std::size_t i = 0; i < 5; ++i)
                s[i] += comp[i];
        return s;
    }
};

// bs_ring_counts — aromatic-ring-centre counts within concentric shells.
// 4 columns: 3 / 5 / 8 / 12 Å (PhysicalConstants RING_COUNT_SHELL).
struct RingCounts {
    double within3A = 0.0;
    double within5A = 0.0;
    double within8A = 0.0;
    double within12A = 0.0;
};

// ── Bond-anisotropy (McConnell) category decomposition ─────────────────
// The mc_category_T2 axis is a 5-value DECOMPOSITION enum — NOT the 8-value
// topology BondCategory in Types.h. Mirror of the SDK _types.py BondCategory:
// the three "total" bins plus the two single-nearest-bond tensors.
enum class McConnellCategory : std::int8_t {
    BackboneTotal = 0,
    SidechainTotal = 1,
    AromaticTotal = 2,
    CONearest = 3,
    CNNearest = 4,
};
constexpr int kMcConnellCategoryCount = 5;

// mc_category_T2 / mopac_mc_category_T2 — 5-component T2 per McConnell
// category (5 x 5 = 25, row-major).
struct PerBondCategoryT2 {
    std::array<std::array<double, 5>, kMcConnellCategoryCount> byCategory = {};
    const std::array<double, 5>& forCategory(McConnellCategory c) const {
        return byCategory[static_cast<std::size_t>(c)];
    }
};

// Electric field gradient — symmetric-traceless, T2 ONLY (5 components,
// m = -2..+2). T0 (trace) and T1 (antisymmetric) are structural zeros after
// the calculator's traceless projection, so this is a 5-component block, NOT
// a SphericalTensor (mirrors the SDK EFGTensor; e3nn irreps "1x2e"). Used by
// the Coulomb / APBS / water / AIMNet2 / MOPAC-Coulomb EFG fields.
struct QtEfg {
    std::array<double, 5> t2 = {};
    double t2Magnitude() const {
        double s = 0.0;
        for (double v : t2)
            s += v * v;
        return std::sqrt(s);
    }
};

// mc_scalars (6) — McConnell summary: the (3cos^2 θ − 1)/r^3 angular sums per
// category (before Δχ) + the nearest-bond distances (Å). Shared by both the
// ff14SB QtMcConnellGroup and QtMopacMcConnellGroup.
//
// nearest_CO_dist / nearest_CN_dist carry the NO_DATA_SENTINEL = 99.0 (Å) when
// no C=O / C–N bond lies within the McConnell cutoff (writers McConnellResult.cpp
// :150-151,255-256 + MopacMcConnellResult.cpp:141-142,244-245 — initialised to
// the sentinel, written through unchanged on a miss; the writer's own "is this
// real" test is `dist < NO_DATA_SENTINEL`). It is NOT a real distance — guard
// with hasNearestCO() / hasNearestCN() before display / ML use. (No atom hits it
// in the 1P9J frame-0 fixture — every atom finds both partners, max 8.93/8.05 Å
// — but isolated / chain-terminal atoms in other frames will.)
struct McConnellScalars {
    static constexpr double kNoDataSentinel = 99.0;  // PhysicalConstants NO_DATA_SENTINEL
    double co_sum = 0.0;
    double cn_sum = 0.0;
    double sidechain_sum = 0.0;
    double aromatic_sum = 0.0;
    double nearest_CO_dist = 0.0;
    double nearest_CN_dist = 0.0;
    bool hasNearestCO() const { return nearest_CO_dist < kNoDataSentinel; }
    bool hasNearestCN() const { return nearest_CN_dist < kNoDataSentinel; }
};

// coulomb_scalars (4) — Coulomb E-field summary (V/Å).
struct CoulombScalars {
    double E_magnitude = 0.0;          // |E_total|
    double E_bond_proj = 0.0;          // E along the primary bond axis (Buckingham σ_iso input)
    double E_backbone_frac = 0.0;      // SIGNED projection E_backbone·Ê_total (V/Å): +aligned/−opposed, |·|≤|E_backbone|; historical name, NOT a 0–1 fraction (CoulombResult.cpp:287)
    double aromatic_E_magnitude = 0.0;
};

// hbond_scalars (4) — H-bond summary.
struct HBondScalars {
    double nearest_dist = 0.0;      // Å to nearest partner
    double inv_d3 = 0.0;            // 1/d^3, d = atom-to-(H-bond-midpoint) distance (SDK/writer name: inv_d3)
    double count_3_5A = 0.0;        // candidates within 3.5 Å
    double mcconnell_scalar = 0.0;  // Σ (3cos^2 θ − 1)/r^3 over contributing H-bonds
};

// ── Shared row -> typed-value decoders (one place; used by every group) ──
// 9-col row -> SphericalTensor, layout [T0, T1[3], T2[5]]. See the T1-basis
// note in Types.h: on the per-frame NPY path the T1 is Cartesian.
inline SphericalTensor UnpackSphericalTensor(const double* r) {
    SphericalTensor st;
    st.T0 = r[0];
    for (std::size_t i = 0; i < 3; ++i)
        st.T1[i] = r[1 + i];
    for (std::size_t i = 0; i < 5; ++i)
        st.T2[i] = r[4 + i];
    return st;
}

// 5-col row -> QtEfg (T2 only).
inline QtEfg UnpackEfg(const double* r) {
    QtEfg e;
    for (std::size_t i = 0; i < 5; ++i)
        e.t2[i] = r[i];
    return e;
}

// ── Scalar / vector family blocks ──────────────────────────────────────
// Mirrors of the SDK's scalar/vector outputs. Several of these arrays have
// NO Python wrapper in _tensors.py (delivered as bare np.ndarray); their
// column semantics live ONLY in the C++ writer's WriteFeatures, so each block
// is grounded against that definitive source (file:line cited), not the
// catalog string. Each exposes a `FromRow(const double*)` decode so the group
// views are uniform one-liners and the column order lives with the type.

// dssp_backbone (N×5) — per-atom, BROADCAST from the atom's residue (every
// atom of a residue shares its residue's DSSP row). Writer DsspResult.cpp:
// 214-224 (libdssp; Joosten 2011 / Kabsch & Sander 1983).
//   phi, psi : backbone dihedrals, RADIANS, in libdssp's NEGATED-IUPAC sign
//              convention (DsspResult.h:22-25) — negate for IUPAC φ/ψ.
//   sasa     : DSSP per-RESIDUE SASA (Å²). DISTINCT from atom_sasa (the
//              Shrake-Rupley per-ATOM SASA in QtSasaGroup) — never conflate.
//   ssHelix  : 1.0 if SS ∈ {H,G,I} (any helix) else 0.0  } binary collapse of
//   ssSheet  : 1.0 if SS ∈ {E,B} (strand/bridge) else 0.0 } the 8-class (full
//              8-class is dssp_ss8 / DsspSs8).
// Residues DSSP did not map are written all-zero (looks like coil, φ=ψ=0); the
// NPY does not carry DsspResidue.observed, so default-coil is indistinguishable
// from real coil at this layer.
struct DsspScalars {
    double phi = 0.0;
    double psi = 0.0;
    double sasa = 0.0;
    double ssHelix = 0.0;
    double ssSheet = 0.0;
    static DsspScalars FromRow(const double* r) {
        return DsspScalars{r[0], r[1], r[2], r[3], r[4]};
    }
};

// dssp_ss8 (N×8) — per-atom one-hot of the DSSP 8-class SS, broadcast from
// residue. Column index == DsspCode ordinal (DsspResult.cpp:234-255):
// H=0,G=1,I=2,E=3,B=4,T=5,S=6,C=7 — identical to Types.h DsspCode, so the hot
// column decodes straight to a DsspCode. Unmapped residue → coil (col 7).
struct DsspSs8 {
    std::array<double, 8> oneHot = {};
    static DsspSs8 FromRow(const double* r) {
        DsspSs8 s;
        for (std::size_t i = 0; i < 8; ++i)
            s.oneHot[i] = r[i];
        return s;
    }
    // The set class (argmax). One-hot, so this is exact; a degenerate all-zero
    // row (shouldn't occur) resolves to AlphaHelix(0).
    DsspCode dominant() const {
        std::size_t best = 0;
        for (std::size_t i = 1; i < 8; ++i)
            if (oneHot[i] > oneHot[best])
                best = i;
        return static_cast<DsspCode>(best);
    }
    bool is(DsspCode c) const {
        auto i = static_cast<std::size_t>(c);
        return i < 8 && oneHot[i] != 0.0;
    }
};

// dssp_hbond_energy (N×4) — per-atom (broadcast from residue), kcal/mol: the
// four DSSP electrostatic H-bond energies, in writer order (DsspResult.cpp:
// 259-274), each sourced from libdssp's per-residue slots (DsspResult.cpp:144-163):
//   acceptor0/1 ← DsspResidue.acceptors[0..1] ← libdssp acceptor(0/1)
//   donor0/1    ← DsspResidue.donors[0..1]    ← libdssp donor(0/1)
// 0.0 = empty slot (no partner). NOTE: dssp.hpp (the libdssp header) does NOT
// document the proton direction of acceptor()/donor() — the classic DSSP
// footgun. The library reads `acceptors` as "a partner C=O accepts this
// residue's N-H, i.e. this residue donates" (DsspResult.h:57-58); treat these
// as libdssp slot labels, NOT an asserted donor/acceptor role, until pinned
// against libdssp's source.
struct DsspHBondEnergy {
    double acceptor0 = 0.0;
    double acceptor1 = 0.0;
    double donor0 = 0.0;
    double donor1 = 0.0;
    static DsspHBondEnergy FromRow(const double* r) {
        return DsspHBondEnergy{r[0], r[1], r[2], r[3]};
    }
};

// dssp_chi (N×12) — per-atom (broadcast from residue) sidechain dihedrals.
// 4 chi angles × 3 cols, INTERLEAVED per angle (DsspResult.cpp:283-321):
//   base = k*3 → [cos χ_{k+1}, sin χ_{k+1}, exists]. (cos,sin) avoids the ±π
// wrap discontinuity. exists==0.0 ⇒ the residue has no χ_{k+1} (GLY/ALA, or k
// beyond its chi count) and cos=sin=0.0 there. NOTE: the writer's own inline
// comment says "NaN" for missing, but the code writes 0.0 (data init 0.0,
// else-branch leaves it) — writer code is definitive over its own comment.
struct DsspChi {
    std::array<double, 12> data = {};
    struct Angle {
        double cos;
        double sin;
        bool exists;
    };
    static DsspChi FromRow(const double* r) {
        DsspChi c;
        for (std::size_t i = 0; i < 12; ++i)
            c.data[i] = r[i];
        return c;
    }
    // k = 0..3 → χ1..χ4.
    Angle forAngle(int k) const {
        const double* a = &data[static_cast<std::size_t>(k) * 3];
        return Angle{a[0], a[1], a[2] != 0.0};
    }
};

// water_shell_counts (N×2) — explicit-solvent water-oxygen counts around the
// atom (WaterFieldResult.cpp:255; shells from WaterFieldResult.h):
//   nFirst  : water O within the first shell (< 3.5 Å)
//   nSecond : water O in the second shell (3.5–5.5 Å)
// Stored as float64 in the NPY (hence double, not int).
struct WaterShellCounts {
    double nFirst = 0.0;
    double nSecond = 0.0;
    static WaterShellCounts FromRow(const double* r) {
        return WaterShellCounts{r[0], r[1]};
    }
};

// hydration_shell (N×4) — explicit-water packing geometry (COM-reference
// variant; HydrationShellResult.cpp:139-142):
//   halfShellAsymmetry : fraction of first-shell waters on the solvent-exposed
//                        side (UCSB half-shell method), 0..1
//   meanWaterDipoleCos : mean cos(water dipole, atom→water) — orientation order
//   nearestIonDist (Å) : distance to nearest ion; +INFINITY when no ion within
//                        the search cutoff (HydrationShellResult.cpp:114-115) —
//                        callers MUST guard (use hasNearestIon()) before
//                        display / ML use.
//   nearestIonCharge   : that ion's charge (e); 0.0 when none in cutoff.
struct HydrationShell {
    double halfShellAsymmetry = 0.0;
    double meanWaterDipoleCos = 0.0;
    double nearestIonDist = 0.0;
    double nearestIonCharge = 0.0;
    bool hasNearestIon() const { return std::isfinite(nearestIonDist); }
    static HydrationShell FromRow(const double* r) {
        return HydrationShell{r[0], r[1], r[2], r[3]};
    }
};

// water_polarization (N×10) — first-shell water orientation in the SASA-normal
// frame (HydrationGeometryResult.cpp:141-150; depends on SasaResult):
//   dipole   (Vec3) : net first-shell water dipole vector (Σ dᵢ)
//   normal   (Vec3) : SASA outward surface normal (reference frame; the SAME
//                     vector as QtSasaGroup::normal, duplicated by design)
//   asymmetry       : half-shell asymmetry in the SASA-normal frame, 0..1
//   alignment       : cos(net dipole, surface normal)
//   coherence       : |Σ dᵢ| / n — ordered (→|d|) vs random (→0)
//   shellCount      : first-shell water count (statistical weight)
struct WaterPolarization {
    Vec3 dipole = Vec3::Zero();
    Vec3 normal = Vec3::Zero();
    double asymmetry = 0.0;
    double alignment = 0.0;
    double coherence = 0.0;
    double shellCount = 0.0;
    static WaterPolarization FromRow(const double* r) {
        WaterPolarization w;
        w.dipole = Vec3(r[0], r[1], r[2]);
        w.normal = Vec3(r[3], r[4], r[5]);
        w.asymmetry = r[6];
        w.alignment = r[7];
        w.coherence = r[8];
        w.shellCount = r[9];
        return w;
    }
};

// bonded_energy (N×7) — per-atom GROMACS bonded-energy decomposition, kJ/mol;
// each interaction's energy is split evenly among its participating atoms
// (BondedEnergyResult.cpp:257-264). Columns:
//   bond, angle, ureyBradley, proper, improper, cmap, total(=Σ of the six).
// The columns are force-field-AGNOSTIC: their values come from whatever bonded
// terms the run's TPR defines. ureyBradley + cmap are zero for force fields
// lacking those terms (AMBER ff14SB carries neither). (BondedEnergyResult.h's
// "CHARMM36m" header line predates the AMBER-only move and is not propagated.)
struct BondedEnergy {
    double bond = 0.0;
    double angle = 0.0;
    double ureyBradley = 0.0;
    double proper = 0.0;
    double improper = 0.0;
    double cmap = 0.0;
    double total = 0.0;
    static BondedEnergy FromRow(const double* r) {
        return BondedEnergy{r[0], r[1], r[2], r[3], r[4], r[5], r[6]};
    }
};

// gromacs_energy (1×43, PROTEIN-axis — one row for the whole frame, NOT
// per-atom) — the GROMACS .edr energy terms for this frame, GROMACS-native
// units. Column order per the writer GromacsEnergyResult.cpp:29-53 (= the
// GromacsEnergy struct in GromacsEnergyResult.h, EXCLUDING its time_ps field).
//
// SCHEMA NOTE (writer-definitive): the writer emits 43 columns; _catalog.py and
// the generated QtFieldCatalog.gen.h declare cols=42 — an off-by-one CONFIRMED
// against the fixture bytes (shape (1,43)). This block decodes 43. Flagged for
// the library team (_catalog.py gromacs_energy 42→43); the loader must take the
// actual NPY shape as truth over the catalog cols and log the mismatch.
struct GromacsEnergy {
    std::array<double, 43> raw = {};
    static GromacsEnergy FromRow(const double* r) {
        GromacsEnergy g;
        for (std::size_t i = 0; i < 43; ++i)
            g.raw[i] = r[i];
        return g;
    }
    // Electrostatic (kJ/mol)
    double coulombShortRange() const { return raw[0]; }   // PME real-space
    double coulombReciprocal() const { return raw[1]; }   // PME reciprocal-space
    double coulomb14() const { return raw[2]; }            // 1-4 intramolecular
    // Bonded (kJ/mol)
    double bond() const { return raw[3]; }
    double angle() const { return raw[4]; }
    double ureyBradley() const { return raw[5]; }
    double properDih() const { return raw[6]; }
    double improperDih() const { return raw[7]; }
    double cmapDih() const { return raw[8]; }
    // Van der Waals (kJ/mol)
    double ljShortRange() const { return raw[9]; }
    double lj14() const { return raw[10]; }
    double dispersionCorrection() const { return raw[11]; }
    // Thermodynamic state
    double potential() const { return raw[12]; }    // kJ/mol
    double kinetic() const { return raw[13]; }       // kJ/mol
    double totalEnergy() const { return raw[14]; }   // kJ/mol
    double enthalpy() const { return raw[15]; }      // kJ/mol
    double temperature() const { return raw[16]; }   // K
    double pressure() const { return raw[17]; }      // bar (scalar)
    double volume() const { return raw[18]; }        // nm³
    double density() const { return raw[19]; }       // kg/m³
    // Box (nm)
    Vec3 box() const { return Vec3(raw[20], raw[21], raw[22]); }
    // Virial tensor (kJ/mol), row-major XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
    Mat3 virial() const { return Mat3RowMajor(&raw[23]); }
    // Pressure tensor (bar), same ordering
    Mat3 pressureTensor() const { return Mat3RowMajor(&raw[32]); }
    // Per-group temperature (K)
    double tProtein() const { return raw[41]; }
    double tNonProtein() const { return raw[42]; }

private:
    static Mat3 Mat3RowMajor(const double* p) {
        Mat3 m;
        m << p[0], p[1], p[2],
             p[3], p[4], p[5],
             p[6], p[7], p[8];
        return m;
    }
};

// aimnet2_aim (N×256) — AIMNet2's per-atom learned electronic-structure
// embedding (the "aim" vector); geometry-dependent, encodes hybridisation,
// charge response, conjugation, charge transfer (AIMNet2Result.h; _tensors.py
// AIMNet2AimEmbedding). FLOAT32 ON DISK — the ONLY float32 array in the format
// (AIMNet2Result.cpp:492 WriteFloat32, native torch precision); the catalog has
// no dtype field, so the loader (#4) reads `<f4` from the NPY header and widens
// to double, and this view then reads doubles.
//
// NON-OWNING view into the snapshot's column buffer (256 contiguous doubles for
// one atom) — avoids a 2 KB copy per access for the GNN / ML-feeding path.
// Valid only while the owning QtConformationSnapshot is alive; do NOT retain a
// view after the source frame is released (same transient-use contract as the
// group views).
struct AIMNet2Embedding {
    static constexpr std::size_t kDim = 256;  // == AIMNET2_AIM_DIMS (library)
    const double* data = nullptr;
    std::size_t size() const { return kDim; }
    double operator[](std::size_t i) const { return data[i]; }
    const double* begin() const { return data; }
    const double* end() const { return data + kDim; }
};

// ── Tripeptide / Larsen reference-shielding family ─────────────────────
// These two families deliver DFT-derived reference shieldings as 9-col
// SphericalTensors (UnpackSphericalTensor) plus a few Vec3 residuals and
// scalars — no bespoke block beyond this method-provenance tag. The two
// families use OPPOSITE per-atom "no value here" sentinels (verified in the
// fixture): Tripeptide writes NaN (check std::isnan / Vec3::hasNaN());
// Larsen writes 0.0 (the per-class tensors are packed unconditionally, so an
// atom outside Larsen's Table-2 dispatch carries a structural ZERO). Both
// distinct from a group view's nullopt = "calculator absent this frame".
//
// tripeptide_bb_method_tag is int8 on disk (widened to double by the loader);
// it records which DFT frame_type produced the matched ProCS15 tripeptide pose
// (TripeptideBackboneShieldingResult.cpp:51 MethodTagFromFrameType;
// ConformationAtom.h:283-291). NoMatch(0) is the in-band integer twin of a NaN
// tripeptide_bb_shielding row (fixture: the 38 tag==0 atoms are exactly the 38
// NaN shielding rows). The OPBE/ORCA-PBE split is load-bearing for the
// methods-mixing caveat (a neighbour term drawing on a SER ASA row mixes
// OPBE + PBE — project_serine_pbe_discontinuity).
enum class TripeptideMethodTag : std::int8_t {
    NoMatch = 0,   // tripeptide_bb_has_match == false; shielding row is NaN
    Opbe    = 1,   // gaussian_standard_orientation (OPBE)
    OrcaPbe = 2,   // orca_input_orientation (PBE)
};

// ── MOPAC (PM7+MOZYME) family blocks ───────────────────────────────────
// The MOPAC family splits into three SDK FieldGroups: MOPACCore (the three
// blocks below), MOPACCoulomb (reuses the Coulomb blocks — SphericalTensor /
// Vec3 / QtEfg / CoulombScalars), and MOPACMcConnell (reuses SphericalTensor /
// PerBondCategoryT2 / McConnellScalars). MOPAC runs ONLY on the FullFat
// `--mopac` trajectory path, where it is the QM charge source reconciled
// against ff14SB; absent elsewhere → the group views return nullopt.

// mopac_scalars (N×4) — per-atom PM7+MOZYME electronic summary
// (MopacResult.cpp:494-503). `charge` DUPLICATES the standalone mopac_charges
// array (both are ConformationAtom.mopac_charge).
//   charge  : Mulliken partial charge (e)  [== mopac_charges]
//   sPop    : Mulliken s-orbital population (electrons)
//   pPop    : Mulliken p-orbital population (electrons)
//   valency : Wiberg valency = Σ of this atom's bond orders
struct MopacScalars {
    double charge = 0.0;
    double sPop = 0.0;
    double pPop = 0.0;
    double valency = 0.0;
    static MopacScalars FromRow(const double* r) {
        return MopacScalars{r[0], r[1], r[2], r[3]};
    }
};

// mopac_bond_orders (B×3, BOND axis) — one row per MOPAC-reported atom pair
// (MopacResult.cpp:506-526). B is MOPAC's OWN unique-pair count
// (bond_order_map_.size(); 896 in the fixture) and is INDEPENDENT of the
// topology bond count — PM7 Wiberg reports pairs the covalent topology may not
// list, so this is its own sparse axis, NOT parallel to QtProtein's bonds. The
// bond-axis ORDINAL is the writer's unordered_map iteration order (arbitrary,
// not stable across runs); the real join back to structure is atomA/atomB.
//   atomA, atomB : the pair's atom indices, atomA < atomB (writer PairKey =
//                  min<<32 | max). Stored as doubles on disk; cast back here.
//   wibergOrder  : Wiberg bond order (continuous, conformation-dependent)
struct MopacBondOrder {
    std::size_t atomA = 0;
    std::size_t atomB = 0;
    double wibergOrder = 0.0;
    static MopacBondOrder FromRow(const double* r) {
        return MopacBondOrder{static_cast<std::size_t>(r[0]),
                              static_cast<std::size_t>(r[1]), r[2]};
    }
};

// mopac_global (4, PROTEIN axis) — one row for the whole frame
// (MopacResult.cpp:529-532). On-disk shape is 1-D (4,), NOT (1,4); the loader
// (#4) interprets it as 1 protein-row × 4 cols (catalog cols=4 authoritative),
// the same protein-axis handling as gromacs_energy.
//   heatOfFormation : PM7 FINAL HEAT OF FORMATION (kcal/mol; MopacResult.cpp:435)
//   dipole (Vec3)   : molecular dipole moment (Debye; MopacResult.cpp:212, the
//                     MOPAC DIPOLE-table SUM line). NOTE: for a net-charged
//                     system the dipole is ORIGIN-DEPENDENT, not a translation-
//                     invariant observable (the 1P9J fixture is net −1 e, |d|≈159
//                     D about MOPAC's internal origin) — compare across frames
//                     only with that caveat in mind.
struct MopacGlobal {
    double heatOfFormation = 0.0;
    Vec3 dipole = Vec3::Zero();
    static MopacGlobal FromRow(const double* r) {
        MopacGlobal g;
        g.heatOfFormation = r[0];
        g.dipole = Vec3(r[1], r[2], r[3]);
        return g;
    }
};

}  // namespace h5reader::model

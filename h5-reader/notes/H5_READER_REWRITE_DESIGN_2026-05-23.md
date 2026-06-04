#H5 - Reader Rewrite Design — 2026 - 05 - 23

> **Historical — not current truth (trued 2026-06-04).** This is the original
> typed-mirror design. The rewrite landed and evolved; current architecture and
> load behavior live in `README.md`, `notes/SCOPE.md`, and
> `notes/UI_STATE_OVERVIEW_2026-06-04.md`.

**Status:** DRAFT, pending user review. No code lands against this design until
approved.

**Provenance:** Synthesis of three parallel investigation-agent runs against
the post-2026-05-13 trajectory.h5 + sidecar format and the library's typed
object model. Source data:

- Agent 1: library typed-object inventory (~100 entities across 4 layers; 8
  strings flagged all projection-only; 0 string-keyed maps).
- Agent 2: library → output projection map (per-source projection tables, gap
  inventory, Snobol-risk flags at output boundaries).
- Agent 3: 1P9J fixture deep-dive contract (`/shared/2026Thesis/1p9j_baseline_runs/run_01_baseline_2026-05-21_sha7d8dbe9/`,
  56 TR groups + 5 sidecar NPYs verified end-to-end).

The user's load-bearing design directive, verbatim:

> "The no strings rule has a reason. If there are _any_ strings then the work
> becomes snobol hell, and you never question the comparison later, and that
> is incredibly bad. We need a vision that is based on types and classes, not
> embedded strings."

> "It would be better to echo the structure of the library and repopulate it."


---

## 1. Vision

The h5-reader's typed object model **structurally mirrors the
nmr-shielding library's typed object model** — same vocabulary, same
hierarchy shape, same const-after-load discipline — repopulated from the
library's serialised state (per-TR-group H5 + 5-NPY topology sidecar +
extraction manifest) rather than recomputed.

The reader is a **trajectory animator** — single protein, T frames,
volumetric BS/HM butterfly re-eval at ring centres (the one place
kernel evaluation lives in the reader). Per-atom chemistry comes from
the wire, never from the atom name string.

**Structural mirror means**: where the library has `Atom` /
`AtomSemanticTable` / `Residue` / `Bond` / `Ring` / `RingTopology` /
`LegacyAmberTopology` / `Protein`, the reader has Qt-prefixed
counterparts with identical typed enum fields and the same predicate
methods (`IsBackbone()`, `IsInAnyRing()`, `IsPolarH()`). Where the
library has `ConformationAtom` with 80+ computed fields written by
~30 ConformationResult singletons, the reader has per-frame view
accessors that read from per-TR-group typed buffers — the library's
ConformationResult subclass hierarchy collapses into the reader's
per-TR buffer family, because the reader consumes results, doesn't
produce them.


---

## 2. The No-Strings Discipline

### The rule

**No string field on any chemistry-typed object dispatches identity or
participates in a comparison that decides physics.** Strings exist as
output projections accessed through an explicit `QtAtomNames` /
`QtResidueNames` projection layer, separated from the typed identity
objects.

The library has 8 strings total across 100+ typed entities (per Agent 1)
and zero string-keyed dispatch maps. The reader mirrors that discipline.

### Concrete consequences

**QtAtom must have no atom-name field**. The library's `Atom` has
`pdb_atom_name` (a single string), used only at the loader boundary and
for display. The reader's atom is fully described by typed substrate
fields (Element, Locant, BranchAddress, DiastereotopicIndex,
BackboneRole, PlanarGroupKind, PolarHKind, ProchiralStereo,
PseudoatomKind, ring positions, etc.). The three atom-name strings in
`atoms_category_info.npy` (amber/iupac/bmrb) belong on a separate
**projection** object (`QtAtomNames`) that the inspector dock and ML-
matching consumer ask for **by name** when they need a display label.

**QtResidue must have no name-string field**. The 3-letter codes and
1-letter code are derivable from `aminoAcid` (an enum) via a typed
helper that mirrors the library's `GetAminoAcidType(AminoAcid).three_letter_code`.
The `chain_id` and `insertion_code` are addressing strings (not
chemistry) — they stay, but behind a `QtChainAddress` typed wrapper so
they can't accidentally be compared like chemistry identity. The
addressing wrapper is `struct QtChainAddress {
    QString chain_id;
    int residue_number;
    QString insertion_code;
}` with deliberately
    - disabled
`operator==` requiring
          an explicit `IsSameAddress()` predicate
              .

                  ** Selections kinds must be
      typed**.The H5's `/trajectory/selections/` subgroups use mangled C
      ++ `type_index` names(`N3nmr34DftPoseCoordinatorTrajectoryResultE` etc.)
              .The reader has a static lookup table mapping these mangled strings → typed
`QtSelectionKind` enum(DftPoseCoordinator, RmsdSpikeSelection, ChiRotamerSelection)
              .The mangled names are read once at the loader boundary;
thereafter the reader holds the enum.

        ** Per
    - record `reason` and `metadata_json` strings** in selection subgroups are parsed at the loader into typed
`QtSelectionEvent` subclass instances.The raw strings stay inaccessible past the loader.

          ** Force
    - field atom type**(`ff_atom_type_string`, S4)is an AMBER ff14SB vocabulary of ~30 - 50 values.Two options : declare a
`QtFfAtomType` enum and map at load;
or leave as opaque QString flagged "PROJECTION — never dispatch".The first is more typed and matches the discipline; the FF type vocabulary is finite and stable
within ff14SB.

### Where strings remain (and why)

| Field | Where stored | Why it stays a string |
|---|---|---|
| `QtAtomNames::amber` (S8) | `QtAtomNames` (projection, separate from QtAtom) | True projection for display + ML matching against BMRB |
| `QtAtomNames::iupac` (S8) | `QtAtomNames` | Genuine atom_nom.tbl projection |
| `QtAtomNames::bmrb` (S8) | `QtAtomNames` | Genuine atom_nom.tbl projection |
| `QtChainAddress::chain_id` (S2) | `QtChainAddress` typed wrapper | Free-form addressing (multi-chain, multimers);
operator== deliberately disabled | | `QtChainAddress::insertion_code` (S1) | `QtChainAddress` | Free - form PDB convention |
    | `QtSelectionEvent::reason_text` | `QtSelectionEvent` (parsed at loader) | Free - form diagnostic text;
never compared | | `extraction_manifest.json` source paths | log - only,
    never on typed objects | Provenance |

        These six string surfaces are **boundary - only **.Once across the boundary,
    the rest of the reader works on typed objects and enum ordinals.


        -- -

        ##3. Mirror Map — Library → Reader

        | Library entity | Reader counterpart | Notes | | -- - | -- - | -- - | | `nmr::Protein` | `QtProtein`
        | Owns atoms / residues / bonds / rings / topology + the trajectory.Non - copyable,
    non - movable(same back - pointer discipline).| | `nmr::Atom` | `QtAtom`
        | Typed substrate fields ONLY.No `pdb_atom_name`.Names → `QtAtomNames`.| | `nmr::AtomSemanticTable`
        | (merged into `QtAtom`) | The atom IS its semantic record on the reader side; library separates them because the runtime needs the substrate-aware lookup `LookupBy(MechanicalIdentity)`. The reader is read-only. |
| `nmr::Residue` | `QtResidue` | Typed fields + `QtChainAddress` for addressing. No `amber_residue_3letter` etc. — derivable via `QtAminoAcidType`. |
| `nmr::Bond` | `QtBond` | Already typed-only. Add `is_rotatable / is_aromatic / is_peptide / is_backbone` flags from bonds.npy. |
| `nmr::Ring` hierarchy | `QtRing` hierarchy | Existing 8 aromatic subclasses + new `QtProPyrrolidineRing` (saturated). Virtual `Intensity()`, `NitrogenCount()`, etc. |
| (new) | `QtRingMembership` | Per-(ring, vertex-atom) typed record. |
| `nmr::LegacyAmberTopology` | `QtTopology` | Holds `QtCovalentTopology` + `QtRingTopology` + `QtAtomSemanticTable` field array. Wrapped under a single `QtTopology` object on the reader side rather than mirroring the "LegacyAmber" name — the reader doesn't carry the variant burden. |
| `nmr::CovalentTopology` | `QtCovalentTopology` | Holds bond list + per-atom bond-index cache. |
| `nmr::RingTopology` | `QtRingTopology` | Aromatic[] + Saturated[] split. |
| `nmr::AminoAcidType` + `GetAminoAcidType()` | `QtAminoAcidType` + `GetQtAminoAcidType(QtAminoAcid)` | Static registry;
reader - side mirror of the library's static lookup. Provides 3-letter / 1-letter codes derivable from enum + `is_aromatic / is_titratable / has_amide_h / chi_angle_count`. |
    | `nmr::ChargeSource` family | NOT MIRRORED | Loader concept;
reader doesn't load charges (they come through TRs if present). | | `nmr::ForceFieldChargeTable` | NOT MIRRORED | Same.|
    | `nmr::NamingApplicator` | NOT MIRRORED | Input - side; reader's input is already canonical. |
| `nmr::ProteinBuildContext` | `QtProteinBuildContext` (optional v2) | Provenance for inspector. v1 stub. |
| `nmr::ProteinConformation` | `QtConformation` | The trajectory. Owns frames + per-TR buffers + Welford rollups + selection bags. |
| `nmr::ConformationAtom` (80+ fields) | `QtFrameAtomView` | Per-(atom, frame) view returning std::optional<> per result-typed field. The library's "one struct holds everything" becomes "compose at access time from per-TR buffers". |
| `nmr::ConformationResult` subclass family (~30) | NOT MIRRORED AS CLASSES | The reader consumes outputs, doesn't compute. Equivalent surface lives in per-TR buffer classes. |
| `nmr::RingNeighbourhood` | `QtRingNeighbourhoodSlot` | Per-(atom, ring-slot) typed record;
sourced from `/ trajectory / ring_neighbourhood_trajectory_stats`.| | `nmr::BondNeighbourhood` | (deferred)
    | The new format doesn't carry per-(atom, bond) per-frame data; reader can re-derive from bond list + positions if needed. |
    | `nmr::TrajectoryProtein` | (merged into `QtProtein`) | The reader's QtProtein IS the trajectory-context protein. No separate class. |
    | `nmr::TrajectoryAtom` | `QtTrajectoryAtomStats` | Per - atom Welford rollup mirror;
one per atom, populated from atom - axis Welford datasets.| | `nmr::TrajectoryResult` subclass family(~40)
                  | NOT MIRRORED AS CLASSES | Same reasoning as ConformationResult — reader consumes via per - TR buffers.|
                  | `nmr::RecordBag<T>` | `QtSelectionBag` | At run - scope only(selections); per-atom events not in the H5 yet (Agent 3). |
| `nmr::SelectionRecord` | `QtSelectionEvent` (typed subclass per kind) | DftPoseSelection / RmsdSpikeSelection / ChiRotamerSelection with typed metadata. |
| `nmr::TrajectoryEnv` (per-frame env) | NOT MIRRORED | Library concept for per-frame stash during Run;
reader has no such phase.| | All enums in `src / Types.h` + `src / SemanticEnums.h` | `QtTypes.h` + `QtSemanticEnums.h`
    | Ordinal - identical mirror.|
    | Value structs(SphericalTensor,
                    FieldValue,
                    BranchAddress,
                    RingMembership,
                    RingPosition,
                    PseudoatomMembership,
                    BondOrderMask,
                    AtomMechanicalIdentity)
    | Same names with Qt prefix | Ordinal / layout - identical mirror.|


    -- -

        ##4. Layer A — Identity
        / Topology Mirror

        ## #4.1 QtProtein

```cpp class QtProtein {
public:
    // Non-copyable, non-movable — same back-pointer discipline as nmr::Protein.
    const QString& proteinId() const;  // from manifest.protein_id

    size_t atomCount() const;
    const QtAtom& atom(size_t i) const;

    size_t residueCount() const;
    const QtResidue& residue(size_t i) const;

    size_t bondCount() const;
    const QtBond& bond(size_t i) const;

    // Ring axis: aromatic first, then saturated. ringCount() == aromatic + saturated.
    size_t ringCount() const;
    size_t aromaticRingCount() const;
    size_t saturatedRingCount() const;
    const QtRing& ring(size_t i) const;  // polymorphic via std::unique_ptr<QtRing>

    size_t ringMembershipCount() const;
    const QtRingMembership& ringMembership(size_t i) const;

    const QtTopology& topology() const;

    // Loader-side optional layer; v1 stub.
    const QtProteinBuildContext& buildContext() const;

    // The trajectory companion. One per QtProtein; owned.
    const QtConformation& trajectory() const;

    // Backbone walk — mirror of library's BackbonePredecessor/BackboneSuccessor;
    // typed int32 walks via residue prev/next links, no chain_id strings.
    std::optional<size_t> backbonePredecessor(size_t residueIdx) const;
    std::optional<size_t> backboneSuccessor(size_t residueIdx) const;

private:
    friend class QtProteinLoader;
    // Fields populated by the loader; const accessors after load.
};
```

    ## #4.2 QtAtom — typed substrate only

```cpp struct QtAtom {
    // Identity indices
    int32_t atomIndex = -1;
    int32_t residueIndex = -1;
    int32_t parentAtomIndex = -1;  // for H atoms

    // Mechanical-identity tuple (mirrors AtomMechanicalIdentity)
    Element element = Element::Unknown;
    Locant locant = Locant::None;
    BranchAddress branch = {};
    DiastereotopicIndex diIndex = DiastereotopicIndex::None;
    BackboneRole backboneRole = BackboneRole::None;

    // Chemistry-substrate fields (mirrors AtomSemanticTable)
    ProchiralStereo prochiral = ProchiralStereo::NotProchiral;
    PlanarGroupKind planarGroup = PlanarGroupKind::None;
    PlanarStereo planarStereo = PlanarStereo::NotApplicable;
    PolarHKind polarH = PolarHKind::NotPolar;

    RingPositionLabel ringPositionPrimary = RingPositionLabel::NotInRing;
    RingPositionLabel ringPositionSecondary = RingPositionLabel::NotInRing;

    PseudoatomKind pseudoatomKind = PseudoatomKind::None;
    bool inSuperGroup = false;
    bool aromatic = false;
    int8_t formalCharge = 0;
    bool isExchangeable = false;
    int8_t equivalenceClass = 0;

    // Force-field atom type — proposed typed enum (see Open Question §11.B)
    QtFfAtomType ffAtomType = QtFfAtomType::Unknown;

    // Predicate methods (mirror AtomSemanticTable's IsBackbone() etc.)
    bool IsBackbone() const { return backboneRole != BackboneRole::None; }
    bool IsBackboneNitrogen() const;
    bool IsBackboneAlphaCarbon() const;
    bool IsBackboneAmideHydrogen() const;
    bool IsAnyAlphaHydrogen() const;
    bool IsSidechainCarboxylateOxygen() const;
    bool IsSidechainAmideOxygen() const;
    bool IsPolarH() const { return polarH != PolarHKind::NotPolar; }
    bool IsInAnyRing() const;
    // Element helpers via free functions on Element enum
};
```

    **Note on what's NOT here**: no `pdb_atom_name`, no `amber_atom_name`, no `iupac_atom_name`,
    no `bmrb_atom_name`.Atom names are accessed through `QtAtomNames` (§4.3).

    ## #4.3 QtAtomNames — separate projection layer

```cpp struct QtAtomNames {
    QString amber;  // ff14SB canonical, S8 from atoms_category_info.npy
    QString iupac;  // atom_nom.tbl projection
    QString bmrb;   // atom_nom.tbl projection
    QtNamingProvenance iupacProvenance = QtNamingProvenance::Match;
    QtNamingProvenance bmrbProvenance = QtNamingProvenance::Match;
};

class QtProtein {
public:
    const QtAtomNames& atomNames(size_t atomIdx) const;  // explicit projection access
};
```

The boundary is explicit: code that wants a name asks the protein for
the names projection. Code that asks `atom.IsBackbone()` operates on
typed enums. The two surfaces never mix.

### 4.4 QtResidue

```cpp
struct QtResidue {
    int32_t residueIndex = -1;
    QtChainAddress address;  // chain_id + insertion_code wrapped

    QtAminoAcid aminoAcid = QtAminoAcid::Unknown;
    int8_t protonationVariantIndex = -1;
    QtTerminalState terminalState = QtTerminalState::Internal;

    // Bond-graph-derived chain links
    int32_t prevResidueIndex = -1;
    int32_t nextResidueIndex = -1;
    QtAminoAcid prevResidueType = QtAminoAcid::Unknown;
    QtAminoAcid nextResidueType = QtAminoAcid::Unknown;

    int32_t atomCount = 0;
    bool isXProContext = false;  // bond-graph derived; one typed flag suffices

    // Atom membership (built from atoms_category_info.residue_index walk)
    std::vector<int32_t> atomIndices;

    // Backbone index cache — built from atoms' typed BackboneRole, no string scan
    static constexpr int32_t NONE = -1;
    int32_t N = NONE, CA = NONE, C = NONE, O = NONE, H = NONE, HA = NONE, CB = NONE;

    // Predicates derived from QtAminoAcid via QtAminoAcidType registry
    bool IsAromatic() const;
    bool IsTitratable() const;
    bool HasAmideH() const;
    bool IsProline() const { return aminoAcid == QtAminoAcid::PRO; }
    bool IsTerminalN() const;
    bool IsTerminalC() const;
};
```

        The five string fields in residues
            .npy(amber_residue_3letter, iupac_residue_3letter, one_letter, plus chain_id / insertion_code inside
`QtChainAddress`) become derivable. `QtResidueNames` returns the 3
        - letter
    and 1 - letter on demand from the AminoAcid enum
            + variant index.

              ## #4.5 QtBond

```cpp struct QtBond {
    int32_t bondIndex = -1;
    int32_t atomIndexA = -1;
    int32_t atomIndexB = -1;
    QtBondOrder order = QtBondOrder::Unknown;
    QtBondCategory category = QtBondCategory::Unknown;
    bool isRotatable = false;
    bool isAromatic = false;
    bool isPeptide = false;
    bool isBackbone = false;
};
```

    Already typed.Per
    - frame geometry(length, midpoint, direction)
returned by `QtFrame::bondLength(bondIdx)` etc.

### 4.6 QtRing hierarchy

Existing 8 aromatic subclasses plus new `QtProPyrrolidineRing`
(saturated, RingAromaticity::None, intensity=0). Base class adds
`ringKind` field and `RingKind::Aromatic` vs `RingKind::Saturated`
discriminator. Virtual `Intensity()`, `LiteratureIntensity()`,
`JBLobeOffset()`, `NitrogenCount()`, `Aromaticity()`, `RingSizeValue()`,
`TypeName()` (TypeName() is the one string surface here, returning a
`const char*` literal — same as library's pattern).

### 4.7 QtRingMembership

```cpp
struct QtRingMembership {
    int32_t ringId = -1;
    int32_t atomIndex = -1;
    int8_t ringAtomOrder = 0;  // 0..ringSize-1, canonical walk order
    bool isVertex = false;
    bool isSubstituent = false;
};
```

### 4.8 QtAminoAcidType + GetQtAminoAcidType

Mirror of library's `AminoAcidType` registry. Static const array indexed
by `QtAminoAcid` enum:

```cpp
struct QtAminoAcidType {
    QtAminoAcid type;
    const char* three_letter_code;  // "ALA", "HIS", ...
    char one_letter_code;           // 'A', 'H', ...
    bool is_aromatic;
    bool is_titratable;
    bool has_amide_H;
    int chi_angle_count;
    // ProtonationVariant table per type (HID/HIE/HIP for HIS, ASH for ASP, etc.)
};

const QtAminoAcidType& GetQtAminoAcidType(QtAminoAcid aa);
```

This is the helper that makes residue name strings derivable. The
3-letter codes are `const char*` literals in the registry — exactly
as the library does it.

### 4.9 QtChainAddress

```cpp
struct QtChainAddress {
    QString chainId;
    int residueNumber = 0;
    QString insertionCode;

    bool IsSameAddress(const QtChainAddress& other) const {
        return chainId == other.chainId && residueNumber == other.residueNumber && insertionCode == other.insertionCode;
    }

    // operator== DELETED — comparison must be explicit via IsSameAddress.
    bool operator==(const QtChainAddress&) const = delete;
};
```

The deleted operator== is intentional. It catches the "snobol-trap"
pattern where someone might compare addresses thinking they're comparing
chemistry.


---

## 5. Layer B+C — Trajectory + Per-Frame Mirror

### 5.1 QtConformation — owns frames + per-TR buffers

```cpp
class QtConformation {
public:
    const QtProtein* protein() const;

    size_t frameCount() const;
    const QtFrame& frame(size_t t) const;

    // Per-TR buffer access — typed by buffer family, sparse-tolerant via optional<>.
    // Each buffer is loaded once at trajectory.h5 open and held for the session.
    std::optional<const QtShieldingTimeSeries&> bsShielding() const;
    std::optional<const QtShieldingTimeSeries&> hmShielding() const;
    std::optional<const QtShieldingTimeSeries&> mcShielding() const;
    std::optional<const QtShieldingTimeSeries&> piQuadShielding() const;
    std::optional<const QtShieldingTimeSeries&> ringChiShielding() const;
    std::optional<const QtShieldingTimeSeries&> dispShielding() const;
    std::optional<const QtShieldingTimeSeries&> hbondShielding() const;
    std::optional<const QtShieldingTimeSeries&> mopacCoulombShielding() const;
    std::optional<const QtShieldingTimeSeries&> mopacMcShielding() const;
    std::optional<const QtShieldingTimeSeries&> mopacVsFf14sbReconciliation() const;
    std::optional<const QtShieldingTimeSeries&> tripeptideBbShielding() const;
    std::optional<const QtShieldingTimeSeries&> tripeptideNeighborShielding() const;
    std::optional<const QtShieldingTimeSeries&> larsenHBond1pHBShielding() const;
    std::optional<const QtShieldingTimeSeries&> larsenHBond1pHaBShielding() const;
    std::optional<const QtShieldingTimeSeries&> larsenHBond2pHBShielding() const;
    std::optional<const QtShieldingTimeSeries&> larsenHBond2pHaBShielding() const;
    std::optional<const QtShieldingTimeSeries&> waterFieldShielding() const;

    std::optional<const QtVec3TimeSeries&> apbsEfield() const;
    std::optional<const QtT2TimeSeries&> apbsEfg() const;
    std::optional<const QtVec3TimeSeries&> tripeptideBbResidualVec() const;
    std::optional<const QtVec3TimeSeries&> tripeptideNeighborResidualVecPrev() const;
    std::optional<const QtVec3TimeSeries&> tripeptideNeighborResidualVecNext() const;
    std::optional<const QtScalarTimeSeries&> aimnet2Charge() const;
    std::optional<const QtVec3TimeSeries&> aimnet2ChargeResponseGradientVec() const;
    std::optional<const QtScalarTimeSeries&> aimnet2ChargeResponseGradientScalar() const;
    std::optional<const QtEmbeddingTimeSeries&> aimnet2Embedding() const;  // (N, T, 256) float32
    std::optional<const QtScalarTimeSeries&> sasa() const;
    std::optional<const QtScalarTimeSeries&> jCoupling() const;
    std::optional<const QtScalarTimeSeries&> larsenHBondCount() const;
    std::optional<const QtScalarTimeSeries&> larsenHBondWaterTerm() const;
    std::optional<const QtScalarTimeSeries&> bondedEnergyTotal() const;
    std::optional<const QtTagTimeSeries&> tripeptideBbMethodTag() const;

    std::optional<const QtPerResidueDihedralTimeSeries&> dihedrals() const;
    std::optional<const QtPerResidueDssp8TimeSeries&> dssp8() const;
    std::optional<const QtPerResidueRingPuckerTimeSeries&> ringPucker() const;
    std::optional<const QtRingNeighbourhoodTimeSeries&> ringNeighbourhood() const;  // (N, T, 11, 4)
    std::optional<const QtPerBondLengthStats&> bondLengthStats() const;

    std::optional<const QtSystemEnergyTimeSeries&> gromacsEnergy() const;
    std::optional<const QtRmsdTracking&> rmsdTracking() const;

    // Welford rollups (atom-axis, no T dimension)
    std::optional<const QtShieldingWelford&> bsWelford() const;
    std::optional<const QtShieldingWelford&> hmWelford() const;
    std::optional<const QtShieldingWelford&> mcWelford() const;
    std::optional<const QtScalarWelford&> sasaWelford() const;
    std::optional<const QtScalarWelford&> eeqWelford() const;
    std::optional<const QtScalarWelford&> hbondCountWelford() const;
    std::optional<const QtScalarWelford&> mopacChargeWelford() const;
    std::optional<const QtBondOrderWelford&> mopacBondOrderWelford() const;
    std::optional<const QtScalarWelford&> waterFieldWelford() const;
    std::optional<const QtHydrationWelford&> hydrationShellWelford() const;
    std::optional<const QtHydrationWelford&> hydrationGeometryWelford() const;
    std::optional<const QtScalarWelford&> aimnet2ChargeResponseGradientWelford() const;
    std::optional<const QtAutocorrelation&> bsT0Autocorrelation() const;

    // Transition trackers
    std::optional<const QtDssp8Transitions&> dssp8Transitions() const;
    std::optional<const QtDihedralBinTransitions&> dihedralBinTransitions() const;

    // Selection bags (run-scope events)
    const QtSelectionBag& selections() const;

    // Frame metadata
    double frameTime(size_t t) const;           // /trajectory/frames/time_ps
    size_t originalFrameIndex(size_t t) const;  // /trajectory/frames/original_index

    // Trajectory source provenance
    const QtTrajectorySource& source() const;  // xtc_path/tpr_path/edr_path
};
```

The `std::optional<const T&>` pattern (technically returned via
pointer-or-similar; the exact API can use `const T*` or
`std::optional<std::reference_wrapper<const T>>` — implementation
detail) makes the sparse-set discipline explicit at every consumer
call site. A TR group's absence is normal, not an error.

### 5.2 QtFrame — per-frame typed view

```cpp
class QtFrame {
public:
    size_t frameIndex() const;
    double timePs() const;
    size_t originalIndex() const;

    // Per-atom position — always available (positions TR is always-attached)
    QtVec3 position(size_t atomIdx) const;

    // Per-atom typed accessors — return std::optional<> when the underlying TR
    // is absent OR when source_attached_per_frame[t] == 0 for that frame.
    std::optional<QtSphericalTensor> bsShielding(size_t atomIdx) const;
    std::optional<QtSphericalTensor> hmShielding(size_t atomIdx) const;
    // ... per shielding family ...

    std::optional<double> sasa(size_t atomIdx) const;
    std::optional<QtVec3> apbsEfield(size_t atomIdx) const;
    // ... etc ...

    // Per-residue accessors
    std::optional<QtDsspCode> dsspCode(size_t residueIdx) const;
    std::optional<double> phi(size_t residueIdx) const;
    std::optional<double> psi(size_t residueIdx) const;
    std::optional<double> omega(size_t residueIdx) const;
    std::optional<double> chi(size_t residueIdx, int k) const;

    // Per-ring geometry (reconstructed from ring_membership + positions)
    QtRingGeometry ringGeometry(size_t ringIdx) const;
    std::vector<QtVec3> ringVertices(size_t ringIdx) const;

    // Composite typed view of one atom at this frame (materialized on demand)
    QtFrameAtomView atomView(size_t atomIdx) const;
};
```

`QtFrameAtomView` is the materialized typed struct for one (atom, frame)
pair — useful for the inspector dock which displays all an atom's
fields at the current frame. It's a **view**, not stored: constructed
on access from the per-TR buffers.

### 5.3 Per-TR buffer family

Six template-shaped buffer classes cover the 56 TR groups in the
fixture:

```cpp
// Per-atom T0+T1+T2 shielding [N, T, 9]
class QtShieldingTimeSeries {
public:
    size_t atomCount() const;
    size_t frameCount() const;
    QtSphericalTensor at(size_t atomIdx, size_t t) const;
    bool sourceAttachedAt(size_t t) const;  // per-frame attach mask
    const QtIrrepLayout& irrepLayout() const;
    const QString& units() const;
};

// Per-atom scalar [N, T]
class QtScalarTimeSeries { ... };

// Per-atom Vec3 [N, T, 3]
class QtVec3TimeSeries { ... };

// Per-atom T2-only [N, T, 5]
class QtT2TimeSeries { ... };

// Per-atom uint8 tag [N, T]
class QtTagTimeSeries { ... };

// Per-atom embedding [N, T, 256] float32
class QtEmbeddingTimeSeries { ... };

// Per-residue family (54 in fixture)
class QtPerResidueDihedralTimeSeries { ... };
class QtPerResidueDssp8TimeSeries { ... };
class QtPerResidueRingPuckerTimeSeries { ... };

// Per-(atom, ring-slot) [N, T, 11, 4]
class QtRingNeighbourhoodTimeSeries { ... };

// Per-bond statistics [B] (no T axis — already reduced)
class QtPerBondLengthStats { ... };

// System-scope [T] / [T, 9]
class QtSystemEnergyTimeSeries { ... };

// Welford rollups (atom-axis [N], no T axis)
class QtShieldingWelford {
public:
    QtWelfordMoments t0(size_t atomIdx) const;      // mean/std/min/max/min_frame/max_frame/m2
    QtWelfordMomentsVec3 t1(size_t atomIdx) const;  // per-component
    QtWelfordMomentsT2 t2(size_t atomIdx) const;    // per-component
    QtWelfordMoments t2magnitude(size_t atomIdx) const;
    QtWelfordMoments t0Delta(size_t atomIdx) const;
    QtWelfordMoments t0AbsDelta(size_t atomIdx) const;
    QtWelfordMoments t0DeltaSquared(size_t atomIdx) const;
    QtWelfordMoments t0Dxdt(size_t atomIdx) const;
    double t0RmsDelta(size_t atomIdx) const;
    size_t nFramesPerAtom(size_t atomIdx) const;
    double meanDtPs() const;
    int ddof() const;
};

class QtScalarWelford { ... };
class QtBondOrderWelford { ... };
class QtHydrationWelford { ... };
class QtAutocorrelation { ... };
```

Each TR-group reader at the io/ layer populates one of these.

### 5.4 QtTrajectoryAtomStats — Welford state per atom

For the inspector dock's per-atom Welford display, materialize a typed
view that mirrors `nmr::TrajectoryAtom`'s substruct shape:

```cpp
struct QtTrajectoryAtomStats {
    QtAtomShieldingWelfordSlice bs;  // moments across BS T0/T1/T2/|T2|/delta variants
    QtAtomShieldingWelfordSlice hm;
    QtAtomShieldingWelfordSlice mc;
    std::optional<QtAtomScalarWelfordSlice> sasa;
    std::optional<QtAtomScalarWelfordSlice> eeq;
    std::optional<QtAtomScalarWelfordSlice> hbondCount;
    // ... etc per Welford TR ...
};

QtTrajectoryAtomStats QtConformation::atomStats(size_t atomIdx) const;
```

### 5.5 QtSelectionBag — typed events

```cpp
enum class QtSelectionKind : int8_t {
    DftPoseCoordinator   = 0,
    RmsdSpikeSelection   = 1,
    ChiRotamerSelection  = 2,
};

// Static lookup table mapping mangled C++ type_index names → typed kind.
// Populated at compile time; demangling not required at runtime.
QtSelectionKind ParseSelectionGroupName(const std::string& mangled);

class QtSelectionBag {
public:
    size_t count() const;
    const QtSelectionEvent& at(size_t i) const;
    std::vector<size_t> indicesByKind(QtSelectionKind k) const;
    std::vector<size_t> indicesInTimeRange(double t_lo_ps, double t_hi_ps) const;
};

struct QtSelectionEvent {
    QtSelectionKind kind;
    size_t frameIdx;
    double timePs;
    // Typed metadata populated per kind by parsing the JSON at load:
    std::variant<QtDftPoseMeta, QtRmsdSpikeMeta, QtChiRotamerMeta, std::monostate> meta;
};
```

The `reason` and `metadata_json` strings in selection subgroups are
parsed at the loader. Past the loader, only typed `meta` variants exist.


---

## 6. Substrate Vocabulary (Layer D)

All enums mirrored from `src/Types.h` and `src/SemanticEnums.h` go into
`h5-reader/src/model/QtTypes.h` (existing — partially done) and a new
`h5-reader/src/model/QtSemanticEnums.h` covering the 14 typed-substrate
fields. Ordinal compatibility is checked at compile time via a
build-time test against the manifest's `enum_vocab` block.

Value structs mirrored: `QtSphericalTensor`, `QtFieldValue`,
`QtBranchAddress`, `QtPseudoatomMembership`, `QtRingMembership`,
`QtRingPosition`, `QtBondOrderMask`, `QtAtomMechanicalIdentity`,
`QtWelfordMoments`.


---

## 7. Loader Pipeline

```
QtProteinLoader::Load(h5_path) ->
  1. Resolve sidecar_dir = dirname(h5_path)
  2. Load extraction_manifest.json -> QtManifest (enum_vocab, axis_sizes, axis_alignment)
  3. Validate QtManifest.axis_sizes against the sidecar NPYs that follow
  4. QtNpyReader -> QtNpyAtomCategoryRow[846]    (atoms_category_info.npy)
                  QtNpyResidueRow[54]            (residues.npy)
                  QtNpyBondRow[862]              (bonds.npy)
                  QtNpyRingRow[16]               (rings.npy)
                  QtNpyRingMembershipRow[96]     (ring_membership.npy)
  5. Decode NPY rows -> QtAtom[], QtResidue[], QtBond[], QtRing[], QtRingMembership[]
     - Strings (amber_atom_name, iupac_atom_name, bmrb_atom_name) ->
       parallel QtAtomNames[] array, not on QtAtom
     - All typed enum int8 columns -> cast to typed enums with
       compile-time-validated ordinals
  6. Construct QtTopology from bonds + rings + ring_membership
  7. Construct QtProtein, finalize identity/topology
  8. Open trajectory.h5 via HighFive:
     - Validate root.protein_id matches sidecar.manifest.protein_id (or log warn)
     - Validate root.n_atoms == manifest.axis_sizes.atom
     - Load /atoms (3 datasets) and cross-check element + residue_index alignment
       with QtAtom[] from sidecar
     - Load /trajectory/frames -> QtFrameMetadata[T]
     - Load /trajectory/source attrs -> QtTrajectorySource
     - For each /trajectory/<name>/ subgroup present:
       a. Classify by axis_attr + dataset shape -> dispatch to per-TR buffer reader
       b. Eager-load full (N, T, K) slab into typed buffer
       c. Load source_attached_per_frame mask if present
       d. Store in QtConformation's per-TR slot
     - For each /trajectory/selections/<mangled-name>/ subgroup:
       a. Look up mangled name -> QtSelectionKind via static table
       b. Parse reason + metadata_json per record into typed QtSelectionEvent
       c. Push into QtSelectionBag
  9. Close HighFive::File. All buffers are now in QtConformation.
 10. Construct QtConformation, register it on QtProtein, return QtLoadResult.

Failure modes (each → ErrorBus log + graceful degradation):
 - Sidecar NPY missing/malformed: typed substrate absent → QtProtein constructed
   from H5 /atoms only with reduced typed fields (Element + residueIndex only)
 - extraction_manifest.json missing: enum_vocab validation skipped, ordinals
   trusted blindly; log Warn
 - axis_size mismatch: hard fail at QtTopologySidecar boundary
 - H5 missing: hard fail with specific error message
 - Single TR group malformed: skip group, log Warn, other groups continue
 - Mangled selection name unknown: skip subgroup, log Warn
```


---

## 8. NPY + H5 Boundary Classes

### 8.1 QtNpyReader (kept from prior work)

Already written. Generic NPY 1.0 reader with `ReadStructured<RecordT>`
template. Cross-platform via QFile + QByteArray. Logs failures via
ErrorBus.

### 8.2 QtNpyRecords (kept; static_asserts already fixed)

Five #pragma pack(1) byte-mirror structs:
- `QtNpyResidueRow` (42 bytes)
- `QtNpyBondRow` (18 bytes)
- `QtNpyRingRow` (24 bytes)
- `QtNpyRingMembershipRow` (12 bytes)
- `QtNpyAtomCategoryRow` (83 bytes)

### 8.3 QtTopologySidecar (new)

Orchestrates the five NPY reads + manifest parse. Produces:
- `QtAtom[]` + `QtAtomNames[]` (parallel arrays)
- `QtResidue[]`
- `QtBond[]`
- `QtRing[]` (unique_ptr, polymorphic)
- `QtRingMembership[]`
- `QtEnumVocab` (from manifest)
- `QtManifest` (axis_sizes + axis_alignment statements for cross-check)

Reads at the loader boundary, never at runtime. After construction, all
typed objects are populated; the NPY row structs are discarded.

### 8.4 QtTrajectoryH5 (new — heavily inspired by ui/src/TrajectoryH5)

Orchestrates the H5 read. Eager-loads all present TR groups into the
typed buffer family (§5.3). Each TR group's reader is a small file-local
function in this .cpp:

```
ReadShieldingTimeSeries(file, "/trajectory/bs_shielding_time_series", N, T)
  -> std::optional<QtShieldingTimeSeries>

ReadScalarTimeSeries(file, "/trajectory/sasa_time_series", "sasa", N, T)
  -> std::optional<QtScalarTimeSeries>

ReadVec3TimeSeries(file, "/trajectory/apbs_efield_time_series", "xyz", N, T)
  -> std::optional<QtVec3TimeSeries>

... etc per TR family ...
```

Constructor takes the H5 path; eager-loads everything; closes the
HighFive handle before returning. Memory footprint on 1P9J fixture:
~1.6 GB (the H5 is 1.8 GB; we drop ~10% format overhead).

### 8.5 QtEnumVocab — loaded from manifest

```cpp
class QtEnumVocab {
public:
    static QtEnumVocab Load(const QString& manifest_path);

    bool ValidateAtCompileTime() const;
    // Asserts each loaded enum's ordinal-vs-name pairs match the
    // reader-side typed enum's static naming. Catches drift between
    // library and reader at load time, not at runtime.

    QString DisplayName(QtEnum::Kind kind, int8_t ordinal) const;
    // For inspector tooltips that want the manifest-sanctioned name
    // rather than a hardcoded reader-side string.
};
```


---

## 9. Polymorphic Accessor — Type-Indexed Result Lookup

The library uses `type_index`-keyed maps for ConformationResult and
TrajectoryResult dispatch. The reader's per-TR buffer access is named
(one method per TR family on QtConformation) rather than type-indexed,
because:

1. The reader's TR vocabulary is fixed and finite (56 TR groups in
   fixture; ~25 will be commonly present).
2. Named methods are more discoverable and IDE-friendly.
3. There's no template polymorphism need — the reader doesn't dispatch
   on result type.

If a future need for `conformation.Result<QtBsShieldingTimeSeries>()`
arises, we can add it as a thin wrapper over the named accessors. Not
in v1.


---

## 10. Disposition of Existing Work

Eight files written before the design pause:

| File | Status | Action |
|---|---|---|
| `model/Types.h` | KEEP, MINOR EDITS | Trim a few derived helpers;
the substrate enums are sound.Move display - only helpers(NameForX)
to a separate QtDisplay.h to keep Types.h chemistry-pure. |
| `model/QtAtom.h` | REWRITE | Remove `amberAtomName`, `iupacAtomName`, `bmrbAtomName`, `ffAtomTypeString` QString fields. Add `ffAtomType` typed enum (or defer). Substrate fields are already correct. |
| `model/QtResidue.h` | REWRITE | Remove `amberResidue3Letter`, `iupacResidue3Letter`, `oneLetter` QString fields (derive via QtAminoAcidType). Replace standalone `chainId` + `insertionCode` strings with `QtChainAddress` typed wrapper. |
| `model/QtBond.h` | KEEP AS IS | Already typed-only. |
| `model/QtRingMembership.h` | KEEP AS IS | Typed-only. |
| `model/QtRing.h` (untouched in this session) | EDIT | Add `QtProPyrrolidineRing` saturated subclass. Add `RingKind` member. |
| `io/QtNpyRecords.h` | KEEP | Byte-mirrors of NPY rows;
strings are intrinsic to the file format.| | `io / QtNpyReader.{
    h, cpp
}
` | KEEP | Generic NPY 1.0 reader; sound. |

**Files to add (the substantial new code):**

| File | Purpose | Lines (est.) |
|---|---|---|
| `model/QtTypes.h` (rename from Types.h or extend) | Core enums (Element, BondOrder, etc.) | (existing) |
| `model/QtSemanticEnums.h` | Substrate-specific enums (separate from QtTypes.h for clarity) | ~400 |
| `model/QtAtomNames.h` | Projection layer for atom names | ~60 |
| `model/QtResidueNames.h` | Derived projections for residue names | ~50 |
| `model/QtChainAddress.h` | Typed wrapper for addressing strings | ~50 |
| `model/QtAminoAcidType.h` + `.cpp` | Static AA registry mirroring library | ~300 |
| `model/QtTopology.h` + `.cpp` | Wraps bonds + rings + atom_semantic | ~200 |
| `model/QtRing.h` extension | Add QtProPyrrolidineRing + RingKind | ~50 |
| `model/QtTimeSeriesBuffers.h` | The per-TR buffer family (10+ types) | ~800 |
| `model/QtWelfordBuffers.h` | Welford rollup buffer types | ~400 |
| `model/QtSelectionBag.h` + `.cpp` | Typed selection events + bag | ~200 |
| `model/QtFrameAtomView.h` | Per-(atom, frame) materialized view | ~100 |
| `io/QtNamingProvenance.h` | Enum + ordinal validation | ~30 |
| `io/QtEnumVocab.h` + `.cpp` | Load + validate manifest enum_vocab | ~200 |
| `io/QtManifest.h` + `.cpp` | Parse extraction_manifest.json | ~150 |
| `io/QtTopologySidecar.h` + `.cpp` | Orchestrate the 5 NPYs + manifest | ~500 |
| `io/QtTrajectoryH5.h` + `.cpp` | Orchestrate the H5 + 56 TR groups | ~1500 |
| `io/QtProteinLoader.{
    h, cpp
}
` (rewrite) | Top - level orchestrator | ~250 |

    **Total** : ~5300 lines of new code over 18 files; ~500 lines of
modifications to 4 existing files. Manageable in 3-4 implementation
sessions following the design.


---

## 11. Open Questions for User Decision

### A. ChainAddress comparison policy

Should `QtChainAddress::operator==` be:
- DELETED (force explicit `IsSameAddress()` calls; protect against
  Snobol-style misuse)
- DEFINED but `[[nodiscard]]` with a comment
- LEFT DEFAULT (regular struct equality)

**Recommendation**: DELETED. The user's discipline is "don't let future
agents accidentally compare strings"; deleting the operator is the
strongest enforcement.

### B. QtFfAtomType — typed enum or QString projection?

Force-field atom type (`ff_atom_type_string`, S4 in NPY) is an AMBER
ff14SB vocabulary of ~30-50 values ("CT", "HC", "N", "N3", etc.). Two
options:

- B1. Declare `enum class QtFfAtomType : int8_t { CT, HC, N, N3, ..., Unknown };
` in QtTypes.h.Loader maps S4 string → enum at boundary; runtime is
  fully typed. Snobol-discipline-clean.

- B2. Store as `QString` flagged "DISPLAY ONLY — never dispatch".
  Comment-discipline-only.

**Recommendation**: B1. The vocabulary is finite and stable within
ff14SB. The library writes only known FF types. A typed enum closes the
last string-dispatch vulnerability.

### C. QtAtomNames — eager arrays or lazy projection map?

Two implementations:
- C1. `QtProtein` holds `std::vector<QtAtomNames> atom_names_` parallel
  to atoms_ (~25 bytes × N atoms = trivial memory cost).
- C2. `QtProtein::atomNames(i)` materializes a `QtAtomNames` on call
  from a single shared blob of S8/S8/S8 columns kept compact.

**Recommendation**: C1. Memory cost is negligible (~25 KB for 1000-atom
protein); call-site clarity wins.

### D. Selection metadata parsing strictness

When a selection's `metadata_json` doesn't parse as the expected typed
metadata for its kind (e.g., a DftPose record has malformed JSON), the
loader should:

- D1. Log Warn, skip that record (lenient).
- D2. Log Error, skip that record but mark the bag "partial" (medium).
- D3. Hard fail the whole load (strict).

**Recommendation**: D1. Selections are diagnostic; one bad record
shouldn't kill loading 700 frames of physics.

### E. Build pace — single design-approved diff or staged?

The full rewrite is ~5800 lines. Should the implementation land:
- E1. Single design-approved branch with the whole rewrite at once.
- E2. Staged across 3-4 PRs (topology + sidecar first; per-TR families
  next; selections + transitions last).

**Recommendation**: E2. The staging matches the natural dependency
order (topology must work before per-TR can be tested) and lets the
user review each layer.

### F. Selection-kind mangled-name strategy

The reader needs to map mangled C++ `type_index` names to
`QtSelectionKind` enum. Two options:

- F1. Static lookup table in `QtSelectionBag.cpp` — hand-coded
  mangled→enum map for the 3 known kinds. Breaks on compiler change.
- F2. Demangle via `abi::__cxa_demangle` (libstdc++) — works on
  Linux/Mac, requires platform shim on Windows.

**Recommendation**: F1 with explicit comment. The 3 kinds are stable;
adding a 4th means adding one line. Cross-platform without
platform-specific code.

### G. Trajectory animation strategy for the volumetric BS/HM butterfly

The h5-reader is a TRAJECTORY animator; the existing volumetric BS/HM
butterfly overlay (`QtFieldGridOverlay`) recomputes a grid per ring per
frame. With 15 aromatic rings × 64³ grid × 7 evaluations/voxel = ~28M
evaluations per frame. At 20 fps that's 560M evals/sec — feasible on a
5090 but might tax the GB10/M3 hardware.

- G1. Re-evaluate per frame (current behaviour).
- G2. Cache "ring-local" grids and translate/rotate per frame.
- G3. Lower resolution + per-ring opacity blending.

**Recommendation**: defer to a separate design pass after the typed-
mirror is in place. Not blocking v1.


---

## 12. Implementation Sequence (post-approval)

**Session 1 — Foundation:**
1. Finish `model/QtTypes.h` + new `model/QtSemanticEnums.h` (substrate
   enums separated for clarity)
2. `model/QtChainAddress.h` (with deleted operator==)
3. `model/QtAminoAcidType.h` + `.cpp` (static AA registry)
4. Rewrite `model/QtAtom.h` (strip strings, add ffAtomType enum)
5. Rewrite `model/QtResidue.h` (strip strings, add QtChainAddress)
6. `model/QtAtomNames.h` (projection layer)
7. `model/QtResidueNames.h` (derivation helpers)
8. `model/QtRing.h` extension (QtProPyrrolidineRing + RingKind)

**Session 2 — Sidecar + topology:**
9. `io/QtManifest.h` + `.cpp` (parse extraction_manifest.json)
10. `io/QtEnumVocab.h` + `.cpp` (validate ordinals at load)
11. `io/QtTopologySidecar.h` + `.cpp` (orchestrate 5 NPYs)
12. `model/QtTopology.h` + `.cpp` (wrap bonds/rings/atoms)

**Session 3 — Per-TR buffers + H5:**
13. `model/QtTimeSeriesBuffers.h` (10+ buffer types)
14. `model/QtWelfordBuffers.h`
15. `io/QtTrajectoryH5.h` + `.cpp` (the big one)

**Session 4 — Selections, conformation, loader, wiring:**
16. `model/QtSelectionBag.h` + `.cpp`
17. `model/QtConformation` rewrite to hold typed buffers
18. `model/QtFrameAtomView`
19. Rewrite `io/QtProteinLoader`
20. Update CMakeLists.txt (drop fileformat dep)
21. Audit + delete dead AnalysisFile-touching code in overlays/docks

**Session 5 — Overlay + dock wire-up:**
22. Touch up overlays (ribbon, ring polygon, field grid, butterfly)
to consume the new buffer API 23. Touch up inspector dock to consume QtFrameAtomView 24. Touch up time
    - series dock to consume the typed buffer family 25. Build + run
    + verify on 1P9J fixture

              ** Session 6 — Cleanup
    : **26. Delete `fileformat
      / analysis_file.{
    h, cpp
}
` and `roundtrip_ *.cpp` 27. Audit JobSpec.h comment + ui / UI_ROADMAP.md doc references 28. Land final dead - format
    - removal commit

        Six sessions,
    each scoped to a coherent layer.After session 4 the project builds;
after session 5 it animates the fixture;
session 6 is the dead - format - deletion that closes out the migration.

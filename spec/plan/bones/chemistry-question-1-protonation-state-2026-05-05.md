# Chemistry Question 1 — Protonation state determination (2026-05-05)

Audit-and-define agent commissioned to articulate the protonation-state
chemistry question rigorously, classify where it currently lives, and
sketch a Phase 1 redesign that does not read `pdb_atom_name` strings on
the runtime side. Phase 2 (the future calculator-walk pass) gets only
flagged considerations, not a solution. This document follows the
2026-05-05 audit's hotspot-1 finding (`Protein::ResolveProtonationStates`).

The framing the user gave is exact: "you code based on name, name is
meaning, and some choices are virtually automatic. It is our job ... to
make sure we get a model where our algorithms can do science, if we
are careful and lucky, not just power to a result." The intent of this
analysis is to surface where the science-versus-string-power decision
sits for protonation, so that Phase 1 fixes the mechanism without
papering over the chemistry.

---

## Summary

`Protein::ResolveProtonationStates` (`src/Protein.cpp:306-414`) reads
`pdb_atom_name` at runtime and *infers* the residue's typed
`protonation_variant_index` from a cookbook of which-H-name-is-present.
This is chemistry inference dressed up as load housekeeping. Eight of
the nine titratable branches in the function consist of a string
disjunction over hydrogen names; the result is a typed enum index that
downstream code (charge resolver, ring detection, AMBER LEAP emitter)
correctly consumes.

The model already has the typed substrate the function should be reading
*from*: `AmberAminoAcidVariantTable` in `src/AminoAcidType.cpp` enumerates
HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN/ARN/TYM with their `force_field_name`,
`formal_charge`, and `description`. Three of the four load paths
(`FullSystemReader.cpp:740-748`, `OrcaRunLoader.cpp:189-198`, and the
trajectory-via-readback path) already produce the typed
`protonation_variant_index` directly from upstream evidence. Only the
PDB load path (`PdbFileReader::BuildFromPdb` via `reduce`) leaves the
field at -1 and depends on `ResolveProtonationStates` to fill it in
from atom-name patterns.

**Phase 1 direction.** Move name-based inference *into the PDB loader*
where the load context is fresh, replace name-equality with
typed-identity-presence (`Element::H, Locant::Delta, BranchAddress{1,0}`
finds HD1) using the `AtomMechanicalIdentity` substrate already landed
at commit `ee1f1b4`. `Protein::ResolveProtonationStates` becomes a
*consistency check* on already-resolved state, not a producer. The four
load paths converge on a single contract: by the time
`FinalizeConstruction` runs, every titratable residue has an authoritative
variant index from its own load source; the function trusts it.

**Phase 2 considerations.** Time-dependent protonation (titration during
MD), tautomer ambiguity in the absence of any H, pKa-environment
coupling under conformational change — all are real chemistry that the
current model does not treat. They are out of scope for the substrate
refactor; they require new calculator/result types and possibly a
ConformationResult (variant-by-frame) layered on Protein-level identity.

---

## 1. The chemistry question

Six of the standard 20 amino acids carry a titratable side chain; their
protonation state at solution pH determines side-chain charge,
side-chain tautomer (for HIS), the chemistry of the affected nitrogen
or oxygen, and — for CYS — whether the residue is part of a covalent
disulfide instead of a free thiol. The protonation state is INVARIANT
under fixed-protonation MD (the project's standard protocol per
2026-04-29 decision O1) and therefore belongs on the chemistry
substrate (`Protein` / `LegacyAmberTopology`), not on
`ProteinConformation`.

Per residue family:

- **HIS**. Three variants — HID (Nδ1-H, Nε2 lone pair, neutral), HIE
  (Nε2-H, Nδ1 lone pair, neutral), HIP (both NH, +1 charge,
  imidazolium). The tautomer choice has chemistry consequences beyond
  the charge: the imidazole's aromaticity character shifts (HIP's
  aromaticity is closer to canonical 6π imidazolium, HID/HIE differ in
  which ring N carries the lone pair), the ring-current intensity is
  variant-dependent (`HidImidazoleRing`, `HieImidazoleRing`,
  `HisImidazoleRing` are distinct typed `Ring` subclasses in this
  project), and the H-bond donor/acceptor topology depends on which N
  is protonated.
- **ASP / GLU**. Charged form (default, -1 on the second carboxylate
  oxygen Oδ2 / Oε2) vs protonated ASH/GLH (neutral, COOH form).
  Protonated form is the minority population at physiological pH but
  occurs in low-pKa active sites and buried environments.
- **CYS**. Free thiol (default, Hγ on Sγ) vs CYX (disulfide-bonded,
  no Hγ, Sγ is in an inter-residue covalent S-S bond) vs CYM
  (deprotonated thiolate, no Hγ, -1 on Sγ). The CYX/CYM distinction
  is critical: both lack the thiol H, but CYX is bonded to another
  CYX's Sγ while CYM is a free anion.
- **LYS**. Charged ammonium (default, +1, NH3+ on Nζ) vs LYN (neutral
  amine, NH2 on Nζ).
- **ARG**. Charged guanidinium (default, +1) vs ARN (neutral, deprotonated
  on Nε; very rare at physiological pH; ff14SB has no canonical ARN
  template per `topology-encoding-dependencies-2026-05-05.md` §B.1).
- **TYR**. Neutral phenol (default, Hη on Oη) vs TYM (deprotonated
  phenolate, no Hη, -1 on Oη). The phenolate carries qualitatively
  different ring-current character than neutral Tyr (per
  `SemanticEnums.h`'s `PlanarGroupKind::AromaticOxide`).

The chemistry question, stated precisely:

> *Given a typed `Protein` constructed from a load source, what
> determines the typed `Residue.protonation_variant_index` (and
> derivatively the variant's force-field name, formal charge,
> aromatic/planar-group classification, and ring chemistry) for each
> titratable residue?*

The answer is **the load source's own evidence**, projected onto the
typed variant index at the load boundary. The answer is NOT "scan the
runtime model's atom names for which H is present and infer."

---

## 2. What the current code does

`Protein::ResolveProtonationStates` runs twice during
`FinalizeConstruction` (once before bond detection, once after, see
`src/Protein.cpp:219` and `:253`). Both invocations execute the same
name-based inference for any residue with `protonation_variant_index < 0`.
Inside the loop (`src/Protein.cpp:325-413`), every titratable residue
follows the same shape:

```cpp
std::map<std::string, size_t> name_to_idx;
bool has_any_H = false;
for (size_t ai : res.atom_indices) {
    const Atom& atom = *atoms_[ai];
    name_to_idx[atom.pdb_atom_name] = ai;          // STRING DISPATCH SETUP
    if (atom.element == Element::H) has_any_H = true;
}

int variant_idx = res.protonation_variant_index;
bool resolved = res.protonation_state_resolved;

if (res.type == AminoAcid::HIS) {
    bool has_HD1 = name_to_idx.find("HD1") != name_to_idx.end();   // STRING
    bool has_HE2 = name_to_idx.find("HE2") != name_to_idx.end();   // STRING
    if (has_HD1 && has_HE2)      variant_idx = 2;  // HIP
    else if (has_HD1)            variant_idx = 0;  // HID
    else if (has_HE2)            variant_idx = 1;  // HIE
    else if (has_any_H)          resolved = true;  // ambiguous; leave -1
}
else if (res.type == AminoAcid::ASP) {
    bool has_HD2 = name_to_idx.find("HD2") != name_to_idx.end();   // STRING
    if (has_HD2)            variant_idx = 0;  // ASH
    else if (has_any_H)     resolved = true;  // charged ASP, leave -1
}
else if (res.type == AminoAcid::GLU) {
    bool has_HE2 = name_to_idx.find("HE2") != name_to_idx.end();   // STRING
    /* same shape */
}
else if (res.type == AminoAcid::CYS) {
    auto sg_it = name_to_idx.find("SG");                            // STRING
    if (sg_it != name_to_idx.end() &&
        disulfide_sg.count(sg_it->second) > 0) {
        variant_idx = 0;  // CYX
        resolved = true;
    } else {
        bool has_HG = name_to_idx.find("HG") != name_to_idx.end();  // STRING
        if (has_HG || has_any_H) resolved = true;  // free CYS default
    }
}
else if (res.type == AminoAcid::LYS) {
    bool has_HZ1 = name_to_idx.find("HZ1") != name_to_idx.end();    // STRING
    bool has_HZ2 = name_to_idx.find("HZ2") != name_to_idx.end();    // STRING
    bool has_HZ3 = name_to_idx.find("HZ3") != name_to_idx.end();    // STRING
    if (has_HZ1 && has_HZ2 && has_HZ3) resolved = true;  // charged LYS
    else if (has_HZ1 || has_HZ2)       variant_idx = 0;  // LYN
    else if (has_any_H)                resolved = true;
}
else if (res.type == AminoAcid::ARG) {
    resolved = true;  // charged ARG default; ARN never inferred
}
else if (res.type == AminoAcid::TYR) {
    bool has_HH = name_to_idx.find("HH") != name_to_idx.end();      // STRING
    if (!has_HH && has_any_H)  variant_idx = 0;  // TYM
    else if (has_HH)           resolved = true;  // neutral TYR
}
```

This is `pdb_atom_name == "HD1"` cosplay as load housekeeping. The
function is ~110 lines and contains nine string-equality decisions
that each *encode chemistry*: the meaning "this residue is HID" lives
in the C++ literal `"HD1"` matched against the atom's PDB name. The
chemistry encoding is fragile against any loader that produces a
slightly off-canonical name (CHARMM's HSD/HSE → HID/HIE alias is
already handled in NamingRegistry; future PDB-with-explicit-charge,
non-AMBER force-field labels, or a CCD-style CCD atom_id like `H51`
for the Trp HE1 would slip through). When the function falls through
without resolving (a HIS with no H, a CYS with no SG bond and no HG),
the variant index stays -1, and a downstream consumer (the AMBER
LEAP emitter at `src/AmberLeapInput.cpp:29-32`) falls back to default
**HIE for ambiguous HIS** — embedding an opinion about the canonical
neutral tautomer at the file-format-emit boundary.

The same name-pattern read is duplicated *as a fallback* inside
`Protein::DetectAromaticRings` at `src/Protein.cpp:484-497`. That site
is gated on `res.protonation_variant_index < 0` so it only fires when
`ResolveProtonationStates` failed to assign a variant — but it
re-implements the same string-cookbook for HIS to choose
`HidImidazoleRing` vs `HieImidazoleRing` vs `HisImidazoleRing`. Two
sites, same fragile encoding.

There is also a near-twin in `FullSystemReader.cpp:74-86`
(`DetectHisVariantFromAtoms`), used for CHARMM-source TPR loads as a
final fallback when the GROMACS readback block didn't supply a
variant. Same string list (`HD1`, `HE2`), three sites total in the
codebase encoding the HID/HIE/HIP-by-H-presence rule.

The audit (`mechanical-identity-model-and-audit-2026-05-05.md` §2.3
hotspot 1) classifies this as **Category E (internal model use that
should not have been allowed)** at HIGH risk, the highest classification
in the audit's rubric. The function's documentation comment at
`src/Protein.cpp:299-304` calls itself "construction-boundary string
interpretation" — but the function runs at every protein construction
on every load path, not at a one-shot import boundary. The boundary
comment marks intent, not status.

---

## 3. Where typed information for this question already lives in the model

The pieces required to answer the chemistry question without reading
strings already exist in the substrate. Listed by where they sit:

### 3.1 The variant table (`src/AminoAcidType.cpp::AMINO_ACID_TYPES`)

Each `AminoAcidType` entry carries a `variants` field of type
`std::vector<AmberAminoAcidVariant>`. Every entry has the canonical
fields:

```cpp
struct AmberAminoAcidVariant {
    const char* name;              // "HID", "HIE", "HIP", "ASH", ...
    const char* description;       // "Nd-protonated (delta)"
    int         formal_charge;     // 0 (HID), +1 (HIP), -1 (CYM), ...
    const char* terminal_state;    // "delta", "epsilon", "doubly", ...
};
```

The variant ordering is asserted at startup by `ValidateVariantIndices()`
(`src/AminoAcidType.cpp:241-268`); a swap is a `std::abort`. Variant
indices 0/1/2 for HIS map to HID/HIE/HIP exactly because that contract
is enforced. This means: once a residue has `protonation_variant_index`,
every chemistry question about the variant ("what's the formal
charge?", "what does AMBER call this?", "how is this tautomer
described?") is a typed lookup, not a string match.

### 3.2 The typed pointer field (`Residue.protonation_variant_index`)

`int protonation_variant_index = -1` plus
`bool protonation_state_resolved = false` (`src/Residue.h:43-44`). The
`int` is effectively a typed enum index; -1 means "default-charged
form, no variant" (so LYS with -1 means charged ammonium; ARG with -1
means charged guanidinium; the ASP/GLU/CYS/TYR default; HIS with -1
means *truly ambiguous*). The boolean is a tri-state bookkeeping flag.

This pair of fields IS the typed chemistry answer surface. The
question is who writes it and when.

### 3.3 The protonation result (`src/ProtonationDetectionResult.{h,cpp}`)

A `ConformationResult` that *reports* the resolved state. Per the
header (`ProtonationDetectionResult.h:6-9`):

> Protonation/variant identity is resolved during Protein construction.
> This result remains as compatibility/reporting plumbing for the existing
> ConformationResult surface; it does not perform atom-name lookup.

The implementation (`ProtonationDetectionResult.cpp:20-43`) is a clean
read of `res.protonation_variant_index` and `res.protonation_state_resolved`
into a per-residue variant name string. It does not look at atom
names. The correct boundary is already drawn by this result type;
`ResolveProtonationStates` is the producer that violates it.

### 3.4 The substrate variant tables (landed 2026-05-05 at `ee1f1b4`)

`src/generated/LegacyAmberSemanticTables.{h,cpp}` carries per-(residue,
variant) tables. HIS has three: `kHisAtoms_HID`, `kHisAtoms_HIE`,
`kHisAtoms_HIP`, indexed at runtime by `LookupBy(AminoAcid::HIS,
variant_idx, identity)` per §H.5 of
`spec/plan/topology-encoding-dependencies-2026-05-05.md`. CYS has
`kCysAtoms_CYS`, `kCysAtoms_CYX`, `kCysAtoms_CYM`. ASP has
`kAspAtoms_ASP`, `kAspAtoms_ASH`. And so on. The tables are correctly
indexed by `variant_idx` — once variant_idx is set typed, all chemistry
follows from `LookupBy`. The lookup key is `AtomMechanicalIdentity`
(`src/SemanticEnums.h:848-863`), a 5-tuple of typed enums; no string
flows through.

### 3.5 The upstream pKa tools (`PropkaProtonator`, `KamlProtonator`)

Both tools (`src/PropkaProtonator.cpp::ApplyHendersonHasselbalch` and
`src/KamlProtonator.cpp::Protonate`) directly produce a typed
`ResidueProtonation` decision with `variant_index` set per residue:

```cpp
// PropkaProtonator.cpp:174-179, HIS branch
case AminoAcid::HIS:
    decision.variant_index = protonated ? 2 : 1;    // HIP if charged, HIE if neutral
    decision.is_charged    = protonated;
    break;
```

Both tools resolve protonation BEFORE the Protein is built, write a
`ProtonationState`, and (if wired in to `BuildFromPdb`, which the
TODO at `PropkaProtonator.cpp:245` notes is not yet done) would set
`Residue.protonation_variant_index` directly. The chemistry decision
is already typed at the tool's output. The model gap: the tools'
output is not currently propagated to `Residue.protonation_variant_index`
at construction time on the PDB load path.

### 3.6 The GROMACS readback block (landed 2026-05-02)

`src/GromacsToAmberReadbackBlock.{h,cpp}` reads `topol.top` rtp
comments produced by `pdb2gmx`, resolves canonical 3-letter codes +
variant indices for every residue in the trajectory, and exposes them
to `FullSystemReader::BuildProtein`. The trajectory load path
(`FullSystemReader.cpp:740-748`) sets `res.protonation_variant_index =
variant_index` from the readback **before** `FinalizeConstruction`
runs. On this path, `ResolveProtonationStates`'s name-based logic is
overridden — the function still runs and still reads names, but the
typed value from upstream takes precedence (per the early-out at
`src/Protein.cpp:321-323`).

### 3.7 The ORCA load path

`src/OrcaRunLoader.cpp:189-198` reads the AMBER-style residue label
from the ORCA run's metadata (e.g., `HID`, `ASH`, `CYX`) and matches
it against `aatype.variants[].name`. If the residue label is a variant
name, `res.protonation_variant_index = vi` is set typed at load. If
the label is the canonical 3-letter code (HIS, CYS, ASP, ...), the
field stays -1 and `ResolveProtonationStates`'s name-based fallback
runs.

### 3.8 What is missing

- **PDB load**: no typed handoff. `BuildFromPdb` calls `reduce` to
  add hydrogens, then `ParsePdb` builds the Protein with all
  `protonation_variant_index = -1`, then `ResolveProtonationStates`
  runs and reads names. This is the only load path where the function's
  name-based logic does the work in production. PROPKA/KaML are
  available but not wired in (`TODO` at
  `PropkaProtonator.cpp:245` and `KamlProtonator.cpp:120`).
- **A typed identity-presence query on Residue**. The substrate's
  `LookupBy(residue, variant_idx, identity)` needs an identity to
  query; identity is currently only in the substrate tables. The
  per-atom `Atom.identity` of type `AtomMechanicalIdentity` proposed
  in the audit (Phase 1) is not yet on `Atom`; without it, the
  typed-identity-presence query has no key to scan. This is the
  single substrate gap. (Once it lands, the typed query
  `res.HasAtomWithIdentity({Element::H, Locant::Delta,
  BranchAddress{1,0}, DiastereotopicIndex::None, BackboneRole::None})`
  returns true if HD1 is present, with no string contact.)

---

## 4. Per-load-path analysis

### 4.1 PDB load (`PdbFileReader::BuildFromPdb`)

**How protonation resolves today.** `BuildFromPdb` runs `reduce` to
protonate the heavy-atom PDB (`PdbFileReader.cpp:192-197`), then
`ParsePdb` to read the protonated PDB and build a `Protein` with
`Residue.protonation_variant_index = -1` for every residue (the field
is never written on this path; default -1). `FinalizeConstruction` then
runs `ResolveProtonationStates`, which reads `pdb_atom_name` strings
to fill the variant index by name pattern. The `BuildContext` records
`protonation_tool = "reduce"` (`PdbFileReader.cpp:210`).

**The model gap.** `reduce` is deterministic about which Hs to add,
but its decisions are recorded only as PDB atom names, not as a
typed variant index. The PDB load path therefore *cannot* skip the
name-based inference unless something else writes the variant index.
The two natural authorities (PROPKA, KaML) are present in source but
not wired in.

**Phase 1 path.** Two layered approaches, both viable:

(a) *Move the inference into the loader, but make it identity-typed.*
After `ParsePdb` constructs atoms (each carrying its identity-typed
`Atom.identity` per the audit's Phase 1 step 1), iterate residues and
use typed-identity-presence to set `protonation_variant_index`. The
chemistry encoding (HD1 → HID, HE2 → HIE, etc.) lives in one place,
inside the loader, expressed as identity literals not strings. Example
shape (illustrative — no code commitment):

```cpp
// In PdbFileReader, after ParsePdb, before FinalizeConstruction:
for (auto& res : protein->MutableResidues()) {
    if (res.type == AminoAcid::HIS) {
        bool has_HD1 = res.HasAtomWithIdentity(MakeIdHd1());
        bool has_HE2 = res.HasAtomWithIdentity(MakeIdHe2());
        if (has_HD1 && has_HE2)      res.protonation_variant_index = 2;
        else if (has_HD1)            res.protonation_variant_index = 0;
        else if (has_HE2)            res.protonation_variant_index = 1;
    }
    // ... ASP/GLU/CYS/LYS/ARG/TYR similarly
}
```

The loader's role is now: "produce a typed variant index from the
load source's evidence, in the load source's own context." The
`MakeIdHd1()` etc. helpers are typed-identity literals defined next
to the loader. No string flows through — the chemistry encoding has
moved from "string equals a literal" to "typed identity equals a
typed literal", which is meaningfully different: the typed-identity
form requires the substrate's identity vocabulary to know what HD1
*is*, not just what its PDB column-13-16 spelling is. The chemistry
is captured by the typed fields (Element::H + Locant::Delta +
BranchAddress{1, 0} + ...), not by the spelling.

(b) *Run PROPKA at load time and set variant_index from its output.*
Wire `PropkaProtonator::Protonate` into `BuildFromPdb` between
`reduce` and `FinalizeConstruction`. The tool already produces a
typed `ResidueProtonation::variant_index`; copy it onto
`Residue.protonation_variant_index`. The Henderson-Hasselbalch decision
becomes the authoritative load-time signal; reduce's H-placement is
decorative. This path skips `ResolveProtonationStates`'s inference
entirely. PROPKA is not deterministic for the HID/HIE choice (it
returns a single pKa per HIS; the tautomer is decided downstream — see
the comment at `PropkaProtonator.cpp:177-179` defaulting to HIE),
which means (b) on its own does not handle HID/HIE; combining (a) and
(b) — PROPKA for the protonated/deprotonated decision plus
identity-presence for tautomer disambiguation when reduce has chosen
a tautomer — is the most informative shape.

The user has explicit guidance at memory entry
`feedback_no_pdbfixer_mid_pipeline`: deterministic geometry or
refuse-and-document, no silent atom-moving. (a) is deterministic
(no atom motion); (b) brings PROPKA's pKa prediction into scope but
doesn't move atoms. Both are admissible; (b) requires the user to
bless the wiring of an external tool that's currently TODO'd out.
**Phase 1 baseline recommendation: (a).** PROPKA wiring (b) is its
own design decision and lives outside this Phase 1 minimum.

### 4.2 AMBER trajectory load (`FullSystemReader::BuildProtein` with readback)

**How protonation resolves today.** The TPR's residue names are
CHARMM-canonical port labels (HIS / HID / HIE / HSD / HSE / ...);
`pdb2gmx` records the full chemistry-resolution trail in `topol.top`
as rtp comments. The `GromacsToAmberReadbackBlock` parses those
comments and supplies `(canonical_three, variant_index)` per residue.
`ResolveResidueTypeAndVariant` (`FullSystemReader.cpp:666-749`) then
sets `res.type` and `res.protonation_variant_index` directly from the
readback. The HIS-by-H-name fallback `DetectHisVariantFromAtoms`
runs *only when* the readback didn't resolve the variant
(`FullSystemReader.cpp:860-864`).

**The model gap on this path.** None for the readback case. When the
readback is present (the production trajectory path post-2026-05-02),
this path already produces a typed variant index from the upstream
authority and `ResolveProtonationStates`'s subsequent run is
no-op-by-early-out at `src/Protein.cpp:321-323` (variant_index already
≥ 0). For the no-readback fallback (legacy CHARMM / non-pdb2gmx TPR),
`DetectHisVariantFromAtoms` runs, which is the same string-list
encoding as `ResolveProtonationStates`. That fallback can be
identity-typed identically to (a) above, and probably should be —
its three string-list encodings (in `Protein.cpp`, `Protein.cpp::DetectAromaticRings`,
and `FullSystemReader.cpp::DetectHisVariantFromAtoms`) are the
"three sites for the same chemistry rule" pattern PATTERNS.md prohibits.

**Phase 1 path.** With readback: nothing changes. Without readback:
typed-identity-presence as in (a). The function
`DetectHisVariantFromAtoms` becomes one call to
`res.HasAtomWithIdentity` per HIS variant, and the same shape
generalises to all titratable residues if the TPR-without-readback
path ever needs a complete fallback (it does not today; the readback
is required by trajectory production).

### 4.3 ORCA run load (`OrcaRunLoader::BuildFromOrcaRun`)

**How protonation resolves today.** ORCA run metadata carries the
AMBER residue label (HID, ASH, CYX, ...); `OrcaRunLoader.cpp:189-198`
matches it against the typed variant table by string equality (`if
(res_labels[ri] == aat.variants[vi].name)`) and sets
`res.protonation_variant_index = vi`. This is a string match against a
project-controlled typed-table contents; it's a Category-B (boundary
projection) string contact, not Category-E. When a protein label is
the canonical 3-letter form (HIS, no specific tautomer), the field
stays -1 and `ResolveProtonationStates`'s string-pattern read kicks in.

**The model gap on this path.** When the ORCA run was produced from
an AMBER-prepared structure (the standard project pipeline), every
residue has its variant label in the metadata, so the variant index
is set typed at load. The string match at line 193 is against
project-controlled literals (the variants table itself); replacing it
with `VariantIndexFromForceFieldName(res.type, res_labels[ri])`
(`AminoAcidType.cpp:309-338`, the typed helper added 2026-05-02 to fix
the wrong-indices bug) is a natural cleanup — this is, strictly, a
typed function call wrapping the string match, but it isolates the
string contact to one helper and per-variant tables. The function
works; it's the right shape.

When the ORCA load is from a 3-letter-code structure, the variant
index stays -1 and the same hotspot logic applies. Phase 1 path: same
as (a).

### 4.4 Upstream tools (PROPKA / KaML) outputs

These tools' `ApplyHendersonHasselbalch` already produces typed
`variant_index`. The "load path" for them is "the Protein has been
built; run the tool; copy variant_index onto each Residue." Today
they are run as standalone callers, not wired into `BuildFromPdb`.
Wiring them in is a Phase 1 *extension*, not a Phase 1 baseline. The
shape exists in source; it is a load-path option, not a missing
component.

### 4.5 Mutant (`--mutant`) and protonated-PDB (`--protonated-pdb`) paths

These reuse `BuildFromPdb` (--mutant) or read a user-supplied PDB
where the protonation variant is encoded in the residue label
(--protonated-pdb). The --protonated-pdb path is a string-matched
ORCA-style label-detection; the --mutant path inherits PDB load.
Phase 1 treatment is identical to the PDB load path.

### 4.6 Summary table

| Load path | Today | Phase 1 path |
|---|---|---|
| PDB (BuildFromPdb) | string-pattern read on H names in `ResolveProtonationStates` | identity-typed presence test in PdbFileReader, before FinalizeConstruction |
| AMBER trajectory with readback | typed variant_index from `GromacsToAmberReadbackBlock` | unchanged (already correct) |
| AMBER trajectory without readback (legacy) | string read via `DetectHisVariantFromAtoms` fallback in FullSystemReader, plus `ResolveProtonationStates` | identity-typed presence in `DetectHisVariantFromAtoms` reformulated as typed query |
| ORCA run | residue label string matched against typed variant table | replace string match with `VariantIndexFromForceFieldName` (already exists; it's a one-line cleanup) |
| --protonated-pdb | residue label string matched | replace string match with `VariantIndexFromForceFieldName` |
| --mutant | inherits PDB | inherits PDB |
| PROPKA / KaML output | typed variant_index produced; not wired into load | optional: wire into PdbFileReader as protonation source |

By the time `Protein::FinalizeConstruction` runs, every load path has
populated `Residue.protonation_variant_index` from its own typed
evidence. `ResolveProtonationStates` becomes a *consistency check*
that asserts the variant index is set for every titratable residue
that should have one (and triggers a diagnostic for any HIS that
truly has no H, the genuine ambiguous case).

---

## 5. Phase 1 design (existing calculators must work)

Sequencing within the audit's seven-phase pivot
(`mechanical-identity-model-and-audit-2026-05-05.md` §2.4): protonation
work depends on the substrate identity surface (Phase 1 of the audit's
sequence — `Atom.identity` populated at every loader). Once that
exists, this protonation refactor is a Phase-2-style task confined to
the `ResolveProtonationStates` site plus the four per-load-path
producers described above. The full protonation cleanup is plausibly
~80 lines of source change, distributed.

### 5.1 What stays

- The contract `Residue.protonation_variant_index` (int, -1 = default
  charged form for residue type, ≥ 0 = variant index per
  `AmberAminoAcidVariantTable`) is unchanged. Existing calculator
  code that consumes it (`AmberLeapInput.cpp`, `AmberChargeResolver.cpp`,
  `ChargeSource.cpp`) is unaffected.
- `Residue.protonation_state_resolved` (bool) stays as the
  bookkeeping flag.
- The `AmberAminoAcidVariantTable` ordering contract enforced by
  `ValidateVariantIndices()` stays. Variant indices 0/1/2 still mean
  HID/HIE/HIP for HIS, etc.
- `ProtonationDetectionResult` stays unchanged. It already only reads
  the typed variant index; it does no name lookup.
- The `HidImidazoleRing` / `HieImidazoleRing` / `HisImidazoleRing` typed
  ring classes stay. The ring-type selection in
  `Protein::DetectAromaticRings` already prefers
  `protonation_variant_index` (the typed branch at
  `Protein.cpp:472-482`); the name-pattern fallback at
  `Protein.cpp:484-497` becomes dead code under the contract that
  variant_index is always set by the time `DetectAromaticRings` runs.
- The AMBER LEAP emitter's HIS-default-to-HIE rule
  (`AmberLeapInput.cpp:29-32`) stays as a *legitimate* fallback for
  the rare case where a HIS truly has no H (crystal-only structure
  with all hydrogens stripped). This is a file-format-emit decision,
  not a chemistry decision the substrate is responsible for; keeping
  it at the emit boundary is the right shape per the audit's
  Crystal Projection Rule.

### 5.2 What changes

**Per-load-path variant_index population.** Each loader becomes
responsible for setting `Residue.protonation_variant_index` from its
own evidence:

- **PdbFileReader::BuildFromPdb**: after `ParsePdb`, before
  `FinalizeConstruction`, run a typed-identity-presence pass that
  encodes the same chemistry rules `ResolveProtonationStates` did
  (HD1 → HID, HE2 → HIE, both → HIP, HD2-on-ASP → ASH, HE2-on-GLU →
  GLH, no-HG-on-CYS-with-no-disulfide → CYM, no-Hζ-on-LYS → LYN,
  no-HH-on-TYR → TYM). The chemistry rules live as identity literals
  next to the loader, not as string literals inside `Protein.cpp`.
- **FullSystemReader (no-readback fallback)**: same pass for HIS in
  `DetectHisVariantFromAtoms`, reformulated as identity-typed.
- **OrcaRunLoader**: replace the string match at
  `OrcaRunLoader.cpp:193` with
  `VariantIndexFromForceFieldName(res.type, res_labels[ri])` (already
  exists in `AminoAcidType.cpp`).
- **--protonated-pdb path**: same as ORCA.
- **PROPKA / KaML wiring** (optional): if BuildFromPdb is enhanced to
  call PROPKA/KaML, copy the tool's `variant_index` onto the residue
  before the identity-typed pass; the identity pass becomes a
  consistency check (PROPKA chose neutral HIS, identity sees HD1
  present → HID is the typed answer).

**Protein::ResolveProtonationStates becomes a consistency check.**
The function still runs at `FinalizeConstruction` (twice, per the
current shape — once before bond detection, once after). Its new body:

```text
For each residue:
    If aatype.is_titratable and aatype.variants is non-empty:
        If protonation_variant_index < 0:
            // For CYS specifically, the disulfide-bond check happens
            // here because CYX detection requires bonds, which is
            // only available after CovalentTopology::Resolve. The
            // first ResolveProtonationStates call (use_covalent=false)
            // does not run this check; the second (use_covalent=true)
            // does, and may bump CYS to CYX = variant 0.
            If type == CYS and use_covalent_topology and SG-is-disulfide-bonded:
                protonation_variant_index = 0  // CYX
                protonation_state_resolved = true
            Else:
                // Diagnostic: this residue's load source did not
                // resolve a variant. Log it; leave -1.
                Log("residue " + ... + ": protonation_variant_index unresolved")
        Set protonation_state_resolved = true
    Else if not titratable:
        Set protonation_state_resolved = true (vacuously)
```

The CYX-from-disulfide branch stays inside the function because it
genuinely requires bond detection to have run first; that's a
consistency result that depends on geometry, and the function's
two-pass invocation pattern (before/after CovalentTopology) is the
right shape for it. Everything else moves out.

**The HIS fallback in DetectAromaticRings becomes dead code.** The
name-pattern read at `Protein.cpp:484-497` is gated on
`res.protonation_variant_index < 0`, which is no longer reachable
under the new contract (every load path sets the variant index, or
flags ambiguous explicitly). Remove the fallback; replace with an
assertion that variant_index is always ≥ 0 for HIS by the time
`DetectAromaticRings` runs, OR fall through to a project-policy
default (HIE, matching the AMBER LEAP emitter's choice).

### 5.3 Acceptance criteria for Phase 1

- All existing tests pass. The chemistry decisions (which residue is
  which variant) on every fixture are identical before and after.
- The `pdb_atom_name` reads inside `Protein::ResolveProtonationStates`
  and `Protein::DetectAromaticRings` are gone (or, for the
  CYX-from-disulfide branch, the read is replaced by the typed
  `LegacyAmber().BondCategoryFor` query).
- `DetectHisVariantFromAtoms` in FullSystemReader is reformulated as
  identity-typed.
- `ProtonationDetectionResult` reports identical variant names per
  residue across the fixture set (the test
  `tests/test_protonation_detection.cpp` if it exists, or a new one
  if not).
- The AMBER LEAP emitter's HIS-default-to-HIE rule continues to fire
  for crystal-only-no-H structures.
- A new diagnostic surface in the consistency-check role: any
  titratable residue with `protonation_variant_index < 0` after
  load is logged with sufficient context for triage.

### 5.4 What this Phase 1 design does not solve

- The `is_disulfide_cys` argument threading in `AmberLeapInput`. The
  emit-side decision "is this CYS in a disulfide?" is currently fed
  from `DetectDisulfides` (`AmberLeapInput.cpp:315-348`), which uses
  geometric SG-SG distance — this is its own audit hotspot 2 (SG
  detection by name). Phase 1 of protonation does not fix this; the
  CYX variant_index propagation is independent (set at TPR-readback
  time on trajectory loads; set by SG-bond-presence in the consistency
  check on PDB loads). The AMBER LEAP emitter remains correct because
  it has a separate disulfide path; cleaning up DetectDisulfides is
  audit hotspot 2's job.
- The HIS-no-H ambiguity. A crystal structure with no Hs on HIS is
  genuinely ambiguous; the model has no way to choose HID/HIE/HIP
  without external evidence. Today the AMBER LEAP emitter defaults
  to HIE; Phase 1 keeps that. A more principled answer (PROPKA, or
  pKa-by-environment) is Phase 2 territory.

---

## 6. Phase 2 considerations

The user's framing was explicit: do not solve Phase 2 now. These are
flags for the future calculator-walk pass.

### 6.1 Time-dependent protonation

The current model treats protonation as substrate-invariant: one
variant_index per residue for the lifetime of the Protein. This is
correct under fixed-protonation MD (the standard ff14SB protocol);
constant-pH MD or QM/MM with explicit proton transfer would violate
it. If the project's calculator set ever adopts constant-pH
trajectories, the substrate's invariant variant_index becomes wrong;
the chemistry would migrate to a per-frame `ProtonationStateResult`
ConformationResult and the typed variant index becomes
ProteinConformation-scope.

### 6.2 Tautomer ambiguity in crystal-only structures

A HIS or LYS deposited without explicit hydrogens (crystal at
moderate resolution) carries no information about its tautomer.
Today's pipeline forces a default (HIE) at the AMBER LEAP emitter;
this is a research limitation but a known one. A future pass could
record "ambiguous" as a typed value and propagate it to calculators
that should refuse to compute (e.g., the variant-dependent ring-current
calculator could emit a `KernelEvaluationFilter` rejection for
ambiguous HIS rather than computing on the default tautomer).

### 6.3 pKa-environment coupling under conformational change

PROPKA and KaML are run once per protein at load. If a residue's
pKa shifts during MD because its electrostatic environment changes
(solvent exposure, salt-bridge formation, conformational rearrangement),
the variant_index assigned at t=0 may be wrong at t=N. The substrate's
single-variant-per-residue contract precludes this from being modelled.
A future pass could expose the pKa as a per-frame computation
(`PkaResult` ConformationResult) that the analysis layer correlates
with shielding.

### 6.4 Aromatic ring chemistry depends on HIS variant

The three HIS variants produce three distinct typed Ring subclasses
(`HidImidazoleRing`, `HieImidazoleRing`, `HisImidazoleRing` for HIP).
The ring-current intensity, lobe offset, and N-protonation pattern
differ. Today's calculators (`BiotSavartResult`, `HaighMallionResult`,
`RingSusceptibilityResult`) use the `Ring::Intensity()` virtual; the
typed dispatch is correct. A Phase 2 walk could consider whether
the variant_index should also drive the `PlanarGroupKind` enum
(`Imidazole` vs new `ImidazoliumPositive` for HIP) and the
`PolarHKind::ImidazoleNH` site selection. The substrate already
encodes per-(residue, variant) tables, so this is a substrate-table
extension rather than a runtime change.

### 6.5 ARN, TYM, CYM as research targets

ARN is uncommon at physiological pH and lacks a canonical ff14SB
template (per `topology-encoding-dependencies-2026-05-05.md` §B.1).
TYM and CYM are accessible but rare. A future pass could examine
whether the project's calibration data set has any examples; if so,
the variant tables' "deprotonated" entries earn calculator coverage;
if not, they remain documentation of the chemistry without runtime
exercise.

### 6.6 Constant-pH downstream MD analysis (relaxation observables)

The 2-protein dense validation plan (memory entry
`project_dense_validation_2proteins_20260424`) explicitly targets
relaxation observables (S², T1, T2, CCR, J(ω)) calibrated against
published data. Some observable interpretations depend on whether
a HIS is HID or HIE (the H-bond donor/acceptor pattern differs).
Phase 2 could add a per-residue protonation-confidence diagnostic
that the observable calculators consult to gate their output.

---

## 7. Choices that need a user decision

### 7.1 Phase 1 baseline: identity-typed inference vs PROPKA wiring

(a) Identity-typed-presence inference inside the loader, encoded as
typed identity literals. Direct, deterministic, no new external tool
dependency, no chemistry change.

(b) Wire PROPKA (and/or KaML) into BuildFromPdb. Brings pKa prediction
into the load path; the H-H decision becomes the authoritative
protonation source. Adds an external-tool dependency. The TODO at
`PropkaProtonator.cpp:245` is currently outside the build wiring.

(c) Both. Use PROPKA for the protonated/deprotonated decision; use
identity-typed presence to disambiguate HID/HIE when reduce has
chosen a tautomer.

**Recommendation: (a) as Phase 1 baseline.** (b) and (c) are research
extensions worth scheduling but should not be inside the substrate
refactor.

### 7.2 What does -1 mean post-Phase 1?

Today: -1 means "default charged form" (LYS = +1, ARG = +1, ASP/GLU/CYS
= charged) AND "ambiguous HIS" depending on residue type. The semantic
overload is fragile.

Two options:
- Keep -1 as polymorphic (residue-type determines its meaning).
- Introduce a typed enum `ProtonationVariant` per residue family,
  with explicit "default" and "ambiguous" values. The `int` field
  becomes a typed enum with named values per residue.

**Recommendation:** keep `int` for Phase 1 (no contract change),
flag this as Phase 2 cleanup. The current contract is correct as
long as ambiguous-HIS is logged at the consistency-check site.

### 7.3 Should ResolveProtonationStates die or survive as consistency check?

(a) Delete the function entirely; loaders are fully responsible.
(b) Keep as a consistency-check + CYX-from-bond-detection site.

The CYX-from-bond-detection branch genuinely needs to run after
CovalentTopology resolves bonds; it can't move into the loader (the
loader doesn't have bonds yet). So (b) is the structurally correct
answer; (a) would force the CYX branch into a separate function with
no obvious home.

**Recommendation: (b).** The function survives but its body shrinks
to ~10 lines: the CYX-from-bond branch plus a per-titratable-residue
diagnostic for unresolved variant indices.

### 7.4 What happens to FullSystemReader::DetectHisVariantFromAtoms?

The function is the no-readback fallback for HIS variant detection
on TPR loads. With production trajectory loads always using a readback,
the function only fires on legacy TPRs.

(a) Keep as identity-typed-presence logic for the legacy path.
(b) Delete; force readback as a hard requirement.

**Recommendation: (a).** The CHARMM-retired memo
(`project_charmm_retired_amber_only_2026-05-02`) makes the
no-readback path quarantined-legacy, but the function might run for
non-pdb2gmx-produced TPRs in the future. Refactoring to identity-typed
takes ~10 lines and removes the third site of the same string-pattern
encoding; keeping it as identity-typed is cheap insurance.

### 7.5 Should the AMBER LEAP emitter's HIS-default-to-HIE be moved?

`AmberLeapInput.cpp:29-32` defaults `AmberResidueNameFor` to "HIE"
when `protonation_variant_index < 0` and the residue type is HIS.
This is a wire-format default at the file-format-emit boundary.

(a) Leave it. The Crystal Projection Rule
(`amber-implementation-plan-2026-04-29.md`) says names are computed
at the emit site, and the default-tautomer choice IS a name decision.
(b) Promote it to a substrate decision: when the consistency check
finds HIS with variant_index < 0, set it to 1 (HIE) and log.

**Recommendation: (a).** Defaults at the file-format boundary are
the right shape per the Crystal Projection Rule. Moving it to the
substrate would imply HIE is the chemistry truth, which it is not
(it's an opinion about the canonical neutral tautomer). Keeping the
opinion at the AMBER-emitter site is honest.

### 7.6 Identity literals: how do they look in calculator code?

The audit notes (closing notes §1) that constructing a full
`AtomMechanicalIdentity` literal in calculator code is verbose
(5 fields, most None for any specific query). The audit recommends
partial-identity accessors on Residue:

```cpp
std::optional<size_t> AtomWithRole(BackboneRole) const;
std::optional<size_t> AtomWithLocant(Locant, BranchAddress = {}) const;
bool HasAtomWithIdentity(const AtomMechanicalIdentity&) const;
```

The protonation Phase 1 work needs `HasAtomWithIdentity` for the
HD1/HE2/HD2/HE2/HZ-presence tests, plus a
`HasAtomWithLocantAndElement(Element, Locant, BranchAddress)` helper
for the (Element::H, Locant::Delta, BranchAddress{1, 0}) query that
distinguishes HD1 from HD2 (different BranchAddress).

**This decision is upstream of protonation Phase 1 — it's part of the
audit's Phase 1 (substrate identity surface).** The protonation
refactor needs to know the partial-identity accessor shape to write
its identity literals. The user should bless the audit's
recommendation (or an alternative shape) before protonation Phase 1
lands.

---

## 8. Risks and unknowns

### 8.1 What could break in Phase 1?

- **Silent variant changes on edge-case fixtures.** Any fixture where
  the current name-pattern read produces a different answer than the
  identity-typed read would change behaviour. The chemistry should be
  identical (identity-typed encodes the same rule), but a fixture
  with a non-canonical PDB atom name (e.g., HSD instead of HID,
  or HE21 misspelled as HE12) might have been silently mis-classified
  by either path. A bit-identical-output gate on the existing
  fixture set is the test for this; the user's mutant-validation gate
  (memory `project_mutant_validation_gate_20260426`) extends it.
- **Two-pass invocation timing**. `ResolveProtonationStates` runs
  twice: before and after `CovalentTopology::Resolve`. The first
  call has no bonds (use_covalent_topology = false); the second
  does. The CYX-from-disulfide branch only fires on the second
  call. If protonation moves into loaders, the first-call work is
  no longer needed; only the CYX-from-bond second-call work
  remains. Verifying that nothing else in `FinalizeConstruction`
  depends on the first-call's variant_index resolution is required;
  in particular, `DetectAromaticRings` runs *between* the two calls
  (`Protein.cpp:220` after first ResolveProtonationStates, before
  CovalentTopology::Resolve at `:225`). DetectAromaticRings consumes
  variant_index for HIS ring-type selection; the loader-side variant
  setting must happen before DetectAromaticRings runs.
- **The `protonation_state_resolved` boolean's semantics**. Today
  it tracks "we tried to resolve this", which collapses three cases
  ("variant_index assigned typed", "variant_index left at -1 by
  policy default", "variant_index left at -1 because ambiguous").
  Phase 1's diagnostic surface needs to disambiguate; the boolean
  may need a typed-enum replacement (Resolved / DefaultByPolicy /
  Ambiguous) for the diagnostic to be useful.
- **CHARMM CHARMM-source fixtures**. The TPR-without-readback
  fallback in `DetectHisVariantFromAtoms` is exercised by old test
  fixtures (the quarantined CHARMM path, `tests/bones/`). If any
  active test still goes through that path, the identity-typed
  reformulation needs the same fixture coverage.

### 8.2 Edge cases not yet covered

- **Crystal HIS with no Hs.** Today: silent -1, AMBER LEAP picks HIE.
  Post-Phase-1: same behaviour, but the diagnostic logs the
  ambiguity. The chemistry hasn't improved; the visibility has.
- **Mixed protonation across copies (e.g., monomer A is HID, monomer
  B is HIE).** Today: each residue resolved independently from its
  own atoms. Same in Phase 1.
- **Non-standard residues (modified amino acids, ligands, cofactors).**
  The substrate is built for the standard 20. A non-standard residue
  with `is_titratable = false` is silent; one with `is_titratable =
  true` and no variants table entry would crash the consistency
  check. Phase 1 doesn't introduce non-standard residues; this is
  inert until the project does.
- **N-terminal and C-terminal protonation states.** `ResidueTerminalState`
  is independent of `protonation_variant_index`. The terminal state
  has its own typed enum (Internal / NTerminus / CTerminus /
  NAndCTerminus) and is resolved at `ResolveResidueTerminalStates`
  (`Protein.cpp:273-295`) by chain-position rather than name. This
  is the right shape and Phase 1 doesn't touch it; flagged here
  because protonation discussions can drift into terminal state
  discussions.
- **Disulfide variant_index timing.** CYX is variant 0; CYM is variant
  1; the variant indices are not co-orderable with "amount of H"
  (CYX has no H, CYM has no H, CYS has Hγ). The consistency check's
  CYX-from-disulfide branch must override variant_index that may
  have been set to CYM or CYS (default) by the loader if a disulfide
  is later detected. The current code's logic at
  `Protein.cpp:374-382` already handles this; preserving it under
  the new shape requires the consistency check to take precedence
  over loader-set variant_index when the loader was wrong about
  CYX (e.g., a PDB without explicit disulfide records but with two
  CYS Sγ within bonding distance).

### 8.3 What the audit would have caught with deeper coverage

The audit named the function as hotspot 1 but did not fully enumerate
the per-load-path producer/consumer matrix. This document is that
enumeration. Specifically:

- The `OrcaRunLoader` and `--protonated-pdb` paths' string match against
  variant table entries (Category B in the audit's classification)
  was named but not flagged as a sequencing dependency on
  `VariantIndexFromForceFieldName` — that is, the cleanup is trivial
  but the audit didn't note it as a Phase 1 task.
- The PROPKA / KaML wiring TODOs are not in the audit at all. They
  may not need to be; the audit was scoped to `pdb_atom_name` sites.
  But the chemistry question around protonation inevitably invokes
  pKa prediction as the alternative authority.
- The semantic overload of `protonation_variant_index = -1`
  (default-charged vs ambiguous) is a design choice the audit didn't
  examine. Whether Phase 1 leaves the polymorphism or introduces a
  typed enum is a user decision.
- `DetectAromaticRings`'s HIS fallback at `Protein.cpp:484-497` is
  the same name-pattern logic in a third location. The audit's
  hotspot 3 names DetectAromaticRings's ring-atom-name match
  (Category E) but does not call out the HIS-protonation fallback
  embedded inside it. This third site should not survive Phase 1;
  with variant_index always set by the time DetectAromaticRings
  runs, the fallback is dead code.

### 8.4 What's not a risk but is a discipline

Per the user's framing, the chemistry encoding (HD1 → HID, HE2 → HIE,
both → HIP, HD2-on-ASP → ASH, no-Hγ-on-CYS → CYM-or-CYX, etc.) is
*physically meaningful*. Moving it from string literals to typed
identity literals does not erase the chemistry; it reifies it. The
typed identity tuple `(Element::H, Locant::Delta, BranchAddress{1,0},
DiastereotopicIndex::None, BackboneRole::None)` is the chemistry
encoding of "HD1" — same information, different vocabulary. The
goal is not to dilute the chemistry; the goal is to put it in a
vocabulary the substrate can answer questions about, instead of a
vocabulary that requires the substrate to be a string-key lookup
table. The chemistry is the same; the structure of how the code
asks about it changes.

---

## End

This document is Phase-1-design-and-Phase-2-considerations only. No
code changes have been proposed in form; the per-load-path text in §5.2
and the snippet in §4.1 are illustrative shape notes, not commit
candidates. The user owns the decisions in §7 before implementation
begins.

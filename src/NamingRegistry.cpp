#include "NamingRegistry.h"
#include <algorithm>
#include <cctype>

namespace nmr {

// ============================================================================
// Helpers
// ============================================================================

static std::string ToUpper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return s;
}

static std::string Trim(std::string s) {
    s.erase(0, s.find_first_not_of(' '));
    s.erase(s.find_last_not_of(' ') + 1);
    return s;
}

bool NamingRegistry::ContextKey::operator<(const ContextKey& o) const {
    if (canonical != o.canonical) return canonical < o.canonical;
    if (context != o.context) return context < o.context;
    return variant < o.variant;
}

bool NamingRegistry::AtomNameKey::operator<(const AtomNameKey& o) const {
    if (atom_name != o.atom_name) return atom_name < o.atom_name;
    if (residue_name != o.residue_name) return residue_name < o.residue_name;
    if (from_context != o.from_context) return from_context < o.from_context;
    return to_context < o.to_context;
}


// ============================================================================
// Construction: populate all naming tables
// ============================================================================

NamingRegistry::NamingRegistry() {
    InitialiseStandardResidues();
    InitialiseAmberContext();
    InitialiseCharmmContext();
    InitialiseAtomNameRules();
}


void NamingRegistry::InitialiseStandardResidues() {
    standard_residues_ = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    };

    // Every standard name maps to itself
    for (const auto& name : standard_residues_)
        to_canonical_[name] = name;

    // Variant -> canonical mappings (accumulated from all tool universes)
    to_canonical_["HID"] = "HIS";  to_canonical_["HIE"] = "HIS";
    to_canonical_["HIP"] = "HIS";
    to_canonical_["HSD"] = "HIS";  to_canonical_["HSE"] = "HIS";
    to_canonical_["HSP"] = "HIS";

    to_canonical_["ASH"] = "ASP";
    to_canonical_["ASPP"] = "ASP";

    to_canonical_["GLH"] = "GLU";
    to_canonical_["GLUP"] = "GLU";

    to_canonical_["CYX"] = "CYS";  to_canonical_["CYM"] = "CYS";
    to_canonical_["CYS2"] = "CYS";

    to_canonical_["LYN"] = "LYS";
    to_canonical_["ARN"] = "ARG";
    to_canonical_["TYM"] = "TYR";
    to_canonical_["MSE"] = "MET";
}


void NamingRegistry::InitialiseAmberContext() {
    auto add = [&](const std::string& canonical, const std::string& variant,
                    const std::string& amber_name) {
        context_map_[{canonical, ToolContext::Amber, variant}] = amber_name;
        // OpenMM uses the same naming as AMBER
        context_map_[{canonical, ToolContext::OpenMM, variant}] = amber_name;
    };

    // HIS tautomers
    add("HIS", "delta",   "HID");
    add("HIS", "epsilon", "HIE");
    add("HIS", "doubly",  "HIP");
    add("HIS", "",        "HIS");

    // Protonated acids
    add("ASP", "protonated", "ASH");
    add("GLU", "protonated", "GLH");

    // Cysteine states
    add("CYS", "disulfide",    "CYX");
    add("CYS", "deprotonated", "CYM");

    // Neutral lysine
    add("LYS", "deprotonated", "LYN");

    // Deprotonated arginine (extremely rare, pKa ~12.5)
    add("ARG", "deprotonated", "ARN");

    // Deprotonated tyrosine
    add("TYR", "deprotonated", "TYM");
}


void NamingRegistry::InitialiseCharmmContext() {
    auto add = [&](const std::string& canonical, const std::string& variant,
                    const std::string& charmm_name) {
        context_map_[{canonical, ToolContext::Charmm, variant}] = charmm_name;
    };

    // HIS tautomers (CHARMM uses HSD/HSE/HSP)
    add("HIS", "delta",   "HSD");
    add("HIS", "epsilon", "HSE");
    add("HIS", "doubly",  "HSP");
    add("HIS", "",        "HIS");

    // Protonated acids (CHARMM uses ASPP/GLUP)
    add("ASP", "protonated", "ASPP");
    add("GLU", "protonated", "GLUP");

    // Cysteine disulfide (CHARMM uses CYS2)
    add("CYS", "disulfide",    "CYS2");
    add("CYS", "deprotonated", "CYM");
}


void NamingRegistry::InitialiseAtomNameRules() {
    // Backbone amide hydrogen: H (PDB/AMBER) <-> HN (CHARMM)
    AddAtomNameRule("H",  "HN", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HN", "H",  "*", ToolContext::Charmm, ToolContext::Standard);

    // Beta methylenes (all residues: HB2,HB3 <-> HB1,HB2)
    AddAtomNameRule("HB2", "HB1", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HB3", "HB2", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HB1", "HB2", "*", ToolContext::Charmm, ToolContext::Standard);
    AddAtomNameRule("HB2", "HB3", "*", ToolContext::Charmm, ToolContext::Standard);

    // Gamma methylenes (all residues: HG2,HG3 <-> HG1,HG2)
    AddAtomNameRule("HG2", "HG1", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HG3", "HG2", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("HG1", "HG2", "*", ToolContext::Charmm, ToolContext::Standard);
    AddAtomNameRule("HG2", "HG3", "*", ToolContext::Charmm, ToolContext::Standard);

    // ========================================================================
    // C-terminal carboxyl atoms — CHARMM port → AMBER canonical
    // ========================================================================
    //
    // CHARMM emits the two C-terminal oxygens as OT1 / OT2; AMBER
    // ff14SB and IUPAC use O (the carbonyl O of the chain) and OXT
    // (the deprotonated carboxyl O). The AMBER kCapCtermDeprotonated
    // / kCapCtermProtonated tables key on canonical O / OXT, so
    // CHARMM-loaded fixtures need to collapse here before substrate
    // composition or the cap lookup misses entirely. Mapping:
    //
    //     CHARMM OT1  <-->  AMBER/IUPAC O   (the chain backbone carbonyl O)
    //     CHARMM OT2  <-->  AMBER/IUPAC OXT (the C-terminal carboxyl O)
    //
    // Verified against the AMBER ff14SB cap convention encoded in
    // src/generated/LegacyAmberSemanticTables.cpp's kCapCterm tables
    // and the ff14SB parameter files at data/ff14sb_params.dat.
    //
    // These rules live in the (Charmm, Standard) keyspace -- they
    // are wire-format-collapse aliases and fire only on CHARMM
    // inputs. They do NOT have (Standard, Amber) counterparts;
    // canonical Standard inputs already use O/OXT.
    AddAtomNameRule("OT1", "O",   "*", ToolContext::Charmm,   ToolContext::Standard);
    AddAtomNameRule("OT2", "OXT", "*", ToolContext::Charmm,   ToolContext::Standard);
    AddAtomNameRule("O",   "OT1", "*", ToolContext::Standard, ToolContext::Charmm);
    AddAtomNameRule("OXT", "OT2", "*", ToolContext::Standard, ToolContext::Charmm);

    // ========================================================================
    // LYN protonation-variant H-on-NZ naming — fleet-vintage variance
    // ========================================================================
    //
    // AMBER ff14SB LYN (deprotonated lysine, neutral NH2 amine on NZ)
    // has TWO H atoms on NZ canonically named HZ2 / HZ3 — preserving
    // the LYS NH3+ numbering convention with HZ1 absent (the proton
    // removed at deprotonation). Some upstream prep flows (notably
    // the 1Z9B fleet input.pdb) carry LYN-chemistry residues with
    // the H atoms named HZ1 / HZ2 (a pre-Markley-1998 numbering).
    // Map those to the canonical HZ2 / HZ3 when the residue is
    // recorded as LYN.
    //
    // RULES ARE RESIDUE-CONTEXT-KEYED. They fire ONLY when the
    // residue's three-letter code is "LYN". For LYS residues
    // (NH3+, HZ1/HZ2/HZ3) the rule does NOT fire and HZ1/HZ2/HZ3
    // pass through canonical.
    //
    // The rules live in the (Standard, Amber) keyspace so they fire
    // in Stage 2 of CanonicaliseAmberAtomName regardless of source
    // wire-format. They do NOT have (Charmm, Standard) counterparts:
    // putting them in both keyspaces would chain incorrectly
    // (HZ1 -> HZ2 -> HZ3 in a single canonicalisation pass).
    //
    // KNOWN GAP: when a fleet PDB labels the residue "LYS" but
    // structurally has LYN chemistry (e.g. only HZ1/HZ2 present),
    // load-time canonicalisation here cannot know to fire. The
    // rename then needs to happen at protonation-detection time
    // downstream (Phase 1 step 13.2 -- PROPKA wiring at the
    // PdbFileReader TODO; not in this commit's scope).
    //
    // Reference: AMBER ff14SB residue templates encoded in
    // data/ff14sb_params.dat lines 433-453; Markley 1998 §2.1.1.
    AddAtomNameRule("HZ1", "HZ2", "LYN", ToolContext::Standard, ToolContext::Amber);
    AddAtomNameRule("HZ2", "HZ3", "LYN", ToolContext::Standard, ToolContext::Amber);

    // ========================================================================
    // Glycine alpha methylene — HA <-> HA2
    // ========================================================================
    //
    // Glycine has two prochiral alpha hydrogens; AMBER ff14SB +
    // IUPAC convention names them HA2 / HA3 (Markley 1998 §2.1.2).
    // Some fixtures collapse the methylene to a single "HA" name on
    // GLY (pre-Markley convention) when the file format only supplies
    // one of the two atoms. Map "HA" -> "HA2" so the substrate's
    // GLY chain table (which keys HA2 / HA3) finds a match. A third
    // pre-Markley atom labelled "HA1" maps to "HA3".
    //
    // RULES ARE RESIDUE-CONTEXT-KEYED to GLY only: non-Gly residues
    // canonically use plain HA, and aliasing HA -> HA2 there would
    // produce non-existent atoms.
    //
    // Like LYN HZ above, these rules live in the (Standard, Amber)
    // keyspace so Stage 2 of CanonicaliseAmberAtomName fires them
    // and Stage 1 does not chain-rewrite.
    //
    // Reference: AminoAcidType.cpp line 14 (BB_GLY definition uses
    // HA2 / HA3); Markley 1998 §2.1.2.
    AddAtomNameRule("HA",  "HA2", "GLY", ToolContext::Standard, ToolContext::Amber);
    AddAtomNameRule("HA1", "HA3", "GLY", ToolContext::Standard, ToolContext::Amber);

    // ========================================================================
    // N-terminal H2N (CHARMM port) — kept as cap-only literal H2N
    // ========================================================================
    //
    // CHARMM occasionally emits a single H2N pseudo-name standing
    // for the entire NTERM amine NH2-block. AMBER ff14SB N-termini
    // resolve to three SEPARATE atoms H1 / H2 / H3 (charged) or
    // two atoms H1 / H2 (neutral). The mapping is therefore
    // 1-to-many and not a clean per-atom rename.
    //
    // The runtime parser's `IsCapOnlyAtomName` recognises H2N as a
    // cap-only atom name and routes it to LookupCap. We do NOT
    // install a NamingRegistry alias for H2N here — collapsing it
    // would lose information (which of H1/H2/H3 it represents). The
    // canonicalisation is intentionally a no-op for H2N; the
    // substrate composition path treats H2N as a cap-only literal,
    // and downstream re-protonation (when it lands) will produce
    // the correct H1/H2/H3 atom set per terminal_state.
    //
    // FOLLOW-UP: when a load path produces an H2N atom that needs
    // expansion to three H atoms, that has to happen at the load
    // boundary (a fan-out, not a rename). Out of scope for this
    // commit; flagged as residual work.

    // ========================================================================
    // CHARMM ↔ IUPAC COVERAGE GAPS + ONE WILDCARD BUG — DEFERRED pending
    // fleet-wide vetting
    // ========================================================================
    //
    // Context
    // -------
    // A BMRB/RefDB experimental-shift audit on the 10-protein calibration
    // set (h5-reader/notes/nmr_forensics/) surfaced eight categories of
    // /atoms/atom_name ↔ BMRB-IUPAC divergence. Six are missing
    // translation rules; one is a library BUG — the wildcard β-methylene
    // rule fires on ALA's 3-H methyl, producing duplicate HB3 labels at
    // every ALA residue. All findings come from the same audit pass;
    // activate together after fleet-wide vetting. Full rationale in
    // spec/ChangesRequiredBeforeProductionH5Run.md.
    //
    // Findings
    // --------
    //
    // Missing rules — H5 carries CHARMM names where IUPAC is wanted:
    //
    //    Position          CHARMM                IUPAC                 Residues
    //    ----------------- --------------------  --------------------  --------
    //    ILE γ-carbon      CD                    CD1                   ILE
    //    ILE δ-methyl      HD1/HD2/HD3           HD11/HD12/HD13        ILE
    //    ILE γ1-methylene  HG11/HG12             HG12/HG13             ILE
    //    δ-methylene       HD1/HD2               HD2/HD3               ARG, LYS, PRO
    //    ε-methylene       HE1/HE2               HE2/HE3               LYS
    //    α-methylene       HA1/HA2               HA2/HA3               GLY
    //
    // Bug — wildcard β-methylene fires on ALA's 3-H β-methyl:
    //
    //    ALA β-methyl is HB1/HB2/HB3 in BOTH CHARMM and IUPAC — no
    //    translation is needed. But the wildcard rule HB1→HB2,
    //    HB2→HB3 (Standard↔CHARMM) over-fires on ALA, producing:
    //
    //        CHARMM HB1 → IUPAC HB2   (rule applied)
    //        CHARMM HB2 → IUPAC HB3   (rule applied)
    //        CHARMM HB3 → IUPAC HB3   (no rule, passes through)
    //
    //    Result: /atoms/atom_name has two atoms named HB3 and one
    //    named HB2 at every ALA residue. The underlying atom data
    //    (positions, charges, shielding tensors, bond topology) is
    //    correct — three distinct atom indices with parent = CB. Only
    //    the display label is corrupted.
    //
    //    Fix: add ALA-specific identity blocker rules that shadow the
    //    wildcard. ALA HB1/HB2/HB3 → themselves in both directions.
    //    The registry's lookup order (residue-specific wins over
    //    wildcard) causes the block to take precedence.
    //
    // Consequence (pre-existing H5 files on the 10-protein set)
    // ---------------------------------------------------------
    // All eight categories surface as CHARMM names (or corrupted-IUPAC
    // for ALA) in the already-generated H5 files under
    // fleet_calibration-{working,stats,backup}. Regeneration is not on
    // the table. Downstream consumers assuming uniform IUPAC (e.g.,
    // BMRB/RefDB shift binding against /atoms/atom_name) fail across
    // these positions.
    //
    // For the 10-protein experimental-shifts work, each case is handled
    // in-place by two typed Python data constants in
    //     h5-reader/notes/nmr_forensics/pack_experimental_shifts.py
    //
    //   - H5_NAME_FROM_IUPAC       dict lookup (residue, IUPAC) → H5 name
    //                              covers the six missing-rule cases.
    //   - BMRB_METHYL_VIA_PARENT   dict lookup (residue, BMRB name) →
    //                              parent heavy atom; binds BMRB row to
    //                              atom set via /topology/parent_atom_index.
    //                              Handles ALA; flags the binding as
    //                              AMBIGUOUS_METHYL_ORDER. No information
    //                              loss — methyl Hs are chemically
    //                              indistinguishable.
    //
    // Why these rules are COMMENTED OUT
    // ---------------------------------
    // The 10-protein probe showed uniform BMRB/RefDB IUPAC conventions
    // with no per-protein drift, but 10 is a small sample. A fleet-wide
    // probe against the full 685-protein set is a prerequisite for
    // activation: if the probe confirms these are the only registry
    // gaps, uncomment the rules below AND apply the ALA blockers; if
    // additional gaps surface, add them to
    // spec/ChangesRequiredBeforeProductionH5Run.md first, then activate
    // everything together in one coordinated library update before
    // production extraction.
    //
    // Partial activation is forbidden: activating some rules without
    // others risks a partial fix that silently hides a different
    // mismatch behind the accepted correction.
    //
    // References
    // ----------
    //   spec/ChangesRequiredBeforeProductionH5Run.md    full decision record
    //   h5-reader/notes/nmr_forensics/SUMMARY.md        empirical probe
    //
    // Activation recipe (uncomment all of the below together)
    // -------------------------------------------------------
    //
    // ALA β-methyl identity blockers — shadow the wildcard β-methylene
    // rule, which over-fires on ALA's 3-H methyl. CHARMM and IUPAC agree
    // on ALA β-methyl (HB1/HB2/HB3 in both), so identity is correct.
    //
    // AddAtomNameRule("HB1", "HB1", "ALA", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HB2", "HB2", "ALA", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HB3", "HB3", "ALA", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HB1", "HB1", "ALA", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HB2", "HB2", "ALA", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HB3", "HB3", "ALA", ToolContext::Charmm, ToolContext::Standard);
    //
    // ILE γ-carbon + δ-methyl (Standard → CHARMM and CHARMM → Standard)
    //
    // AddAtomNameRule("CD1",  "CD",   "ILE", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HD11", "HD1",  "ILE", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HD12", "HD2",  "ILE", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HD13", "HD3",  "ILE", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("CD",   "CD1",  "ILE", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HD1",  "HD11", "ILE", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HD2",  "HD12", "ILE", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HD3",  "HD13", "ILE", ToolContext::Charmm, ToolContext::Standard);
    //
    // ILE γ1-methylene (CHARMM HG11/HG12 ↔ IUPAC HG12/HG13)
    //
    // AddAtomNameRule("HG12", "HG11", "ILE", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HG13", "HG12", "ILE", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HG11", "HG12", "ILE", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HG12", "HG13", "ILE", ToolContext::Charmm, ToolContext::Standard);
    //
    // δ-methylene — ARG, LYS, PRO (CHARMM HD1/HD2 ↔ IUPAC HD2/HD3)
    //
    // AddAtomNameRule("HD2", "HD1", "ARG", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HD3", "HD2", "ARG", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HD1", "HD2", "ARG", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HD2", "HD3", "ARG", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HD2", "HD1", "LYS", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HD3", "HD2", "LYS", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HD1", "HD2", "LYS", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HD2", "HD3", "LYS", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HD2", "HD1", "PRO", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HD3", "HD2", "PRO", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HD1", "HD2", "PRO", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HD2", "HD3", "PRO", ToolContext::Charmm, ToolContext::Standard);
    //
    // ε-methylene — LYS (CHARMM HE1/HE2 ↔ IUPAC HE2/HE3)
    //
    // AddAtomNameRule("HE2", "HE1", "LYS", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HE3", "HE2", "LYS", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HE1", "HE2", "LYS", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HE2", "HE3", "LYS", ToolContext::Charmm, ToolContext::Standard);
    //
    // α-methylene — GLY (CHARMM HA1/HA2 ↔ IUPAC HA2/HA3)
    //
    // AddAtomNameRule("HA2", "HA1", "GLY", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HA3", "HA2", "GLY", ToolContext::Standard, ToolContext::Charmm);
    // AddAtomNameRule("HA1", "HA2", "GLY", ToolContext::Charmm, ToolContext::Standard);
    // AddAtomNameRule("HA2", "HA3", "GLY", ToolContext::Charmm, ToolContext::Standard);
}


void NamingRegistry::AddAtomNameRule(const std::string& from_name,
                                      const std::string& to_name,
                                      const std::string& residue_name,
                                      ToolContext from_context,
                                      ToolContext to_context) {
    atom_name_map_[{from_name, residue_name, from_context, to_context}] = to_name;
}


// ============================================================================
// Residue name operations
// ============================================================================

bool NamingRegistry::IsKnownResidueName(const std::string& name) const {
    std::string upper = ToUpper(Trim(name));
    return to_canonical_.count(upper) > 0;
}

std::string NamingRegistry::ToCanonical(const std::string& name) const {
    std::string upper = ToUpper(Trim(name));
    auto it = to_canonical_.find(upper);
    if (it != to_canonical_.end())
        return it->second;
    return "";  // unknown -- return empty, caller checks
}

std::string NamingRegistry::ResolveForTool(const std::string& canonical,
                                            ToolContext context,
                                            const std::string& variant) const {
    // Exact (canonical, context, variant) match
    auto it = context_map_.find({canonical, context, variant});
    if (it != context_map_.end())
        return it->second;

    // Fall back to default variant
    if (!variant.empty()) {
        it = context_map_.find({canonical, context, ""});
        if (it != context_map_.end())
            return it->second;
    }

    // No context-specific mapping: return canonical name unchanged
    return canonical;
}


// ============================================================================
// Atom name operations
// ============================================================================

std::string NamingRegistry::TranslateAtomName(const std::string& atom_name,
                                               const std::string& residue_name,
                                               ToolContext from_context,
                                               ToolContext to_context) const {
    if (from_context == to_context) return atom_name;

    // Try residue-specific rule first
    auto it = atom_name_map_.find({atom_name, residue_name, from_context, to_context});
    if (it != atom_name_map_.end())
        return it->second;

    // Try wildcard residue
    it = atom_name_map_.find({atom_name, "*", from_context, to_context});
    if (it != atom_name_map_.end())
        return it->second;

    // No rule: return unchanged
    return atom_name;
}


std::string
NamingRegistry::CanonicaliseAmberAtomName(const std::string& atom_name,
                                           const std::string& residue_name,
                                           ToolContext source_context) const {
    // Stage 1: collapse the wire-format universe to Standard (canonical
    // AMBER ff14SB / IUPAC). For Standard inputs this is a no-op via
    // TranslateAtomName's same-context fast path; for CHARMM inputs
    // this fires the H<->HN, OT1->O, OT2->OXT, HA->HA2 (Gly) etc. rules.
    std::string after_stage1 = (source_context == ToolContext::Standard)
            ? atom_name
            : TranslateAtomName(atom_name, residue_name,
                                 source_context, ToolContext::Standard);

    // Stage 2: residue-context-keyed canonicalisation that applies
    // even on Standard inputs (e.g. LYN HZ1->HZ2 / HZ2->HZ3 fleet
    // variance, GLY HA->HA2 collapsed-methylene). Implemented as
    // (Standard -> Amber) rules so the existing residue-keyed
    // dispatch in TranslateAtomName fires; the result IS canonical
    // AMBER ff14SB.
    return TranslateAtomName(after_stage1, residue_name,
                              ToolContext::Standard, ToolContext::Amber);
}


// ============================================================================
// Global singleton
// ============================================================================

NamingRegistry& GlobalNamingRegistry() {
    static NamingRegistry instance;
    return instance;
}

}  // namespace nmr

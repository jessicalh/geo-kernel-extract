# V2 Library Bugs Found During Unification (2026-03-31)

## Bug 1: PdbFileReader destroys histidine protonation state (CRITICAL)

**Symptom**: All histidine rings classified as HIS. HID/HIE never produced.
Old v1 output: 168K atoms with HIS contributions, 20K with HIE, 0 HID.
New v2 output: 188K atoms with HIS (= 168K + 20K merged), 0 HIE, 0 HID.

**Impact**: Every v2 calculation on histidine-containing proteins is wrong:
- Ring current intensities differ for HIS vs HID vs HIE
- Feature extraction produces wrong per-type sums
- Heuristic tier classification affected
- Viewer displays incorrect ring types

**Location**: src/Protein/PdbFileReader.cpp -- histidine residue name normalization.
The v2 reader normalizes HID/HIE/HIP to HIS. The v1 reader preserved them.

**Fix needed**: Preserve the PDB residue name (HID/HIE/HIS) through to AromaticRing::ParentAminoAcid and RingTypeIndex. The RingType table already has separate entries for HIS/HID/HIE.

## Bug 2: MolecularGraph BFS distances differ (HIGH)

**Symptom**: graph_dist_ring differs by up to 111 between v1 and new port.
graph_dist_N differs by 99, graph_dist_O by 97.

**Likely cause**: Bond detection threshold, hydrogen handling, or ring atom source set differs between the ported MolecularGraph and v1 original. The ported code may use different covalent radii or a different tolerance factor.

**Location**: src/Features/V1Physics.cpp -- MolecularGraph::detect_bonds() and multi_source_bfs().

**Fix needed**: Compare bond detection line-by-line against v1 model/MolecularGraph.cpp.

## Bug 3: Coulomb E-field / charge loading differs (MEDIUM)

**Symptom**: coulomb_E_mag max_diff=0.33, coulomb_E_backbone_frac max_diff=5.0.

**Likely cause**: Different partial charges loaded. v1 loaded from ff14sb_params.dat with its own matching logic. v2 uses ChargeAssigner or Conformation charges. The V1Adapter charge conversion may not match.

**Location**: src/Features/V1Adapter.cpp charge handling.

## Bug 4: T2 target normalization (FIXED)

**Symptom**: T2 off-diagonal components differed by sqrt(3/2). T0 and T2 m=0 were correct.

**Cause**: New code read T2 from JSON "spherical" field (unnormalized). Old binary re-decomposed from delta_tensor with CG normalization.

**Fix**: Changed to re-decompose using v1::decompose_tensor. COMMITTED.

## Verification data
- Old features at: data/features/ (produced by old v1 binary, now in bones)
- New features at: data/features-unified/ (produced by new unified binary)
- Comparison: 63 features perfect, 50 close (<1e-3), 23 significantly different

# Illustrative test peptides — build instructions

Per `spec/PLANNED_CALCULATORS_2026-04-22.md` Amendment 2026-05-08(c).
Two reference structures used by the calculator-walkthrough Mathematica
workbook chapters: a folded mini-protein (Trp-cage) for citation-friendly
illustrations and a designed all-canonical peptide (synthetic 22-mer)
for chemistry features Trp-cage doesn't reach.

## Files

- `1l2y.pdb` — raw NMR ensemble from RCSB, 38 models. Source of
  truth; preserved for ensemble-aware analyses if needed later.
- `trp_cage_1l2y_model1.pdb` — first NMR model only, extracted from
  `1l2y.pdb`. The single-conformation file the workbooks load.
- `synthetic_22mer.pdb` — built via `tleap` from the script
  `build_synthetic.tleap`. 22 residues, sequence `ACDEFGHIKLMNPQRSTVWYCV`
  (one of each canonical AA + Cys-21 paired with Cys-2 via SS bond
  + Val-22 for QG super-group illustration). 349 atoms in extended
  conformation.
- `build_synthetic.tleap` — the deterministic tleap input that
  produces `synthetic_22mer.pdb`. ff14SB force field, CYX residues
  for the disulfide pair, HIE for histidine.
- `aimnet2_requires_grad_check.py` — pre-flight investigation for
  the planned AIMNet2PolarisabilityResult calculator slice (per
  PLANNED_CALCULATORS Amendment 2026-05-08(b)). Verifies that
  `data/models/aimnet2_wb97m_0.jpt` propagates gradients through
  the coordinate input tensor. Result on 2026-05-09: **PASSED**
  (autograd path is viable; the polarisability calculator can be
  built without re-exporting the model). See script docstring for
  full details.

## Regenerating

### Trp-cage

```bash
curl -s -o 1l2y.pdb 'https://files.rcsb.org/download/1L2Y.pdb'
awk '/^MODEL / && model_count++ > 0 { exit } { print }' 1l2y.pdb > trp_cage_1l2y_model1.pdb
```

Deterministic — RCSB serves the same PDB from the canonical entry.
The first MODEL block is what becomes the single-conformation file.

### Synthetic 22-mer

```bash
$AMBERHOME/bin/tleap -f build_synthetic.tleap
# or
/home/jessica/micromamba/envs/mm/bin/tleap -f build_synthetic.tleap
```

Deterministic — tleap applies ff14SB ideal angles and the standard
extended-chain phi/psi to the sequence. Re-running produces a
byte-identical PDB modulo timestamp comments.

### AIMNet2 requires_grad check

```bash
python3 aimnet2_requires_grad_check.py
```

Exit code 0 = pass (autograd path viable); 1 = .jpt blocks gradients
(need re-export); 2 = setup error (missing model or load failure);
3 = forward pass error.

## Caveats

- The synthetic 22-mer is a **chemistry-coverage substrate**, not a
  fold. It does not represent a stable conformation; the workbook
  chapters consume the geometry as-is to illustrate calculator
  per-atom outputs. If you need a folded reference, use Trp-cage.
- The disulfide bond between Cys-2 and Cys-21 is wired by the
  `bond mol.2.SG mol.21.SG` line in the tleap script; it requires
  `CYX` in the sequence (not `CYS`) so the residues have the
  disulfide-bonded atom types (S not SH; no HG). Switching to `CYS`
  produces tleap errors about missing SH-SH parameters.
- Variants (HID/HIE/HIP, ASH, GLH, LYN, ARN, TYM) are **not** all
  represented in the synthetic peptide — only HIE for histidine,
  CYX for cysteine, neutral defaults for everything else. To
  illustrate variant chemistry, generate a parallel peptide with the
  swapped variant and compare side-by-side. The categorical record
  in `atoms_category_info.npy` makes per-variant joins trivial.

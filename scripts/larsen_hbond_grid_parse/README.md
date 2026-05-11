# Larsen H-bond DFT grid parser

Parses Larsen 2015 ProCS15 H-bond DFT grid logs from the ERDA archive
into 6 (NPZ + JSON) pairs under `data/larsen_hbond_grids/` consumed
by `LarsenHBondShieldingResult` at runtime.

## Inputs

- `/mnt/expansion/larsen_archive/hydrogenbondnmrlogs.tar` — outer tar
  containing 6 nested `.tar.bz2` archives (NMANMA, NMACOH, NMACOO,
  ALANMA, ALACOH, ALACOO). Each archive holds Gaussian 09 NMR logs
  (GIAO OPBE/6-311++G(2d,p) on PM6-frozen monomers) scanned over
  (rOH, θ, ρ) grid points. Total ~21K logs, ~500 MB compressed.

## Outputs

Per archive `<STEM>` (NMANMA, NMACOH, NMACOO, ALANMA, ALACOH, ALACOO):

- `data/larsen_hbond_grids/<STEM>_grid.npz` — numpy zip with:
  - `r`, `theta`, `rho` (shape (N,)) — actual measured geometry per grid point.
  - `donor_<X>` for X ∈ {N, CA, CB, C, HA, HN} (shape (N, 3, 3)) —
    Δσ tensor on donor atom X, expressed in canonical donor frame.
  - `acceptor_<X>` for X ∈ {N, HN, HA} (NMA acceptor only, shape (N, 3, 3))
    — Δσ tensors representing residue i+1's backbone atoms.
- `data/larsen_hbond_grids/<STEM>_meta.json` — axis ranges, readout
  atom names, canonical atom indices (1-based, matches Gaussian log
  numbering).

## Canonical donor frame

For each grid point's tensors, the parser rotates from Gaussian's
Standard orientation frame into a canonical donor frame defined by:

- Origin: donor H (Hα for ALA donor, amide H for NMA donor).
- z-axis: donor anchor → donor H direction. For ALA: Cα→Hα bond.
  For NMA: N→H bond.
- x-axis: donor third atom (Cα-side N for ALA, carbonyl C for NMA)
  projected orthogonal to z.
- y-axis: z × x.

At runtime, the calculator computes the same canonical frame from
protein atom positions, looks up the σ tensor in the grid, rotates
back to protein frame via the protein-side rotation.

## Reference subtraction (proxy)

The Larsen H-bond archive does NOT contain free-monomer reference σ
calculations. Larsen's pipeline performed those separately and stored
the resulting Δσ in pre-interpolated `.npy` files (probably in
`predictions.tar.bz2` or `proteinnmrlogs.tar.bz2`, both unfetched
from ERDA).

This parser uses the **r-max grid edge as a proxy reference**: per
readout atom, reference σ = mean over grid points where r is within
one grid step of r_max. Subtracted from every grid point's σ to give
approximate Δσ. Bias is small (~0.1 ppm; H-bond effect at r ≥ 3.8 Å
for ALA donor and r ≥ 2.875 Å for NMA donor is near-zero per the
Larsen Table 1 RMSD scale) but documented.

To replace with actual free-monomer σ later: run Gaussian 09 GIAO
OPBE/6-311++G(2d,p) on isolated PM6-optimised NMA + Ac-A-NMe, swap
the parser's `compute_reference_and_subtract` to subtract those
values instead.

## Filename ↔ geometry caveats

Verified empirically across 10 sampled logs (2026-05-11):

- **NMA donor archives**: r_filename ≈ r_actual (offset ≤ 0.001 Å).
  θ_filename = θ_actual at the O vertex (Hα–O–C(=O) angle).
  ρ_filename = −ρ_actual (sign-flip; Larsen's dihedral convention
  opposite to standard IUPAC).
- **ALA donor archives**: r_filename = r_actual + 0.2 Å (consistent
  offset of one grid step — likely off-by-one in Larsen's input-file
  labeling). θ matches as above. ρ sign-flipped as above.

Parser keys grid on **actual measured (r, θ, ρ)**, not filename
values. Filenames are labels only.

## Known limitations

1. **Scattered grid**: PM6-optimisation of monomers introduces small
   geometric scatter, so the rho axis has ~100-300 distinct values
   rather than the nominal 24-25. The C++ grid loader bins to nearest
   nominal axis value at load time, or does scattered nearest-N
   interpolation.

2. **Close-approach extremes**: at r near 1.5 Å (NMA donor minimum)
   some grid points have large |Δσ| (up to ~600 ppm on N) where DFT
   produces anomalous values due to steric strain. Larsen 2015
   smoothed these via interpolation. The calculator may want to clip
   or de-weight near-collision points.

3. **Larsen-failed DFT logs**: Some grid points failed Gaussian
   optimization ("Error termination via Lnk1e ... Small interatomic
   distances encountered"); these are skipped at parse time and
   filled by runtime trilinear interpolation from neighbouring grid
   points.

4. **Duplicate filenames**: Larsen's archive contains some duplicate
   geometries (same r, θ, ρ from different filename labels). Parser
   detects and skips duplicates by rounded geometry key.

## Usage

```bash
# Parse all 6 archives, default output location
python3 scripts/larsen_hbond_grid_parse/parse_larsen_hbond_grids.py

# Specific archive
python3 scripts/larsen_hbond_grid_parse/parse_larsen_hbond_grids.py --archive ALANMA

# Smoke test (N logs)
python3 scripts/larsen_hbond_grid_parse/parse_larsen_hbond_grids.py --archive ALANMA --max-logs 5

# Custom output
python3 scripts/larsen_hbond_grid_parse/parse_larsen_hbond_grids.py --out /custom/path/
```

Full run: ~35 seconds for all 6 archives on a workstation.

## References

- Larsen, Bratholm, Christensen, Hamelryck, Boomsma & Jensen,
  *PeerJ* 3:e1344 (2015), DOI 10.7717/peerj.1344. §H-bond scans
  paragraph defines the parameterization; Table 2 specifies which
  atom types receive which contributions.
- `references/larsen-2015-procs15-dft-chemical-shift-predictor.pdf`
  for the local copy.
- `spec/plan/larsen-hbond-shielding-design-2026-05-11.md` for the
  authoritative design document.

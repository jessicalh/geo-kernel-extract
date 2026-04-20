# RefDB snapshot — 2016 vintage

The Re-referenced Protein Chemical Shift Database (RefDB, Zhang/Neal/Wishart
2003, *J. Biomol. NMR* 25:173-195) provides SHIFTX/SHIFTCOR-derived corrections
to BMRB chemical shift depositions.

## Provenance

- **Source:** `https://refdb.wishartlab.com/`
- **Fetched:** 2026-04-20 (this directory)
- **Server state at fetch time:** Apache 2.2.15 / CentOS; expired SSL cert;
  homepage `Last-Modified: Thu, 14 Jul 2016 16:15:20 GMT`. **RefDB has not
  been updated since 2016.** The chemical-shift corrections reflect the
  SHIFTX prediction model and BMRB data as of that snapshot date.
- **Mapping file `bmrpdbnew.txt`:** `Last-Modified: Thu, 15 Dec 2011` —
  the BMRB→PDB index is from 2011; later BMRB entries (post-2011) will
  not appear in the mapping, though some later entries may still have
  corrdata files in the tarball.

## Contents

- `RefDB_files.tar.gz` (30 MB) — full 2016 RefDB snapshot containing
  `data/www/lib/RefDB/RefDB/corrdata/bmr<id>.str.corr` files for
  2427 entries (NMR-STAR 2.1.1 format).
- `bmrpdbnew.txt` (350 KB) — tab-delimited per-entry summary with
  columns: `BMRB_ID`, `PDB_ID`, `RESOLUTION`, `SEQ_IDENTITY`,
  `SEQ_LENGTH`, `SEQ_OFFSET`, nucleus-present flags (`HA`, `CA`, `CB`,
  `CO`, `N`, `HN`), average correlation, per-nucleus correlation
  (`*_CORR`), RMSD (`*_RMSD`), chemical-shift delta (`*_CSDIFF`).

## Hashes

See the matching `provenance.toml` blocks in each per-protein
`experimental_shifts/` bundle for the tarball SHA256 at extraction
time. Present hashes at snapshot time in this directory's own
provenance record (to be written alongside the per-protein
extractions).

## Usage

The snapshot is the authoritative RefDB source for this project.
Per-protein `bmr<id>.str.corr` files are extracted on demand by
`parse_shifts.py` (sibling directory). The full tarball is not
committed to git; re-download via `curl -ksL -o RefDB_files.tar.gz
https://refdb.wishartlab.com/RefDB_files.tar.gz` if missing.
`bmrpdbnew.txt` is small enough to commit.

## Caveats for the Huxley-discipline reader

RefDB itself drops some BMRB rows that SHIFTCOR could not reconcile
against the structural model (e.g., for 1DV0/4757: BMRB raw has 309
H shifts; RefDB has 261). **RefDB is a curated re-referenced subset,
not a complete mirror of BMRB.** Downstream work needs both sources
separately: BMRB raw for complete experimental record; RefDB for the
re-referencing correction metadata.

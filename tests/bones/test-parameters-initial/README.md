# Baseline Parameter Extraction — Pre-TOML

Captured 2026-04-07 at commit 6a39328 (after doc reframing, before
TOML-ising constants). All parameters are literature defaults hardcoded
in ring type classes and calculator .cpp files.

## Protein: A0A7C5FAR6

543 atoms (WT), 501 atoms (ALA), 35 residues. Run via:

    build/nmr_extract --mutant \
      --wt-xyz .../A0A7C5FAR6_WT.xyz --wt-prmtop .../A0A7C5FAR6_WT.prmtop \
      --wt-nmr .../A0A7C5FAR6_WT_20260327_120049_nmr.out \
      --ala-xyz .../A0A7C5FAR6_ALA.xyz --ala-prmtop .../A0A7C5FAR6_ALA.prmtop \
      --ala-nmr .../A0A7C5FAR6_ALA_20260327_122543_nmr.out \
      --output test-parameters-initial/npy

## Contents

- npy/              — 50 NPY files (full WT mutant comparison output)
- manifest.json     — per-file shape, dtype, md5, min/max
- extraction_log.jsonl — 212 JSON log lines (OperationLog)
- extraction_stdout.txt — non-JSON console output

## How to verify after TOML changes

```python
import numpy as np, json

with open('test-parameters-initial/manifest.json') as f:
    baseline = json.load(f)

for name, info in baseline['files'].items():
    old = np.load(f'test-parameters-initial/npy/{name}')
    new = np.load(f'new-output-dir/{name}')
    assert old.shape == new.shape, f'{name}: shape mismatch'
    if old.dtype.kind == 'f':
        assert np.allclose(old, new, atol=1e-12), f'{name}: values differ'
    else:
        assert np.array_equal(old, new), f'{name}: values differ'
print('All 50 files match.')
```
